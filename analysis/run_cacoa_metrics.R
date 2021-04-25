library(conos)
library(tidyverse)
library(splatter)
library(muscat)
library(SingleCellExperiment)
library(optparse)
library(parallel)
devtools::load_all('/home/larsc/SecretUtils')

option_list = list(
  make_option(c("--input", "-i"), type="character", default="/home/larsc/data/cacoa_objects_type_de.rds",
              help="output folder, [default= %default]", metavar="path")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

cacoa_objects <- readRDS(opt$input)

# just a wrapper for generating the expression shift results
getEESM <- function(cao.obj, min.cells=10, n.cells=1e3, dist='cor', n.subsamples=50){
  res <- cao.obj$estimateExpressionShiftMagnitudes(min.cells=min.cells, n.cells=n.cells, dist=dist, n.subsamples=n.subsamples)
  dist.per.type <- res$dist.df %$% split(value, Type)
  return(dist.per.type)
}

# convenience function to loop over the datasets contained in the list of datasets
lapplypersublist <- function(items.per.seed, function_name, ...){
  results <- items.per.seed %>% lapply(function(x){
    print(names(x))
    get(function_name)(x$cao, ...)
  })
}

eesms <- cacoa_objects %>% lapply(lapplypersublist, 'getEESM')

# get the medians from each result per seed to compare
extractmedianeesm <- function(eesm.result, varied.factor){
  df <- eesm.result %>% as_tibble()
  summary.medians <- df %>% names %>% lapply(function(col.name){
    median.result <- df[[col.name]] %>% median()
    de.level <- gsub('.*_', '', col.name)
    assign(varied.factor, gsub('_.*', '', col.name))
    result.list <- list(median.result, de.level, get(varied.factor))
    names(result.list) <- c('median.result', 'de.level', varied.factor)
    return(result.list)
  })
  df.out <- dplyr::bind_rows(summary.medians)
  df.out <- df.out %>% mutate(de.level=as.factor(as.double(de.level)))
  df.out[[varied.factor]] <- as.factor(as.integer(df.out[[varied.factor]]))
  return(df.out)
}

eesm_ncell_medians <- eesms %>% lapply(function(x){
  extractmedianeesm(x$ncell, 'ncell')
}) %>% dplyr::bind_rows()

eesm_ngenes_medians <- eesms %>% lapply(function(x){
  extractmedianeesm(x$ngenes, 'ngenes')
}) %>% dplyr::bind_rows()

eesm_ncell_medians %>% ggplot(aes(x=ncell, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

eesm_ngenes_medians %>% ggplot(aes(x=ngenes, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

eesms_ncell_all <- eesms %>% lapply(function(x){x$ncell}) %>% dplyr::bind_rows() %>%
  pivot_longer(cols=everything(), names_to='ncell_defactor', values_to='EESM') %>%
  mutate(ncell_defactor=factor(ncell_defactor, levels=stringr::str_sort(unique(ncell_defactor), numeric=TRUE)))
eesms_ncell_all %>% ggplot(aes(x=ncell_defactor, y=EESM))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90, hjust = 1))
