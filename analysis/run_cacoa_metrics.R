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
getESM <- function(cao.obj, min.cells=10, n.cells=1e3, dist='cor', n.subsamples=50){
  res <- cao.obj$estimateExpressionShiftMagnitudes(min.cells=min.cells, n.cells=n.cells, dist=dist, n.subsamples=n.subsamples)
  dist.per.type <- res$dist.df %$% split(value, Type)
  return(dist.per.type)
}

getCESM <- function(cao.obj){
  res <- cao.obj$estimateCommonExpressionShiftMagnitudes(n.cores=20)
  cell.types <- names(res[[1]])
  dist.per.type <- lapply(setNames(cell.types, cell.types), function(n) {
    lapply(res, `[[`, n) %>% do.call(rbind, .) %>% colMeans()
  })
  return(dist.per.type)
}

getCFES <- function(cao.obj, min.expr.frac=0.01){
  res <- cao.obj$estimateClusterFreeExpressionShifts(n.top.genes=1000)
  dist.per.type <- res %>% split(cao.obj$sample.groups[names(.)])
  return(dist.per.type)
}

#test_cesm <- getCESM(cacoa_objects$`687`$ncell$cao)
#test_esm <- getESM(cacoa_objects$`687`$ncell$cao)

# convenience function to loop over the datasets contained in the list of datasets
lapplypersublist <- function(items.per.seed, function_name, ...){
  results <- items.per.seed %>% lapply(function(x){
    #print(names(x))
    get(function_name)(x$cao, ...)
  })
}

# we only need to calculate this value once per cacoa object
lapplyCFZ <- function(items.per.seed){
  items.per.seed %>% lapply(function(x){
    x$cao$estimateClusterFreeZScores(min.expr.frac=0.01)
  })
  return(NULL)
}

cacoa_objects %>% lapply(lapplyCFZ)

esms <- cacoa_objects %>% lapply(lapplypersublist, 'getESM')
cesms <- cacoa_objects %>% lapply(lapplypersublist, 'getCESM')
cfes <- cacoa_objects %>% lapply(lapplypersublist, 'getCFES')

test_cao <- cacoa_objects$`687`$ncell$cao
test_cao$estimateClusterFreeZScores(min.expr.frac=0.01) # This one you need just once. Then if you change parameters of expression shifts, you don't need to re-estimate z-scores
res_test <- test_cao$estimateClusterFreeExpressionShifts(n.top.genes=1000)
test_dist_per_type <- res_test %>% split(test_cao$sample.groups[names(.)])

# get the medians from each result per seed to compare
extractmedianesm <- function(esm.result, varied.factor){
  df <- esm.result %>% as_tibble()
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


### expression shift magnitudes
esm_ncell_medians <- esms %>% lapply(function(x){
  extractmedianesm(x$ncell, 'ncell')
}) %>% dplyr::bind_rows()

esm_ngenes_medians <- esms %>% lapply(function(x){
  extractmedianesm(x$ngenes, 'ngenes')
}) %>% dplyr::bind_rows()

esm_ncell_medians %>% ggplot(aes(x=ncell, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

esm_ngenes_medians %>% ggplot(aes(x=ngenes, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

esms_ncell_all <- esms %>% lapply(function(x){x$ncell}) %>% dplyr::bind_rows() %>%
  pivot_longer(cols=everything(), names_to='ncell_defactor', values_to='esm') %>%
  mutate(ncell_defactor=factor(ncell_defactor, levels=stringr::str_sort(unique(ncell_defactor), numeric=TRUE)))
esms_ncell_all %>% ggplot(aes(x=ncell_defactor, y=esm))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90, hjust = 1))


### common expression shift magnitude

cesm_ncell_medians <- cesms %>% lapply(function(x){
  extractmedianesm(x$ncell, 'ncell')
}) %>% dplyr::bind_rows()

cesm_ngenes_medians <- cesms %>% lapply(function(x){
  extractmedianesm(x$ngenes, 'ngenes')
}) %>% dplyr::bind_rows()

cesm_ncell_medians %>% ggplot(aes(x=ncell, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

cesm_ngenes_medians %>% ggplot(aes(x=ngenes, y=median.result, col=de.level))+geom_point()+theme(legend.position="top")

cesms_ncell_all <- cesms %>% lapply(function(x){x$ncell}) %>% dplyr::bind_rows() %>%
  pivot_longer(cols=everything(), names_to='ncell_defactor', values_to='cesm') %>%
  mutate(ncell_defactor=factor(ncell_defactor, levels=stringr::str_sort(unique(ncell_defactor), numeric=TRUE)))
cesms_ncell_all %>% ggplot(aes(x=ncell_defactor, y=cesm))+geom_boxplot()+
  theme(axis.text.x = element_text(angle = -90, hjust = 1))


###
