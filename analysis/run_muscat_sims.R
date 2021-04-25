library(conos)
library(tidyverse)
library(splatter)
library(muscat)
library(SingleCellExperiment)
library(optparse)
library(parallel)
devtools::load_all('/home/larsc/SecretUtils')

option_list = list(
  make_option(c("--de_type"), type="character", default="de",
              help="which type of de to use, must be either de, dp, dm or db, [default= %default]", metavar="character"),
  make_option(c("--n_steps"), type="integer", default=5,
              help="how many DE steps, [default= %default]", metavar="integer"),
  make_option(c("--n_samples"), type="integer", default=2,
              help="number of samples for simulation, [default= %default]", metavar="integer"),
  make_option(c("--output", "-o"), type="character", default="/home/larsc/data",
              help="output folder, [default= %default]", metavar="path")

)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ref <- readRDS(file.path('/home/larsc/data/lamp5_prepsim.rds'))

de_type_to_idx <- c(3,4,5,6) %>% setNames(c('de', 'dp', 'dm', 'db'))
de_type_idx <- de_type_to_idx[opt$de_type]
n_steps <- opt$n_steps


# we will keep EP at 0 for now, and make simulations with increasing de proportion
generateDEvectors <- function(de.type.idx, n.steps){
  de.fractions <- 1:n.steps %>% sapply(function(x){x*1e-1})
  de.vectors <- de.fractions %>% lapply(function(de.fraction){
    init.vec <- diag(6)[1, ]
    init.vec[1] <- 1-de.fraction
    init.vec[de.type.idx] <- de.fraction
    return(init.vec)
  })
  return(de.vectors)
}

de_proportions <- generateDEvectors(de_type_idx, n_steps)


generateSims <- function(sce.obj, ncell, ngenes, nreplicates, varied.factor, de.proportions, n.cores){
  sims <- de.proportions %>% mclapply(function(prop.vec){
    sim.obj <- muscat::simData(
      sce.obj, p_dd = prop.vec, nk = 1, ns = nreplicates, nc = ncell*2, ng = ngenes, force = TRUE)
    sim.obj@colData$de <- prop.vec[3]
    sim.obj@colData[[varied.factor]] <- get(varied.factor)
    sim.obj@colData$cellid <- paste(rownames(sim.obj@colData), prop.vec[3], get(varied.factor), sep='_')
    colnames(sim.obj@assays@data$counts) <- sim.obj@colData$cellid
    return(sim.obj)
  }, mc.cores=n.cores)
  sim.names <- de.proportions %>% lapply(function(x){x[x>0][2]}) %>% unlist
  names(sims) <- sim.names
  return(sims)
}

extractSimItems <- function(sims.by.de){
  cm.bound <- sims.by.de %>% lapply(function(sim.obj){
    Matrix::t(sim.obj@assays@data$counts)
  }) %>% do.call(rbind, .)
  annot.bound <- sims.by.de %>% lapply(function(sim.obj){
    df <- data.frame(sim.obj@colData)
    rownames(df) <- NULL
    return(df)
  }) %>% dplyr::bind_rows()
  gene.infos <- sims.by.de %>% lapply(function(sim.obj){
    sim.obj@metadata$gene_info
  })
  return(list(cm = cm.bound, annot = annot.bound, gene.infos = gene.infos))
}

SimAndExtract <- function(sce.obj, ncell, ngenes, nreplicates, varied.factor, de.proportions, n.cores){
  if (varied.factor=='ncell'){
    sims.per.factor <- ncell %>% lapply(function(ncell.val){
      sim.items <- generateSims(sce.obj, ncell.val, ngenes, nreplicates, varied.factor,
                                de.proportions, n.cores) %>% extractSimItems()
    })
    cm <- sims.per.factor %>% lapply(function(x){x$cm}) %>% do.call(rbind, .)
  } else if (varied.factor=='ngenes'){
    sims.per.factor <- ngenes %>% lapply(function(ngene.val){
      sim.items <- generateSims(sce.obj, ncell, ngene.val, nreplicates, varied.factor,
                                de.proportions, n.cores) %>% extractSimItems()
    })
    cms <- sims.per.factor %>% lapply(function(x){x$cm})
    most.genes.idx <- which(unlist(lapply(cms, ncol))==max(unlist(lapply(cms, ncol))))
    total.genes <- cms[[most.genes.idx]] %>% colnames()
    cm <- cms %>% lapply(Matrix::t) %>% lapply(SecretUtils::PadGenesRows, total.genes) %>%
      lapply(Matrix::t) %>% do.call(rbind, .)
  }
  annot <- sims.per.factor %>% lapply(function(x){x$annot}) %>% dplyr::bind_rows()
  gene.infos <- sims.per.factor %>% lapply(function(x){x$gene.infos})
  return(list(cm = cm, annot = annot, gene.infos = gene.infos))
}


ncell_values <- c(50, 100, 200, 500)
ngenes_values <- c(500, 1000, 2000, 5000)
seeds <- c(687, 788, 894, 273, 795)

# choose 500 cells and 5000 genes as default values for now
datasets_per_seed <- seeds %>% lapply(function(seed, n.samples=opt$n_samples){
  print(paste0('creating datasets for seed: ', as.character(seed)))
  set.seed(seed)
  sims.ncell <- SimAndExtract(ref, ncell_values, 5000, n.samples, 'ncell', de_proportions, 5)
  sims.ngenes <- SimAndExtract(ref, 500, ngenes_values, n.samples, 'ngenes', de_proportions, 5)
  return(list(ncell=sims.ncell, ngenes=sims.ngenes))
})


createCacoa <- function(simulated.dataset, varied.factor){
  sim.annot <- simulated.dataset$annot
  sim.annot %<>% mutate(celltype=paste0(get(varied.factor), '_', de))
  sim.cm <- simulated.dataset$cm
  cms.by.sample <- sim.annot$cellid %>% split(sim.annot$sample_id) %>% lapply(function(x){sim.cm[x, ]})
  sim.p2s <- lapply(cms.by.sample %>% lapply(Matrix::t), pagoda2::basicP2proc, n.cores=5, min.cells.per.gene=0, n.odgenes=2e3,
                    get.largevis=FALSE, make.geneknn=FALSE, get.tsne=FALSE, min.transcripts.per.cell=1)
  sim.conos <- Conos$new(sim.p2s)
  sim.conos$buildGraph(n.odgenes=2e3) # crashes if the value is higher than actual genes present
  sim.conos$embedGraph(method='largeVis')
  sim.samples.con <- setNames(sim.annot$sample_id, sim.annot$cellid)
  sim.condition.con <- setNames(sim.annot$group_id, sim.annot$cellid)
  sim.celltype.con <- setNames(sim.annot$celltype, sim.annot$cellid)
  conos.sample.plot <- sim.conos$plotGraph(groups = sim.samples.con, show.legend=T)
  conos.celltype.plot <- sim.conos$plotGraph(groups = sim.celltype.con, show.legend=T)
  conos.condition.plot <- sim.conos$plotGraph(groups = sim.condition.con, show.legend=T)
  # cacoa
  sim.sample.grps <- sim.annot %>% split(sim.annot$sample_id) %>% lapply(function(x){x$group_id %>% unique}) %>% unlist()
  sim.cell.grps <- setNames(sim.annot$celltype, sim.annot$cellid)
  cao <- cacoa::Cacoa$new(sim.conos, sample.groups=sim.sample.grps, cell.groups=sim.cell.grps,
                          n.cores=20, target.level='B', ref.level='A')
  return(list(conos.sample.plot=conos.sample.plot, conos.celltype.plot=conos.celltype.plot,
              conos.condition.plot=conos.condition.plot, cao=cao))
}

cacoas_per_seed <- datasets_per_seed %>% lapply(function(seed){
  cacoa.objs <- names(seed) %>% lapply(function(varied.factor){
    createCacoa(seed[[varied.factor]], varied.factor)
  })
  names(cacoa.objs) <- names(seed)
  return(cacoa.objs)
})
names(cacoas_per_seed) <- seeds
saveRDS(cacoas_per_seed, paste0(opt$output, '/cacoa_objects_type_', opt$de_type,'.rds'))


### for testing/sanity checking
#cacoas_per_seed <- readRDS('/d0/home/larsc/data/cacoa_objects_type_de.rds')
#res <- cacoas_per_seed$`687`$ncell$cao$estimateExpressionShiftMagnitudes(min.cells=5, dist="cor", n.subsamples=50)
#dist.per.type <- res$dist.df %$% split(value, Type)
#dist.df <- dist.per.type %>% data.frame
#names(dist.df) <- names(dist.per.type) # stupid appending of characters when casting to df
#dist.df <- dist.df %>% tidyr::pivot_longer(cols=everything(), names_to='celltype', values_to='distance')
#ncell_de_list <- strsplit(dist.df$celltype, '_')

#dist.df %<>% mutate(ncell = as.factor(as.numeric(ncell_de_list %>% lapply(function(x){x[[1]]}) %>% unlist)), de = ncell_de_list %>% lapply(function(x){x[[2]]}) %>% unlist,)
#dist.df %>% ggplot(aes_string(x='ncell', y='distance', col='de'))+
#  geom_jitter()+ theme(legend.position="top")

#test_cm <- cacoas_per_seed$`687`$ncell$cao$getJointCountMatrix(raw=TRUE)
#test_cm2 <- datasets_per_seed[[1]]$ncell$cm
#all.equal(colSums(test_cm2), colSums(test_cm))
