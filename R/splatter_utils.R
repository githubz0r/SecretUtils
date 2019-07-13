Selectsimcons <- function(connectivity.mat, de.prob.vec, n.batch.vec) {
  #browser()
  connectivity.mat %<>% as.matrix()
  n.bv <- length(n.batch.vec)
  n.pv <- length(de.prob.vec)
  ind.intervals <- 1:n.bv %>% lapply(function(x){((x-1)*n.pv+1):(x*n.pv)})
  n.intervals <- de.prob.vec %>% unique %>% length
  n.de.reps <- table(de.prob.vec) %>% unique
  interval.start.inds <- seq(1, n.pv, n.de.reps)
  interval.end.inds <- c(interval.start.inds[-1]-1, n.pv)
  interval.all.inds <- Map(function(x,y){x:y}, interval.start.inds, interval.end.inds)
  obtainsubintervals <- function(interval){
    reference <- interval[interval.all.inds[[1]]]
    comps <- 1:n.intervals %>% lapply(function(i){interval[interval.all.inds[[i]]]})
    submats <- comps %>% lapply(function(x){connectivity.mat[reference,x]})
  }
  comp.sections <- ind.intervals %>% lapply(obtainsubintervals)
  comp.sections <- comp.sections %>% lapply(function(x){names(x)<-unique(de.prob.vec); return(x)})
  names(comp.sections) <- n.batch.vec
  meltMats <- function(list.of.mats){
    diag(list.of.mats[[1]]) <- NA
    list.of.mats[[1]][upper.tri(list.of.mats[[1]])] <- NA
    molten.mats <- list.of.mats %>% lapply(reshape2::melt) %>% lapply(na.omit)
    molten.dfs <- 1:length(molten.mats) %>% lapply(function(i){y <- molten.mats[[i]]; y$de.prob=names(list.of.mats)[[i]];
    return(y)})
    return(molten.dfs)
  }
  comp.dfs <- comp.sections %>% lapply(function(x){meltMats(x)})
  return(comp.dfs)
}

SimulateGroups <- function(splatter.params, n.cells, de.probabilities = de_prob, group.prob = group_prob,
                           n.genes, lib.loc, mean.shape=0.414, seed=22071) {
  sim.result <- splatSimulateGroups(splatter.params, group.prob = group.prob, dropout.type='experiment', mean.shape=mean.shape,
                                    de.prob = de.probabilities, nGenes=n.genes, batchCells=n.cells*10, lib.loc=lib.loc,
                                    verbose = FALSE, seed=seed)
  #browser()
  sim.annot <- sim.result@colData %>% as.data.frame
  nrowann <- dim(sim.annot)[1]
  sim.annot %<>% mutate(ncell=rep(n.cells, nrowann))
  sim.annot %<>% mutate(ngenes=rep(n.genes, nrowann), lib.loc=rep(lib.loc, nrowann))
  sim.annot %<>% mutate(cellid = paste0(sim.annot$Cell, '-', n.cells, '.', n.genes, ' ', lib.loc))
  sim.cm <- counts(sim.result)
  sim.cm <- sim.cm[, as.character(sim.annot$Cell)]
  colnames(sim.cm) <- sim.annot$cellid
  sim.genes <- rownames(sim.cm)
  return(list(cm=t(sim.cm), sim.annot=sim.annot, whole.result=sim.result))
}

lapplyCells <- function(n.cell.vec, seed, lib.loc=8, n.genes=10000) {

  results <- n.cell.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=x, n.genes=n.genes, lib.loc=lib.loc)
    return(a.sim)
  })
  boundcm <- do.call(rbind, results %>% lapply(function(x){x$cm}))
  boundannot <- do.call(rbind, results %>% lapply(function(x){x$sim.annot}))
  boundannot %<>% mutate(group=as.numeric(gsub('Group', '', Group)))
  #names(results) <- n.cell.vec %>% as.character()
  return(return(list(cm=boundcm, annot=boundannot)))
}

lapplyGenes <- function(n.gene.vec, seed, n.cells=400, lib.loc=8) {
  results <- n.gene.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=x, lib.loc=lib.loc)
    return(a.sim)
  })
  names(results) <- n.gene.vec %>% as.character()
  return(results)
}

lapplyCover <- function(lib.loc.vec, mean.shape.vec, seed, n.genes=10000, n.cells=400) {
  doSimulateGroups <- function(lib.loc, mean.shape, seed, n.genes, n.cells){
    return(SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=n.genes, lib.loc=lib.loc, mean.shape=mean.shape))
  }
  results <- Map(doSimulateGroups, lib.loc.vec, mean.shape.vec, MoreArgs=list(seed, n.genes, n.cells))
  boundcm <- do.call(rbind, results %>% lapply(function(x){x$cm}))
  boundannot <- do.call(rbind, results %>% lapply(function(x){x$sim.annot}))
  boundannot %<>% mutate(group=as.numeric(gsub('Group', '', Group)))
  return(return(list(cm=boundcm, annot=boundannot)))
}

annotSubtypeCondition <- function(simannot, arrange.factor){
  if (arrange.factor=='ncell') {
    simannot %<>% arrange(ncell, group)
  } else if (arrange.factor=='ngenes') {
    simannot %<>% arrange(ngenes, group)
  } else if (arrange.factor=='cover') {
    simannot %<>% arrange(lib.loc, group)
  }
  #simannot %<>% arrange(ncell, group)
  simannot$Cell <- as.character(simannot$Cell)
  de.test.level <- c(0, 0, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5)
  condition.level <- rep(c('healthy', 'diseased'), 5)
  de.switch <- setNames(de.test.level, simannot$Group %>% unique)
  condition.switch <- setNames(condition.level, simannot$Group %>% unique)
  simannot %<>% mutate(de.level=de.switch[Group], condition=condition.switch[Group])
  if (arrange.factor=='ncell') {
    simannot %<>% mutate(subtype=paste(ncell, de.level))
  } else if (arrange.factor=='ngenes') {
    simannot %<>% mutate(subtype=paste(ngenes, de.level))
  } else if (arrange.factor=='cover') {
    simannot %<>% mutate(subtype=paste(lib.loc, de.level))
  }
  return(simannot)
}

MakeSimRepPaga <- function(rep.annot, rep.cm, varied.factor) {
  rep.cm <- rep.cm[rep.annot$cellid, ]
  rep.annot$Cell <- rep.annot$Cell %>% as.character
  batch.p2 <- NeuronalMaturation::GetPagoda(Matrix::t(rep.cm), n.odgenes=2000)
  batch.knn <- igraph::as_adjacency_matrix(batch.p2$graphs$PCA, attr='weight')[rep.annot$cellid, rep.annot$cellid]
  paga.batch <- GeneratePagaItems(batch.knn, rep.annot$subtype, rep.annot$condition, by.subtypes.condition=T)
  n.cellmatch <- sum(rep.annot$cellid %in% rownames(batch.knn))
  if(n.cellmatch != nrow(rep.annot)) {
    stop('dimensions do not match')
  }
  paga.df.batch <- paga.batch$paga.df
  if (varied.factor=='ncells') {
    paga.df.batch %<>% mutate(ncell=gsub(' .*', '', subtype), de.prob = gsub('.* ', '', subtype))
    paga.df.batch$ncell <- paga.df.batch$ncell %>% as.numeric %>% as.factor
  } else if (varied.factor=='ngenes') {
    paga.df.batch %<>% mutate(ngenes=gsub(' .*', '', subtype), de.prob = gsub('.* ', '', subtype))
    paga.df.batch$ngenes <- paga.df.batch$ngenes %>% as.numeric %>% as.factor
  } else if (varied.factor=='cover') {
    paga.df.batch %<>% mutate(lib.loc=gsub(' .*', '', subtype), de.prob = gsub('.* ', '', subtype))
    paga.df.batch$lib.loc <- paga.df.batch$lib.loc %>% as.numeric %>% as.factor
  }
  return(paga.df.batch)
}
