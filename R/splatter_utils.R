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

SimulateGroups <- function(splatter.params, n.cells, de.probabilities = c(0, 0, 0, 0.2, 0, 0.3, 0, 0.4, 0, 0.5),
                           group.prob = c(rep(1/10, 10)), n.genes, lib.loc, mean.shape=0.414, seed=22071, dropout.type='none') {
  sim.result <- splatSimulateGroups(splatter.params, group.prob = group.prob, dropout.type=dropout.type, mean.shape=mean.shape,
                                    de.prob = de.probabilities, nGenes=n.genes, batchCells=n.cells*length(group.prob),
                                    lib.loc=lib.loc, verbose = FALSE, seed=seed)
  #browser()
  sim.annot <- sim.result@colData %>% as.data.frame
  nrowann <- dim(sim.annot)[1]
  sim.annot %<>% mutate(ncell=rep(n.cells, nrowann), mean.shape=as.numeric(rep(mean.shape, nrowann)))
  sim.annot %<>% mutate(ngenes=rep(n.genes, nrowann), lib.loc=rep(lib.loc, nrowann))
  sim.annot %<>% mutate(cellid = paste0(sim.annot$Cell, '-', n.cells, '.', n.genes, ' ', lib.loc))
  sim.cm <- counts(sim.result)
  sim.cm <- sim.cm[, as.character(sim.annot$Cell)]
  colnames(sim.cm) <- sim.annot$cellid
  sim.genes <- rownames(sim.cm)
  return(list(cm=t(sim.cm), sim.annot=sim.annot, whole.result=sim.result))
}

lapplyCells <- function(n.cell.vec, seed, lib.loc=8, n.genes=10000,  de.probs, bind=F) {
  group.probs <- rep(1/length(de.probs), length(de.probs))
  results <- n.cell.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=x, n.genes=n.genes, lib.loc=lib.loc, de.probabilities=de.probs,
                            group.prob=group.probs)
    return(a.sim)
  })
  if (bind) {
    boundcm <- do.call(rbind, results %>% lapply(function(x){x$cm}))
    boundannot <- do.call(rbind, results %>% lapply(function(x){x$sim.annot}))
    boundannot %<>% mutate(group=as.numeric(gsub('Group', '', Group)))
    return.value <- list(cm=boundcm, annot=boundannot)
  } else {
    cms <- results %>% lapply(function(x){x$cm})
    annots <- results %>% lapply(function(x){x$sim.annot})
    annots <- annots %>% lapply(function(x){x %<>% mutate(group=as.numeric(gsub('Group', '', Group)))})
    names(cms) <- n.cell.vec -> names(annots)
    return.value <- list(cms=cms, annots=annots)
  }
  return(return.value)
}

lapplyLibloc <- function(lib.loc.vec, seed, n.cells=500, n.genes=10000,  de.probs, bind=F) {
  group.probs <- rep(1/length(de.probs), length(de.probs))
  results <- lib.loc.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=n.genes, lib.loc=x, de.probabilities=de.probs,
                            group.prob=group.probs)
    return(a.sim)
  })
  if (bind) {
    boundcm <- do.call(rbind, results %>% lapply(function(x){x$cm}))
    boundannot <- do.call(rbind, results %>% lapply(function(x){x$sim.annot}))
    boundannot %<>% mutate(group=as.numeric(gsub('Group', '', Group)))
    return.value <- list(cm=boundcm, annot=boundannot)
  } else {
    cms <- results %>% lapply(function(x){x$cm})
    annots <- results %>% lapply(function(x){x$sim.annot})
    annots <- annots %>% lapply(function(x){x %<>% mutate(group=as.numeric(gsub('Group', '', Group)))})
    names(cms) <- lib.loc.vec -> names(annots)
    return.value <- list(cms=cms, annots=annots)
  }
  return(return.value)
}

lapplyGenes <- function(n.gene.vec, seed, n.cells=500, lib.loc=8,  de.probs, bind=F) {
  #browser()
  group.probs <- rep(1/length(de.probs), length(de.probs))
  results <- n.gene.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=x, lib.loc=lib.loc,
                            group.prob=group.probs, de.probabilities=de.probs)
    return(a.sim)
  })
  if (bind) {
    boundcm <- do.call(rbind, results %>% lapply(function(x){x$cm}))
    boundannot <- do.call(rbind, results %>% lapply(function(x){x$sim.annot}))
    boundannot %<>% mutate(group=as.numeric(gsub('Group', '', Group)))
    return.value <- list(cm=boundcm, annot=boundannot)
  } else {
    cms <- results %>% lapply(function(x){x$cm})
    annots <- results %>% lapply(function(x){x$sim.annot})
    annots <- annots %>% lapply(function(x){x %<>% mutate(group=as.numeric(gsub('Group', '', Group)))})
    names(cms) <- n.gene.vec -> names(annots)
    return.value <- list(cms=cms, annots=annots)
  }
  return(return.value)
}

lapplyGenesOld <- function(n.gene.vec, seed, n.cells=500, lib.loc=8,  de.probs, bind=F) {
  results <- n.gene.vec %>% lapply(function(x){
    a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=x, lib.loc=lib.loc,
                            group.prob=group.probs)
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
  simannot$ngenes <- as.numeric(simannot$ngenes)
  simannot$group <- as.numeric(simannot$group)
  simannot$lib.loc <- as.numeric(simannot$lib.loc)
  simannot$ncell <- as.numeric(simannot$ncell)
  if (arrange.factor=='ncell') {
    simannot %<>% arrange(ncell, group)
  } else if (arrange.factor=='ngenes') {
    simannot %<>% arrange(ngenes, group)
  } else if (arrange.factor=='cover') {
    simannot %<>% arrange(lib.loc, group)
  }
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
    simannot %<>% mutate(subtype=paste0(lib.loc, '_', mean.shape,  ' ', de.level))
  }
  return(simannot)
}

MakeSimRepPaga <- function(rep.annot, rep.cm, varied.factor, n.odgenes=3000, embedding.type=NULL) {
  #browser()
  if (class(rep.cm)!='Pagoda2') {
    rep.cm <- rep.cm[rep.annot$cellid, ]
    rep.annot$Cell <- rep.annot$Cell %>% as.character
    batch.p2 <- NeuronalMaturation::GetPagoda(Matrix::t(rep.cm), n.odgenes=n.odgenes, embeding.type=embedding.type)
    batch.knn <- igraph::as_adjacency_matrix(batch.p2$graphs$PCA, attr='weight')[rep.annot$cellid, rep.annot$cellid]
  } else {
    batch.knn <- igraph::as_adjacency_matrix(rep.cm$graphs$PCA, attr='weight')[rep.annot$cellid, rep.annot$cellid]
  }
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
    paga.df.batch %<>% mutate(cover=gsub(' .*', '', subtype), de.prob = gsub('.* ', '', subtype))
    paga.df.batch %<>% mutate(cover=gsub(' .*', '', subtype), de.prob = gsub('.* ', '', subtype))
    paga.df.batch %<>% mutate(lib.loc=gsub('_.*', '', cover), mean.shape=gsub('.*_', '', cover))
    paga.df.batch$lib.loc <- paga.df.batch$lib.loc %>% as.numeric %>% as.factor
  }
  return(paga.df.batch)
}

ProcessGeneSim <- function(genesim, gene.name.vec){
  #browser()
  annot <- do.call(rbind, genesim %>% lapply(function(x){x$sim.annot}))
  cms <- genesim %>% lapply(function(x){x$cm})
  cms <- cms %>% lapply(Matrix::t) %>% lapply(PadGenesRows, gene.name.vec) %>% lapply(Matrix::t)
  annot %<>% mutate(group=gsub('Group', '', Group))
  cmbound <- do.call(rbind, cms)
  return(list(cm=cmbound, annot=annot))
}

performDistance <- function(subtype.vector, cellid.vector, condition.vector, count.matrix, genes){
  sub.mats <- GetSubMatrices(subtype.vector, cellid.vector, condition.vector, count.matrix, genes, avg=T)
  corr.dists <- Map(function(x,y){1-cor(x,y)}, sub.mats[[1]], sub.mats[[2]])
  names(corr.dists) <- names(sub.mats[[1]])
  ncells <- table(subtype.vector)[names(sub.mats[[1]])]
  return(dplyr::bind_cols(correlation.distance=unlist(corr.dists), subtype=names(sub.mats[[1]]), ncells=ncells))
}

ExtendFactorAnnot <- function(factor.annot, factor.name, factor.class){
  de.switch <- setNames(c('ref', 0, 0.2, 0.3, 0.4, 0.5), (factor.annot$Group %>% unique %>% sort))
  factor.annot %<>% dplyr::mutate(de.level=de.switch[factor.annot$Group],
                                  factor.name=as.numeric(rep(factor.name, nrow(factor.annot))))
  factor.annot %<>% dplyr::rename(!!factor.class:=factor.name)
  return(factor.annot)
}

SimPaga3 <- function(factor.cm, factor.annot.extended, factor.iteration, factor.identity, static.factors=NULL, return.p2=F){
  #browser()
  factor.annot <- factor.annot.extended
  if (class(factor.cm)!='Pagoda2') {
    p2 <- NeuronalMaturation::GetPagoda(Matrix::t(factor.cm), n.odgenes=3000, embeding.type = NULL)
    graph <- igraph::as_adjacency_matrix(p2$graphs$PCA, attr='weight')[factor.annot$cellid, factor.annot$cellid]
  } else {
    graph <- igraph::as_adjacency_matrix(factor.cm$graphs$PCA, attr='weight')[factor.annot$cellid, factor.annot$cellid]
  }
  sample.connectivity <- GeneratePagaItems(graph, sample.vector=factor.annot.extended$de.level, by.sample=T)
  connectivity.mat <- sample.connectivity$connectivities
  de.levels <- (paste0(factor.annot$de.level) %>% unique)[order(paste0(factor.annot$de.level) %>% unique)]
  n.levels <- length(de.levels)
  comparisons <- connectivity.mat[n.levels, ] %>% setNames(de.levels)
  ncell.delevels <- factor.annot$cellid %>% split(factor.annot$de.level) %>% lapply(length) %>% unlist
  ncell.true <- ncell.delevels+ncell.delevels[n.levels]
  df <- dplyr::bind_cols(paga.connectivity.value=comparisons, de.levels=de.levels, ncell.comparison=ncell.true,
                         which.factor=rep(factor.iteration, n.levels))
  df %<>% dplyr::rename(!!factor.identity:=which.factor)
  if (return.p2) {
    return(list(p2=p2, paga.df=df))
  } else {
    return(df)
  }
}

SimPaga4 <- function(connectivity.mat, factor.annot.extended, factor.iteration, factor.identity){
  factor.annot <- factor.annot.extended
  de.levels <- (paste0(factor.annot$de.level) %>% unique)[order(paste0(factor.annot$de.level) %>% unique)]
  n.levels <- length(de.levels)
  comparisons <- connectivity.mat[n.levels, ] %>% setNames(de.levels)
  ncell.delevels <- factor.annot$cellid %>% split(factor.annot$de.level) %>% lapply(length) %>% unlist
  ncell.true <- ncell.delevels+ncell.delevels[n.levels]
  df <- dplyr::bind_cols(paga.connectivity.value=comparisons, de.levels=de.levels, ncell.comparison=ncell.true,
                         which.factor=rep(factor.iteration, n.levels))
  df %<>% dplyr::rename(!!factor.identity:=which.factor)
  return(df)
}

SimPagaFactor <- function(factor.list, factor.class, return.p2=F){
  #browser()
  if (names(factor.list)[1]=='p2s') {
    names(factor.list)[1]='cms'
  }
  factor.names <- names(factor.list$cms)
  factor.class.vector <- rep(factor.class, length(factor.names))
  extended.annots <- Map(ExtendFactorAnnot, factor.list$annots, factor.names, factor.class.vector)
  paga.results <- Map(SimPaga3, factor.list$cms, extended.annots, factor.names, factor.class.vector,
                      MoreArgs=list(return.p2=return.p2))
  if (return.p2==F) {
    paga.df <- dplyr::bind_rows(paga.results)
    paga.df[[factor.class]] <- paga.df[[factor.class]] %>% as.numeric %>% as.factor
    return(paga.df)
  } else {
    p2s <- paga.results %>% lapply(function(x){x$p2})
    paga.dfs <- paga.results %>% lapply(function(x){x$paga.df})
    paga.df <- dplyr::bind_rows(paga.dfs)
    paga.df[[factor.class]] <- paga.df[[factor.class]] %>% as.numeric %>% as.factor
    return(list(paga.df=paga.df, p2s=p2s))
  }
}

SimCor <- function(cm, annot, factor.iteration, factor.class){
  subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]}) %>% lapply(Matrix::colMeans)
  n.de.levels <- length(subvecs)
  ref.mat <- subvecs[[n.de.levels]]
  cor.dists <- subvecs %>% lapply(function(x){1-cor(x,ref.mat)}) %>% unlist
  df <- dplyr::bind_cols(correlation.distance=cor.dists, de.levels=names(subvecs))
  df %<>% dplyr::mutate(which.factor = rep(factor.iteration, n.de.levels))
  df %<>% dplyr::rename(!!factor.class:=which.factor)
  return(df)
}

SimDist <- function(cm, annot, factor.iteration, factor.class, distance, pseudo.prob=1e-6){
  subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]}) %>% lapply(Matrix::colMeans)
  n.de.levels <- length(subvecs)
  ref.mat <- subvecs[[n.de.levels]]
  if(distance=='correlation.distance') {
    dist.function <- function(x){1-cor(x, ref.mat)}
  } else if (distance=='jensen_shannon'){
    #browser()
    renormalize <- function(a.vec){
      a.vec <- a.vec/sum(a.vec); a.vec <- a.vec+pseudo.prob; a.vec <- a.vec/sum(a.vec); return(a.vec)}
    subvecs <- subvecs %>% lapply(renormalize)
    ref.mat <- subvecs[[n.de.levels]]
    dist.function <- function(x){JensenShannon(x, ref.mat)}
  }
  cor.dists <- subvecs %>% lapply(dist.function) %>% unlist
  df <- dplyr::bind_cols(correlation.distance=cor.dists, de.levels=names(subvecs))
  df %<>% dplyr::mutate(which.factor = rep(factor.iteration, n.de.levels))
  df %<>% dplyr::rename(!!factor.class:=which.factor)
  df %<>% dplyr::rename(!!distance:=correlation.distance)
  return(df)
}

doSimCor <- function(cms, annots, factor.class){
  extended.annots <- Map(ExtendFactorAnnot, annots, names(annots), factor.class)
  factor.iterations <- names(annots)
  dfs <- Map(SimCor, cms, extended.annots, factor.iterations, MoreArgs=list(factor.class=factor.class))
  df <- dplyr::bind_rows(dfs)
  df[[factor.class]] <- df[[factor.class]] %>% as.numeric %>% as.factor
  return(df)
}

doSimDist <- function(cms, annots, factor.class, distance){
  extended.annots <- Map(ExtendFactorAnnot, annots, names(annots), factor.class)
  factor.iterations <- names(annots)
  dfs <- Map(SimDist, cms, extended.annots, factor.iterations, MoreArgs=list(factor.class=factor.class, distance=distance))
  df <- dplyr::bind_rows(dfs)
  df[[factor.class]] <- df[[factor.class]] %>% as.numeric %>% as.factor
  return(df)
}

PairwiseComparisonsFullMat <- function(p2.list.or.pca.cm, annot.list, factor.identity, return.p2=T){
  some.object <- p2.list.or.pca.cm
  if (class(some.object)=='list'){
    if (some.object %>% lapply(class) %>% unique=='Pagoda2'){
      if (factor.identity=='ngenes'){
        mat.n.cols <- some.object %>% lapply(function(x){x$misc$rawCounts %>% ncol})
        max.ind <- which(mat.n.cols %>% unlist==max(mat.n.cols %>% unlist))
        gene.vec <- some.object[[max.ind]]$misc$rawCounts %>% colnames
        cm.raw <- some.object %>% lapply(function(x){x$misc$rawCounts}) %>% lapply(Matrix::t) %>%
          lapply(PadGenesRows, gene.vec) %>% lapply(Matrix::t) %>% do.call(rbind,.)
      } else {
        cm.raw <- some.object %>% lapply(function(x){x$misc$rawCounts}) %>% do.call(rbind,.)
      }
      p2 <- NeuronalMaturation::GetPagoda(Matrix::t(cm.raw), n.odgenes=3000, embeding.type=NULL)
      cm.pca <- p2$reductions$PCA
    } else {
      if (factor.identity=='ngenes'){
        mat.n.cols <- some.object %>% lapply(function(x){x %>% ncol})
        max.ind <- which(mat.n.cols %>% unlist==max(mat.n.cols %>% unlist))
        gene.vec <- some.object[[max.ind]]%>% colnames
        cm.raw <- some.object %>% lapply(function(x){x}) %>% lapply(Matrix::t) %>%
          lapply(PadGenesRows, gene.vec) %>% lapply(Matrix::t) %>% do.call(rbind,.)
      } else {
        cm.raw <- some.object %>% do.call(rbind,.)
      }
      p2 <- NeuronalMaturation::GetPagoda(Matrix::t(cm.raw), n.odgenes=3000, embeding.type=NULL)
      cm.pca <- p2$reductions$PCA
    }
  } else {
    if (class(some.object=='Pagoda2')){
      cm.pca <- some.object$reductions$PCA
    } else {
      cm.pca <- some.object
    }
  }
  factor.levels <- (annot.list %>% lapply(function(x){x[[factor.identity]] %>% unique})) %>% unlist
  annots.extended <- Map(ExtendFactorAnnot, annot.list, factor.levels, MoreArgs=(list(factor.identity)))
  bound.annot <- annots.extended %>% dplyr::bind_rows()
  partition.levels <- paste0(bound.annot[[factor.identity]], '_', bound.annot$de.level) %>% str_sort(numeric=TRUE) %>% unique
  bound.annot$partitions <- factor(paste0(bound.annot$ncell, '_', bound.annot$de.level), levels=partition.levels)
  if (return.p2) {
    return(list(p2=p2, bound.annot=bound.annot))
  } else {
    return(list(cm.pca=cm.pca, bound.annot=bound.annot))
  }
}

PagaForBound <- function(p2, annot, varied.factor){
  graph <- igraph::as_adjacency_matrix(p2$graphs$PCA, attr='weight')[annot$cellid, annot$cellid]
  mem.vec <- annot$partitions %>% as.numeric
  part.connects <- GetPagaMatrix(graph, mem.vec)
  inds <- 1:30 %>% split(seq(1,30, by=5)) %>% as.data.frame %>% t %>% as.matrix %>% as.data.frame
  sub.mats <- inds %>% lapply(function(x){part.connects[x, x]})
  factor.levels <- annot[[varied.factor]] %>% unique %>% sort
  names(sub.mats) <- factor.levels
  annot.split <- annot %>% split(annot[[varied.factor]])
  cons <- Map(SimPaga4, sub.mats, annot.split, factor.levels, MoreArgs=list(varied.factor))
  names(cons) <- factor.levels
  return(cons)
}

GetProbDistPerSeed <- function(p2s.anns, factor.class, distance='correlation.distance'){
  #browser()
  if (distance=='correlation.distance'){
    cms <- p2s.anns$p2s %>% lapply(function(x){x$reductions$PCA})
  } else if (distance=='jensen_shannon'){
    cms <- p2s.anns$p2s %>% lapply(function(x){x$counts})
  }
  dist.df <- doSimDist(cms, p2s.anns$annots, factor.class, distance)
}
