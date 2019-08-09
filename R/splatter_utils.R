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

SimulateGroups <- function(splatter.params, n.cells, de.probabilities,
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

MakeSimPerFactor <- function(factor.vec, n.cells=500, seed, lib.loc=8, n.genes=10000,  de.probs, factor.identity,
                             n.cl.sim=NULL) {
  group.probs <- rep(1/length(de.probs), length(de.probs))
  if (factor.identity=='ncell'){
    results <- factor.vec %>% pbapply::pblapply(function(x){
      a.sim <- SimulateGroups(params, seed=seed, n.cells=x, n.genes=n.genes, lib.loc=lib.loc, de.probabilities=de.probs,
                            group.prob=group.probs); return(a.sim)}, cl=n.cl.sim)
  } else if (factor.identity=='ngenes') {
    results <- factor.vec %>% pbapply::pblapply(function(x){
      a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=x, lib.loc=lib.loc, de.probabilities=de.probs,
                              group.prob=group.probs); return(a.sim)}, cl=n.cl.sim)
  } else if (factor.identity=='lib.loc') {
    results <-factor.vec %>% pbapply::pblapply(function(x){
      a.sim <- SimulateGroups(params, seed=seed, n.cells=n.cells, n.genes=n.genes, lib.loc=x, de.probabilities=de.probs,
                              group.prob=group.probs); return(a.sim)}, cl=n.cl.sim)
  }

  cms <- results %>% lapply(function(x){x$cm})
  annots <- results %>% lapply(function(x){x$sim.annot})
  annots <- annots %>% lapply(function(x){x %<>% mutate(group=as.numeric(gsub('Group', '', Group)))})
  annots <- annots %>% lapply(appendDelevel, de.probs=de.probs)
  names(cms) <- factor.vec -> names(annots)
  return.value <- list(cms=cms, annots=annots)
  return(return.value)
}

appendDelevel <- function(annot, de.probs){
  if (is.null(names(de.probs))){
    de.level.vec <- c('ref', de.probs[2:length(de.probs)])
  } else {
    de.level.vec <- names(de.probs)
  }
  #browser()
  de.switch <- setNames(de.level.vec, annot$Group %>% unique %>% sort)
  de.levels <- de.switch[annot$Group]
  annot %<>% mutate(de.level=de.levels)
  return(annot)
}

MakeSimsPerSeed <- function(seed, factor.vec, de.probs, factor.identity, n.cl.sim){
  sims <- MakeSimPerFactor(factor.vec, de.probs=de.probs, seed=seed, factor.identity=factor.identity, n.cl.sim=n.cl.sim)
  return(sims)
}

MakeSimsAllSeeds <- function(seed.vec, factor.vec, de.probs, factor.identity, make.p2=T, n.cl.tsne=30, n.cl.sim=NULL){
  sims.per.seed <- seed.vec %>% lapply(MakeSimsPerSeed, factor.vec, de.probs, factor.identity, n.cl.sim=n.cl.sim)
  names(sims.per.seed) <- seed.vec
  if (make.p2) {
    sims.per.seed <- sims.per.seed %>% lapply(MakeP2PerSeed, n.cl.tsne=n.cl.tsne)
  }
  return(sims.per.seed)
}

MakeP2PerSeed <- function(cms.anns.list, n.cl.tsne){
  p2s <- cms.anns.list$cms %>% lapply(Matrix::t) %>% lapply(NeuronalMaturation::GetPagoda, n.odgenes=3000, n.cores=n.cl.tsne)
  return(list(p2s=p2s, annots=cms.anns.list$annots))
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

ExtendFactorAnnot <- function(factor.annot, factor.name=NULL, factor.class=NULL, de.level.vec=c('ref', 0, 0.2, 0.3, 0.4, 0.5)){
  de.switch <- setNames(de.level.vec, (factor.annot$Group %>% unique %>% sort))
  factor.annot %<>% dplyr::mutate(de.level=de.switch[factor.annot$Group])
  #factor.annot %<>% dplyr::mutate(de.level=de.switch[factor.annot$Group],
                                  #factor.name=as.numeric(rep(factor.name, nrow(factor.annot))))
  #factor.annot %<>% dplyr::rename(!!factor.class:=factor.name)
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
  df <- dplyr::bind_cols(paga.connectivity.value=comparisons, de.level=de.levels, ncell.comparison=ncell.true,
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
  df <- dplyr::bind_cols(paga.connectivity.value=comparisons, de.level=de.levels, ncell.comparison=ncell.true,
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
  if (is.null(factor.list$annots[[1]]$de.level)) {
    extended.annots <- Map(ExtendFactorAnnot, factor.list$annots, factor.names, factor.class.vector)
  } else {
    extended.annots <- factor.list$annots
  }
  paga.results <- Map(SimPaga3, factor.list$cms, extended.annots, factor.names, factor.class.vector,
                      MoreArgs=list(return.p2=return.p2))
  if (return.p2==F) {
    paga.df <- dplyr::bind_rows(paga.results)
    paga.df[[factor.class]] <- paga.df[[factor.class]] %>% as.numeric #%>% as.factor
    paga.df <- paga.df[, c(1,2,4,3)]
    return(paga.df)
  } else {
    p2s <- paga.results %>% lapply(function(x){x$p2})
    paga.dfs <- paga.results %>% lapply(function(x){x$paga.df})
    paga.df <- dplyr::bind_rows(paga.dfs)
    paga.df[[factor.class]] <- paga.df[[factor.class]] %>% as.numeric #%>% as.factor
    paga.df <- paga.df[, c(1,2,4,3)]
    return(list(paga.df=paga.df, p2s=p2s))
  }
}

SimCor <- function(cm, annot, factor.iteration, factor.class){
  subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]}) %>% lapply(Matrix::colMeans)
  n.de.levels <- length(subvecs)
  ref.mat <- subvecs[[n.de.levels]]
  cor.dists <- subvecs %>% lapply(function(x){1-cor(x,ref.mat)}) %>% unlist
  df <- dplyr::bind_cols(correlation.distance=cor.dists, de.level=names(subvecs))
  df %<>% dplyr::mutate(which.factor = rep(factor.iteration, n.de.levels))
  df %<>% dplyr::rename(!!factor.class:=which.factor)
  return(df)
}

SimDist <- function(cm, annot, factor.iteration, factor.class, distance, pseudo.prob=1e-6, col.medians=F){
  if (col.medians){
    subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]}) %>% lapply(function(x){x %>%
        as.matrix %>% matrixStats::colMedians})
  } else {
    subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]}) %>% lapply(Matrix::colMeans)
  }
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
  } else if (distance=='euclidean'){
    dist.function <- function(x){sqrt(sum((x-ref.mat)^2))}
  } else if (distance=='CMD') {
    subvecs <- annot$cellid %>% split(annot$de.level) %>% lapply(function(x){cm[x, ]})
    ref.mat <- subvecs[[n.de.levels]]
    dist.function <- function(x){
      cor.ref <- cor(ref.mat)
      cor.dis <- cor(x)
      return(CMD(cor.ref, cor.dis))
    }
  }
  cor.dists <- subvecs %>% lapply(dist.function) %>% unlist
  df <- dplyr::bind_cols(correlation.distance=cor.dists, de.level=names(subvecs))
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
  if (is.null(annots[[1]]$de.level)) {
    extended.annots <- Map(ExtendFactorAnnot, annots, names(annots), factor.class)
  } else {
    extended.annots <- annots
  }
  factor.iterations <- names(annots)
  dfs <- Map(SimDist, cms, extended.annots, factor.iterations, MoreArgs=list(factor.class=factor.class, distance=distance))
  df <- dplyr::bind_rows(dfs)
  df[[factor.class]] <- df[[factor.class]] %>% as.numeric #%>% as.factor
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
  if (is.null(annot.list[[1]]$de.level)) {
    annots.extended <- Map(ExtendFactorAnnot, annot.list, factor.levels, MoreArgs=(list(factor.identity)))
  } else {
    annots.extended <- annot.list
  }
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
  n.inds <- part.connects %>% nrow()
  #browser()
  n.de.levels <- length(annot$de.level %>% unique)
  n.split.by <- n.inds/n.de.levels
  inds <- 1:n.inds %>% split(seq(1,n.inds, by=n.split.by)) %>% as.data.frame %>% t %>% as.matrix %>% as.data.frame
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
  } else if (distance=='euclidean' | distance=='CMD'){
    cms <- p2s.anns$p2s %>% lapply(function(x){x$reductions$PCA})
  }
  dist.df <- doSimDist(cms, p2s.anns$annots, factor.class, distance)
}

KnnCor <- function(p2, annot, avg.meds=F){
  #browser()
  if(is.null(annot$de.level)){
    annot <- annot %>% ExtendFactorAnnot
  }
  de.split <- annot %>% split(annot$de.level)
  ref.annot <- de.split$ref
  ref.cellid <- ref.annot$cellid
  non.ref <- de.split[names(de.split)[names(de.split)!='ref']]
  #ref.annot %<>% mutate(condition='healthy')
  #de.split[non.ref] <- de.split[non.ref] %>% lapply(function(x){x %<>% mutate(condition='diseased')})
  getConcAnn <- function(x){
    diseased.cellids <- x$cellid
    per.type <- list(diseased=diseased.cellids, healthy=ref.cellid)
    return(per.type)
  }
  cell.ids.per.clust.per.type <- non.ref %>% lapply(getConcAnn)
  cor.dists.pairwise <- lapply(cell.ids.per.clust.per.type, scConditionDifference::corDistPairwise, p2$reductions$PCA,
                               n.closest=15)

  cor.dists.disease <- lapply(cell.ids.per.clust.per.type, scConditionDifference::corDistPairwise, p2$reductions$PCA,
                              n.closest=15, group.ids=c(1, 1))

  cor.dists.healthy <- lapply(cell.ids.per.clust.per.type, scConditionDifference::corDistPairwise, p2$reductions$PCA,
                              n.closest=15, group.ids=c(2, 2))
  cor.dists.pairwise.avg <- lapply(cor.dists.pairwise, rowMeans)
  cor.dists.disease.avg <- lapply(cor.dists.disease, rowMeans)
  #cor.dists.disease.avg <- cor.dists.disease.avg %>% lapply(function(x){x[!duplicated(names(x))]})
  cor.dists.healthy.avg <- lapply(cor.dists.healthy, rowMeans)
  #cor.dists.healthy.avg <- cor.dists.healthy.avg %>% lapply(function(x){x[!duplicated(names(x))]})
  plot.df <- scConditionDifference::distancesToDataFrame(cor.dists.pairwise.avg, trim.level=NULL, type="between") %>%
    rbind(scConditionDifference::distancesToDataFrame(cor.dists.healthy.avg, trim.level=NULL, type="healthy")) %>%
    rbind(scConditionDifference::distancesToDataFrame(cor.dists.disease.avg, trim.level=NULL, type="diseased"))
  triple.plot <- ggplot(plot.df) + geom_violin(aes(x=Cluster, y=Distance, fill=Type)) +
    theme(axis.text.x = element_text(angle = 30, vjust=1, hjust=1), legend.justification=c(1, 1),
          legend.background=element_blank())
  if(avg.meds) {
    mean.dist <- mapply(function(v1, v2) {med1 <- median(v1); med2 <- median(v2); return(mean(c(med1,med2)))},
                        cor.dists.healthy.avg[names(cor.dists.disease.avg)], cor.dists.disease.avg)
    mean.std <- mapply(function(v1, v2) {mad1 <- mad(v1); mad2 <- mad(v2); return(mean(c(mad1,mad2)))},
                       cor.dists.healthy.avg[names(cor.dists.disease.avg)], cor.dists.disease.avg)
  } else {
    mean.dist <- mapply(function(v1, v2) median(c(v1, v2)),
                        cor.dists.healthy.avg[names(cor.dists.disease.avg)], cor.dists.disease.avg)
    mean.std <- mapply(function(v1, v2) mad(c(v1, v2)),
                       cor.dists.healthy.avg[names(cor.dists.disease.avg)], cor.dists.disease.avg)
  }
  cor.dists.pairwise.avg.centered <- mapply(`-`, cor.dists.pairwise.avg[names(mean.dist)], mean.dist)
  cor.dists.pairwise.avg.norm <- mapply(`/`, cor.dists.pairwise.avg.centered[names(mean.std)], mean.std)
  z.df <- distancesToDataFrame(cor.dists.pairwise.avg.norm, trim=NULL)
  z.df %<>% dplyr::rename(de.level=Cluster, knncor.z.score=Distance)
  z.plot <- ggplot(z.df) +
    geom_violin(aes(x=Cluster, y=knncor.z.score), fill=scales::hue_pal()(10)[1]) +
    theme(axis.text.x = element_text(angle = 30, vjust=1, hjust=1), legend.position="none")
  return(list(plot.df=plot.df, z.df=z.df, z.plot=z.plot, triple.plot=triple.plot))
}

doKnnCor <- function(list.p2.anns, factor.identity, avg.meds=T, get.medians=F){
  p2s <- list.p2.anns$p2s
  anns <- list.p2.anns$annots
  knn.cor.results <- Map(KnnCor, p2s, anns, MoreArgs=list(avg.meds=avg.meds))
  z.dfs <- knn.cor.results %>% lapply(function(x){x$z.df})
  appendFactor <- function(df, factor.iteration){
    df$some.variable <- rep(factor.iteration, nrow(df))
    df %<>% dplyr::rename(!!factor.identity:=some.variable)
    return(df)
  }
  factor.levels <- names(anns)
  z.dfs <- Map(appendFactor, z.dfs, factor.levels)
  if (get.medians){
    z.dfs <- z.dfs %>% lapply(medianofz, factor.identity)
    #z.dfs <- z.dfs %>% lapply(function(x){x[[1]]<-x[[1]] %>% as.numeric; return(x)})
    #return(z.dfs)
    z.df.bound <- dplyr::bind_rows(z.dfs)
    z.df.bound[[factor.identity]] <- z.df.bound[[factor.identity]] %>% as.numeric
    z.df.bound <- z.df.bound[, c(3,2,1,4)]
  } else {
    z.df.bound <- dplyr::bind_rows(z.dfs)
    z.df.bound[[factor.identity]] <- z.df.bound[[factor.identity]] %>% as.numeric %>% as.factor
    z.df.bound <- z.df.bound[, c(1,2,4,3)]
  }
  return(z.df.bound)
}

NormalizedRelativeEntropy_bk <- function(p2, annot, leiden.resolution, factor.identity, cl.subgraph=NULL){
  factor.iteration <- annot[[factor.identity]] %>% unique
  annot <- ExtendFactorAnnot(annot, factor.iteration, factor.identity)
  annot.splits <- annot %>% split(annot$de.level)
  conc.annots <- annot.splits[names(annot.splits)!='ref'] %>% lapply(function(x){dplyr::bind_rows(x, annot.splits$ref)})
  conc.cellids <- conc.annots %>% lapply(function(x){x$cellid})
  raw.cm <- p2$misc$rawCounts
  sub.cms <- conc.cellids %>% lapply(function(x){raw.cm[x, ]})
  #browser()
  getsubGraph <- function(sub.cm){
    p2.sub <- sub.cm %>% GetPagoda(n.odgenes=3000, embeding.type = NULL, verbose=F, clustering.type = NULL)
    graph <- p2.sub$graphs$PCA
    return(graph)
  }
  sub.graphs <- sub.cms %>% lapply(Matrix::t) %>% pbapply::pblapply(getsubGraph, cl=cl.subgraph)
  getLeid <- function(graph){
    mem <- conos:::leiden.community(graph, resolution=leiden.resolution)$membership
  }
  mem.per.level <- sub.graphs %>% lapply(getLeid)
  appendLeid <- function(annot, leid.mem.vec){
    mem <- leid.mem.vec[annot$cellid]
    annot$membership <- mem
    return(annot)
  }
  conc.annots2 <- conc.annots %>% Map(appendLeid, ., mem.per.level)
  calcEntropy <- function(conc.annot.2){
    base.count <- conc.annot.2$de.level %>% table
    n.samples <- length(base.count)
    #other.factor <- names(base.count)[names(base.count)!='ref']
    split.by.mem <- conc.annot.2 %>% split(conc.annot.2$membership)
    cluster.counts <- split.by.mem %>% lapply(function(x){x$de.level %>% table})
    ent.numerators <- cluster.counts %>% lapply(doEntKLD, base.count=base.count)
    norm.rel.ent <- sum(ent.numerators %>% unlist)/(log(n.samples)*sum(cluster.counts %>% unlist))
    return(1 - norm.rel.ent)
  }
  doEntKLD <- function(cluster.count, base.count){
    #browser()
    if (length(cluster.count) < 2) {
      orig.name <- names(cluster.count)
      if (orig.name=='ref'){
        cluster.count <- c(1e-3, cluster.count%>% as.numeric)
      } else {
        cluster.count <- c(cluster.count %>% as.numeric, 1e-3)
      }
    }
    c.c.norm <- cluster.count/sum(cluster.count)
    b.c.norm <- base.count/sum(base.count)
    ent.val <- sum(cluster.count) * KLD(c.c.norm, b.c.norm)
    return(ent.val)
  }
  rel.ents <- conc.annots2 %>% lapply(calcEntropy)
  df <- dplyr::bind_cols(norm.rel.entropy=rel.ents %>% unlist, de.level=names(conc.annots2))
  nrowdf <- nrow(df)
  df %<>% mutate(leiden.resolution=rep(leiden.resolution, nrowdf), some.factor=rep(factor.iteration, nrowdf))
  df %<>% dplyr::rename(!!factor.identity:=some.factor)
  return(df)
}

NormalizedRelativeEntropy <- function(p2, annot, leiden.resolutions, factor.identity, cl.subgraph=NULL){
  factor.iteration <- annot[[factor.identity]] %>% unique
  annot <- ExtendFactorAnnot(annot, factor.iteration, factor.identity)
  annot.splits <- annot %>% split(annot$de.level)
  conc.annots <- annot.splits[names(annot.splits)!='ref'] %>% lapply(function(x){dplyr::bind_rows(x, annot.splits$ref)})
  conc.cellids <- conc.annots %>% lapply(function(x){x$cellid})
  raw.cm <- p2$misc$rawCounts
  sub.cms <- conc.cellids %>% lapply(function(x){raw.cm[x, ]})
  #browser()
  getsubGraph <- function(sub.cm){
    p2.sub <- sub.cm %>% GetPagoda(n.odgenes=3000, embeding.type = NULL, verbose=F, clustering.type = NULL, n.pcs = 100)
    graph <- p2.sub$graphs$PCA
    return(graph)
  }
  sub.graphs <- sub.cms %>% lapply(Matrix::t) %>% pbapply::pblapply(getsubGraph, cl=cl.subgraph)
  #browser()
  getLeid <- function(graph, leiden.resolutions){
    #browser()
    mems <- leiden.resolutions %>% lapply(function(x){
      conos:::leiden.community(graph, resolution=x)$membership
  })}
  mems.per.level <- sub.graphs %>% lapply(function(x){getLeid(x,leiden.resolutions)})
  appendLeid <- function(annot, leid.mem.vec){
    #browser()
    mems <- leid.mem.vec %>% lapply(function(x){x[annot$cellid]})
    annots.ext <- mems %>% lapply(function(x){z<-annot; z$membership <- x;return(z)})
    return(annots.ext)
  }
  #browser()
  conc.annots2 <- conc.annots %>% Map(appendLeid, ., mems.per.level)
  calcEntropy <- function(conc.annot.2){
    #browser()
    base.count <- conc.annot.2$de.level %>% table
    n.samples <- length(base.count)
    #other.factor <- names(base.count)[names(base.count)!='ref']
    split.by.mem <- conc.annot.2 %>% split(conc.annot.2$membership)
    cluster.counts <- split.by.mem %>% lapply(function(x){x$de.level %>% table})
    ent.numerators <- cluster.counts %>% lapply(doEntKLD, base.count=base.count)
    norm.rel.ent <- sum(ent.numerators %>% unlist)/(log(n.samples)*sum(cluster.counts %>% unlist))
    return(1 - norm.rel.ent)
  }
  doEntKLD <- function(cluster.count, base.count){
    #browser()
    if (length(cluster.count) < 2) {
      orig.name <- names(cluster.count)
      if (orig.name=='ref'){
        cluster.count <- c(1e-3, cluster.count%>% as.numeric)
      } else {
        cluster.count <- c(cluster.count %>% as.numeric, 1e-3)
      }
    }
    c.c.norm <- cluster.count/sum(cluster.count)
    b.c.norm <- base.count/sum(base.count)
    ent.val <- sum(cluster.count) * KLD(c.c.norm, b.c.norm)
    return(ent.val)
  }
  #browser()
  rel.ents.list <- conc.annots2 %>% lapply(function(x){x %>% lapply(calcEntropy)})
  inds.leid <- 1:length(leiden.resolutions)
  switchEntDe <- function(leid.ind){rel.ents.list %>% lapply(function(x){x[[leid.ind]]})}
  rel.ents.switched <- inds.leid %>% lapply(switchEntDe)
  names(rel.ents.switched) <- leiden.resolutions
  makeRes <- function(ent.res, leiden.resolution){
    df <- dplyr::bind_cols(norm.rel.entropy=ent.res %>% unlist, de.level=names(conc.annots2))
    nrowdf <- nrow(df)
    df %<>% mutate(leiden.resolution=rep(leiden.resolution, nrowdf), some.factor=rep(factor.iteration, nrowdf))
    df %<>% dplyr::rename(!!factor.identity:=some.factor)
    return(df)
  }
  dfs <- rel.ents.switched %>% Map(makeRes, ., leiden.resolutions)
  return(dfs)
}


doEntropy <- function(seed, factor.identity, resolutions, cl.subgraph=NULL){
  #browser()
  factor.levels <- names(seed$p2s)
  makeEntDfs <- function(resolutions, p2s, annots, factor.identity=factor.identity){
    #browser()
    ent.dfs <- Map(SecretUtils::NormalizedRelativeEntropy, p2s, annots, MoreArgs=list(leiden.resolutions=resolutions, #remove s if bk
                                                                                      factor.identity=factor.identity,
                                                                                      cl.subgraph=cl.subgraph))
    #browser()
    #bound.df <- ent.dfs %>% dplyr::bind_rows()
    #return(bound.df)
    return(ent.dfs)
  }
  dfs.per.res <- makeEntDfs(resolutions, seed$p2s, seed$annots, factor.identity)
  #browser()
  dfs.bound <- dfs.per.res %>% lapply(function(x){x %>% dplyr::bind_rows()}) %>% dplyr::bind_rows()
  #return(dfs.per.res %>% dplyr::bind_rows())
  return(dfs.bound)
}

doEntropy_bk <- function(seed, factor.identity, resolutions, cl.subgraph=NULL){
  #browser()
  factor.levels <- names(seed$p2s)
  makeEntDfs <- function(resolution, p2s, annots, factor.identity=factor.identity){
    #browser()
    ent.dfs <- Map(SecretUtils::NormalizedRelativeEntropy, p2s, annots, MoreArgs=list(leiden.resolution=resolution,
                                                                                      factor.identity=factor.identity,
                                                                                      cl.subgraph=cl.subgraph))
    bound.df <- ent.dfs %>% dplyr::bind_rows()
    return(bound.df)
  }
  dfs.per.res <- resolutions %>% lapply(makeEntDfs, seed$p2s, seed$annots, factor.identity)
  return(dfs.per.res %>% dplyr::bind_rows())
}

doBoundPaga <- function(p2s.anns.list, factor.identity){
  boundp2s.anns <- PairwiseComparisonsFullMat(p2s.anns.list$p2s, p2s.anns.list$annots, factor.identity)
  bound.paga.res <- PagaForBound(boundp2s.anns$p2, boundp2s.anns$bound.annot, factor.identity) %>% bind_rows
  return(bound.paga.res)
}

doSomeDistance <- function(x,y,distance){
  df <- x %>% lapply(GetProbDistPerSeed, y, distance) %>% bind_rows; return(df)
}

dodoKnnCor <- function(x,y,...){
  df <- x %>% lapply(doKnnCor, y, ...) %>% bind_rows; return(df)
}

dodoEntropy <- function(x,y,...){
  df <- x %>% lapply(doEntropy, y, ...) %>% bind_rows; return(df)
}

doKnnCorMed <- function(x,y,...){
  df <- x %>% lapply(doKnnCor, y, get.medians=T) %>% bind_rows; return(df)
}

dfsPerDistance <- function(p2.ann.per.seed.per.factor, factor.identities, distance, bind.paga=F,...){
  items.per.seed <- p2.ann.per.seed.per.factor
  if (distance=='paga'){
    if (!bind.paga){
      doPaga <- function(x,y){
        df <- x %>% lapply(SimPagaFactor, y) %>% bind_rows; return(df)
      }
      dfs.per.factor <- Map(doPaga, items.per.seed, factor.identities)
    } else {
      getBoundPaga <- function(x,y){
        df <- x %>% lapply(doBoundPaga, y) %>% bind_rows; return(df)
      }
      dfs.per.factor <- Map(getBoundPaga, items.per.seed, factor.identities)
    }
  } else if (distance=='correlation.distance'){
    dfs.per.factor <- Map(doSomeDistance, items.per.seed, factor.identities, MoreArgs=list(distance=distance))
  } else if (distance=='jensen_shannon'){
    dfs.per.factor <- Map(doSomeDistance, items.per.seed, factor.identities, MoreArgs=list(distance=distance))
  } else if (distance=='euclidean'){
    dfs.per.factor <- Map(doSomeDistance, items.per.seed, factor.identities, MoreArgs=list(distance=distance))
  } else if (distance=='CMD'){
    dfs.per.factor <- Map(doSomeDistance, items.per.seed, factor.identities, MoreArgs=list(distance=distance))
  } else if (distance=='knncor.z'){
    dfs.per.factor <- Map(dodoKnnCor, items.per.seed, factor.identities, MoreArgs=list(...))
  } else if (distance=='entropy'){
    dfs.per.factor <- Map(dodoEntropy, items.per.seed, factor.identities, MoreArgs=list(...))
  } else if (distance=='knncor.z.med'){
    dfs.per.factor <- Map(doKnnCorMed, items.per.seed, factor.identities)
  }
  names(dfs.per.factor) <- factor.identities
  return(dfs.per.factor)
}

dfsPerDist <- function(distance, items.per.seed.per.factor, factor.identities, bind.paga=F, avg.meds=T,
                       leiden.resolutions=NULL){
  if (distance=='knncor.z'){
    dfs <- dfsPerDistance(items.per.seed.per.factor, factor.identities, distance=distance, avg.meds=avg.meds)
  } else if (distance=='entropy'){
    dfs <- dfsPerDistance(items.per.seed.per.factor, factor.identities, distance=distance, resolutions=leiden.resolutions)
  } else if (distance=='paga'){
    dfs <- dfsPerDistance(items.per.seed.per.factor, factor.identities, distance=distance, bind.paga=bind.paga)
  } else {
    dfs <- dfsPerDistance(items.per.seed.per.factor, factor.identities, distance=distance)
  }
}

AllDistsDfs <- function(items.per.seed.per.factor, factor.identities, distances, ...){
  all.dfs <- distances %>% lapply(dfsPerDist, items.per.seed.per.factor, factor.identities, ...)
  names(all.dfs) <- distances
  return(all.dfs)
}

PlotsPerFactor_old <- function(dfs.per.factor, alpha, jitter=F, geom.smooth=T, is.entropy=F, is.knncor=F, is.knncor.med=F){
  factor.identities <- dfs.per.factor %>% names
  if (!is.entropy & !is.knncor &!is.knncor.med){
    plots <- dfs.per.factor %>% lapply(function(x){
      vars <- colnames(x)
      if (jitter){
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_jitter(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_point(alpha=alpha) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    })
  } else if (is.entropy) {
    split.by.resolution <- dfs.per.factor %>% lapply(function(x){split(x, x$leiden.resolution)})
    EntropyPlot <- function(x){
      vars <- colnames(x)
      which.res <- x[3] %>% unique
      if (jitter){
        p <- x %>% ggplot(aes_string(x=vars[4], y=vars[1], col=vars[2]), shape=vars[3])+
          geom_jitter(alpha=alpha)+ggtitle(paste('leiden resolution: ', which.res)) + theme(legend.position="top")
      } else {
        p <- x %>% ggplot(aes_string(x=vars[4], y=vars[1], col=vars[2]), shape=vars[3])+
          geom_point(alpha=alpha)+ggtitle(paste('leiden resolution: ', which.res)) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    }
    plots <- split.by.resolution %>% lapply(function(x){
      plots.for.a.res <- x %>% lapply(EntropyPlot)
    })
  } else if (is.knncor) {
    dfs.per.factor <- dfs.per.factor %>% lapply(function(x){x[[3]] <- x[[3]] %>% as.factor; return(x)})
    plots <- dfs.per.factor %>% lapply(function(x){
      vars <- colnames(x)
      if (!jitter){
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_boxplot(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_violin(alpha=alpha) + theme(legend.position="top")+stat_summary(fun.y=median, geom="point",
                                                                               position=position_dodge(width=0.9))
      }
      if (geom.smooth) {
        #p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
        NULL
      }
      return(p)
    })
  } else if (is.knncor.med) {
    plots <- dfs.per.factor %>% lapply(function(x){
      vars <- colnames(x)
      if (jitter){
        p <- x %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_jitter(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% ggplot(aes_string(x=vars[3], y=vars[1], col=vars[2]))+
          geom_point(alpha=alpha) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    })
  }
  return(plots)
}

PlotsPerFactor <- function(dfs.per.factor, alpha, jitter=F, geom.smooth=T, is.entropy=F, is.knncor=F, is.knncor.med=F){
  factor.identities <- dfs.per.factor %>% names
  if (!is.entropy & !is.knncor &!is.knncor.med){
    plots <- dfs.per.factor %>% Map(function(x,factor.identity){
      vars <- colnames(x)
      if (jitter){
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_jitter(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_point(alpha=alpha) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    }, ., factor.identities)
  } else if (is.entropy) {
    split.by.resolution <- dfs.per.factor %>% lapply(function(x){split(x, x$leiden.resolution)})
    idents.for.resolutions <- factor.identities %>% lapply(rep,4)
    EntropyPlot <- function(x, factor.identity){
      vars <- colnames(x)
      which.res <- x$leiden.resolution %>% unique
      if (jitter){
        p <- x %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_jitter(alpha=alpha)+ggtitle(paste('leiden resolution: ', which.res)) + theme(legend.position="top")
      } else {
        p <- x %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_point(alpha=alpha)+ggtitle(paste('leiden resolution: ', which.res)) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    }
    plots <- split.by.resolution %>% Map(function(x, idents){
      #browser()
      plots.for.a.res <- x %>% Map(EntropyPlot, ., idents)
    }, ., idents.for.resolutions)
  } else if (is.knncor) {
    dfs.per.factor <- dfs.per.factor %>% Map(function(x, factor.identity){x[[factor.identity]] <- x[[factor.identity]] %>% as.factor; return(x)},
                                             ., factor.identities)
    plots <- dfs.per.factor %>% Map(function(x, factor.identity){
      vars <- colnames(x)
      if (!jitter){
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_boxplot(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% filter(de.level!='ref') %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_violin(alpha=alpha) + theme(legend.position="top")+stat_summary(fun.y=median, geom="point",
                                                                               position=position_dodge(width=0.9))
      }
      if (geom.smooth) {
        #p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
        NULL
      }
      return(p)
    }, ., factor.identities)
  } else if (is.knncor.med) {
    plots <- dfs.per.factor %>% Map(function(x, factor.identity){
      vars <- colnames(x)
      if (jitter){
        p <- x %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_jitter(alpha=alpha) + theme(legend.position="top")
      } else {
        p <- x %>% ggplot(aes_string(x=factor.identity, y=vars[1], col='de.level'))+
          geom_point(alpha=alpha) + theme(legend.position="top")
      }
      if (geom.smooth) {
        p <- p+stat_smooth(geom='line', se=F, alpha=alpha)
      }
      return(p)
    }, ., factor.identities)
  }
  return(plots)
}

doPlotsPerFactor <- function(dfs.per.distance, alpha=0.6, jitter=F, geom.smooth=F, use.old=F){
  if(use.old){
    PlotsPerFactor <- PlotsPerFactor_old
  }
  plots.per.distance <- names(dfs.per.distance) %>% lapply(function(x){
    if (x=='entropy'){
      plots <- PlotsPerFactor(dfs.per.distance[[x]], alpha=alpha, jitter=jitter, geom.smooth=geom.smooth, is.entropy=T)
      plots <- plots %>% unlist(recursive=F)
    } else if (x=='knncor.z'){
      plots <- PlotsPerFactor(dfs.per.distance[[x]], alpha=alpha, jitter=jitter, geom.smooth=geom.smooth, is.knncor=T)
    } else if (x=='knncor.z.med'){
      plots <- PlotsPerFactor(dfs.per.distance[[x]], alpha=alpha, jitter=jitter, geom.smooth=geom.smooth, is.knncor.med=T)
    } else {
      plots <- PlotsPerFactor(dfs.per.distance[[x]], alpha=alpha, jitter=jitter, geom.smooth=geom.smooth)
    }
  })
  names(plots.per.distance) <- names(dfs.per.distance)
  return(plots.per.distance)
}

medianofz <- function(z.df, factor.identity){
  z.df %<>% dplyr::rename(some.factor:=factor.identity)
  summarised.df <- z.df %>% group_by(some.factor, de.level) %>% summarise(distance.median=median(knncor.z.score),
                                                                         distance.mad=mad(knncor.z.score))
  summarised.df %<>% dplyr::rename(!!factor.identity:=some.factor)
  return(summarised.df)
}

CreateGrid <- function(plots.per.factor, leiden.resolutions, row.names=NULL){
  if (is.null(row.names)){
    factor.names <- plots.per.factor %>% names
  } else {
    factor.names <- row.names
  }
  n.factors <- plots.per.factor[[1]] %>% length
  if (length(plots.per.factor$entropy)>n.factors){
    extra.ents <- rep('entropy', (length(leiden.resolutions)-1)) %>% as.list
    factor.names <- list(factor.names, extra.ents) %>% unlist(recursive=F)
  }
  GroupEnts <- function(ent.plots){
    #browser()
    n.ent <- length(plots.per.factor$entropy)
    split.inds <- 1:n.ent %>% split(seq(1,n.ent, by=n.factors))
    ent.groups <- split.inds %>% lapply(function(x){plots.per.factor$entropy[x]})
    names(ent.groups) <- leiden.resolutions
    return(ent.groups)
  }
  #browser()
  grouped.ents <- GroupEnts(plots.per.factor$entropy)
  plots.per.factor <- plots.per.factor[names(plots.per.factor)!='entropy']
  plots.per.factor <- list(plots.per.factor, grouped.ents) %>% unlist(recursive=F)
  sub.grids <- plots.per.factor %>% Map(function(x,y){plt.grid <- cowplot::plot_grid(plotlist=x, labels=y,ncol=3)},.,factor.names)
  grid.of.grids <- cowplot::plot_grid(plotlist=sub.grids, nrow=length(sub.grids))
}
