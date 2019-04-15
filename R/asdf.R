#' Turn a count matrix into a panel list.
#'
#' @param count.mat count matrix
#' @param annotation annotation file
#' @param cat.col column number specifying sample identity in annotation
#' @return A panel, i.e. a list of count matrices for each sample.
Panelize <- function(count.mat, annotation, cat.col) {
  samples <- unique(annotation[,cat.col])
  # might need to sort annot and countmat if order of cellids differ between annotation and countmat
  sample.inds <- samples %>% lapply(function(cat, annotation, cat.col){return(annotation[,cat.col]==cat)}, annotation, cat.col)
  panel <- sample.inds %>% lapply(function(inds, count.mat){return(count.mat[, inds])}, count.mat)
  names(panel) = samples
  return(panel)
}

#' Get unique samples for a category in annotation
#'
#' @param cat category
#' @param annotation annotation file
#' @param sample.col the number of the sample column
#' @param cat.col the number of the category column
#' @return unique samples for a category
GetUniqSamp <- function(cat,annotation, sample.col, cat.col) {
  samples <- unique(annotation[annotation[ ,cat.col]==cat,][,sample.col])
  return(as.character(as.matrix(samples)))
}

CreateSampleGroups <- function(annotation, cat.col, sample.col) {
  sample.cats <- as.matrix(unique(annotation[ ,cat.col]))
  sample.groups <- sample.cats %>% lapply(GetUniqSamp, annotation, sample.col, cat.col)
  names(sample.groups) <- sample.cats
  return(sample.groups)
}

GetMostvar <- function(de.object,nr) {
  mostvar.names <- de.object[!is.na(de.object)] %>%
    lapply(function(df) rownames(df$res)[order(df$res$pvalue)][1:nr]) %>%
    Reduce(union, .)
}

GetClusterIDs <- function(con.object, str) {
  # will prolly change this later since very specific for clusters
  cluster.cellid <- con.object[['clusters']][[str]][['result']] %$% split(names, membership)
  return(cluster.cellid)
}

RbindPanel <- function(con.panel) {
  count.list <- con.panel$samples %>% lapply(function(x){return(x$counts)}) #%>% Reduce(rbind,.)
  return(do.call(rbind, count.list))
  #return(count.list)
}

GetColMeans <- function(cluster.array) {
  if (dim(cluster.array) %>%  is.null) {
    meaned <- cluster.array
  } else {
    meaned <- cluster.array %>% apply(2, mean)
  }
  return(meaned)
}
CategoryExper <- function(cellids, rbound, genes) {
  ReturnExps <- function(cellids) {
    if (length(cellids) < 2) {
      gene.exps <- rbound[cellids,] %>% as.matrix %>% t %>% as.data.frame
    } else {
      gene.exps <- rbound[cellids,]
    }
    return(gene.exps)
  }
  #cluster.exps <- cellids %>% lapply(function(x){return(rbound[x,])})
  cluster.exps <- cellids %>% lapply(ReturnExps)
  cluster.exps <- cluster.exps %>% lapply(function(x){return(x[,genes])})
  cluster.exps <- as.data.frame(lapply(cluster.exps,GetColMeans)) # gives genes*cells, might wanna change in future
  return(cluster.exps)
}

GetClusterExp <- function(con.panel, genes, cluster.str) {
  cell.ids <- GetClusterIDs(con.panel, cluster.str)
  rbound.panel <- RbindPanel(con.panel)
  cluster.exps <- CategoryExper(cell.ids, rbound.panel, genes)
  colnames(cluster.exps) <- names(cell.ids)
  return(cluster.exps)
}

GetSubtypeIDs <- function(annotation, sub.col, id.col) {
  annotation <- as.data.frame(annotation)
  split.annot <- annotation %$% split(annotation[,id.col],annotation[,sub.col])
  return(split.annot)
}

GetSubtypeExp <- function(con.panel, genes, annotation, sub.col, id.col) {
  cell.ids <- GetSubtypeIDs(annotation, sub.col, id.col)
  rbound.panel <- RbindPanel(con.panel)
  subtype.exps <- CategoryExper(cell.ids, rbound.panel, genes)
  colnames(subtype.exps) <- names(cell.ids)
  return(subtype.exps)
}

GetIntersect <- function(list.of.vecs) {
  return(list.of.vecs %>% lapply(function(x){x}) %>% Reduce(intersect, .))
}

FilterAnnot <- function(con.object, annotation, cellid.col) {
  count.names <- con.object$samples %>% lapply(function(x){return(rownames(x$counts))}) %>% unlist
  filt.inds <- as.data.frame(annotation)[,cellid.col] %in% count.names
  return(annotation[filt.inds,])
}

# input: annotation corresponding to one of the states, splits annotation into list of subtypes containing list of samples
SplitCells <- function(annotation, samp.col, sub.col) {
  sub.split <- split(annotation, annotation[, sub.col], drop=T)
  samp.split <- sub.split %>% lapply(function(x){split(x, x[, samp.col], drop=T)})
}

SelectCellProbs <- function(annotation, rbound.panel, cellid.col, nr.cell, genes, pseudo.count){
  if (nr.cell > dim(annotation)[1]) {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=TRUE)
  } else {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=FALSE)
  }
  annot.sampled <- annotation[inds,]
  exps <- rbound.panel[annot.sampled[, cellid.col], genes]
  probs <- exps/Matrix::rowSums(exps)
  probs <- probs+pseudo.count
  probs <- probs/Matrix::rowSums(probs)
  probs.list <- 1:nr.cell %>% lapply(function(i){return(probs[i,])})
  names(probs.list) <- rownames(probs)
  return(probs.list)
}

SelectCellProbsAggregated <- function(annotation, genes, rbound.panel, cellid.col, nr.cell, pseudo.count){
  if (nr.cell > dim(annotation)[1]) {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=TRUE)
  } else {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=FALSE)
  }
  annot.sampled <- annotation[inds,]
  exps <- rbound.panel[annot.sampled[, cellid.col], genes]
  prob.dist <- exps %>% GetColMeans
  prob.dist <- prob.dist/sum(prob.dist)
  prob.dist <- prob.dist+pseudo.count
  prob.dist <- prob.dist/sum(prob.dist)
  return(prob.dist)
}

IndividualCellProbs <- function(annotation, rbound.panel, cellid.col, sub.col, nr.cell, genes, pseudo.count=0){
  sub.split <- split(annotation, annotation[, sub.col], drop=T)
  cell.distributions <- sub.split %>% lapply(fuck::SelectCellProbs, rbound.panel, cellid.col, nr.cell, genes, pseudo.count)
}


GetSampProbs <- function(list.of.samp, rbound.panel, genes, cellid.col, pseudo.count){
  exps.list <- list.of.samp %>% lapply(function(x){rbound.panel[x[,cellid.col], genes]})
  exps.list <- exps.list %>% lapply(fuck::GetColMeans)
  exps.list <- exps.list %>% lapply(function(x){x+pseudo.count})
  probs.list <- exps.list %>% lapply(function(vector){return(vector/sum(vector))})
}

GetAllProbs <- function(annotation, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count){
  first.split <- fuck::SplitCells(annotation, samp.col, sub.col)
  prob.dist <- first.split %>% lapply(GetSampProbs, rbound.panel, genes, cellid.col, pseudo.count)
}

GetSubMatrices <- function(list.of.cats, rbound.panel, genes, cellid.col) {
  exps.list <- list.of.cats %>% lapply(function(x){rbound.panel[x[,cellid.col], genes]})
}

ObtainProbabilities <- function(annotation.list, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count=0){
  #rbound.panel <- RbindPanel(con.object)
  prob.list <- annotation.list %>% lapply(GetAllProbs, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count)
  return(prob.list)
}

CalculateJSD <- function(x, a.list) {
  a.list %>% lapply(JensenShannon,x)
}

CalculateAllJSD <- function(list1, list2) {
  alldists <- list1 %>% lapply(CalculateJSD, list2)
  return(unlist(alldists))
}


KLD <- function(P,Q) {
  divergence <- sum(P*log(P/Q))
  return(divergence)
}

JensenShannon <- function(P,Q){
  M <- (P+Q)/2
  kbl.pm <- KLD(P,M)
  kbl.qm <- KLD(Q,M)
  return((kbl.pm+kbl.qm)/2)
}

Bhattacharyya <- function(x,y){
  -log(sum(sqrt(x*y)))
}

CalculateBhat <- function(x, a.list) {
  a.list %>% lapply(Bhattacharyya,x)
}

CalculateAllBhat <- function(list1, list2) {
  alldists <- list1 %>% lapply(CalculateBhat, list2)
  return(unlist(alldists))
}

KSMatrix <- function(X,Y, pseudo.count=0){
  X <- X+pseudo.count
  Y <- Y+pseudo.count
  a_range <- 1:dim(X)[2]
  ks_stats <- a_range %>% sapply(function(i){ks.test(X[,i], Y[,i])$statistic})
}


#' Wrapper for SelectCellProbs that rearranges input order for convenience
#'
#' @param sub.annot annotation subset
#' @param some.genes some genes
#' @param rbound.panel matrix of all expression values post-pagoda'ing
#' @param cellid.col cellid column number
#' @param nr.cell number of cells
#' @param pseudo.count pseudo probability to add
#' @return A number of cells turned into probability distributions (of the genes)
ObtainDistributions <- function(sub.annot, some.genes, rbound.panel, cellid.col, nr.cell, pseudo.count){
  return(SelectCellProbs(sub.annot, rbound.panel, cellid.col, nr.cell, some.genes, pseudo.count))
}

IndividualCellProbsList <- function(annotation, rbound.panel, cellid.col, sub.col, nr.cell, genes.list, pseudo.count=0){
  sub.split <- split(annotation, annotation[, sub.col], drop=T)
  cell.distributions <- mapply(ObtainDistributions, sub.split, genes.list, MoreArgs=list(rbound.panel, cellid.col, nr.cell, pseudo.count), SIMPLIFY=F)
  return(cell.distributions)
}

ProbsToJSD <- function(state.split, rbound.panel, cellid.col, sub.col, nr.cell, genes.list, pseudo.count){
  list.of.lists <- state.split %>% lapply(IndividualCellProbsList, rbound.panel, cellid.col, sub.col,
                                          nr.cell, genes.list, pseudo.count)
  all.sc.dists.de <- Map(CalculateAllJSD, list.of.lists$healthy, list.of.lists$epilepsy)
  return(all.sc.dists.de %>% as_tibble)
} # specific for this data but saves some lines

IndividualCellProbsAgg <- function(annotation, rbound.panel, cellid.col, sub.col, nr.cell, genes.list, pseudo.count=0){
  sub.split <- split(annotation, annotation[, sub.col], drop=T)
  cell.distributions <- mapply(SelectCellProbsAggregated, sub.split, genes.list, MoreArgs=list(rbound.panel, cellid.col, nr.cell, pseudo.count), SIMPLIFY=F)
  return(cell.distributions)
}
