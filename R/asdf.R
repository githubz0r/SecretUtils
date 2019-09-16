#' Turn a count matrix into a panel list.
#'
#' @param count.mat count matrix
#' @param annotation annotation file
#' @param cat.col column number specifying sample identity in annotation
#' @return A panel, i.e. a list of count matrices for each sample.
PanelizeOld <- function(count.mat, annotation, cat.col) {
  samples <- unique(annotation[,cat.col])
  # might need to sort annot and countmat if order of cellids differ between annotation and countmat
  sample.inds <- samples %>% lapply(function(cat, annotation, cat.col){return(annotation[,cat.col]==cat)}, annotation, cat.col)
  panel <- sample.inds %>% lapply(function(inds, count.mat){return(count.mat[, inds])}, count.mat)
  names(panel) = samples
  return(panel)
}

Panelize <- function(a.matrix, cell.ids, sample.ids) { # replace old panelize with this
  inds <- split(cell.ids, sample.ids)
  submats <- inds %>% lapply(function(x){a.matrix[x,]})
  return(submats)
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

RbindRaw <- function(con.object){
  count.list <- con.object$samples %>% lapply(function(x){return(x$misc$rawCounts)})
  return(do.call(rbind, count.list))
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

SelectCellProbs <- function(annotation, rbound.panel, cellid.col, nr.cell, genes, pseudo.count, norm.var=FALSE, vec=FALSE){
  if (nr.cell > dim(annotation)[1]) {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=TRUE)
  } else {
    inds <- sample(1:dim(annotation)[1], nr.cell, replace=FALSE)
  }
  annot.sampled <- annotation[inds,]
  exps <- rbound.panel[annot.sampled[, cellid.col], genes]
  if (norm.var) {
    sds <- apply(exps, 2, sd)
    sds <- sds+1e-8
    exps <- sweep(exps, 2, sds, '/')
  }
  if (vec){
    exps.list <- 1:nr.cell %>% lapply(function(i){return(exps[i,])})
    names(exps.list) <- rownames(exps)
    return(exps.list)
  } else {
    probs <- exps/Matrix::rowSums(exps)
    probs <- probs+pseudo.count
    probs <- probs/Matrix::rowSums(probs)
    probs.list <- 1:nr.cell %>% lapply(function(i){return(probs[i,])})
    names(probs.list) <- rownames(probs)
    return(probs.list)
  }
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

IndividualCellProbs <- function(annotation, rbound.panel, cellid.col, sub.col, nr.cell, genes,
                                pseudo.count=0, norm.var=FALSE, vec=FALSE){
  sub.split <- split(annotation, annotation[, sub.col], drop=T)
  cell.distributions <- sub.split %>% lapply(SelectCellProbs, rbound.panel, cellid.col,
                                             nr.cell, genes, pseudo.count, norm.var, vec)
}


GetSampProbs <- function(list.of.samp, rbound.panel, genes, cellid.col, pseudo.count){
  exps.list <- list.of.samp %>% lapply(function(x){rbound.panel[x[,cellid.col], genes]})
  exps.list <- exps.list %>% lapply(GetColMeans)
  exps.list <- exps.list %>% lapply(function(x){
    x<-x/sum(x)
    x<-x+pseudo.count}) # missing the first 'normalization' step
  probs.list <- exps.list %>% lapply(function(vector){return(vector/sum(vector))})
}

GetAllProbs <- function(annotation, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count){
  first.split <- SplitCells(annotation, samp.col, sub.col)
  prob.dist <- first.split %>% lapply(GetSampProbs, rbound.panel, genes, cellid.col, pseudo.count)
}

GetSubMatricesOld <- function(list.of.cats, rbound.panel, genes, cellid.col, avg=F) {
  exps.list <- list.of.cats %>% lapply(function(x){rbound.panel[x[,cellid.col], genes]})
  if (avg){
    return(exps.list %>% lapply(GetColMeans))
  } else {
    return(exps.list)
  }
}

GetSubMatrices <- function(subtype.vector, cellid.vector, condition.vector, count.matrix, genes, avg=T) {
  sub.cell.df <- dplyr::bind_cols(subtype.vector=subtype.vector, cellid.vector=cellid.vector)
  condition.split <- split(sub.cell.df, condition.vector)
  subtype.splits <- condition.split %>% lapply(function(x){split(x, x$subtype.vector)})
  obtainSubMats <- function(subtype.cellid.annot, count.matrix, genes, avg){
    sub.cond.vals <- count.matrix[subtype.cellid.annot$cellid.vector, genes]
    if (avg) {
      sub.cond.vals <- sub.cond.vals %>% Matrix::colMeans()
    }
    return(sub.cond.vals)
  }
  sub.mats <- subtype.splits %>% lapply(function(x){x %>% lapply(obtainSubMats, count.matrix, genes, avg)})
  return(sub.mats)
}

GetSubSampMats <- function(annotation, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count=0){
  first.split <- SplitCells(annotation, samp.col, sub.col)
  GetSampMats <- function(list.of.samp, rbound.panel, genes, cellid.col, pseudo.count){
    exps.list <- list.of.samp %>% lapply(function(x){rbound.panel[x[,cellid.col], genes]})
    exps.list <- exps.list %>% lapply(GetColMeans)
  }
  prob.dist <- first.split %>% lapply(GetSampMats, rbound.panel, genes, cellid.col, pseudo.count)
}
ObtainSubSampleMats <- function(annotation.list, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count=0){
  prob.list <- annotation.list %>% lapply(GetSubSampMats, rbound.panel, samp.col, sub.col, cellid.col, genes, pseudo.count)
  return(prob.list)
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

CalculateAllCor <- function(list1, list2) {
  CalculateCor <- function(x, a.list) {
    a.list %>% lapply(function(a,b){return(1-cor(a,b))},x)
  }
  alldists<- list1 %>% lapply(CalculateCor, list2)
  return(unlist(alldists))
}

CalculateAllBhat <- function(list1, list2) {
  CalculateBhat <- function(x, a.list) {
    a.list %>% lapply(Bhattacharyya,x)
  }
  alldists <- list1 %>% lapply(CalculateBhat, list2)
  return(unlist(alldists))
}

KLD <- function(P,Q) {
  divergence <- sum(P*log(P/Q))
  return(divergence)
}

JensenShannon <- function(P,Q){
  M <- (P+Q)/2
  kbl.pm <- sum(P*log(P/M))
  kbl.qm <- sum(Q*log(Q/M))
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

FractionalChanges <- function(annotation, patient.col, subtype.col, condition.col){
  sub.split <- split(annotation, annotation[, subtype.col])
  sub.split <- sub.split %>% lapply(function(x){x <- mutate(x, pat.cond = paste(x[, patient.col], x[, condition.col], sep='-'))})
  CountPatConds <- function(annotation) {
    patcond.col <- dim(annotation)[2]
    combs.df <- table(annotation[, patcond.col]) %>% as.data.frame
    colnames(combs.df)[1] <- 'pat.cond'
    colnames(combs.df)[2] <- 'count'
    combs.df <- combs.df %>% mutate(patient = gsub("-.*","", pat.cond))
    combs.df <- combs.df %>% mutate(condition = gsub(".*-","", pat.cond))
    combs.df <- combs.df %>% mutate(subtype = annotation[, subtype.col][1:nrow(combs.df)])
  }
  sub.split.counts <- sub.split %>% lapply(CountPatConds)
  sub.split.df <- do.call(rbind, sub.split.counts)
  rownames(sub.split.df) <- NULL
  patconds.split <- split(sub.split.df, sub.split.df$pat.cond)
  NormalizeCounts <- function(x){
    x$count <- x$count/sum(x$count)
    return(x)
  }
  patconds.split <- patconds.split %>% lapply(NormalizeCounts)
  freq.df <- do.call(rbind, patconds.split)
  freq.df <- freq.df %>% mutate(freq=count, count=NULL)
  return(freq.df)
}

GetSubMatsList <- function(list.of.annots, rbound.panel, list.of.genes, cellid.col) {
  GetMat <- function(sub.annot, genes, rbound.panel, cellid.col) {
    cellids <- sub.annot[, cellid.col]
    return(rbound.panel[cellids, genes])
  }
  exps.list <- mapply(GetMat, list.of.annots, list.of.genes, MoreArgs=list(rbound.panel, cellid.col))
  return(exps.list)
}

CalculateJBLD <- function(cov1, cov2) {
  JBLD <- log(det((cov1+cov2)/2)) -0.5*log(det(cov1%*%cov2))
  return(JBLD)
}

AddPseudo <- function(matrix, pseudo.count=1e-4){
  matrix <- as.matrix(matrix)
  return(matrix+pseudo.count)
}

countzerocols <- function(mat){
  mat <- as.matrix(mat)
  n.zerocol <- sum((mat %>% apply(2, sum))==0)
  return(n.zerocol)
}

removezerocols <- function(mat) {
  mat <- as.matrix(mat)
  mat <- mat[, !apply(mat, 2, sum)==0]
  return(mat)
}

CMD <- function(cor1, cor2) {
  vec1 <- as.vector(cor1)
  vec2 <- as.vector(cor2)
  l2norm <- function(x){
    return(sqrt(sum(x^2)))
  }
  cmd <- 1 - (vec1 %*% vec2)/(l2norm(vec1)*l2norm(vec2))
  return(cmd)
}

ComputeCor <- function(mat1, mat2) { # makes sure both matrices have the same columns
  commoncols <- intersect(colnames(mat1), colnames(mat2))
  mat1 <- mat1[, commoncols]
  mat2 <- mat2[, commoncols]
  return(list(mat1, mat2) %>% lapply(cor))
}

PadGenesRows <- function(a.matrix, gene.vector) {
  N <- gene.vector %>% length
  Nold <- dim(a.matrix)[1]
  M <- dim(a.matrix)[2]
  new.mat <- matrix(0L, nrow = N, ncol = M)
  new.mat[1:Nold, 1:M] <- a.matrix %>% as.vector
  missing.genes<- setdiff(gene.vector, rownames(a.matrix))
  name.vector <- c(rownames(a.matrix), missing.genes)
  colnames(new.mat) <- colnames(a.matrix); rownames(new.mat) <- name.vector
  return(new.mat[gene.vector,])
}

isNeg <- function(x){
  if (x<0){return ('Neg')} else {return('Pos')}
}

GetSubMats <- function(count.matrix, cellid.vector, subtype.vector, condition.vector, genes=NULL, avg=T, sum=F,
                       normalize=F, pseudo.prob=0){
  bound.annot <- dplyr::bind_cols(cellid=cellid.vector, subtype=subtype.vector, condition=condition.vector)
  condition.split <- bound.annot %>% split(bound.annot$condition)
  sub.cms.split <- condition.split %>% lapply(function(x){
    sub.split <- x %>% split(x$subtype)
    if (!is.null(genes)){
      cms <- sub.split %>% lapply(function(x){count.matrix[x$cellid, genes]})
    } else {
      cms <- sub.split %>% lapply(function(x){count.matrix[x$cellid, ]})
    }
    if (avg){
      values <- cms %>% lapply(Matrix::colMeans)
    } else if (sum) {
      values <- cms %>% lapply(Matrix::colSums)
    } else {
      values <- cms
    }
    if (normalize){
      values <- values %>% lapply(function(x){
        x <- x/sum(x)
        if(pseudo.prob>0){x<-x+pseudo.prob; x<-x/sum(x)}
        return(x)
      })
    }
    return(values)
  })
}
