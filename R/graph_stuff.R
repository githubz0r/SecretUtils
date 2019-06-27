GetPagaMatrix <- function(dst.matrix, membership.vector, scale=F, linearize=T) {
  if (class(dst.matrix)!='dsTMatrix'){
    dst.matrix %<>% as('dgTMatrix') %>% as('symmetricMatrix')
  }
  dst.matrix@Dimnames <- list(NULL, NULL)
  ones <- dst.matrix
  ones@x <- rep(1, length(ones@x))
  g <- graph_from_adjacency_matrix(ones, mode='directed')
  vc <- igraph::make_clusters(g, membership = membership.vector, algorithm = 'conos',
                              merges = NULL, modularity = F)
  ns <- vc %>% igraph::sizes()
  n <- sum(ns)
  getEsCount <- function(index, commu.obj, graph){
    comm.sub <- commu.obj[[index]]
    sub.graph <- igraph::induced_subgraph(graph, comm.sub, impl = 'copy_and_delete')
    return(2*(sub.graph %>% gsize))
  }
  es.inner.cluster <- 1:length(vc) %>% lapply(getEsCount, vc, g) %>% unlist
  cg1 <- igraph::contract(g, vc$membership, vertex.attr.comb = 'sum')
  cg2 <- igraph::simplify(cg1, remove.multiple = F, remove.loops = T,
                  edge.attr.comb = igraph_opt("sum"))
  inter.es <- as_adj(cg2)
  inter.es2 <- inter.es + Matrix::t(inter.es)
  es <- es.inner.cluster + inter.es2 %>% apply(1, sum)
  inter.es3 <- inter.es2 + Matrix::t(inter.es2)
  connectivities <- inter.es3
  inter.es4 <- inter.es3 %>% as('dgTMatrix')
  expected.random.null.vals <- (es[inter.es4@i + 1]*ns[inter.es4@j + 1] + es[inter.es4@j + 1]*ns[inter.es4@i + 1])/(n - 1)
  if (!linearize) {
    sd.random.null.vals <- (es[inter.es4@i + 1]*ns[inter.es4@j + 1]*(n-ns[inter.es4@j + 1]-1) +
                              es[inter.es4@j + 1]*ns[inter.es4@i + 1])*(n-ns[inter.es4@i + 1]-1) / ((n - 1)^2)
    scaled.values <- ifelse(sd.random.null.vals !=0, (inter.es4@x - expected.random.null.vals)/sd.random.null.vals, 0) %>% as.numeric
  } else {
    scaled.values <- ifelse(expected.random.null.vals != 0, inter.es4@x / expected.random.null.vals, 1) %>% as.numeric
    scale <- scale %>% rep(length(scaled.values))
    scaled.values <- ifelse(scale, scaled.values[scaled.values>1]<-1, scaled.values)
  }
  connectivities <- inter.es4
  connectivities@x <- scaled.values
  return(connectivities)
}

GenerateFactorVectors <- function(subtype.vector, sample.vector, condition.vector) {
  conc <- paste0(subtype.vector, "-;;-", sample.vector, ";__;" , condition.vector) %>% as.factor %>% levels
  subtypes <- gsub("-;;-.*", "", conc)
  samples <- gsub(".*-;;-", "", conc)
  samples <- gsub(';__;.*', "", samples)
  condition <- gsub('.*;__;', "", conc)
  return(dplyr::bind_cols(subtypes=subtypes, samples=samples, condition=condition, concatenated=conc))
}

MeltMatrix <- function(x, symmetric){
  if (symmetric) {
    x[lower.tri(x)] <- NA; diag(x) <- NA
    df <- na.omit(reshape2::melt(as.matrix(x)))
  } else {
    df <- reshape2::melt(as.matrix(x))
  }
  df <- dplyr::bind_cols(value=df$value, comparison=paste0(df$Var1, '-', df$Var2))
  return(df)
}

MeltAndAppend <- function(mat.list, factor.identity, symmetric=TRUE) {
  molten.mats <- mat.list %>% lapply(MeltMatrix, symmetric)
  AppendCols <- function(df, subtype.name, factor.identity){
    df$subtype = subtype.name
    df$condition = factor.identity
    return(df)
  }
  extended.dfs <- Map(AppendCols, molten.mats, names(mat.list), MoreArgs=list(factor.identity))
  return(extended.dfs)
}

GeneratePagaSubSampDFOld <- function(paga.connectivities, subtype.vector, sample.vector, condition.vector) {
  factor.vectors <- GenerateFactorVectors(subtype.vector, sample.vector, condition.vector)
  sub.cond.indices <- as.factor(factor.vectors$concatenated) %>% as.numeric %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.cond.indices <- sub.cond.indices[order(sub.cond.indices %>% names)]

  sub.samp.factor <- as.factor(factor.vectors$samples) %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.samp.factor <- sub.samp.factor[order(sub.samp.factor %>% names)]

  sub.cond.factor <- as.factor(factor.vectors$condition) %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.cond.factor <- sub.cond.factor[order(sub.cond.factor %>% names)]

  GetSubConnectivity <- function(indices1, indices2, connectivity.matrix){
    return(connectivity.matrix[indices1, indices2])
  }

  factor1.mats <- seq(1, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i]], paga.connectivities);
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i]]; return(sub.mat)})
  names(factor1.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  factor2.mats <- seq(2, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i]], paga.connectivities);
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i]]; return(sub.mat)})
  names(factor2.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  between.mats <- seq(1, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i+1]], paga.connectivities);
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i+1]]; return(sub.mat)})
  names(between.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  factor1.identity <- sub.cond.factor[[1]] %>% unique %>% as.character
  factor2.identity <- sub.cond.factor[[2]] %>% unique %>% as.character

  factor1.dfs <- MeltAndAppend(factor1.mats, factor1.identity)
  factor2.dfs <- MeltAndAppend(factor2.mats, factor2.identity)
  between.dfs <- MeltAndAppend(between.mats, 'between', symmetric = FALSE)
  return(dplyr::bind_rows(factor1.dfs, factor2.dfs, between.dfs))
}

GeneratePagaSubSampDF <- function(paga.connectivities, subtype.vector, sample.vector, condition.vector) {
  factor.vectors <- GenerateFactorVectors(subtype.vector, sample.vector, condition.vector)
  sub.cond.indices <- as.factor(factor.vectors$concatenated) %>% as.numeric %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.cond.indices <- sub.cond.indices[order(sub.cond.indices %>% names)]

  sub.samp.factor <- as.factor(factor.vectors$samples) %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.samp.factor <- sub.samp.factor[order(sub.samp.factor %>% names)]

  sub.cond.factor <- as.factor(factor.vectors$condition) %>%
    split(list(factor.vectors$subtypes, factor.vectors$condition))
  sub.cond.factor <- sub.cond.factor[order(sub.cond.factor %>% names)]

  GetSubConnectivity <- function(indices1, indices2, connectivity.matrix){
    return(connectivity.matrix[indices1, indices2])
  }

  factor1.mats <- seq(1, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i]], paga.connectivities) %>% as.matrix;
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i]]; return(sub.mat)})
  names(factor1.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  factor2.mats <- seq(2, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i]], paga.connectivities) %>% as.matrix;
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i]]; return(sub.mat)})
  names(factor2.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  between.mats <- seq(1, length(sub.cond.indices), 2) %>%
    lapply(function(i){sub.mat <- GetSubConnectivity(sub.cond.indices[[i]], sub.cond.indices[[i+1]], paga.connectivities) %>% as.matrix;
    rownames(sub.mat) <- sub.samp.factor[[i]]; colnames(sub.mat) <- sub.samp.factor[[i+1]]; return(sub.mat)})
  names(between.mats) <- factor.vectors$subtypes %>% as.factor %>% levels

  factor1.count <- sub.cond.factor[[1]] %>% length
  factor2.count <- sub.cond.factor[[2]] %>%  length
  factor1.identity <- sub.cond.factor[[1]] %>% unique %>% as.character
  factor2.identity <- sub.cond.factor[[2]] %>% unique %>% as.character

  if (factor1.count > 1) {
    factor1.dfs <- MeltAndAppend(factor1.mats, factor1.identity)
  } else {
    factor1.dfs <- NULL
  }
  if (factor2.count > 1) {
    factor2.dfs <- MeltAndAppend(factor2.mats, factor2.identity)
  } else {
    factor2.dfs <- NULL
  }
  between.dfs <- MeltAndAppend(between.mats, 'between', symmetric = FALSE)
  return(dplyr::bind_rows(factor1.dfs, factor2.dfs, between.dfs))
}

GenerateUnalignedAdjOld <- function(raw.count.mat, cellid.vector, k=15){
  p2 <- Pagoda2$new(Matrix::t(raw.count.mat), log.scale=F)
  p2$adjustVariance(gam.k=10)
  p2$calculatePcaReduction(nPcs=100,n.odgenes=3e3)
  p2$makeKnnGraph(k=k,type='PCA',center=T,distance='angular')
  unaligned.graph.adj <- igraph::as_adjacency_matrix(p2$graphs$PCA, attr="weight")[cellid.vector, cellid.vector]
  return(unaligned.graph.adj)
}

GenerateUnalignedAdj <- function(raw.count.mat, cellid.vector, k=15){
  p2 <- Matrix::t(raw.count.mat) %>% basicP2proc(n.cores=1, min.cells.per.gene=0, n.odgenes=3e3,
                                                 get.largevis=FALSE, make.geneknn=FALSE, get.tsne=FALSE)
  p2$makeKnnGraph(k=k,type='PCA',center=T,distance='angular')
  unaligned.graph.adj <- igraph::as_adjacency_matrix(p2$graphs$PCA, attr="weight")[cellid.vector, cellid.vector]
  return(unaligned.graph.adj)
}

SubsetPAGABySubtype <- function(connectivities, membership.levels, sample.order, subtype.order) {
  membership.split <- membership.levels %>% as.factor %>% as.numeric %>% split(subtype.order)
  sample.split <- sample.order %>% as.factor %>% split(subtype.order)
  subsetConnectivities <- function(sub.name){
    sub.inds <- membership.split[[sub.name]]
    submat <- connectivities[sub.inds, sub.inds]
    colnames(submat) <- sample.split[[sub.name]]
    rownames(submat) <- sample.split[[sub.name]]
    return(submat)
  }
  subtype.connectivities <- subtype.order %>% unique %>% lapply(subsetConnectivities)
  names(subtype.connectivities) <- subtype.order %>% unique
  return(subtype.connectivities)
}

PadMatZeroes <- function(a.matrix, sample.list) {
  N <- sample.list %>% length
  Nold <- dim(a.matrix)[1]
  new.mat <- matrix(0L, nrow = N, ncol = N)
  new.mat[1:Nold, 1:Nold] <- a.matrix %>% as.vector
  missing.samples <- setdiff(sample.list, colnames(a.matrix))
  name.vector <- c(colnames(a.matrix), missing.samples)
  colnames(new.mat) <- name.vector; rownames(new.mat) <- name.vector
  return(new.mat[sample.list, sample.list])
}




GeneratePagaItems <- function(graph.adj, subtype.vector=NULL, condition.vector=NULL, sample.vector=NULL,
                              by.subtypes.condition=FALSE, by.subtypes.samples=FALSE, by.samples=FALSE,
                              linearize=T, log.scale=F, pseudo.connectivity=1e-3){
  if (by.subtypes.condition) {
    subtype.condition <- paste0(subtype.vector, '-', condition.vector)
    membership.vector <- as.numeric(factor(subtype.condition))
    subtype.order <- (paste0(subtype.vector) %>% unique)[order(paste0(subtype.vector) %>% unique)]
    connectivities <- GetPagaMatrix(graph.adj, membership.vector, scale=F, linearize=linearize)
    statistics <- seq(1, dim(connectivities)[1], 2) %>% sapply(function(i){connectivities[i,i+1]})
    if (log.scale){
      statistics <- log(statistics+pseudo.connectivity)
    }
    paga.df <- dplyr::bind_cols(paga.connectivity.value=statistics, subtype=subtype.order)
    paga.df$ncells <- table(subtype.vector)[subtype.order] %>% as.numeric
    p <- ggplot(paga.df, aes(y=paga.connectivity.value, x=subtype)) + geom_point()+
      theme(axis.text.x = element_text(angle = -90, hjust = 1))
    return(list(connectivities=connectivities, paga.df=paga.df, scatter.plot=p))
  } else if (by.subtypes.samples) {
    subtype.samples <- paste0(subtype.vector, '-', sample.vector)
    membership.vec.subsamp <- as.numeric(factor(subtype.samples))
    connectivities <- GetPagaMatrix(graph.adj, membership.vec.subsamp, scale=F, linearize=linearize)
    paga.df <- GeneratePagaSubSampDF(connectivities, subtype.vector, sample.vector, condition.vector)
    paga.df <- paga.df %>% dplyr::rename(paga.connectivity.value=value)
    if (log.scale){
      paga.df$paga.connectivity.value <- log(paga.df$paga.connectivity.value+pseudo.connectivity)
    }
    p <- do.call(SecretUtils::GeneratePagaPlot, list(paga.df, subset=NULL))
    return(list(connectivities=connectivities, paga.df=paga.df, sub.cond.plot=p))
  } else if (by.samples) {
    membership.vector <- as.numeric(factor(sample.vector))
    sample.order <- (paste0(sample.vector) %>% unique)[order(paste0(sample.vector) %>% unique)]
    connectivities <- GetPagaMatrix(graph.adj, membership.vector, scale=F, linearize=linearize) %>% as.matrix
    rownames(connectivities) <- sample.order; colnames(connectivities) <- sample.order
    if (log.scale){
      connectivities <- log(connectivities+pseudo.connectivity)
    }
    return(list(connectivities=connectivities))
  }
}

GeneratePagaPlot <- function(paga.df, subset=NULL, geom.boxplot=T, geom.jitter=T) {
  if (!is.null(subset)){
    paga.df %<>% filter(condition==subset)
    p <- paga.df %>% ggplot(aes(x=subtype, y=paga.connectivity.value)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5))+theme(legend.position="top")
    if (geom.boxplot) {
      p <- p+geom_boxplot(notch=FALSE, outlier.shape=NA)
    }
    if (geom.jitter) {
      p <- p + geom_jitter(aes(col=comparison))
    }

  } else {
    p <- paga.df %>% ggplot(aes(x=subtype, y=paga.connectivity.value, dodge=condition, fill=condition)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(angle = 90, hjust = 0.5))+theme(legend.position="top") +
      geom_boxplot(notch=F)

  }
  return(p)
}
