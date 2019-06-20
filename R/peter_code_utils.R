FractionalPlot <- function(patient.vec, subtype.vec, condition.vec, fraction.palette=NULL, return.plot=TRUE){
  annotation <- bind_cols(list(patient.vec, subtype.vec, condition.vec)) %>% as.data.frame
  sub.split <- split(annotation, annotation[, 2])
  sub.split <- sub.split %>% lapply(function(x){x <- mutate(x, pat.cond = paste(x[, 1], x[, 3], sep='-'))})
  CountPatConds <- function(annotation) {
    patcond.col <- dim(annotation)[2]
    combs.df <- table(annotation[, patcond.col]) %>% as.data.frame
    colnames(combs.df)[1] <- 'pat.cond'
    colnames(combs.df)[2] <- 'count'
    combs.df <- combs.df %>% mutate(patient = gsub("-.*","", pat.cond))
    combs.df <- combs.df %>% mutate(condition = gsub(".*-","", pat.cond))
    combs.df <- combs.df %>% mutate(subtype = annotation[, 2][1:nrow(combs.df)])
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

  p <- ggplot(na.omit(freq.df),aes(x=subtype,y=freq,dodge=condition,fill=condition))+
    geom_boxplot(notch=FALSE,outlier.shape=NA) +
    geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=condition),size=0.8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    xlab("") +ylab("fraction of total cells")+ theme(legend.position="top")
  if (fraction.palette %>% is.null==FALSE) {
    p <- p+scale_fill_manual(values=fraction.palette)
  }
  if (return.plot==TRUE) {
    return(p)
  } else {
    return(freq.df)
  }
}

Makesubdistmat <- function(con.object, sample.vec, subtype.vec, cellid.vec) {
  type.factor <- setNames(subtype.vec, cellid.vec) %>% as.factor
  sample.factor <- setNames(sample.vec, cellid.vec) %>% as.factor
  CM <- con.object$getClusterCountMatrices(groups=type.factor)
  type.factor.count <- table(type.factor, sample.factor)

  sub.dist.mat <- parallel::mclapply(setNames(colnames(CM[[1]]), colnames(CM[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(CM,function(x) x[,ct]))
    tcm <- t(tcm/pmax(1,rowSums(tcm)))
    tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm)); # calls a c function
    # calculate how many cells there are
    attr(tcd,'cc') <- type.factor.count[ct,colnames(tcm)]
    tcd
  },mc.cores=1)

  return(sub.dist.mat)
}

PlotDistanceMatRed <- function(sub.dist.mat.object, sample.vec, subtype.vec, patient.vec,
                               cellid.vec, condition.vec, perplexity, max_iter=1e3, get.mat=FALSE, by.subtype=TRUE) {
  annotation <- bind_cols(list(sample.vec, subtype.vec, patient.vec, cellid.vec, condition.vec)) %>% as.data.frame
  type.factor <- setNames(subtype.vec, cellid.vec) %>% as.factor
  samp.factor <- setNames(sample.vec, cellid.vec) %>% as.factor
  type.factor.count <- table(type.factor, samp.factor)
  if (by.subtype) {
  x <- abind::abind(lapply(sub.dist.mat.object,function(x) {
    nc <- attr(x,'cc')
    wm <- sqrt(outer(nc,nc,FUN='pmin'))
    return( x*wm )
  }),along=3)
  y <- abind::abind(lapply(sub.dist.mat.object,function(x) {
    nc <- attr(x,'cc')
    sqrt(outer(nc,nc,FUN='pmin'))
  }),along=3)
  agg.dist.mat <- apply(x,c(1,2),sum)/apply(y,c(1,2),sum)
  } else {
    agg.dist.mat <- sub.dist.mat.object
  }

  sample.split <- split(annotation, sample.vec)
  sample.conds <- sample.split %>% lapply(function(x){x[, 5][1]})
  sample.pats <- sample.split %>% lapply(function(x){x[, 3][1]})

  agg.dist.mat.tsne <- Rtsne::Rtsne(agg.dist.mat,is_distance=TRUE, perplexity=perplexity,max_iter=max_iter)$Y
  df <- data.frame(agg.dist.mat.tsne); colnames(df) <- c("x","y")
  df <- df %>% mutate(samples=rownames(agg.dist.mat))
  df <- df[order(df$samples),]
  sample.conds <- sample.conds[order(sample.conds %>% names)]
  df <- df %>% mutate(condition=unlist(sample.conds), patient=unlist(sample.pats))
  df <- df %>% mutate(ncells=colSums(type.factor.count)[df$samples])

  p <- ggplot(df,aes(x,y,color=patient,shape=condition,size=log10(ncells))) + geom_point() +
    theme_bw() +xlab('tsne 1') + ylab('tsne 2') + guides(color=guide_legend(ncol=2))
  if (get.mat){
    return(agg.dist.mat)
  } else {
    return(p)
  }
}

ConditionDistanceDensity <- function(sub.dist.mat.object, sample.vec, subtype.vec, patient.vec,
                               cellid.vec, condition.vec, notch=FALSE, fraction.palette=NULL, plot=TRUE, by.subtype=TRUE) {
  type.factor <- setNames(subtype.vec, cellid.vec) %>% as.factor
  samp.factor <- setNames(sample.vec, cellid.vec) %>% as.factor

  type.sample.count <- table(type.factor, samp.factor)
  if (by.subtype) {
  x <- abind::abind(lapply(sub.dist.mat.object, function(x) {
    nc <- attr(x,'cc')
    wm <- sqrt(outer(nc,nc,FUN='pmin'))
    return( x*wm )
  }),along=3)
  y <- abind::abind(lapply(sub.dist.mat.object, function(x) {
    nc <- attr(x,'cc')
    sqrt(outer(nc,nc,FUN='pmin'))
  }),along=3)
  agg.dist.mat <- apply(x,c(1,2),sum)/apply(y,c(1,2),sum)
  } else {
    agg.dist.mat <- sub.dist.mat.object
  }
  x <- agg.dist.mat; x[upper.tri(x)] <- NA; diag(x) <- NA
  df.jsd <- na.omit(melt(x))

  annotation <- bind_cols(list(sample.vec, patient.vec, condition.vec)) %>% as.data.frame
  sample.split <- split(annotation, sample.vec)
  GetPat <- function(x){
    samp.sub <- annotation[annotation[, 1]==x, ]
    pat <- samp.sub[, 2][1]
  }
  GetCond <- function(x){
    samp.sub <- annotation[annotation[, 1]==x, ]
    pat <- samp.sub[, 3][1]
  }
  df.jsd$patient1 <- df.jsd$Var1 %>% lapply(GetPat) %>% unlist
  df.jsd$patient2 <- df.jsd$Var2 %>% lapply(GetPat) %>% unlist
  df.jsd$condition1 <- df.jsd$Var1 %>% lapply(GetCond) %>% unlist
  df.jsd$condition2 <- df.jsd$Var2 %>% lapply(GetCond) %>% unlist

  df.jsd$samePatient <- df.jsd$patient1==df.jsd$patient2
  df.jsd$sameCondition <- df.jsd$condition1==df.jsd$condition2
  df.jsd <- df.jsd[df.jsd$sameCondition==TRUE,]

  ts.within.fraction <- ggplot(na.omit(df.jsd),aes(x=condition1,y=value))+geom_boxplot(notch=notch,outlier.shape=NA,aes(fill=condition1))+
    geom_jitter(position=position_jitter(0.2), color=adjustcolor('black', alpha=0.2)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    guides(fill=FALSE) + xlab('') + ylab('inter-patient distance')
  if(fraction.palette %>% is.null()==FALSE) {
    ts.within.fraction <- ts.within.fraction + scale_fill_manual(values=fraction.palette)
  }
  if(plot){
    return(ts.within.fraction)
  } else {
    return(df.jsd)
  }
}

PlotCellTypeDists <- function(sub.dist.mat.object, min.cells, sample.vector, condition.vector, plot=T, cell.count.attr=TRUE) {
  look.up <- split(condition.vector, sample.vector) %>% lapply(unique)
  InterPatDF <- function(cell.type, min.cells=0, look.up, cell.count.attr) {
    factors <- look.up %>% unlist %>% unique
    x <- sub.dist.mat.object[[cell.type]] # jsdists for subtype xn
    x[upper.tri(x)] <- NA; diag(x) <- NA
    df2 <- melt(x)
    if (cell.count.attr){
      nc <- attr(x,'cc') # count of cells
      wm <- outer(nc,nc,FUN='pmin') # min cells for the comparison
      wm[upper.tri(wm)] <- NA; diag(wm) <- NA
      df2$ncells <- melt(wm)$value
      df2 <- na.omit(df2)
      df2 <- df2[df2$ncells>=min.cells,]
    }
    df2 <- na.omit(df2)

    df2$patient1 <- df2$Var1
    df2$patient2 <- df2$Var2
    inds.factor1.var1 <- which((df2$Var1 %>% lapply(function(x)look.up[[as.character(x)]]) %>% unlist)==factors[1])
    inds.factor1.var2 <- which((df2$Var2 %>% lapply(function(x)look.up[[as.character(x)]]) %>% unlist)==factors[1])
    inds.factor2.var1 <- which((df2$Var1 %>% lapply(function(x)look.up[[as.character(x)]]) %>% unlist)==factors[2])
    inds.factor2.var2 <- which((df2$Var2 %>% lapply(function(x)look.up[[as.character(x)]]) %>% unlist)==factors[2])
    factor1.nooverlap <- setdiff(union(inds.factor1.var1, inds.factor1.var2), intersect(inds.factor1.var1, inds.factor1.var2))
    factor2.nooverlap <- setdiff(union(inds.factor2.var1, inds.factor2.var2), intersect(inds.factor2.var1, inds.factor2.var2))
    inds.mixed <- intersect(factor1.nooverlap, factor2.nooverlap)

    df2 <- df2[inds.mixed,]
    df2$cell <- cell.type
    df2
  }
  xl <- do.call(rbind,lapply(names(sub.dist.mat.object), InterPatDF, min.cells, look.up, cell.count.attr=cell.count.attr))
  xl$comparison <- paste(xl$patient1, xl$patient2, sep='.')
  p <- ggplot(na.omit(xl),aes(x=cell,y=value))+geom_boxplot(notch=FALSE)+geom_jitter(aes(col=comparison))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    xlab("") +ylab("inter-patient expression distance")+ theme(legend.position="top")
  if (plot) {
    return(p)
  } else {
    return(xl)
  }
}
