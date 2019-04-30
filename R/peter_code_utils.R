FractionalPlot <- function(annotation, patient.col, subtype.col, condition.col, fraction.palette=NULL){
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

  p <- ggplot(na.omit(freq.df),aes(x=subtype,y=freq,dodge=condition,fill=condition))+
    geom_boxplot(notch=FALSE,outlier.shape=NA) +
    geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=condition),size=0.8) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    xlab("") +ylab("fraction of total cells")+ theme(legend.position="top")
  if (fraction.palette %>% is.null==FALSE) {
    p <- p+scale_fill_manual(values=fraction.palette)
  }
  return(p)
}

Makectdm <- function(con.object, annotation, sample.col, subtype.col, cellid.col) {
  type.factor <- setNames(annotation[, subtype.col], annotation[, cellid.col]) %>% as.factor
  sample.factor <- setNames(annotation[, sample.col], annotation[, cellid.col]) %>% as.factor

  CM <- con.object$getClusterCountMatrices(groups=type.factor)
  CCT <- table(type.factor, sample.factor)

  ctdm <- mclapply(setNames(colnames(CM[[1]]), colnames(CM[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(CM,function(x) x[,ct]))
    tcm <- t(tcm/pmax(1,rowSums(tcm)))
    tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm)); # calls a c function
    # calculate how many cells there are
    attr(tcd,'cc') <- CCT[ct,colnames(tcm)]
    tcd
  },mc.cores=1)

  return(ctdm)
}

PlotDistanceMatRed <- function(ctdm.object, annotation, sample.col, subtype.col, patient.col,
                               cellid.col, condition.col, perplexity=30, max_iter=1e3) {
  type.factor <- setNames(annotation[, subtype.col], annotation[, cellid.col]) %>% as.factor
  sample.factor <- setNames(annotation[, sample.col], annotation[, cellid.col]) %>% as.factor
  CCT <- table(type.factor, sample.factor)
  x <- abind(lapply(ctdm.object,function(x) {
    nc <- attr(x,'cc');
    wm <- sqrt(outer(nc,nc,FUN='pmin'))
    return( x*wm )
  }),along=3)
  y <- abind(lapply(ctdm.object,function(x) {
    nc <- attr(x,'cc');
    sqrt(outer(nc,nc,FUN='pmin'))
  }),along=3)
  XD <- apply(x,c(1,2),sum)/apply(y,c(1,2),sum)

  sample.split <- split(annotation, annotation[, sample.col])
  sample.conds <- sample.split %>% lapply(function(x){x[, condition.col][1]})
  sample.pats <- sample.split %>% lapply(function(x){x[, patient.col][1]})

  XDE <- Rtsne::Rtsne(XD,is_distance=TRUE, perplexity=perplexity,max_iter=max_iter)$Y
  df <- data.frame(XDE); colnames(df) <- c("x","y")
  df <- df %>% mutate(samples=rownames(XD))
  df <- df[order(df$samples),]
  sample.conds <- sample.conds[order(sample.conds %>% names)]
  df <- df %>% mutate(condition=unlist(sample.conds), patient=unlist(sample.pats))
  df <- df %>% mutate(ncells=colSums(CCT)[df$samples])

  p <- ggplot(df,aes(x,y,color=patient,shape=condition,size=log10(ncells))) + geom_point() +
    theme_bw() + xlab("") + ylab("") +
    theme(axis.title=element_blank(),  axis.text=element_blank(), axis.ticks=element_blank()) +
    guides(color=guide_legend(ncol=2))
  return(p)
}

PlotCellTypeDists <- function(ctdm.object, min.cells, between.conditions=TRUE) {
  InterPatDF <- function(cell.type, min.cells=0, between.conditions) {
    x <- ctdm.object[[cell.type]] # jsdists for subtype xn
    nc <- attr(x,'cc'); # count of cells
    wm <- outer(nc,nc,FUN='pmin') # min cells for the comparison

    x[upper.tri(x)] <- NA; diag(x) <- NA;
    wm[upper.tri(wm)] <- NA; diag(wm) <- NA;
    df2 <- melt(x);
    df2$ncells <- melt(wm)$value
    df2 <- na.omit(df2)
    df2 <- df2[df2$ncells>=min.cells,];

    df2$patient1 <- df2$Var1
    df2$patient2 <- df2$Var2
    inds.c <- grep('c_',df2$Var1)
    inds.e <- grep('ep_',df2$Var2)
    if (between.conditions){
      inds.mixed <- intersect(inds.c, inds.e)
    } else {
      inds.mixed <- setdiff(union(inds.c, inds.e), intersect(inds.c, inds.e))
    }
    df2 <- df2[inds.mixed,]
    df2$cell <- cell.type
    df2
  }
  xl <- do.call(rbind,lapply(names(ctdm.object), InterPatDF, min.cells, between.conditions=TRUE))
  xl$comparison <- paste(xl$patient1, xl$patient2, sep='.')
  p <- ggplot(na.omit(xl),aes(x=cell,y=value))+geom_boxplot(notch=FALSE)+geom_jitter(aes(col=comparison))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
    xlab("") +ylab("inter-patient expression distance")+ theme(legend.position="top")
  return(p)
}
