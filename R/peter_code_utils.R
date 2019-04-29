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

  return(list(ctdm, CCT))
}

WeighDistances <- function(ctdm.object) {
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
  return(XD)
}

MakeXl <- function(ctdm.object, min.cells, between.conditions=TRUE) {
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
  xl <- do.call(rbind,lapply(names(ctdm.object), InterPatDF, min_cells, between.conditions=TRUE))
  xl$comparison <- paste(xl$patient1, xl$patient2, sep='.')
  return(xl)
}
