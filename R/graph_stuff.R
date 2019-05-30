GetPagaMatrix <- function(dst.matrix, membership.vector, scale=F) {
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
  cg1 <- contract(g, vc$membership, vertex.attr.comb = 'sum')
  cg2 <- simplify(cg1, remove.multiple = F, remove.loops = T,
                  edge.attr.comb = igraph_opt("sum"))
  inter.es <- as_adj(cg2)
  inter.es2 <- inter.es + Matrix::t(inter.es)
  es <- es.inner.cluster + inter.es2 %>% apply(1, sum)
  inter.es3 <- inter.es2 + Matrix::t(inter.es2)
  connectivities <- inter.es3
  inter.es4 <- inter.es3 %>% as('dgTMatrix')
  expected.random.null.vals <- (es[inter.es4@i + 1]*ns[inter.es4@j + 1] + es[inter.es4@j + 1]*ns[inter.es4@i + 1])/(n - 1)
  scaled.values <- ifelse(expected.random.null.vals != 0, inter.es4@x / expected.random.null.vals, 1) %>% as.numeric
  scale <- scale %>% rep(length(scaled.values))
  scaled.values <- ifelse(scale, scaled.values[scaled.values>1]<-1, scaled.values)
  connectivities <- inter.es4
  connectivities@x <- scaled.values
  return(connectivities)
}
