PPI2Mat <- function(ppi)
{
  net<-graph.data.frame(ppi,directed=F)
  stopifnot(is.connected(net))
  net <- simplify(net)
  mat<-get.adjacency(net)
  Matrix(mat,sparse=T)
}