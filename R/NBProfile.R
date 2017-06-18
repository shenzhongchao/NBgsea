NBProfile <- function(geneSet,RWmatrix)
{
  allGenes <- rownames(mat)
  geneSet <- intersect(geneSet,allGenes)
  v <- rep(0,length(allGenes))
  names(v) <- allGenes
  v[geneSet] <- 1
  v <- v/sum(v)
  geneP <- v%*%mat
  geneP <- as.vector(geneP)
  names(geneP) <- allGenes
  geneP
}