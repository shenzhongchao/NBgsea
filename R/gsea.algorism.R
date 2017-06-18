GSEA.EnrichmentScore2 <- function(gene.list, gene.set,correl.vector) {  
  #   Inputs:
  #   gene.list: The ordered gene list (e.g. integers indicating the original position in the input dataset)  
  #   gene.set: A gene set (e.g. integers indicating the location of those genes in the input dataset) 
  #   correl.vector: A vector with the coorelations (e.g. signal to noise scores) corresponding to the genes in the gene list 
  
  N <- length(gene.list) 
  Nh <- length(gene.set) 
  Nm <-  N - Nh 
  
  loc.vector <- vector(length=N, mode="numeric")
  peak.res.vector <- vector(length=Nh, mode="numeric")
  tag.correl.vector <- vector(length=Nh, mode="numeric")
  tag.diff.vector <- vector(length=Nh, mode="numeric")
  tag.loc.vector <- vector(length=Nh, mode="numeric")
  loc.vector[gene.list] <- seq(1, N)
  tag.loc.vector <- loc.vector[gene.set]
  tag.loc.vector <- sort(tag.loc.vector, decreasing = F)
  tag.correl.vector <- correl.vector[tag.loc.vector]
  #tag.correl.vector <- rep(1, Nh)
  norm.tag <- 1.0/sum(tag.correl.vector)
  tag.correl.vector <- tag.correl.vector * norm.tag
  norm.no.tag <- 1.0/Nm
  tag.diff.vector[1] <- (tag.loc.vector[1] - 1) 
  tag.diff.vector[2:Nh] <- tag.loc.vector[2:Nh] - tag.loc.vector[1:(Nh - 1)] - 1
  tag.diff.vector <- tag.diff.vector * norm.no.tag
  peak.res.vector <- cumsum(tag.correl.vector - tag.diff.vector)
  max.ES <- max(peak.res.vector)
  max.ES
}

maxESs<-function(L.ordered,ord,genesets)
{
  genesetNum<-length(genesets)
  maxESs<-rep(0,genesetNum)
  for(i in 1:genesetNum)
  {
    gset<-genesets[[i]]
    maxESs[i]<-GSEA.EnrichmentScore2(ord,gset,L.ordered)
    
  }
  maxESs
}