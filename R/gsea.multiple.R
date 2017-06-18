gsea.multiple<-function(Ls,genesets,n=1000)
{
  result<-matrix(0,ncol(Ls),length(genesets))
  for(i in 1:ncol(Ls))
  {
    t<-gsea.one(Ls[,i],genesets,n)
    pvalue<-rowSums(t[,1]<t[,2:n+1])/n
    result[i,]<-pvalue
  }
  rownames(result)<-colnames(Ls)
  colnames(result)<-names(genesets)
  result
}