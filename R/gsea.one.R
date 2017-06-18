gsea.one<-function(L,genesets,n=1000)
{ 
  ord<-order(L,decreasing=T)
  L.ordered<-sort(L,decreasing=T)
  names(ord)<-names(L.ordered)
  genesets2<-sapply(genesets,function(e) ord[e])
  randomESset<-matrix(0,length(genesets),n)#逐列计算
  ESset<-maxESs(L.ordered,ord,genesets2)
  for(i in 1:n)#permutation
  {
    artificial_genesets<-sapply(genesets2,function(e) 
    {
      len<-length(e)
      sample(ord,len)
    }
    )
    randomESset[,i]<-maxESs(L.ordered,ord,artificial_genesets)
    print(i)
  }
  result<-cbind(ESset,randomESset)
  rownames(result)<-names(genesets)
  rowMeans(result[,2:(n+1)]>result[,1])
}