randomwalk <- function(mat,run=1,stay=NULL)
{
  if(!is(mat,"sparseMatrix"))
    warning("'mat' is better to be sparse Matrix that can be constructed by Matrix package.")
  if(is.null(stay))
  {
    diag(mat)<-1
    mat<-mat/rowSums(mat)
  }else
  {
    if(stay<0|stay>1) stop("stay should be in [0,1]")
    mat<-mat/rowSums(mat)*(1-stay)
    diag(mat)<-stay
  }
  Ls <- list()
  for(i in 1:run)
    Ls[[i]] <- mat
  Reduce("%*%",Ls)
}