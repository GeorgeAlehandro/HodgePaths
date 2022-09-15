#' Title
#'
#' @param cmplx Creates boundary matrix of the filtration
#'
#' @return list of boundary matricies with out[[n]] being boundary matrix of n+1 simplicies
#' @export
#'
#' @examples
complex_to_boundaryF<-function(cmplx){
  LL<-sapply(cmplx,length)
  i<-rep(1:length(cmplx),times=LL)
  j<-unlist(cmplx[LL>1])
  x<-(-1)^sequence(LL[LL>1])
  B<-Matrix::sparseMatrix(i=j,j=i,x=x,dims=c(length(cmplx),length(cmplx)))
  LL[LL==0]<-1
  out<-list()[1:(max(LL)-1)]
  for(i in 1:(max(LL)-1)){
    out[[i]]<-B[LL==i,LL==(i+1)]
  }
  out

}
