projPV<-function(p,v,v0) {
  p<-p-v0
  v<-v-v0
  v*as.numeric(((t(v)%*%p)/t(v)%*%v))+v0
}


gradientphi<-function(A,B,C, normalize=FALSE){

  g1<-A-projPV(A,B,C)
  g1<-g1/as.numeric(t(g1)%*%g1)
  g2<-B-projPV(B,A,C)
  g2<-g2/as.numeric(t(g2)%*%g2)
  g3<-C-projPV(C,A,B)
  g3<-g3/as.numeric(t(g3)%*%g3)
  cbind(g1,g2,g3)

}

#' Title
#'
#' @param cmplx simplicial complex
#' @param flt filtration complex
#' @param X Coordinates of points
#' @param Rahm Scalar mapping to 1-simplicies
#'
#' @return Returns oneforms and points of their origin
#' @export
#'
#' @examples
WhitneyMap1f_BC <- function(cmplx,flt,X,Rahm) {
  pnul<-nrow(X) #pocet 0simp - tedy Xty index 1simp v cmplx bude (X-pnul)ty index v BB[[1]], za predpokladu serazeni cmplx
  dvaS<-which(sapply(cmplx, length)==3)
  matice<-stredy<-matrix(NA,ncol=dim(X)[2],nrow=length(dvaS))
  oneS<-which(sapply(cmplx,length)==2)
  map<-rep(-1,length(cmplx))
  map[oneS]<-1:length(oneS)

  j<-0
  for(index in dvaS) {
    j<-j+1
    ABC<-X[flt[[index]],]
    g<-gradientphi(ABC[1,],ABC[2,],ABC[3,])
    hranice=cmplx[[index]] #hranice 2simp
    hrana1=cmplx[[hranice[[1]]]] #hranice  1 simp
    s1<-which(flt[[index]] %in% hrana1)
    hrana2=cmplx[[hranice[[2]]]]
    s2<-which(flt[[index]] %in% hrana2)
    hrana3=cmplx[[hranice[[3]]]]
    s3<-which(flt[[index]] %in% hrana3)


    vysledek=Rahm[map[hranice[[1]]]]*(g[,s1[2]]-g[,s1[1]]) + Rahm[map[hranice[[2]]]]*(g[,s2[2]]-g[,s2[1]]) + Rahm[map[hranice[[3]]]]*(g[,s3[2]]-g[,s3[1]])


    matice[j,]<-vysledek/3 ##1/3 - barycenter
    stredy[j,]<-colMeans(ABC)


  }


  return(list(oneform=matice,BC=stredy))
}
