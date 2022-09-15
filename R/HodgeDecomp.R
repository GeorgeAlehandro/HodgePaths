#' Calculation of Hodge Decomposition of a vector field on simplicial complex
#'
#' @param BB List of boundary matricies (return of complex_to_boundaryF(...))
#' @param derahm DeRahm map on 1-simplicies
#' @param alg Algorithm for lsqr possibilities: "pinv", "minres" and "bicg"
#'
#' @return Returns a list where 1st entry is gradient field, 2nd is curl field and 3rd is harmonic the component
#' @export
#'
#' @examples

HodgeDecomp <- function(BB,derahm,alg="minres") {
  dif=Matrix::t(BB[[1]])
  codif=Matrix::t(BB[[2]])
  adjdif=BB[[1]]
  adjcodif=BB[[2]]

  A1<-adjdif%*%dif
  b1<-adjdif%*%derahm

  A2<-codif%*%adjcodif
  b2<-codif%*%derahm


  if (alg=="minres") {
    res1<-minres(A1,as.numeric(b1))
    yg<-dif%*%res1$x

    res2<-minres(A2,as.numeric(b2))
    yc<-adjcodif%*%res2$x

    yh<-derahm-yg-yc
  }

  else if (alg=="pinv"){
    s<-pracma::pinv(as.matrix(A1))%*%(b1)
    yg<-dif%*%s
    z<-pracma::pinv(as.matrix(A2))%*%(b2)
    yc<-adjcodif%*%z
    yh<-derahm-yg-yc
  }

  else if (alg=="bicg") {
    res1<-bicgSparse(A1,as.numeric(b1))
    yg<-dif%*%res1$x

    res2<-bicgSparse(A2,as.numeric(b2))
    yc<-adjcodif%*%res2$x

    yh<-derahm-yg-yc
  }
  else {
    stop("Incorrect algorithm for calculation of lsqr entered")
  }



  return(list(yg,yc,yh))
}


