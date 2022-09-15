#' Calculation of Hodge Decomposition of a vector field on simplicial complex
#'
#' @param BB List of boundary matricies (return of complex_to_boundaryF(...))
#' @param derahm DeRahm map on 1-simplicies
#'
#' @return Returns a list where 1st entry is gradient field, 2nd is curl field and 3rd is harmonic component
#' @export
#'
#' @examples
HodgeDecomp <- function(BB,derahm) {
  dif=Matrix::t(BB[[1]])
  codif=Matrix::t(BB[[2]])
  adjdif=BB[[1]]
  adjcodif=BB[[2]]
  A=pracma::pinv(as.matrix(adjdif%*%dif))
  s=A%*%(adjdif%*%derahm)
  yg=dif%*%s
  A2=pracma::pinv(as.matrix(codif%*%adjcodif))
  z=A2%*%(codif%*%derahm)
  yc=adjcodif%*%z
  yh=derahm-yg-yc
  return(list(yg,yc,yh))
}
