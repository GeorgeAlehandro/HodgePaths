#' Do DeRahm map from 0-simplicies to 1-simplicies
#'
#' @param B Boundary matrix of 1-simplicies
#' @param X Coordinates of the points
#' @param V Vectors associated with the each point
#'
#' @return DeRahm map
#' @export
#'
#' @examples
deRahmMap1f <- function(B,X,V) { #deRahmovo zobrazeni
  MV <- abs(Matrix::t(B))%*%V/2
  Msigma <- -Matrix::t(X)%*%B
  Vysledek <- Matrix::colSums(Matrix::t(MV)*Msigma)
  return(Vysledek)
}
