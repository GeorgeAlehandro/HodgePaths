#' Title
#'
#' @param complex filtration complex e.g. (flt$cmplx) where flt=TDA::alphaComplexFiltration(X)
#'
#' @return BB
#' @export
#'
#' @examples
build_boundary_CuR <- function(complex){
  return(build_boundary_Cu(complex))
}
