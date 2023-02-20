#' Create graph using transition matrix
#'
#' @param tv1 tviblindi object
#' @param n number of neighbors of the transition matrix
#'
#' @return list object
#' @export
#'
#' @examples
create_graph_from_transition_matrix <- function(tv1, n) {
  transition_matrix <- transition.matrix(tv1, n)
  G <- graph_from_adjacency_matrix(transition_matrix,weighted = TRUE,mode = "undirected")
  listOfEdges <- split(as_edgelist(G), seq(nrow(as_edgelist(G))))
  names(listOfEdges)<-NULL
  return(list(transition_matrix,G, listOfEdges))
}



#' Create the boundary matrix using the list of nodes and edges
#'
#' @param flt filter of nodes and edges
#'
#' @return boundary matrix
#' @export
#'
#' @examples
calculate_boundary_matrix <- function(flt) {
  cmplx<-build_boundary_CuR(flt)
  BB<-complex_to_boundaryF(cmplx = cmplx)
  B<-BB[[1]]
  B
}



#' Calculate the adjoint matrix
#'
#' @param tv1 tv
#' @param n number of neighbors of the transition matrix
#'
#' @return list object
#' @export
#'
#' @examples
calculate_adjoint_matrix <- function(B, D0, D1) {
  adjoint_matrix<-D0%*%B%*%D1
  return(adjoint_matrix)
}



#' Get the laplacean matrix using the boundary matrix, transition matrix and weights of the edges
#'
#' @param B Boundary matrix
#' @param D0 Diagonal matrix of the sums of all the transition possibilities of the nodes e.g: D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
#' @param D1 Diagonal matrix of the weights associated to the edges e.g: D1<-Matrix::Diagonal(x=E(G)$weight)
#'
#' @return laplacean matrix
#' @export
#'
#' @examples
calculate_Lap <- function(B, D0, D1) {
  Lap<-D0%*%B%*%D1%*%t(B)
  return(Lap)
}



#' Create the symmetrical laplacean using the boundary matrix, transition matrix and weights of the edges
#'
#' @param B Boundary matrix
#' @param D0 Diagonal matrix of the sums of all the transition possibilities of the nodes e.g: D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
#' @param D1 Diagonal matrix of the weights associated to the edges e.g: D1<-Matrix::Diagonal(x=E(G)$weight)
#'
#' @return laplacean symmetrical matrix
#' @export
#'
#' @examples
calculate_Lap_sym <- function(B, D0, D1) {
  Dm <- Matrix::Diagonal(x=diag(D0)^(1/2))
  Lap_sym<-Dm%*%B%*%D1%*%t(B)%*%Dm
  return(Lap_sym)
}



#' Calculate pseudotime with root cell selected using symmetrical laplacean
#'
#' @param B Boundary matrix
#' @param D0 Diagonal matrix of the sums of all the transition possibilities of the nodes e.g: D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
#' @param D1 Diagonal matrix of the weights associated to the edges e.g: D1<-Matrix::Diagonal(x=E(G)$weight)
#'
#' @return laplacean symmetrical matrix
#' @export
#'
#' @examples
calculate_alpha <- function(Lap, adjoint_matrix, VD, IndexOfOriginCell) {
  Dm_minus <- Matrix::Diagonal(x=diag(D0)^(-1/2))
  Dm_minus_inverse <- Matrix::Diagonal(x=diag(Dm_minus)^(-1))
  left_hand <- Lap[-IndexOfOriginCell,-IndexOfOriginCell]
  right_hand <- -(Dm_minus%*%adjoint_matrix%*%VD)[-IndexOfOriginCell]
  epsilon<-HodgePaths:::bicgSparse(left_hand,as.numeric(right_hand), nb_iter = 1500)
  alpha <- Dm_minus_inverse[-IndexOfOriginCell,-IndexOfOriginCell]%*%epsilon$x
  return(alpha)
}



#' Calculate pseudotime without root cell selected using symmetrical laplacean
#'
#' @param Lap Symmetrical laplacean matrix
#' @param D0 Diagonal matrix of the sums of all the transition possibilities of the nodes e.g: D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
#' @param D1 Diagonal matrix of the weights associated to the edges e.g: D1<-Matrix::Diagonal(x=E(G)$weight)
#'
#' @return pseudotime solution
#' @export
#'
#' @examples
calculate_alpha_no_root_selected <- function(Lap, adjoint_matrix, D0, VD) {
  Dm_minus <- Matrix::Diagonal(x=diag(D0)^(-1/2))
  Dm_minus_inverse <- Matrix::Diagonal(x=diag(Dm_minus)^(-1))
  left_hand <- Lap
  right_hand <- -(Dm_minus%*%adjoint_matrix%*%VD)
  epsilon<-HodgePaths:::bicgSparse(left_hand,as.numeric(right_hand), nb_iter = 1500)
  alpha <- Dm_minus_inverse%*%epsilon$x
  return(alpha)
}



#' Calculate pseudotime using RNA Velocity
#' This uses the symmetrical Laplacean for faster performance
#' @param tv1 tviblindi object
#' @param nearest_neighbour_number Number of the nearest neighbor nodes in the graph
#' @param IndexOfRootCell (Optional) Index of the origin cell
#' @param MatrixOfVelocity Velocity matrix (calculated using RNA Velocity tools)
#'
#' @return pseudotime
#' @export
#'
#' @examples
get_pseudotime_from_velocity <- function(tv1, nearest_neighbour_number=30, IndexOfRootCell=NULL, MatrixOfVelocity = velocity) {
  # create graph and extract edges from transition matrix
  graph <- create_graph_from_transition_matrix(tv1, nearest_neighbour_number)
  transition_matrix <- graph[[1]]
  G <- graph[[2]]
  listOfEdges <- graph[[3]]
  # calculate complex and boundary
  flt<-c(1:nrow(tv1$data),listOfEdges)
  B<-calculate_boundary_matrix(flt)

  # load data
  V <- as.matrix(MatrixOfVelocity)
  X<-tv1$data
  VD<-deRahmMap1f(B = B,X = X,V = V)
  D0<-Matrix::Diagonal(x=rowSums(transition_matrix)^-1)
  D1<-Matrix::Diagonal(x=E(G)$weight)
  adjoint_matrix <- calculate_adjoint_matrix(B,D0,D1)
  # calculate Lap and adjoint matrix
  Lap <- calculate_Lap_sym(B, D0, D1)
  # calculate alpha depending if the root cell is selected or not
  if (is.null(IndexOfRootCell)){
    alpha <- calculate_alpha_no_root_selected(Lap, adjoint_matrix, VD, D0=D0)
    return(alpha)
  }
  else{
    alpha_without_origin <- calculate_alpha(Lap, adjoint_matrix, VD, IndexOfRootCell)
    alpha <-rep(0,nrow(tv1$data))
    alpha[-IndexOfRootCell]<-as.numeric(alpha_without_origin)
    return(alpha)
  }

}
