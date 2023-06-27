similarity <- function(data){
  dist <- as.matrix(dist(data) )
  sigma <- stats::median(dist)
  dist <- exp(-dist^2/(2*sigma^2))

  similarity <- dist
  return(similarity)
}

#' Title
#'
#' @param tnsr the adjacency tensor
#' @param rank the number of columns of the output matrix U
#' @param embedding_type SumAdj for ‘Node’ and M3SC for ‘Layer’
#'
#' @return The embeddding result can be applied in cluster methods like kmeans.
#' @importFrom rTensor modeSum
#' @importFrom Matrix sparseMatrix
#' @importFrom geigen geigen
#' @export
#'
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' emb_result = SpecClustering(tnsr,3)
#'
SpecClustering <- function(tnsr, rank, embedding_type = "Layer")
{
  if(embedding_type == "Layer")
  {
    k = k_unfold(tnsr, m = 3)
    matrx = k@data
  }
  if(embedding_type == "Nodes")
  {
    k = modeSum(tnsr, m = 3, drop = TRUE)
    matrx = k@data
  }
  # Given n by n similarity this function first calculate the Laplacian
  # matrix L then generate n by ncol matrix U of top ncol eigenvectors of L.
  #
  # Args:
  #     matrx: an n by n matrix
  #
  #     rank: number of columns of the output matrix U
  #
  # Returns:
  #    U: n by ncol numeric matrix that contains the ncol tops
  #       eigenvectors of Laplacian matrix as column
  #
  #calculate degree matrix
  similarity <- similarity(matrx)
  diag <- apply(similarity, 1 , sum)
  l <- length(diag)
  D <- sparseMatrix(1:l,1:l,x=diag)
  L <- D - similarity
  # Compute Normalized Laplacian
  L <- as.matrix(L)
  D <- as.matrix(D)
  eig <- geigen(L,D, symmetric=T)
  rm(L,D)
  U <- eig$vectors[,1:rank]
  return(U)

}
