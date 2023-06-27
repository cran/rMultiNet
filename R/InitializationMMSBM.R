#' Title  A function for initialization
#'
#' @param m the number of types of the network
#' @param k the number of groups of vertices
#' @param tnsr  the tensor of network
#' @param rank the rank of the core tensor calculated by the equation
#'
#' @return U_list a list including the core tensor Z, network embedding and node embedding
#' @export
#' @importFrom rTensor rs_unfold
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
InitializationMMSBM<-function(tnsr, m, k, rank=NULL)
{
  if (is.null(rank))
  {
    rank = m*k-m+1

  }
  ranks = c(rank,rank,m)
  num_modes <- tnsr@num_modes
  U_list <- vector("list", num_modes)
  temp_mat <-matrix(0,ncol=tnsr@modes[1], nrow=tnsr@modes[2])
  for(i in 1:tnsr@modes[3])
  {
    temp_mat=temp_mat+tnsr@data[,,i]
  }
  U_list[[1]] <- eigen(temp_mat,symmetric = T)$vector[,c(1:ranks[1])]
  U_list[[2]] <- U_list[[1]]
  outer_production=NULL
  for(i in 1:dim(U_list[[1]])[1])
  {
    row_matrix=NULL
    for (j in 1:dim(U_list[[1]])[2])
    {
      temp=U_list[[1]][i,j]*U_list[[2]]
      row_matrix=cbind(row_matrix,temp)
    }
    outer_production=rbind(outer_production,row_matrix)
  }
  temp_mat <- rs_unfold(tnsr, m = 3)@data %*% outer_production
  U_list[[3]] <- svd(temp_mat, nu = ranks[3])$u

  return(U_list)
}
