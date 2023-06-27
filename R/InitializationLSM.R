## Intialization for U,W
Init_UW <- function(tensor_A, n, r, M, Utrue, Wtrue, perturb = 0.1, int_type) {
  L <- tensor_A@modes[3]
  M1A <- k_unfold(tensor_A, m = 1)@data
  M3A <- k_unfold(tensor_A, m = 3)@data
  if(int_type == "spec")
  {
    U_int <- svd(M1A, nu = r)[["u"]]
    W_int <- svd(M3A, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
  if (int_type == "rand") {
    # random Initialization
    U_int <- generate_U(0, n, r) / sqrt(n)
    W_int <- svd(M3A, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
  if (int_type == "warm") {
    # warm Initialization
    per_U <- matrix(runif(n * n, -0.1, 0.1), ncol = n)
    per_W <- matrix(runif(L * L, -perturb, perturb), ncol = L)
    U_int <- svd(Utrue %*% t(Utrue) / n + per_U, nu = r)[["u"]]
    W_int <- svd(Wtrue %*% t(Wtrue) + per_W, nu = M)[["u"]]
    return(list("U_int" = U_int, "W_int" = W_int))
  }
}


#' Title
#'
#' @param gen_list a list including the adjacency tensor and the parameter of the mixture multilayer network
#' @param n the number of nodes
#' @param rank rank of U
#' @param m the number of network types
#' @param perturb the upper bound of Uniform distribution
#' @param k the number of groups of vertices
#' @param int_type the method to initialize U and W ( ‘spec’, ‘rand’ or ‘warm’)
#'
#' @importFrom rTensor tucker
#' @importFrom rTensor k_unfold
#' @return a list including the adjacency tensor, U0, W0 and tuning parameters
#' @export
#'
#' @examples
#' gen_list = GenerateMMLSM(200,3,10,2,d=NULL)
#' InitializationLSM(gen_list,200,3,2)
InitializationLSM <- function(gen_list, n, m, k, rank=NULL,perturb = 0.1, int_type = "warm")
{

  if (is.null(rank))
  {
    rank = m*k-m+1
  }
  adj_tensor <- gen_list[[1]]
  theta <- gen_list[[2]]
  tuc_res <- tucker(theta, ranks = c(rank,rank,m), max_iter = 25, tol = 1e-05)
  Utrue_scaled <- tuc_res$U[[1]]
  Wtrue <- tuc_res$U[[3]]
  delta_1 <- sqrt(max(rowSums(Utrue_scaled^2)))
  delta_3 <- sqrt(max(rowSums(Wtrue^2)))
  delta <- c(delta_1,delta_1,delta_3)
  initial_UW <- Init_UW(adj_tensor,n,rank,m, Utrue_scaled, Wtrue, perturb, int_type)
  U_int <- initial_UW$U_int
  W_int <- initial_UW$W_int;
  return(list(adj_tensor,U_int, W_int, delta))
}
