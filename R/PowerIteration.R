
norm_vec <- function(x) sqrt(sum(x^2))
reg_vec <- function(x,delta) min(delta,norm_vec(x))/norm_vec(x)*x

#' Title
#'
#' @param tnsr the adjacency tensor of the network
#' @param type specifies the iterative algorithm to run ‘TWIST’ or ‘Tucker’
#' @param U_0_list InitializationMMSBM outputs
#' @param delta1 tuning parameters for regularization in mode1
#' @param delta2 tuning parameters for regularization in mode2
#' @param max_iter the max times of iteration
#' @param tol the convergence tolerance
#' @param m the number of types of the network
#' @param k the number of groups of vertices
#' @param rank the rank of the core tensor calculated by the equation
#'
#' @importFrom rTensor fnorm
#' @importFrom rTensor ttl
#' @importFrom rTensor ttm
#' @importFrom utils setTxtProgressBar
#' @importFrom utils tail
#' @importFrom utils txtProgressBar
#' @return a list including the core tensor Z, network embedding and node embedding
#' @export
#'
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
#' embed_list = PowerIteration(tnsr,3,2,rank=NULL,type="TUCKER",U_0_list=U_list)
PowerIteration<- function(tnsr, m, k, rank=NULL, type="TWIST", U_0_list, delta1=1000, delta2=1000, max_iter = 5, tol = 1e-05)
{

  if (is.null(rank))
  {
    rank = m*k-m+1

  }
  ranks = c(rank,rank,m)
  if (sum(ranks > tnsr@modes) != 0)
    stop("ranks must be smaller than the corresponding mode")
  if (sum(ranks <= 0) != 0)
    stop("ranks must be positive")
  if(type == "TWIST")
  {
    num_modes <- tnsr@num_modes
    U_list <- U_0_list
    tnsr_norm <- fnorm(tnsr)
    curr_iter <- 1
    converged <- FALSE
    fnorm_resid <- rep(0, max_iter)
    CHECK_CONV <- function(Z, U_list) {
      est <- ttl(Z, U_list, ms = 1:num_modes)
      curr_resid <- fnorm(tnsr - est)
      fnorm_resid[curr_iter] <<- curr_resid
      if (curr_iter == 1)
        return(FALSE)
      if (abs(curr_resid - fnorm_resid[curr_iter - 1])/tnsr_norm <
          tol)
        return(TRUE)
      else {
        return(FALSE)
      }
    }
    pb <- txtProgressBar(min = 0, max = max_iter, style = 3)
    while ((curr_iter < max_iter) && (!converged)) {
      message("iteration", curr_iter, "\n")
      setTxtProgressBar(pb, curr_iter)
      modes <- tnsr@modes
      modes_seq <- 1:num_modes

      ##Regularization
      U_list_reg = U_list
      for(m in modes_seq)
      {
        if(m == 1 | m == 2)
        {
          U_list_reg[[m]] = t(apply(U_list_reg[[m]],1,reg_vec, delta=delta1))
        }
        if(m == 3)
        {
          U_list_reg[[m]] = t(apply(U_list_reg[[m]],1,reg_vec, delta=delta2))
        }
      }

      ##Iterate
      for (m in modes_seq) {
        X <- ttl(tnsr, lapply(U_list_reg[-m], t), ms = modes_seq[-m])
        U_list[[m]] <- svd(rs_unfold(X, m = m)@data, nu = ranks[m])$u
      }

      Z <- ttm(X, mat = t(U_list[[num_modes]]), m = num_modes)

      if (CHECK_CONV(Z, U_list)) {
        converged <- TRUE
        setTxtProgressBar(pb, max_iter)
      }
      else {
        curr_iter <- curr_iter + 1
      }
    }
    close(pb)
    fnorm_resid <- fnorm_resid[fnorm_resid != 0]
    norm_percent <- (1 - (tail(fnorm_resid, 1)/tnsr_norm)) *
      100
    est <- ttl(Z, U_list, ms = 1:num_modes)
    invisible(list(Z = Z, U = U_list, conv = converged, est = est,
                   norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid,
                                                                   1), all_resids = fnorm_resid))
    network_embedding <- U_list[[3]]
    node_embedding <- U_list[[1]]
    return(list(Z, network_embedding, node_embedding))
  }
  if(type == "TUCKER")
  {
    decomp=tucker(tnsr,ranks,max_iter = 10000,tol=1e-05)
    node_embedding=decomp[["U"]][[1]] #nodes' embedding
    network_embedding=decomp[["U"]][[3]] #layers' embedding
    Z = decomp[["Z"]]
    return(list(Z, network_embedding, node_embedding))
  }

}
