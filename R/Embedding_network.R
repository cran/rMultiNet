#' Title
#'
#' @param network_membership the number of types of the network or the number of groups of vertices
#' @param L the number of layers
#' @param paxis the number of eigenvectors to use in the plot
#'
#' @importFrom graphics text
#' @return a plot table If the number of eigenvectors is more than two or plot the image
#' @export
#'
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
#' embed_list = PowerIteration(tnsr,3,2,rank=NULL,type="TUCKER",U_0_list=U_list)
#' Embedding_network(embed_list[[2]],10,2)
Embedding_network <- function(network_membership,L, paxis = 2)
{
  included = rep(NA,L)
  for(i in 1:L)
  {
    included[i] = i
  }
  #network_membership = decomp[["U"]][[3]]
  d_2 = dim(network_membership)
  if(paxis>2)
  {
    data <- as.data.frame(network_membership[,2:paxis])
    return(data)
  }
  else{
    plot(x = network_membership[,2],y = network_membership[,3], xlab="second eigen vector", ylab="third eigen vector")
    text(network_membership[,2],network_membership[,3],as.character(included),cex = 1,col = "red")
  }
}

