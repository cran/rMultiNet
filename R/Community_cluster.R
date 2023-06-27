#' Title
#'
#' @param embedding the embedding results from different methods
#' @param type node embedding ‘n’ or network embedding ‘N’
#' @param cluster_number the number of clusters for Kmeans
#' @importFrom stats kmeans
#' @importFrom dbscan dbscan
#' @importFrom plotly plot_ly
#' @return the embedding results
#' @export
#'
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
#' embed_list = PowerIteration(tnsr,3,2,rank=NULL,type="TUCKER",U_0_list=U_list)
#' em = embed_list[[2]]
#' Community_cluster_km(em,"N",5)
Community_cluster_km <- function(embedding,type,cluster_number)
{
  if (type == "N")
  {
    re = kmeans(embedding[,c(2,3)],cluster_number,iter.max = 100)
    plot(embedding[,c(2,3)], col=as.factor(re$cluster),main = "cluster result", xlab = "second eigen vector", ylab = "third eigen vector")
  }

  if(type == "n")
  {
    re = kmeans(embedding[,c(2,3,4)],cluster_number, iter.max = 100)
    plot_ly(x = embedding[,2], y =embedding[,3], z =embedding[,4],size=2,color=as.character(re$cluster))
  }

}


#' Title
#'
#' @param embedding the embedding results from different methods
#' @param type node embedding ‘n’ or network embedding ‘N’
#' @param eps_value parameters for DBSCAN
#' @param pts_value parameters for DBSCAN
#'
#' @return the embedding results
#' @export
#'
#' @examples
#' tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#' U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
#' embed_list = PowerIteration(tnsr,3,2,rank=NULL,type="TUCKER",U_0_list=U_list)
#' em = embed_list[[2]]
#' Community_cluster_dbscan(em,"N")
Community_cluster_dbscan <- function(embedding,type,eps_value =.05,pts_value=5)
{
  if(type == "N")
  {
    #embedding = decomp[["U"]][[3]]
    dbscan_tensor=dbscan(embedding[,2:3],eps=eps_value, minPts =pts_value)
    #png("cluster_network.png")
    plot(embedding[,2:3], col=as.factor(dbscan_tensor$cluster),main = "cluster result", xlab = "second eigen vector", ylab = "third eigen vector")
    #dev.off()
  }
  if(type == "n")
  {
    #embedding = decomp[["U"]][[1]]
    #library(igraph)
    dbscan_tensor=dbscan(embedding[,c(2,3,4)],eps=eps_value,minPts = eps_value)
    plot_ly( x = embedding[,2], y =embedding[,3], z =embedding[,4],size=2
             ,color=as.character(dbscan_tensor$cluster))

  }

}
