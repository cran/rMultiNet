generate_sbm<-function(n,K,r,d,membership)
{
  p=d*K/(n*(1+(K-1)*r))
  q=p*r
  B=matrix(q,K,K)
  diag(B)=p

  Z=matrix(0,nrow=n,ncol=K)
  for(i in 1:n)
  {
    Z[i,membership[i]]=1
  }

  P=Z%*%B%*%t(Z)

  A=matrix(NA,ncol=n,nrow=n)
  diag(A)=0
  for (i in 1:(n-1))
  {
    for (j in (i+1):n)
    {
      A[i,j]=rbinom(1,1,P[i,j])
      A[j,i]=A[i,j]
    }
  }
  return(A)
}


generate_tensor<-function(n,m,K,L,r,d,Pi,Pi_network_type)
{
  ###Generate Membership matrix n by m
  local_membership=matrix(NA,nrow=n,ncol=m)
  colnames(local_membership)=c(1:m)
  for(i in 1:m)
  {
    local_membership[,i]=sample(x = K, size = n, replace = TRUE, prob = Pi)
  }

  ###Assign network type for each layer
  network_type=sample(x = m, size = L, replace = TRUE, prob = Pi_network_type)

  mylist=list()

  for (i in 1:L)
  {
    mylist[[i]]=generate_sbm(n=n,K=K,r=r[i],d=d[i],membership=local_membership[,network_type[i]])
  }

  global_membership=rep(0,n)
  for(i in 1:m)
  {
    global_membership=global_membership+K^(i-1)*(local_membership[,i]-1)
  }
  indices <- c(n,n,L)
  arr<-array(as.numeric(unlist(mylist)), dim=indices)
  arrT <- as.tensor(arr)
  output<-list(arrT,global_membership,local_membership,network_type)
  return(output)
}


#' Title
#'
#' @param n the number of vertices
#' @param m the number of types of the network
#' @param L the number of layers
#' @param K the number of groups of vertices
#' @param d the average degree of the network
#' @param r the out-in ratio in each layer
#'
#' @return a list including an adjacency tensor and the generating parameters
#' @export
#' @importFrom rTensor as.tensor
#' @importFrom stats rnorm
#' @importFrom stats rbinom
#' @examples
#' GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
#

GenerateMMSBM <- function(n,m,L,K,d=NULL,r=NULL)
{

  ###Initialize the average degree list
  if (is.null(d))
  {
    d_list = abs(rnorm(L,5,5))
  }
  else
    d_list = rep(d,L)

  ###Initialize the out-in ratio of each layer list
  if (is.null(r))
  {
    r_list = rep(0.4,L)
  }
  else
    r_list = rep(r,L)
  Pi_network_type = rep(1/m,m)
  Pi_network_community = rep(1/K,K)
  temp = generate_tensor(n=n,m=m,K=K,L=L,r=r_list,d=d_list,Pi=Pi_network_community,Pi_network_type=Pi_network_type)
  arrT = temp[[1]]
  Global_node_type_power = temp[[2]]
  number_gloab_community = length(unique(Global_node_type_power))
  return(arrT)

}


