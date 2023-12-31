% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Community_cluster.R
\name{Community_cluster_dbscan}
\alias{Community_cluster_dbscan}
\title{Title}
\usage{
Community_cluster_dbscan(embedding, type, eps_value = 0.05, pts_value = 5)
}
\arguments{
\item{embedding}{the embedding results from different methods}

\item{type}{node embedding ‘n’ or network embedding ‘N’}

\item{eps_value}{parameters for DBSCAN}

\item{pts_value}{parameters for DBSCAN}
}
\value{
the embedding results
}
\description{
Title
}
\examples{
tnsr = GenerateMMSBM(200, 3, 10, 2, d = NULL, r = NULL)
U_list = InitializationMMSBM(tnsr, 3, 2, rank = NULL)
embed_list = PowerIteration(tnsr,3,2,rank=NULL,type="TUCKER",U_0_list=U_list)
em = embed_list[[2]]
Community_cluster_dbscan(em,"N")
}
