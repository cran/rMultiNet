% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitializationLSM.R
\name{InitializationLSM}
\alias{InitializationLSM}
\title{Title}
\usage{
InitializationLSM(
  gen_list,
  n,
  m,
  k,
  rank = NULL,
  perturb = 0.1,
  int_type = "warm"
)
}
\arguments{
\item{gen_list}{a list including the adjacency tensor and the parameter of the mixture multilayer network}

\item{n}{the number of nodes}

\item{m}{the number of network types}

\item{k}{the number of groups of vertices}

\item{rank}{rank of U}

\item{perturb}{the upper bound of Uniform distribution}

\item{int_type}{the method to initialize U and W ( ‘spec’, ‘rand’ or ‘warm’)}
}
\value{
a list including the adjacency tensor, U0, W0 and tuning parameters
}
\description{
Title
}
\examples{
gen_list = GenerateMMLSM(200,3,10,2,d=NULL)
InitializationLSM(gen_list,200,3,2)
}
