% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhplot.R
\name{lh_eta_dist}
\alias{lh_eta_dist}
\title{DISTRIBUTION OF ETA}
\usage{
lh_eta_dist(
  data = eta1,
  lst.eta = c("Ka", "F", "Vc"),
  ncol = 3,
  nrow = 2,
  fancy = "yes"
)
}
\arguments{
\item{data}{Data frame, merged ETA and COVAR data}

\item{lst.eta}{List of ETA names}
}
\description{
Generate GOF
}
\examples{
p1<-lh_con_cov(data=cateta,lst.eta=keta,lst.cov=cat,save.path=NULL)
}
\keyword{lh_con_cov}
