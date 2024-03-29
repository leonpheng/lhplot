% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhplot.R
\name{lh_cwres}
\alias{lh_cwres}
\title{Plot cwres vs x}
\usage{
lh_cwres(
  data,
  y = "CWRES",
  x = "TAD",
  type = "log",
  scale = c(0.1, 100),
  xtit = "Individual Predicted Concentration (ng/mL)",
  ytit = "Conditional Weighted Residuals",
  sortby = NULL,
  col.obs = "#A6CEE3",
  col.ident = "#1F78B4",
  brew.col = NULL
)
}
\arguments{
\item{data}{data frame}
}
\description{
Generate plot with sort color of observe
}
\examples{
lh_cwres(...)
}
\keyword{lh_cwres}
