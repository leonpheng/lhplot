% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lhplot.R
\name{forest.dat}
\alias{forest.dat}
\title{Prepare dataset for Forest plot in Shiny App}
\usage{
forest.dat(
  data = t,
  parameter = c("Cmax..ng.mL.", "AUCtau..ng.h.mL."),
  catcov = c("Cohort", "WT_group"),
  stats = c("quantile(x,0.5)=mid", "quantile(x,0.05)=lower", "quantile(x,0.95)=upper",
    "length(x)=n"),
  N = T
)
}
\arguments{
\item{data}{data frame (prepare categorical covariates before hand)}

\item{parameter}{list of parameters}

\item{catcov}{list categorical covariates.}

\item{stats}{define statistics for mid, lower and upper.}

\item{N}{number of subjects to be included}
}
\description{
Generate dataset for coveffectsplot
}
\examples{
dat<-forest.dat(data=t,parameter=c("Cmax..ng.mL.","AUCtau..ng.h.mL."),catcov=c("Cohort","WT_group"),stats=c("quantile(x,0.5)=mid","quantile(x,0.05)=lower","quantile(x,0.95)=upper","length(x)=n"))
Save dataset as CSV then open shiny APP using coveffectsplot::run_interactiveforestplot(). The package could be installed from Github: devtools::install_github('smouksassi/coveffectsplot')
}
\keyword{forest.dat}
