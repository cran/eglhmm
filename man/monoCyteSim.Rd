\name{monoCyteSim}
\alias{monoCyteSim}
\alias{bivarSim}
\alias{ccSim}
\docType{data}
\title{
Simulated monocyte counts and psychosis symptoms.
}
\description{
Discretised values of monocyte counts, and ratings of level of
psychosis simulated from a model fitted to a data set consisting of
observations made on a number of patients from the Northland District
Health Board system.  The real data must be kept confidential due
to ethics constraints.
}
\section{Notes}{
The data sets are named \code{bivarSim} and \code{cSim}.

These data are \bold{not} immediately available in the \code{eglhmm}
package.  Their presence would cause the size of the \code{data}
directory to exceed 4.5 Mb., which is unacceptably large.
Consequently these data sets have been placed in a separate
\dQuote{data only} package called \code{monoCyteSim}, which is
available from \code{github}.  This package may be obtained by
executing the command:
\preformatted{
    install.packages("monoCyteSim",repos="https://rolfturner.r-universe.dev")
}
After having installed the \code{monoCyteSim} package, you may load
it via \code{library(monoCyteSim)} and then access the data sets
in the usual way, e.g. \code{X <- ccSim}.

Alternatively (after having installed the \code{monoCyteSim}
package) you may use the \code{::} syntax to access a single data
set, e.g. \code{X <- monoCyteSim::ccSim}.

You can access the documentation via, e.g., \code{?monoCyteSim::ccSim}.
}
\keyword{datasets}
