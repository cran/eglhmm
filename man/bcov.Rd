\name{bcov}
\alias{bcov}
\title{ Bootstrap covariance. }
\description{
  Creates an estimate of the covariance matrix of the parameter
  estimates for an extended generalised linear hidden Markov model via
  parametric bootstrapping.
}
\usage{
   bcov(object, nsim = 50, itmax = 500, verbose = TRUE)
}
\arguments{
  \item{object}{ An object of class \code{eglhmm} as produced by
  \code{\link{eglhmm}()}. }
  \item{nsim}{ The number of data sets to simulate, from which
  to estimate parameters.  From each data set a vector of parameters
  is estimated; the estimated covariance matrix is the empirical
  covariance matrix of these \code{nsim} vectors.}
  \item{itmax}{ The maximum number of iterations to be used in attempting
  to achieve convergence when fitting models to the simulated data sets.
  Note that if convergence is not achieved, the simulated data set being
  used is discarded (i.e. it ``doesn't count'') and a replacement data set
  is simulated. }
  \item{verbose}{Logical scalar.  Should a \dQuote{progress report}
  be printed out at each step of the fitting procedure?}
}
\section{Remarks}{
  Although this documentation refers to \dQuote{extended
  generalised linear models}, the only such models currently
  (\today) available are the Gaussian model with the identity link,
  the Poisson model, with the log link, the Binomial model with
  the logit link, the Dbd (discretised beta distribution model),
  and the Multinom model.  The latter two are generalised linear
  models only in the \dQuote{extended} sense.  Other models may be
  added at a future date.

  When \code{eglhmm()} is called by \code{bcov()} the argument
  \code{checkDecrLL} is set equal to \code{FALSE}.  This has an
  effect only when the \code{method} used in fitting the models
  is \code{"em"}.  In this case a decrease in the log likelihood
  is treated as meaning that the algorithm has converged.   Setting
  \code{checkDecrLL} equal to \code{FALSE} is done so as to decrease
  the number of discarded data sets and thereby speed up the rate
  at which the iterations proceed.
}
\value{
  A list with components:
  \item{C_hat}{The parametric bootstrap estimate of the covariance
  matrix of the parameter estimates.}
  \item{nc.count}{A count of the total number of times that the
  algorithm failed to converge during the bootstrapping procedure.}
  \item{an.count}{A count of the \dQuote{anomalies} that occurred,
  i.e. the number of times that there was a decrease in the log
  likelihood.  Present only if the \code{method} used in fitting
  the models is \code{"em"}.
}
}
\references{
  See the help for \code{\link{eglhmm}()} for references.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\seealso{ \code{\link{fitted.eglhmm}()} \code{\link{reglhmm}()}
          \code{\link{reglhmm.default}()} \code{\link{reglhmm.eglhmm}()} }
\examples{
    \dontrun{ # Takes too long.
        fitP <- eglhmm(y~locn+depth,data=SydColCount,distr="P",cells=c("locn","depth"),
                     K=2,contr="sum",verb=TRUE,itmax=300)
	cvrP <- bcov(fitP)
        fitD <- eglhmm(y~locn+depth,data=SydColCount,distr="D",cells=c("locn","depth"),
                     K=2,nbot=0,ntop=11,contr="sum",verb=TRUE)
	cvrD <- bcov(fitD)
        fitM <- eglhmm(y~locn+depth,data=SydColDisc,distr="M",cells=c("locn","depth"),
                     K=2,contr="sum",verb=TRUE)
	cvrM <- bcov(fitM)
    }
}
\keyword{ models }
