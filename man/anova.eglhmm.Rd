\name{anova.eglhmm}
\alias{anova.eglhmm}
\title{
Test for a difference between two fitted \code{eglhmm} models.
}
\description{
Apply a likelihood ratio test to determine whether the difference,
between the log likelihood statistics of two fitted \code{eglhmm} models,
is statistically significant.
}
\usage{
\method{anova}{eglhmm}(object, ...)
}
\arguments{
  \item{object}{
  An object of class \code{"eglhmm"} as returned by \code{\link{eglhmm}()}.
}
  \item{\dots}{
  Precisely \emph{one} more object of class \code{"eglhmm"}, to be
  compared with \code{object}.
}
}
\details{
  This anova method handles only comparisons between \emph{two}
  models.)  The order of the arguments (i.e. which object is passed
  as \dQuote{\code{object}} and which is passed as the sole entry
  of the \dots argument) is immaterial.
}
\value{
A list with components
\itemize{
    \item \code{stat} the likelihood ratio statistic, i.e. the difference
    between the log likelihoods of the two models.  That for the model
    with the smaller number of parameters is subtracted from that for
    the model with the larger number.
    \item \code{df} the degrees of freedom of the likelihood ratio statistic,
    i.e. the difference between the number of parameters of the respective
    models.  The smaller number is subtracted from the larger.
    \item \code{pvalue} the \eqn{p}-value of the test as given by
    \code{pchisq(stat, df, lower.tail = FALSE)}.
}
This list has an attribute \code{"details"} consisting of a numeric
vector of length four with entries \code{ll1} (the smaller of the
log likelihoods), \code{ll2} (the larger of the log likelihoods),
\code{np1} (the smaller of the parameter counts) and \code{np2}
(the larger of the parameter counts).
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\examples{
fit1 <- eglhmm(formula=y~locn+depth,data=SydColCount,
               cells=c("locn","depth"),distr="P",K=2,
               method="em",verb=TRUE)
fit2 <- eglhmm(formula=y~locn+depth+ma.com+nh.com+bo.com,data=SydColCount,
               cells=c("locn","depth"),distr="P",K=2,
               method="em",verb=TRUE)
anova(fit1,fit2)
}

\keyword{ htest }
