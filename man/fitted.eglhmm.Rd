\name{fitted.eglhmm}
\alias{fitted.eglhmm}
\title{ Predict method for extended generalised linear hidden Markov models. }
\description{
  Predicted values based on an extended generalised linear hidden Markov model
  object. }
\usage{
\method{fitted}{eglhmm}(object, ...)
}
\arguments{
  \item{object}{ An object of class \code{eglhmm} as returned
  by \code{\link{eglhmm}()}. }
  \item{\dots}{ Not used. }
}
\section{Remark}{
  Although this documentation refers to \dQuote{generalised linear
  models}, the only such models currently (\today) available are
  the Gaussian model with the identity link, the Poisson model,
  with the log link, and the Binomial model with the logit link.
  Other models may be added at a future date.
}
\value{
  A vector of fitted values of the same length as that of the
  observed values (i.e. length equal to the row dimension of the
  data frame to which the model was fitted.  This data frame is
  equal to \code{object$data} but with repeated rows corresponding to
  different states collapsed to a single row.  The row dimension of
  this data frame is thus \code{nrow(object$data)/K} where \code{K}
  is the number of states in the model.  This data frame, with
  columns \code{cf} and \code{state} omitted, is returned as an
  attribute \code{data} of the vector of fitted values.
}
\references{
  See the help for \code{\link{eglhmm}()} for references.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\seealso{ \code{\link{reglhmm}()} \code{\link{reglhmm.default}()}
          \code{\link{reglhmm.eglhmm}()} \code{\link{bcov}()} }
\examples{
    loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
    SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,] 
    SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
    rownames(SCC4) <- 1:nrow(SCC4)
    fit <- eglhmm(y~locn+depth,data=SCC4,cells=c("locn","depth"),
                 K=2,distr="P",contr="sum",verb=TRUE)
    fv  <- fitted(fit)
    with(attr(fv,"data"),plot(y[locn=="BondiOff" & depth=="40"],
             xlab="time",ylab="count"))
    with(attr(fv,"data"),lines(fv[locn=="BondiOff" & depth=="40"]))
}
\keyword{ models }
