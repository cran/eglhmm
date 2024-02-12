\encoding{UTF-8}
\name{postHocGradHess}
\alias{postHocGradHess}
\title{
   Obtain gradient and Hessian, post hoc.
}
\description{
   Calculates the gradient and Hessian of the log likelihood of
   an extended generalised hidden Markov model, from the components
   of a \code{"eglhmm"} object, or in certain circumstances, simply
   extracts these quantities from that object.
}
\usage{
postHocGradHess(object,inclTau=TRUE)
}
\arguments{
  \item{object}{
  An object of class \code{"eglhmm"} as returned by \code{\link{eglhmm}()}.
}
  \item{inclTau}{
  Logical scalar; should the vector of \dQuote{\code{tau}} parameters
  be included in the parameters under consideration?
}
}
\details{
  If \code{object} is the result of fitting a bivariate model (i.e.
  if it inherits from \code{"eglhmm.bivariate"} then an error
  is thrown.  (No gradient or Hessian is available in this case.)

  If \code{object$method} is \code{"lm"} and if \code{inclTau}
  matches the value used in fitting \code{method}, then the
  appropriate gradient and Hessian have already been calculated
  and are simply extracted from \code{object}.  If \code{inclTau}
  does not match the value used in fitting \code{method}, then
  the gradient and Hessian are recalculated (by the undocumented
  function \code{getHgl()}) with the value of \code{inclTau} being
  that specified by the function argument.

  If \code{object$method} is \code{"em"} or \code{"bf"}, then the
  gradient and Hessian are calculated (by the undocumented function
  \code{getHgl()}) with the value of \code{inclTau} being that
  specified by the function argument.

  If \code{object$method} is \code{"bf"} then the gradient has
  been calculated numerically and the Hessian \emph{may} have
  been calculated numerically (f the argument \code{hessian}
  of \code{eglhmm()} was set equal to \code{TRUE}).  The
  corresponding value of the gradient will comprise a component,
  named \code{"numGrad"}, of the list returned by this function.
  The corresponding value of the Hessian, if this was indeed
  calculated, will comprise a component, named \code{"numHess"},
  of the list returned by this function.
}
\value{
  A list with components
  \itemize{
  \item \code{gradient}: The gradient of the log likelihood.
  \item \code{Hessian}: The Hessian of the log likelihood.
  \item \code{numGrad}: The numerically calculated gradient of the log likelihood.
  Present only if \code{object$method} is \code{"bf"}.
  \item \code{numHess}: The numerically calculated Hessian of the log likelihood.
  Present only if \code{object$method} is \code{"bf"} and if argument
  \code{hessian} was set equal to \code{TRUE} in the call to \code{eglhmm()}
  that produced \code{object}.
  }
}

\references{
R. Nazim Khan, (2002). \emph{Statistical modelling and analysis of ion
channel data based on hidden Markov models and the EM algorithm}.
Ph.D. thesis, the University of Western Australia, Crawley, WA 6009.

David Oakes, Direct calculation of the information matrix via the EM
algorithm (1999).  \emph{Journal of the Royal Statistical Society},
series B, \bold{61}, pp. 479 -- 482.

Theodore C. Lystig and James P. Hughes (2002). Exact computation
of the observed information matrix for hidden Markov models.
\emph{Journal of Computational and Graphical Statistics}, \bold{11}
(3), pp. 678 -- 689.

Olivier Cappé and Eric Moulines (July 2005).  Recursive computation
of the score and observed information matrix in hidden Markov
models, \emph{IEEE Workshop on Statistical Signal Processing},
Bordeaux.
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}

\seealso{
\code{\link{eglhmm}()}
}
\examples{
fit.em <- eglhmm(y~locn+depth,data=SydColCount,distr="P",
                 cells=c("locn","depth"),K=2,method="em",verb=TRUE)
gh.em  <- postHocGradHess(fit.em) # Calculates using inclTau=TRUE.
\dontrun{
    gh.em.noTau <- postHocGradHess(fit.em,inclTau=FALSE)
    fit.lm <- eglhmm(y~locn+depth,data=SydColCount,distr="P",
                     cells=c("locn","depth"),K=2,verb=TRUE)
    gh.lm  <- postHocGradHess(fit.lm) # Just extracts the relevant components.
    gh.lm.noTau  <- postHocGradHess(fit.lm,inclTau=FALSE)
    fit.bf <- eglhmm(y~locn+depth,data=SydColCount,distr="P",
                     cells=c("locn","depth"),K=2,method="bf",verb=TRUE,
                     hessian=TRUE)
    gh.bf  <- postHocGradHess(fit.bf) # Calculates using inclTau=TRUE; also
                                      # extracts numerically computed quantities.
    gh.bf.noTau <- postHocGradHess(fit.bf,inclTau=FALSE) # Calculates; also
                                                         # extracts numerically
                                                         # computed quantities.
}
}
\keyword{ utilities }
