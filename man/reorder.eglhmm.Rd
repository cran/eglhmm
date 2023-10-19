\name{reorder.eglhmm}
\alias{reorder.eglhmm}
\title{
    Reorder the states of a fitted eglhmm model
}
\description{
    Reorder the states of a fitted eglhmm model so that the
    state effects are in decreasing order.
}
\usage{
\method{reorder}{eglhmm}(x, ...)
}
\arguments{
  \item{x}{
  An object of class \code{"eglhmm"} as returned by the function
  \code{\link{eglhmm}()}.
}
  \item{\dots}{
  Not used.
}
}
\details{
   The states of a fitted hidden Markov model are usually in a
   rather arbitrary order, which can sometimes make it difficult
   to compare different fits.  This function reorders the states
   so that the state corresponding to the \dQuote{largest state
   effect} comes first (and so on down the line). What is meant by
   \dQuote{largest state effect} depends on whether the distribution
   used in the model is \dQuote{Dbd}.  If the distribution is
   \emph{not} \dQuote{Dbd}, then what is meant is simply the largest
   of those entries of \code{phi} which correspond to \code{state}.
   (The vector \code{phi} is the vector of coefficients of the
   linear predictor in the model.  Note that, since the formula
   for the model is constructed as \code{y~0+state+\dots}, the
   \dQuote{state} coefficients are unconstrained and there are as
   many of them as there are states.)

   If the distribution in question \emph{is} \dQuote{Dbd} then things
   are a bit more complicated.  We calculate the theoretical expected
   values for \dQuote{Dbd}s with parameters \eqn{\alpha =}{alpha=}
   \code{alpha[k]} and \eqn{\beta =}{beta=} \code{beta[k]} where
   \code{alpha[k]} and \code{beta[k]} are the parameter values
   corresponding to the \code{k}th state.   The states are then
   ordered according to the decreasing order of these expected
   values.  These expected values are the expected values of the
   emissions given that all predictors other than the \emph{state}
   predictors are zero.
}
\value{
    An object of class \code{c("eglhmm","reordered")} which
    is identical to the argument \code{x} in most respects.
    The components which (may) differ are:
    \itemize{
    \item \code{tpm}
    \item \code{ispd}
    \item \code{phi}
    \item \code{theta}
    \item \code{Hessian}
    \item \code{gradient}
    \item \code{mean} and \code{sd}, or \code{lambda} or \code{p},
    or \code{alpha} and \code{beta} (depending on which distribution
    is being used)
    \item \code{fy}
    }
    The entries of these components will have the same numerical
    values as before but, given that the ordering of the states has
    actually changed, will have different orderings, corresponding
    to the new ordering of the states.

    Note that the \emph{attribute} \code{preSpecSigma} of the
    component \code{theta}, may differ from what it was in \code{x}.

    The returned value will also have a component \code{neworder}
    which is an integer vector providing the indices of the
    reordering of the states.  It also currently (\today) has a
    component \code{newlog.like}.  This should (if there is any
    justice in the world --- but there isn't!) have the same value
    as the component \code{log.like}.  Once I am confident that
    everything is working as it should, the \code{newlog.like}
    component will be removed.
} \author{Rolf Turner
  \email{rolfturner@posteo.net}
}

\seealso{
    \code{\link{eglhmm}()}
}
\examples{
    loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
    SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,]
    SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
    rownames(SCC4) <- 1:nrow(SCC4)
    fit  <- eglhmm(y~locn+depth,data=SCC4,cells=c("locn","depth"),
                  K=2,distr="P",verb=TRUE)
    ofit <- reorder(fit)
}
