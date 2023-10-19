\name{reglhmm}
\alias{reglhmm}
\alias{reglhmm.default}
\alias{reglhmm.eglhmm}
\title{ Simulate data from a hidden generalised linear Markov
        model.}
\description{
  Takes a specification of the model and simulates the data from
  that model.  The model may be specified in terms of the individual
  components of that model (the default method).  The components
  include a data frame that provides the predictor variables,
  and various parameters of the model.  For the \code{"eglhmm"}
  method the model is specified as a fitted model, an object of
  class \code{"eglhmm"}.
}
\usage{
reglhmm(x,...)
\method{reglhmm}{default}(x, formula, response, cells=NULL, data=NULL, nobs=NULL,
                         distr=c("Gaussian","Poisson","Binomial","Dbd","Multinom"),
                         phi, Rho, sigma, size, ispd=NULL, ntop=NULL, zeta=NULL,
                         missFrac = 0, fep=NULL,
                         contrast=c("treatment","sum","helmert"),...)
\method{reglhmm}{eglhmm}(x, missFrac = NULL, ...)
}
\arguments{
  \item{x}{
  For the default method, the transition probability matrix of
  the hidden Markov chain.  For the \code{"eglhmm"} method,
  an object of class \code{"eglhmm"} as returned by the function
  \code{\link{eglhmm}()}.
}
  \item{formula}{
  The formula specifying the generalised linear model from which data
  are to be simulated.  Note that the predictor variables in
  this formula must include a factor \code{state}, which specifies
  the state of the hidden Markov chain.  Note also that this formula
  must determine a design matrix having a number of columns equal to
  the length of the vector \code{phi} of model coefficients provided
  in \code{object} (and to the length of \code{psi} in the case of
  the Gaussian distribution).  If this condition is not satisfied,
  an error is thrown.

  It is advisable to use a formula specified in the manner
  \code{y~0+state+\dots} where \code{\dots} represents the predictors
  in the model other than \code{state}.  Of course \code{phi} must
  be supplied in a manner that is consistent with this structure.
}
  \item{response}{A character vector of length 2, specifying
  the names of the responses.  Ignored unless \code{distr} is
  \code{"Multinom"}.  If \code{distr} is \code{"Multinom"} and if
  \code{response} is provided appropriately, then the simulated
  data are bivariate multinomial.
}
  \item{cells}{
  A character vector specifying the names of the factors which
  determine the ``cells'' of the model.  These factors must be
  columns of the data frame \code{data}.  (See below.)  Each cell
  corresponds to a time series of (simulated) observations.
  If \code{cells} is not supplied (left equal to \code{NULL})
  then the model is taken to have a single cell, i.e. data from a
  \dQuote{simple} hidden Markov model is generated.   The parameters
  of that model may be time-varying, and still depend on the
  predictors specified by \code{formula}.
}
  \item{data}{
  A data frame containing the predictor variables referred to by
  \code{formula}, i.e. the predictors for the model from which
  data are to be simulated.  If \code{data} is not specified,
  the \code{nobs} (see below) must be.  If \code{data} is not
  specified then \code{formula} must have the structure \code{y ~
  state} or preferably \code{y ~ 0 + state}.  Of course \code{phi}
  must be specified in a consistent manner.
}
  \item{nobs}{
  Integer scalar.  The number of observations to be generated in
  the setting in which the generalised linear model in question is
  vacuous.  Ignored if \code{data} is supplied.
}

  \item{distr}{
  Character string specifying the distribution of the
  \dQuote{emissions} from the model, i.e., of the observations.
  This distribution determines \dQuote{emission probabilities}.
}
  \item{phi}{
  A numeric vector specifying  the coefficients of the linear
  predictor of the generalised linear model.  The length of
  \code{phi} must be equal to the number of columns of the
  design matrix determined by \code{formula} and \code{data}.
  The entries of \code{phi} must match up appropriately with
  the columns of the design matrix.
}
  \item{Rho}{A matrix, or a list of two matrices or a three dimensional
  array specifying the emissions probabilities for a multinomial
  distribution.  Ignored unless \code{distr} is \code{"Multinomial"}.
}
  \item{sigma}{
  A numeric vector of length equal to the number of states.
  Its \eqn{i}th entry is the standard deviation of the (Gaussian)
  distribution corresponding to the \eqn{i}th state.  Ignored unless
  \code{distr} is \code{"Gaussian"}.
}
  \item{size}{
  Integer scalar.  The number of trials (sample size) from which
  the number of \dQuote{successes} are counted, in the context of
  the binomial distribution.  (I.e. the \code{size} parameter of
  \code{rbinom()}.)  Ignored unless \code{distr} is \code{"Binomial"}.
}

  \item{ispd}{
  An optional numeric vector specifying the initial state probability
  distribution of the model. If \code{ispd} is not provided then it
  is taken to be the stationary/steady state distribution determined
  by the transition probability matrix \code{x}.  If specified,
  \code{ispd} must be a \emph{probability} vector of length equal
  to the number of rows (equivalently the number of columns)
  of \code{x}.
}

 \item{ntop}{
  Integer scalar, strictly greater than 1.  The maximum possible
  value of the db distribution.  See \code{\link[dbd]{db}()}.
  Used only if \code{distr} is \code{"Dbd"}.
}
  \item{zeta}{
  Logical scalar.  Should zero origin indexing be used?
  I.e. should the range of values of the db distribution be taken to
  be \code{{0,1,2,...,ntop}} rather than \code{{1,2,...,ntop}}?
  Used only if \code{distr} is \code{"Dbd"}.
}

  \item{missFrac}{
  A non-negative scalar, less than 1.  Data will be randomly set
  equal to \code{NA} with probability \code{miss.frac}.  Note that
  for the \code{"eglhmm"} method, if \code{"miss.frac"} is not
  supplied then it is extracted from \code{object}
}

  \item{fep}{A list of length 1 or 2.  The first entry of this
  list is a logical scalar.  If this is \code{TRUE}, then the first
  entry of the simulated emissions (or at least one entry of the first
  pair of simulated emissions) is forced to be \dQuote{present},
  i.e. non-missing.  The second entry of \code{fep}, if present, is
  a numeric scalar, between 0 and 1 (i.e. a probability).  It
  is equal to the probability that both entries of the first
  pair of emissions are present.  It is ignored if the emissions
  are univariate.  If the emissions are bivariate but the second
  entry of \code{fep} is not provided, then this second entry defaults
  to the \dQuote{overall} probability that both entries of a pair of
  emission are present, given that at least on is present.
  This probability is calculated from \code{nafrac}.
}

  \item{contrast}{
  A character string, one of ``treatment'', ``helmert'' or ``sum'',
  specifying what contrast (for unordered factors) to use in
  constructing the design matrix.  (The contrast for ordered factors,
  which is has no relevance in this context, is left at it default
  value of \code{"contr.poly"}.)  Note that the meaning of the
  coefficient vector \code{phi} depends on the contrast specified,
  so make sure that the contrast is the same as what you had in
  mind when you specified phi!!!  Note that for the \code{"eglhmm"}
  method, \code{contrast} is extracted from \code{x}.
}
  \item{\dots}{
  Not used.
}
}
\section{Remark}{
  Although this documentation refers to \dQuote{generalised linear
  models}, the only such models currently (\today) available are
  the Gaussian model with the identity link, the Poisson model,
  with the log link, and the Binomial model with the logit link.
  The Multinomial model, which is also available, is not exactly
  a generalised linear model; it might be thought of as an
  \dQuote{extended} generalised linear model.  Other models may be
  added at a future date.
}
\value{
  A data frame with the same columns as those of \code{data}
  and an added column, whose name is determined from \code{formula},
  containing the simulated \emph{response}
}
\references{
  T. Rolf Turner, Murray A. Cameron, and Peter J. Thomson
  (1998).  Hidden Markov chains in generalized linear models.
  \emph{Canadian Journal of Statististics} \bold{26}, pp. 107 --
  125, DOI: https://doi.org/10.2307/3315677.

  Rolf Turner (2008).  Direct maximization of the likelihood of a
  hidden Markov model. \emph{Computational Statistics
  and Data Analysis} \bold{52}, pp. 4147 -- 4160, DOI:
  https://doi.org/10.1016/j.csda.2008.01.029
}
\author{Rolf Turner
  \email{rolfturner@posteo.net}
}
\seealso{
    \code{\link{fitted.eglhmm}()}
    \code{\link{bcov}()}
}
\examples{
    loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
    SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,] 
    SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
    rownames(SCC4) <- 1:nrow(SCC4)
    Tpm   <- matrix(c(0.91,0.09,0.36,0.64),byrow=TRUE,ncol=2)
    Phi   <- c(0,log(5),-0.34,0.03,-0.32,0.14,-0.05,-0.14)
    # The "state effects" are 1 and 5.
    Dat   <- SCC4[,1:3]
    fmla  <- y~0+state+locn+depth
    cells <- c("locn","depth")
# The default method.
    X     <- reglhmm(Tpm,formula=fmla,cells=cells,data=Dat,distr="P",phi=Phi,
                    miss.frac=0.75,contrast="sum")
# The "eglhmm" method.
    fit <- eglhmm(y~locn+depth,data=SCC4,cells=cells,K=2,
                 verb=TRUE,distr="P")
    Y   <- reglhmm(fit)
# Vacuous generalised linear model.
    Z   <- reglhmm(Tpm,formula=y~0+state,nobs=300,distr="P",phi=log(c(2,7)))
    # The "state effects" are 2 and 7.
}
\keyword{ models }
\keyword{ datagen }
