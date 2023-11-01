\name{eglhmm}
\alias{eglhmm}
\title{ Fit (extended) generalised linear hidden Markov models.}
\description{
  Fits an (extended) generalised linear model to a data set where the
  response in each \dQuote{cell} of the model consists of a
  time series whose serial dependence is modelled by a
  hidden Markov model.
}
\usage{
eglhmm(formula = NULL, response = NULL, data,
      distr = c("Gaussian", "Poisson", "Binomial", "Dbd", "Multinom", "discnp"),
      inclTau=TRUE,preSpecSigma=NULL, indep = NULL, size = NULL, nbot = NULL, ntop = NULL,
      cells = NULL, K = NULL, par0 = NULL, randStart = NULL,
      method = c("lm", "em", "bf"), optimiser = c("optim", "nlm"),
      optimMethod = "BFGS", nlmWarn = FALSE, lmc = 10, tolerance = NULL,
      digits = NULL, verbose = FALSE, itmax = 200,
      contrast = c("treatment", "sum", "helmert"),
      crit = c("CLL", "L2", "Linf", "ABSGRD"), breaks = NULL, hessian = FALSE,
      useAnalGrad = FALSE, ca = FALSE, checkDecrLL=TRUE)
}
\arguments{
  \item{formula}{
  A model formula specifying the linear predictor for the model. The
  formula should \emph{not} include \code{state} as a predictor
  variable.  The variable \code{state} gets added to the formula
  automatically.  Ignored if the model is bivariate, i.e. if the
  length of \code{response} is 2.
}
  \item{response}{
  A character scalar or a length-2 vector of such scalars,
  specifying the name or names of the response(s).
  If \code{response} is not specified (i.e. if it is left as
  \code{NULL}) then \code{formula} (see below) \emph{must} be
  specfied and \code{response} is taken to be the left hand side
  of \code{formula}.  (In this case, it is of course univariate.)
}
  \item{data}{
  A data frame with columns providing the response(s) and the predictor
  variables in the model.
}
\item{distr}{
  Character string specifying the distribution of the response(s)
  (\dQuote{emissions}) variable(s).  Currently (\today) the only
  distributions accommodated are Gaussian, Poisson, Binomial, Dbd,
  and Multinom.  Note that \code{"discnp"} is just an alternative
  expression for \code{"Multinom"}.  Ignored if the response is
  bivariate, in which case \code{distr} is forcibly set equal
  to \code{"Multinom"}.  I.e.  bivariate models are, currently,
  fitted only to data in which the emissions have the \code{"Multinom"}
  distribution.
}
\item{inclTau}{
  Logical scalar. Should the transition probability matrix parameters
  \dQuote{\code{tau}} be included in those that are estimated via the
  Hessian/gradient pardigm? In this case, they are included in the
  set of parameters to which the gradient and Hessian are applicable.
  If not, they are estimated via the method of moments as is done
  when the EM algorithm is used.  In this latter case the dimensions
  of the Hessian are reduced (by a substantial amount if \code{K}
  is \dQuote{large}).
}
\item{preSpecSigma}{
  Numeric vector of length \code{K} (see below) with strictly
  positive entries.  Ignored if \code{distr} is not equal to
  \code{"Gaussian"}. This vector provides \dQuote{pre-specified}
  values of the standard deviations \code{sigma} of the Gaussian
  distribution associated with each state.  If \code{preSpecSigma}
  is specified, then it is used as the value of \code{sigma}
  throughout the fitting process, and \code{sigma} is \emph{not}
  estimated from the data. If \code{distr} is \code{"Gaussian"}
  and \code{preSpecSigma} is specified, then an error will be thrown
  if the length of \code{preSpecSigma} is not equal to \code{K}, or
  if any entries of \code{preSpecSigma} fail to be strictly positive.
}

\item{indep}{
   Logical scalar; should the components of a bivariate
   model be considered to be independent?  Ignored unless the
   model is bivariate (i.e. unless \code{response} is of length 2.
   If the model is bivariate and \code{indep} is not specified,
   an error is thrown.
}
  \item{size}{
  Scalar integer specifying the number of trials in the experiment generating
  binomial data (see the \code{size} argument of \code{\link{dbinom}()}).
  Ignored unless \code{distr} is equal to \code{"Binomial"}.
}
  \item{nbot}{
  Scalar integer specifying the lower end (0 or 1) of the range
  of values of the discretised Beta distribution.  Ignored unless
  \code{distr} is \code{"Dbd"}.
}
  \item{ntop}{
  Scalar integer specifying the upper end of the range of values
  of the discretised Beta distribution.  Ignored unless \code{distr}
  is \code{"Dbd"}.
}
  \item{cells}{A character vector giving the names of the factors
  (columns of the \code{data} data frame) which determine what the
  \dQuote{cells} of the model are considered to be.  The cells correspond
  to the combinations of levels of the factors named by \code{cells}.
}
  \item{K}{
  Scalar integer specifying the number of states of the hidden Markov
  model in question.  If \code{K} is not specified and \code{par0}
  (see below) \emph{is} specified, and has a component \code{tpm},
  then \code{K} is set equal to \code{nrow(tpm)}.  If \code{par0}
  does not have a \code{tpm} component, an error is thrown.
  An error is also thrown in this setting if \code{K} is specified
  to a value different from \code{nrow(tpm)}.
}
  \item{par0}{
  A list comprising starting values parameter estimates, to be used
  by the various methods.  (See \code{method} below.)  This list
  may have components \code{tpm} (an estimate of the transition
  probability matrix), \code{phi} (a vector of estimates of
  the coefficients in the linear predictor in the generalised
  linear model) and \code{Rho} (a matrix, a list of two matrices,
  or a three dimensional array) that specifyies the emission probabilities
  when \code{distr} is \code{"Multinomial"}.
}
  \item{randStart}{
  Either a logical scalar or a list of three logical scalars
  named \code{tpm}, \code{phi}, and \code{Rho}.  If the former,
  it is converted internally into a list with entries named
  \code{tpm}, \code{phi} and \code{Rho}, all having the same
  value as the original argument.  If \code{tpm} is TRUE then the
  (undocumented) function \code{inititialise()} chooses entries
  for the starting value of \code{tpm} at random; likewise for
  \code{phi} and \code{Rho}.  If left \code{NULL}, this argument
  defaults to \code{list(tpm=FALSE,phi=FALSE,Rho=FALSE)}.
}
  \item{method}{
  Character string specifying the method used to fit the model. This
  may be \code{"lm"} (Levenberg-Marquardt algorithm), \code{"em"}
  (EM algorithm) or \code{"bf"} (\dQuote{brute force}).  The latter
  calls upon \code{\link{optim}()} or \code{\link{nlm}()} to do the
  heavy lifting).  If the response is bivariate, then \code{method}
  is forcibly (and silently) set equal to \code{"em"}.
}
  \item{optimiser}{
  Character string specifying which of \code{optim()} or \code{nlm()}
  should be used when \code{method} is \code{"bf"}.  Ignored unless
  \code{method} is \code{"bf"}.
}
  \item{optimMethod}{
  Character string specifying the optimisation method to be used
  by \code{optim()}.  See \code{\link{optim}()} for details.
  Ignored unless \code{method} is \code{"bf"} and \code{optimiser}
  is \code{"optim"}.
}
  \item{nlmWarn}{
  The \code{nlm()} function sometimes produces, in the first few
  iterations,  warnings to the effect \dQuote{NA/Inf
  replaced by maximum positive value}.  These warnings are almost
  surely irrelevant and are annoying.  If \code{nlmWarn} is \code{FALSE}
  (the default) then these warnings are suppressed.  This argument
  is provided to allow for the remote possibilty that the user might
  want to see these warnings.
}
\item{lmc}{
  Positive numeric scalar.  The initial \dQuote{Levenberg-Marquardt
  constant}.  Ignored unless \code{method} is \code{"lm"}.
}
  \item{tolerance}{
  Positive numeric scalar.  The convergence tolerance to be used.
  What this value actually \emph{means} depends upon \code{method}.
  If left as \code{NULL} it defaults to \code{1e-6} for the bivariate
  methods, to \code{sqrt(.Machine$double.eps} for the \code{"em"}
  and \code{"lm"} methods, and to the default value of \code{reltol}
  used by \code{optim()} when \code{method} is \code{"bf"} and
  \code{optimiser} is \code{"optim"}.  It is ignored if \code{method}
  is \code{"bf"} and \code{optimiser} is \code{"nlm"}.
}
  \item{digits}{
  Integer scalar. The number of digits to which \dQuote{progress
  reports} are printed when \code{verbose} (see below) is
  \code{TRUE}.  There is a \dQuote{sensible} default which is
  calculated in terms of \code{tolerance}.  This argument is ignored
  if \code{method} is \code{"bf"}.
}
  \item{verbose}{
  Logical scalar; if \code{TRUE}, rudimentary \dQuote{progress
  reports} are printed out at appropriate points during the
  iteration process.  The nature of these \dQuote{reports} varies
  with \code{method}.
}
  \item{itmax}{
  Integer scalar.  The maximum number of iterative steps to
  take.  Has a somewhat different meaning when \code{method} is
  \code{"bf"}, in which case the meaning depends on \code{optimiser}.
  For methods \code{"em"} and \code{"lm"}, if convergence is
  not achieved by \code{itmax} steps, the function gives up,
  prints a message to this effect, and returns a value with a
  component \code{converged=FALSE}.  This returned value may be
  used as a starting (the value of the argument \code{par0}) so
  that the iterations may be continued from where they left off.
  Unfortunately this facility is not available when \code{method}
  is \code{"bf"}.
}
  \item{contrast}{
  Text string specifying the contrast (in respect of
  unordered factors) (see \code{\link{contrasts}()} and
  \code{\link{options}()}) that will be used when the design matrix
  is constructed from the model formula.  May be abbreviated
  (e.g. to "t", "s" or "h").
}
  \item{crit}{
  Text string specifying the stopping criterion to be used.
  Possible values are \dQuote{CLL} (scaled change in log likelihood),
  \dQuote{L2} (scaled square root of the sum of squares of the
  changes in the parameter estimates), \dQuote{Linf} (scaled maximum
  of the absolute value of the changes in the parameter estimates),
  and \dQuote{ABSGRD} (scaled maximum of the absolute values of
  the entries of the gradient vector).  The latter only makes sense
  for the Levenberg-Marquardt algorithm.

  This argument is ignored if \code{method} is \code{"bf"}.
  It seems that the \code{"bf"} method effectively uses \dQuote{CLL}
  when \code{optimiser} is \code{"optim"}.  When \code{optimiser}
  is \code{"nlm"} it seems that a combination of (something like)
  \dQuote{ABSGRD} and \dQuote{CLL} is used.
}
  \item{breaks}{
  A vector of \code{K+1} values used to construct a set of guesses
  at the states corresponding to each observation.  These are
  in turn used to calculate an initial estimate of the transition
  probability matrix.  There is a \dQuote{sensible} default (produced
  by the undocumented function \code{breaker()}.
}
  \item{hessian}{
  Logical scalar; should a Hessian matrix obtained by numerical
  differentiation be returned?  Ignored unless \code{method}
  is \code{"bf"}.
}
  \item{useAnalGrad}{
  Logical scalar; should \dQuote{analytical} calculation of the
  gradient be conducted?  This argument is ignored unless the method
  is \code{"bf"}.
}
  \item{ca}{
  Logical scalar; \dQuote{check analyticals}.  Used only when
  the method is \code{"bf"} and \code{optimiser} is \code{"nlm"},
  and is passed on to \code{\link{nlm}()}.
}
  \item{checkDecrLL}{
  Logical scalar; \dQuote{check for a \emph{decrease} in the log likelihood}.
  Ignored unless the \code{method} is \code{"em"}.  Should the software
  check for a decrease in the log likelihood after an EM step?  See
  the \bold{Remarks} for further discussion.
}
}
\section{Remarks}{
  \describe{
  \item{\bold{Available models:}}{Although this documentation refers to
  (extended) \dQuote{generalised linear models}, the only such
  models currently (\today) available are the Gaussian model
  with the identity link, the Poisson model, with the log link,
  and the Binomial model with the logit link.  When \code{distr}
  is \code{"Dbd"} or \code{"Multinom"} the model fitted is is
  a generalised linear model only in a rather extended sense.
  Even the Gaussian model is not strictly speaking a generalised
  linear model, since the (state dependent) standard deviations
  are estimated by a method separate from the generalised linear
  model paradigm.  Other models may be added at a future date.}

  \item{\bold{Decrease in the log likelihood:}}{
  If \code{method} is equal to \code{"EM"} there may be a
  \emph{decrease} (!!!) in the log likelihood at some EM step.
  This is \dQuote{theoretically impossible} but can occur in
  practice due to an intricacy in the way that the EM algorithm
  treats \code{ispd} when \code{stationary} is \code{TRUE}.
  It turns out to be effectively impossible to maximise the expected
  log likelihood unless the term in that quantity corresponding
  to \code{ispd} is ignored (whence it \emph{is} ignored).
  Ignoring this term is \dQuote{asymptotically negligible} but
  can have the unfortunate effect of occasionally leading to a
  decrease in the log likelihood.  If \code{method} is equal
  to \code{"em"}, then the object returned by \code{eglhmm()}
  has a component \code{anomaly} which is \code{TRUE} if such
  a decrease in the log likelihood was detected, and \code{FALSE}
  otherwise.

  If such a decrease/anomaly is detected, then (provided that
  \code{checkDecrLL} is \code{TRUE}) the algorithm terminates and
  the \code{converged} component of the returned value is set equal
  to \code{NA}.  The algorithm issues a message to the effect that
  the decrease occurred.  The message suggests that another method
  be used and that perhaps the results from the penultimate EM step
  (which are returned by this function) be used as starting values.
  This of course is not possible if the response is bivariate,
  in which case only the EM algorithm is applicable.

  Note that if \code{checkDecrLL} is \code{FALSE}, then the algorithm
  proceeds \dQuote{normally}.  That is, it treats the decrease in
  the log likelihood to mean that the \dQuote{increase} in the
  log likeihood is less than \code{tolerance} and deems convergence
  to be achieved.

  The value of \code{checkDecrLL} is set to \code{FALSE} in the
  function \code{\link{bcov}()} so as to speed up the rate at which
  the iterations proceed.  In other circumstances it is probably
  judicious to leave it at its default value of \code{TRUE}.}
  }
}
\value{
  An object of class \code{"eglhmm"}, consisting of a list with
  components:

  \item{call}{The call by which this object was created.  Present
  so that \code{update()} can be applied to objects returned by
  \code{eglhmm()}.}

  \item{tpm}{The estimated transition probability matrix.}

  \item{ispd}{The estimated initial state probability distribution.}

  \item{phi}{Except for the \code{"Multinom"} distribution this
   is the vector of estimated coefficients of the linear predictor
   in the generalised linear model.  For the \code{"Multinom"}
   distribution it consists of the entries of \code{Rho} (see below)
   with the final all-zero column remove.  In this case \code{phi}
   is of course redundant.}

  \item{theta}{The vector of parameter estimates that the estimation
   procedure actually works with. It consists of the catenation of
   the non-redundant parameterization of the transition probability
   matrix and the vector \code{phi}. It is redundant in the case of
   the \code{"Multinom"} distribution.}

  \item{Rho}{A matrix, or a list of two matrices or a three dimensional
   array specifying the emissions probabilities for a multinomial
   distribution.  Present only if \code{distr} is \code{"Multinom"}.}

  \item{log.like}{The value of the log likelihood of the model
   evaluated at the parameter estimates, i.e. the (approximately)
   maximal value of the log likelihood.}

  \item{gradient}{(Not present for the \code{"em"} method.) The gradient
   vector of the log likelihood at the final parameter estimates;
   it \emph{should} be effectively the zero vector.}

  \item{numHess}{(Present only if \code{method} is \code{"bf"}
   and only if the argument \code{hessian} is \code{TRUE}.)  A value
   of the Hessian matrix (see below), obtained by means of numerical
   differentiation.}

  \item{Hessian}{(Present only if \code{method} is \code{"lm"}.) The
   Hessian matrix, i.e. the matrix of second partial derivatives of
   the log likelihood, evaluated at the final parameter estimates.
   The inverse of the negative of this matrix constitutes an estimate
   of the covariance matrix of the parameter estimates.}

  \item{mu}{A data frame with \code{npred+1} columns where \code{npred}
   is the number of predictors in the model.  The rows contain, in their
   first \code{npred} entries, all possible combinations of the predictor
   values.  The last (\code{npred+1}) entry of each row is the fitted mean
   of the Gaussian distribution, as determined by that combination
   of predictors.  Present only if \code{distr} is \code{"Gaussian"}.}

  \item{sigma}{Numeric vector of length \code{K} whose entries
   consist of the fitted standard deviations for the underlying
   Gaussian distribution, corresponding to each of the states.
   Present only if \code{distr} is \code{"Gaussian"} and
   \code{preSpecSigma} is \emph{not} supplied.}

  \item{preSpecSigma}{Numeric vector equal to the \code{preSpecSigma}
   argument, with names \code{"sigma1"}, \code{"sigma2"}, \dots,
   \code{"sigmaK"} added.  Present only if \code{distr} is
   \code{"Gaussian"} and \code{preSpecSigma} is supplied.}

  \item{stopCritVal}{Numeric scalar equal to the value, assumed by the
   stopping criterion specified by the argument \code{crit}, at
   the termination of the algorithm.  If the algorithm converged
   then \code{stopCritVal} will be less than \code{tolerance}.
   Not present if \code{method} is \code{"bf"}.  If \code{converged}
   (see below) is \code{NA} then \code{stopCritVal} is \code{NA}
   also.}

   \item{anomaly}{Logical scalar.  Did an \dQuote{anomaly} occur
   in an application of the EM algorithm?  (See \bold{Remarks}.) `
   Present only if \code{method} was equal to \code{"em"}.  This
   entry of the returned value is provided mainly for use by the
   \code{bcov()} function. Note that \code{anomaly} is added to the
   returned object, irrespective of the value of \code{checkDecrLL}.
   When \code{checkDecrLL} is \code{TRUE}, \code{anomaly} is
   somewhat redundant, since it will be \code{TRUE} if aand only if
   \code{converged} is \code{NA}.  However when \code{checkDecrLL}
   is \code{FALSE}, \code{anomaly} is informtive, since it is not
   possible to tell from other entries of the returned value when
   an anomaly has occurred.}

  \item{converged}{A logical scalar.  For the \code{"lm"}, and
   \code{"em"} methods it is \code{TRUE} if convergence is achieved
   within \code{itmax} iterations and \code{FALSE} otherwise.
   For the \code{"em"} method, if \code{checkDecrLL} is \code{TRUE},
   then \code{converged} may be \code{NA}.  See \bold{Remarks}
   for some discussion.

   For the \code{"bf"} method \code{converged} is \code{TRUE}
   if the \code{convergence} component of the object returned by
   \code{\link{optim}()} is equal to 0 or if the \code{code}
   component of the object returned by \code{\link{nlm}()} is
   less than or equal to 2, and is \code{FALSE} otherwise.  When
   \code{nlm()} is used, the value of \code{converged} has an attribute
   \code{"code"} equal to the actual value of the \code{code}
   component.}

  \item{nstep}{The number of steps (iterations) actually used by
   the algorithm.  For the \code{"lm"} and \code{"em"} methods
   this is the number of Levenberg-Marquardt steps, or EM steps,
   respectively, taken by the algorithm.  For the \code{"bf"}
   method it is the \code{counts} component of the object returned
   by \code{\link{optim}()} when \code{optimiser} is \code{"optim"}
   and it is the \code{iterations} component of the object returned by
   \code{\link{nlm}()} when \code{optimiser} is \code{"nlm"}.}

  \item{mean}{A vector of the fitted mean values underlying each
   combination of observed predictors and state (i.e. corresponding to each
   entry of \code{y} in the data frame used to fit the model.
   See the description of \code{data} below.  Present only if
   \code{distr} is \code{"Gaussian"}.}

  \item{sd}{A vector of the fitted values of the standard
   deviations underlying each combination of observed predictors and state,
   i.e. corresponding to each entry of \code{y} in the data frame
   used to fit the model.  See the description of \code{data} below.
   Present only if \code{distr} is \code{"Gaussian"}.}

  \item{lambda}{A vector of estimated values of the Poisson
   parameter associated with each combination of observed predictors and state,
   i.e. corresponding to each entry of \code{y} in the data frame
   used to fit the model.  See the description of \code{data} below.
   Present only if \code{distr} is \code{"Poisson"}.}

  \item{p}{A vector of estimated values of the \dQuote{success}
   probabilities associated with each combination of observed
   predictors and state, i.e. corresponding to each entry of
   \code{y} in the data frame used to fit the model.  See the
   description of \code{data} below.  Present only if \code{distr}
   is \code{"Binomial"}.}

  \item{alpha}{A numeric vector of the fitted \dQuote{alpha}
   parameters, of the discretised Beta distribution, corresponding to
   each observation.  Present only if \code{distr} is \code{"Dbd"}.}
   
  \item{beta}{A numeric vector of the fitted \dQuote{beta}
   parameters, of the discretised Beta distribution, corresponding to
   each observation.  Present only if \code{distr} is \code{"Dbd"}.}

  \item{fy}{The values of the \dQuote{emission probability (density)}
   function, calculated at each observed value, for each state
   (i.e. at each entry of \code{y} in \code{data}.  See below.)
   These values are calculated using the (final) fitted parameters.}

  \item{message}{A (long) text string that is produced if the EM
   algorithm encounters the anomaly of a decrease in the log
   likelihood after an EM step.  It warns the user that this has
   occurred and suggests consulting the help file for an explanation.
   Present only if \code{method=="em"}, the anomaly referred to
   has occurred, and \code{checkDecrLL} is \code{TRUE}.}

  \item{par0}{The starting values used in the estimation procedure.
   Either those provided by the argument \code{par0} or those created
   by the (undocumented) function \code{initialise}.}

  \item{cells}{A character vector indicating the names of the
   factors specifying the \dQuote{cells} of the model.  (Equal to
   the \code{cells} argument.)}

  \item{formula}{The formula for the model that was fitted; equal
   to the \code{formula} argument, augmented by \code{state}.}

  \item{distr}{Text string specifying the distribution of the
   response variable.  Equal to the \code{distr} argument of this
   function.}

  \item{nbot}{Integer scalar. The lower endpoint of the range of
   values of the discretised beta distribution.  Equal to the value
   of the \code{nbot} argument of this function.  Present only if
   \code{distr} is \code{"Dbd"}.}

  \item{ntop}{Integer scalar. The upper endpoint of the range of
   values of the discretised beta distribution.  Equal to the value
   of the \code{nbot} argument of this function.  Present only if
   \code{distr} is \code{"Dbd"}.}

  \item{size}{Scalar integer equal to the number of trials in the
   \dQuote{experiments} generating the data.  Equal to the \code{size}
   argument of this function.  Present only if \code{distr} is
   \code{"Binomial"}.}

  \item{tolerance}{The convergence tolerance used to fit the model.
   Equal to the \code{tolerance} argument.}

  \item{crit}{Character scalar specifying the stopping criterion
   that was used. Equal to the \code{crit} argument of this function.
   Not present if \code{method} is \code{"bf"}.}

  \item{contrast}{Text string specifying the contrast for unordered
   factors that was used in fitting the model.  Equal to the
   \code{contrast} argument of this function.}

  \item{method}{The method (\code{"lm"}, \code{"em"}, or \code{"bf"}
   used to fit the model.  Equal to the \code{method} argument.}

  \item{stationary}{Logical scalar.  Was a stationary Markov
   chain fitted?  Currently (\today) \code{stationary} is always
   \code{TRUE}.}

  \item{data}{The data frame to which the model was fitted.
   It is a rearrangement of the \code{data} argument, with rows of
   that argument replicated \code{K} times (once for each state).
   A \code{state} column (factor) has been added, as has a column
   \code{cf} (\dQuote{cell factor}), which indicates, by means
   of a single factor, which cell of the model a given row of
   \code{data} corresponds to, is added.  The aforementioned
   rearrangement consists of ordering the cells in the order of
   the levels of \code{cf}.  When \code{distr} is \code{"Multinom"}
   the \code{"response"} variables are coerced into factors. }

  \item{bicm}{Numerical scalar.  The number by which \code{npar}
   is multiplied to form the \code{BIC} criterion.  It is
   essentially the log of the number of observations.  See the
   code of \code{eglhmm()} for details.}

  \item{AIC}{Numerical scalar.  The Akaike Information criterion,
   calculated as \code{-2*ll + 2*npar} where \code{ll} is the log
   likelihood of the fitted model and \code{npar} is the number of
   fitted parameters.}

  \item{BIC}{Numerical scalar.  The Bayesian Information criterion,
   calculated as \code{-2*ll + bicm*npar} where \code{ll} is the log
   likelihood of the fitted model, \code{npar} is the number of
   fitted parameters, and \code{bicm} is the log of the number of
   observations.}

  \item{missFrac}{The fraction or proportion of missing values in
   the observations.}
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

\seealso{ \code{\link{fitted.eglhmm}()} \code{\link{reglhmm.default}()}
          \code{\link{reglhmm.eglhmm}()} \code{\link{bcov}()} }

\examples{
    loc4 <- c("LngRf","BondiE","BondiOff","MlbrOff")
    SCC4 <- SydColCount[SydColCount$locn \%in\% loc4,] 
    SCC4$locn <- factor(SCC4$locn) # Get rid of unused levels.
    rownames(SCC4) <- 1:nrow(SCC4)
    fitP.em <- eglhmm(y~locn+depth,data=SCC4,distr="P",cells=c("locn","depth"),
                    K=2,method="em",verb=TRUE)
    \donttest{
        fitP.lm <- eglhmm(y~locn+depth,data=SCC4,distr="P",cells=c("locn","depth"),
                        K=2,verb=TRUE)
        fitD.lm <- eglhmm(formula=y~ma.com+nh.com+bo.com,data=SCC4,nbot=0,ntop=11,
                      cells=c("locn","depth"),distr="Dbd",K=2,method="lm",verb=TRUE,
                      tolerance=NULL)
        SCD4 <- SydColDisc[SydColDisc$locn \%in\% loc4,] 
        SCD4$locn <- factor(SCD4$locn) # Get rid of unused levels.
        fitM.lm  <- eglhmm(formula=y~ma.com+nh.com+bo.com,data=SCD4,
                      cells=c("locn","depth"),distr="Multinom",K=2,
                      verb=TRUE)
        xxx <- split(SCD4,f=SCD4$locn)
        X   <- with(xxx,data.frame(y.LngRf=LngRf$y,y.BondiE=BondiE$y,depth=LngRf$depth))
        fitBiv <- eglhmm(response=c("y.LngRf","y.BondiE"),data=X,K=2,cells="depth",
                         indep=FALSE,verb=TRUE)
    }
# See the help for ionChannelData for more examples involving the
# ion channel data.
}
\keyword{ models }
