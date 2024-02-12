\name{crossval}
\alias{crossval}
\title{
    Cross validate an extended generalised linear hidden Markov model.
}
\description{
    Calculates a number of cross validation log likelihood values
    for a hidden Markov model (usually one fitted by the \code{eglhmm()}
    function).
}
\usage{
    crossval(model,data,nrep,frac=0.8,type,id="id",minNcomp=100,
             seed=NULL,crossVerb=FALSE,lastPar=NULL,...)
}
\arguments{
  \item{model}{
  A list having components with names selected from those
  of objects returned by \code{eglhmm()} e.g. \code{distr},
  \code{theta}, \code{formula} etc.  Typically \code{model} will
  be of class \code{"eglhmm"} and will have been returned by
  the \code{\link{eglhmm}()} function.
}
  \item{data}{
  A data frame containing the observations to which the cross
  validation are to be fitted.

  It is of course up to the user to ensure that a specified
  value of \code{data} makes sense, i.e. is consistent with the
  other arguments.
}
  \item{nrep}{
  Positive nteger scalar; the number of replications, i.e. the
  number of cross validation calculations undertaken.  If argument
  \code{lastPar} (see below) is supplied then \code{nrep} is ignored
  and is silently set equal to 1.  If \code{lastPar} is \code{NULL}
  then \code{nrep} must be supplied, otherwise an error is thrown.
}
  \item{frac}{Postive scaler, less than 1.  The fraction of the
  data randomly selected to be used as \code{training data} on each
  replication.  (The remaining data, i.e. those data \emph{not}
  used as training data, are used as \code{validation} data.
}
  \item{type}{
  Integer scalar, equal to either 1 or 2, which determines the nature
  of the sampling used to produce the training and validation data.
  If \code{type=1} then these data sets are obtained by sampling
  data points individually.  The training data are obtained by
  setting a fraction \code{1 - frac} of the observed emission values
  (those which are \emph{not} missing already) equal to \code{NA}.
  The validation data are the complement of the training data.

  If \code{type=2} then the \emph{components} of \code{data} are
  randomly sampled (and used in their entirety, either for training
  or for validation).  The components are determined from the column
  of \code{data} the name of which is specified by the argument
  \code{id}; if there is no such column, then an error is thrown.
  Sampling the components means sampling the levels of (the factor)
  \code{data[,id]}.

  Obviously it is sensible to use \code{type=2} \emph{only}
  if \code{data} has a \emph{large} number of components.  By
  default this number is required to be at least 100.  (See
  \code{minNcomp} below.)
}
  \item{id}{
  Character scalar specifying the column to be used to determine the
  individual independent time series that make up the data.  Ignored
  unless \code{type=2}.
}
  \item{minNcomp}{
  Integer scalar specifying the minimum number of components (number
  of levels of \code{data["id"]}) that \code{data} must have in order
  for \code{type=2} to used.  Ignored if \code{type} is equal to 1.

  If the number of components is less than the default value of
  \code{minNcomp} (i.e. 100) then it is strongly recommended that
  \code{type=1} be used instead.
}
  \item{seed}{
  Positive integer scalar to be used as a seed for the random
  number generator (so as to enable reproducibility).  If not
  supplied, it is randomly chosen from the sequence \code{1:1e6}.
  Note that if \code{nrep > 1} then after this seed is set, a vector
  \code{SEEDS} of \dQuote{auxiliary} seeds, of length \code{nrep},
  is chosen from the sequence \code{1:1e6} and the seed is set
  from the corresponding entry of this vector at the start of each
  replication.  If \code{nrep==1} then the sampling for the single
  replication that occurs is determined by \code{seed}.
}
  \item{crossVerb}{
  Logical scalar.  Should brief \dQuote{progress reports} (letting
  the user know what is happening with respect to replicate \code{repl},
  for each \code{repl}) be produced?
}
  \item{lastPar}{
  The last values of the (relevant) fitted paramenters, provided as
  an attribute of a component of the list returned by the function
  currently under consideration (i.e. \code{crossval()}, whenever
  the process of fitting the model in question to the training data
  did not converge.  These values can be used as starting values so
  as to carry on with the fitting process from where the previous
  attempt left off.
}
  \item{\dots}{
  Possible additional arguments to be passed to \code{eglhmm()} via
  \code{update()}, e.g. \code{itmax}, \code{tolerance}, \dots .
}
}

\section{Note}{
  If the function fails to fit the model, obtained from the training
  data, to the validation data, then the value returned is \code{NA}.
  This value will have an attribute \code{lastPar}.  This attribute
  will in turn have attributes, \code{trnDat} and \code{valDat},
  the training data and validation data which were being used
  in the failed fitting procedure.  Supplying an appropriate
  value of \code{lastPar} enables the continuation of the fitting
  procedure, starting from where the procedure previously left off.
  See \bold{Details} for a little more information.
}

\details{
   On each replication a random subset comprising \code{frac} of
   the data is selected to serve as training data.  The complement
   of this subset is used as validation data.  The model specified
   by \code{model} is fitted to the training data.  It is possible
   to over-ride some of the details of the specifications producing
   \code{model}, via the \code{\dots} argument of \code{crossval()}.
   After the model is fitted to the training data, the log likelihood
   of the validation data is calculated on the basis of that
   fitted model.

   If the procedure for fitting a model to the training data fails
   to converge, then the corresponding entry of the list returned by
   this function is \code{NA}.  In this case, the entry is assigned an
   attribute \code{lastPar}, (the estimates of the model parameters
   that were current when the fitting algorithm terminated)
   which will in turn have attributes \code{trnDat}, \code{valDat}
   (the training data in question and the corresponding validation
   data), and \code{seed} (the value of the seed that was set before
   the sampling that determined \code{trnDat} and \code{valDat}
   took place).  The value of \code{seed} is \code{SEEDS[i]}
   (if \code{nrep>1} and the entry in question was the \code{i}th
   entry of the returned list) or the value of the \code{seed}
   argument of this function or its random replacement if this
   argument was not supplied (if \code{nrep==1}).

   The attribute \code{lastPar} enables the user to continue the
   procedure for fitting a model to the training data, starting
   from where the procedure, that failed to converge, left off.
   Continuing the procedure is easily effected by calling
   \code{crossval()} with argument \code{par0} set equal to the
   \code{lastPar} attribute of the relevant entry of the list that
   was previously returned by this function.

   If \code{type==1} then the training and validation data are
   created in a somewhat subtle manner.  The procedure necessitates
   referring to the \dQuote{original} data.  The data frame that
   is passed to \code{eglhmm()} is a \dQuote{replicated} version
   of the original data, with each row of the original data being
   repeated once for each level of \code{state} (and a \code{"state"}
   column --- factor --- being added to the resulting data frame).
   If \code{type==2} the procedure is conceptually simpler.
   The procedure in the \code{type==1} setting is as follows:
   \itemize{
   \item Let \eqn{S} be the set of indices of all non-missing values
   in the column of the original data that contains the emissions.
   \item Select at random a subset \eqn{V} of \eqn{S} so that
   the size of \eqn{V} is the fraction \code{frac} of the size
   of \eqn{S}.
   \item Let \eqn{T} be the complement in \eqn{S} of \eqn{V}.
   \item The training data are formed by replacing by \code{NA} all
   those values of the emissions column in the original data whose
   indices are in \eqn{T}.
   \item The validation data are formed
   by replacing by \code{NA} all those values of the emissions
   column in the original data, whose indices are in \eqn{V}.
   \item Then replicate both the training and validation data in
   the manner described above.

   If \code{type==2} then the training data are formed by selecting
   at random a fraction \code{frac} of the levels of the column of
   \code{data} named \code{"id"}.  (If there is no such column, an
   error is thrown.)  The training data then consist of those rows
   of \code{data} corresponding to the selected levels of \code{id}.
   The validation data then consist of those rows of \code{data}
   which are not in the training data.
   }
}
\value{
   If \code{nrep>1} the returned value is list of length \code{nrep}.
   The \code{i}th entry of this list is the log likelihood of the
   validation data with respect to the model fitted to the training
   data, for the \code{i}th random selection of these two subsets.
   This entry will be \code{NA} if the attempt to fit a model
   to the training data was unsuccessful.  The \code{i}th entry
   has an attribute \code{seed} (singular) which is the value of
   the seed that was set prior to the random sampling that chose
   the training and validation data.  If the \code{i}th entry is
   \code{NA} it will also have an attribute \code{lastPar} which
   in turn will have attributes \code{trnDat} and \code{valDat}.
   See \bold{Details}.

   If \code{nrep>1} then the returned value also has an attribute
   \code{seeds} (plural) which is vector of length \code{nrep+1},
   consisting of the \dQuote{auxiliary} seed vector \code{SEEDS}
   (see the argument \code{seed}) together with the \dQuote{over all}
   seed (possibly equal to the \code{seed} argument) that was set for
   the random number generator before any sampling was undertaken.
   Note that the \code{i}th entry of this \code{seeds} attribute
   is the same as the \code{seed} attribute of the \code{i}th entry
   of the returned valuel

   If \code{nrep==1} then the returned value is a single numeric
   scalar which is the log likelihood of the validation data or
   \code{NA} if the fitting procedure did not converge for the
   training data.  It has an attribute \code{seed} which is equal
   to the \code{seed} argument or its random replacement.  If the
   value is \code{NA} then it has a further attribute \code{lastPar}.
   (See above.)
}

\references{
Celeux, Gilles and Durand, Jean-Baptiste (2008).  Selecting
hidden Markov model state number with cross-validated likelihood.
\emph{Computational Statistics} \bold{23} 541--564, \code{DOI
10.1007/s00180-007-0097-1}.

Smyth, Padhraic (2000).  Model selection for probabilistic clustering
using cross-validated likelihood.  \emph{Statistics and Computing}
\bold{9} 63--72.
}

\author{Rolf Turner
  \email{r.turner@auckland.ac.nz}\cr
}

\seealso{
   \code{\link{eglhmm}()}
}
\examples{
\dontrun{
ids   <- paste0("s",1001:1101)
cc  <- ccSim[ccSim$id %in% ids,]
cc$id <- factor(cc$id)
cvll1 <- vector("list",9)
set.seed(42)
SEEDS <- sample(1:1e6,9)
for(k in 1:9) {
    cat("k =",k,"started\n")
    fit  <- eglhmm(categMC ~ 1, distr="M", method="em", data=cc, K=k,
                   itmax=1500,cells="id",verb=TRUE)
    cvll1[[k]] <- crossval(fit,nrep=5,type=2,seed=SEEDS[k],tolerance=1e-4,
                           verbose=FALSE,crossVerb=TRUE)
    cat("k =",k,"finished\n")
}
}
fit   <- eglhmm(y ~ 1, data=Downloads,K=4,distr="P",verb=TRUE,cf="singlecell")
# Use artifically low value of itmax so that crossval() fails to
# fit the model to the training data
cvll2 <- crossval(fit,nrep=5,type=1,verbose=TRUE,seed=322596,itmax=5)
cvll3 <- crossval(fit,type=1,verbose=TRUE,lastPar=attr(cvll2[[1]],"lastPar"))
# So cvll3 carried on, in one instance, from where the first
# attempted fit in cvll2 gave up.

}
