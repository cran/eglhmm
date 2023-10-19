eglhmmBI <- function(data,par0,K,response,tolerance,digits,verbose,itmax,crit) {
#
# Function eglhmmBI.  To conduct the fitting of a Hidden Markov
# model when the observations are bivariate and assumed to be
# conditionally independent.  Called by eglhmm().  Not intended to
# be called directly.  Note that the application of this function
# is appropriate only when the data are discrete (i.e. only when
# the distribution involved is "Multinom".

# Pick out the index of the stopping criterion:
icrit <- match(crit,c('CLL','L2','Linf'))
if(is.na(icrit)) stop(paste("Stopping criterion",crit,"not recognized.\n"))

# Set the tolerance if it is NULL.
if(is.null(tolerance)) tolerance <- 1e-6

# Perform initial setting-up.
tpm   <- par0$tpm
ispd  <- reviseIspd(tpm)
Rho   <- par0$Rho
m1    <- nrow(Rho[[1]])
m2    <- nrow(Rho[[2]])

# Set the number of digits with which to print out
# "progress reports".
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

# Set the level below which probabilities are considered to
# be noise:

old.theta <- c(as.vector(tpm[,-K]),as.vector(Rho[[1]][1:(m1-1),]),
                                   as.vector(Rho[[2]][1:(m2-1),]))
fy        <- ffun(data,fmla=NULL,response,Rho,type=2)
fy[is.na(fy)] <- 1
fy <- split(fy,f=data[["cf"]])
rp        <- recurse(fy,tpm)
old.ll    <- getLL(rp)
msge      <- ""

# Ready to go.
    if(verbose){
        cat("\n      Initial set-up completed ...\n")
        cat("\n      Initial log-likelihood: ",
        format(round(old.ll,digits)),"\n\n",sep="")
    }

# Revise:
em.step <- 1
if(verbose) cat('Repeating ...\n\n')
scrit <- numeric(3)
repeat{
    if(verbose) cat(paste('EM step ',em.step,':\n',sep=''))

# Calculate the parameters.
    tpm  <- reviseTpm(rp$xisum,mixture=FALSE)
    ispd <- reviseIspd(tpm)
    waits <- as.vector(rp$gamma)
    data$weights <- waits
    Rho   <- reviseRho(data=data,response=response,fmla=NULL,type=2)

# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
    fy <- ffun(data,fmla=NULL,response=response,Rho=Rho,type=2)
    fy[is.na(fy)] <- 1
    fy  <- split(fy,f=data[["cf"]])
    rp  <- recurse(fy,tpm)
    ll  <-  getLL(rp)
    scrit[1]  <- llDiff(ll,old.ll,tolerance)
    if(ll < old.ll) {
        decll  <- round(old.ll - ll, digits)
        stpnum <- ordinal(em.step)
        whinge <- paste0("\n  At the ",stpnum," step, there was a DECREASE in the\n",
                         "  log likelihood in the amount of ",decll,". The EM\n",
                         "  algorithm has apparently encountered the problem\n",
                         "  induced by the \"ispd\" parameters which is explained\n",
                         "  in the help. It is possible, although not certain, that\n",
                         "  the results from the penultimate EM step are indeed a\n",
                         "  local  maximum, or close to it.  Ideally one would\n",
                         "  check this by applying a different method, but\n",
                         "  unfortunately for bivariate data only the EM method\n",
                         "  is available.  Hence this function returns the\n",
                         "  results from the penultimate EM step, that being\n",
                         "  the best that can be done.\n")
        message(whinge)
        ll <- old.ll
        converged <- NA
        nstep <- em.step-1
        msge  <- paste0("These results may be unreliable.  There was an\n",
                        "\"impossible\", decrease in the log likelihood at\n",
                        "the ",stpnum," step.  See the help for a possible\n",
                        "explanation.\n")
        class(msge) <- "kitty"
        break
    }

# Test for convergence:
    new.theta <- c(as.vector(tpm[,-K]),as.vector(Rho[[1]][1:(m1-1),]),
                                           as.vector(Rho[[2]][1:(m2-1),]))
    scrit[2]  <- sqrt(sum((old.theta-new.theta)^2))/
                     (sqrt(sum(new.theta^2)) + tolerance)
    scrit[3]  <- max(abs(old.theta-new.theta))/
                     (max(abs(new.theta)) + tolerance)
    if(verbose){
        cat('     Log-likelihood: ',
                format(round(ll,digits)),'\n',sep='')
        cat('     Scaled percent increase in log-likelihood: ',
            format(round(scrit[1],digits)),'\n',sep='')
        cat('     Scaled root-SS of change in coef.: ',
            format(round(scrit[2],digits)),'\n',sep='')
        cat('     Scaled max. abs. change in coef.: ',
            format(round(scrit[3],digits)),'\n',sep='')
    }

    if(scrit[icrit] < tolerance) {
            converged <- TRUE
            nstep <- em.step
            break
        }

    if(em.step >= itmax) {
        cat('Failed to converge in ',itmax,' EM steps.\n',sep='')
        converged <- FALSE
        nstep <- em.step
        break
    }

# Replace the ``old'' parameter and log likelihood values
# by the new ones.
    old.theta <- new.theta
    old.ll    <- ll
# Increment the step number.
    em.step   <- em.step + 1
}

npar <- K*(K-1) + K*(sum(sapply(Rho,nrow))-2)
rslt <- list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,par0=par0,npar=npar,
             crit=crit,tolerance=tolerance,stopCritVal=scrit[icrit],
             converged=converged,nstep=nstep,message=msge,distr="Multinom",
             response=response)
rslt
}
