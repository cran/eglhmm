eglhmmEm <- function(formula,data,distr,preSpecSigma,size,tau,zeta,phi,
                     mixture=FALSE,itmax=200,crit="CLL",tolerance=NULL,
                     digits=NULL,verbose=FALSE,checkDecrLL=TRUE) {
#
# Function eglhmmEm.  Fit an extended generalised linear hidden
# Markov model using the EM algorithm.
#

# Check that distr is NOT Dbd.
if(distr=="Dbd") {
   whinge <- paste0("The eglhmmEm() function is not applicable when\n",
                    "  \"distr\" is \"Dbd\".\n")
   stop(whinge)
}

# Pick out the index of the stopping criterion:
if(crit=="ABSGRD")
    stop(paste("Stopping criterion \"ABSGRD\" cannot be",
               "used with EM algorithm.\n"))
icrit <- match(crit,c("CLL","L2","Linf"))
if(is.na(icrit)) stop(paste("Stopping criterion",crit,"not recognized.\n"))

# Do some initial housekeeping.
fmla  <- formula
K     <- length(levels(data$state))
npro  <- K*(K-1)
X     <- model.matrix(fmla[-2],data)
nmst  <- names(tau)
nmsz  <- if(is.null(preSpecSigma)) names(zeta) else NULL
nmsp  <- colnames(X)
nms   <- c(nmst,nmsz,nmsp)
y     <- data[,as.character(fmla)[2]]
cf    <- data$cf
state <- data$state
theta <- list(tau=tau,zeta=zeta,phi=phi)
attr(theta,"inclTau")      <- TRUE
attr(theta,"preSpecSigma") <- preSpecSigma
oldtvec <- unlist(theta)
fmt     <- paste0("%.",digits,"f")
scrit   <- numeric(4)

# Set the tolerance if it is NULL.
if(is.null(tolerance)) tolerance <- sqrt(.Machine$double.eps)

# Set the number of digits with which to print out
# "progress reports".
if(is.null(digits)) digits <- ceiling(abs(log10(tolerance)))
fmt       <- paste0("%.",digits,"f")

# Build the first go at the parameter vectors and fy:
X    <- model.matrix(fmla[-2],data=data)
ynm  <- as.character(fmla)[2]
y    <- data[[ynm]]
switch(EXPR=distr,
    Gaussian = {
        gmu      <- X%*%phi
        if(is.null(preSpecSigma)) {
            sigma <- exp(zeta)
        } else {
            sigma <- preSpecSigma
        }
        names(sigma) <- paste0("sigma",1:K)
        sd       <- sigma[data$state]
        fy       <- dnorm(y,mean=gmu,sd=sd)
        modSpecs <- list(phi=phi,sigma=sigma,gmu=gmu,sd=sd)
    },
    Poisson = {
        lam      <- exp(X%*%phi)
        fy       <- dpois(y,lam)
        modSpecs <- list(phi=phi,lambda=lam)
    },
    Binomial = {
        p        <- logistic(X%*%phi)
        fy       <- dbinom(y,prob=p,size=size)
        modSpecs <- list(phi=phi,p=p,size=size)
    },
    Dbd = {
# Can't happen.
    },
    Multinom = {
        Rho      <- phi2Rho(phi,K=K,rhovals=levels(data[[ynm]]),preds=colnames(X))
        fy       <- ffun(data,fmla,response=NULL,Rho,type=1)
        modSpecs <- list(phi=phi,Rho=Rho) # Of course phi is redundant,
                                          # but this makes the syntax
                                          # more consistent than it
                                          # might otherwise be.

    }
)
fy  <- split(fy,f=data[["cf"]])
fy  <- lapply(fy,function(x){x[is.na(x)] <- 1; x})
modSpecs$fy <- fy

# Get the initial log likelihood.
tpm <- getTpm(theta,K)
rp  <- recurse(fy,tpm)
if(any(rp$llc <= 0)) {
    oldll <- -Inf
} else {
    oldll <- sum(log(rp$llc))
}
msge   <- ""

if(verbose){
    cat("\n      Initial set-up completed ...\n")
    cat("\n      Initial log-likelihood: ",
    sprintf(fmt,oldll),"\n\n",sep="")
}

# Update:
em.step <- 1
if(verbose) cat("Repeating ...\n\n")

repeat{ # Swear again; i.e. recurse.
    if(verbose) cat(paste0("EM step ",em.step,":\n"))

# Back up modSpecs and tpm
    modSpecs.bak <- modSpecs
    tpm.bak <- tpm

# Update tpm (and tau).  Note that tau is needed only for
# assessment of the change in the model coefficients.
    tpm <- reviseTpm(rp$xisum,mixture)
    tau <- fixTau(tpm)

# Update the model.
    waits <- rp[["gamma"]]
    if(distr=="Gaussian") waits <- waits/modSpecs$sigma^2
    data$weights <- as.vector(waits)
    modSpecs <- reviseModel(fmla,data,distr,preSpecSigma,size)
# Note that we use log(sigma) in the following construction
# of tvec in order to keep tvec consistent with its use in
# eglhmmBf() and eglhmmLm() in the Gaussian setting.  The
# vector tvec is used only in assessing the amount of change
# in the model since the previous iteration.
    tvec <- if(distr=="Gaussian") {
        if(is.null(preSpecSigma)) {
            zeta <- log(modSpecs$sigma)
            names(zeta) <- paste0("zeta",1:K)
            with(modSpecs,c(tau,zeta,phi))
        } else {
            with(modSpecs,c(tau,phi))
        }
    } else {
        with(modSpecs,c(tau,phi))
    }

# Update rp.
    rp <- recurse(modSpecs$fy,tpm)

# Test for convergence:
    if(any(rp$llc <= 0)) {
        ll <- -Inf
    } else {
        ll    <-  sum(log(rp$llc))
    }
    if(ll < oldll & checkDecrLL) {
        decll  <- sprintf("%.6f",oldll - ll)
        stpnum <- ordinal(em.step)
        whinge <- paste0("\n  At the ",stpnum," step, the log likelihood was ",
                         sprintf("%.6f",ll),"\n",
                         "  which means that there was a DECREASE in the\n",
                         "  log likelihood in the amount of ",decll,". The EM\n",
                         "  algorithm has apparently encountered the problem\n",
                         "  induced by the \"ispd\" parameters which is explained\n",
                         "  in the help. It is possible, although not certain, that\n",
                         "  the results from the penultimate EM step are indeed a\n",
                         "  local maximum, or close to it. It is advisable to\n",
                         "  check on this by applying one of the other methods\n",
                         "  using the results from the penultimate EM step (as\n",
                         "  returned by this function) as starting values.\n")
        message(whinge)
        ll <- oldll
        modSpecs <- modSpecs.bak
        tpm <- tpm.bak
        converged <- NA
        scrit[icrit] <- NA
        nstep <- em.step-1
        msge  <- paste0("These results may be unreliable.  There was an\n",
                        "\"impossible\", decrease in the log likelihood at\n",
                        "the ",stpnum," step.  See the help for a possible\n",
                        "explanation.\n")
        class(msge) <- "kitty"
        break
    }

    scrit[1] <- if(oldll > -Inf) {
                    (ll - oldll)/(abs(ll) + tolerance)
    } else {
        Inf
    }
    scrit[2] <- sqrt(sum((oldtvec-tvec)^2))/(sqrt(sum(tvec^2)) + tolerance)
    scrit[3] <- max(abs(oldtvec-tvec))/(max(abs(tvec)) + tolerance)
    if(verbose){
        cat("     Log-likelihood: ",
                  sprintf(fmt,ll),"\n",sep="")
                  #format(round(ll,digits)),"\n",sep="")
        cat("     Scaled increase in log-likelihood: ",
                  sprintf(fmt,scrit[1]),"\n",sep="")
                  #format(round(scrit[1],digits)),"\n",sep="")
        cat("     Scaled root-SS of change in coef.: ",
                  sprintf(fmt,scrit[2]),"\n",sep="")
                  #format(round(scrit[2],digits)),"\n",sep="")
        cat("     Scaled max. abs. change in coef.: ",
                  sprintf(fmt,scrit[3]),"\n",sep="")
                  #format(round(scrit[3],digits)),"\n",sep="")
        }

    if(scrit[icrit] < tolerance) {
        converged <- TRUE
        nstep     <- em.step
        break
    }

    if(em.step >= itmax) {
        cat("Failed to converge in ",itmax," EM steps.\n",sep="")
        converged <- FALSE
        nstep     <- em.step
        break
    }
    oldtvec <- tvec
    oldll   <- ll
    em.step <- em.step + 1
}

ispd         <- reviseIspd(tpm)
fy           <- modSpecs$fy
data$weights <- NULL

if(checkDecrLL) {
    anomaly <- is.na(converged)
} else {
    anomaly <- scrit[1] < 0
}
th.new <- list(tau=tau,zeta=zeta,phi=modSpecs$phi)
attr(th.new,"inclTau") <- TRUE
attr(th.new,"preSpecSigma") <- preSpecSigma

part1 <- list(tpm=tpm,ispd=ispd,phi=modSpecs$phi,theta=th.new,mixture=mixture,
              log.like=ll,crit=crit,tolerance=tolerance,
              stopCritVal=scrit[icrit],anomaly=anomaly,converged=converged,
              nstep=em.step,message=msge,formula=fmla,distr=distr,fy=fy)
part2 <- switch(EXPR=distr,
    Gaussian = {
        X   <- model.matrix(fmla,data=data)
        gmu <- X%*%modSpecs$phi
        mu  <- getMu(gmu,data,fmla)
        if(is.null(preSpecSigma)) {
            list(mean=modSpecs$gmu,sd=modSpecs$sd,mu=mu,sigma=modSpecs$sigma)
        } else {
            list(mean=modSpecs$gmu,sd=modSpecs$sd,mu=mu,preSpecSigma=modSpecs$sigma)
        }
    },
    Poisson = {
        list(lambda=modSpecs$lambda)
    },
    Binomial = {
        list(p=modSpecs$p,size=size)
    },
    Multinom = {
        list(Rho=modSpecs$Rho)
    }
)
c(part1,part2)
}
