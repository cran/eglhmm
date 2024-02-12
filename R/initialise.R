initialise <- local({

initGauss <- function(fmla,data,phi,sigma,preSpecSigma) {
    if(is.null(phi) | is.null(sigma)) {
        rslt <- glm(formula=fmla,data=data,family=gaussian,na.action=na.exclude)
    }
    if(is.null(phi)) phi <- coef(rslt)
    if(is.null(preSpecSigma)) {
        if(is.null(sigma)) { # Here's where the rubber hits the road!
            ynm   <- as.character(fmla[2])
            y     <- data[[ynm]]
            y     <- y[!is.na(y)]
            state <- data[["state"]]
            state <- state[!is.na(state)]
            X     <- model.matrix(fmla,data=data)
            sigma <- revSigUnwtd(phi=phi,X=X,y=y,state=state)
        }
        rslt <- list(phi=phi,sigma=sigma)
   } else {
        rslt <- list(phi=phi)
   }
   rslt
}

initDbd <- function(y,K,nbot,ntop) {
# Form the overall mle of phi.
phi0 <- dbd::mleDb(unlist(y),ntop=ntop,zeta=nbot==0)
rslt <- vector("list",2)
names(rslt) <- c("alpha","beta")
# Perturb phi0 to get a different value for each state.
for(i in 1:2) {
    ri <- sort(phi0[i]*c(0.25,0.75))
    if(diff(ri) < 0.1) {
        incr <- (0.1 - diff(ri))/2
        ri   <- ri + incr*c(-1,1)
    }
    si <- seq(ri[1],ri[2],length=K)
    rslt[[i]] <- numeric(K)
    for(k in 1:K) {
        rslt[[i]][k] <- si[k]
    }
}
rslt
}

function(distr,data,formula,rsplvls,par0,K,preSpecSigma,size,nbot,ntop,
         breaks,randStart,indep) {
#
# Set fmla (for parallelism with eglhmm{Bf,Em,Lm}.R).
fmla <- formula
# Make randStart if necessary.  Note that randStart must either be
# NULL, or a logical scalar, or a list of three logical scalars
# named "tpm", "phi", and "Rho".
if(is.null(randStart)) {
    randStart <- list(tpm=FALSE,phi=FALSE,Rho=FALSE)
} else if(length(randStart)==1 && is.logical(randStart)) {
    randStart <- list(tpm=randStart,phi=randStart,Rho=randStart)
} else if(inherits(randStart,"list")) {
    nms <- sort(names(randStart))
    if(!identical(nms,c("phi","Rho","tpm"))) {
        stop("The names of argument \"randStart\" are inappropriate.\n")
    }
    ok <- sapply(randStart,function(x){is.logical(x) & length(x)==1})
    if(!ok) {
        whinge <- paste0("Some components of argument \"randStart\" are",
                         " either non-logical or non-scalar.\n")
        stop(whinge)
    }
} else {
    whinge <- paste0("Argument \"randStart\" must either be NULL, a logical\n",
                     "  scalar or a list of three logical scalars.\n")
    stop(whinge)
} 

# Basic checks:
if(is.null(par0)) {
    if(is.null(K)) {
        stop("One of \"par0\" or \"K\" must be supplied.\n")
    } else {
        if(randStart$tpm)
            tpm <- matrix(runif(K*K),K,K)
        else tpm <- matrix(1/K,K,K) + 1.5*diag(K)
# I have written a function est.tpm() to produce an initial value
# for tpm, using "fake states".  (See below.)  The value returned by
# this function could possibly be used instead of the ad hoc value
# produced by the preceding line.  Would it be worth the effort?
        tpm <- tpm/apply(tpm,1,sum)

# Note: "phi" is taken care of later on.
    }
} else {
    tpm <- par0$tpm
    ok  <- checkTpm(tpm) # "ok" will always be TRUE; if it weren't
                         # an error would have been thrown by checkTpm().
    if(is.null(K)) {
        K <- nrow(tpm)
    }  else if(K != nrow(tpm)) {
        stop("The values of \"K\" and \"par0$tpm\" are inconsistent.\n")
    }
    if(!is.null(par0$distr) && distr!=par0$distr)
        stop("The values of \"distr\" and \"par0$distr\" are inconsistent.\n")
}

# Check whether the response is univariate.
if(inherits(rsplvls,"list")) {
    if(length(rsplvls) == 2) {
        univar <- FALSE
    } else { 
        stop("If \"rsplvls\" is a list, it must be of length 2.\n")
    }
} else univar <- TRUE

# If univariate, add "fake" states to the data.
if(univar & distr != "Multinom") {
    data <- fakeStates(data=data,K=K,fmla=fmla)
}

# Calculate the exponential form "tau" of the tpm entries:
tau <- if(K > 1) fixTau(tpm) else NULL

# Allow for random start for phi
if(is.null(par0$phi) & randStart[["phi"]]) {
    dumX   <- model.matrix(fmla,data=data[1,])
    nphi   <- ncol(dumX)
    par0$phi <- switch(EXPR=distr,
        Gaussian = rnorm(nphi),
        Poisson  = rnorm(nphi),
        Binomial = rnorm(nphi),
        Dbd      = rnorm(2*nphi),
        Multinom = NA
    )
}

# Run through the possible distributions.
switch(EXPR=distr,
    Gaussian = {
        phi   <- par0$phi
        sigma <- par0$sigma

# The possibility of par0$phi and/or par0$sigma being NULL
# is taken care of in initGauss().
        if(is.null(preSpecSigma)) {
            rslt <- initGauss(fmla,data,phi,sigma=sigma,preSpecSigma=NULL)
        } else {
            rslt <- initGauss(fmla,data,phi,sigma=sigma,preSpecSigma=preSpecSigma)
        }
# Note that if preSpecSigma is provided (is not NULL) then pars has
# no sigma component, i.e. pars$sigma is NULL.
        pars  <- c(list(tpm=tpm,tau=tau),rslt)
    },
    Poisson = {
        if(is.null(par0$phi)) {
            rslt  <- glm(fmla,data=data,family=poisson,na.action=na.omit)
            phi   <- coef(rslt)
        } else phi <- par0$phi
        pars <- list(tpm=tpm,tau=tau,phi=phi)
    },
    Binomial = {
        if(is.null(size)) {
            if(is.null(par0$size)) {
                stop("When \"distr\" is \"Binomial\", \"size\" must be supplied.\n")
            }
            size <- par0$size
        } else if(!is.null(par0$size) && par0$size != size) {
            stop("The values of \"size\" and \"par0$size\" are inconsistent.\n")
        }
        if(is.null(par0$phi)) {
            bfmla <- binForm(fmla)
            rslt  <- glm(formula=bfmla,data=cbind(data,size=size),family=binomial,
                        na.action=na.omit)
            phi  <- coef(rslt)
        } else phi <- par0$phi
        pars <- list(tpm=tpm,tau=tau,phi=phi,size=size)
    },
    Dbd = {
        if(is.null(nbot)) {
            if(is.null(par0$nbot)) {
                stop("When \"distr\" is \"Dbd\", \"nbot\" must be supplied.\n")
            }
            nbot <- par0$nbot
        } else if(!is.null(par0$nbot) && par0$nbot != nbot) {
            stop("The values of \"nbot\" and \"par0$nbot\" are inconsistent.\n")
        }
        if(is.null(ntop)) {
            if(is.null(par0$ntop)) {
                stop("When \"distr\" is \"Dbd\", \"ntop\" must be supplied.\n")
            }
            ntop <- par0$ntop
        } else if(!is.null(par0$ntop) && par0$ntop != ntop) {
            stop("The values of \"ntop\" and \"par0$ntop\" are inconsistent.\n")
        }
        dumX  <- model.matrix(fmla,data=data[1,])
        nms   <- dimnames(dumX)[[2]]
        if(is.null(par0$phi)) {
            ynm   <- as.character(fmla)[2]
            rslt  <- initDbd(data[[ynm]],K,nbot,ntop)
            npad  <- length(grep("state",nms,invert=TRUE))
            zpad  <- rep(0,npad)
            phi   <- c(rslt[[1]],zpad,rslt[[2]],zpad)
        } else phi <- par0$phi
        names(phi) <- c(paste0("alpha.",nms),paste0("beta.",nms))
        pars <- list(tpm=tpm,tau=tau,phi=phi,nbot=nbot,ntop=ntop)
    },
    Multinom = {
# In the univariate case Rho should be a k x m matrix. The first
# dimension k is the number of predictors, i.e. the number of columns
# of the design matrix X. Note that if there are no "auxilliary"
# predictors, then "state" is the only predictor, so that k = K,
# the number of states.  The second dimenson m is the number of
# possible values y_1, ..., y_m of "Y", the "emissions" variate.
# If M = X%*%Rho, then the i-th row of M is a vector ("expressed
# in exponential form") of coefficients determining Pr(Y = y_j)
# (given the i-th vector of predictors), j = 1, ..., m.  Note that
# the m-th column of Rho, and hence of M, is identically 0.
#
# In the bivariate independent case, Rho is a list of two matrices
# like unto the matrix described above for the univariate case.
# (Do/should we allow "auxilliary" predictors in the bivariate
# independent case?)
#
# In the bivariate dependent case Rho is a three dimentional array
# of probabilities, Rho[i,j,k] = Pr(Y_1 = y_i, Y_2 = y_j | state k).
# No predictors are allowed in the bivariate dependent case.

        if(!univar & is.null(indep)) {
            stop("In the bivariate setting \"indep\" must be specified.\n")
        }

        if(is.null(par0$Rho)) {
            Rho <- initRho(data,K,fmla,randStart=randStart[["Rho"]],
                           indep=indep,rsplvls=rsplvls)
        } else {
            Rho <- par0$Rho
            ok  <- checkRho(par0$Rho,K,rsplvls,indep) # "ok" will always be TRUE;
                                                   # if it weren't, an error
                                                   # would have been thrown.
        }
        pars <- list(tpm=tpm,tau=tau,Rho=Rho)
    }
)
pars
}
})
