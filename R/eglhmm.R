eglhmm <- function(formula=NULL, response=NULL, data,
                  distr=c("Gaussian","Poisson","Binomial","Dbd","Multinom","discnp"),
                  inclTau=TRUE,preSpecSigma=NULL, indep=NULL, size=NULL, nbot=NULL,
                  ntop=NULL, cells=NULL, cf="singlecell", K=NULL, par0=NULL,
                  randStart=NULL, method=c("lm","em","bf"),
                  optimiser=c("optim","nlm"), optimMethod="BFGS",
                  nlmWarn=FALSE,lmc=10, tolerance=NULL, digits=NULL,
                  verbose=FALSE,itmax=200, contrast=c("treatment","sum","helmert"),
                  crit=c("CLL", "L2", "Linf","ABSGRD"), breaks=NULL,
                  hessian=FALSE,useAnalGrad=FALSE, ca=FALSE,
                  checkDecrLL=TRUE) {
#
# Function eglhmm.  To carry out the fitting of an extended
# generalised hidden Markov model structure given data from a
# distribution the parameter(s) of which depends in part upon the
# unobserved state of an associated Markov chain.
#

# Get the call (so that update() can be used).
ccc <- match.call()

# Check on presence of K.
if(is.null(K)) {
    if(is.null(par0)) stop("One of \"K\" or \"par0\" must be supplied.\n")
    K <- nrow(par0[["tpm"]])
}

# Check that "response" is of the right length.
nresp <- length(response)
if(!(nresp %in% 0:2))
    stop("Argument \"response\" is of the wrong length.\n")
bivar <- length(response) == 2

if(is.null(formula)) {
    if(is.null(response))
        stop("One of \"formula\" or \"response\" must be specified.\n")
    if(!bivar) {
        formula <- if(K>1) {
            as.formula(paste0(response," ~ 0+state"))
        } else {
            as.formula(paste0(response," ~ 1"))
        }
    }
} else {
    if(bivar) {
        formula <- NULL
    } else {
        cform <- as.character(formula)
        if(cform[3]=="1") {
            if(K > 1) {
                cform[3] <- "0+state"
            }
        } else {
            if(grepl("state",cform[3])) {
                whinge <- paste0("The model formula should NOT contain \"state\".\n",
                                 "  The \"state\" predictor is added internally.\n")
                stop(whinge)
            }
            if(K > 1) cform[3] <- paste0("0+state+",cform[3])
        }
        formula <- as.formula(paste(cform[2],cform[3],sep=" ~ "))
        if(is.null(response)) response <- cform[2]
    }
}

# Get the distribution.
if(bivar) {
    distr <- "Multinom"
} else {
    distr <- match.arg(distr)
    if(distr == "discnp") distr <- "Multinom"
}
if(distr == "Dbd") {
    if(!requireNamespace("dbd"))
        stop("Package dbd is needed and seems not to be available.\n")
}

# If distr is "Gaussian" then check that preSpecSigma is of the right
# length and that it has strictly positive entries.  Otherwise set
# preSpecSigma equal to NA.
if(distr == "Gaussian") {
    if(!is.null(preSpecSigma)) {
        if(length(preSpecSigma) != K) {
            stop("Argument \"preSpecSigma\" must be of length K.\n")
        }
        if(!all(preSpecSigma > 0)) {
            if(K > 1) {
                stop("All entries of \"preSpecSigma\" must be strictly positive.\n")
            } else {
                stop("Argument \"preSpecSigma\" must be strictly positive.\n")
            }
        }
    }
} else {
    preSpecSigma <- NA
}

# For the Multinom distribution, coerce the response(s) to be
# factor valued.  But don't lose any fucking levels if the response(s)
# is/are already factors!
if(distr=="Multinom") {
    if(bivar) {
        for(i in 1:2) {
            data[[response[i]]] <- as.factor(data[[response[i]]])
        }
    } else {
        data[[response]] <- as.factor(data[[response]])
    }
}

# Get "rsplvls" (the possible values of the response) in the context
# of the Multinom distribution.
if(distr == "Multinom") {
    if(bivar) {
        rsplvls <- list(levels(data[[response[1]]]),levels(data[[response[2]]]))
    } else {
        rsplvls <- levels(data[[response]])
    }
} else {
    if(bivar) {
        whinge <- paste0("A bivariate response is permitted only for the ",
                         "Multinom distribution.\n")
        stop(whinge)
    }
    rsplvls <- NULL
}

# Get the method.
    if(bivar) {
        method <- "em"
    } else {
        method <- match.arg(method)
    }
    if(distr=="Dbd") {
        if(method=="em") {
            stop("Cannot use the \"em\" method when \"distr\" is \"Dbd\".\n")
        }
    }

# Check on the availability of "indep".
if(bivar & is.null(indep))
    stop("When the data are bivariate \"indep\" must be specified.\n")

# Set up the multiplier for BIC.
if(bivar) {
    nobs <- with(data,sum(!is.na(get(response[1]))) + sum(!is.na(get(response[2]))))/2
} else {
    nobs <- with(data,sum(!is.na(get(response))))
}
bicm <- log(nobs)

# Calculate the "missing fraction".
nafrac <- nafracCalc(data,response)

# Set the contrast:
contrast <- match.arg(contrast)
contr    <- c(paste("contr",contrast,sep="."),"contr.poly")
OC       <- options(contrasts=contr)
on.exit(options(OC))

# Do the K=1, case --- iid data.
if(K==1) {
    attr(data,"rsplvls") <- rsplvls
    rslt <- doKeq1(data=data,fmla=formula,distr=distr,preSpecSigma=preSpecSigma,
                   response=response,indep,size,nbot,ntop,bicm,nafrac)
} else {

# Add the cells factor "cf" to "data".
    if(is.null(cells)) {
        if(inherits(cf,"character")) {
            ssc <- identical(cf,"singlecell") # ssc = set single cell
        } else if(is.null(cf)) {
            ssc <- TRUE
            whinge <- paste0("Neither \"cells\" nor \"cf\" has been specified.\n",
                             "  Therefore \"cf\" has been set equal to\n",
                             "  factor(rep(1,NR)) where NR is the number of rows\n",
                             "  of the data frame \"data\".  This results in the\n",
                             "  observation sequence being treated as a single\n",
                             "  time series.  If you have something different in\n",
                             "  mind, then specify either \"cells\" or \"cf\".\n")
            warning(whinge)
        } else ssc <- FALSE
        if(ssc) cf <- factor(rep(1,nrow(data)))
        data[,"cf"] <- cf # Use either the argument "cf" or the default.
   } else {
       if(!all(cells %in% names(data))) {
            stop("Some of the \"cells\" names are not in names(data).\n")
       }
       data[,"cf"] <- interaction(data[,cells])
   }
    
    # Build the replicated data frame; first a temporary version.
    M         <- nrow(data)
    tmp       <- as.data.frame(lapply(data,function(x,K,M){rep(x,rep(K,M))},K=K,M=M))
    tmp$state <- factor(rep(1:K,M))
    
    # Now tidy the cell structure of the temporary version.  (This is
    # a *CRUCIAL* step, previously --- before 25/03/2023 --- missing).
    # Note that this step puts the cells of the data into contiguous rows,
    # although these rows are still in the same order *within* the
    # cells as they were in the original data.  This *has* to be so;
    # since the observations within each cell form a *time series* whence
    # the order within cells is of fundamental importance.  The cells
    # are independent of each other so their order is of no theoretical
    # importance.  However consideration of their order is of vital
    # importance in terms of the proper functioning of the code.  It is
    # imperitive that the observations and the derived quantities "fy"
    # and "gamma" be properly aligned.  Note that the cells now occur in
    # the order of the levels of the cell factor "cf".
    if(length(levels(tmp[["cf"]])) > 1) {
        stmp  <- split(tmp,f=tmp[["cf"]])
        rd.df <- do.call(rbind,stmp)
    } else {
        rd.df <- tmp
    }
    rownames(rd.df) <- 1:nrow(rd.df)
    attr(rd.df,"rsplvls") <- rsplvls
    
    # Initialize the parameters.
    parz <- initialise(distr=distr,data=data,formula=formula,rsplvls=rsplvls,
                       par0=par0,K=K,preSpecSigma=preSpecSigma,size=size,
                       nbot=nbot,ntop=ntop,breaks=breaks,randStart=randStart,
                       indep=indep)
    
    # Get the stopping criterion.
    crit <- match.arg(crit)
    if(crit=="ABSGRD" & method != "lm")
        stop("The \"ABSGRD\" criterion can only be used when \"method\" = \"lm\".\n")
    
    # Dig out the parameters to be optimised.
    # Fix up the following:
    if(bivar) {
        tau <- zeta <- phi <- NULL # No actual effect; just cosmetic code!
    } else {
        tau  <- parz$tau
        if(distr=="Multinom") {
            zeta <- NULL
            Rho  <- parz$Rho
            phi  <- rho2Phi(Rho)
        } else if(distr=="Gaussian") {
            phi <- parz$phi
            if(is.null(preSpecSigma)) {
                zeta <- log(parz$sigma)
                names(zeta) <- paste0("zeta",1:K)
            } else {
                zeta <- NULL
            }
        } else {
            zeta <- NULL
            phi  <- parz$phi
        }
    }
    
    if(bivar) {
    # BI <--> "bivariate independent".
    # BD <--> "bivariate dependent".
        if(indep) {
            rslt <- eglhmmBI(data=rd.df,par0=parz,K=K,response=response,
                          tolerance=tolerance,digits=digits,verbose=verbose,itmax=itmax,
                          crit=crit)
        } else {
            rslt <- eglhmmBD(data=rd.df,par0=parz,K=K,response=response,
                          tolerance=tolerance,digits=digits,verbose=verbose,itmax=itmax,
                          crit=crit)
        }
    } else {
    # Note: should tidy up the order of the argument lists in the following calls.
        switch(EXPR=method,
            lm = {
    	    rslt <- eglhmmLm(formula=formula,data=rd.df,distr=distr,
                                 inclTau=inclTau,preSpecSigma=preSpecSigma,size=size,
                                 nbot=nbot,ntop=ntop,tau=tau,zeta=zeta,phi=phi,
                                 lmc=lmc,itmax=itmax,crit=crit,tolerance=tolerance,
                                 digits=digits,verbose=verbose)
                 },
            em = {
    	    rslt <- eglhmmEm(formula=formula,data=rd.df,distr=distr,
                                 preSpecSigma=preSpecSigma,size=size,
                                 tau=tau,zeta=zeta,phi=phi,itmax=itmax,
                                 crit=crit,tolerance=tolerance,
                                 digits=digits,verbose=verbose,
                                 checkDecrLL=checkDecrLL)
                 },
            bf = {
                optimiser <- match.arg(optimiser)
                theta <- list(tau=tau,zeta=zeta,phi=phi)
                attr(theta,"inclTau")      <- inclTau
                attr(theta,"preSpecSigma") <- preSpecSigma
    	        rslt <- eglhmmBf(formula=formula,data=rd.df,distr=distr,
                                 theta=theta,size=size,nbot=nbot,ntop=ntop,
                                 optimiser=optimiser,optimMethod=optimMethod,
                                 nlmWarn=nlmWarn,hessian=hessian,
                                 useAnalGrad=useAnalGrad,ca=ca,itmax=itmax,
                                 tolerance=tolerance,verbose=verbose)
                 },
                 {
    	         stop(paste("Unrecognized method ",method,".\n",sep=""))
                 }
        )
    }
    # Set AIC and BIC.
    npar <- if(bivar) {
        Rho <- parz$Rho
        if(indep) {
            K*(K-1) + K*(sum(sapply(Rho,ncol))-2)
        } else {
            K*(K-1) + prod(dim(Rho))-K
        }
    } else {
        length(tau)*inclTau + length(zeta) + length(phi)
    }
    ll   <- rslt$log.like
    AIC  <- -2*ll + 2*npar
    BIC  <- -2*ll + bicm*npar
    rslt$contrast   <- if(is.null(formula)) NULL else contrast
    rslt$par0       <- parz[names(parz) != "tau"]
    rslt$method     <- method
    rslt$data       <- rd.df
    rslt$bicm       <- bicm
    rslt$npar       <- npar
    rslt$AIC        <- AIC
    rslt$BIC        <- BIC
    rslt$cells      <- cells
    rslt$missFrac   <- nafrac
    msge            <- rslt$message
    if(!is.null(msge) && msge == "") rslt$message <- NULL
    if(method == "em") rslt$checkDecrLL <- checkDecrLL
}
rslt$call       <- ccc
    
# Get the components of the returned object to be always in
# the same order; "nmsiro" <--> "names in right order".
nmsiro <- c("call","tpm","ispd","phi","theta","Rho","log.like","gradient","numHess",
            "Hessian","mu","sigma","preSpecSigma","stopCritVal",
            "anomaly","converged","nstep","mean","sd","lambda","p","alpha",
            "beta","fy","message","par0","cells","formula","response",
            "distr","nbot","ntop","size","tolerance","crit","mixture",
            "contrast","method","checkDecrLL","data","bicm","npar","AIC",
            "BIC","missFrac")
iii    <- match(nmsiro,names(rslt),nomatch=0)
rslt   <- rslt[iii]
class(rslt) <- if(bivar) {
                   c("eglhmm","eglhmm.bivariate")
               } else {
                   "eglhmm"
               }
rslt
}
