reglhmm.default <-function(x,formula,response,cells=NULL,data=NULL,nobs=NULL,
                          distr=c("Gaussian","Poisson","Binomial","Dbd","Multinom"),
                          phi,Rho,sigma,size,ispd=NULL,ntop=NULL,zeta=NULL,
                          missFrac=0,fep=NULL,
                          contrast=c("treatment","sum","helmert"),...) {
#
# Function reglhmm (r for random) to simulate data from a hidden
# generalized linear Markov model.
#

# Make sure that the kontrasts are konsistent.  The interpretation
# of the coefficients of the linear predictor contained in object$phi
# depends on the contrast used.  So the contrast must be the
# same as what you had in mind when you specified phi!!!
# Set the contrasts.
contrast <- match.arg(contrast)
contr    <- c(paste("contr",contrast,sep="."),"contr.poly")
old.con  <- options(contrasts=contr)
on.exit(options(old.con))

distr <- match.arg(distr)
tpm   <- x
if(is.null(ispd)) ispd <- reviseIspd(tpm)

if(is.null(data)) {
    if(is.null(nobs)) {
        stop("If \"data\" is not supplied, then \"nobs\" must be.\n")
    }
    data <- data.frame(cf=factor(rep(1,nobs)))
} else {
# Add the "cells" factor to "data".
    if(is.null(cells)) {
        data[,"cf"] <- factor(rep(1,nrow(data)))
    } else {
        if(!all(cells %in% names(data))) {
            stop("Some of the \"cells\" names are not in names(data).\n")
        }
        data[,"cf"] <- interaction(data[,cells])
    }
# Get the cells in the right fucking order!!!
    sdata <- split(data,f=data$cf)
    data  <- do.call(rbind,sdata)
}

# Build the replicated data frame.
K    <- nrow(tpm)
M    <- nrow(data)
data <- as.data.frame(lapply(data,function(x,K,M){rep(x,rep(K,M))},K=K,M=M))
data$state <- factor(rep(1:K,M))

simMlt(formula,response,distr,data,ispd,tpm,phi,Rho,sigma,size,ntop,zeta,missFrac,fep)
}
