#
# Functions that are used to convert parameter structures from
# one form to another.
#

p2expForm <- function (x) {
# Convert a vector of probabilities (summing to 1) to
# a vector of parameters for the logistic style
# parametrisation of these probabilities.  Notice that
# the last entry of this vector is constrained to be 0.
# There has to be a constraint of course, corresponding
# to the original constraint that the probabilities
# sum to 1.
    z  <- log(x)
    nok <- z==-Inf
    if(any(nok)) {
        btm <- min(z[!nok])
        z[nok] <- min(btm,-300)
    }
    z - z[length(z)]
}

expForm2p <- function(x){
# Convert a vector of parameters for the logistic or exponential style
# parameterisation to a vector of probabilities summing to 1.
    m  <- max(x)
    xr <- exp(x-m)
    xr/sum(xr)
}

getIspd <- function(theta,K) {
    return(expForm2p(c(theta[1:(K-1)],0)))
}

getTpm <- function(theta,K) {
# Note that theta is now (30/07/2023) a *list* with components
# tau, zeta and phi --- and includes possibly redundant parameters
# --- rather than being a shaganappi vector of (already updated)
# non-redundant parameters.  Get the transition probability matrix
# from the tau component of theta.
    if(K > 1) {
        rawmat <- cbind(matrix(theta$tau,nrow=K),0)
        tpm    <- t(apply(rawmat,1,expForm2p))
    } else tpm <- NA
    tpm
}

rho2Phi <- function(Rho) {
    phi     <- t(Rho)[-ncol(Rho),,drop=FALSE]
    nms     <- if(nrow(Rho) > 1) {
                   outer(rownames(phi),colnames(phi),
                         function(a,b){paste(a,b,sep=".")})
               } else {
                   colnames(Rho)[-ncol(Rho)]
               }
    phi     <- as.vector(phi)
    names(phi) <- as.vector(nms)
    phi
}

phi2Rho <- function(phi,K,rhovals,preds) {
    m      <- length(rhovals)
    Rho    <- matrix(phi,ncol=m-1,byrow=TRUE)
    if(nrow(Rho) != length(preds)) {
        stop("Wrong length for \"phi\".\n")
    }
    Rho <- cbind(Rho,0)
    rownames(Rho) <- preds
    colnames(Rho) <- rhovals
    Rho
}

cnvrtRho <- function(Rho) {
    if(inherits(Rho,"RhoExpForm")) {
#
# Rho is a K x m matrix (K = number of states, m = number of possible
# y values), with its last --- mth --- column identically zero.
# Columns correspond to y values, rows to states.
# Pr(Y = y_j | state = i) = exp(Rho[i,j])/sum_ell(exp(Rho[ell,j])).
#
        ok <- isTRUE(all.equal(Rho[,ncol(Rho)],rep(0,nrow(Rho)),
                               check.attributes=FALSE))
        if(!ok)
            stop("The last column of \"Rho\" must be all zeroes.\n")
        newRho <- apply(Rho,1,expForm2p)
        class(newRho) <- c(class(newRho),"RhoProbForm")
        return(newRho)
    }
    if(inherits(Rho,"RhoProbForm")) {
#
# Rho is an m x K matrix with with non-negative entries,
# all columns summing to 1.
# Columns correspond to states, rows to y values.
# Pr(Y = y_i | state = j) = Rho[i,j].
#
        ok <- all.equal(apply(Rho,2,sum),rep(1,ncol(Rho)),
                        check.attributes=FALSE)
        if(!ok)
            stop("The columns of \"Rho\" must be all sum to 1.\n")
        newRho <- t(apply(Rho,2,p2expForm))
        class(newRho) <- c(class(newRho),"RhoExpForm")
        return(newRho)
    }
    stop("Argument \"Rho\" is not of an appropriate class.\n")
}    
