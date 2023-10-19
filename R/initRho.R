initRho <- local({

grnp <- function(K,randStart,lvls) {
# grnp <--> "generic Rho, no predictors" (other than state)
nval <- length(lvls)
if(randStart) {
    Rho <- matrix(runif(K*nval),K,nval)
} else {
    Rho <- matrix(1:(nval*K),K,nval)
}
Rho           <- t(Rho/apply(Rho,1,sum))
rownames(Rho) <- lvls
colnames(Rho) <- paste0("state",1:K)
class(Rho)    <- c(class(Rho),"RhoProbForm")
Rho           <- cnvrtRho(Rho)
Rho
}

function(data,K,fmla,randStart,indep,lvls) {
#
# To get the probabilities of the observations (given the predictor
# values) from this initial Rho or its updates, we form
# M <- X%*%Rho, where X is the model matrix, and then calculate
# P <- apply(M,1,expForm2p).  The matrix P is m x N, where m is the
# number of possible y values (= ncol(Rho)) and N = n*K is the "total"
# number of observations, each of the "actual" observations (of which
# there are n)  having been replicated K times (once for each state).
# Note that
#               P[i,j] = Pr(Y = y_i | the jth observation).
# To get the probablities of the observed values of y, given the
# predictors we form P[cbind(y,1:N)].
#
# In the future I might try using "fake states" to produce initial
# estimates of Rho.  However it seems probable that le jeu n'en
# vaut pas la chandelle.

# This should never happen, but just in case ....
if(K==1) return(NA)

# Now do Rho ....
nval   <- if(inherits(lvls,"list")) sapply(lvls,length) else length(lvls)
univar <- length(nval) == 1
if(univar) { # Univariate.
    dumDat <- cbind(data[1,],state=factor(1,levels=1:K))
    dumX   <- model.matrix(fmla,data=dumDat)
    ncX    <- ncol(dumX)
    if(ncX == K) {
        Rho <- grnp(K,randStart,lvls)
    } else {
        if(randStart) {
            Rho <- matrix(rnorm((nval-1)*ncX),nrow=ncX)
            Rho <- cbind(Rho,0)
        } else {
            Rho1 <- grnp(K,randStart,lvls) # K x nval
            Rho2 <- matrix(0,ncol=nval,nrow=ncX-K)
            Rho  <- rbind(Rho1,Rho2)
        }
        class(Rho)    <- "RhoExpForm"
        rownames(Rho) <- colnames(dumX)
        colnames(Rho) <- lvls
    }

} else { # Bivariate.
    if(indep) {
        Rho <- vector("list",2)
        for(i in 1:2) {
# This gives Rho in "exponential form".  Is this what we want?
# Yes it is; now!!!
            Rho[[i]] <- grnp(K,randStart,lvls[[i]])
        }
    } else {
        if(randStart) {
             Rho <- array(runif(K*prod(nval)),dim=c(nval[1],nval[2],K))
        } else {
             Rho <- array(1:(prod(nval)*K),dim=c(nval[1],nval[2],K))
        }
        div <- apply(Rho,3,sum)
        Rho <- aperm(Rho,c(3,1,2))/div
        Rho <- aperm(Rho,c(2,3,1))
        dimnames(Rho) <- c(lvls,list(1:K))
    }
}
Rho
}
})
