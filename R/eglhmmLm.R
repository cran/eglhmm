eglhmmLm <- function(formula,data,distr,inclTau,preSpecSigma,
                     size,nbot,ntop,tau,zeta,phi,lmc=10,mixture=FALSE,
                     itmax=200,crit,tolerance=NULL,digits=NULL,verbose=FALSE) {
#
# Function eglhmmLm.  Fit an extended generalised hidden Markov
# model using the Levenberg-Marquardt algorithm.
#

if(mixture)
    stop(paste("Sorry; can't do mixtures by the",
               "Levenberg-Marquardt algorithm yet.\n"))

# Pick out the index of the stopping criterion:
icrit <- match(crit, c("CLL", "L2", "Linf","ABSGRD"))
if(is.na(icrit)) stop('Stopping criterion not recognized.\n')

# Set the tolerance if it is NULL.
if(is.null(tolerance)) tolerance <- sqrt(.Machine$double.eps)

# Set the number of digits with which to print out
# "progress reports".
if(is.null(digits)) digits <- 2+ceiling(abs(log10(tolerance)))

# Do some initial housekeeping.
fmla  <- formula
K     <- length(levels(data$state))
npro  <- K*(K-1)
X     <- model.matrix(fmla[-2],data)
nmst  <- if(inclTau) names(tau) else NULL
nmsz  <- if(is.null(preSpecSigma)) names(zeta) else NULL
nmsp  <- names(phi)
nms   <- c(nmst,nmsz,nmsp)
y     <- data[,as.character(fmla)[2]]
cf    <- data$cf
state <- data$state
theta <- list(tau=tau,zeta=zeta,phi=phi)
attr(theta,"inclTau")      <- inclTau
attr(theta,"preSpecSigma") <- preSpecSigma
oldll   <- getHgl(nd=0,distr,theta,data,fmla,size,nbot,ntop)$ll
oldtvec <- unlist(theta)
fmt     <- paste0("%.",digits,"f")
scrit   <- numeric(4)

# Ready to go.
if(verbose){
    cat("\n      Initial set-up completed ...\n")
    cat("\n      Initial log-likelihood: ",
    sprintf(fmt,oldll),"\n\n",sep="")
}

# Update:
if(verbose) cat('Repeating ...\n\n')
lm.step <- 1

th.new <- theta
repeat { # Swear again; i.e. recurse.
    if(verbose) cat(paste('Lev.-Marq. step ',lm.step,':\n',sep=''))
    xxx <- lmstep(th.new,data,fmla,distr,size,nbot,ntop,lmc)
    th.new  <- xxx$theta
    tpm <- getTpm(th.new,K)
    if(identical(all.equal(tpm,diag(K)),TRUE)) {
        stop("The \"tpm\" estimate has converged to the identity.\n")
    }
    tvec <- unlist(th.new)
    ll   <- xxx$ll
    grad <- xxx$gradient
    scrit[1] <- if(oldll > -Inf)
                        (ll - oldll)/(abs(ll) + tolerance)
                    else Inf
    scrit[2] <- sqrt(sum((tvec-oldtvec)^2))/(sqrt(sum(tvec^2)) + tolerance)
    scrit[3] <- max(abs(tvec-oldtvec))/(max(abs(tvec)) + tolerance)
    scrit[4] <- max(abs(grad))/(1+abs(ll))
    if(verbose){
            cat('     Log-likelihood: ',
            sprintf(fmt,ll),"\n",sep="")
            cat('     Scaled increase in log-likelihood: ',
            sprintf(fmt,scrit[1]),"\n",sep="")
            cat('     Scaled root-SS of change in coef.: ',
            sprintf(fmt,scrit[2]),"\n",sep="")
            cat('     Scaled max. abs. change in coef.: ',
            sprintf(fmt,scrit[3]),"\n",sep="")
            cat('     Scaled max. abs. val. of grad. vector.: ',
            sprintf(fmt,scrit[4]),"\n",sep="")
            cat('     Used steepest ascent? ',
                      attr(xxx,"usedSteep"),"\n",sep="")
    }

    if(scrit[icrit] < tolerance) {
        converged <- TRUE
        nstep     <- lm.step
        break
    }
    if(lm.step >= itmax) {
        cat(paste("Failed to converge in ",itmax,
                  " Levenberg-Marquardt steps.\n",sep=""))
        converged <- FALSE
        nstep     <- lm.step
        break
    }
    oldtvec <- tvec
    oldll   <- ll
    lm.step <- lm.step + 1
    lmc     <- 0.1*xxx$lmc
}

theta <- xxx$theta
xxx   <- getHgl(nd=2,distr,theta,data,fmla,size,nbot,ntop)
Hess  <- xxx$Hessian
grad  <- xxx$grad
tpm   <- getTpm(theta,K)
ispd  <- reviseIspd(tpm)
phi   <- theta$phi
switch(EXPR=distr,
    Gaussian = {
        gmu <- X%*%phi
        if(is.null(preSpecSigma)) {
            zeta  <- theta$zeta
            sigma <- exp(zeta)
        } else {
            sigma <- preSpecSigma
        }
        names(sigma) <- paste0("sigma",1:K)
        sd <- sigma[state]
        fy <- dnorm(y,mean=gmu,sd=sd)
    },
    Poisson = {
        lambda <- exp(X%*%phi)
        fy     <- dpois(y,lambda=lambda)
    },
    Binomial = {
        p  <- logistic(X%*%phi)
        fy <- dbinom(y,prob=p,size=size)
    },
    Dbd = {
        np <- length(phi)
        if(np%%2 != 0)
            stop("There must be an even number of phi coefficients.\n")
        kp     <- np/2
        phia   <- phi[1:kp]
        phib   <- phi[(kp+1):np]
        alpha  <- X%*%phia
        beta   <- X%*%phib
        fy     <- dbd::ddb(y,alpha,beta,ntop,nbot==0)
    },
    Multinom = {
        Rho <- phi2Rho(phi,K,rhovals=levels(y),preds=colnames(X))
        Q   <- X%*%Rho
        P   <- apply(Q,1,expForm2p)
        fy  <- P[cbind(y,1:nrow(X))]
    }
)
fy[is.na(fy)] <- 1	
fy    <- split(fy,f=data[["cf"]])
nmsfy <- unique(data$cf)
fy    <- fy[nmsfy] # Keeps "weights" properly aligned with observations.
rownames(Hess) <- colnames(Hess) <- nms
names(grad)    <- nms

part1 <- list(tpm=tpm,ispd=ispd,phi=phi,theta=theta,mixture=FALSE,
              Hessian=Hess,gradient=grad,log.like=ll,crit=crit,
              tolerance=tolerance,converged=converged,nstep=lm.step,formula=fmla,
              distr=distr,stopCritVal=scrit[icrit],fy=fy)
part2 <- switch(EXPR=distr,
    Gaussian = {
        mu <- getMu(gmu,data,fmla)
        if(is.null(preSpecSigma)) {
            list(mean=gmu,sd=sd,mu=mu,sigma=sigma)
        } else {
            list(mean=gmu,sd=sd,mu=mu,preSpecSigma=sigma)
        }
    },
    Poisson = {
        list(lambda=lambda)
    },
    Binomial = {
        list(p=p,size=size)
    },
    Dbd = {
        list(alpha=alpha,beta=beta,nbot=nbot,ntop=ntop)
    },
    Multinom = {
        list(Rho=Rho)
    }
)
c(part1,part2)
}
