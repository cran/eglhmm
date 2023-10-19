eglhmmBf <- function(formula,data,distr,theta,size,nbot,ntop,
                     optimiser,optimMethod,nlmWarn,hessian,useAnalGrad,ca,
                     itmax,tolerance,verbose) {
# Function eglhmmBf.  Fit an extended generalised hidden Markov
# model, using brute force.

# Dig out the crucial attributes.
inclTau      <- attr(theta,"inclTau")
preSpecSigma <- attr(theta,"preSpecSigma")

# Do some initial housekeeping.
fmla <- formula
K    <- length(levels(data$state))
npro <- K*(K-1)
tau  <- if(inclTau) theta$tau else NULL
zeta <- if(is.null(preSpecSigma)) theta$zeta else NULL
phi  <- theta$phi
nll  <- function(tvec,data,fmla,distr,theta,size,nbot,ntop,useAnalGrad) {
            theta <- reviseTheta(tvec,theta,distr,fmla,data,size,nbot,ntop)
            if(useAnalGrad) {
                xxx <- getHgl(nd=1,distr,theta,data,fmla,size,nbot,ntop) 
                rslt <- -xxx$ll
                attr(rslt,"gradient") <- -xxx$gradient
            } else {
               rslt <- -getHgl(nd=0,distr,theta,data,fmla,size,nbot,ntop)$ll
            }
         rslt
}

# Do it to it.
tvec <- c(tau,zeta,phi)
attr(tvec,"inclTau")      <- inclTau
attr(tvec,"preSpecSigma") <- preSpecSigma
if(optimiser=="optim") {
    if(useAnalGrad) {
        gradFun <- function(theta,data,fmla,distr,size,nbot,ntop,...) {
                       -getHgl(nd=1,distr,theta,data,fmla,size,nbot,ntop)$gradient
                   }
    } else gradFun <- NULL

# Set the control.
    ctrl <- if(verbose) list(trace=6,maxit=itmax) else list(maxit=itmax)
    if(!is.null(tolerance)) ctrl <- c(ctrl,list(reltol=tolerance))
        rslt   <- optim(tvec,nll,gr=gradFun,method=optimMethod,
                        hessian=hessian,control=ctrl,
                        fmla=fmla,theta=theta,data=data,distr=distr,size=size,
                        nbot=nbot,ntop=ntop,useAnalGrad=FALSE)
} else {
    PL      <- if(verbose) 2 else 0
    savenms <- names(tvec) # Check that this works!!!
    warnFn  <- if(nlmWarn) function(x){x} else suppressWarnings
    rslt    <- warnFn(nlm(f=nll,p=tvec,hessian=hessian,print.level=PL,iterlim=itmax,
                          data=data,fmla=fmla,distr=distr,theta=theta,size=size,
                          nbot=nbot,ntop=ntop,useAnalGrad=useAnalGrad,
                          check.analyticals=ca))
}

if(optimiser=="optim") {
    tvec  <- rslt$par
    ll    <- -rslt$value
    kount <- rslt$counts
    conv  <- rslt$convergence==0
} else {
    tvec  <- rslt$estimate
    ll    <- -rslt$minimum
    kount <- rslt$iterations
    conv  <- rslt$code==1
    names(tvec)       <- savenms
    attr(conv,"code") <- rslt$code
}
if(hessian) numHess <- -rslt$hessian

# Create the updated theta list.
theta <- reviseTheta(tvec,theta,distr,fmla,data,size,nbot,ntop)
zeta  <- theta$zeta
phi   <- theta$phi
tpm   <- getTpm(theta,K)
ispd  <- reviseIspd(tpm)
X     <- model.matrix(fmla[-2],data)
y     <- data[,as.character(fmla)[2]]
state <- data[,"state"]
switch(EXPR=distr,
    Gaussian = {
        mean <- try(X%*%phi)
        if(inherits(mean,"try-error")) browser()
        mu   <- getMu(mean,data,fmla)
        if(is.null(preSpecSigma)) {
            sigma <- exp(zeta)
        } else {
            sigma <- preSpecSigma
        }
        names(sigma) <- paste0("sigma",1:K)
        sd    <- sigma[state]
        fy    <- dnorm(y,mean=mean,sd=sd)
    },
    Poisson = {
        lam <- exp(X%*%phi)
        fy  <- dpois(y,lam)
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
fy   <- split(fy,f=data[["cf"]])
nms  <- unique(data$cf)
fy   <- fy[nms] # Keeps "weights" properly aligned with observations.
xxx  <- getHgl(nd=2,distr,theta,data,fmla,size,nbot,ntop)
grad <- xxx$gradient
names(grad) <- names(tvec)
if(hessian) rownames(numHess) <- colnames(numHess) <- names(tvec)

part1 <- list(tpm=tpm,ispd=ispd,phi=phi,theta=theta,
              gradient=grad,log.like=ll,converged=conv,
              nstep=kount,formula=fmla,distr=distr,fy=fy)
part2 <- switch(EXPR=distr,
    Gaussian = {
        if(is.null(preSpecSigma)) {
            list(mean=mean,sd=sd,mu=mu,sigma=sigma)
        } else {
            list(mean=mean,sd=sd,mu=mu,preSpecSigma=sigma)
        }
    },
    Poisson = {
        list(lambda=lam)
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
part3 <- if(hessian) list(numHess=numHess) else NULL
c(part1,part2,part3)
}
