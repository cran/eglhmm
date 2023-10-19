getHgl <- function(nd,distr,theta,data,fmla,size,nbot,ntop) {
#
# Function getHgl --- get Hessian, gradient, log likelihood.
#
# Note: nd == number of derivatives to calculate.
# 0 <--> just calculate log likelihood (previously done by getLl())
# 1 <--> calculate gradient and log likelihood (previously done by getGl())
# 2 <--> calculate Hessian, gradient and log likelihood
#
dists <-c("Gaussian","Poisson","Binomial","Dbd","Multinom")
if(!(distr %in% dists)) {
    whinge <- paste0("Distribution ",distr," not recognised.\n")
    stop(whinge)
}
ndistr <- match(distr,dists)
K      <- length(levels(data[["state"]]))
if(is.null(size)) size <- 1

# Determine whether tau is to be revised/estimated in the LM step.
inclTau <- attr(theta,"inclTau")

# Dig out tau, zeta and phi from theta.
tau  <- theta$tau
zeta <- theta$zeta
phi  <- theta$phi
npt  <- length(tau) + length(zeta) + length(phi)
npar <- if(inclTau) npt else npt - length(tau) 
nms  <- c(if(inclTau) names(tau) else  NULL,names(zeta),names(phi))
#theta$tau <- reviseTau(distr,fmla,data,theta,size,nbot,ntop)

ynm  <- as.character(fmla[2])
y    <- data[[ynm]]
X    <- model.matrix(fmla[-2],data=data)
nxc  <- ncol(X)
if(distr %in% c("Gaussian","Poisson","Binomial")) {
    if(length(phi) != nxc) {
        stop("Alles upgefucken ist!\n")
    }
}

# Check that y is a factor when distr is "Multinom" (and if so,
# convert it to a numeric vector) and is numeric otherwise.  In the
# "Multinom" case, we must dig out the levels of y before doing the
# conversion to numeric mode, i.e. before buggering up y.
if(distr=="Multinom") {
    if(inherits(y,"factor")) {
        nyv   <- length(levels(y))
        y     <- as.numeric(y)
    } else {
        stop("For the Multinom distribution, \"y\" must be a factor.\n")
    }
} else {
    if(!is.numeric(y)) {
        whinge <- paste0("For the Gaussian, Poisson, Binomial and Dbd\n",
                         "  distributions, \"y\" must be numeric.\n")
        stop(whinge)
    }
    nyv <- 1
}

# Get the model components.
mc <- getModComps(distr,fmla,data,theta,size,nbot,ntop)
fy <- mc$fy

# The following components will be dummies/space-holders except when
# "distr" is the appropriate one.
gmu    <- mc$gmu    # Dummy except when distr is Gaussian.
sd     <- mc$sd     # Dummy except when distr is Gaussian.
lambda <- mc$lambda # Dummy except when distr is Poisson.
p      <- mc$p      # Dummy except when distr is Binomial.
ashp   <- mc$ashp   # Dummy except when distr is Dbd.
bshp   <- mc$bshp   # Dummy except when distr is Dbd.
Rho    <- mc$Rho    # Dummy except when distr is Multinom.


# Get tpm, and ispd
tpm  <- getTpm(theta,K)
ispd <- reviseIspd(tpm)

d12p  <- derivp(npt,K,tau,TRUE)
d1p   <- d12p$d1p
d2p   <- d12p$d2p
d12pi <- derivpi(ispd,tpm,npt,d12p)
d1pi  <- d12pi$d1pi
d2pi  <- d12pi$d2pi
nc    <- length(fy)
lcf   <- names(fy)

xll        <- numeric(nc)
names(xll) <- lcf
if(nd >= 1) {
    alist <- vector("list",nc)
    names(alist) <- lcf
}
if(nd == 2) {
    hlist        <- vector("list",nc)
    names(hlist) <- lcf
}

# Run through the cells:
for(xl in lcf) {
    ind   <- (data[["cf"]]==xl)
    Xxl   <- X[ind,,drop=FALSE]
    nind  <- nrow(Xxl)
    fyxl  <- fy[[xl]]
    yxl   <- matrix(y[ind],ncol=K,byrow=TRUE)[,1] # All columns are identical.
    ymiss <- is.na(yxl)
    ny    <- length(yxl)
    xxx   <- .Fortran(
        "gethgl",
        NAOK=TRUE,
        nd=as.integer(nd),
        ndistr=as.integer(ndistr),
        fy=as.double(fyxl),
        gmu=as.double(gmu[ind]),
        sd=as.double(sd[ind]),
        lambda=as.double(lambda[ind]),
        p=as.double(p[ind]),
        ashp=as.double(ashp[ind]),
        bshp=as.double(bshp[ind]),
        phimat=as.double(Rho),
        y=as.double(yxl),
        ymiss=as.integer(ymiss),
        tdm=as.double(t(Xxl)),
        nind=as.integer(nind),
        n=as.integer(ny),
        tpm=as.double(tpm),
        xispd=as.double(ispd),
        size=as.integer(size),
        nbot=as.integer(nbot),
        ntop=as.integer(ntop),
        d1pi=as.double(d1pi),
        d2pi=as.double(d2pi),
        kstate=as.integer(K),
        npar=as.integer(npar),
        npt=as.integer(npt),
        nyv=as.integer(nyv),
        nxc=as.integer(nxc),
        d1p=as.double(d1p),
        d2p=as.double(d2p),
        d1f=double(K*npar),
        d2f=double(K*npar*npar),
        d1a=double(K),
        d1b=double(K),
        d2aa=double(K),
        d2ab=double(K),
        d2bb=double(K),
        alpha=double(K),
        alphw=double(K),
        a=double(K*npar),
        b=double(K*npar*npar),
        aw=double(K*npar),
        bw=double(K*npar*npar),
        xlc=double(ny),
        hess=double(npar*npar),
        PACKAGE="eglhmm"
        )
    xlc <- xxx$xlc
    if(nd >= 1) alist[[xl]]  <- matrix(xxx$a,K,npar)/xlc[ny]
    xll[xl] <- sum(log(xlc))
    if(nd == 2) hlist[[xl]]  <- matrix(xxx$hess,npar,npar)
} 
ll <- sum(xll)
if(nd == 0) return(list(ll=ll))
a           <- array(unlist(alist),c(K,npar,nc))
grad        <- apply(a,2,sum)
names(grad) <- nms
if(nd == 1) return(list(ll=ll,gradient=grad))
hpre           <- array(unlist(hlist),c(npar,npar,nc))
Hessian        <- apply(hpre,c(1,2),sum)
dimnames(Hessian) <- list(nms,nms)
list(ll=ll,gradient=grad,Hessian=Hessian)
}
