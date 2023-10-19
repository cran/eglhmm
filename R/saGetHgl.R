saGetHgl <- function(nd,theta,data,fmla,inclTau=TRUE,preSpecSigma=NULL) {
#
# Function saGetHgl --- get Hessian, gradient, log likelihood,
# using stand-alone R code; Gaussian distribution only.
#
# Note: nd == number of derivatives to calculate.
# 0 <--> just calculate log likelihood (previously done by getLl())
# 1 <--> calculate gradient and log likelihood (previously done by getGl())
# 2 <--> calculate Hessian, gradient and log likelihood
#

dyn.load("sag.so")
on.exit(dyn.unload("sag.so"))

# Adjust the attributes (and possibly the zeta component) of theta.
if(!missing(inclTau)) {
    attr(theta,"inclTau") <- inclTau
}
if(!missing(preSpecSigma)) {
    attr(theta,"preSpecSigma") <- preSpecSigma
    if(!is.null(preSpecSigma)) theta$zeta <- NULL
}

# Dig out K.
K <- length(levels(data$state))

# Determine whether tau is to be revised/estimated in the LM step.
inclTau <- attr(theta,"inclTau")

# Dig out tau, zeta and phi from theta.
tau  <- theta$tau
zeta <- theta$zeta
phi  <- theta$phi
npt  <- length(tau) + length(zeta) + length(phi)
npar <- if(inclTau) npt else npt - length(tau) 
nms  <- c(if(inclTau) names(tau) else  NULL,names(zeta),names(phi))

ynm  <- as.character(fmla[2])
y    <- data[[ynm]]
X    <- model.matrix(fmla[-2],data=data)
nxc  <- ncol(X)
if(length(phi) != nxc) {
    stop("Alles upgefucken ist!\n")
}

if(!is.numeric(y)) {
    whinge <- paste0("For the Gaussian distribution, \"y\" must be numeric.\n")
    stop(whinge)
}

# Get the model components.
mc <- getModComps("Gaussian",fmla,data,theta,size=NULL,nbot=NULL,ntop=NULL)
fy <- mc$fy

# Assign values to gmu and sd.
gmu    <- mc$gmu
sd     <- mc$sd


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
    ind  <- (data[["cf"]]==xl)
    Xxl  <- X[ind,,drop=FALSE]
    fyxl <- fy[[xl]]
    yxl  <- matrix(y[ind],ncol=K,byrow=TRUE)[,1] # All columns are identical.
    ny   <- length(yxl)
    xxx  <- saSubGetHgl(nd=nd,fy=fyxl,gmu=gmu[ind],sd=sd[ind],
                      y=yxl,tdm=t(Xxl),tpm=tpm,xispd=ispd,
                      d1pi=d1pi,d2pi=d2pi,kstate=K,npar=npar,npt=npt,
                      nxc=nxc,d1p=d1p,d2p=d2p,alpha=numeric(K),
                      alphw=numeric(K),a=numeric(K*npar),
                      b=numeric(K*npar*npar),aw=numeric(K*npar*npar),
                      bw=numeric(K*npar*npar),xlc=numeric(ny),
                      hess=numeric(npar*npar),xl)

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
Hess           <- apply(hpre,c(1,2),sum)
dimnames(Hess) <- list(nms,nms)
list(ll=ll,gradient=grad,Hessian=Hess)
}
