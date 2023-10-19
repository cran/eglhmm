reglhmm.eglhmm <-function(x,missFrac=NULL,...) {
#
# The "eglhmm" method for reglhmm().  Simulates data from a fitted
# hidden generalised linear Markov model.
#

# Keep the kontrast konsistent. (Note that simMlt() uses
# the contrast in forming the model matrix.)
contrast <- x$contrast
contr    <- c(paste("contr",contrast,sep="."),"contr.poly")
oldcon   <- options(contrasts=contr)
on.exit(options(oldcon))

ispd  <- x$ispd
tpm   <- x$tpm
phi   <- x$phi
Rho   <- x$Rho
sigma <- x$sigma
fmla  <- x$formula
resp  <- x$response
data  <- x$data
distr <- x$distr
size  <- x$size
ntop  <- x$ntop
nbot  <- x$nbot
if(distr=="Dbd") {
    zeta  <- nbot==0
} else zeta <- NULL

if(is.null(missFrac)) missFrac <- x$missFrac
fep <- list(...)[["fep"]]

simMlt(fmla,resp,distr,data,ispd,tpm,phi,Rho,sigma,size,ntop,zeta,missFrac,fep)
}
