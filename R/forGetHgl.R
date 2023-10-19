forGetHgl <- function(nd,theta,data,fmla,inclTau=TRUE,preSpecSigma=NULL) {
#
# Function forGetHgl --- get/revise the Hessian, gradient, log
# likelihood, using the Fortran code.  Gaussian distribution only.
#
# Note: nd == number of derivatives to calculate.
# 0 <--> just calculate log likelihood (previously done by getLl())
# 1 <--> calculate gradient and log likelihood (previously done by getGl())
# 2 <--> calculate Hessian, gradient and log likelihood
#

requireNamespace("eglhmm")

# Adjust the attributes (and possibly the zeta component) of theta.
if(!missing(inclTau)) {
    attr(theta,"inclTau") <- inclTau
}
if(!missing(preSpecSigma)) {
    attr(theta,"preSpecSigma") <- preSpecSigma
    if(!is.null(preSpecSigma)) theta$zeta <- NULL
}

xxx  <- getHgl(nd,distr="Gaussian",theta,data,fmla,size=NULL,nbot=NULL,ntop=NULL)
xxx
}
