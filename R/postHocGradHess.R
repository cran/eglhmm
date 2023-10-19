postHocGradHess <- function(object,inclTau=TRUE) {
if(inherits(object,"eglhmm.bivariate")) {
    stop("Gradient and Hessian are not available for bivariate models.\n")
}
if(!inherits(object,"eglhmm")) {
    stop("Argument \"object\" must be of class \"eglhmm\".\n")
}
# Method = "em" --- no gradient or Hessian present in "object"
# Method = "lm" --- gradient and Hessian present in "object";
#                   these are calculated analytically.
# Method = "bf" --- gradient present in "object", calculated numerically;
#                   Hessian *may* be present; if so it is calculated numerically.

method <- object$method
revGH  <- method %in% c("em","bf") || attr(object$theta,"inclTau") != inclTau
if(revGH) {
    attr(object$theta,"inclTau") <- inclTau
    xxx <- with(object,getHgl(nd=2,distr=distr,theta,
                data=data,fmla=formula,size=object$size,
                nbot=object$nbot,ntop=object$ntop))
}
switch(EXPR=method,
    em = {
        grad <- xxx$grad
        hess <- xxx$Hess
        numGrad <- NULL
        numHess <- NULL
    },
    lm = {
        if(revGH) {
            grad <- xxx$grad
            hess <- xxx$Hess
        } else {
            grad <- object$grad
            hess <- object$Hess
        }
        numGrad <- NULL
        numHess <- NULL
    },
    bf = {
        grad <- xxx$grad
        hess <- xxx$Hess
        numGrad <- object$grad
        numHess <- object$numHess
    }
)
rslt <- list(gradient=grad,Hessian=hess)
if(!is.null(numGrad)) rslt <- c(rslt,list(numGrad=numGrad))
if(!is.null(numHess)) rslt <- c(rslt,list(numHess=numHess))
rslt
}
