fitted.eglhmm <- function(object,...) {
	fy      <- object$fy
        nms     <- unique(object$data$cf)
        fy      <- fy[nms] # To align gamma and the means correctly.
	tpm     <- object$tpm
        distr   <- object$distr
	gamma   <- recurse(fy,tpm)$gamma
        parnm   <- switch(EXPR=distr,Gaussian="mean",Poisson="lambda",
                                     Binomial="p",Dbd="mean")
	fit     <- tapply(object[[parnm]],object$data$state,function(x){x})
	fit     <- do.call("cbind",fit)
	yhat    <- apply(t(gamma)*fit,1,sum)
	data    <- object$data
	data    <- data[data$state==levels(data$state)[1],]
        jjj     <- match(c("cf","state"),names(data))
	data    <- data[,-jjj]
	attr(yhat,"data") <- data
	class(yhat) <- "fitted.eglhmm"
	yhat
}
