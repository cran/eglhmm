bcov <-function(object,nsim=50,itmax=500,verbose=TRUE) {
# Function bcov (bootstrap covariance) to create a
# version of the covariance matrix of the parameters of a
# generalized linear hidden Markov model through  parametric bootstrapping.

missFrac <- object$missFrac
rslt     <- list()
i        <- 0
nc.count <- 0
emUsed   <- object$method == "em"
if(emUsed) an.count <- 0
repeat{
	dat.sim <- reglhmm(object,missFrac=missFrac)
        p <- try(update(object,data=dat.sim,checkDecrLL=FALSE,
                      par0=object,itmax=itmax,verbose=FALSE))
        if(inherits(p,"try-error")) browser()
	if(p$converged) {
		i <- i+1
		rslt[[i]] <- p$theta
                if(emUsed) {
                    if(p$anomaly) an.count <- an.count + 1
                }
		if(verbose) cat(i,"")
                if(verbose & i%%10 == 0) cat("\n")
	}
	else nc.count <- nc.count+1
	if(i>=nsim) break
}
if(verbose & i%%10 != 0) cat("\n")

# Put things together and return.
m <- matrix(unlist(rslt),byrow=TRUE,ncol=length(rslt[[1]]))
C_hat <- var(m)
nms   <- names(object$theta)
dimnames(C_hat) <- list(nms,nms)
rslt <- if(emUsed) {
    list(C_hat=C_hat,nc.count=nc.count,an.count=an.count)
} else {
    list(C_hat=C_hat,nc.count=nc.count)
}
rslt
}
