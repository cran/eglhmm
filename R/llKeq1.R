llKeq1 <- function(data,object) {
    bivar <- length(object$response) == 2
    if(bivar) {
        type <- if(object$indep) 2 else 3
        fy   <- ffun(data,fmla=NULL,response=object$response,Rho=Rho,type=type)
    } else {
        ynm <- as.character(object$formula[2])
        y   <- data[[ynm]]
        switch(EXPR=object$distr,
            Gaussian = {
                fy  <- dnorm(y,mean=object$mean,sd=object$sigma)
            },
            Poisson = {
                fy <- dpois(y,object$lambda)
            },
            Binomial = {
                fy <- dbinom(y,size=object$size,prob=object[["p"]])
            },
            Dbd = {
                fy    <- dbd::ddb(y,alpha=object$alpha,beta=object$beta,
                                  ntop=object$ntop,object$nbot==0)
            },
            Multinom = {
                data$state <- factor(1)
                fy   <- ffun(data,fmla=object$formula,response=NULL,Rho=object$Rho,type=1)
            }
        )
    }
    if(object$distr=="Multinom") {
        Rho <- object$Rho
        if(bivar) {
            if(object$indep) {
                npar <- sum(sapply(Rho,nrow))-2
            } else {
                npar <- prod(dim(Rho))-1
            }
        } else {
            npar <- ncol(Rho)-1
        }
    } else {
        npar <- length(object$phi)
    }
    fy[is.na(fy)] <- 1
    ll   <- sum(log(fy))
    AIC  <- -2*ll+2*npar
    BIC  <- -2*ll+(object$bicm)*npar
    rslt <- list(ll=ll)
    attr(rslt,"AIC") <- AIC
    attr(rslt,"BIC") <- BIC
    rslt
}
