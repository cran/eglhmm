reviseModel <- function(formula,data,distr,preSpecSigma,size) {
#
# Function reviseModel.  Calculates updated value of the
# model parameters, given a data frame which includes
# weights (== the gamma's) for each state, at each time.
#
# Note that the third term of the "auxiliary function"
# Q(theta,theta'), which is being maximised, is the weighted log
# likelihood of the model.  Except when distr is "Multinom" this
# is a generalised linear model.  This generalised linear model may
# be fitted thereby using glm(), with weights equal to the gamma's.
# The resulting coefficients are the updated estimate of phi.
# When distr is "Multinom" things work differently.
#
ynm <- as.character(formula)[2]
switch(EXPR=distr,
    Gaussian = {
        rslt <- lm(formula, data=data, weights=weights,na.action=na.exclude)
        phi <- coef(rslt)
        nms <- names(phi)
        nms[nms=="(Intercept)"] <- "mu"
        names(phi) <- nms
        gmu   <- fitted(rslt)
        if(is.null(preSpecSigma)) {
            rrr   <- resid(rslt)
            sigma <- reviseSigma(rrr,data$weights,data$state)
            if(any(sigma < sqrt(.Machine$double.eps))) {
                whinge <- paste0("Some entries of \"sigma\" appear to have converged\n",
                                 "  to zero.  This probably indicates that a Gaussian\n",
                                 "  model is inappropriate for your data.\n")
                stop(whinge)
            }
        } else {
            sigma <- preSpecSigma
        }
        names(sigma) <- paste0("sigma",levels(data[["state"]]))
        sd    <- sigma[data$state]
        fy    <- dnorm(data[,ynm],mean=gmu,sd=sd)
        modSpecs <- list(phi=phi,sigma=sigma,gmu=gmu,sd=sd)
    },
    Poisson = {
        rslt <- glm(formula, data=data, weights=weights, family=poisson,
                    na.action=na.exclude)
        phi <- coef(rslt)
        nms <- names(phi)
        nms[nms=="(Intercept)"] <- "mu"
        names(phi) <- nms
        lam  <- predict(rslt,newdata=data,type="response")
        fy   <- dpois(data[,ynm],lam)
        modSpecs <- list(phi=phi,lambda=lam)
    },
    Binomial = {
        bfmla <- binForm(formula)
        rslt <- glm(bfmla, data=cbind(data,size=size), weights=weights,
                    family=binomial,na.action=na.exclude)
        phi <- coef(rslt)
        nms <- names(phi)
        nms[nms=="(Intercept)"] <- "mu"
        names(phi) <- nms
        p <- predict(rslt,newdata=data,type="response")
        fy <- dbinom(data[,ynm],size=size,prob=p)
        modSpecs <- list(phi=phi,p=p,size=size)
    },
    Multinom = {
        Rho     <- reviseRho(data,response=ynm,formula,type=1)
        fy      <- ffun(data,formula,response=NULL,Rho=Rho,type=1)
        phi     <- rho2Phi(Rho)
        modSpecs <- list(phi=phi,Rho=Rho) # Of course phi is redundant,
                                          # but this makes the syntax
                                          # more consistent.
    }

)
fy[is.na(fy)] <- 1	
fy  <- split(fy,f=data[["cf"]])
nms <- unique(data$cf)
fy  <- fy[nms] # Keeps "weights" properly aligned with observations.
modSpecs$fy <- fy
modSpecs
}
