getModComps <- function(distr,fmla,data,theta,size,nbot,ntop) {
X    <- model.matrix(fmla[-2],data=data)
K    <- length(levels(data$state))
ynm  <- as.character(fmla[2])
y    <- data[[ynm]]
lvls <- attr(data,"lvls")
nxr  <- nrow(X)
phi  <- theta$phi
zeta <- theta$zeta
cmplst <- list(gmu=numeric(nxr),sd=numeric(nxr),
               lambda=numeric(nxr),
               p=numeric(nxr),
               ashp=numeric(nxr),bshp=numeric(nxr))
switch(EXPR=distr,
    Gaussian = {
        gmu    <- X%*%phi
        preSpecSigma <- attr(theta,"preSpecSigma")
        if(is.null(preSpecSigma)) {
            sigma <- exp(zeta)
        } else {
            sigma <- preSpecSigma
        }
        names(sigma) <- paste0("sigma",1:K)
        sd           <- sigma[data[["state"]]]
        cmplst$gmu   <- gmu
        cmplst$sd    <- sd
        fy           <- dnorm(y,mean=gmu,sd=sd)
    },
    Poisson = {
        lambda        <- exp(X%*%phi)
	cmplst$lambda <- lambda
        fy            <- dpois(y,lambda)
    },
    Binomial = {
        p        <- logistic(X%*%phi)
	cmplst$p <- p
        fy       <- dbinom(y,size,p)
    },
    Dbd = {
        np <- length(phi)
        if(np%%2 != 0)
            stop("There must be an even number of phi coefficients.\n")
        kp     <- np/2
        if(kp != ncol(X)) stop("Alles upgefucken ist!\n")
        phia   <- phi[1:kp]
        phib   <- phi[(kp+1):np]
        ashp   <- X%*%phia
        bshp   <- X%*%phib
        cmplst$ashp <- ashp
        cmplst$bshp <- bshp
        fy     <- dbd::ddb(y,ashp,bshp,ntop,nbot==0)
    },
    Multinom = {
        Rho    <- phi2Rho(phi,K,rhovals=lvls,preds=colnames(X))
        cmplst <- list(Rho=Rho)
        fy     <- ffun(data,fmla=fmla,response=ynm,Rho=Rho,type=1)
    }
)

fy[is.na(fy)] <- 1
fy            <- split(fy,f=data[["cf"]])
nms           <- unique(data$cf)
fy            <- fy[nms] # Keeps "weights" properly aligned with observations.
cmplst$fy     <- fy
cmplst
}
