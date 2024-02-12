doKeq1 <- local({

    nll  <- function(phi,X,y,nbot,ntop) {
                np   <- length(phi)
                kp   <- np/2
                phia   <- phi[1:kp]
                phib   <- phi[(kp+1):np]
                ashp   <- X%*%phia
                bshp   <- X%*%phib
                rslt <- -dbd::ddb(y,ashp,bshp,ntop,nbot==0,log=TRUE)
                rslt[is.na(rslt)] <- 0
                sum(rslt)
    }

function(data,fmla,distr,preSpecSigma,response,indep,size,nbot,ntop,bicm,nafrac) {
# When K=1 fit an i.i.d. model.
    rsplvls <- attr(data,"rsplvls")
    bivar <- length(response)==2
    if(bivar) {
        y1 <- response[1]
        y2 <- response[2]
        if(indep) {
            X  <- factor(data[[y1]],levels=rsplvls[[1]])
            Y  <- factor(data[[y2]],levels=rsplvls[[2]])
            Rho <- vector("list",2)
            Rho[[1]] <- t(as.matrix(table(X)))
            Rho[[1]] <- Rho[[1]]/sum(Rho[[1]])
            Rho[[2]] <- t(as.matrix(table(Y)))
            Rho[[2]] <- Rho[[2]]/sum(Rho[[2]])
            fy  <- ffun(data,response=response,Rho=Rho,type=2)
        } else {
            X  <- factor(data[[y1]],levels=c(rsplvls[[1]],NA),exclude=NULL)
            Y  <- factor(data[[y2]],levels=c(rsplvls[[2]],NA),exclude=NULL)
            G    <- table(X,Y,useNA="always")
            m    <- nrow(G)
            n    <- ncol(G)
            Rho0 <- G[-m,-n]
            Rho0 <- Rho0/sum(Rho0)
            Rho0 <- array(Rho0,dim=c(dim(Rho0),1))
            dnms <- c(rsplvls,list("1"))
            dimnames(Rho0) <- dnms
            G    <- array(G,dim=c(dim(G),1))
            Rho  <- msRho(Rho0,G)
            dimnames(Rho) <- dnms
            fy   <- ffun(data,fmla=NULL,response=response,Rho=Rho,type=3)
        }
        rslt <- list(Rho=Rho)
    } else {
        ynm  <- as.character(fmla[2])
        y    <- data[[ynm]]
        if(distr=="Binomial") {
            fmla <- binForm(fmla)
            data <- cbind(data,size=size)
        }
        if(distr %in% c("Gaussian","Poisson","Binomial")) {
            fam   <- switch(EXPR=distr,Gaussian="gaussian",Poisson="poisson",
                            Binomial="binomial")
            fit   <- glm(fmla,data,family=fam,na.action=na.exclude)
            phi   <- coef(fit)
        }
        switch(EXPR=distr,
            Gaussian = {
                gmu <- as.matrix(fitted(fit))
                if(is.null(preSpecSigma)) {
                    sigma <- sqrt(summary(fit)[["deviance"]]/summary(fit)[["df"]][2])
                } else {
                    sigma <- preSpecSigma
                }
                fy  <- dnorm(y,mean=gmu,sd=sigma)
                mu  <- getMu(gmu,data,fmla)
                rslt <- list(mu=mu,sigma=sigma,mean=gmu,sd=rep(sigma,length(gmu)),
                             phi=phi)
            },
            Poisson = {
                lambda <- fitted(fit)
                fy <- dpois(y,lambda)
                rslt <- list(lambda=lambda,phi=phi)
            },
            Binomial = {
                p  <- fitted(fit)
                fy <- dbinom(y,size=size,prob=p)
                rslt <- list(p=p,phi=phi)
            },
            Dbd = {
                X     <- model.matrix(fmla[-2],data=data)
                kp    <- ncol(X)
                np    <- 2*kp
                phi0  <- rep(0,np)
                fit   <- optim(phi0,nll,method="BFGS",X=X,y=y,nbot=nbot,ntop=ntop)
                phi   <- fit$par
                alpha <- X%*%phi[1:kp]
                beta  <- X%*%phi[(kp+1):np]
                fy    <- dbd::ddb(y,alpha,beta,ntop,nbot==0)
                rslt  <- list(alpha=alpha,phi=phi)
            },
            Multinom = {
                data$state   <- factor(1)
                data$weights <- 1
                Rho  <- reviseRho(data,response,fmla,type=1)
                phi  <- rho2Phi(Rho)
                fy   <- ffun(data,fmla=fmla,response=NULL,Rho=Rho,type=1)
                rslt <- list(Rho=Rho,phi=phi)
            }
        )
    }
    if(distr=="Multinom") {
        if(bivar) {
            if(indep) {
                npar <- sum(sapply(Rho,nrow))-2
            } else {
                npar <- prod(dim(Rho))-1
            }
        } else {
            npar <- ncol(Rho)-1
        }
    } else {
        npar <- length(phi)
    }
    fy[is.na(fy)] <- 1
    ll   <- sum(log(fy))
    AIC  <- -2*ll+2*npar
    BIC  <- -2*ll+bicm*npar
    theta <- list(tau=NULL,zeta=if(distr=="G") log(sigma) else NULL,phi=phi)
    attr(theta,"inclTau") <- FALSE
    rslt <- c(rslt,list(theta=theta,log.like=ll,fy=fy,
                        formula=fmla,distr=distr,bicm=bicm,npar=npar,
                        AIC=AIC,BIC=BIC,data=data,missFrac=nafrac))

    rslt <- c(rslt,list(tpm=NA,ispd=NA,converged=NA,nstep=NA,par0=NA,
                        tolerance=NA,crit=NA,mixture=NA,method=NA))
    rslt
}
})
