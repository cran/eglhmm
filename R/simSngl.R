simSngl <- function(distr,tpm,ispd,nt,mlp,Rho,yvals,datxl,fmla,response,
                          sig,size,ntop,zeta,mf,fep) {
# Function simSngl.  Simulates data from a single cell of a hidden
# Markov model.  Meaning of arguments:
# mlp  <--> linear predictor for the mean (mu, lambda, p).
# sig  <--> standard deviations (Gaussian only).
# ispd <--> initial state probability distribution.
# tpm  <--> transition probability matrix.
# mf   <--> missFrac == fraction of missing values.
# fep  <--> first emission present
#
    K   <- length(ispd)
    ynm <- as.character(fmla)[2]
    switch(EXPR=distr,
        Gaussian = {
            mu <- numeric(nt) 
            sd <- numeric(nt) 
            st <- sample(1:K,1,prob=ispd)
            mu[1] <- mlp[1,st]
            sd[1] <- sig[st]
            for(t in 2:nt) { 
                st <- sample(1:K,1,prob=tpm[st,])
                mu[t] <- mlp[t,st]
                sd[t] <- sig[st]
            }
            y <- data.frame(y=rnorm(nt,mean=mu,sd=sd))
            names(y) <- ynm
        },
        Poisson = {
            lambda <- numeric(nt)
            st     <- sample(1:K,1,prob=ispd)
            lambda[1] <- exp(mlp[1,st])
            for(t in 2:nt) {
                st <- sample(1:K,1,prob=tpm[st,])
                lambda[t] <- exp(mlp[t,st])
            }
            y <- data.frame(y=rpois(nt,lambda))
            names(y) <- ynm
        },
        Binomial = {
            p <- numeric(nt)
            st<- sample(1:K,1,prob=ispd)
            p[1] <- logistic(mlp[1,st])
            for(t in 2:nt) {
                st   <- sample(1:K,1,prob=tpm[st,])
                p[t] <- logistic(mlp[t,st])
            }
            y <- data.frame(y=rbinom(nt,size,p))
            names(y) <- ynm
        },
        Dbd = {
            mlpa  <- mlp$mlpa
            mlpb  <- mlp$mlpb
            alpha <- numeric(nt)
            beta  <- numeric(nt)
            st    <- sample(1:K,1,prob=ispd)
            alpha[1] <- mlpa[1,st]
            beta[1]  <- mlpb[1,st]
            for(t in 2:nt) {
                st <- sample(1:K,1,prob=tpm[st,])
                alpha[t] <- mlpa[t,st]
                beta[t]  <- mlpb[t,st]
            }
            if(requireNamespace("dbd")) {
                y <- data.frame(y=dbd::rdb(nt,alpha,beta,ntop,zeta))
                names(y) <- ynm
            } else {
                stop("Required package \"dbd\" is not available.\n")
            }
        },
# This function now handles univariate Multinom, bivariate
# independent and bivariate dependent data.

        Multinom = {
            sts <- numeric(nt)
            sts[1]  <- sample(1:K,1,prob=ispd)
            for(t in 2:nt) {
                sts[t]  <- sample(1:K,1,prob=tpm[sts[t-1],])
            }
            if(!inherits(yvals,"list")) { # Univariate.
                y   <- numeric(nt)
                nyv   <- length(yvals)
                ok    <- rep(sts,each=K) == rep(1:K,nt)
                datxl <- datxl[ok,]
                nr    <- nrow(datxl)
                datxl.rep <- as.data.frame(lapply(datxl,
                                 function(x,nyv,nr){rep(x,rep(nyv,nr))},
                                          nyv=nyv,nr=nr))
                datxl.rep[[ynm]] <- factor(rep(yvals,nr),levels=yvals)
                fy <- ffun(datxl.rep,fmla,response=NULL,Rho,type=1)
                fym <- matrix(fy,ncol=nyv,byrow=TRUE)
                y   <- apply(fym,1,function(x,yvals){sample(yvals,1,prob=x)},yvals=yvals)
                y   <- factor(y,levels=yvals)
                y   <- data.frame(y)
                names(y) <- ynm
            } else { # Bivariate
                if(inherits(Rho,"list")) { # Bivariate independent.
                    type <- 2
                } else if(inherits(Rho,"array")) { # Bivariate dependent
                    if(length(dim(Rho)) == 3) {
                        type <- 3
                    } else {
                        stop("Argument \"Rho\" is of the wrong dimension.\n")
                    }
                } else {
                    stop("Argument \"Rho\" has the wrong structure.\n")
                }
                c1 <- factor(rep(yvals[[1]],length(yvals[[2]])),levels=yvals[[1]])
                c2 <- factor(rep(yvals[[2]],each=length(yvals[[1]])),levels=yvals[[2]])
                nudeAt <- data.frame(c1,c2)
                names(nudeAt) <- response
                M     <- nrow(nudeAt)
                rnudeAt <- as.data.frame(lapply(nudeAt,
                                        function(x,K,M){rep(x,rep(K,M))},K=K,M=M))
                rnudeAt$state <- factor(rep(1:K,M))
                fy <- ffun(rnudeAt,fmla=NULL,response=response,Rho=Rho,type=type)
                fy <- matrix(fy,nrow=K)
                iresp <- sapply(sts,function(k,fy) {
                             sample(1:ncol(fy),1,prob=fy[k,])
                             },fy=fy)
                y <- nudeAt[iresp,]
                rownames(y) <- 1:nrow(y)
            }
        }
    )
    misstify(y,response=names(y),nafrac=mf,fep)
}
