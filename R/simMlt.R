simMlt <- function(formula,response,distr,data,ispd,tpm,phi,Rho,
                   sigma,size,ntop,zeta,missFrac,fep) {
# Function simMlt.  Simulates data from all cells, in turn,
# of a Hidden Markov Model.

cf   <- data$cf
lcf  <- levels(cf)
K    <- nrow(tpm)
ynm  <- if(is.null(formula)) NULL else as.character(formula)[2]
rslt <- vector("list",length(lcf))
names(rslt) <- lcf
if(distr == "Multinom") {
    if(is.null(ynm)) {
        yvals <- lapply(data[response],levels)
    } else {
        yvals <- levels(data[[ynm]])
    }
    mlp   <- NULL
} else {
    yvals <- NULL
}
for(xl in lcf) {
    nt <- sum(cf==xl)/K
    if(distr == "Multinom") {
        datxl <- data[data$cf == xl,]
    } else {
        datxl <- NULL
        X     <- model.matrix(formula[-2],data[cf==xl,])
        if(distr=="Dbd") {
            np <- length(phi)
            if((np %% 2) != 0) {
                whinge <- paste0("For the \"Dbd\" case, \"phi\" must",
                                 " have an even length.\n")
                stop(whinge)
            }
            npo2 <- np/2
            phia <- phi[1:npo2]
            phib <- phi[-(1:npo2)]
	        mlpa <- matrix(X%*%phia,ncol=K,byrow=TRUE)
	        mlpb <- matrix(X%*%phib,ncol=K,byrow=TRUE)
            mlp  <- list(mlpa=mlpa,mlpb=mlpb)
        } else {
	        mlp <- matrix(X%*%phi,ncol=K,byrow=TRUE)
        }
        sig <- if(distr=="Gaussian") rep(sigma,length=K) else NULL
    }
    rslt[[xl]] <- simSngl(distr,tpm,ispd,nt,mlp,Rho,yvals,datxl,fmla=formula,
                          response,sig,size,ntop,zeta,mf=missFrac,fep)
}
y     <- do.call(rbind,rslt)
K     <- length(ispd)
n     <- nrow(data)
ind   <- K*(1:(n/K))
norep <- data[ind,]
norep[["state"]]   <- NULL
norep[["cf"]]      <- NULL
norep[["weights"]] <- NULL
norep[names(y)]    <- y
row.names(norep)   <- 1:nrow(norep)
norep
}
