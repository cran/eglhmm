reviseRho <- function(data,response,fmla,type) {

lvls  <- lapply(data[response],levels)
if(is.null(data$state)) { # Handle the K=1 case.
    noWts <- TRUE
} else {
    K     <- length(levels(data[["state"]]))
    gamma <- matrix(data[["weights"]],nrow=K)
    noWts <- FALSE
}

if(type==1) {
    lvls  <- lvls[[1]]
    preds <- attr(terms(fmla),"term.labels")
    if(identical(preds,"state")) { # The simple case.
        f    <- data[["state"]]
        ynm  <- as.character(fmla[2])
        y    <- split(data[[ynm]],f=f)[[1]]
        Rho  <- apply(gamma,1,function(x,y){xtabs(x ~ y)},y=y)
        Rho  <- t(t(Rho)/apply(Rho,2,sum))
        rownames(Rho) <- lvls
        colnames(Rho) <- paste0("state",1:K)
        class(Rho)    <- c(class(Rho),"RhoProbForm")
        Rho           <- cnvrtRho(Rho)
    } else {
        mnfit <- suppressWarnings(
                     if(noWts) {
                         multinom(fmla,data=data,trace=FALSE)
                     } else {
                         multinom(fmla,data=cbind(data,weights=gamma),
                                  weights=weights,trace=FALSE)
                     }
        )

# Make provision for there being zero actual appearances in the
# data of some levels of data$y:
        mc  <- coef(mnfit)
        btm <- min(mc)
        z   <- -300 - (abs(btm) - btm)/2
        M   <- matrix(z,nrow=length(lvls),ncol=ncol(mc))
        rownames(M) <- lvls
        M[rownames(mc),] <- mc
        M[1,] <- 0

# The nnet::multinom() function uses the convention that the *first*
# exponent is 0, rather than the last, whereas the convention that
# is used elsewhere in this package is that the *last* exponent is
# set to 0.  The nnet::multinom() convention messes up the procedure
# applied by the getRho() function.
#
# We can rectify the situation by subtracting the last row of
# the matrix produced, by nnet::multinom(), from all rows (including
# itself), making the last row equal to the zero vector.  In the
# simple intercept-only setting this could be done by:
#            x <- c(0,coef(mnfit))
#            x <- x - x[length(x)]
# It is only slightly more complicated in the general case.
# Note that coef(mnfit) is always a matrix, even if there is
# only a single column.
        Rho <- t(M) - M[nrow(M),]
        dumX          <- model.matrix(fmla[-2],data=data[1,])
        rownames(Rho) <- colnames(dumX)
    }
}

if(type==2) {
    Rho <- vector("list",2)
    for(j in 1:2) {
        ynm  <- response[j]
	yj   <- data[[ynm]]
        yj   <- split(data[[ynm]],f=data[["state"]])[[1]]
        Rhoj <- apply(gamma,1,function(x,y){xtabs(x ~ y)},y=yj)
	Rhoj <- t(t(Rhoj)/apply(Rhoj,2,sum))
        rownames(Rhoj) <- lvls[[j]]
        colnames(Rhoj) <- paste0("state",1:K)
        class(Rhoj) <- c(class(Rhoj),"RhoProbForm")
        Rho[[j]]    <- cnvrtRho(Rhoj)
    }
}

if(type==3) {

# Here Rho must be a 3 dimensional array.  We shall set it up
# so that the third dimension ("layers") corresponds to "state".
# Each layer is a matrix whose (i,j)th entry is the estimated
# probability that X = x_i and Y = y_j where X and Y are the
# two variables that are emitted and where x_i and y_j are the
# possible values of X and Y respectfully.
# The tricky bit is handling what the contribution to Rho[i,j,k]
# should be when there are missing values in the obsereved responses.
    s1   <- data$state == 1
    sdat <- split(data[s1,response],f=data[s1,"cf"])
    ym   <- as.matrix(do.call(rbind,sdat))
    aNA  <- any(is.na(as.vector(ym)))
    Rho0 <- array(0,dim=c(sapply(lvls,length),K))
    G    <- array(0,dim=dim(Rho0)+c(1,1,0))
    m    <- dim(G)[1]    
    n    <- dim(G)[2]    
    X    <- factor(ym[,1],levels=c(lvls[[1]],NA),exclude=NULL)
    Y    <- factor(ym[,2],levels=c(lvls[[2]],NA),exclude=NULL)
    for(k in 1:K) {
        G[,,k] <- xtabs(gamma[k,] ~ X + Y,exclude=NULL)
        Rho0[,,k] <- G[-m,-n,k]
        Rho0[,,k] <- Rho0[,,k]/sum(Rho0[,,k])
    }
    if(aNA) {

# Call upon optim() to find the optimum value of Rho,
# Using Rho0 as starting values.
        Rho <- msRho(Rho0,G)
    } else {
        Rho <- Rho0
    }
    dimnames(Rho) <- c(lvls,list(paste0("state",1:K)))
}
Rho
}
