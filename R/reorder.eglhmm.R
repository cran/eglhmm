reorder.eglhmm <- function(x,...) {
#
# First we need to order the state effects.  For the Gaussian, Poisson
# and Binomial distributions this is reasonably straightforward.  Just
# extract the state effect coefficients (the first K entries of phi)
# and determine the decreasing order.
#
# For the Db and Multinom distributions things are more complicated.
#
# For the Dbd there are state effects for each of the two shape parameters,
# alpha and beta.  We determine the expected value of the Dbd for each
# pair of state effects and find the decreasing order of these expected
# values.
# For the Multinom distribution there is a seperate tabular distribution
# for each state.  The y-values are interpreted as being numeric (if
# this makes any sense) or replaced by their indices, and the mean of
# the result is calculated for each state.  The decreasing order of
# these means is then determined.

if(length(x$response)) stop("Cannot reorder bivariate models.\n")

phi   <- x$phi
tpm   <- x$tpm
ispd  <- x$ispd
theta <- x$theta
inclTau      <- attr(x$theta,"inclTau")
preSpecSigma <- attr(x$theta,"preSpecSigma")
K     <- nrow(tpm)
if(x$distr=="Dbd") {
    nbot <- x$nbot
    ntop <- x$ntop
    ia   <- grep("alpha\\.state",names(phi))
    ib   <- grep("beta\\.state",names(phi))
    ast.eff <- phi[ia]
    bst.eff <- phi[ib]
    phi0 <- cbind(ast.eff,bst.eff)
    icpt <- apply(phi0,1,function(x,nbot,ntop){
                  dbd::expValDb(ao=x[1],beta=x[2],ntop=ntop,zeta=(nbot==0))
                  },nbot=nbot,ntop=ntop)
    o    <- order(icpt,decreasing=TRUE)
    ast.eff <- ast.eff[o]
    bst.eff <- bst.eff[o]
    phi[ia] <- ast.eff
    phi[ib] <- bst.eff
} else if(x$distr=="Multinom") {
    ynm  <- as.character(x$formula[2])
    rsplvls <- levels(x$data[[ynm]])
    nmbs <- suppressWarnings(as.numeric(rsplvls))
    if(any(is.na(nmbs))) nmbs <- 1:length(rsplvls)
    pmat <- apply(x$Rho[1:K,],1,expForm2p)
    mns  <- apply(nmbs*pmat,2,sum)
    o    <- order(mns,decreasing=TRUE)
    x$Rho[1:K,] <- x$Rho[1:K,][o,]
    phi <- rho2Phi(x$Rho)
} else {
    st.eff <- phi[1:K]
    o      <- order(st.eff,decreasing=TRUE)
    phi[1:K] <- phi[1:K][o]
}

# We now have the order "o" of the state effects.  Adjust the relevant
# components of "x" according to this order.
tpm    <- tpm[o,o]
ispd   <- ispd[o]
tau    <- fixTau(tpm)
if(x$distr=="Gaussian") {
    if(is.null(preSpecSigma)) {
        zeta  <- log(x$sigma[o])
        names(zeta) <- paste0("zeta",1:K)
    } else {
        zeta <- NULL
        preSpecSigma <- preSpecSigma[o]
    }
} else {
    zeta <- NULL
}
rtheta <- list(tau=tau,zeta=zeta,phi=phi)
attr(rtheta,"inclTau")      <- inclTau
attr(rtheta,"preSpecSigma") <- preSpecSigma
if(is.null(x[["Hessian"]])) {
    if(is.null(x[["gradient"]])) {
        nd <- 0
    } else {
        nd <- 1
    }
} else {
    nd <- 2
}
xxx <- getHgl(nd=nd,distr=x$distr,theta=rtheta,data=x$data,fmla=x$formula,
              size=x$size,nbot=x$nbot,ntop=x$ntop)
x$newlog.like <- xxx$ll # Get rid of this when sure the function is working.
x$phi  <- phi
x$tpm  <- tpm
x$ispd <- ispd
x$theta    <- rtheta
x$Hessian  <- xxx$Hess
x$gradient <- xxx$gradient
X <- with(x,model.matrix(formula[-2],data=data))
y <- with(x,data[,as.character(formula)[2]])

# Some of the components of "x" that need to be adjusted depend on the
# distribution that is in use.

switch(EXPR=x$distr,
    Gaussian = {
        x$mean  <- X%*%phi
        if(is.null(preSpecSigma)) {
            x$sigma <- x$sigma[o]
            sd      <- x$sigma[x$data$state]
        } else {
            x$preSpecSigma <- x$preSpecSigma[o]
            sd             <- x$preSpecSigma[x$data$state]
        }
        x$sd <- sd
        x$mu <- getMu(x$mean,x$data,x$formula)
        fy   <- dnorm(y,mean=x$mean,sd=sd)
    },
    Poisson = {
        x$lambda <- exp(X%*%phi)
        fy       <- dpois(y,x$lam)
    },
    Binomial = {
        x$p <- logistic(X%*%phi)
        fy <- dbinom(y,prob=x$p,size=x$size)
    },
    Dbd = {
        np <- length(phi)
        if(np%%2 != 0)
            stop("There must be an even number of phi coefficients.\n")
        kp      <- np/2
        phia    <- phi[1:kp]
        phib    <- phi[(kp+1):np]
        x$alpha <- X%*%phia
        x$beta  <- X%*%phib
        fy      <- dbd::ddb(y,x$alpha,x$beta,x$ntop,x$nbot==0)
    },
    Multinom  = {
        fy    <- ffun(x$data,fmla=x$formula,response=ynm,Rho=x$Rho,type=1)
    }
)
fy[is.na(y)] <- 1
fy   <- split(fy,f=x$data$cf)
nms  <- unique(x$data$cf)
fy   <- fy[nms] # Keeps "weights" properly aligned with observations.
x$fy <- fy
x$neworder <- o
class(x)   <- c(class(x),"reordered")
x
}
