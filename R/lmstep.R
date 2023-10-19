lmstep <- function(theta,data,fmla,distr,size,nbot,ntop,lmc) {
#
# Function lmstep to do one Levenberg-Marquardt step.
# Note that theta is a *list* with components tau, possibly zeta,
# and phi; the zeta component may be NULL.  Some or all of these
# components will later be concatenated to form a vector.
#
K    <- length(levels(data[["state"]]))
xxx  <- getHgl(nd=2,distr,theta,data,fmla,size,nbot,ntop)
Hess <- xxx$Hess
ev   <- eigen(Hess)$values
ok   <- all(Im(ev) == 0)
if(!ok) {
    miev <- max(abs(Im(ev)))
    if(miev <= sqrt(.Machine$double.eps)) {
        ev <- Re(ev)
    } else {
        stop("Hessian has complex eigenvalues.\n")
    }
}
osag    <- sum(abs(xxx$gradient))

# Now form the vector tvec ("theta vector") by selecting the
# appropriate components of theta and catenating them. to create
# a vector.  The vector tvec gets updated by the Levenberg-Marquardt
# step.
inclTau      <- attr(theta,"inclTau")
preSpecSigma <- attr(theta,"preSpecSigma")
tau  <- if(inclTau) theta$tau else NULL
zeta <- theta$zeta
phi  <- theta$phi
tvec <- c(tau,zeta,phi)
npar <- length(tvec)
old.ll  <- xxx$ll
steepit <- FALSE
repeat {
    lll     <- lmc + max(0,ev)
    HessMod <- Hess - lll*diag(npar)
    tvec.new  <- try(tvec - solve(HessMod,xxx$gradient),TRUE)
    if(inherits(tvec.new,"try-error")) {
        steepit <- TRUE
        break
    }

# Now update theta by extracting the relevant bits from tvec.new.
# Note that if inclTau is FALSE, then tau is not updated by the LM step.
# Instead it gets updated via the method of moments inside reviseTheta().
    th.new <- reviseTheta(tvec.new,theta,distr,fmla,data,size,nbot,ntop)
    yyy    <- try(getHgl(nd=1,distr,th.new,data,fmla,size,nbot,ntop),TRUE)
    if(inherits(yyy,"try-error")) {
        steepit <- TRUE
        break
    }
    new.ll <- yyy$ll
    nsag   <- sum(abs(yyy$gradient))
    ok     <- {all(is.finite(c(old.ll,new.ll,osag,nsag))) &&
                      (new.ll > old.ll & nsag < osag)}
    if(ok) break
    lmc <- 10*lmc
    if(1/lmc < sqrt(.Machine$double.eps)) {
        steepit <- TRUE
        break
    }
}
if(steepit) {
    attr(tvec,"inclTau")      <- inclTau
    attr(tvec,"preSpecSigma") <- preSpecSigma
    tvec.new <- steepest(tvec,theta,data,fmla,distr,size,nbot,ntop)
    th.new <- reviseTheta(tvec.new,theta,distr,fmla,data,size,nbot,ntop)
    yyy    <- getHgl(nd=1,distr,th.new,data,fmla,size,nbot,ntop)
    lmc    <- 1
}

rslt  <- list(theta=th.new,ll=yyy$ll,gradient=yyy$gradient,lmc=lmc)
attr(rslt,"usedSteep") <- steepit
rslt
}
