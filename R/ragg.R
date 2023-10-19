ragg <- function(theta,data,fmla,distr,size=NULL,
                 nbot=NULL,ntop=NULL,delta=0.01){
# Rough as guts gradient.  Used for debugging purposes only.
f0   <- getHgl(nd=0,distr=distr,theta=theta,data=data,fmla=fmla,
               size=size,nbot=nbot,ntop=ntop)
npar <- length(theta)
m    <- diag(npar)
gr   <- numeric(npar)
for(i in 1:npar) {
    t1    <- theta + delta*m[i,]
    f1 <- getHgl(nd=0,distr=distr,theta=t1,data=data,fmla=fmla,
                 size=size,nbot=nbot,ntop=ntop)
    gr[i] <- (f1-f0)/delta
}
gr
}
