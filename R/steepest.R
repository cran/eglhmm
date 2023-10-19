steepest <- function(tvec,theta,data,fmla,distr,size,nbot,ntop) {
K             <- length(levels(data$state))
foo <- function(x,data,fmla,tvec,theta,gv,K) { 
                tvecp <- tvec+x*gv # "tvec perturbed"
                th.new <- reviseTheta(tvecp,theta,distr,fmla,data,size,nbot,ntop)
		getHgl(nd=0,distr,th.new,data,fmla,size,nbot,ntop)$ll
        }
th.new <- reviseTheta(tvec,theta,distr,fmla,data,size,nbot,ntop)
gv     <- getHgl(nd=1,distr,th.new,data,fmla,size,nbot,ntop)$grad
gv     <- gv/sqrt(sum(gv^2))
fooMax <- optimize(foo,c(0,1),maximum=TRUE,data=data,fmla=fmla,
                   tvec=tvec,theta=theta,gv=gv,K=K)
# Numerical inaccuracy can result in a "maximum" that is
# a miniscule amount *less* than than the log likelihood value
# at the original tvec.  (In such circumstances the original
# tvec is a local maximum.)  Using the maximum thus produced
# by optimize() would result in a (small) *decrease* in the
# log likelihood.  We guard against this as follows:
lwbnd <- foo(0,data,fmla,tvec,theta,gv,K)
con   <- if(fooMax$objective < lwbnd) 0 else fooMax$maximum
tvec.new <- tvec + con*gv
attr(tvec.new,"inclTau")      <- attr(tvec,"inclTau")
attr(tvec.new,"preSpecSigma") <- attr(tvec,"preSpecSigma")
tvec.new
}
