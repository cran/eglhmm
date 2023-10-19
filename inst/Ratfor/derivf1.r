subroutine derivf1(y,gmu,sd,fy,tdm,kstate,npar,npt,nxc,nd,d1f,d2f,kt)
# 1 <--> Gaussian.
# nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
dimension gmu(kstate), sd(kstate)
dimension fy(kstate), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)
logical sigfix # Sigma fixed, i.e. sigma was pre-specified.

# The argument npt is the total number of parameters, always
# including "tau", the transition probabilities, which *may*
# be estimated via the method of moments, as is done in the
# EM algorithm.  The argument npar is the number of parameters
# estimated using the gradient and Hessian and may *not* include
# "tau". If npar == npt then eglhmm() was called with inclTau equal
# to TRUE. Otherwise not.
if(npt > npar) {
    ignore = 0
} else {
    ignore = kstate*(kstate-1)
}
sigfix = (npt == kstate*(kstate-1) + nxc)
if(sigfix) {
    jincr = 0
} else {
    jincr = kstate
}

zero  = 0.d0
one   = 1.d0
three = 3.d0

# If inclTau was set to TRUE then the first "ignore" columns of
# d1f() correspond to the transition probabilities.  The underlying
# pdf/pmf of the model (conceptually denoted "f()") does not depend
# on these probabiities and so the corresponding entries of d1f()
# are left as 0.  If inclTau was set to FALSE then there are no such
# columns to ignore (whence "ignore" is 0.)  Likewise for the first
# "ignore" rows and columns of d2f().
do i = 1,kstate {
    z        = (y-gmu(i))/sd(i)
    dfdmu    = fy(i)*z/sd(i)
    if(!sigfix) {
        dfdzeta = fy(i)*(z**2 - one)
        dfdsigma = dfdzeta/sd(i)
        d1f(i,ignore+i) = dfdzeta
    } else {
        dfdsigma = zero
    }
    do j = 1,nxc {
        d1f(i,ignore+jincr+j) = dfdmu*tdm(j,i)
    }
    if(nd > 1) {
#
# M = d2f[i,,] = 2 x 2 array of matrices
#    --             --
#    | M[1,1] M[1,2] |
#    | M[2,1] M[2,2] |
#    --             --
#
#
        d2fdmu2 = dfdsigma/sd(i)
# M[2,2]
        do j = 1,nxc {
            do k = 1,nxc {
                d2f(i,ignore+jincr+j,ignore+jincr+k) = d2fdmu2*tdm(j,i)*tdm(k,i)
            }
        }
        if(sigfix) next
        d2fdsigma2   = fy(i)*((z**2 - one)**2 + one - three*z**2)/sd(i)**2
        d2fdzeta2    = sd(i)*(dfdsigma + d2fdsigma2*sd(i))
        d2fdmudzeta  = fy(i)*(z**2 - three)*z/sd(i)
# M[1,1]
        d2f(i,ignore+i,ignore+i) = d2fdzeta2
# M[1,2]
        do k = 1,nxc {
            d2f(i,ignore+i,ignore+kstate+k) = d2fdmudzeta*tdm(k,i)
        }
# M[2,1]
        do j = 1,nxc {
            d2f(i,ignore+kstate+j,ignore+i) = d2fdmudzeta*tdm(j,i)
        }
    }
}
return
end
