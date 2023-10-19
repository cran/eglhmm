subroutine derivf2(y,lambda,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
# 2 <--> Poisson.
# nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
double precision lambda(kstate)
dimension fy(kstate), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)


one  = 1.d0
if(npar == nxc) {
    npro = 0
} else {
    npro = kstate*(kstate-1)
}
# Got here.
do i = 1,kstate {
    fctr   = y/lambda(i) - one
    dfdlam = fy(i)*fctr
    do j = 1,nxc {
        d1f(i,npro+j)    = dfdlam*lambda(i)*tdm(j,i)
        if(nd > 1) {
            d2fdlam2     = fy(i)*(fctr**2 - y/lambda(i))
            do k = 1,nxc {
                xxt                  = tdm(j,i)*tdm(k,i)
                d2f(i,npro+j,npro+k) = (dfdlam*lambda(i) + d2fdlam2*lambda(i)**2)*xxt
            }
        }
    }
}
return
end
