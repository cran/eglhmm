subroutine derivf3(y,p,size,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
# 3 <--> Binomial.
# nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
double precision p(kstate)
integer size
dimension fy(kstate), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)

one   = 1.d0
if(npar == nxc) {
    npro = 0
} else {
    npro = kstate*(kstate-1)
}
do i = 1,kstate {
    fctr = (y/p(i) - (size - y)/(one - p(i)))
    dfdp = fy(i)*fctr
    u    = dlog(p(i)/(one-p(i)))
    emu  = dexp(-u)
    h1   = emu/(one+emu)**2
    do j = 1,nxc {
        d1f(i,npro+j) = dfdp*h1*tdm(j,i)
        if(nd > 1) {
            h2 = h1*(emu-one)/(emu+one)
            d2fdp2 = fy(i)*(fctr**2 - y/p(i)**2 - (size - y)/(one - p(i))**2)
            do k = 1,nxc {
                d2f(i,npro+j,npro+k) = (dfdp*h2 + d2fdp*h1**2)*tdm(j,i)*tdm(k,i)
            }
        }
    }
}
return
end
