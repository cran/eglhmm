subroutine derivf5(y,phimat,tdm,kstate,npar,nxc,nyv,nd,
                   d1a,d1b,d2aa,d2ab,d2bb,d1f,d2f)
# 5 <--> Multinom.
# nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
dimension phimat(nxc,nyv), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)
dimension d1a(kstate), d1b(kstate)
dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
integer ell, r, s, t, u
iy = idnint(y)
if(npar == nxc) {
    npro = 0
} else {
    npro = kstate*(kstate-1)
}
nyvm1 = nyv - 1 # phi_{nyv,j} == 0 for all j, so derivs. are 0.
do k = 1,kstate {
    call pmf(iy,tdm(1,k),phimat,nyv,nxc,pmfy)
    do r = 1,nxc {
        do s = 1,nyvm1 {
            j = npro + (r-1)*nyvm1 + s
            call pmf(s,tdm(1,k),phimat,nyv,nxc,pmfs)
            call delta(iy,s,iysd)
            part2 =  iysd-pmfs
            xr = tdm(r,k)
            d1f(k,j) = pmfy*part2*xr
            if(nd>1) {
                do t = 1, nxc {
                    do u = 1,nyvm1 {
                        ell = npro + (t-1)*nyvm1 + u
                        call pmf(u,tdm(1,k),phimat,nyv,nxc,pmfu)
                        call delta(s,u,isud)
                        call delta(iy,u,iyud)
                        part3 = pmfs*pmfu - pmfs*isud
                        part4 = (iysd - pmfs)*(iyud - pmfu)
                        xt = tdm(t,k)
                        d2f(k,j,ell) = pmfy*(part3 + part4)*xr*xt
                    }
                }
            }
        }
    }
}
return
end
