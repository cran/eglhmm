subroutine prephi(ndistr,y,fy,kstate,zeta,sigma,size,d1zeta,d2zeta,d1u,d2u,nd)

# Argument ndistr determines the distribution being used.
# 1 <--> Gaussian.
# 2 <--> Poisson.
# 3 <--> Binomial.
# 4 <--> Dbd. # Does not arise in prephi.
# 5 <--> Multinom.
#
implicit double precision(a-h,o-z)
integer size
dimension fy(kstate), zeta(kstate), sigma(kstate)
dimension d1zeta(kstate), d2zeta(kstate), d1u(kstate),d2u(kstate)
zero = 0.d0
one  = 1.d0
half = 0.5d0

if(ndistr==1) {
# zeta = mu (mean)
    do k = 1,kstate {
        d1zeta(k) = fy(k)*(y-zeta(k))/sigma(k)**2
        d1u(k) = one
        if(nd > 1) {
            #d2zeta(k) = (fy(k)/sigma(k)**2)*(((y-zeta(k))/sigma(k))**2 - one)
            d2zeta(k) = (fy(k)/sigma(k)**2)*((y-zeta(k))**2/sigma(k)**2 - one)
            d2u(k) = zero
        }
    }
} else if(ndistr==2) {
# zeta = lambda
    do k = 1,kstate {
        d1zeta(k) = fy(k)*(y/zeta(k) - one)
        d1u(k)     = zeta(k)
        if(nd > 1) {
            d2zeta(k) = fy(k)*((y/zeta(k) - one)**2 - y/zeta(k)**2)
            d2u(k)    = zeta(k)
        }
    }
} else if(ndistr==3) {
# zeta = p
    do k = 1,kstate {
        d1zeta(k) = fy(k)*(y/zeta(k) - (size-y)/(1-zeta(k)))
        u          = dlog(zeta(k)/(1-zeta(k)))
        emu        = dexp(-u)
        d1u(k)     = emu/(1+emu)**2
        if(nd > 1) {
            d2zeta(k) = fy(k)*((y/zeta(k) - (size-y)/(1-zeta(k)))**2 -
                                (size-y)/(1-zeta(k))**2 - y/zeta(k)**2)
            d2u(k)     = emu*(emu-1)/(1+emu)**3 
        }
   }
} else if(ndistr==5) {
# There is no zeta.
    do k = 1,kstate {
        d1u(k) = 0.d0
        d2u(k) = 0.d0
        d1zeta(k) = 0.d0
        d2zeta(k) = 0.d0
    }
}
return
end

