subroutine derivfdbd(fy,zeta,size,y,ymiss,tdm,kstate,npar,nphi,
                     nxc,nbot,ntop,d1f,d2f,nd,d1a,d1b,d2aa,d2ab,d2bb)
# Derivatives of the likelihood in the "special"  case of
# the Dbd distribution.
# Note that nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
integer ymiss
integer size
dimension fy(kstate), zeta(kstate,2), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)
dimension d1a(kstate), d1b(kstate)
dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)

zero = 0.d0
do i = 1,kstate {
    do j = 1,npar {
        d1f(i,j) = zero
        do k = 1,npar {
            d2f(i,j,k) = zero
        }
    }
}
if(ymiss > 0) return

# Calculate the derivatives of the likelihood w.r.t. alpha and beta.
call derivfab(y,fy,kstate,zeta,nbot,ntop,nd,d1a,d1b,d2aa,d2ab,d2bb)
npro     = kstate*(kstate-1)
nphio2   = nphi/2
nphio2p1 = nphio2+1

# Do the first derivatives.
do i = 1,kstate {
    do j = 1,nxc {
        d1f(i,npro+j) = d1a(i)*tdm(j,i)
        d1f(i,npro+nxc+j) = d1b(i)*tdm(j,i)
    }
}

# Do the second derivatives.
if(nd > 1) {
    do i = 1,kstate {
        do j = 1,nxc{
            do k = 1,nxc {
                d2f(i,npro+j,npro+k) = d2aa(i)*tdm(j,i)*tdm(k,i)
                d2f(i,npro+j,npro+nxc+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
                d2f(i,npro+nxc+j,npro+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
                d2f(i,npro+nxc+j,npro+nxc+k) = d2bb(i)*tdm(j,i)*tdm(k,i)
            }
        }
    }
}
return
end
