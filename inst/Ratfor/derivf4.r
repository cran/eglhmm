subroutine derivf4(y,ashp,bshp,nbot,ntop,fy,tdm,kstate,npar,
                     nxc,nd,d1f,d2f,d1a,d1b,d2aa,d2ab,d2bb)
# 4 <--> Dbd
# nd = number of derivatives to calculate.
#
implicit double precision(a-h,o-z)
dimension fy(kstate), ashp(kstate), bshp(kstate), tdm(nxc,kstate)
dimension d1f(kstate,npar), d2f(kstate,npar,npar)
dimension d1a(kstate), d1b(kstate)
dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)

# Calculate the derivatives of the likelihood w.r.t. alpha and beta.
# This is messy and is hence isolated inside the subroutine derivfab().
call derivfab(y,fy,kstate,ashp,bshp,nbot,ntop,nd,d1a,d1b,d2aa,d2ab,d2bb)

# Now carry on with calculating the derivatives of f with respect
# to phi.
if(npar == 2*nxc) {
    npro = 0
} else {
    npro = kstate*(kstate-1)
}

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
