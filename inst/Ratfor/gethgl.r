subroutine gethgl(nd,ndistr,fy,gmu,sd,lambda,p,ashp,bshp,phimat,y,ymiss,
                  tdm,nind,n,tpm,xispd,size,nbot,ntop,d1pi,d2pi,kstate,
                  npar,npt,nyv,nxc,d1p,d2p,d1f,d2f,d1a,d1b,d2aa,d2ab,d2bb,
                  alpha,alphw,a,b,aw,bw,xlc,hess)
#
# Note: nd is the number of derivatives to calculate; may be 0, 1 or 2.
# 0 <--> just calculate the log likelihood
# 1 <--> calculate the log likelihood and the gradient
# 2 <--> calculate the log likelihood, the gradient and the hessian
#
# Note: nind is the number of rows of the data corresponding
# to the cell that is currently being dealt with, *after*
# replication.  This data frame was obtained by replicating each
# row of the original data frame "kstate" times, one copy of
# the row for each state.
#
# In contrast n is the number of rows in the original data,
# i.e. the number of observations (in the given cell) without
# replication.  Thus nind = n*kstate.
#
# Note: tdm = transposed design matrix.
#
implicit double precision(a-h,o-z)
integer ymiss(n), size
dimension fy(kstate,n),gmu(kstate,n),sd(kstate,n)
double precision lambda(kstate,n)
dimension p(kstate,n)
dimension ashp(kstate,n),bshp(kstate,n),phimat(nxc,nyv)
dimension y(n),tdm(nxc,nind)
dimension tpm(kstate,kstate),xispd(kstate)
dimension d1pi(kstate,npt),d2pi(kstate,npt,npt)
dimension d1p(kstate,kstate,npt),d2p(kstate,kstate,npt,npt)
dimension d1f(kstate,npar),d2f(kstate,npar,npar)
dimension d1a(kstate),d1b(kstate)
dimension d2aa(kstate),d2ab(kstate),d2bb(kstate)
dimension alpha(kstate),alphw(kstate)
dimension a(kstate,npar),b(kstate,npar,npar)
dimension aw(kstate,npar),bw(kstate,npar,npar)
dimension xlc(n),hess(npar,npar)

# Set zero.
zero = 0.d0

kt = 1
if(nd >= 1) {
    kstart = 1
    call derivf(ndistr,y(kt),ymiss(kt),fy(1,kt),phimat,tdm(1,kstart),
                gmu(1,kt),sd(1,kt),lambda(1,kt),p(1,kt),ashp(1,kt),
                bshp(1,kt),kstate,npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b,
                d2aa,d2ab,d2bb,nd,d1f,d2f)
}
sxlc = zero

kstpr = npt - npar
do j = 1,kstate {
    alpha(j) = xispd(j)*fy(j,kt)
    sxlc = sxlc + alpha(j)
    if(nd >= 1) {
        do k1 = 1,npar {
            if(ymiss(1) == 1) {
                d1fx1 = zero
            } else {
                d1fx1 = d1f(j,k1)
            }
            a(j,k1) = xispd(j)*d1fx1 + fy(j,kt)*d1pi(j,kstpr+k1)
            if(nd == 2) {
                do k2 = 1,npar {
                    if(ymiss(1) == 1) {
                        d1fx2 = zero
                    } else {
                        d1fx2 = d1f(j,k2)
                    }
                    if(ymiss(1) == 1) {
                        d2fx = zero
                    } else {
                        d2fx = d2f(j,k1,k2)
                    }
                    b(j,k1,k2) = xispd(j)*d2fx + d1pi(j,kstpr+k1)*d1fx2 +
                                                 d1pi(j,kstpr+k2)*d1fx1 +
                                                 fy(j,kt)*d2pi(j,kstpr+k1,kstpr+k2)
                }
            }
        }
    }
}
xlc(1) = sxlc

do j = 1,kstate {
    alpha(j) = alpha(j)/sxlc
}

if(n>1) {
    do kt = 2,n {
        if(nd >= 1) {
            kstart = 1+(kt-1)*kstate
            call derivf(ndistr,y(kt),ymiss(kt),fy(1,kt),phimat,tdm(1,kstart),
                        gmu(1,kt),sd(1,kt),lambda(1,kt),p(1,kt),ashp(1,kt),
                        bshp(1,kt),kstate,npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b,
                        d2aa,d2ab,d2bb,nd,d1f,d2f)
        }
        if(nd == 2) {
# Do the b's:
            do j = 1,kstate {
                do k1 = 1,npar {
                    if(ymiss(kt) == 1) {
                        d1fx1 = zero
                    } else {
                        d1fx1 = d1f(j,k1)
                    }
                    do k2 = 1,npar {
                        if(ymiss(kt) == 1) {
                            d1fx2 = zero
                            d2fx  = zero
                        } else {
                            d1fx2 = d1f(j,k2)
                            d2fx  = d2f(j,k1,k2)
                        }
                        vvv = zero
                        xxx = zero
                        yy1 = zero
                        yy2 = zero
                        zz1 = zero
                        zz2 = zero
                        www = zero
                        do i = 1,kstate {
                            vvv = vvv + alpha(i)*d2p(i,j,kstpr+k1,kstpr+k2)
                            xxx = (xxx + a(i,k1)*d1p(i,j,kstpr+k2) +
                                         a(i,k2)*d1p(i,j,kstpr+k1) +
                                         b(i,k1,k2)*tpm(i,j))
                            yy1 = yy1 + alpha(i)*d1p(i,j,kstpr+k2)
                            yy2 = yy2 + a(i,k2)*tpm(i,j)
                            zz1 = zz1 + alpha(i)*d1p(i,j,kstpr+k1)
                            zz2 = zz2 + a(i,k1)*tpm(i,j)
                            www = www + alpha(i)*tpm(i,j)
                        }
                        vvv = fy(j,kt)*vvv
                        xxx = fy(j,kt)*xxx/sxlc
                        yyy = d1fx1*(yy1 + yy2/sxlc)
                        zzz = d1fx2*(zz1 + zz2/sxlc)
                        www = d2fx*www
                        bw(j,k1,k2) = vvv + xxx + yyy + zzz + www
                    }
                }
            }
            do j = 1,kstate {
                do k1 = 1,npar {
                    do k2 = 1,npar {
                        b(j,k1,k2) = bw(j,k1,k2)
                    }
                }
            }
        }
        if(nd >= 1) {
# Do the a's:
            do j = 1,kstate {
                do k = 1,npar {
                    if(ymiss(kt) == 1) {
                        d1fx = zero
                    } else {
                        d1fx = d1f(j,k)
                    }
                    xxx = zero
                    yyy = zero
                    zzz = zero
                    do i = 1, kstate {
                        xxx = xxx + alpha(i)*d1p(i,j,kstpr+k)
                        yyy = yyy + a(i,k)*tpm(i,j)
                        zzz = zzz + alpha(i)*tpm(i,j)
                    }
                    aw(j,k) = fy(j,kt)*(xxx + yyy/sxlc) + d1fx*zzz
                }
            }
            do j = 1,kstate {
                do k = 1,npar {
                    a(j,k) = aw(j,k)
                }
            }
        }

# Update the alpha's:
        sxlc = zero
        do j = 1,kstate {
            alphw(j) = zero
            do i = 1,kstate {
                alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
            }
            alphw(j) = fy(j,kt)*alphw(j)
            sxlc = sxlc + alphw(j)
        }
        xlc(kt) = sxlc
        do j = 1,kstate {
            alpha(j) = alphw(j)/sxlc
        }
    }
}

if(nd == 2) {
# Build the Hessian increment.
    do k1 = 1,npar {
        do k2 = 1,npar {
            xxx = zero
            yyy = zero
            zzz = zero
            do i = 1,kstate {
                xxx = xxx + b(i,k1,k2)
                yyy = yyy + a(i,k1)
                zzz = zzz + a(i,k2)
            }
            hess(k1,k2) = (xxx - yyy*zzz/sxlc)/sxlc
        }
    }
}

return
end
