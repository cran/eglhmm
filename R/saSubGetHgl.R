saSubGetHgl <- function(nd,fy,gmu,sd,y,tdm,tpm,xispd,
                        d1pi,d2pi,kstate,npar,npt,nxc,d1p,d2p,
                        alpha,alphw,a,b,aw,bw,xlc,hess,xl) {
#
# Note: nd is the number of derivatives to calculate; may be 0, 1 or 2.
# 0 <--> just calculate the log likelihood
# 1 <--> calculate the log likelihood and the gradient
# 2 <--> calculate the log likelihood, the gradient and the hessian
#
# Note: tdm = transposed design matrix.
#

# Turn fy, gmu, sd, a, aw and hess into matrices.
fy   <- matrix(fy,nrow=kstate)
gmu  <- matrix(gmu,nrow=kstate)
sd   <- matrix(sd,nrow=kstate)
a    <- matrix(a,nrow=kstate)
aw   <- matrix(aw,nrow=kstate)
hess <- matrix(hess,nrow=npar)

# Turn b and bw into arrays.
b  <- array(b,dim=c(kstate,npar,npar))
bw <- array(bw,dim=c(kstate,npar,npar))

# Set zero.
zero <- 0

kt <- 1
if(nd >= 1) {
    jtdm  <- (1+(kt-1)*kstate):(kt*kstate)
    rslt1 <- sasubrf1(y[kt],gmu[1:kstate,kt],sd[1:kstate,kt],fy[1:kstate,kt],
                       tdm[1:nxc,jtdm],kstate,npar,npt,nxc,nd)
    d1f <-rslt1$d1f
    d2f <-rslt1$d2f
}
sxlc <- zero

kstart <- npt - npar
for(j in 1:kstate) {
    alpha[j] <- xispd[j]*fy[j,kt]
    sxlc <- sxlc + alpha[j]
    if(nd >= 1) {
        for(k1 in 1:npar) {
            if(is.na(y[1]) == 1) {
                d1fx1 <- zero
            } else {
                d1fx1 <- d1f[j,k1]
            }
            a[j,k1] <- xispd[j]*d1fx1 + fy[j,kt]*d1pi[j,kstart+k1]
            #cat(j,k1,d1pi[j,kstart+k1],"\n")
            if(nd == 2) {
                for(k2 in 1:npar) {
                    if(is.na(y[kt]) == 1) {
                        d1fx2 <- zero
                    } else {
                        d1fx2 <- d1f[j,k2]
                    }
                    if(is.na(y[1]) == 1) {
                        d2fx <- zero
                    } else {
                        d2fx <- d2f[j,k1,k2]
                    }
                    b[j,k1,k2] <- xispd[j]*d2fx + d1pi[j,kstart+k1]*d1fx2 +
                                                  d1pi[j,kstart+k2]*d1fx1 +
                                                  fy[j,kt]*d2pi[j,kstart+k1,kstart+k2]
                }
            }
        }
    }
}
xlc[1] <- sxlc

for(j in 1:kstate) {
    alpha[j] <- alpha[j]/sxlc
}

n <- length(y)
if(n>1) {
    for(kt in 2:n) {
        if(nd >= 1) {
            jtdm  <- (1+(kt-1)*kstate):(kt*kstate)
            rslt2 <- sasubrf1(y[kt],gmu[1:kstate,kt],sd[1:kstate,kt],fy[1:kstate,kt],
                              tdm[1:nxc,jtdm],kstate,npar,npt,nxc,nd)
            d1f <-rslt2$d1f
            d2f <-rslt2$d2f
        }
        if(nd == 2) {
# Do the b's:
            for(j in 1:kstate) {
                for(k1 in 1:npar) {
                    if(is.na(y[kt]) == 1) {
                        d1fx1 <- zero
                    } else {
                        d1fx1 <- d1f[j,k1]
                    }
                    for(k2 in 1:npar) {
                        if(is.na(y[kt]) == 1) {
                            d1fx2 <- zero
                            d2fx  <- zero
                        } else {
                            d1fx2 <- d1f[j,k2]
                            d2fx  <- d2f[j,k1,k2]
                        }
                        vvv <- zero
                        xxx <- zero
                        yy1 <- zero
                        yy2 <- zero
                        zz1 <- zero
                        zz2 <- zero
                        www <- zero
                        for(i in 1:kstate) {
                            vvv <- vvv + alpha[i]*d2p[i,j,kstart+k1,kstart+k2]
                            xxx <- (xxx + a[i,k1]*d1p[i,j,kstart+k2] +
                                          a[i,k2]*d1p[i,j,kstart+k1] +
                                          b[i,k1,k2]*tpm[i,j])
                            yy1 <- yy1 + alpha[i]*d1p[i,j,kstart+k2]
                            yy2 <- yy2 + a[i,k2]*tpm[i,j]
                            zz1 <- zz1 + alpha[i]*d1p[i,j,kstart+k1]
                            zz2 <- zz2 + a[i,k1]*tpm[i,j]
                            www <- www + alpha[i]*tpm[i,j]
                        }
                        vvv <- fy[j,kt]*vvv
                        xxx <- fy[j,kt]*xxx/sxlc
                        yyy <- d1fx1*(yy1 + yy2/sxlc)
                        zzz <- d1fx2*(zz1 + zz2/sxlc)
                        www <- d2fx*www
                        bw[j,k1,k2] <- vvv + xxx + yyy + zzz + www
                    }
                }
            }
            for(j in 1:kstate) {
                for(k1 in 1:npar) {
                    for(k2 in 1:npar) {
                        b[j,k1,k2] <- bw[j,k1,k2]
                    }
                }
            }
        }
        if(nd >= 1) {
# Do the a's:
            for(j in 1:kstate) {
                for(k in 1:npar) {
                    if(is.na(y[kt]) == 1) {
                        d1fx <- zero
                    } else {
                        d1fx <- d1f[j,k]
                    }
                    xxx <- zero
                    yyy <- zero
                    zzz <- zero
                    for(i in 1:kstate) {
                        xxx <- xxx + alpha[i]*d1p[i,j,kstart+k]
                        yyy <- yyy + a[i,k]*tpm[i,j]
                        zzz <- zzz + alpha[i]*tpm[i,j]
                    }
                    aw[j,k] <- fy[j,kt]*(xxx + yyy/sxlc) + d1fx*zzz
                }
            }
            for(j in 1:kstate) {
                for(k in 1:npar) {
                    a[j,k] <- aw[j,k]
                }
            }
        }

# Update the alpha's:
        sxlc <- zero
        for(j in 1:kstate) {
            alphw[j] <- zero
            for(i in 1:kstate) {
                alphw[j] <- alphw[j] + alpha[i]*tpm[i,j]
            }
            alphw[j] <- fy[j,kt]*alphw[j]
            sxlc <- sxlc + alphw[j]
        }
        xlc[kt] <- sxlc
        for(j in 1:kstate) {
            alpha[j] <- alphw[j]/sxlc
        }
    }
}

if(nd == 2) {
# Build the Hessian increment.
    for(k1 in 1:npar) {
        for(k2 in 1:npar) {
            xxx <- zero
            yyy <- zero
            zzz <- zero
            for(i in 1:kstate) {
                xxx <- xxx + b[i,k1,k2]
                yyy <- yyy + a[i,k1]
                zzz <- zzz + a[i,k2]
            }
            hess[k1,k2] <- (xxx - yyy*zzz/sxlc)/sxlc
        }
    }
}
list(xlc=xlc,a=a,hess=hess)
}
