dotzp <- function(tvec,K,inclTau,preSpecSigma) {
# The letters "do" indicate "dig out".  I.e.  dig out the values of
# tau, zeta and phi from the vector tvec ("theta vector").  Note that
# either or both of tau and zeta may be absent.  Thus either or
# both of the tau and zeta in the list created below may be NULL.
# However phi is always non-NULL
npro  <- if(inclTau) K*(K-1) else 0
tau   <- if(inclTau) tvec[1:npro] else NULL
zeta  <- if(is.null(preSpecSigma)) tvec[(npro+1):(npro+K)] else NULL
nzeta <- length(zeta)
ntz   <- npro + nzeta
phi   <- tvec[(ntz+1):length(tvec)]
list(tau=tau,zeta=zeta,phi=phi)
}
