reviseTau <- function(distr,fmla,data,theta,size,nbot,ntop) {
    fy  <- getModComps(distr,fmla,data,theta,size,nbot,ntop)$fy
    K   <- length(levels(data$state))
    tpm <- getTpm(theta,K)
    rp  <- recurse(fy,tpm)
    tpm <- reviseTpm(rp$xisum,mixture=FALSE)
    fixTau(tpm)
}
