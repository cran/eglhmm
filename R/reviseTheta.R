reviseTheta <- function(tvec,theta,distr,fmla,data,size,nbot,ntop) {
# Put together a new theta list, using the tvec ("theta vector") that
# has been updated by the LM step. Revise tau, if inclTau is FALSE,
# and put the revised tau, (however it was obtainded), the revised
# version of zeta (if this was indeed revised), and the revised
# version of phi into the new theta list.
    inclTau      <- attr(theta,"inclTau")
    preSpecSigma <- attr(theta,"preSpecSigma")
    K   <- length(levels(data$state))
    ttt <- dotzp(tvec,K,inclTau,preSpecSigma)
    if(!inclTau) ttt$tau <- theta$tau
    theta <- ttt
    attr(theta,"inclTau")      <- inclTau
    attr(theta,"preSpecSigma") <- preSpecSigma
    if(!inclTau) theta$tau     <- reviseTau(distr,fmla,data,theta,size,nbot,ntop)
    theta
}
