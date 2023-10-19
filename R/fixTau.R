fixTau <- function(tpm) {
    K   <- nrow(tpm)
    A   <- t(apply(tpm,1,p2expForm))
    tau <- as.vector(A[,-K])
    names(tau) <- paste0("tau",row(tpm)[,-K],col(tpm)[,-K])
    tau
}
