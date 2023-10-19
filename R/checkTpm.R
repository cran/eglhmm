checkTpm <- function(tpm) {
if(!inherits(tpm,"matrix"))
    stop("Argument \"tpm\" must be a matrix.\n")
if(nrow(tpm) != ncol(tpm))
    stop("Argument \"tpm\" must be a square matrix.\n")
if(!all(tpm >= 0))
    stop("All entries of tpm must be non-negative.\n")
if(!isTRUE(all.equal(apply(tpm,1,sum),rep(1,nrow(tpm)))))
    stop("All rows of \"tpm\" must sum to 1.\n")
TRUE
}
