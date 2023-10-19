reviseIspd <- function(tpm) {
# Function reviseIspd.  To revise the initial state probability
# distribution.

# Ispd is taken to be the steady state distribution.
    ok <- checkTpm(tpm)
    if(!ok)
        stop("Argument \"tpm\" is not of the correct form.\n")
    K <- nrow(tpm)
    if(isTRUE(all.equal(tpm,diag(K)))) {
        whinge <- paste0("The transition probability matrix is equal to the",
                " identity.\n  The initial state probability distribution is",
                " undefined.\n")
        stop(whinge)
    } else { 
        eee <- eigen(t(tpm))
        k   <- which(as.logical(match(round(eee$values,6),1,nomatch=0)))
        if(length(k) == 1) {
        v <- Re(eee$vectors[,k])
            ispd <- v/sum(v)
        } else {
            whinge <- paste0("The ispd is undefined.  The",
                             " eigenvalues of\n","  t(tpm) are ",
                             paste(eee$values,collapse=" "),".\n")
            stop(whinge)
        }
    }
    ispd
}
