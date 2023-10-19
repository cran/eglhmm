reviseTpm <- function(xisum,mixture) {
    if(mixture)  {
        matrix(apply(xisum,2,sum)/sum(xisum),byrow=TRUE,
               nrow=nrow(xisum),ncol=ncol(xisum))
    } else {
        den <- apply(xisum,1,sum)
        if(any(den == 0))
            stop("At least one row of new tpm would be all zeroes.\n")
        xisum/den
    }
}
