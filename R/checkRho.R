checkRho <- function(Rho,K,rsplvls,indep) {
if(inherits(rsplvls,"list")) { # Bivariate.
    if(indep) {
        if(!inherits(Rho,"list"))
            stop("In the bivariate independent setting, \"Rho\" must be a list.\n")
        for(i in 1:2) {
            if(!inherits(Rho[[i]],"matrix") || ncol(Rho[[i]]) != length(rsplvls[[i]])) {
                whinge <- paste0("Argument \"Rho\"[[",i,"]] is not of the right form.\n")
                stop(whinge)
            }
        }
    } else {
        if(!inherits(Rho,"array"))
            stop("In the bivariate dependent setting, \"Rho\" must be an array.\n")
        m1 <- length(rsplvls[[1]])
        m2 <- length(rsplvls[[2]])
        if(!isTRUE(all.equal(unname(dim(Rho)),c(m1,m2,K))))
            stop("Argument \"Rho\" is of the wrong dimension.\n")
        if(!all(Rho>=0))
            stop("All entries of \"Rho\" must be non-negative.\n")
        some <- unname(apply(Rho,3,sum))
        if(!isTRUE(all.equal(some,rep(1,K))))
            stop("Each \"layer\" of \"Rho\" must sum to 1.\n")
    }
} else { # Univariate
    if(!inherits(Rho,"matrix") || ncol(Rho) != length(rsplvls))
        stop("Argument \"Rho\" is not of the right form.\n")
}
TRUE
}
