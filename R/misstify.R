misstify <- local({

fix1 <- function(x,nafrac,present){
# Take an array (possibly with only 1 column) and randomly
# replace entries with NA, with probability given by the
# appropriate entries of nafrac.  If "present" is TRUE,
# make all entries of the first row *non*-missing.
    pad <- if(present) 0 else NULL
    s   <- if(present) 1 else 0
    ina1 <- c(pad,rbinom(nrow(x)-s,1,nafrac[1]))
    x[ina1 == 1,1] <- NA
    if(ncol(x)==2) {
        ina2 <- c(pad,rbinom(nrow(x)-s,1,nafrac[2]))
        x[ina2 == 1,2] <- NA
    }
    x
}

fix2 <- function(x,nafrac,p2){
# Used only if "present" is TRUE and only if the data are bivariate.
# Adjust the first row of x so that it is not necessarily the case
# that *both* entries are non-missing.  Note that at least one
# entry of the result will be non-missing.  The argument p2 is
# the probability both entries of the first row of the result are
# *non*-missing.  It defaults to the probability that both entries
# of an arbitrary row of the input data are non-missing, given that
# at least one is non-missing.

if(ncol(x)==1) return(x)
    pq <- nafrac/sum(nafrac)
    x1 <- x[1,]
    j  <- sample(1:2,1,prob=pq)
    ind <- rbinom(1,1,p2)
    if(ind==0) x1[j] <- NA
    x[1,] <- x1
    x
}

function(y,response,nafrac,fep=NULL) {
    nafrac <- rep(nafrac,length=length(response))
    if(!all(nafrac >= 0 & nafrac < 1)) {
        whinge <- paste0("All entries of \"nafrac\" must be probabilities",
                         " strictly less than 1.\n")
        stop(whinge)
    }
    x <- as.matrix(y[,response])
    bivar <- length(nafrac) == 2
    if(is.null(fep)) {
        fep <- list(present=TRUE)
    }
    if(bivar) {
        if(length(fep) == 1) {
            fep <- c(fep,list(p2=prod(1-nafrac)/(1 - prod(nafrac))))
        } else if(fep[[2]] < 0 | fep[[2]] > 1) {
            whinge <- paste0("Component \"p2\" of \"fep\" must be a probability.\n")
            stop(whinge)
        }
    }
    x <- fix1(x,nafrac=nafrac,present=fep[[1]])
    if(bivar) {
        if(fep[[1]] & fep[[2]] < 1) {
           x <- fix2(x,nafrac=nafrac,p2=fep[[2]])
        }
    }
    x <- as.data.frame(x)
    for(j in 1:ncol(x)) {
        if(inherits(y[,response[j]],"factor")) {
             x[,j] <- factor(x[,j],levels=levels(y[,response[j]]))
        }
    }
    y[,response] <- x
    y
}})
