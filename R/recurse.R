recurse <- function(fy,tpm,level=2L) {
#
# Function recurse to revise the ``recursive probabilities'',
# given the transition probability matrix and the function values.
#

# Check on the value of "level".
if(!level %in% (1L:2L))
    stop("The value of \"level\" must be either 1 or 2.\n")

# Set epsilon.
#epsilon <- sqrt(.Machine$double.eps)
epsilon <- .Machine$double.eps

# Get the initial state probability distribution.
ispd <- reviseIspd(tpm)

# Set a bunch of constants:
K   <- nrow(tpm)
lcf <- names(fy)
nc  <- length(lcf)

# Recursive probabilities:
llc   <- vector("list",nc)
gamma <- vector("list",nc)
xi    <- vector("list",nc)
xisum <- vector("list",nc)
names(llc)   <- lcf
names(gamma) <- lcf
names(xi)    <- lcf
names(xisum) <- lcf
for(xl in lcf) {
        xfy <- fy[[xl]]
	ny  <- length(xfy)/K
        sto <- .Fortran(
                'recurse',
		fy=as.double(xfy),
                xispd=as.double(ispd),
                tpm=as.double(tpm),
                epsilon=as.double(epsilon),
                kstate=as.integer(K),
		n=as.integer(ny),
                wrk=double(K*K),
                xlc=double(ny),
                alpha=double(K*ny),
                beta=double(K*ny),
                gamma=double(K*ny),
                xi=double(K*K*ny),
                xisum=double(K*K),
                level=as.integer(level),
		PACKAGE="eglhmm"
        )
	llc[[xl]]   <- sto$xlc
	gamma[[xl]] <- matrix(sto$gamma,nrow=K)
        if(level==2L) xisum[[xl]] <- matrix(sto$xisum,nrow=K)
}
gamma <- do.call(cbind,gamma)
llc   <- unlist(llc)
if(level==1L) {
   return(list(gamma=gamma,llc=llc))
}
xisum <- Reduce("+",xisum)
list(gamma=gamma,xisum=xisum,llc=llc)
}
