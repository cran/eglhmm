derivp <- function(npar,K,tau=NULL,expo=FALSE) {

d1p <- array(0,c(K,K,npar))
d2p <- array(0,c(K,K,npar,npar))

if(expo) {
	if(is.null(tau)) stop("When expo is TRUE tau must be supplied.\n")
	E <- exp(cbind(matrix(tau,nrow=K,byrow=TRUE),0))
	Id <- diag(K)
	for(i in 1:K) {
		den <- sum(E[i,])
		for(j in 1:K) {
			for(k in 1:(K-1)) {
				m <- (K-1)*(i-1) + k
				s <- Id[j,k]*den - E[i,k]
				d1p[i,j,m] <- E[i,j]*s/den^2
				for(l in 1:(K-1)) {
					n  <- (K-1)*(i-1) + l 
					ds <- Id[j,k]*E[i,l] - Id[k,l]*E[i,l]
					a  <- (ds + Id[j,l]*s)*den
					b  <- 2*E[i,l]*s
					d2p[i,j,m,n] <- E[i,j]*(a-b)/den^3
				}
			}
		}
	}
}
else {
	for(i in 1:K) {
		for(j in 1:K) {
			for(k in 1:(K-1)) {
				m <- (K-1)*(i-1) + k
				d1p[i,j,m] <- if(j<K) 1 else -1
			}
		}
	}
}
list(d1p=d1p,d2p=d2p)
}
