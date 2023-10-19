reviseSigma <- function(r,weights,state) {
    K     <- length(levels(state))
    sr    <- split(r,f=state)
    wr    <- split(weights,f=state)
    sigma <- sapply(1:K,function(k,x,w){
                  ok <- !is.na(x[[k]])
                  xok <- x[[k]][ok]
                  wok <- w[[k]][ok]
                  sqrt(sum(wok*xok^2/sum(wok)))
              },x=sr,w=wr)
    names(sigma) <- paste0("sigma",1:K)
    sigma
}
