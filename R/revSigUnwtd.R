revSigUnwtd <- function(phi,X,y,state) {
#
# Rough as guts.
#
    yhat <- X%*%phi
    rrr  <- y-yhat
    srr  <- split(rrr,f=state)
    sigma <- sapply(srr,sd,na.rm=TRUE)
    names(sigma) <- paste0("sigma",1:length(sigma))
    sigma
}
