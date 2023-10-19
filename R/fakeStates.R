fakeStates <- function(data,K,fmla) {
ynm <- as.character(fmla[2])
y   <- data[[ynm]]
o   <- order(y[!is.na(y)])
b   <- seq(0,1+length(o),length=K+1)
u   <- cut(order(o),breaks=b,labels=1:K)
s   <- rep(NA,nrow(data))
s[!is.na(y)] <- u
s   <- factor(s,levels=1:K)
data[["state"]] <- s
ok <- isTRUE(all.equal(is.na(data[[ynm]]), is.na(data[["state"]])))
if(!ok) {
    whinge <- paste0("The missing values in ",ynm," do not match up with\n",
                     "  the missing values in the fake states that were\n",
                     "  created in order to produce starting values for\n",
                     "  the parameters of the model.\n")
    stop(whinge)
}
data
}
