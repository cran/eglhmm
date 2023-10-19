getMu <- function(gmu,data,fmla) {
if(as.character(fmla[3]) == "1") { # Happens only if K=1.
    mu <- unique(gmu)
    colnames(mu) <- "mu"
    return(mu)
}
preds <- attr(terms(fmla),"term.labels")
numc  <- sapply(data[,preds],is.numeric)
ifac  <- lapply(preds,function(x,ddd){ddd[[x]]},ddd=data)
nms   <- as.character(with(data,interaction(ifac,sep="_")))
mu1   <- unique(gmu)
unms  <- nms[as.numeric(rownames(mu1))]
unms  <- matrix(unlist(strsplit(unms,split="_")),byrow=TRUE,ncol=length(preds))
rownames(unms) <- 1:nrow(unms)
unms <- as.data.frame(unms)
unms <- as.data.frame(lapply(1:ncol(unms),function(k,doit,u) {
                      if(doit[k]) {
                          sprintf("%f",as.numeric(u[,k]))
                      } else u[,k] },doit=numc,u=unms))
colnames(unms) <- preds
rslt <- cbind(unms,mu=round(mu1,6))
rownames(rslt) <- NULL
rslt
}
