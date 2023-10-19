binForm <- function(fmla) {
    ynm      <- as.character(fmla)[2]
    wise     <- paste0("cbind(",ynm,",size-",ynm,")")
    cform    <- as.character(fmla)
    cform[2] <- wise
    as.formula(paste(cform[c(2,1,3)],collapse=" "))
}
