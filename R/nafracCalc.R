nafracCalc <- function(data,response) {
#
if(!(is.character(response) & length(response) %in% 1:2)) {
    stop("Argument \"response\" is of an incorrect form.\n")
}
if(!all(response %in% names(data))) {
     stop("Response(s) not found in \"data\".\n")
}

n <- nrow(data)
if(length(response) == 2) {
    y1   <- response[1]
    y2   <- response[2]
    nay1 <- sum(is.na(data[[y1]]))
    nay2 <- sum(is.na(data[[y2]]))
    rslt <- c(nafrac1=nay1/n,nafrac2=nay2/n)
} else {
    y    <- response
    nay  <- sum(is.na(data[[y]]))
    rslt <- c(nafrac=nay/n)
}
rslt
}
