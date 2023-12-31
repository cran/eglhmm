checkYval <- function(yval,Rho,type,warn=TRUE) {
fname <- as.character(sys.call(-1))[1]
if(is.na(fname)) fname <- "call from the command line"

# Univariate.
if(type==1) {
    rn <- levels(Rho$y)
    if(is.null(rn))
        stop("In ",fname," something is wrong; \"Rho$y\" is not a factor.\n")
    yval <- as.character(yval)
    if(all(yval %in% rn)) return(Rho)
    whinge <- paste0("In ",fname," some y values do not match ",
                     "the levels of \"Rho$y\".\n")
    stop(whinge,call.=FALSE)
}

# Bivariate independent.
if(type==2) {
    yval <- lapply(yval,as.character)
    for(j in 1:2) {
        rn <- rownames(Rho[[j]])
        yv <- yval[[j]]
        wstrng <- paste0("\"Rho[[",j,"]]\"")
        if(!is.null(rn)) {
            if(!all(yv %in% rn))
                stop(paste0("In ",fname," some y values do not match the row names\n",
                           "of ",wstrng,".\n"),call.=FALSE)
        } else {
            if(length(yval[[j]]) != nrow(Rho[[j]]))
                stop(paste0("In ",fname," wrong number of rows in ",
                            "\"Rho[[",j,"]]\".\n",call.=FALSE))
            if(warn) {
                whinge <- paste(wstrng,"has no row names.  I am assuming",
                                "that the\n rows of",wstrng,"correspond to the",
                                "sorted unique values of the",nmbr,
                                "component of \"y\".\n")
                warning(whinge)
            }
            nv <- as.numeric(yv)
            yv <- if(!any(is.na(nv))) yv[order(nv)] else sort(yv)
            rownames(Rho[[j]]) <- yv
        }
    }
    return(Rho)
}
# Bivariate dependent.
if(type==3) {
    for(j in 1:2) {
        nmbr <- if(j==1) "first" else "second"
        rn <- dimnames(Rho)[[j]]
        yv <- yval[[j]]
        if(!is.null(rn)) {
            if(!all(yv %in% rn))
                stop(paste0("In ",fname," some y values do not match the\n",
                           "names of the ",nmbr," dimension of \"Rho\".\n"),
                            call.=FALSE)
        } else {
            if(length(yv) != dim(Rho)[j])
                stop(paste0("In ",fname," the ",nmbr,
                           " dimension of \"Rho\" is wrong.\n",sep=""),call.=FALSE)
            if(warn) {
                whinge <- paste0("\"Rho\" has no",nmbr,"dimension names.\n",
                                "  I am assuming that the names of this",
                                "dimension correspond to the sorted unique values",
                                "of the",nmbr,"component of \"y\".\n")
                warning(whinge)
            }
            nv <- as.numeric(yv)
            yv <- if(!any(is.na(nv))) yv[order(nv)] else sort(yv)
            dimnames(Rho)[[j]] <- yv
        }
    }
    return(Rho)
}
stop(paste("The value",type,"of \"type\" is not recognised.\n"))
}
