crossval <- local({

crossvalEngine <- function(trnDat,valDat,model,par0=NULL,crossVerb,repl,...) {
    if(crossVerb) {
        part1 <- "Trying to fit model to training data"
        part2 <- if(is.null(repl)) ".\n" else paste0(", replicate ",repl,".\n")
        cat(paste0(part1,part2))
    }
    if(is.null(par0)) par0 <- model[names(model$par0)]
    aarghs  <- c(list(object=model,par0=par0,data=origData(trnDat)),list(...))
    fit.trn <- do.call(update,aarghs)
    OK      <- if(is.na(fit.trn$converged)) TRUE else fit.trn$converged
    if(OK){
        if(crossVerb) {
            part1 <- "Model fitted to training data"
            part2 <- if(is.null(repl)) ".\n" else paste0(", replicate ",repl,".\n")
            cat(paste0(part1,part2))
            cat("Calculating log likelihood of validation data.\n")
        }
        if(length(levels(fit.trn$data$state)) <= 1) {
            rslt <- llKeq1(data=valDat,object=model)
        } else {
            rslt <- with(fit.trn,getHgl(nd=0,distr,theta,data=valDat,
                         fmla=formula,size=model$size,nbot=model$nbot,
                         ntop=model$ntop))
            rp <- recurse(fit.trn$fy,fit.trn$tpm)
        }
    } else {
        lastPar <- fit.trn[names(fit.trn$par0)]
        attr(lastPar,"trnDat") <- trnDat
        attr(lastPar,"valDat") <- valDat
        rslt <- NA
        if(crossVerb) {
            part1 <- "Failed to fit model to training data"
            part2 <- if(is.null(repl)) ".\n" else paste0(", replicate ",repl,".\n")
            cat(paste0(part1,part2))
            whinge <- paste0("If convergence failure seems to be due to a",
                             " need to\n","  increase itmax, then it is",
                             " probably advisable to do so!!!\n")
            cat(whinge)
        }
        attr(rslt,"lastPar") <- lastPar
    }
    if(crossVerb & !is.null(repl)) cat("Replicate",repl,"finished.\n")
    rslt
}

function(model,data=NULL,nrep,frac=0.8,type,id="id",minNcomp=100,
         seed=NULL,crossVerb=FALSE,lastPar=NULL,...) {

if(is.null(lastPar)) {
    if(missing(nrep)) {
        stop("When \"lastPar\" is NULL, \"nrep\" must be specified.\n")
    }
} else {
    nrep <- 1
}

# Set the seeds.
if(is.null(seed)) seed <- sample(1:1e6,1)
set.seed(seed)
if(nrep>1) SEEDS <- sample(1:1e6,nrep)

# If lastPar is supplied, use it.
if(!is.null(lastPar)) {
   trnDat <- attr(lastPar,"trnDat")
   valDat <- attr(lastPar,"valDat")
   cvll   <- crossvalEngine(trnDat=trnDat,valDat=valDat,model=model,par0=lastPar,
                              crossVerb=crossVerb,repl=NULL,...)
   attr(cvll,"seed") <- seed
   return(cvll)
} 

# Else .... get the starting values of the parameters.
if(identical(model$tpm,NA)) { # K = 1.
    par0 <- NULL
} else {
    par0 <- model[names(model$par0)]
}

# Get the observations.
if(is.null(data)) data <- model[["data"]]
if(is.null(data)) stop("No observations supplied.\n")

# Do the cross validation calculations.
cvll <- vector("list",nrep)

# For type=1, the training data are formed by removing "all but"
# a random fraction frac of the "complete" (original, not replicated
# over state) data set.  Here the complete data set consists of
# those rows of the original data where the emissions are *NOT*
# missing (not equal to NA).  Let the number of such rows be ntot.
# "Removing" data means setting emissions (which were not previously
# NA) equal to NA.  The number of such emissions set equal to NA is
# nlvo = (1 - frac)*ntot (rounded).  The number of training data
# is ntrn = ntot - nlvo.  The validation data are the complement
# (in the "complete data set") of the training data.  To form the
# validation data we take the complement of the indices of the
# training data and set the corresponding emissions equal to NA.
# After the training and validation date are thus formed, they are
# replicated over state, and the results are passed to crossvalEngine()
#
# For type=2, the training data are formed by selecting a random
# fraction "frac" of the levels of "id" where id is a factor
# identifying the individual components of the data. These
# components are independent time series, each modelled by the
# model in question.  Thus if NC = the length of the levels of the
# id factor, we take ntrn = frac*NC (rounded).  The training data
# consists of those components of the data for which the id level
# is in the selected fraction.  The number of such components is ntrn.
# The validation data consist of all components that are *NOT*
# in the training data.  The number of such components is NC - ntrn.
#
# For type=1, ntrn is the count of the number of observations in
# the training data, for each state.  For type=2, ntrn is the count
# of the number of components (time series) in the training data.

if(length(levels(data$state)) < 1) {
    data$state <- factor(rep(1,nrow(data)))
}
if(type==1) {
    ynm  <- as.character(model$formula[2])
    yvec <- data[[ynm]][data$state==1]
    ina  <- which(is.na(yvec))  # Is NA.
    inna <- which(!is.na(yvec)) # Is not NA.
    nlvo <- round((1-frac)*length(inna)) # lvo == leave out
    sdat <- split(data,f=data$state)
    K    <- length(levels(data$state))
    ii   <- as.vector(matrix(1:nrow(data),nrow=K,byrow=TRUE))
} else {
    idnm <- id
    id   <- data[[idnm]]
    if(is.null(id)) {
        whinge <- paste0("When \"type\" is 2, \"data\" must have a component",
                         " with the name specified\n  by argument \"id\".\n")
        stop(whinge)
    }
    NC <- length(levels(id))
    if(NC < minNcomp)
        stop("Too few observation components to use type=2.\n")
    ntrn <- round(frac*NC)
}

for(repl in 1:nrep) {
    if(nrep>1) set.seed(SEEDS[repl])
    if(type==1) {
        tlvo <- sample(inna,nlvo)
        vlvo <- setdiff(inna,tlvo)
        tyvec <- yvec
        tyvec[tlvo] <- NA
        tsdat  <- lapply(sdat,function(x,ynm){x[,ynm] <- tyvec; x},ynm=ynm)
        trnDat <- do.call(rbind,tsdat)
        trnDat <- trnDat[ii,]
        rownames(trnDat) <- 1:nrow(trnDat)
        vyvec <- yvec
        vyvec[vlvo] <- NA
        vsdat  <- lapply(sdat,function(x,ynm){x[,ynm] <- vyvec; x},ynm=ynm)
        valDat <- do.call(rbind,vsdat)
        valDat <- valDat[ii,]
        rownames(valDat) <- 1:nrow(valDat)
    } else {
        itrn   <- sample(1:NC,ntrn)
        ltrn   <- levels(id)[itrn]
        ival   <- setdiff(1:NC,itrn)
        lval   <- levels(id)[ival]
        trnDat <- data[id%in%ltrn,]
        valDat <- data[id%in%lval,]
        trnDat[,idnm] <- factor(trnDat[,idnm])
        valDat[,idnm] <- factor(valDat[,idnm])
    }
    cvll[[repl]] <- crossvalEngine(trnDat,valDat,model,crossVerb=crossVerb,
                                   repl=if(nrep==1) NULL else repl,...)
    if(nrep > 1) attr(cvll[[repl]],"seed") <- SEEDS[repl]
}
if(nrep==1) {
    cvll <- cvll[[1]]
    attr(cvll,"seed") <- seed
} else {
    seeds <- c(SEEDS,seed)
    names(seeds) <- c(paste0("seed.",1:nrep),"overAllSeed")
    attr(cvll,"seeds") <- seeds 
}
cvll
}
})
