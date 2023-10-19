plot.eglhmm <- local({

rbox <- function() {
    if(is.null(dev.list()))
        stop("No device open.\n")
    uuu <- par()$usr
    lines(c(uuu[1],uuu[1]),c(0,uuu[4]),lwd=2)
    lines(c(uuu[1],uuu[2]),c(uuu[4],uuu[4]),lwd=2)
    lines(c(uuu[2],uuu[2]),c(0,uuu[4]),lwd=2)
    lines(c(uuu[1],uuu[2]),c(0,0),lwd=1)
    invisible()
}

function(x,...,wcells=NULL,col="red",nrnc=NULL,ntop=NULL,xlab=NULL,
         ylab=NULL,xlim=NULL,ylim=NULL,main=NULL,cex.main=1.5) {
if(is.null(xlab)) xlab <- "x"
K  <- nrow(x$tpm)
cf <- x$data$cf
if(is.null(wcells)) {
    wcells <- levels(cf)
}
if(!is.null(main)) main <- rep(main,length=length(wcells))
if(is.null(nrnc)) {
    nr <- ceiling(sqrt(K))
    nc <- ceiling(K/nr)
    nrnc <- c(nr,nc)
} else {
    if(nrnc[1]*nrnc[2] < K)
       stop("Specified dimensions for the array of plots is too small.\n")
}
OP <- par(mfrow=nrnc,oma=c(0,0,3,0))
on.exit(par(OP))
nwc    <- length(wcells)
devint <- dev.interactive() | (interactive() & is.null(dev.list()))
distr  <- x$distr
size   <- x$size
switch(EXPR=distr,
    Gaussian = {
        mu <- split(x$mu,cf) # I think this needs fixed.
        if(any(sapply(mu,function(x){length(unique(x)) != K}))) {
            whinge <- paste0("The number of unique values of \"mu\" is not",
                             " equal to the number of states for some cells.\n")
            stop(whinge)
        }
        mu    <- lapply(mu,unique)
        sigma <- x$sigma
        tipe  <- "l"
        if(is.null(ylab)) ylab <- "probability density"
    },
    Poisson = {
        lambda <- split(x$lambda,cf)
        if(any(sapply(lambda,function(x){length(unique(x)) != K}))) {
            whinge <- paste0("The number of unique values of \"lambda\" is not",
                             " equal to the number of states for some cells.\n")
            stop(whinge)
        }
        lambda <- lapply(lambda,unique)
        tipe  <- "h"
        if(is.null(ntop)) {
            lamall <- unlist(lambda)
            ntop <- max(sapply(lamall,function(x){qpois(1e-7,x,lower.tail=FALSE)}))
        }
        xi <- 0:ntop
        if(is.null(xlim)) {
            xlim <- c(0,ntop)
        } else {
           xi <- xi[xlim[1] <= xi & xi <= xlim[2]]
        }
        if(is.null(ylab)) ylab <- "probability"
    },
    Binomial = {
        p <- split(x$p,cf)
        if(any(sapply(p,function(x){length(unique(x)) != K}))) {
            whinge <- paste0("The number of unique values of \"p\" is not",
                             " equal to the number of states for some cells.\n")
            stop(whinge)
        }
        p <- lapply(p,unique)
        tipe  <- "h"
        xi <- 0:size
        if(is.null(xlim)) {
            xlim <- c(0,size)
        } else {
           xi <- xi[xlim[1] <= xi & xi <= xlim[2]]
        }
        if(is.null(ylab)) ylab <- "probability"
    }
)

fval  <- vector("list",K*nwc)
ik <- 0
for(i in 1:nwc) {
    for(k in 1:K) {
        ik <- ik+1
        fval[[ik]] <- switch(EXPR=distr,
            Gaussian = {
                muik <- mu[[i]][k]
                sdk  <- sigma[k]
                a    <- muik - 3*sdk
                b    <- muik + 3*sdk
                xi   <- seq(a,b,length=100)
                dnorm(xi,mean=mu[[i]][k],sd=sigma[k])
            },
            Poisson = {
                dpois(xi,lambda[[i]][k])
            },
            Binomial = {
                dbinom(xi,p[[i]][k],size=size)
            }
        )
    }
}
if(is.null(ylim)) {
    ylim <- c(0,max(unlist(fval)))
}
ik <- 0
for(i in 1:nwc) {
    for(k in 1:K) {
        ik <- ik+1
        fk <- fval[[ik]]
        plot(xi,fk,type=tipe,xlim=xlim,ylim=ylim,main=paste0("state ",k),
             col=col,xlab=xlab,ylab=ylab,axes=FALSE)
        axis(side=2,lwd=0,lwd.ticks=1)
        axis(side=1,lwd=0,pos=0)
        rbox()
    }
    if(is.null(main)) {
        mane <- paste0("Cell ",i)
    } else {
        mane <- main[i]
    }
    if(mane!="") mtext(side=3,outer=TRUE,text=mane,cex=cex.main)
    if(devint & i < nwc) readline("Go? ")
}
}
})
