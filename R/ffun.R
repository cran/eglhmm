ffun <- local({
# This function calculates the probabilty of each observation,
# given each corresponding state, when the distributions of the
# observations are Multinom (discnp).
# 
# Note that "data" is a data frame with column(s) constituting the
# response(s) (with name(s) given by "response" or "ynm" (see below)
# and other columns which include the predictors.  These predictors
# will include, at least, those determined by the state factor.

f1 <- function(data,fmla,response,Rho) {
# Univariate.
# Argument "response" is ignored; the response name "ynm" is
# extracted from "fmla".
    X   <- model.matrix(fmla[-2],data=data)
    M   <- X%*%Rho
    P   <- apply(M,1,expForm2p)
    ynm <- as.character(fmla[2])
    P[cbind(data[[ynm]],1:nrow(data))]
}

f2 <- function(data,fmla,response,Rho) {
# Bivariate independent.
# Here data is a data frame with columns comprising the bivariate
# response.  No auxiliary predictors are allowed, so if K > 1 then
# the formula must have ~0+state as its right hand side, and if K==1
# it must simply have ~1 as its right hand side.  Argument Rho is a
# list with two entries, one for each response.

    K <- nrow(Rho[[1]]) # Equal to nrow(Rho[[2]]) also!
    if(K > 1) {
        fmla1 <- as.formula(paste0(response[1],"~0+state"))
        fmla2 <- as.formula(paste0(response[2],"~0+state"))
    } else {
        fmla1 <- as.formula(paste0(response[1],"~1"))
        fmla2 <- as.formula(paste0(response[2],"~1"))
    }
    poot1 <- f1(data,fmla=fmla1,response=NULL,Rho=Rho[[1]])
    poot1[is.na(poot1)] <- 1
    poot2 <- f1(data,fmla=fmla2,response=NULL,Rho=Rho[[2]])
    poot2[is.na(poot2)] <- 1
    poot1*poot2
}

f3 <- function(data,fmla,response,Rho) {
# Bivariate independent.
# Here data is a data frame with columns comprising the bivariate
# response.  Argument Rho is a 3-dimensional array.
# Argument "fmla" is ignored.
K    <- dim(Rho)[3]
iii <- if(K > 1) which(data$state == 1) else 1:nrow(data)
data <- as.matrix(data[,response])[iii,] # Character matrix
lll  <- apply(data,1,is.na)
if(!any(lll)) {
    xxx <- lapply(1:nrow(data),function(i,Rho,y){Rho[y[i,1],y[i,2],]},
                  Rho=Rho,y=data)
} else {
    xxx <- lapply(1:nrow(data),function(j,lll,Rho,y){
        u <- lll[,j]
        v <- y[j,]
        nat <- switch(1+sum(u),1,1+which(u),4)
        switch(nat,
            Rho[v[1],v[2],],
            apply(Rho[,v[2],,drop=FALSE],3,sum),
            apply(Rho[v[1],,,drop=FALSE],3,sum),
            rep(1,dim(Rho)[3])
        )
    },lll=lll,Rho=Rho,y=data)
}
as.vector(matrix(unlist(xxx),nrow=dim(Rho)[3]))
}

f  <- list(f1,f2,f3)

function(data,fmla,response,Rho,type) {
#
# Function ffun to calculate
#        f(y) = Pr(Y=y | the model parameters and predictors)
# when the distributions of the observations are Multinom (discnp).
# In the univariate setting data has a response column which
# is a factor constituting the observations.  The name of this
# column is determined by as.character(fmla[2]).  Other columns
# constititute predictors. In the bivariate setting, data has
# columns, with names given by "response", which are factors
# constituting the (bivariate) observations.  In the bivariate
# setting there can be no auxiliary predictors.
#
# The returned result, fy, is a numeric vector (of probabilities).
#
# The "type" argument:
# type = 1 <--> univariate             -- f1
# type = 2 <--> bivariate, independent -- f2
# type = 3 <--> bivariate, dependent   -- f3

fy <- f[[type]](data=data,fmla=fmla,response,Rho=Rho)
fy
}
})
