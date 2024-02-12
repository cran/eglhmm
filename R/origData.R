origData <- function(rData) {
#
# Function to recover the original data from the replicated data.
# If "distr" was "Multinom" in the call to eglhmm() that created
# rData, then the response variable(s) may gave been coerced into
# factors.  This coersion is not undone, so the object ("oD") that
# is returned may be not *quite* exactly the same as the original
# data.
#
state <- rData$state
if(is.null(state)) return(rData) # K = 1; no relication.
oD           <- rData[rData$state==1,]
oD$state     <- NULL
oD[["cf"]]   <- NULL
rownames(oD) <- 1:nrow(oD)
oD
}
