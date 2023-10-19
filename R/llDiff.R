llDiff <- function(new.ll,old.ll,tolerance) {
iii <- c(is.finite(old.ll),is.finite(new.ll))
jjj <- 1+sum((1:2)*iii)
switch(EXPR=jjj,NA,{stop("Infinite decrease in log likelihood.\n")},Inf,
                 (new.ll - old.ll)/(abs(old.ll) + tolerance))
}
