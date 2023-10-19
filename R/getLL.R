getLL <- function(rp) {
    if(any(rp$llc <= 0)) -Inf else sum(log(rp$llc))
}
