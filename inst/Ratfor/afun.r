subroutine afun(fy,xispd,tpm,epsilon,n,kstate,wrk,xlc,alpha)
implicit double precision(a-h,o-z)
dimension wrk(kstate), xispd(kstate), xlc(n)
dimension fy(kstate,n), tpm(kstate,kstate), alpha(kstate,n)

# Set some constants
one  = 1.d0
zero = 0.d0

# Set the value to give to the ``log-likelihood constant'', xlc(...)
# if this is indeterminate --- i.e. less than epsilon.
# Possible choices: -1, 1, or epsilon.
dummy = -one

# Update the initial alpha.
tsum = zero
do j = 1,kstate {
	wrk(j) =  fy(j,1)*xispd(j)
	tsum = tsum + wrk(j)
}

if(tsum < epsilon) {
	xlc(1) = dummy
	do j = 1,kstate {
		alpha(j,1) = one/kstate
	}
}
else {
	xlc(1) = tsum
	do j = 1,kstate {
		alpha(j,1) = wrk(j)/tsum
	}
}

# Run through the remaining n-1 of the alphas (recursing!).
do kt = 2,n {
	tsum = zero
	ktm = kt - 1
	do j = 1,kstate {
		wrk(j) = zero
		do i = 1,kstate {
			wrk(j) = wrk(j) + alpha(i,ktm)*tpm(i,j)
#if(kt==7247 & j==7) {
#    call intpr("i:",-1,i,1)
#    call dblepr("wrk(7):",-1,wrk(7),1)
#    call dblepr("alpha:",-1,alpha(i,ktm),1)
#    call dblepr("tpm:",-1,tpm(i,j),1)
#}
		}
		wrk(j) = fy(j,kt)*wrk(j)
		tsum = tsum + wrk(j)
	}
#if(kt==7247) call rexit("Oh fuck!")
	if(tsum < epsilon) {
		xlc(kt) = dummy
		do j = 1,kstate {
			alpha(j,kt) = one/kstate
		}
	}
	else {
		xlc(kt) = tsum
		do j = 1,kstate {
			alpha(j,kt) = wrk(j)/tsum
		}
	}
#if(kt >= 7246 & 7247 >= kt) {
#    call dblepr("tsum:",-1,tsum,1)
#    call dblepr("epsilon:",-1,epsilon,1)
#}
#if(kt==7247) call rexit("Oh fuck!")
#if(tsum > one + 1d-8) {
#    call intpr1("like. const. > 1 at t =",-1,kt)
#    call dblepr1("like. const. is:",-1,tsum)
#    call rexit("Oh fuck!")
#}
}

return
end
