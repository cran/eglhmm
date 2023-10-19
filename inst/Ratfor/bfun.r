subroutine bfun(fy,xispd,tpm,epsilon,n,kstate,wrk,beta)
implicit double precision(a-h,o-z)
dimension wrk(kstate), xispd(kstate)
dimension fy(kstate,n), tpm(kstate,kstate), beta(kstate,n)

# Set some constants.
one  = 1.d0
zero = 0.d0

# Set the last beta's.

do j = 1,kstate {
	beta(j,n) = one # Doesn't matter that the beta_n's are
}                       # not rescaled.

# Run through the remaining n-1 of the betas (recursing!), backwards.
do ktb = 2,n {
	kt  = n - ktb + 1
	ktp = kt + 1
	tsum = zero
	do i = 1,kstate {
		wrk(i) = zero
		do j = 1,kstate {
			wrk(i) = wrk(i) + tpm(i,j)*beta(j,ktp)*fy(j,ktp)
		}
		tsum = tsum + wrk(i)
	}
	if(tsum < epsilon) {
		do j = 1,kstate {
			beta(j,kt) = one/kstate
		}
	}
	else {
		do j = 1,kstate {
			beta(j,kt) = wrk(j)/tsum
		}
	}
}

return
end
