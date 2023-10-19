subroutine gfun(alpha,beta,epsilon,n,kstate,wrk,gamma)
implicit double precision(a-h,o-z)
dimension alpha(kstate,n), beta(kstate,n), gamma(kstate,n)
dimension wrk(kstate)

zero = 0.d0
ook  = 1.d0/dble(kstate)

do kt = 1,n {
	tsum = zero
	do i = 1,kstate {
		wrk(i) = alpha(i,kt)*beta(i,kt)
		tsum = tsum + wrk(i)
	}
	if(tsum<epsilon) {
		do i = 1,kstate {
			gamma(i,kt) = ook
		}
	}
	else {
		do i = 1,kstate {
			gamma(i,kt) = wrk(i)/tsum
		}
	}
}

return
end
