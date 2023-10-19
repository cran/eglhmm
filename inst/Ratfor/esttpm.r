subroutine esttpm(ns,n,k,tpm,mixture,wrk)
implicit double precision(a-h,o-z)
dimension ns(n), tpm(k,k), wrk(k)

zero = 0.d0
one  = 1.d0
ook  = one/dble(k)

do i = 1,k {
	do j = 1,k {
		tpm(i,j) = zero
	}
}

do nt = 2,n {
	nb = nt-1
	do i = 1,k {
		do j = 1,k {
			if(ns(nb)==i&ns(nt)==j) tpm(i,j) = tpm(i,j)+one
		}
	}
}

if(mixture > 0) {
	den = zero
	do j = 1,k {
		wrk(j) = zero
		do i = 1,k {
			den = den + tpm(i,j)
			wrk(j) = wrk(j) + tpm(i,j)
		}
	}
	do i = 1,k {
		do j = 1,k {
			tpm(i,j) = wrk(j)/den
		}
	}
}

else {
	do i = 1,k {
		den = zero
		do j = 1,k {
			den = den + tpm(i,j)
		}
		if(den>=one) {
			do j = 1,k {
				tpm(i,j) = tpm(i,j)/den
			}
		}
		else {
			do j = 1,k {
				tpm(i,j) = ook
			}
		}
	}
}

return
end
