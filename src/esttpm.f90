! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine esttpm(ns,n,k,tpm,mixture,wrk)
implicit double precision(a-h,o-z)
dimension :: ns(n), tpm(k,k), wrk(k)

zero = 0.d0
one  = 1.d0
ook  = one/dble(k)

do i = 1,k
    do j = 1,k
        tpm(i,j) = zero
    enddo
enddo

do nt = 2,n
    nb = nt-1
    do i = 1,k
        do j = 1,k
            if(ns(nb)==i .and. ns(nt)==j) tpm(i,j) = tpm(i,j)+one
        enddo
    enddo
enddo

if(mixture > 0) then
    den = zero
    do j = 1,k
        wrk(j) = zero
        do i = 1,k
            den = den + tpm(i,j)
            wrk(j) = wrk(j) + tpm(i,j)
        enddo
    enddo
    do i = 1,k
        do j = 1,k
            tpm(i,j) = wrk(j)/den
        enddo
    enddo
else
    do i = 1,k
        den = zero
        do j = 1,k
            den = den + tpm(i,j)
        enddo
        if(den>=one) then
            do j = 1,k
                tpm(i,j) = tpm(i,j)/den
            enddo
        else
            do j = 1,k
                tpm(i,j) = ook
            enddo
        endif
    enddo
endif

end subroutine esttpm
