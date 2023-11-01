! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine xfun(alpha,beta,fy,tpm,epsilon,n,kstate,wrk,xi,xisum)
implicit double precision(a-h,o-z)
dimension :: alpha(kstate,n), beta(kstate,n), fy(kstate,n)
dimension :: tpm(kstate,kstate), wrk(kstate,kstate), xi(kstate,kstate,n)
dimension :: xisum(kstate,kstate)

one  = 1.d0
zero = 0.d0
dns2 = dble(kstate*kstate)

do ktp = 2,n
    kt = ktp - 1
    tsum = zero
    do i = 1,kstate
        do j = 1, kstate
            wrk(i,j) = alpha(i,kt)*fy(j,ktp)*beta(j,ktp)*tpm(i,j)
            tsum = tsum + wrk(i,j)
        enddo
    enddo
    if(tsum<epsilon) then
        do i = 1,kstate
            do j = 1,kstate
                xi(i,j,kt) = one/dns2
            enddo
        enddo
    else
        do i = 1,kstate
            do j = 1,kstate
                xi(i,j,kt) = wrk(i,j)/tsum
            enddo
        enddo
    endif
enddo

! Sum up the xi's.
nm1 = n - 1
do i = 1,kstate 
    do j = 1,kstate
        xisum(i,j) = zero
        do k = 1,nm1
            xisum(i,j) = xisum(i,j) + xi(i,j,k)
        enddo
    enddo
enddo

end subroutine xfun
