! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine gfun(alpha,beta,epsilon,n,kstate,wrk,gamma)
implicit double precision(a-h,o-z)
dimension :: alpha(kstate,n), beta(kstate,n), gamma(kstate,n)
dimension :: wrk(kstate)

zero = 0.d0
ook  = 1.d0/dble(kstate)

do kt = 1,n
    tsum = zero
    do i = 1,kstate
        wrk(i) = alpha(i,kt)*beta(i,kt)
        tsum = tsum + wrk(i)
    enddo
    if(tsum<epsilon) then
        do i = 1,kstate
            gamma(i,kt) = ook
        enddo
    else
        do i = 1,kstate
            gamma(i,kt) = wrk(i)/tsum
        enddo
    endif
enddo

end subroutine gfun
