! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine afun(fy,xispd,tpm,epsilon,n,kstate,wrk,xlc,alpha)
implicit double precision(a-h,o-z)
dimension :: wrk(kstate), xispd(kstate), xlc(n)
dimension :: fy(kstate,n), tpm(kstate,kstate), alpha(kstate,n)

! Set some constants
one  = 1.d0
zero = 0.d0

! Set the value to give to the ``log-likelihood constant'', xlc(...)
! if this is indeterminate --- i.e. less than epsilon.
! Possible choices: -1, 1, or epsilon.
dummy = -one

! Update the initial alpha.
tsum = zero
do j = 1,kstate 
    wrk(j) =  fy(j,1)*xispd(j)
    tsum = tsum + wrk(j)
enddo

if(tsum < epsilon) then
    xlc(1) = dummy
    do j = 1,kstate
        alpha(j,1) = one/kstate
    enddo
else
    xlc(1) = tsum
    do j = 1,kstate
        alpha(j,1) = wrk(j)/tsum
    enddo
endif

! Run through the remaining n-1 of the alphas (recursing!).
do kt = 2,n
    tsum = zero
    ktm = kt - 1
    do j = 1,kstate
        wrk(j) = zero
        do i = 1,kstate
            wrk(j) = wrk(j) + alpha(i,ktm)*tpm(i,j)
        enddo
        wrk(j) = fy(j,kt)*wrk(j)
        tsum = tsum + wrk(j)
    enddo
    if(tsum < epsilon) then
        xlc(kt) = dummy
        do j = 1,kstate
            alpha(j,kt) = one/kstate
        enddo
    else
        xlc(kt) = tsum
        do j = 1,kstate
            alpha(j,kt) = wrk(j)/tsum
        enddo
    endif
enddo

end subroutine afun
