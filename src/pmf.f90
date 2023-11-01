! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine pmf(y,x,phi,nyv,nxc,pmfy)
implicit double precision(a-h,o-z)
integer :: y
dimension :: x(nxc), phi(nxc,nyv)
if(y < 1 .or. y > nyv) then
    call intpr1("The value of y is:",-1,y)
    call rexit("This value is out of bounds.\n")
endif
nyvm1 = nyv-1 ! All values of phi(i,nyv), i = 1, ..., nxc, are 0.
zed = 1.d0
top = 1.d0
do j = 1,nyvm1
    some = 0.d0
    do i = 1,nxc
        some = some + x(i)*phi(i,j)
    enddo
    esome = exp(some)
    if(j==y) top = esome
    zed = zed + esome
enddo
pmfy = top/zed
end subroutine pmf
