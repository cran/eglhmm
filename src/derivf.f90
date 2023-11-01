! Manually recoded from Ratfor to Fortran 90, 29/10/2013.
subroutine derivf(ndistr,y,ymiss,fy,phimat,tdm,gmu,sd,lambda,p,ashp,bshp,kstate,&
                   npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b,d2aa,d2ab,d2bb,&
                   nd,d1f,d2f)
! Argument ndistr determines the distribution being used.
! 1 <--> Gaussian.
! 2 <--> Poisson.
! 3 <--> Binomial.
! 4 <--> Dbd.
! 5 <--> Multinom.
!
! nd = number of derivatives to calculate.
!
implicit double precision(a-h,o-z)
dimension :: fy(kstate), phimat(nxc,nyv), tdm(nxc,kstate)
dimension :: gmu(kstate), sd(kstate)
double precision :: lambda(kstate)
dimension :: p(kstate), ashp(kstate), bshp(kstate)
dimension :: d1a(kstate), d1b(kstate)
dimension :: d2aa(kstate), d2ab(kstate), d2bb(kstate)
dimension :: d1f(kstate,npar), d2f(kstate,npar,npar)
integer :: ymiss, size

do i = 1, kstate
    do j = 1, npar
        d1f(i,j) = 0.d0
        do k = 1, npar
            d2f(i,j,k) = 0.d0
        enddo
    enddo
enddo

if(ymiss > 0) return

if(ndistr==1) then
    call derivf1(y,gmu,sd,fy,tdm,kstate,npar,npt,nxc,nd,d1f,d2f)
    return
endif
if(ndistr==2) then
    call derivf2(y,lambda,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
    return
endif
if(ndistr==3) then
    call derivf3(y,p,size,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
    return
endif
if(ndistr==4) then
    call derivf4(y,ashp,bshp,nbot,ntop,fy,tdm,kstate,npar,&
                 nxc,nd,d1f,d2f,d1a,d1b,d2aa,d2ab,d2bb)
    return
endif
if(ndistr==5) then
    call derivf5(y,phimat,tdm,kstate,npar,nxc,nyv,nd,d1f,d2f)
    return
endif
call intpr1("The value of ndistr is",-1,ndistr)
call rexit("This value should be between 1 and 5 inclusive; bailing out.")
end subroutine derivf
