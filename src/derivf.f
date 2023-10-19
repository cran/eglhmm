C Output from Public domain Ratfor, version 1.03
      subroutine derivf(ndistr,y,ymiss,fy,phimat,tdm,gmu,sd,lambda,p,ash
     *p,bshp,kstate, npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b,d2aa,d2ab,d
     *2bb, nd,d1f,d2f)
      implicit double precision(a-h,o-z)
      dimension fy(kstate), phimat(nxc,nyv), tdm(nxc,kstate)
      dimension gmu(kstate), sd(kstate)
      double precision lambda(kstate)
      dimension p(kstate), ashp(kstate), bshp(kstate)
      dimension d1a(kstate), d1b(kstate)
      dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      integer ymiss
      integer size
      do23000 i = 1, kstate 
      do23002 j = 1, npar 
      d1f(i,j) = 0.d0
      do23004 k = 1, npar 
      d2f(i,j,k) = 0.d0
23004 continue
23005 continue
23002 continue
23003 continue
23000 continue
23001 continue
      if(ymiss .gt. 0)then
      return
      endif
      if(ndistr.eq.1)then
      call derivf1(y,gmu,sd,fy,tdm,kstate,npar,npt,nxc,nd,d1f,d2f)
      return
      endif
      if(ndistr.eq.2)then
      call derivf2(y,lambda,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
      return
      endif
      if(ndistr.eq.3)then
      call derivf3(y,p,size,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
      return
      endif
      if(ndistr.eq.4)then
      call derivf4(y,ashp,bshp,nbot,ntop,fy,tdm,kstate,npar, nxc,nd,d1f,
     *d2f,d1a,d1b,d2aa,d2ab,d2bb)
      return
      endif
      if(ndistr.eq.5)then
      call derivf5(y,phimat,tdm,kstate,npar,nxc,nyv,nd, d1a,d1b,d2aa,d2a
     *b,d2bb,d1f,d2f)
      return
      endif
      call intpr1("The value of ndistr is",-1,ndistr)
      call rexit("This value should be between 1 and 5 inclusive; bailin
     *g out.")
      end
