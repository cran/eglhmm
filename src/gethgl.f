C Output from Public domain Ratfor, version 1.03
      subroutine gethgl(nd,ndistr,fy,gmu,sd,lambda,p,ashp,bshp,phimat,y,
     *ymiss, tdm,nind,n,tpm,xispd,size,nbot,ntop,d1pi,d2pi,kstate, npar,
     *npt,nyv,nxc,d1p,d2p,d1f,d2f,d1a,d1b,d2aa,d2ab,d2bb, alpha,alphw,a,
     *b,aw,bw,xlc,hess)
      implicit double precision(a-h,o-z)
      integer ymiss(n), size
      dimension fy(kstate,n),gmu(kstate,n),sd(kstate,n)
      double precision lambda(kstate,n)
      dimension p(kstate,n)
      dimension ashp(kstate,n),bshp(kstate,n),phimat(nxc,nyv)
      dimension y(n),tdm(nxc,nind)
      dimension tpm(kstate,kstate),xispd(kstate)
      dimension d1pi(kstate,npt),d2pi(kstate,npt,npt)
      dimension d1p(kstate,kstate,npt),d2p(kstate,kstate,npt,npt)
      dimension d1f(kstate,npar),d2f(kstate,npar,npar)
      dimension d1a(kstate),d1b(kstate)
      dimension d2aa(kstate),d2ab(kstate),d2bb(kstate)
      dimension alpha(kstate),alphw(kstate)
      dimension a(kstate,npar),b(kstate,npar,npar)
      dimension aw(kstate,npar),bw(kstate,npar,npar)
      dimension xlc(n),hess(npar,npar)
      zero = 0.d0
      kt = 1
      if(nd .ge. 1)then
      kstart = 1
      call derivf(ndistr,y(kt),ymiss(kt),fy(1,kt),phimat,tdm(1,kstart), 
     *gmu(1,kt),sd(1,kt),lambda(1,kt),p(1,kt),ashp(1,kt), bshp(1,kt),kst
     *ate,npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b, d2aa,d2ab,d2bb,nd,d1f
     *,d2f)
      endif
      sxlc = zero
      kstpr = npt - npar
      do23002 j = 1,kstate 
      alpha(j) = xispd(j)*fy(j,kt)
      sxlc = sxlc + alpha(j)
      if(nd .ge. 1)then
      do23006 k1 = 1,npar 
      if(ymiss(1) .eq. 1)then
      d1fx1 = zero
      else
      d1fx1 = d1f(j,k1)
      endif
      a(j,k1) = xispd(j)*d1fx1 + fy(j,kt)*d1pi(j,kstpr+k1)
      if(nd .eq. 2)then
      do23012 k2 = 1,npar 
      if(ymiss(1) .eq. 1)then
      d1fx2 = zero
      else
      d1fx2 = d1f(j,k2)
      endif
      if(ymiss(1) .eq. 1)then
      d2fx = zero
      else
      d2fx = d2f(j,k1,k2)
      endif
      b(j,k1,k2) = xispd(j)*d2fx + d1pi(j,kstpr+k1)*d1fx2 + d1pi(j,kstpr
     *+k2)*d1fx1 + fy(j,kt)*d2pi(j,kstpr+k1,kstpr+k2)
23012 continue
23013 continue
      endif
23006 continue
23007 continue
      endif
23002 continue
23003 continue
      xlc(1) = sxlc
      do23018 j = 1,kstate 
      alpha(j) = alpha(j)/sxlc
23018 continue
23019 continue
      if(n.gt.1)then
      do23022 kt = 2,n 
      if(nd .ge. 1)then
      kstart = 1+(kt-1)*kstate
      call derivf(ndistr,y(kt),ymiss(kt),fy(1,kt),phimat,tdm(1,kstart), 
     *gmu(1,kt),sd(1,kt),lambda(1,kt),p(1,kt),ashp(1,kt), bshp(1,kt),kst
     *ate,npar,npt,nyv,nxc,size,nbot,ntop,d1a,d1b, d2aa,d2ab,d2bb,nd,d1f
     *,d2f)
      endif
      if(nd .eq. 2)then
      do23028 j = 1,kstate 
      do23030 k1 = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx1 = zero
      else
      d1fx1 = d1f(j,k1)
      endif
      do23034 k2 = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx2 = zero
      d2fx = zero
      else
      d1fx2 = d1f(j,k2)
      d2fx = d2f(j,k1,k2)
      endif
      vvv = zero
      xxx = zero
      yy1 = zero
      yy2 = zero
      zz1 = zero
      zz2 = zero
      www = zero
      do23038 i = 1,kstate 
      vvv = vvv + alpha(i)*d2p(i,j,kstpr+k1,kstpr+k2)
      xxx = (xxx + a(i,k1)*d1p(i,j,kstpr+k2) + a(i,k2)*d1p(i,j,kstpr+k1)
     * + b(i,k1,k2)*tpm(i,j))
      yy1 = yy1 + alpha(i)*d1p(i,j,kstpr+k2)
      yy2 = yy2 + a(i,k2)*tpm(i,j)
      zz1 = zz1 + alpha(i)*d1p(i,j,kstpr+k1)
      zz2 = zz2 + a(i,k1)*tpm(i,j)
      www = www + alpha(i)*tpm(i,j)
23038 continue
23039 continue
      vvv = fy(j,kt)*vvv
      xxx = fy(j,kt)*xxx/sxlc
      yyy = d1fx1*(yy1 + yy2/sxlc)
      zzz = d1fx2*(zz1 + zz2/sxlc)
      www = d2fx*www
      bw(j,k1,k2) = vvv + xxx + yyy + zzz + www
23034 continue
23035 continue
23030 continue
23031 continue
23028 continue
23029 continue
      do23040 j = 1,kstate 
      do23042 k1 = 1,npar 
      do23044 k2 = 1,npar 
      b(j,k1,k2) = bw(j,k1,k2)
23044 continue
23045 continue
23042 continue
23043 continue
23040 continue
23041 continue
      endif
      if(nd .ge. 1)then
      do23048 j = 1,kstate 
      do23050 k = 1,npar 
      if(ymiss(kt) .eq. 1)then
      d1fx = zero
      else
      d1fx = d1f(j,k)
      endif
      xxx = zero
      yyy = zero
      zzz = zero
      do23054 i = 1, kstate 
      xxx = xxx + alpha(i)*d1p(i,j,kstpr+k)
      yyy = yyy + a(i,k)*tpm(i,j)
      zzz = zzz + alpha(i)*tpm(i,j)
23054 continue
23055 continue
      aw(j,k) = fy(j,kt)*(xxx + yyy/sxlc) + d1fx*zzz
23050 continue
23051 continue
23048 continue
23049 continue
      do23056 j = 1,kstate 
      do23058 k = 1,npar 
      a(j,k) = aw(j,k)
23058 continue
23059 continue
23056 continue
23057 continue
      endif
      sxlc = zero
      do23060 j = 1,kstate 
      alphw(j) = zero
      do23062 i = 1,kstate 
      alphw(j) = alphw(j) + alpha(i)*tpm(i,j)
23062 continue
23063 continue
      alphw(j) = fy(j,kt)*alphw(j)
      sxlc = sxlc + alphw(j)
23060 continue
23061 continue
      xlc(kt) = sxlc
      do23064 j = 1,kstate 
      alpha(j) = alphw(j)/sxlc
23064 continue
23065 continue
23022 continue
23023 continue
      endif
      if(nd .eq. 2)then
      do23068 k1 = 1,npar 
      do23070 k2 = 1,npar 
      xxx = zero
      yyy = zero
      zzz = zero
      do23072 i = 1,kstate 
      xxx = xxx + b(i,k1,k2)
      yyy = yyy + a(i,k1)
      zzz = zzz + a(i,k2)
23072 continue
23073 continue
      hess(k1,k2) = (xxx - yyy*zzz/sxlc)/sxlc
23070 continue
23071 continue
23068 continue
23069 continue
      endif
      return
      end
