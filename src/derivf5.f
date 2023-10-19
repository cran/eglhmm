C Output from Public domain Ratfor, version 1.03
      subroutine derivf5(y,phimat,tdm,kstate,npar,nxc,nyv,nd, d1a,d1b,d2
     *aa,d2ab,d2bb,d1f,d2f)
      implicit double precision(a-h,o-z)
      dimension phimat(nxc,nyv), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      dimension d1a(kstate), d1b(kstate)
      dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
      integer ell, r, s, t, u
      iy = idnint(y)
      if(npar .eq. nxc)then
      npro = 0
      else
      npro = kstate*(kstate-1)
      endif
      nyvm1 = nyv - 1
      do23002 k = 1,kstate 
      call pmf(iy,tdm(1,k),phimat,nyv,nxc,pmfy)
      do23004 r = 1,nxc 
      do23006 s = 1,nyvm1 
      j = npro + (r-1)*nyvm1 + s
      call pmf(s,tdm(1,k),phimat,nyv,nxc,pmfs)
      call delta(iy,s,iysd)
      part2 = iysd-pmfs
      xr = tdm(r,k)
      d1f(k,j) = pmfy*part2*xr
      if(nd.gt.1)then
      do23010 t = 1, nxc 
      do23012 u = 1,nyvm1 
      ell = npro + (t-1)*nyvm1 + u
      call pmf(u,tdm(1,k),phimat,nyv,nxc,pmfu)
      call delta(s,u,isud)
      call delta(iy,u,iyud)
      part3 = pmfs*pmfu - pmfs*isud
      part4 = (iysd - pmfs)*(iyud - pmfu)
      xt = tdm(t,k)
      d2f(k,j,ell) = pmfy*(part3 + part4)*xr*xt
23012 continue
23013 continue
23010 continue
23011 continue
      endif
23006 continue
23007 continue
23004 continue
23005 continue
23002 continue
23003 continue
      return
      end
