C Output from Public domain Ratfor, version 1.03
      subroutine derivf2(y,lambda,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
      implicit double precision(a-h,o-z)
      double precision lambda(kstate)
      dimension fy(kstate), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      one = 1.d0
      if(npar .eq. nxc)then
      npro = 0
      else
      npro = kstate*(kstate-1)
      endif
      do23002 i = 1,kstate 
      fctr = y/lambda(i) - one
      dfdlam = fy(i)*fctr
      do23004 j = 1,nxc 
      d1f(i,npro+j) = dfdlam*lambda(i)*tdm(j,i)
      if(nd .gt. 1)then
      d2fdlam2 = fy(i)*(fctr**2 - y/lambda(i))
      do23008 k = 1,nxc 
      xxt = tdm(j,i)*tdm(k,i)
      d2f(i,npro+j,npro+k) = (dfdlam*lambda(i) + d2fdlam2*lambda(i)**2)*
     *xxt
23008 continue
23009 continue
      endif
23004 continue
23005 continue
23002 continue
23003 continue
      return
      end
