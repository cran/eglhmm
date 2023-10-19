C Output from Public domain Ratfor, version 1.03
      subroutine derivf3(y,p,size,fy,tdm,kstate,npar,nxc,nd,d1f,d2f)
      implicit double precision(a-h,o-z)
      double precision p(kstate)
      integer size
      dimension fy(kstate), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      one = 1.d0
      if(npar .eq. nxc)then
      npro = 0
      else
      npro = kstate*(kstate-1)
      endif
      do23002 i = 1,kstate 
      fctr = (y/p(i) - (size - y)/(one - p(i)))
      dfdp = fy(i)*fctr
      u = dlog(p(i)/(one-p(i)))
      emu = dexp(-u)
      h1 = emu/(one+emu)**2
      do23004 j = 1,nxc 
      d1f(i,npro+j) = dfdp*h1*tdm(j,i)
      if(nd .gt. 1)then
      h2 = h1*(emu-one)/(emu+one)
      d2fdp2 = fy(i)*(fctr**2 - y/p(i)**2 - (size - y)/(one - p(i))**2)
      do23008 k = 1,nxc 
      d2f(i,npro+j,npro+k) = (dfdp*h2 + d2fdp*h1**2)*tdm(j,i)*tdm(k,i)
23008 continue
23009 continue
      endif
23004 continue
23005 continue
23002 continue
23003 continue
      return
      end
