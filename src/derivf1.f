C Output from Public domain Ratfor, version 1.03
      subroutine derivf1(y,gmu,sd,fy,tdm,kstate,npar,npt,nxc,nd,d1f,d2f,
     *kt)
      implicit double precision(a-h,o-z)
      dimension gmu(kstate), sd(kstate)
      dimension fy(kstate), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      logical sigfix
      if(npt .gt. npar)then
      ignore = 0
      else
      ignore = kstate*(kstate-1)
      endif
      sigfix = (npt .eq. kstate*(kstate-1) + nxc)
      if(sigfix)then
      jincr = 0
      else
      jincr = kstate
      endif
      zero = 0.d0
      one = 1.d0
      three = 3.d0
      do23004 i = 1,kstate 
      z = (y-gmu(i))/sd(i)
      dfdmu = fy(i)*z/sd(i)
      if(.not.sigfix)then
      dfdzeta = fy(i)*(z**2 - one)
      dfdsigma = dfdzeta/sd(i)
      d1f(i,ignore+i) = dfdzeta
      else
      dfdsigma = zero
      endif
      do23008 j = 1,nxc 
      d1f(i,ignore+jincr+j) = dfdmu*tdm(j,i)
23008 continue
23009 continue
      if(nd .gt. 1)then
      d2fdmu2 = dfdsigma/sd(i)
      do23012 j = 1,nxc 
      do23014 k = 1,nxc 
      d2f(i,ignore+jincr+j,ignore+jincr+k) = d2fdmu2*tdm(j,i)*tdm(k,i)
23014 continue
23015 continue
23012 continue
23013 continue
      if(sigfix)then
      goto 23004
      endif
      d2fdsigma2 = fy(i)*((z**2 - one)**2 + one - three*z**2)/sd(i)**2
      d2fdzeta2 = sd(i)*(dfdsigma + d2fdsigma2*sd(i))
      d2fdmudzeta = fy(i)*(z**2 - three)*z/sd(i)
      d2f(i,ignore+i,ignore+i) = d2fdzeta2
      do23018 k = 1,nxc 
      d2f(i,ignore+i,ignore+kstate+k) = d2fdmudzeta*tdm(k,i)
23018 continue
23019 continue
      do23020 j = 1,nxc 
      d2f(i,ignore+kstate+j,ignore+i) = d2fdmudzeta*tdm(j,i)
23020 continue
23021 continue
      endif
23004 continue
23005 continue
      return
      end
