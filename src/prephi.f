C Output from Public domain Ratfor, version 1.03
      subroutine prephi(ndistr,y,fy,kstate,zeta,sigma,size,d1zeta,d2zeta
     *,d1u,d2u,nd)
      implicit double precision(a-h,o-z)
      integer size
      dimension fy(kstate), zeta(kstate), sigma(kstate)
      dimension d1zeta(kstate), d2zeta(kstate), d1u(kstate),d2u(kstate)
      zero = 0.d0
      one = 1.d0
      half = 0.5d0
      if(ndistr.eq.1)then
      do23002 k = 1,kstate 
      d1zeta(k) = fy(k)*(y-zeta(k))/sigma(k)**2
      d1u(k) = one
      if(nd .gt. 1)then
      d2zeta(k) = (fy(k)/sigma(k)**2)*((y-zeta(k))**2/sigma(k)**2 - one)
      d2u(k) = zero
      endif
23002 continue
23003 continue
      else
      if(ndistr.eq.2)then
      do23008 k = 1,kstate 
      d1zeta(k) = fy(k)*(y/zeta(k) - one)
      d1u(k) = zeta(k)
      if(nd .gt. 1)then
      d2zeta(k) = fy(k)*((y/zeta(k) - one)**2 - y/zeta(k)**2)
      d2u(k) = zeta(k)
      endif
23008 continue
23009 continue
      else
      if(ndistr.eq.3)then
      do23014 k = 1,kstate 
      d1zeta(k) = fy(k)*(y/zeta(k) - (size-y)/(1-zeta(k)))
      u = dlog(zeta(k)/(1-zeta(k)))
      emu = dexp(-u)
      d1u(k) = emu/(1+emu)**2
      if(nd .gt. 1)then
      d2zeta(k) = fy(k)*((y/zeta(k) - (size-y)/(1-zeta(k)))**2 - (size-y
     *)/(1-zeta(k))**2 - y/zeta(k)**2)
      d2u(k) = emu*(emu-1)/(1+emu)**3
      endif
23014 continue
23015 continue
      else
      if(ndistr.eq.5)then
      do23020 k = 1,kstate 
      d1u(k) = 0.d0
      d2u(k) = 0.d0
      d1zeta(k) = 0.d0
      d2zeta(k) = 0.d0
23020 continue
23021 continue
      endif
      endif
      endif
      endif
      return
      end
