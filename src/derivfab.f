C Output from Public domain Ratfor, version 1.03
      subroutine derivfab(y,fy,kstate,ashp,bshp,nbot,ntop,nd,d1a,d1b,d2a
     *a,d2ab,d2bb)
      implicit double precision(a-h,o-z)
      dimension fy(kstate)
      dimension ashp(kstate), bshp(kstate)
      dimension d1a(kstate), d1b(kstate)
      dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      xtop = dble(ntop)
      xbot = dble(nbot)
      sy = (y - xbot + one)/(xtop - xbot + two)
      t1 = log(sy)
      t2 = log(1-sy)
      s1 = (two - xbot)/(xtop - xbot + two)
      tt1 = log(s1)
      tt2 = log(1-s1)
      do23000 i = 1,kstate 
      alpha = ashp(i)
      beta = bshp(i)
      xam = alpha*tt1 + beta*tt2
      eee = zero
      d1ea = zero
      d1eb = zero
      if(nd.eq.2)then
      d2eaa = zero
      d2eab = zero
      d2ebb = zero
      endif
      do23004 k = nbot,ntop 
      sk = (dble(k)- xbot + one)/(xtop - xbot + two)
      tk1 = log(sk)
      tk2 = log(1-sk)
      xk = alpha*tk1 + beta*tk2
      if(xk .gt. xam)then
      xam = xk
      endif
23004 continue
23005 continue
      do23008 k = nbot,ntop 
      sk = (dble(k)- xbot + one)/(xtop - xbot + two)
      tk1 = log(sk)
      tk2 = log(1-sk)
      aitchk = one/(sk*(one - sk))
      eee = eee + aitchk * exp(alpha*tk1 + beta*tk2 - xam)
      d1ea = d1ea + aitchk * tk1*exp(alpha*tk1 + beta*tk2 - xam)
      d1eb = d1eb + aitchk * tk2*exp(alpha*tk1 + beta*tk2 - xam)
      if(nd.eq.2)then
      d2eaa = d2eaa + aitchk * tk1*tk1*exp(alpha*tk1 + beta*tk2 - xam)
      d2eab = d2eab + aitchk * tk1*tk2*exp(alpha*tk1 + beta*tk2 - xam)
      d2ebb = d2ebb + aitchk * tk2*tk2*exp(alpha*tk1 + beta*tk2 - xam)
      endif
23008 continue
23009 continue
      einv = one/eee
      pred1a = einv*d1ea
      pred1b = einv*d1eb
      d1a(i) = fy(i)*(t1 - pred1a)
      d1b(i) = fy(i)*(t2 - pred1b)
      if(nd.eq.2)then
      einv2 = einv**2
      pred2aa = einv*d2eaa - einv2*(d1ea**2)
      pred2ab = einv*d2eab - einv2*d1ea*d1eb
      pred2bb = einv*d2ebb - einv2*(d1eb**2)
      d2aa(i) = fy(i)*((t1 - pred1a)**2 - pred2aa)
      d2ab(i) = fy(i)*((t1 - pred1a)*(t2-pred1b) - pred2ab)
      d2bb(i) = fy(i)*((t2 - pred1b)**2 - pred2bb)
      endif
23000 continue
23001 continue
      return
      end
