C Output from Public domain Ratfor, version 1.03
      subroutine derivf4(y,ashp,bshp,nbot,ntop,fy,tdm,kstate,npar, nxc,n
     *d,d1f,d2f,d1a,d1b,d2aa,d2ab,d2bb)
      implicit double precision(a-h,o-z)
      dimension fy(kstate), ashp(kstate), bshp(kstate), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      dimension d1a(kstate), d1b(kstate)
      dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
      call derivfab(y,fy,kstate,ashp,bshp,nbot,ntop,nd,d1a,d1b,d2aa,d2ab
     *,d2bb)
      if(npar .eq. 2*nxc)then
      npro = 0
      else
      npro = kstate*(kstate-1)
      endif
      do23002 i = 1,kstate 
      do23004 j = 1,nxc 
      d1f(i,npro+j) = d1a(i)*tdm(j,i)
      d1f(i,npro+nxc+j) = d1b(i)*tdm(j,i)
23004 continue
23005 continue
23002 continue
23003 continue
      if(nd .gt. 1)then
      do23008 i = 1,kstate 
      do23010 j = 1,nxc
      do23012 k = 1,nxc 
      d2f(i,npro+j,npro+k) = d2aa(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+j,npro+nxc+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+nxc+j,npro+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+nxc+j,npro+nxc+k) = d2bb(i)*tdm(j,i)*tdm(k,i)
23012 continue
23013 continue
23010 continue
23011 continue
23008 continue
23009 continue
      endif
      return
      end
