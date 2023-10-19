C Output from Public domain Ratfor, version 1.03
      subroutine derivfdbd(fy,zeta,size,y,ymiss,tdm,kstate,npar,nphi, nx
     *c,nbot,ntop,d1f,d2f,nd,d1a,d1b,d2aa,d2ab,d2bb)
      implicit double precision(a-h,o-z)
      integer ymiss
      integer size
      dimension fy(kstate), zeta(kstate,2), tdm(nxc,kstate)
      dimension d1f(kstate,npar), d2f(kstate,npar,npar)
      dimension d1a(kstate), d1b(kstate)
      dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)
      zero = 0.d0
      do23000 i = 1,kstate 
      do23002 j = 1,npar 
      d1f(i,j) = zero
      do23004 k = 1,npar 
      d2f(i,j,k) = zero
23004 continue
23005 continue
23002 continue
23003 continue
23000 continue
23001 continue
      if(ymiss .gt. 0)then
      return
      endif
      call derivfab(y,fy,kstate,zeta,nbot,ntop,nd,d1a,d1b,d2aa,d2ab,d2bb
     *)
      npro = kstate*(kstate-1)
      nphio2 = nphi/2
      nphio2p1 = nphio2+1
      do23008 i = 1,kstate 
      do23010 j = 1,nxc 
      d1f(i,npro+j) = d1a(i)*tdm(j,i)
      d1f(i,npro+nxc+j) = d1b(i)*tdm(j,i)
23010 continue
23011 continue
23008 continue
23009 continue
      if(nd .gt. 1)then
      do23014 i = 1,kstate 
      do23016 j = 1,nxc
      do23018 k = 1,nxc 
      d2f(i,npro+j,npro+k) = d2aa(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+j,npro+nxc+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+nxc+j,npro+k) = d2ab(i)*tdm(j,i)*tdm(k,i)
      d2f(i,npro+nxc+j,npro+nxc+k) = d2bb(i)*tdm(j,i)*tdm(k,i)
23018 continue
23019 continue
23016 continue
23017 continue
23014 continue
23015 continue
      endif
      return
      end
