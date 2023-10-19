C Output from Public domain Ratfor, version 1.03
      subroutine pmf(y,x,phi,nyv,nxc,pmfy)
      implicit double precision(a-h,o-z)
      integer y
      dimension x(nxc), phi(nxc,nyv)
      if(y .lt. 1 .or. y .gt. nyv)then
      call intpr1("The value of y is:",-1,y)
      call rexit("This value is out of bounds.\n")
      endif
      nyvm1 = nyv-1
      zed = 1.d0
      top = 1.d0
      do23002 j = 1,nyvm1 
      some = 0.d0
      do23004 i = 1,nxc 
      some = some + x(i)*phi(i,j)
23004 continue
23005 continue
      esome = exp(some)
      if(j.eq.y)then
      top = esome
      endif
      zed = zed + esome
23002 continue
23003 continue
      pmfy = top/zed
      end
