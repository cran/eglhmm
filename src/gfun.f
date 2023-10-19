C Output from Public domain Ratfor, version 1.03
      subroutine gfun(alpha,beta,epsilon,n,kstate,wrk,gamma)
      implicit double precision(a-h,o-z)
      dimension alpha(kstate,n), beta(kstate,n), gamma(kstate,n)
      dimension wrk(kstate)
      zero = 0.d0
      ook = 1.d0/dble(kstate)
      do23000 kt = 1,n 
      tsum = zero
      do23002 i = 1,kstate 
      wrk(i) = alpha(i,kt)*beta(i,kt)
      tsum = tsum + wrk(i)
23002 continue
23003 continue
      if(tsum.lt.epsilon)then
      do23006 i = 1,kstate 
      gamma(i,kt) = ook
23006 continue
23007 continue
      else
      do23008 i = 1,kstate 
      gamma(i,kt) = wrk(i)/tsum
23008 continue
23009 continue
      endif
23000 continue
23001 continue
      return
      end
