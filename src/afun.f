C Output from Public domain Ratfor, version 1.03
      subroutine afun(fy,xispd,tpm,epsilon,n,kstate,wrk,xlc,alpha)
      implicit double precision(a-h,o-z)
      dimension wrk(kstate), xispd(kstate), xlc(n)
      dimension fy(kstate,n), tpm(kstate,kstate), alpha(kstate,n)
      one = 1.d0
      zero = 0.d0
      dummy = -one
      tsum = zero
      do23000 j = 1,kstate 
      wrk(j) = fy(j,1)*xispd(j)
      tsum = tsum + wrk(j)
23000 continue
23001 continue
      if(tsum .lt. epsilon)then
      xlc(1) = dummy
      do23004 j = 1,kstate 
      alpha(j,1) = one/kstate
23004 continue
23005 continue
      else
      xlc(1) = tsum
      do23006 j = 1,kstate 
      alpha(j,1) = wrk(j)/tsum
23006 continue
23007 continue
      endif
      do23008 kt = 2,n 
      tsum = zero
      ktm = kt - 1
      do23010 j = 1,kstate 
      wrk(j) = zero
      do23012 i = 1,kstate 
      wrk(j) = wrk(j) + alpha(i,ktm)*tpm(i,j)
23012 continue
23013 continue
      wrk(j) = fy(j,kt)*wrk(j)
      tsum = tsum + wrk(j)
23010 continue
23011 continue
      if(tsum .lt. epsilon)then
      xlc(kt) = dummy
      do23016 j = 1,kstate 
      alpha(j,kt) = one/kstate
23016 continue
23017 continue
      else
      xlc(kt) = tsum
      do23018 j = 1,kstate 
      alpha(j,kt) = wrk(j)/tsum
23018 continue
23019 continue
      endif
23008 continue
23009 continue
      return
      end