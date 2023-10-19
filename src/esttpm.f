C Output from Public domain Ratfor, version 1.03
      subroutine esttpm(ns,n,k,tpm,mixture,wrk)
      implicit double precision(a-h,o-z)
      dimension ns(n), tpm(k,k), wrk(k)
      zero = 0.d0
      one = 1.d0
      ook = one/dble(k)
      do23000 i = 1,k 
      do23002 j = 1,k 
      tpm(i,j) = zero
23002 continue
23003 continue
23000 continue
23001 continue
      do23004 nt = 2,n 
      nb = nt-1
      do23006 i = 1,k 
      do23008 j = 1,k 
      if(ns(nb).eq.i.and.ns(nt).eq.j)then
      tpm(i,j) = tpm(i,j)+one
      endif
23008 continue
23009 continue
23006 continue
23007 continue
23004 continue
23005 continue
      if(mixture .gt. 0)then
      den = zero
      do23014 j = 1,k 
      wrk(j) = zero
      do23016 i = 1,k 
      den = den + tpm(i,j)
      wrk(j) = wrk(j) + tpm(i,j)
23016 continue
23017 continue
23014 continue
23015 continue
      do23018 i = 1,k 
      do23020 j = 1,k 
      tpm(i,j) = wrk(j)/den
23020 continue
23021 continue
23018 continue
23019 continue
      else
      do23022 i = 1,k 
      den = zero
      do23024 j = 1,k 
      den = den + tpm(i,j)
23024 continue
23025 continue
      if(den.ge.one)then
      do23028 j = 1,k 
      tpm(i,j) = tpm(i,j)/den
23028 continue
23029 continue
      else
      do23030 j = 1,k 
      tpm(i,j) = ook
23030 continue
23031 continue
      endif
23022 continue
23023 continue
      endif
      return
      end
