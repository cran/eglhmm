C Output from Public domain Ratfor, version 1.03
      subroutine xfun(alpha,beta,fy,tpm,epsilon,n,kstate,wrk,xi,xisum)
      implicit double precision(a-h,o-z)
      dimension alpha(kstate,n), beta(kstate,n), fy(kstate,n)
      dimension tpm(kstate,kstate), wrk(kstate,kstate), xi(kstate,kstate
     *,n)
      dimension xisum(kstate,kstate)
      one = 1.d0
      zero = 0.d0
      dns2 = dble(kstate*kstate)
      do23000 ktp = 2,n 
      kt = ktp - 1
      tsum = zero
      do23002 i = 1,kstate 
      do23004 j = 1, kstate 
      wrk(i,j) = alpha(i,kt)*fy(j,ktp)*beta(j,ktp)*tpm(i,j)
      tsum = tsum + wrk(i,j)
23004 continue
23005 continue
23002 continue
23003 continue
      if(tsum.lt.epsilon)then
      do23008 i = 1,kstate 
      do23010 j = 1,kstate 
      xi(i,j,kt) = one/dns2
23010 continue
23011 continue
23008 continue
23009 continue
      else
      do23012 i = 1,kstate 
      do23014 j = 1,kstate 
      xi(i,j,kt) = wrk(i,j)/tsum
23014 continue
23015 continue
23012 continue
23013 continue
      endif
23000 continue
23001 continue
      nm1 = n - 1
      do23016 i = 1,kstate 
      do23018 j = 1,kstate 
      xisum(i,j) = zero
      do23020 k = 1,nm1 
      xisum(i,j) = xisum(i,j) + xi(i,j,k)
23020 continue
23021 continue
23018 continue
23019 continue
23016 continue
23017 continue
      return
      end
