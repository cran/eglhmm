C Output from Public domain Ratfor, version 1.03
      subroutine recurse(fy,xispd,tpm,epsilon,kstate,n,wrk,xlc, alpha,be
     *ta,gamma,xi,xisum,level)
      implicit double precision(a-h,o-z)
      call afun(fy,xispd,tpm,epsilon,n,kstate,wrk,xlc,alpha)
      call bfun(fy,xispd,tpm,epsilon,n,kstate,wrk,beta)
      call gfun(alpha,beta,epsilon,n,kstate,wrk,gamma)
      if(level.eq.1)then
      return
      endif
      call xfun(alpha,beta,fy,tpm,epsilon,n,kstate,wrk,xi,xisum)
      return
      end
