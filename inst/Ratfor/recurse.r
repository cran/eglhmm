subroutine recurse(fy,xispd,tpm,epsilon,kstate,n,wrk,xlc,
                   alpha,beta,gamma,xi,xisum,level)
#
implicit double precision(a-h,o-z)

# Realised that there is no need to dimension variables that
# are just being passed on to other subroutines so that their
# indices are not used in the current subroutine!!!

# Update the alpha's.
	call afun(fy,xispd,tpm,epsilon,n,kstate,wrk,xlc,alpha)

# Update the beta's.
	call bfun(fy,xispd,tpm,epsilon,n,kstate,wrk,beta)

# Update the gamma's.
	call gfun(alpha,beta,epsilon,n,kstate,wrk,gamma)
if(level==1) return

# Update the xi's.
	call xfun(alpha,beta,fy,tpm,epsilon,n,kstate,wrk,xi,xisum)

return
end
