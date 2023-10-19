subroutine derivfab(y,fy,kstate,ashp,bshp,nbot,ntop,nd,d1a,d1b,d2aa,d2ab,d2bb)
# Derivatives of f w.r.t. alpha and beta.
implicit double precision(a-h,o-z)
dimension fy(kstate)
dimension ashp(kstate), bshp(kstate)
dimension d1a(kstate), d1b(kstate)
dimension d2aa(kstate), d2ab(kstate), d2bb(kstate)

zero = 0.d0
one  = 1.d0
two  = 2.d0
xtop = dble(ntop)
xbot = dble(nbot)
sy   = (y - xbot + one)/(xtop - xbot + two)
t1   = log(sy)
t2   = log(1-sy)
s1   = (two - xbot)/(xtop - xbot + two)
tt1  = log(s1)
tt2  = log(1-s1)
do i = 1,kstate {
    alpha = ashp(i)
    beta  = bshp(i)
    xam   = alpha*tt1 + beta*tt2
    eee   = zero
    d1ea  = zero
    d1eb  = zero
    if(nd==2) {
        d2eaa = zero
        d2eab = zero
        d2ebb = zero
    }
    do k = nbot,ntop {
        sk     = (dble(k)- xbot + one)/(xtop - xbot + two)
        tk1    = log(sk)
        tk2    = log(1-sk)
        xk     = alpha*tk1 + beta*tk2
        if(xk > xam) {
            xam = xk
        }
    }
    do k = nbot,ntop {
        sk     = (dble(k)- xbot + one)/(xtop - xbot + two)
        tk1    = log(sk)
        tk2    = log(1-sk)
        aitchk = one/(sk*(one - sk))
        eee    = eee + aitchk * exp(alpha*tk1 + beta*tk2 - xam)
        d1ea   = d1ea + aitchk * tk1*exp(alpha*tk1 + beta*tk2 - xam)
        d1eb   = d1eb + aitchk * tk2*exp(alpha*tk1 + beta*tk2 - xam)
        if(nd==2) {
            d2eaa  = d2eaa + aitchk * tk1*tk1*exp(alpha*tk1 + beta*tk2 - xam)
            d2eab  = d2eab + aitchk * tk1*tk2*exp(alpha*tk1 + beta*tk2 - xam)
            d2ebb  = d2ebb + aitchk * tk2*tk2*exp(alpha*tk1 + beta*tk2 - xam)
        }
    }
    einv    = one/eee
    pred1a  = einv*d1ea
    pred1b  = einv*d1eb
    d1a(i)  = fy(i)*(t1 - pred1a)
    d1b(i)  = fy(i)*(t2 - pred1b)
    if(nd==2) {
        einv2   = einv**2
        pred2aa = einv*d2eaa - einv2*(d1ea**2)
        pred2ab = einv*d2eab - einv2*d1ea*d1eb
        pred2bb = einv*d2ebb - einv2*(d1eb**2)
        d2aa(i) = fy(i)*((t1 - pred1a)**2 - pred2aa)
        d2ab(i) = fy(i)*((t1 - pred1a)*(t2-pred1b) - pred2ab)
        d2bb(i) = fy(i)*((t2 - pred1b)**2 - pred2bb)
    }
}
return
end
