! SUBROUTINE polint.f90 (Numerical Recipes, page 103)

SUBROUTINE polint(xa,ya,n,x,y,dy)
!Given arrays xa and ya, each of length n, and given a value x, this routine returns a
!value y, and an error estimate dy. If P(x) is the polynomial of degree N − 1 such that
!P(xai ) = yai, i = 1,..., n, then the returned value y = P(x).
INTEGER n,NMAX
REAL(8) :: dy,x,y,xa(n),ya(n)
PARAMETER (NMAX=10)
INTEGER i,m,ns
REAL(8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
ns=1
dif=abs(x-xa(1))
do 11 i=1,n
  dift=abs(x-xa(i))
  if (dift.lt.dif) then
    ns=i
    dif=dift
  endif
  c(i)=ya(i)
  d(i)=ya(i)
11    CONTINUE
y=ya(ns)
ns=ns-1
do 13 m=1,n-1
  do 12 i=1,n-m
    ho=xa(i)-x
    hp=xa(i+m)-x
    w=c(i+1)-d(i)
    den=ho-hp
    if(den.eq.0.d0) STOP 'failure in polint'
    den=w/den
    d(i)=hp*den
    c(i)=ho*den
12      CONTINUE
  if (2*ns.lt.n-m)then
    dy=c(ns+1)
  else
    dy=d(ns)
    ns=ns-1
  endif
  y=y+dy
13    CONTINUE
return
END SUBROUTINE
