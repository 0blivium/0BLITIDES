! SUBROUTINE banbks.f90 (Numerical Recipes, page 46)

SUBROUTINE banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
! DESCRIPTION:
! Given the arrays a, al, and indx as returned from BANDEC, and given a right-hand side
! vector b(1:n), solves the band diagonal linear equations A·x = b. The solution vector x
! overwrites b(1:n). The other input arrays are not modified, and can be left in place for
! successive calls with different right-hand sides.
INTEGER m1,m2,mp,mpl,n,np,indx(n)
DOUBLE PRECISION a(np,mp),al(np,mpl),b(n)
INTEGER i,k,l,mm
DOUBLE PRECISION dum
mm=m1+m2+1
IF(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) STOP 'bad args in banbks'
l=m1
DO k=1,n
  i=indx(k)
  IF(i.ne.k) THEN
    dum=b(k)
    b(k)=b(i)
    b(i)=dum
  ENDIF
  IF(l.lt.n)l=l+1
  DO i=k+1,l
    b(i)=b(i)-al(k,i-k)*b(k)
  ENDDO
ENDDO
l=1
DO i=n,1,-1
  dum=b(i)
  DO k=2,l
    dum=dum-a(i,k)*b(k+i-1)
  ENDDO
  b(i)=dum/a(i,1)
  IF(l.lt.mm) l=l+1
ENDDO
RETURN
END SUBROUTINE