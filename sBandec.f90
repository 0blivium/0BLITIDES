! SUBROUTINE bandec.f90 (Numerical Recipes, page 45)

SUBROUTINE bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
! DESCRIPTION:
! Given an n × n band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal
! rows, compactly stored in the array a(1:n,1:m1+m2+1), this routine constructs 
! an LU decomposition of a rowwise permutation of A. 
! The upper triangular matrix replaces a, while the lower triangular matrix is returned
! in al(1:n,1:m1). indx(1:n) is an output vector which records the row permutation
! effected by the partial pivoting; d is output as ±1 depending on whether the number of
! row interchanges was even or odd, respectively.
INTEGER m1,m2,mp,mpl,n,np,indx(n)
DOUBLE PRECISION d,a(np,mp),al(np,mpl),TINY
PARAMETER (TINY=1.d-20)
INTEGER i,j,k,l,mm
DOUBLE PRECISION dum
mm=m1+m2+1
IF(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) STOP 'bad args in bandec'
l=m1
DO i=1,m1
  DO j=m1+2-i,mm
    a(i,j-l)=a(i,j)
  ENDDO
  l=l-1
  DO j=mm-l,mm
    a(i,j)=0.d0
  ENDDO
ENDDO
d=1.d0
l=m1
DO k=1,n
  dum=a(k,1)
  i=k
  IF(l.lt.n)l=l+1
  DO j=k+1,l
    IF(abs(a(j,1)).gt.abs(dum)) THEN
      dum=a(j,1)
      i=j
    ENDIF
  ENDDO
  indx(k)=i
  IF(dum.eq.0.d0) a(k,1)=TINY
  IF(i.ne.k) THEN
    d=-d
    DO j=1,mm
      dum=a(k,j)
      a(k,j)=a(i,j)
      a(i,j)=dum
    ENDDO
  ENDIF
  DO i=k+1,l
    dum=a(i,1)/a(k,1)
    al(k,i-k)=dum
    DO j=2,mm
      a(i,j-1)=a(i,j)-dum*a(k,j)
    ENDDO
    a(i,mm)=0.d0
  ENDDO
ENDDO
RETURN
END SUBROUTINE