! MODULE: Auxiliery subroutines implementing spectrally decomposed PDEs into a matrix (row by row).

MODULE mEqs
USE mPar
IMPLICIT NONE

    ! ORIENTATION: matrix goes from the bottom (core) to the top (surface)
    ! SHAPE: compactly stored in a (nn) x (m1+m2+1) matrix
    ! BOUNDARY CONDITION: dynamic topography by default (bottom bc overridden by no-slip in main.f90 if necessary)
    ! INTERFACES: continuity equation, rheology equations, boundary conditions
    ! CENTERS: equations of motion

    CONTAINS
        ! CONTINUITY EQUATION
        SUBROUTINE con_eq(A,cmat,i,k,dr,r,rhoF,gF,lmd,mu)
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(N),INTENT(IN)            :: r,gF,rhoF,lmd,mu
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od+0)=cmat(1,1)/(2d0*r(k))-cmat(1,3)/dr  ! u--
                A(i,od+1)=cmat(1,2)/(2d0*r(k))-cmat(1,4)/dr  ! u-+
                A(i,od+6)=cmat(1,1)/(2d0*r(k))+cmat(1,3)/dr  ! u+-
                A(i,od+7)=cmat(1,2)/(2d0*r(k))+cmat(1,4)/dr  ! u++
                IF(COMPRESS) THEN !% compressibility term
                    A(i,od+2)=A(i,od+2)+cmat(1,5)                                              ! t0
                    A(i,od+0)=A(i,od+0)*(lmd(k)+(2d0/3d0)*mu(k))-(cmat(1,3)/2d0)*rhoF(k)*gF(k) ! u--
                    A(i,od+1)=A(i,od+1)*(lmd(k)+(2d0/3d0)*mu(k))-(cmat(1,4)/2d0)*rhoF(k)*gF(k) ! u-+
                    A(i,od+6)=A(i,od+6)*(lmd(k)+(2d0/3d0)*mu(k))-(cmat(1,3)/2d0)*rhoF(k)*gF(k) ! u+-
                    A(i,od+7)=A(i,od+7)*(lmd(k)+(2d0/3d0)*mu(k))-(cmat(1,4)/2d0)*rhoF(k)*gF(k) ! u++
                ENDIF
                i=i+1
        END SUBROUTINE con_eq   

        ! RHEOLOGY EQUATION Y_j-2,2
        SUBROUTINE rheo_eq_jm2(A,cmat,i,k,dr,r,mu)    
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(n),INTENT(IN)            :: r,mu
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od-3)=cmat(2,1)/(2d0*r(k))-cmat(2,2)/dr  ! u--
                A(i,od+0)=(-1d0)/(2d0*mu(k))                 ! t-2
                A(i,od+3)=cmat(2,1)/(2d0*r(k))+cmat(2,2)/dr  ! u-+
                i=i+1
        END SUBROUTINE rheo_eq_jm2

        ! RHEOLOGY EQUATION Y_j+2,2
        SUBROUTINE rheo_eq_jp2(A,cmat,i,k,dr,r,mu)   
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(out)                        :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(n),INTENT(IN)            :: r,mu
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od-4)=cmat(2,7)/(2d0*r(k))-cmat(2,8)/dr  ! u+-
                A(i,od+0)=(-1d0)/(2d0*mu(k))                 ! t+2
                A(i,od+2)=cmat(2,7)/(2d0*r(k))+cmat(2,8)/dr  ! u++
                i=i+1
        END SUBROUTINE rheo_eq_jp2

        ! RHEOLOGY EQUATION Y_j,2
        SUBROUTINE rheo_eq_j(A,cmat,i,k,dr,r,mu)
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(n),INTENT(IN)            :: r,mu
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od-4)=cmat(2,3)/(2d0*r(k))-cmat(2,5)/dr  ! u--
                A(i,od-3)=cmat(2,4)/(2d0*r(k))-cmat(2,6)/dr  ! u-+
                A(i,od+0)=(-1d0)/(2d0*mu(k))                 ! t2
                A(i,od+2)=cmat(2,3)/(2d0*r(k))+cmat(2,5)/dr  ! u+-
                A(i,od+3)=cmat(2,4)/(2d0*r(k))+cmat(2,6)/dr  ! u++
                i=i+1
        END SUBROUTINE rheo_eq_j

        ! EQUATION OF MOTION Y_j-1
        SUBROUTINE mot_eq_jm1(A,cmat,i,k,dr,r,rhoC,rhoF,gC,mu,lmd)
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(N),INTENT(IN)            :: r,rhoF,lmd,mu
            REAL(8),DIMENSION(N+1),INTENT(IN)          :: rhoC,gC
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od-6)=cmat(3,1)/(2d0*(r(k-1)+dr/2d0))-cmat(3,4)/dr  ! t0
                A(i,od-5)=cmat(3,2)/(2d0*(r(k-1)+dr/2d0))-cmat(3,5)/dr  ! t-2
                A(i,od-4)=cmat(3,3)/(2d0*(r(k-1)+dr/2d0))-cmat(3,6)/dr  ! t2
                A(i,od+0)=cmat(3,1)/(2d0*(r(k-1)+dr/2d0))+cmat(3,4)/dr  ! t0
                A(i,od+1)=cmat(3,2)/(2d0*(r(k-1)+dr/2d0))+cmat(3,5)/dr  ! t-2
                A(i,od+2)=cmat(3,3)/(2d0*(r(k-1)+dr/2d0))+cmat(3,6)/dr  ! t2
                IF(COMPRESS) THEN  !% compressibility term
                    A(i,od-6)=A(i,od-6)-(cmat(1,5)/2d0)*cmat(1,3)*rhoC(k)*gC(k)/(lmd(k)+(2d0/3d0)*mu(k)) ! t0
                    A(i,od+0)=A(i,od+0)-(cmat(1,5)/2d0)*cmat(1,3)*rhoC(k)*gC(k)/(lmd(k)+(2d0/3d0)*mu(k)) ! t0
                    A(i,od-2)=A(i,od-2)+cmat(9,2)*((rhoC(k)*gC(k))**2)/(lmd(k)+(2d0/3d0)*mu(k))          ! u-
                    A(i,od-1)=A(i,od-1)+cmat(9,1)*((rhoC(k)*gC(k))**2)/(lmd(k)+(2d0/3d0)*mu(k))          ! u+
                ENDIF
                IF(abs(rhoF(k+1)-rhoF(k)).gt.JTINY.and.k.lt.N) THEN !% density jumps
                    A(i+6,od-2)=A(i+6,od-2)+cmat(9,2)*((rhoF(k+1)-rhoF(k))/dr)*gC(k+1) ! u-
                    A(i+6,od-1)=A(i+6,od-1)+cmat(9,1)*((rhoF(k+1)-rhoF(k))/dr)*gC(k+1) ! u+
                ENDIF
                i=i+1
        END SUBROUTINE mot_eq_jm1

        ! EQUATION OF MOTION Y_j+1
        SUBROUTINE mot_eq_jp1(A,cmat,i,k,dr,r,rhoC,rhoF,gC,mu,lmd)
            INTEGER,INTENT(IN)                         :: k
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),INTENT(IN)                         :: dr
            REAL(8),DIMENSION(N),INTENT(IN)            :: r,rhoF,lmd,mu
            REAL(8),DIMENSION(N+1),INTENT(IN)          :: rhoC,gC
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
                A(i,od-5)=cmat(4,1)/(2d0*(r(k-1)+dr/2d0))-cmat(4,4)/dr  ! t0
                A(i,od-3)=cmat(4,2)/(2d0*(r(k-1)+dr/2d0))-cmat(4,5)/dr  ! t2
                A(i,od-2)=cmat(4,3)/(2d0*(r(k-1)+dr/2d0))-cmat(4,6)/dr  ! t+2
                A(i,od+1)=cmat(4,1)/(2d0*(r(k-1)+dr/2d0))+cmat(4,4)/dr  ! t0
                A(i,od+3)=cmat(4,2)/(2d0*(r(k-1)+dr/2d0))+cmat(4,5)/dr  ! t2
                A(i,od+4)=cmat(4,3)/(2d0*(r(k-1)+dr/2d0))+cmat(4,6)/dr  ! t+2
                IF(COMPRESS) THEN   !% compressibility term
                    A(i,od-5)=A(i,od-5)-cmat(1,5)*(cmat(1,4)/2d0)*rhoC(k)*gC(k)/(lmd(k)+(2d0/3d0)*mu(k)) ! t0
                    A(i,od+1)=A(i,od+1)-cmat(1,5)*(cmat(1,4)/2d0)*rhoC(k)*gC(k)/(lmd(k)+(2d0/3d0)*mu(k)) ! t0
                    A(i,od-1)=A(i,od-1)+cmat(9,1)*((rhoC(k)*gC(k))**2)/(lmd(k)+(2d0/3d0)*mu(k))          ! u-
                    A(i,od+0)=A(i,od+0)+cmat(10,1)*((rhoC(k)*gC(k))**2)/(lmd(k)+(2d0/3d0)*mu(k))         ! u+
                ENDIF
                IF(abs(rhoF(k+1)-rhoF(k)).gt.JTINY.and.k.lt.N) THEN !% density jumps
                    A(i+6,od-1)=A(i+6,od-1)+cmat(10,2)*((rhoF(k+1)-rhoF(k))/dr)*gC(k+1) ! u-   
                    A(i+6,od+0)=A(i+6,od+0)+cmat(10,1)*((rhoF(k+1)-rhoF(k))/dr)*gC(k+1) ! u+
                ENDIF
                i=i+1
        END SUBROUTINE mot_eq_jp1

        ! TOP BOUNDARY CONDITION Y_j-1
        SUBROUTINE top_boundary_jm1(A,cmat,i)
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
                A(i,od-6)=cmat(5,1)/2d0     ! u--
                A(i,od-5)=cmat(5,2)/2d0     ! u-+
                A(i,od-4)=cmat(5,3)         ! t0
                A(i,od-3)=cmat(5,4)         ! t-2
                A(i,od-2)=cmat(5,5)         ! t2
                A(i,od+0)=cmat(5,1)/2d0     ! u+-
                A(i,od+1)=cmat(5,2)/2d0     ! u++
                i=i+1
        END SUBROUTINE top_boundary_jm1

        ! TOP BOUNDARY CONDITION Y_j+1
        SUBROUTINE top_boundary_jp1(A,cmat,i)
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
                A(i,od-7)=cmat(6,1)/2d0     ! u--
                A(i,od-6)=cmat(6,2)/2d0     ! u-+
                A(i,od-5)=cmat(6,3)         ! t0
                A(i,od-3)=cmat(6,4)         ! t2
                A(i,od-2)=cmat(6,5)         ! t+2
                A(i,od-1)=cmat(6,1)/2d0     ! u+-
                A(i,od+0)=cmat(6,2)/2d0     ! u++
                i=i+1
        END SUBROUTINE top_boundary_jp1

        ! BOTTOM BOUNDARY CONDITION Y_j-1
        SUBROUTINE bot_boundary_jm1(A,cmat,i)         
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
                A(i,od-1)=cmat(7,1)/2d0     ! u--
                A(i,od+0)=cmat(7,2)/2d0     ! u-+
                A(i,od+1)=cmat(7,3)         ! t0
                A(i,od+2)=cmat(7,4)         ! t-2
                A(i,od+3)=cmat(7,5)         ! t2
                A(i,od+5)=cmat(7,1)/2d0     ! u+-
                A(i,od+6)=cmat(7,2)/2d0     ! u++
                i=i+1
        END SUBROUTINE bot_boundary_jm1

        ! BOTTOM BOUNDARY CONDITION Y_j+1
        SUBROUTINE bot_boundary_jp1(A,cmat,i)
            INTEGER,INTENT(INOUT)                      :: i
            REAL(8),DIMENSION(nn,m),INTENT(INOUT)      :: A
            REAL(8),DIMENSION(10,10),INTENT(IN)        :: cmat
                A(i,od-2)=cmat(8,1)/2d0     ! u--
                A(i,od-1)=cmat(8,2)/2d0     ! u-+
                A(i,od+0)=cmat(8,3)         ! t0
                A(i,od+2)=cmat(8,4)         ! t-2
                A(i,od+3)=cmat(8,5)         ! t2
                A(i,od+4)=cmat(8,1)/2d0     ! u+-
                A(i,od+5)=cmat(8,2)/2d0     ! u++
                i=i+1
        END SUBROUTINE bot_boundary_jp1
END MODULE mEqs