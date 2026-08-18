! MODULE: Coefficients of spectrally decomposed PDEs.

MODULE mCoeff
USE mPar
USE mFunc
IMPLICIT NONE
    REAL(8) :: cmat(10,10)=0d0
    CONTAINS
        SUBROUTINE init_coeff(j,g,rho)
            REAL(8),INTENT(IN) :: j
            REAL(8),INTENT(IN) :: g(N+1),rho(N+1)
                ! continuity equation:
                cmat(1,1)=-sqrt(j/(2.*j+1.))*(j-1.)                                             ! = a2(j)
                cmat(1,2)=-sqrt((j+1.)/(2.*j+1.))*(j+2)                                         ! = b2(j)
                cmat(1,3)=sqrt(j/(2.*j+1.))                                                     ! = a1(j)  
                cmat(1,4)=-sqrt((j+1.)/(2.*j+1.))                                               ! = b1(j)
                cmat(1,5)=(1./sqrt(3.))                                                         ! = п1(j)
                cmat(1,6)=sqrt(1./(j*(j+1.)))                                                   ! = п2(j)
                ! rheology equation:
                cmat(2,1)=sqrt((j-1.)/(2.*j-1.))*j                                              ! = p2(j)
                cmat(2,2)=sqrt((j-1.)/(2.*j-1.))                                                ! = p1(j)
                cmat(2,3)=sqrt(((j+1.)*(2.*j+3.))/(6.*(2.*j-1.)*(2.*j+1.)))*(j-1.)              ! = q2(j)
                cmat(2,4)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3.)))*(j+2.)                   ! = r2(j)
                cmat(2,5)=-sqrt(((j+1.)*(2.*j+3.))/(6.*(2.*j-1.)*(2.*j+1.)))                    ! = q1(j)
                cmat(2,6)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3.)))                          ! = r1(j)
                cmat(2,7)=sqrt((j+2.)/(2.*j+3.))*(j+1.)                                         ! = s2(j)
                cmat(2,8)=-sqrt((j+2.)/(2.*j+3.))                                               ! = s1(j)
                ! equation of motion: Y_j-1:
                cmat(3,1)=-sqrt(j/(3.*(2.*j+1.)))*(j+1.)                                        ! = c2(j)
                cmat(3,2)=-sqrt((j-1.)/(2.*j-1.))*(j-2.)                                        ! = d2(j)
                cmat(3,3)=-sqrt(((j+1.)*(2.*j+3.))/(6.*(2.*j-1.)*(2.*j+1.)))*(j+1.)             ! = e2(j)
                cmat(3,4)=-sqrt(j/(3.*(2.*j+1.)))                                               ! = c1(j)
                cmat(3,5)=sqrt((j-1.)/(2.*j-1.))                                                ! = d1(j)
                cmat(3,6)=-sqrt(((j+1.)*(2.*j+3.))/(6.*(2.*j-1)*(2.*j+1.)))                     ! = e1(j)
                ! equation of motion: Y_j+1:
                cmat(4,1)=-sqrt((j+1.)/(3.*(2.*j+1.)))*j                                        ! = f2(j)
                cmat(4,2)=-sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3.)))*j                       ! = g2(j)
                cmat(4,3)=-sqrt((j+2.)/(2.*j+3.))*(j+3.)                                        ! = h2(j)
                cmat(4,4)=sqrt((j+1.)/(3.*(2.*j+1.)))                                           ! = f1(j)
                cmat(4,5)=sqrt(j*(2.*j-1.)/(6.*(2.*j+1.)*(2.*j+3.)))                            ! = g1(j)
                cmat(4,6)=-sqrt((j+2)/(2*j+3))                                                  ! = h1(j)
                ! top boundary condition: Y_j-1:
                cmat(5,1)=rho(N+1)*g(N+1)*j/(2.*j+1.)                                           ! = z1(j)
                cmat(5,2)=-rho(N+1)*g(N+1)*sqrt(j*(j+1.))/(2.*j+1.)                             ! = z2(j)
                cmat(5,3)=-sqrt(j/(3.*(2.*j+1.)))                                               ! = z3(j)
                cmat(5,4)=sqrt((j-1.)/(2.*j-1.))                                                ! = z4(j)
                cmat(5,5)=-sqrt(((j+1.)*(2.*j+3))/(6.*(2.*j+1.)*(2.*j-1.)))                     ! = z5(j)
                ! top boundary condition: Y_j+1:
                cmat(6,1)=-rho(N+1)*g(N+1)*sqrt(j*(j+1.))/(2.*j+1.)                             ! = w1(j)
                cmat(6,2)=rho(N+1)*g(N+1)*(j+1.)/(2.*j+1.)                                      ! = w2(j)
                cmat(6,3)=sqrt((j+1.)/(3.*(2.*j+1.)))                                           ! = w3(j)
                cmat(6,4)=sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3.)))                          ! = w4(j)
                cmat(6,5)=-sqrt((j+2)/(2*j+3))                                                  ! = w5(j)
                ! bottom boundary condition: Y_j-1:
                cmat(7,1)=(rho(1)-rho(2))*g(1)*j/(2.*j+1.)                                      ! = x1(j)
                cmat(7,2)=-(rho(1)-rho(2))*g(1)*sqrt(j*(j+1.))/(2.*j+1.)                        ! = x2(j)
                cmat(7,3)=sqrt(j/(3.*(2.*j+1.)))                                                ! = x3(j)
                cmat(7,4)=-sqrt((j-1.)/(2.*j-1.))                                               ! = x4(j)
                cmat(7,5)=sqrt(((j+1.)*(2.*j+3.))/(6.*(2.*j+1.)*(2.*j-1.)))                     ! = x5(j)
                ! bottom boundary condition: Y_j+1:
                cmat(8,1)=-(rho(1)-rho(2))*g(1)*sqrt(j*(j+1.))/(2.*j+1.)                        ! = y1(j)
                cmat(8,2)=(rho(1)-rho(2))*g(1)*(j+1.)/(2.*j+1.)                                 ! = y2(j)
                cmat(8,3)=-sqrt((j+1.)/(3.*(2.*j+1.)))                                          ! = y3(j)
                cmat(8,4)=-sqrt((j*(2.*j-1.))/(6.*(2.*j+1.)*(2.*j+3.)))                         ! = y4(j)
                cmat(8,5)=sqrt((j+2.)/(2.*j+3.))                                                ! = y5(j)
                ! gradient of density and compressible equation of motion: Y_j-1:
                cmat(9,1)=-sqrt(j*(j+1.))/(2.*j+1.)                                             ! = v1(j) 
                cmat(9,2)=j/(2.*j+1.)                                                           ! = v2(j)
                ! gradient of density and compressible equation of motion: Y_j+1:
                cmat(10,1)=(j+1.)/(2.*j+1.)                                                     ! = v3(j)
                cmat(10,2)=-sqrt(j*(j+1.))/(2.*j+1.)                                            ! = v4(j)
        END SUBROUTINE init_coeff
END MODULE mCoeff
