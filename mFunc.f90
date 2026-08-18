! MODULE: Procedures utilized in main.f90.

MODULE mFunc
USE mPar 
USE mEqs
    CONTAINS
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%% FUNCTIONS %%

        ! TIDES (OBLIQUITY & ECCENTRIC)
        FUNCTION obltidepot(t,r) RESULT(Vj2)                 !% obliquity tidal potential
            REAL(8)     :: t,r
            COMPLEX(8)  :: Vj2
                Vj2=(3d0/2d0)*sqrt((2d0*PI)/15d0)*(OMEGA**2)*(r**2)*THETA0*PI/180d0*(cos(OMEGA*t)+complex(0,1)*sin(OMEGA*t))
        END FUNCTION obltidepot

        FUNCTION obltidef(t,r) RESULT(fjm1)                  !% obliquity tidal force, Y_j-1
            REAL(8)     :: t,r  
            COMPLEX(8)  :: fjm1
                fjm1=sqrt(3d0*PI)*THETA0*PI/180d0*(OMEGA**2)*r*(cos(OMEGA*t)+complex(0,1)*sin(OMEGA*t))
        END FUNCTION obltidef

        ! SELFGRAVITY
        FUNCTION sgext(r,ur,j,rbar,drho) RESULT (sgtopj)     !% self-gravitational potential (exterior)
            REAL(8) :: j,r,ur,rbar,drho,sgtopj
                sgtopj=((-4d0*ur*Gconst*PI*r*drho)/(2d0*j+1d0))*(rbar/r)**(j+2d0)
        END FUNCTION sgext

        FUNCTION sgint(r,ur,j,rbar,drho) RESULT (sgbotj)     !% self-gravitational potential (interior) 
            REAL(8) :: j,r,ur,rbar,drho,sgbotj
                sgbotj=((-4d0*ur*Gconst*PI*r*drho))/(2d0*j+1d0)*(r/rbar)**(j-1d0)
        END FUNCTION sgint

        FUNCTION sgfext(r,ur,j,rbar,drho) RESULT (fjp1)      ! self-gravitational fj-1 force
            REAL(8) :: j,r,ur,rbar,drho,fjp1
                fjp1=4d0*ur*PI*Gconst*drho*sqrt((j+1d0)/(2d0*j+1d0))*(rbar/r)**(j+2d0)
        END FUNCTION sgfext

        FUNCTION sgfint(r,ur,j,rbar,drho) RESULT (fjm1)      ! self-gravitational fj+1 force
            REAL(8) :: j,r,ur,rbar,drho,fjm1
                fjm1=4d0*ur*PI*Gconst*drho*sqrt(j/(2d0*j+1d0))*(r/rbar)**(j-1d0)
        END FUNCTION sgfint
    
        ! AUXILIARY ROUNDING FUNCTION:
        FUNCTION round(num) RESULT(idx)
            INTEGER            :: idx
            REAL(8),INTENT(IN) :: num
                IF(abs((nint(num))+1d0-num).lt.1) THEN
                    idx=nint(num)+1
                ELSE; idx=nint(num); ENDIF;
        END FUNCTION round

        ! DISSIPATION FOR RADIALLY DEPENDENT VISCOSITY:
        FUNCTION calc_diss(bRe,bIm,dr,etaFaces,rFaces) RESULT(diss)
            INTEGER                               :: i
            REAL(8)                               :: diss,tmp
            REAL(8),INTENT(IN)                    :: dr
            REAL(8),DIMENSION(NN+2),INTENT(IN)    :: bRe,bIm 
            REAL(8),DIMENSION(N),INTENT(IN)       :: etaFaces,rFaces
                diss=0d0
                DO i=1,N !% iterate through the layers
                    tmp=(bRe(6*i+2)**2+bIm(6*i+2)**2)+(bRe(6*i+3)**2+bIm(6*i+3)**2)+(bRe(6*i+4)**2+bIm(6*i+4)**2)
                    tmp=tmp/etaFaces(i)
                    IF(i.eq.0.or.i.eq.N-1) THEN
                        tmp=tmp*(rFaces(i)**2)*dr/2d0
                    ELSE
                        tmp=tmp*(rFaces(i)**2)*dr
                    ENDIF
                    diss=diss+tmp
                ENDDO
        END FUNCTION calc_diss

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%% SUBROUTINES %%

        ! INITIALIZATION OF THE BANDED MATRIX:
        SUBROUTINE init_matrix(A,c,li,mi,dr,r,rhoC,rhoF,gC,gF,mu,lmd)
            !% variable order in the matrix: u-,u+,t0,t-2,t2,t+2
            INTEGER                                :: l
            INTEGER,INTENT(INOUT)                  :: li,mi
            REAL(8),INTENT(IN)                     :: dr
            REAL(8),DIMENSION(N),INTENT(IN)        :: r,gF,mu,lmd
            REAL(8),DIMENSION(N),INTENT(IN)        :: rhoF
            REAL(8),DIMENSION(N+1),INTENT(IN)      :: rhoC,gC    
            REAL(8),DIMENSION(10,10),INTENT(IN)    :: c
            REAL(8),DIMENSION(NN,m),INTENT(INOUT)  :: A 
                !% 1st block & bottom boundary conditions:
                CALL con_eq(A,c,mi,li,dr,r,rhoF,gF,lmd,mu)
                CALL bot_boundary_jm1(A,c,mi)
                CALL bot_boundary_jp1(A,c,mi)
                CALL rheo_eq_jm2(A,c,mi,li,dr,r,mu)
                CALL rheo_eq_j(A,c,mi,li,dr,r,mu)
                CALL rheo_eq_jp2(A,c,mi,li,dr,r,mu)
                !% inner blocks:
                LAYER_LOOP: DO l=1,N-1 !% iteration over N-1 layers!
                li=li+1
                CALL con_eq(A,c,mi,li,dr,r,rhoF,gF,lmd,mu)
                CALL mot_eq_jp1(A,c,mi,li,dr,r,rhoC,rhoF,gC,mu,lmd)
                CALL mot_eq_jm1(A,c,mi,li,dr,r,rhoC,rhoF,gC,mu,lmd)
                CALL rheo_eq_jm2(A,c,mi,li,dr,r,mu)
                CALL rheo_eq_j(A,c,mi,li,dr,r,mu)
                CALL rheo_eq_jp2(A,c,mi,li,dr,r,mu)
                END DO LAYER_LOOP
                ! top boundary conditions:
                CALL top_boundary_jm1(A,c,mi)
                CALL top_boundary_jp1(A,c,mi)
        END SUBROUTINE init_matrix
    
        ! PRINT MATRIX INTO A TEXT FILE:
        SUBROUTINE print_mat(A,nrow,ncol)
            INTEGER,INTENT(IN)                       :: nrow,ncol
            REAL(8),DIMENSION(nrow,ncol), INTENT(IN) :: A
            INTEGER                                  :: j1,j2
                OPEN(11,file='run/banded_matrix.txt')
                DO j1=1,nrow
                    WRITE(11,*) (A(j1,j2), j2=1,ncol)
                END DO
                CLOSE(11)        
        END SUBROUTINE print_mat

        ! AUXILIARY PROCEDURE FOR ARRAY REVERSING:
        SUBROUTINE arr_reverse(array) 
            REAL(8),DIMENSION(:),INTENT(INOUT) :: array
            INTEGER                            :: indx
                array=[(array(size(array)+1-indx),indx=1,size(array))]
        END SUBROUTINE
    
        ! INITIALIZE THE SPATIAL DISCRETIZATION:
        SUBROUTINE init_rad(dr,rCenters,rFaces)
            INTEGER                               :: l
            REAL(8),INTENT(INOUT)                 :: dr      
            REAL(8),DIMENSION(N),INTENT(INOUT)    :: rFaces
            REAL(8),DIMENSION(N+1),INTENT(INOUT)  :: rCenters
                dr=(RTOP-RBOT)/(N-1)
                DO l=1,N
                    rFaces(l)=RBOT+dr*(l-1)
                    IF (l.eq.1) THEN; rCenters(l)=rFaces(l); ELSEIF (l.eq.2) THEN;&
                        rCenters(l)=rFaces(l)-dr/2d0; ELSEIF (l.eq.N) THEN;&
                        rCenters(l)=rFaces(l)-dr/2d0; ELSE; rCenters(l)=rCenters(l-1)+dr;&
                    ENDIF
                END DO 
                rCenters(N+1)=RTOP
        END SUBROUTINE init_rad
    
        ! INITIALIZE THE HOMOGENEOUS MODEL:
        SUBROUTINE init_homog_model(rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces)
            REAL(8),DIMENSION(N),INTENT(INOUT)   :: rhoFaces,muFaces,lmdFaces,etaFaces
            REAL(8),DIMENSION(N+1),INTENT(INOUT) :: rhoCenters
                muFaces=muHmg               
                lmdFaces=lmdHmg
                IF(.not.ETAEXP) THEN; etaFaces=etaHmg; ENDIF;
                rhoCenters=rhoHmg           
                rhoCenters(1)=rhoCore
                rhoFaces=rhoHmg
        END SUBROUTINE init_homog_model

        ! INITIALIZE THE LAYERED MODEL:
        SUBROUTINE init_layered_model(rCenters,rFaces,rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces,jindx,dr)
            INTEGER                                 :: i
            INTEGER,INTENT(INOUT)                   :: jindx(noj)
            REAL(8),INTENT(IN)                      :: dr
            REAL(8),DIMENSION(N),INTENT(INOUT)      :: rFaces,rhoFaces,muFaces,etaFaces,lmdFaces
            REAL(8),DIMENSION(N+1),INTENT(INOUT)    :: rCenters,rhoCenters
                DO i=1,noj
                    WHERE(rCenters.gt.rJumps(i).and.rCenters.lt.rJumps(i+1)) 
                        rhoCenters=rhoJumps(i+1)
                    END WHERE
                    WHERE(rFaces.gt.rJumps(i).and.rFaces.lt.rJumps(i+1)) 
                        muFaces=muJumps(i+1)
                        etaFaces=etaJumps(i+1)
                        rhoFaces=rhoJumps(i+1)
                        lmdFaces=lmdJumps(i+1)
                    END WHERE
                    jindx(i)=round((rJumps(i+1)-rJumps(1))/dr)
                END DO
                !% the last chunk:
                rhoCenters(jindx(noj)+2:)=rhoJumps(noj+2)
                rhoFaces(jindx(noj)+1:)=rhoJumps(noj+2)
                muFaces(jindx(noj)+1:)=muJumps(noj+2)
                etaFaces(jindx(noj)+1:)=etaJumps(noj+2)
                lmdFaces(jindx(noj)+1:)=lmdJumps(noj+2)
                !% core:
                rhoFaces(1)=rhoJumps(2)
                rhoCenters(1)=rhoJumps(1) 
                etaFaces(1)=etaJumps(2)
                muFaces(1)=muJumps(2)
                lmdFaces(1)=lmdJumps(2)
        END SUBROUTINE

        ! INITIALIZE REFERENCE GRAVITY g0:
        SUBROUTINE init_g0m0(rCenters,rhoCenters,gCenters,rFaces,rhoFaces,gFaces,mFaces,mCenters)
            INTEGER                                 :: ig,jcntr=1
            REAL(8)                                 :: mShell,rhoShell,tmpdr
            REAL(8),DIMENSION(N),INTENT(INOUT)      :: rFaces,rhoFaces
            REAL(8),DIMENSION(N),INTENT(INOUT)      :: gFaces,mFaces 
            REAL(8),DIMENSION(N+1),INTENT(INOUT)    :: gCenters,mCenters
            REAL(8),DIMENSION(N+1),INTENT(INOUT)    :: rCenters,rhoCenters
                !% Analytical formula for layered and homogeneous models
                IF(MODEL.eq.1.or.MODEL.eq.2) THEN
                    ! gravitational acceleration evaluated at centers:
                    mCenters(1)=(4d0/3d0)*PI*rhoCenters(1)*rCenters(1)**3 
                    gCenters(1)=Gconst*mCenters(1)/rCenters(1)**2
                    DO ig=2,N+1
                        ! density jumps:
                        IF(abs(rhoCenters(ig)-rhoCenters(ig-1)).gt.JTINY.and.ig.ne.2) THEN
                            tmpdr=rJumps(1+jcntr)-rCenters(ig-1)
                            mShell=(4d0/3d0)*PI*rhoCenters(ig-1)*((tmpdr+rCenters(ig-1))**3-rCenters(ig-1)**3) 
                            mShell=mShell+(4d0/3d0)*PI*rhoCenters(ig)*((rCenters(ig)**3-(rCenters(ig-1)+tmpdr)**3))
                            mCenters(ig)=mShell+mCenters(ig-1) 
                            gCenters(ig)=Gconst*mCenters(ig)/(rCenters(ig)**2)
                            jcntr=jcntr+1
                        ! homogeneous chunks:
                        ELSE 
                            mShell=(4d0/3d0)*PI*rhoCenters(ig)*(rCenters(ig)**3-rCenters(ig-1)**3)
                            mCenters(ig)=mShell+mCenters(ig-1)
                            gCenters(ig)=Gconst*mCenters(ig)/rCenters(ig)**2
                        ENDIF  
                    ENDDO
                    jcntr=1
                    ! gravitational acceleration evaluated at interfaces:
                    mFaces(1)=(4d0/3d0)*PI*rhoCore*rFaces(1)**3
                    gFaces(1)=Gconst*mFaces(1)/rFaces(1)**2
                    DO ig=2,N 
                        ! density jump
                        IF(abs(rhoFaces(ig)-rhoFaces(ig-1)).gt.JTINY.and.ig.ne.2) THEN
                            tmpdr=rJumps(1+jcntr)-rFaces(ig-1)
                            mShell=(4d0/3d0)*PI*rhoFaces(ig-1)*((tmpdr+rFaces(ig-1))**3-rFaces(ig-1)**3)
                            mShell=mShell+(4d0/3d0)*PI*rhoFaces(ig)*((rFaces(ig)**3-(rFaces(ig-1)+tmpdr)**3))
                            mFaces(ig)=mShell+mFaces(ig-1)
                            gFaces(ig)=Gconst*mFaces(ig)/rFaces(ig)**2
                            jcntr=jcntr+1
                        ! homogeneous chunks:
                        ELSE
                            mShell=(4d0/3d0)*PI*rhoFaces(ig)*(rFaces(ig)**3-rFaces(ig-1)**3)
                            mFaces(ig)=mShell+mFaces(ig-1)
                            gFaces(ig)=Gconst*mFaces(ig)/rFaces(ig)**2
                        ENDIF
                    ENDDO
                !% Numerical integration for PREM
                ELSEIF(MODEL.eq.3) THEN
                    mCenters(1)=(4d0/3d0)*PI*rhoCenters(1)*rCenters(1)**3 
                    mFaces(1)=(4d0/3d0)*PI*rhoCore*rCenters(1)**3 
                    gCenters(1)=Gconst*mCenters(1)/rCenters(1)**2
                    gFaces(1)=Gconst*mFaces(1)/rFaces(1)**2
                    DO ig=2,N+1
                        rhoShell=(rhoCenters(ig)+rhoCenters(ig-1))/2d0
                        mShell=(4d0/3d0)*PI*rhoShell*(rCenters(ig)**3-rCenters(ig-1)**3)
                        mCenters(ig)=mShell+mCenters(ig-1)
                        gCenters(ig)=Gconst*mCenters(ig)/rCenters(ig)**2
                    ENDDO
                    DO ig=2,N
                        rhoShell=(rhoFaces(ig)-rhoFaces(ig-1))/2d0
                        mShell=(4d0/3d0)*PI*(rFaces(ig)**3)*rhoShell-(4d0/3d0)*PI*(rFaces(ig-1)**3)*rhoShell
                        mFaces(ig)=mShell+mFaces(ig-1)
                        gFaces(ig)=Gconst*mFaces(ig)/rFaces(ig)**2
                    ENDDO
                ENDIF
        END SUBROUTINE init_g0m0
    
        ! INITIALIZE PREM MODEL:
        SUBROUTINE init_PREM(rCenters,rFaces,rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces)
            REAL(8),DIMENSION(N),INTENT(IN)      :: rFaces
            REAL(8),DIMENSION(N+1),INTENT(IN)    :: rCenters
            REAL(8),DIMENSION(N+1),INTENT(INOUT) :: rhoCenters
            REAL(8),DIMENSION(N),INTENT(INOUT)   :: rhoFaces,muFaces,etaFaces,lmdFaces
            REAL(8),DIMENSION(6)                 :: tmp
            REAL(8),DIMENSION(:),ALLOCATABLE     :: rPrem(:),rhoPrem(:),muPrem(:),KPrem(:),etaPrem(:),lmdPrem(:)
            REAL(8)                              :: core,rhoFerr,rhoCerr,muerr,etaerr,lmderr
            LOGICAL                              :: findcore=.true.
            INTEGER                              :: lcntr=0,i=1,IERR=1,ir=2,irho=3,imu=4,iK=5,ieta=6
            INTEGER                              :: j,k,from,to
                !% PREM.dat: 1.num. interface 2.radius, 3.density, 4.shear modulus, 5.viscosity, 6.bulk modulus
                OPEN(55,file='PREM.dat')
                IERR=1
                DO WHILE(IERR.ge.0)
                    !% end of the file, IERR changes from 0 to -1
                    READ(55,fmt=*,iostat=IERR) tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6)
                    IF(tmp(2).ge.RBOT) THEN
                    IF(findcore) THEN; findcore=.false.; core=tmp(1); ENDIF
                    lcntr=lcntr+1
                    ENDIF
                ENDDO
                REWIND 55
                ALLOCATE(rPrem(lcntr-1),rhoPrem(lcntr-1),muPrem(lcntr-1),KPrem(lcntr-1),etaPrem(lcntr-1),&
                lmdPrem(lcntr-1),SOURCE=0d0)
                IERR=1
                DO WHILE(IERR.ge.0) 
                    READ(55,fmt=*,iostat=IERR) tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6)
                    IF(tmp(1).ge.core) THEN 
                        rPrem(i)=tmp(ir)
                        rhoPrem(i)=tmp(irho)     
                        muPrem(i)=tmp(imu)
                        KPrem(i)=tmp(iK)
                        etaPrem(i)=tmp(ieta)
                        i=i+1
                    ENDIF
                ENDDO    
                i=i-1
                print *, i
                lmdPrem=KPrem-(2d0/3d0)*muPrem
                !% Using polint to interpolate prem.DAT quantities onto the computational grid
                IF(N.lt.size(rPrem)) STOP "Denser grid required."
                IF(rFaces(N).gt.rPrem(lcntr-1)) STOP "Invalid grid."
                DO j=2,N-1 ! grid
                    DO k=2,i; 
                      IF(rFaces(j).lt.rPrem(k)) EXIT
                    ENDDO; ! finds corrensonding place on the PREM grid
                    print *, k, j
                    IF(k.lt.i/2) THEN
                        from=k-1; to=k+NPOLINT-2
                    ELSE
                        from=k-NPOLINT+1; to=k
                    ENDIF
                    !% NPOLINT=2 gives linear interpolation
                    CALL polint(rPrem(from:to),rhoPrem(from:to),NPOLINT,rFaces(j),rhoFaces(j),rhoFerr)
                    CALL polint(rPrem(from:to),rhoPrem(from:to),NPOLINT,rCenters(j),rhoCenters(j),rhoCerr)
                    CALL polint(rPrem(from:to),muPrem(from:to),NPOLINT,rFaces(j),muFaces(j),muerr)
                    CALL polint(rprem(from:to),etaPrem(from:to),NPOLINT,rFaces(j),etaFaces(j),etaerr)
                    CALL polint(rPrem(from:to),lmdPrem(from:to),NPOLINT,rFaces(j),lmdFaces(j),lmderr)
                ENDDO
                !% Fixing boundary points & HOMOGENEOUS CORE
                rhoFaces(1)=rhoFaces(2); rhoFaces(N)=rhoFaces(N-1)
                muFaces(1)=muFaces(2); muFaces(N)=muFaces(N-1)
                etaFaces(1)=etaFaces(2); etaFaces(N)=etaFaces(N-1)
                lmdFaces(1)=lmdFaces(2); lmdFaces(N)=lmdFaces(N-1)
                rhoCenters(1)=rhoCore; rhoCenters(N+1)=rhoCenters(N-1); rhoCenters(N)=rhoCenters(N-1)
                CLOSE(55)
        END SUBROUTINE init_PREM

        ! SELF-GRAVITY ITERATION SUBROUTINE:
        SUBROUTINE sg_iter(A,al,indx,b,rhs,ur_in,ut2_in,urtop,urbot,r,rho,drhobot,drhotop,cmat,trhs,jindx)
            INTEGER,DIMENSION(noj),INTENT(IN)      :: jindx
            INTEGER,DIMENSION(nn),INTENT(IN)       :: indx 
            REAL(8),DIMENSION(nn,m),INTENT(IN)     :: A
            REAL(8),DIMENSION(nn,m1),INTENT(IN)    :: al
            REAL(8),DIMENSION(10,10),INTENT(IN)    :: cmat 
            REAL(8),DIMENSION(N+1),INTENT(INOUT)   :: ur_in
            REAL(8),DIMENSION(N),INTENT(INOUT)     :: ut2_in
            REAL(8),DIMENSION(nn),INTENT(INOUT)    :: b,rhs
            REAL(8),DIMENSION(nn),INTENT(INOUT)    :: r,rho
            REAL(8),POINTER,INTENT(IN)             :: trhs
            REAL(8),INTENT(IN)                     :: drhotop,drhobot
            REAL(8),INTENT(INOUT)                  :: urtop,urbot
            REAL(8),DIMENSION(N+1)                 :: ur
            REAL(8),DIMENSION(N)                   :: ur2,ut2
            INTEGER                                :: i,l,k
            REAL(8)                                :: sgcorr,tmp1,tmp2,urjump,tmpdrho
                !% todo model 3
                ur=ur_in; ut2=ut2_in;
                DO WHILE(abs(urbot-tmp1).gt.EPS.and.abs(urtop-tmp2).gt.EPS) !% iteration
                    tmp1=urbot
                    tmp2=urtop
                    IF(SGPOT) THEN !% sg potential correction
                        sgcorr=sgint(RBOT,urtop,JMAX,RTOP,drhotop)+sgext(RBOT,urbot,JMAX,RBOT,drhobot) 
                        IF(SGJUMP.and.MODEL.eq.2) THEN 
                            DO l=1,noj !% loop over jumps
                                urjump=ur(jindx(l)+1) 
                                tmpdrho=rhojumps(l+1)-rhojumps(l+2) ! >0
                                sgcorr=sgcorr+sgint(RBOT,urjump,JMAX,r(jindx(l)+1),tmpdrho)                               
                            END DO
                        ELSEIF(SGJUMP.and.MODEL.eq.3) THEN
                            DO i=2,N !% centers of layers
                                IF((rho(i+1)-rho(i)).gt.JTINY) THEN
                                    urjump=ur(i)
                                    tmpdrho=rho(i+1)-rho(i)
                                    sgcorr=sgcorr+sgint(RBOT,urjump,JMAX,r(i),tmpdrho)
                                ENDIF
                            ENDDO
                        ENDIF
                        IF(CORE_PRESSURE) THEN
                            rhs(2)=-cmat(1,3)*rhoCore*real(obltidepot(trhs,RBOT)+sgcorr)
                            rhs(3)=-cmat(1,4)*rhoCore*real(obltidepot(trhs,RBOT)+sgcorr)
                        ENDIF
                    ENDIF
                    IF(SGFORCE) THEN
                        DO k=2,N   !% top and bottom topography contribution
                            rhs((k-1)*6+2)=-rho(k)*sgfext(r(k),urbot,JMAX,RBOT,drhobot) ! Y_j+1
                            rhs((k-1)*6+3)=rho(k)*real(obltidef(trhs,r(k)))&
                            -rho(k)*sgfint(r(k),urtop,JMAX,RTOP,drhotop) ! Y_j-1
                        END DO 
                        IF(SGJUMPF.and.MODEL.eq.2) THEN
                            DO l=1,noj  ! loop over jumps
                                urjump=ur(jindx(l)+1) 
                                tmpdrho=rhojumps(l+1)-rhojumps(l+2) ! >0
                                DO k=2,N     ! loop over centers of layers
                                    IF(k.eq.jindx(l)+1) THEN
                                        rhs((k-1)*6+2)=rhs((k-1)*6+2)-rho(k)*sgfext(r(k),urjump,&
                                        JMAX,r(jindx(l)+1),tmpdrho) ! Y_j+1 <---- QUESTIONABLE!!
                                    ELSEIF(r(k).lt.r(jindx(l)+1)) THEN
                                        rhs((k-1)*6+3)=rhs((k-1)*6+3)-rho(k)*sgfint(r(k),urjump,&
                                        JMAX,r(jindx(l)+1),tmpdrho) ! Y_j-1 
                                    ELSEIF(r(k).gt.r(jindx(l)+1)) THEN
                                        rhs((k-1)*6+2)=rhs((k-1)*6+2)-rho(k)*sgfext(r(k),urjump,&
                                        JMAX,r(jindx(l)+1),tmpdrho) ! Y_j+1
                                    ELSE; STOP "Unexpected error in sg_iter (model 2).";
                                    ENDIF
                                END DO
                            END DO
                        ELSEIF(SGJUMPF.and.MODEL.eq.3) THEN
                            DO i=2,N !% centers of layers
                                IF((rho(i+1)-rho(i)).gt.JTINY) THEN
                                    urjump=ur(i)
                                    tmpdrho=rho(i+1)-rho(i)
                                    DO k=2,N     ! loop over centers of layers
                                        IF(k.eq.i) THEN
                                            rhs((k-1)*6+2)=rhs((k-1)*6+2)-rho(k)*sgfext(r(k),urjump,&
                                            JMAX,r(i),tmpdrho) ! Y_j+1 <---- QUESTIONABLE!!
                                        ELSEIF(r(k).lt.r(i)) THEN
                                            rhs((k-1)*6+3)=rhs((k-1)*6+3)-rho(k)*sgfint(r(k),urjump,&
                                            JMAX,r(i),tmpdrho) ! Y_j-1 
                                        ELSEIF(r(k).gt.r(i)) THEN
                                            rhs((k-1)*6+2)=rhs((k-1)*6+2)-rho(k)*sgfext(r(k),urjump,&
                                            JMAX,r(i),tmpdrho) ! Y_j+1
                                        ELSE; STOP "Unexpected error in sg_iter (model 3).";
                                        ENDIF
                                    END DO
                                ENDIF
                            ENDDO
                        ENDIF
                    ENDIF !% end of self-gravity
                    b=rhs   !% beware, banbks rewrite this variable!
                    CALL banbks(A,nn,m1,m2,np,mp,al,mpl,indx,b)
                    ur(:)=cmat(1,3)*b(1::6)+cmat(1,4)*b(2::6)
                    ur2(:)=(ur(1:N)+ur(2:N+1))/2d0
                    ut2(:)=-cmat(1,4)*(b(1:6*N-5:6)+b(7::6))/2d0+cmat(1,3)*(b(2:6*N-4:6)+b(8::6))/2d0
                    urbot=ur2(1)
                    urtop=ur2(N)
                END DO !% end of the while loop
                ur_in=ur
                ut2_in=ut2
        END SUBROUTINE sg_iter
END MODULE mFunc   
