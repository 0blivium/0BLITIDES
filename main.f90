! AUTHOR: Bc. Vít Beran
! CLASS: Mechanics of Continuum 2 (NGEO069), Charles University, Prague.
! PROGRAM: Response of a compressible viscoleastic shell to the tidal potential.
! COMPILATION: gfortran -O3 main.f90 mCoef.f90 mEqs.f90 mFunc.f90 mPar.f90 sPolint.f90 sBan* 

PROGRAM main
    USE mPar
    USE mFunc 
    USE mCoeff
    USE mEqs
    IMPLICIT NONE
        EXTERNAL sBandec,sBanbks,sPolint   !% numerical recipes   
        TYPE(LoveT)                         :: LnumT
        REAL(8)                             :: dt=100.*P
        REAL(8),DIMENSION(:,:),ALLOCATABLE  :: A(:,:),Asvd(:,:),al(:,:)
        REAL(8),DIMENSION(:),ALLOCATABLE    :: rFaces(:),rhoFaces(:),muFaces(:),etaFaces(:),gFaces(:),mFaces(:),lmdFaces(:)
        REAL(8),DIMENSION(:),ALLOCATABLE    :: rCenters(:),rhoCenters(:),gCenters(:),mCenters(:)
        REAL(8),DIMENSION(:),ALLOCATABLE    :: rhsRe(:),rhsIm(:),bRe(:),bIm(:),ur(:),urIm(:),ur2(:),ur2Im(:),ut2(:),mxwlTime(:)
        INTEGER,DIMENSION(:),ALLOCATABLE    :: indx(:)
        INTEGER,DIMENSION(noj)              :: jindx
        REAL(8),DIMENSION(noj+2)            :: drhojumps=0d0
        REAL(8),POINTER                     :: trhs
        REAL(8),TARGET                      :: t0=0d0,t=0d0
        REAL(8)                             :: c(N+1)=0d0,f(N)=0d0
        REAL(8)                             :: urbot,urtop,urbotIm,urtopIm,uttop,utbot
        REAL(8)                             :: dr,d,drpsg=0d0,Vsurf,drhotop,drhobot,urjump,num,tmpdrho
        REAL(8)                             :: c1,c2,tmp,tmpeta,Tmid
        REAL(8)                             :: dSol,dur,urold,mxwlMax,mxwlMin,dtnext
        INTEGER                             :: simpcntr=0,steps=1,mindx=1,lindx=1
        INTEGER                             :: i,j,l,k1,k2,rl
        LOGICAL                             :: step1=.true.

        !% Commence with the run
        CALL SYSTEM('mkdir -p run/')

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%% SPATIAL DISCRETIZATION AND MATRIX CONSTRUCTION %%
        ALLOCATE(rCenters,rhoCenters,gCenters,mCenters,ur,urIm,SOURCE=c)
        ALLOCATE(rFaces,rhoFaces,gFaces,mFaces,muFaces,etaFaces,lmdFaces,mxwlTime,SOURCE=f)
        ALLOCATE(ur2,ut2,ur2Im,SOURCE=f)
        ALLOCATE(A(NN,m),al(NN,m1),Asvd(NN,m),rhsRe(NN),rhsIm(NN),bRe(NN),bIm(NN),SOURCE=0d0)
        ALLOCATE(indx(NN),SOURCE=0)
    
        SELECT CASE(MODEL)
            CASE(1) !% Homogeneous model
                CALL init_rad(dr,rCenters,rFaces)
                CALL init_homog_model(rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces)
                CALL init_g0m0(rCenters,rhoCenters,gCenters,rFaces,rhoFaces,gFaces,mFaces,mCenters)
                ! overrides init_g0m0
                IF(GHOMO) gCenters=gHmg; gFaces=gHmg; 
                IF(ETAEXP) THEN  ! Arrhenius law implemented for the homogeneous model, overrides constant viscosity
                    IF(TMODE.eq.1) THEN ! conductive profile
                        c1=(RBOT*RTOP*(TBOT-TTOP))/(RTOP-RBOT)
                        c2=TBOT-c1/RBOT
                        DO i=1,N
                            tmp=c1/rFaces(i)+c2
                            tmpeta=ETA_BASE*exp((E/(R*TBOT))*(TBOT/tmp-1d0))
                            IF(tmpeta.lt.eta_cut) THEN; etaFaces(i)=tmpeta; ELSE; etaFaces(i:)=ETA_CUT; EXIT; ENDIF;
                        END DO
                    ELSEIF(TMODE.eq.2) THEN ! convective profile
                        k1=nint(0.1*N)
                        k2=nint(2d0/3d0*N)
                        Tmid=250d0
                        !IF(Tmid.gt.TTOP.or.Tmid.lt.TBOT) THEN; STOP "Adjust Tmid."; ENDIF;
                        c1=(Tmid*RBOT-TBOT*rFaces(k1))/(RBOT-rFaces(k1))
                        c2=(Tmid-c1)/rFaces(k1)
                        DO i=1,k1
                            tmp=c2*rFaces(i)+c1
                            tmpeta=ETA_BASE*exp((E/(R*TBOT))*(TBOT/tmp-1d0))
                            IF(tmpeta.lt.eta_cut) THEN; etaFaces(i)=tmpeta; ELSE; etaFaces(i:)=ETA_CUT; EXIT; ENDIF;
                        ENDDO
                        etaFaces(k1+1:k2)=etaFaces(k1)
                        c1=(rFaces(k2)*RTOP*(Tmid-TTOP))/(rFaces(k2)-RBOT)
                        c2=Tmid-c1/rFaces(k2)
                        DO i=k2,N 
                            tmp=c1/rFaces(i)+c2
                            tmpeta=etaFaces(k1)*exp((E/(R*Tmid))*(Tmid/tmp-1d0))
                            IF(tmpeta.lt.eta_cut) THEN; etaFaces(i)=tmpeta; ELSE; etaFaces(i:)=ETA_CUT; EXIT; ENDIF;
                        ENDDO
                    ENDIF
                ENDIF
            CASE(2) !% Layered model
                CALL init_rad(dr,rCenters,rFaces)
                CALL init_layered_model(rCenters,rFaces,rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces,jindx,dr)
                CALL init_g0m0(rCenters,rhoCenters,gCenters,rFaces,rhoFaces,gFaces,mFaces,mCenters)
                !% if a density jump occurs, the value is overwritten by an average value to capture the density gradient more accuratelly
                DO i=1,noj; rhoCenters(jindx(i)+1)=(rhoJumps(i+1)+rhoJumps(i+2))/2d0; ENDDO
                DO l=1,noj+1; drhojumps(l)=rhojumps(l)-rhojumps(l+1); ENDDO
                drhojumps(noj+2)=rhojumps(noj+2)
            CASE(3) !% PREM model - simplified model with homogeneous liquid core 
                CALL init_rad(dr,rCenters,rFaces)
                CALL init_PREM(rCenters,rFaces,rhoCenters,rhoFaces,muFaces,etaFaces,lmdFaces)
                CALL init_g0m0(rCenters,rhoCenters,gCenters,rFaces,rhoFaces,gFaces,mFaces,mCenters)
            CASE DEFAULT
                STOP 'Invalid model selected.'
        END SELECT

        OPEN(44,file='run/output_evol.dat')
        drhobot=rhoCenters(1)-rhoCenters(2) ! density jump on the bottom
        drhotop=rhoCenters(N+1)             ! density jump on the surface

        !% Constant amplitude vs periodic tides
        SELECT CASE(DEFORM)
            CASE(1); trhs=>t0;
            CASE(2); trhs=>t;
            CASE DEFAULT; STOP 'Invalid deformation response.'
        END SELECT

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%% SPECTRAL PROBLEM %%
        CALL init_coeff(JMAX,gCenters,rhoCenters)
        CALL init_matrix(A,cmat,lindx,mindx,dr,rFaces,rhoCenters,rhoFaces,gCenters,gFaces,muFaces,lmdFaces)
        IF(MATRIX_OUTPUT) THEN; CALL print_mat(A,nn,m); ENDIF;
    
        !% boundary condition (by default dynamic topography)
        IF(NOSLIP.eqv..true..and.CORE_PRESSURE.eqv..true.) STOP 'Undetermined bottom boundary condition.'
        IF(NOSLIP) THEN 
            A(2,:)=0d0; A(3,:)=0d0 ! overrides BC
            A(2,od-1)=1d0/2d0     
            A(2,od+5)=1d0/2d0           
            A(3,od-1)=1d0/2d0
            A(3,od+5)=1d0/2d0 
        ELSEIF(CORE_PRESSURE) THEN  ! bottom boundary condition (dynamic topography & core pressure)
            !% real component
            rhsRe(2)=-cmat(1,3)*rhoCore*real(obltidepot(0d0,RBOT))
            rhsRe(3)=-cmat(1,4)*rhoCore*real(obltidepot(0d0,RBOT))
            !% complex component
            rhsIm(2)=-cmat(1,4)*rhoCore*aimag(obltidepot(0d0,RBOT))
            rhsIm(3)=-cmat(1,4)*rhoCore*aimag(obltidepot(0d0,RBOT))
        ENDIF

        print *, rCenters(N), obltidef(0d0,rCenters(N))

        !%% Initialize the right-hand side 
        DO i=2,N ! Y_j-1
            rhsRe((i-1)*6+3)=rhoCenters(i)*real(obltidef(0d0,rCenters(i)))     !% real component 
            rhsIm((i-1)*6+3)=rhoCenters(i)*aimag(obltidef(0d0,rCenters(i)))    !% imaginary component 
        END DO 
        
        bRe=rhsRe; bIm=rhsIm; Asvd=A;
        ! LU decomposition:
        CALL bandec(A,nn,m1,m2,np,mp,al,mpl,indx,d)
        !% real system & complex system
        CALL banbks(A,nn,m1,m2,np,mp,al,mpl,indx,bRe)
        CALL banbks(A,nn,m1,m2,np,mp,al,mpl,indx,bIm)
        !%% Radial and Vertical displacement %%
        ur(:)=cmat(1,3)*bRe(1::6)+cmat(1,4)*bRe(2::6)
        ur2(:)=(ur(1:N)+ur(2:N+1))/2d0
        urIm(:)=cmat(1,3)*bIm(1::6)+cmat(1,4)*bIm(2::6)
        ur2Im(:)=(urIm(1:N)+urIm(2:N+1))/2d0
        ut2(:)=cmat(1,3)*(bRe(2:6*N-4:6)+bRe(8::6))/2d0-cmat(1,4)*(bRe(1:6*N-5:6)+bRe(7::6))/2d0
        urbot=ur2(1); urbotIm=ur2Im(1); utbot=ut2(1); 
        urtop=ur2(N); urtopIm=ur2Im(N); uttop=ut2(N);
        !%% Iterate the solution for self-gravity until the desired precision is achieved
        CALL sg_iter(A,al,indx,bRe,rhsRe,ur,ut2,urtop,urbot,rCenters,rhoCenters,drhobot,drhotop,cmat,trhs,jindx)
        !%% Surface deformation
        uttop=ut2(N); utbot=ut2(1);
        Vsurf=sgext(RTOP,urbot,JMAX,RBOT,drhobot)+sgint(RTOP,urtop,JMAX,RTOP,drhotop)  
        IF(MODEL.eq.2) THEN
            DO l=1,noj ! loop over density jumps
                urjump=ur(jindx(l)+1)     
                tmpdrho=rhojumps(l+1)-rhojumps(l+2)
                Vsurf=Vsurf+sgext(RTOP,urjump,JMAX,rCenters(jindx(l)+1),tmpdrho)
            ENDDO
        ELSEIF(MODEL.eq.3) THEN
            DO i=2,N !% centers of layers
                IF((rhoCenters(i+1)-rhoCenters(i)).gt.JTINY) THEN
                    urjump=ur(i)
                    tmpdrho=rhoCenters(i+1)-rhoCenters(i)
                    Vsurf=Vsurf+sgext(RTOP,urjump,JMAX,rCenters(i),tmpdrho)
                ENDIF
            ENDDO
        ENDIF
        IF(DEFORM.eq.1) THEN
            !% Tidal Love numbers & driving+sg potential [m]
            LnumT%k2T=real(Vsurf/(obltidepot(trhs,RTOP)))
            LnumT%h2T=real((gCenters(N+1)*urtop))/(abs(obltidepot(trhs,RTOP)))
            LnumT%l2T=real((cmat(1,6)*gCenters(N+1)*uttop)/(abs(obltidepot(trhs,RTOP))))
            drpsg=real(obltidepot(trhs,RTOP)+Vsurf)/gCenters(N+1)
            PRINT '(a60,f18.14)','Elastic limit of degree 2 tidal Love number k2e:',LnumT%k2T
            PRINT '(a60,f18.14)','Elastic limit of degree 2 tidal Love number h2e:',LnumT%h2T
            PRINT '(a60,f18.14)','Elastic limit of degree 2 tidal Love number l2e:',LnumT%l2T
        ENDIF

        !% 1st row: elastic solution at t=0
        IF(TIMEEVOL) WRITE(44,"(10e30.15)") 0d0, 0d0, urtop, urtopIm, urbot, urbotIm, LnumT, drpsg

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%  TIME EVOLUTION  %%    
        IF(TIMEEVOL) THEN 
            IF(RULE.eq."simpson".and.TSTEP.eq.2) STOP 'Simpson rule is not implemented for non-constant time step!'
            mxwlMin=minval(etaFaces/muFaces)/YTOSEC
            mxwlMax=maxval(etaFaces/muFaces)/YTOSEC
            IF(dt/YTOSEC.gt.mxwlMin) STOP 'Your time step is too large!'
            ! PRINT *, mxwlMin/YTOSEC, mxwlMax/YTOSEC
            IF(VISCOEL) THEN
                A=Asvd
                !% modification of the matrix
                SELECT CASE(RULE)
                    CASE('trapez')
                        num=4d0
                    CASE('simpson')
                        num=6d0
                    CASE DEFAULT 
                        STOP 'Unrecognized integration rule.'
                END SELECT
                DO i=1,N  
                    A((i-1)*6+4,od+0)=A((i-1)*6+4,od+0)-dt/(num*etaFaces(i))
                    A((i-1)*6+5,od+0)=A((i-1)*6+5,od+0)-dt/(num*etaFaces(i))
                    A((i-1)*6+6,od+0)=A((i-1)*6+6,od+0)-dt/(num*etaFaces(i))
                ENDDO
                Asvd=A
                CALL bandec(A,nn,m1,m2,np,mp,al,mpl,indx,d) !% -> new A,al,indx
            ENDIF
            !%% Main time loop %%
            DO WHILE(t.lt.TMAX.and.steps.lt.MAXSTEPS)                 
                steps=steps+1
                dt=dt/YTOSEC ! [sec]->[yr]
                SELECT CASE(TSTEP)
                    CASE(1) ! constant time step
                        t=t+dt
                    CASE(2) ! setting length of the time step (with respect to urtop)
                        IF(mod(steps,10).eq.2.and.steps.gt.9) THEN
                            dSol=abs(1.-dur*dt/(urtop-urold))
                            IF(dSol.lt.UPTOL) THEN
                                dtnext=min(dt*up,mxwlMin)
                            ELSEIF(dSol.gt.DOWNTOL) THEN
                                dtnext=dt/DOWN
                            ELSE; dtnext=dt;
                            ENDIF
                            dur=(urtop-urold)/dt
                        ELSE
                            dtnext=dt
                            IF(steps.eq.2) dur=(urtop-urold)/dt
                        ENDIF
                        dt=dtnext
                        urold=urtop
                        t=t+dt
                    CASE DEFAULT
                        STOP 'Invalid time step adjustment method.'
                END SELECT
                dt=dt*YTOSEC ! [yr] -> [sec]
                DO i=1,N-1  !% update of the right-hand-side!
                    rhsRe(i*6+3)=rhoCenters(i+1)*real(obltidef(trhs,rCenters(i+1))) 
                    rhsIm(i*6+3)=rhoCenters(i+1)*aimag(obltidef(trhs,rCenters(i+1)))
                    IF(VISCOEL) THEN
                        !% real part
                        rhsRe((i-1)*6+4)=rhsRe((i-1)*6+4)+(bRe((i-1)*6+4)*dt)/(num*etaFaces(i)) ! tj-2,2
                        rhsRe((i-1)*6+5)=rhsRe((i-1)*6+5)+(bRe((i-1)*6+5)*dt)/(num*etaFaces(i)) ! tj,2
                        rhsRe((i-1)*6+6)=rhsRe((i-1)*6+6)+(bRe((i-1)*6+6)*dt)/(num*etaFaces(i)) ! tj+2,2
                        !% complex part
                        rhsIm((i-1)*6+4)=rhsIm((i-1)*6+4)+(bIm((i-1)*6+4)*dt)/(num*etaFaces(i)) ! tj-2,2
                        rhsIm((i-1)*6+5)=rhsIm((i-1)*6+5)+(bIm((i-1)*6+5)*dt)/(num*etaFaces(i)) ! tj,2
                        rhsIm((i-1)*6+6)=rhsIm((i-1)*6+6)+(bIm((i-1)*6+6)*dt)/(num*etaFaces(i)) ! tj+2,2
                    ELSEIF(CORE_PRESSURE) THEN
                        rhsRe(2)=-cmat(1,3)*rhoCore*real(obltidepot(trhs,RBOT)) 
                        rhsRe(3)=-cmat(1,4)*rhoCore*real(obltidepot(trhs,RBOT))
                        rhsIm(2)=-cmat(1,4)*rhoCore*aimag(obltidepot(trhs,RBOT))
                        rhsIm(3)=-cmat(1,4)*rhoCore*aimag(obltidepot(trhs,RBOT))
                    ENDIF
                ENDDO
                IF(RULE.eq.'trapez') THEN
                    IF(step1) THEN
                        num=2d0
                        step1=.false.
                    ENDIF
                ELSEIF(RULE.eq.'simpson') THEN
                    simpcntr=simpcntr+1
                    IF(mod(simpcntr,2).eq.0) THEN
                        num=(3d0/1d0) !% even terms
                    ELSE
                        num=(3d0/2d0) !% odd terms
                    ENDIF
                ELSE
                    STOP 'Unrecognized integration rule.'
                ENDIF
                bRe=rhsRe
                bIm=rhsIm
                CALL banbks(A,nn,m1,m2,np,mp,al,mpl,indx,bRe)
                !% Skip the complex problem for calculating response to Heviside
                IF(DEFORM.eq.2) THEN
                    CALL banbks(A,nn,m1,m2,np,mp,al,mpl,indx,bIm)
                ENDIF
                ur(:)=cmat(1,3)*bRe(1::6)+cmat(1,4)*bRe(2::6)
                ur2(:)=(ur(1:N)+ur(2:N+1))/2d0
                urIm(:)=cmat(1,3)*bIm(1::6)+cmat(1,4)*bIm(2::6)
                ur2Im(:)=(urIm(1:N)+urIm(2:N+1))/2d0
                ut2(:)=cmat(1,3)*(bRe(2:6*N-4:6)+bRe(8::6))/2d0-cmat(1,4)*(bRe(1:6*N-5:6)+bRe(7::6))/2d0
                urbot=ur2(1); urbotIm=ur2Im(1); utbot=ut2(1); 
                urtop=ur2(N); urtopIm=ur2Im(N); uttop=ut2(N);
                !%% Iterate the solution until a desired precision is achieved
                CALL sg_iter(A,al,indx,bRe,rhsRe,ur,ut2,urtop,urbot,rCenters,rhoCenters,drhobot,drhotop,cmat,trhs,jindx)
                !% Update vertical displacement at the top/bottom
                utbot=ut2(1); uttop=ut2(N);
                Vsurf=sgext(RTOP,urbot,JMAX,RBOT,drhobot)+sgint(RTOP,urtop,JMAX,RTOP,drhotop)  
                IF(MODEL.eq.2) THEN
                    DO l=1,noj !% loop over density jumps
                        urjump=ur(jindx(l)+1)     
                        tmpdrho=rhojumps(l+1)-rhojumps(l+2)
                        Vsurf=Vsurf+sgext(RTOP,urjump,JMAX,rCenters(jindx(l)+1),tmpdrho)
                    ENDDO 
                ELSEIF(MODEL.eq.3) THEN
                    DO i=2,N !% centers of layers
                        IF((rhoCenters(i+1)-rhoCenters(i)).gt.JTINY) THEN
                            urjump=ur(i)
                            tmpdrho=rhoCenters(i+1)-rhoCenters(i)
                            Vsurf=Vsurf+sgext(RTOP,urjump,JMAX,rCenters(i),tmpdrho)
                        ENDIF
                    ENDDO
                ENDIF
                IF(DEFORM.eq.1) THEN
                    !%% Tidal Love numbers & driving+sg potential [m]
                    LnumT%k2T=real(Vsurf/(obltidepot(trhs,RTOP)))
                    LnumT%h2T=real((gCenters(N+1)*urtop)/abs(obltidepot(trhs,RTOP)))
                    LnumT%l2T=real((cmat(1,6)*gCenters(N+1)*uttop)/(abs(obltidepot(trhs,RTOP))))
                    drpsg=real(obltidepot(trhs,RTOP)+Vsurf)/gCenters(N+1)
                ENDIF
                PRINT '(a20,f10.5,a1)', 'Time evolution:',(t/TMAX)*1.e2,'%'
                WRITE (44,"(10e30.15)") t, dt/YTOSEC, urtop, urtopIm, urbot, urbotIm, LnumT, drpsg
            END DO
            IF(steps.gt.MAXSTEPS) STOP 'The maximum number of steps has been exceeded.'
        ENDIF !% end of time evolution branch
        CLOSE(44) 
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%% MODEL OUTPUT %%
    
        IF(MODEL_OUTPUT) THEN
            ! reverse the arrays for benchmarking with reversed models
            IF(MODEL_REVERSE) THEN
                CALL arr_reverse(rFaces)
                CALL arr_reverse(muFaces)
                CALL arr_reverse(etaFaces)
                CALL arr_reverse(lmdFaces)
                CALL arr_reverse(rhoFaces)
                CALL arr_reverse(gFaces)
                CALL arr_reverse(rCenters)
                CALL arr_reverse(rhoCenters)
                CALL arr_reverse(gCenters)
                CALL arr_reverse(mCenters)
            ENDIF
            IF(MODEL.eq.2) THEN
                OPEN(22,file='run/log.dat')
                WRITE (22,*) 'Layers with jumps:',jindx
                WRITE (22,*) 'Radii with density jumps:',rjumps
                WRITE (22,*) 'Density jumps:',drhojumps
                CLOSE(22)
            ENDIF
            OPEN(1,file='run/gridfaces.txt')
            OPEN(2,file='run/gridcenters.txt')
            DO rl=1,N
                ! 1: radius 2: density 3: mass 4: gravity 5: material parameters
                WRITE (1,'(7e30.17)') (rFaces(rl),rhoFaces(rl),mFaces(rl),gFaces(rl),muFaces(rl),etaFaces(rl),lmdFaces(rl),j=1,1)
                WRITE (2,'(4e30.17)') (rCenters(rl),rhoCenters(rl),mCenters(rl),gCenters(rl),j=1,1)
            END DO
            WRITE (2,'(4e30.17)') rCenters(N+1),rhoCenters(N+1),mCenters(N+1),gCenters(N+1)
            CLOSE(1)
            CLOSE(2)
        ENDIF  
        PRINT *, 'The program terminated successfully.'
END PROGRAM main
