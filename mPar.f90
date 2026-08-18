! MODULE: Global constants, parameters, and attributes (on/off) of the model.

MODULE mPar
    IMPLICIT NONE

    ! ADDITIONAL DATA TYPES:
    TYPE LoveT; REAL(8) :: k2T=0d0,h2T=0d0,l2T=0d0; END TYPE LoveT;

    ! MATHEMATICAL AND PHYSICAL CONSTANTS:
    REAL(8),PARAMETER :: PI            = atan(1d0)*4d0           
    REAL(8),PARAMETER :: Gconst        = 6.6732d-11                                                      ! gravitational constant
    REAL(8),PARAMETER :: EPS           = 1d-9                                                            ! 8-bit real variables comparison
    REAL(8),PARAMETER :: JTINY         = 1d-1                                                            ! 8-bit real parameter jumps identifier (pretty much arbitrary)                                                    
    REAL(8),PARAMETER :: YTOSEC        = 31556926d0                                                      ! [s] conversion variable
    REAL(8),PARAMETER :: UP            = 1.01                                                            ! parameters for setting the length of the time step
    REAL(8),PARAMETER :: DOWN          = 1.0
    REAL(8),PARAMETER :: UPTOL         = 0.04
    REAL(8),PARAMETER :: DOWNTOL       = 0.07                                                            
    INTEGER,PARAMETER :: MAXSTEPS      = 1000000                                                         ! maximum number of steps for the time evolution   
    INTEGER,PARAMETER :: NPOLINT       = 2                                                               ! =2 (linear interpolation)

    ! GLOBAL PARAMETERS OF THE PLANETARY BODY:
    INTEGER,PARAMETER :: N             = 500                                                             ! number of interfaces (i.e., (N-1) layers)
    REAL(8),PARAMETER :: JMAX          = 2.                                                              ! spectral degree of the problem
    REAL(8),PARAMETER :: P             = 86164.1                                                         ! [s] spin period
    REAL(8),PARAMETER :: OMEGA         = (2.*PI)/P                                                       ! [s^-1] angular velocity
    REAL(8),PARAMETER :: RTOP          = 6371d3                                                          ! [m] radius of the top boundary
    REAL(8),PARAMETER :: RBOT          = 3480d3                                                          ! [m] radius of the bottom boundary
    REAL(8),PARAMETER :: THETA0        = -0.35                                                           ! [°] obliquity of the body 
                                                                                                         ! ^NOTE: VALUE FOR TRITON!
    ! ATTRIBUTES OF THE MODEL:
    INTEGER,PARAMETER :: MODEL         = 3                                                               ! mode type 1 (HOMOGENEOUS); 2 (LAYERED), 3 (PREM)
    INTEGER,PARAMETER :: DEFORM        = 1                                                               ! mode type 1 (Heviside at 0); 2 (periodic tides)
    INTEGER,PARAMETER :: TSTEP         = 2                                                               ! mode type 1 (constant step); 2 (gradient)
    LOGICAL,PARAMETER :: NOSLIP        = .false.                                                         ! on/off - no slip on the bottom
    LOGICAL,PARAMETER :: CORE_PRESSURE = .true.                                                          ! on/off - pressure in the liquid core, induced by tides and SG
    LOGICAL,PARAMETER :: SGPOT         = .false.                                                          ! on/off - self-gravity potential (top and bottom topography) 
    LOGICAL,PARAMETER :: SGFORCE       = .false.                                                          ! on/off - self-gravity force (top and bottom topography) 
    LOGICAL,PARAMETER :: SGJUMP        = .false.                                                          ! on/off - self-gravity potential (density jumps)
    LOGICAL,PARAMETER :: SGJUMPF       = .false.                                                          ! on/off - self-gravity force (density jumps)
    LOGICAL,PARAMETER :: TIMEEVOL      = .false.                                                          ! on/off - time evolution
    LOGICAL,PARAMETER :: VISCOEL       = .false.                                                          ! on/off - viscoelasticity 
    LOGICAL,PARAMETER :: GHOMO         = .false.                                                         ! on/off - homogeneous gravitational acceleration
    LOGICAL,PARAMETER :: ETAEXP        = .false.                                                         ! on/off - Arrhenius law (homogeneous model)
    LOGICAL,PARAMETER :: COMPRESS      = .false.                                                          ! on/off - compressible CE & EOM

    ! MODEL 1 (HOMOGENEOUS):
    REAL(8),PARAMETER :: gHmg          = 10d0                                                            ! [m/s^2] homogeneous gravitational acceleration
    REAL(8),PARAMETER :: rhoCore       = 10750d0                                                         ! [kg/m^3] core density
    REAL(8),PARAMETER :: rhoHmg        = 3037d0                                                          ! [kg/m^3] homogeneous model density
    REAL(8),PARAMETER :: muHmg         = 0.70363d11                                                      ! [Pa] homogeneous shear modulus
    REAL(8),PARAMETER :: etaHmg        = 1d21                                                            ! [Pa·s] homogeneous viscosity
    REAL(8),PARAMETER :: lmdHmg        = 1d11                                                            ! [Pa] Lamé elastic parameter

    !LAYERED MODEL (adapted from Spada et al. 2011):
    INTEGER,PARAMETER                  :: noj=3                                                                     ! number of density jumps
    REAL(8),PARAMETER,DIMENSION(noj+2) :: rJumps=(/RBOT,5701d3,5951d3,6301d3,RTOP/)                                 ! radii corresponding to density jumps 
    REAL(8),PARAMETER,DIMENSION(noj+2) :: rhoJumps=(/rhoCore,4978d0,3871d0,3438d0,3037d0/)                          ! densities corresponding to density jumps
    REAL(8),PARAMETER,DIMENSION(noj+2) :: muJumps=(/0d0,2.28340d11,1.05490d11,0.70363d11,0.50605d11/)               ! shear modulus corresponding to jumps
    REAL(8),PARAMETER,DIMENSION(noj+2) :: etaJumps=(/0d0,2d21,1d21,1d21,1d30/)                                      ! viscosity corresponding to jumps
    REAL(8),PARAMETER,DIMENSION(noj+2) :: lmdJumps=(/1d11,1d11,1d11,1d11,1d11/)                                     ! Lamé elastic constant

    ! TIME EVOLUTION:
    CHARACTER(15),PARAMETER :: RULE    = "trapez"                                                        ! = "trapez", = "simpson"
    REAL(8),PARAMETER       :: TMAX    = 1d6                                                             ! [yrs] maximum time of the evolution

    ! ARRHENIUS LAW FOR HOMOGENEOUS MODEL:
    INTEGER,PARAMETER :: TMODE         = 1                                                               ! 1 - conductive profile, 2 - convective profile                                                   
    REAL(8),PARAMETER :: TTOP          = 38d0                                                            ! [K] surface temperature
    REAL(8),PARAMETER :: TBOT          = 273d0                                                           ! [K] bottom temperature   
    REAL(8),PARAMETER :: E             = 59000d0                                                         ! [J/mol]
    REAL(8),PARAMETER :: R             = 8.3144598d0                                                       ! [J/(K·mol)]
    REAL(8),PARAMETER :: ETA_BASE      = 1d13                                                            ! [Pa·s]
    REAL(8),PARAMETER :: ETA_CUT       = 1d20                                                            ! [Pa·s]

    ! LINEAR ALGEBRA AND NUMERICAL RECIPES:
    INTEGER,PARAMETER :: m1            = 7                                                               ! number of subdiagonal rows
    INTEGER,PARAMETER :: m2            = 7                                                               ! number of superdiagonal rows
    INTEGER,PARAMETER :: od            = m1+1                                                            ! original diagonal of the square matrix
    INTEGER,PARAMETER :: m             = m1+m2+1                                                         ! main diagonal & subdiagonals (m1) & superdiagonals (m2)
    INTEGER,PARAMETER :: nn            = 6*N+2                                                           ! size of the matrix
    INTEGER,PARAMETER :: mpl           = m1                                                              ! auxiliary parameter for bandec & banbks
    INTEGER,PARAMETER :: np            = nn                                                              ! auxiliary parameter for bandec & banbks
    INTEGER,PARAMETER :: mp            = m                                                               ! auxiliary parameter for bandec & banbks

    ! MODEL OUTPUT:
    LOGICAL,PARAMETER :: MATRIX_OUTPUT = .true.                                           
    LOGICAL,PARAMETER :: MODEL_OUTPUT  = .true.
    LOGICAL,PARAMETER :: MODEL_REVERSE = .true.

END MODULE mPar
