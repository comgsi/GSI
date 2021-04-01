MODULE radaremul_iface
! initialize variables required for radaremul library
! required for direct reflectivity DA capabilities
CONTAINS
SUBROUTINE init_mphyopt(mype,iret)

! use mpeu_util,only: die
! use mpimod, only: mype

  use radaremul_cst
  use dualpolepara, only: qgh_opt, grpl_ON, hl_ON, lambda, calcConstants

  IMPLICIT NONE

  integer, intent(in   ) :: mype
  integer, intent(inout) :: iret

  INTEGER, EXTERNAL :: get_qgh_opt

  IF (mype == 0) &
    WRITE(6,*) "INIT_MPHYOPT: Initializing radar emulator ...  (mphyopt=",mphyopt,")"

  P_QC = 0; P_QR = 0; P_QI = 0; P_QS = 0; P_QH = 0; P_QG = 0
  P_NC = 0; P_NR = 0; P_NI = 0; P_NS = 0; P_NH = 0; P_NG = 0
            P_ZR = 0; P_ZI = 0; P_ZS = 0; P_ZH = 0; P_ZG = 0

  alpharain = 0.0
  alphasnow = 0.0
  alphagrpl = 0.0
  alphahail = 0.0
  alphaice = 0.0

  IF (mphyopt == 1) THEN
    nscalar = 2
    P_QC = 1; P_QR = 2

    graupel_ON = 0
    hail_ON = 1
  ELSE IF ( mphyopt == 2 .OR. mphyopt == 3 .OR. mphyopt == 4 ) THEN
    nscalar = 5
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QH = 5

    graupel_ON = 0
    hail_ON = 1
  ELSE IF ( mphyopt == 5 .OR. mphyopt == 6 .OR. mphyopt == 7 ) THEN
    nscalar = 5
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5

    graupel_ON = 1
    hail_ON = 0
  ELSE IF ( mphyopt == 8 ) THEN
    nscalar = 6
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5; P_QH = 6;

    graupel_ON = 1
    hail_ON = 1
  ELSE IF ( mphyopt == 9 .OR. mphyopt == 10 .OR. mphyopt == 12) THEN
    nscalar = 12
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS =  4; P_QG =  5; P_QH =  6;
    P_NC = 7; P_NR = 8; P_NI = 9; P_NS = 10; P_NG = 11; P_NH = 12;

    graupel_ON = 1
    hail_ON = 1
  ELSE IF ( mphyopt == 11) THEN
    nscalar = 17
    P_QC = 1; P_QR =  2; P_QI =  3; P_QS =  4; P_QG =  5; P_QH =  6;
    P_NC = 7; P_NR =  8; P_NI =  9; P_NS = 10; P_NG = 11; P_NH = 12;
              P_ZR = 13; P_ZI = 14; P_ZS = 15; P_ZG = 16; P_ZH = 17;

    graupel_ON = 1
    hail_ON = 1
  ELSE IF (mphyopt == 100) THEN      ! passiveqv
    nscalar = 0

    graupel_ON = 0
    hail_ON = 0
  ELSE IF (mphyopt == 101 .OR. mphyopt == 103) THEN  ! kesslerscheme, wsm3scheme
    P_QC = 1; P_QR = 2
    nscalar = 2

    graupel_ON = 0
    hail_ON = 0
  ELSE IF (mphyopt == 102 .OR. mphyopt == 106) THEN  ! linscheme, wsm6scheme
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
    nscalar  = 5

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mphyopt == 104 .OR. mphyopt == 199) THEN  ! wsm5scheme,ncepcloud5
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4
    nscalar  = 4

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mphyopt == 105) THEN                      ! etampnew
    P_QC = 1; P_QR = 2; P_QS = 3
    nscalar  = 3

    graupel_ON = 0
    hail_ON = 0
  ELSE IF (mphyopt == 108) THEN                      ! thompson
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
              P_NR = 6; P_NI = 7
    nscalar = 7

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mphyopt == 110) THEN                      ! Morrison DB
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
              P_NR = 6; P_NI = 7; P_NS = 8; P_NG = 9
    nscalar  = 9

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mphyopt == 114) THEN                      ! DM-5
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4
    P_NC = 5; P_NR = 6
    nscalar  = 6

    graupel_ON = 0
    hail_ON = 0
  ELSE IF (mphyopt == 116) THEN                      ! DM-6
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
    P_NC = 6; P_NR = 7
    nscalar  = 7

    graupel_ON = 1
    hail_ON = 0
  ELSE IF (mphyopt == 198) THEN                      ! Old Thompson (07)
    P_QC = 1; P_QR = 2; P_QI = 3; P_QS = 4; P_QG = 5
    P_NI = 6
    nscalar  = 6

    graupel_ON = 1
    hail_ON = 0
  END IF

  grpl_ON = graupel_ON
  hl_ON = hail_ON
  qgh_opt = get_qgh_opt(graupel_ON, hail_ON)

! JP: assign default values for radaremul library 
! (used in TM forward operator for EnKF) !
  IF (mphyopt == 108 ) then
    n0rain=8.0e6
    n0snow=3.0e6
    n0hail=4.0e4
    n0grpl=4.0e5
    rhosnow=100.0
    rhohail=913.0
    rhogrpl=400.0
    rfopt=1
    wavelen=107.0
  END IF

  lambda = wavelen
  call calcConstants()

! iret = 0
! IF (mphyopt /= 2 .AND. mphyopt /= 108) THEN
!!      CALL die('init_mphyopt', 'Invalid microphysics option', mphyopt) 
!     iret=-1
! END IF

  iret = -1
  SELECT CASE (mphyopt)
  CASE(2,3,4)             ! Lin
      iret = 0
  CASE(5,6,7)             ! WSM6
      iret = 0
  CASE(108)               ! Thompson
      iret = 0
  CASE DEFAULT            ! not ready for dbz operator
      iret = -1
  END SELECT



END SUBROUTINE init_mphyopt

END MODULE
