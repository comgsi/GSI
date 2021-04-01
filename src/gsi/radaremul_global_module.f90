! Adopted from RADREMUL library of ARPS system
! Used for direct reflectity DA capability
!########################################################################
!########################################################################
!#########                                                      #########
!#########                 module radaremul_cst                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE radaremul_cst
  IMPLICIT NONE

  INTEGER :: mphyopt
  INTEGER :: hail_ON, graupel_ON
  INTEGER :: nscalar
  INTEGER :: MFflg = 2

  INTEGER :: P_QC, P_QR, P_QI, P_QS, P_QH, P_QG
  INTEGER :: P_NC, P_NR, P_NI, P_NS, P_NH, P_NG
  INTEGER ::       P_ZR, P_ZI, P_ZS, P_ZH, P_ZG

  REAL :: n0rain, n0snow, n0hail, n0grpl
  REAL :: rhosnow, rhohail, rhogrpl
  REAL :: alpharain,alphasnow,alphagrpl,alphahail,alphaice
  REAL :: alpha_dsd(6)

  INTEGER :: dsdparaopt = 0

  INTEGER :: nen = 1

  INTEGER :: rfopt
  REAL :: wavelen

END MODULE radaremul_cst


!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module global_paraest               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE global_paraest

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declare variables to be used in DSD parameter retrieval
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/4/2007
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! Declare variable arrays
! admisbl_ub  Upper bound of admissible range of DSD parameters (in log)
! admisbl_lb  Lower bound of admissible range of DSD parameters (in log)
!-----------------------------------------------------------------------
  REAL, allocatable :: para(:)
  REAL, allocatable :: paraen(:,:)
  REAL, allocatable :: paraprim(:,:)
  REAL, allocatable :: paraspd(:)
  REAL, DIMENSION(5) :: admisbl_ub = (/79.03,80.00,66.02,26.02,29.60/)
  REAL, DIMENSION(5) :: admisbl_lb = (/64.77,56.99,26.02,13.01,26.02/)

  INTEGER, ALLOCATABLE, TARGET :: corflgrdr(:,:,:,:,:)
  INTEGER, POINTER     :: corflg(:,:,:,:)
  REAL,    ALLOCATABLE :: corprdr(:,:,:,:,:)
  INTEGER, ALLOCATABLE :: histogz(:,:)
  REAL,    ALLOCATABLE :: rmscorprdr(:,:)

  CONTAINS

  SUBROUTINE allocate_paraestArray(nen,paranum)

    INTEGER :: nen,paranum,istatus

    ALLOCATE(para(paranum),stat=istatus)
    ALLOCATE(paraen(nen,paranum),stat=istatus)

    para   = 0.
    paraen = 0.

  END SUBROUTINE allocate_paraestArray

  SUBROUTINE allocate_paraprimArray(nen,paranum)

    INTEGER :: nen,paranum,istatus

    ALLOCATE(paraprim(nen,paranum),stat=istatus)
    ALLOCATE(paraspd(paranum),stat=istatus)

    paraprim = 0.
    paraspd = 0.

  END SUBROUTINE allocate_paraprimArray

  SUBROUTINE deallocate_paraestArray()
    DEALLOCATE(para,paraen,paraprim,paraspd)
  END SUBROUTINE deallocate_paraestArray

  SUBROUTINE allocate_covArray(nx,ny,nv,paranum,numvar,nrank)

    INTEGER :: nx,ny,nv,paranum,nrank,numvar,istatus

    ALLOCATE(corprdr(nx,ny,nv,paranum,numvar),stat=istatus)
    ALLOCATE(histogz(nrank,paranum),stat=istatus)
    ALLOCATE(rmscorprdr(paranum,numvar),stat=istatus)

    corprdr = 0.0
    histogz = 0
    rmscorprdr = 0.0

  END SUBROUTINE allocate_covArray

END MODULE global_paraest


!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module mod_reflec_ferrier           #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE mod_reflec_ferrier

!-----------------------------------------------------------------------
!
! PURPOSE:
!    Arrays to store precalculated exponentail values to be used in
!    the calculation of radar reflectivity factors.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 4/17/2008
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! Variable declaration
!-----------------------------------------------------------------------

  REAL, ALLOCATABLE :: pwrrqrpowf(:)
  REAL, ALLOCATABLE :: pwrrqsnpowf(:)
  REAL, ALLOCATABLE :: pwrrqsppowf(:)
  REAL, ALLOCATABLE :: pwrrqhnpowf(:)
  REAL, ALLOCATABLE :: pwrrqhppowf(:)

END MODULE mod_reflec_ferrier

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module rsa_table                    #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE rsa_table

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine read in the precalculated scattering amplitude
! from the tables, which are calcualted using T-matrix method.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 4/17/2008
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! PARAMETER
! dsr : rain drop size,  rsa: scattering amplitude for rain drop
! dss : snow drop size,  ssa: scattering amplitude for snow aggregate
! dsh : hail drop size,  hsa: scattering amplitude for hailstone
! dsg : grpl drop size,  gsa: scattering amplitude for graupel
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: nd = 100, ns = 8, nfw = 21
  REAL, PRIVATE, ALLOCATABLE :: rsa(:,:)
  REAL, PRIVATE, ALLOCATABLE :: ssa(:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: hsa(:,:,:)
  REAL, PRIVATE, ALLOCATABLE :: gsa(:,:,:)
  REAL, ALLOCATABLE :: dsr(:), dss(:), dsh(:), dsg(:)

  COMPLEX :: far_b(nd), fbr_b(nd), far_f(nd), fbr_f(nd)
  COMPLEX :: fas_b(nd,nfw), fbs_b(nd,nfw), fas_f(nd,nfw), fbs_f(nd,nfw)
  COMPLEX :: fah_b(nd,nfw), fbh_b(nd,nfw), fah_f(nd,nfw), fbh_f(nd,nfw)
  COMPLEX :: fag_b(nd,nfw), fbg_b(nd,nfw), fag_f(nd,nfw), fbg_f(nd,nfw)

  CONTAINS

  SUBROUTINE read_table (rsafndir)
!-----------------------------------------------------------------------
! Store radar scattering amplitudes table calculated using T-matrix
! method.
!-----------------------------------------------------------------------

    INTEGER :: istatus, i, j, k
    CHARACTER (LEN=256) :: rsafndir
    CHARACTER (LEN=256) :: rsafn
    CHARACTER (LEN=3), DIMENSION(nfw) :: extn = (/'000','005','010',    &
           '015','020','025','030','035','040','045','050','055','060', &
           '065','070','075','080','085','090','095','100'/)
    CHARACTER (LEN=256) :: head

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------
    ALLOCATE(dsr      (nd),stat=istatus)
    ALLOCATE(dss      (nd),stat=istatus)
    ALLOCATE(dsh      (nd),stat=istatus)
    ALLOCATE(dsg      (nd),stat=istatus)
    ALLOCATE(rsa   (ns,nd),stat=istatus)
    ALLOCATE(ssa(ns,nd,nfw),stat=istatus)
    ALLOCATE(hsa(ns,nd,nfw),stat=istatus)
    ALLOCATE(gsa(ns,nd,nfw),stat=istatus)

!-----------------------------------------------------------------------
!  Read rain
!-----------------------------------------------------------------------
    rsafn = TRIM(rsafndir)//'/SCTT_RAIN_fw100.dat'
    OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
    READ(51,*) head

    DO j=1,100
      READ(51,'(f5.2,8e13.5)') dsr(j), (rsa(i,j),i=1,8)
    ENDDO
    CLOSE(51)

    far_b = CMPLX(rsa(1,:),rsa(2,:))
    fbr_b = CMPLX(rsa(3,:),rsa(4,:))
    far_f = CMPLX(rsa(5,:),rsa(6,:))
    fbr_f = CMPLX(rsa(7,:),rsa(8,:))

!-----------------------------------------------------------------------
!  Read snow
!-----------------------------------------------------------------------
    DO k=1,nfw
      rsafn = TRIM(rsafndir)//'/SCTT_SNOW_fw'//extn(k)//'.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,100
        READ(51,'(f5.2,8e13.5)') dss(j), (ssa(i,j,k),i=1,8)
      ENDDO
      CLOSE(51)
    ENDDO

    fas_b = CMPLX(ssa(1,:,:),ssa(2,:,:))
    fbs_b = CMPLX(ssa(3,:,:),ssa(4,:,:))
    fas_f = CMPLX(ssa(5,:,:),ssa(6,:,:))
    fbs_f = CMPLX(ssa(7,:,:),ssa(8,:,:))

!-----------------------------------------------------------------------
!  Read hail
!-----------------------------------------------------------------------
    DO k=1, nfw
      rsafn = TRIM(rsafndir)//'/SCTT_HAIL_fw'//extn(k)//'.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,100
        READ(51,'(f5.2,8e13.5)') dsh(j), (hsa(i,j,k),i=1,8)
      ENDDO
      CLOSE(51)
    ENDDO

    fah_b = CMPLX(hsa(1,:,:),hsa(2,:,:))
    fbh_b = CMPLX(hsa(3,:,:),hsa(4,:,:))
    fah_f = CMPLX(hsa(5,:,:),hsa(6,:,:))
    fbh_f = CMPLX(hsa(7,:,:),hsa(8,:,:))

!-----------------------------------------------------------------------
!  Read graupel
!-----------------------------------------------------------------------
    DO k=1, nfw
      rsafn = TRIM(rsafndir)//'/SCTT_GRPL_fw'//extn(k)//'.dat'
      OPEN(UNIT=51,FILE=TRIM(rsafn),STATUS='old',FORM='formatted')
      READ(51,*) head

      DO j=1,100
        READ(51,'(f5.2,8e13.5)') dsg(j), (gsa(i,j,k),i=1,8)
      ENDDO
      CLOSE(51)
    ENDDO

    fag_b = CMPLX(gsa(1,:,:),gsa(2,:,:))
    fbg_b = CMPLX(gsa(3,:,:),gsa(4,:,:))
    fag_f = CMPLX(gsa(5,:,:),gsa(6,:,:))
    fbg_f = CMPLX(gsa(7,:,:),gsa(8,:,:))

    deallocate(rsa,ssa,hsa,gsa)

  END SUBROUTINE read_table

END MODULE rsa_table

