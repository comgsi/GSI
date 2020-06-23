
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

!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module scatt_table                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE scatt_table

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This module creates a look-up table of values - the portion of the simulated 
! radar variable  calculation the uses scattering amplitudes, slope, and 
! shape parameters - using a range of slope and shape parameter values
! along with the rfopt2 T-matrix method for faster use in data
! assimilation. 
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 2/1/2015
!
!-----------------------------------------------------------------------

  USE rsa_table 

  IMPLICIT NONE
  SAVE

!-----------------------------------------------------------------------
! PARAMETER
! rlss: scattering sum with lamda for rain drop
! slss: scattering sum with lamda for snow aggregate
! hlss: scattering sum with lamda for hailstone
! glss: scattering sum with lamda for graupel
!-----------------------------------------------------------------------

  INTEGER :: n_lam, n_alpha 
  REAL :: intvlrl,intvlra,intvlsl,intvlsa,intvlhl,intvlha,intvlgl,intvlga 
  REAL, ALLOCATABLE :: rlss(:,:,:)
  REAL, ALLOCATABLE :: slss(:,:,:,:)
  REAL, ALLOCATABLE :: hglss(:,:,:,:,:)
  REAL, ALLOCATABLE :: lam_list(:,:)
  REAL, ALLOCATABLE :: alf_list(:,:) 

  REAL, ALLOCATABLE :: scatth_r(:,:), scattv_r(:,:), scatthv_r(:,:), scattk_r(:,:)
  REAL, ALLOCATABLE :: scatth_s(:,:,:),scattv_s(:,:,:),scatthv_s(:,:,:),scattk_s(:,:,:)
  REAL, ALLOCATABLE :: scatth_ha(:,:,:),scattv_ha(:,:,:),scatthv_ha(:,:,:),scattk_ha(:,:,:)
  REAL, ALLOCATABLE :: scatth_g(:,:,:),scattv_g(:,:,:),scatthv_g(:,:,:),scattk_g(:,:,:)

  CONTAINS

  SUBROUTINE calc_scatt_table (on_hail,on_grpl,nlam,nalpha,l_min,l_max,   &
                                alf_min,alf_max)

!-----------------------------------------------------------------------
! Store radar scattering amplitudes table calculated using T-matrix
! method.
!-----------------------------------------------------------------------

    INTEGER, INTENT(IN) :: on_hail,on_grpl 

    INTEGER, INTENT(IN) :: nlam, nalpha 

    REAL, INTENT(IN) :: l_min(4),l_max(4),alf_min(4),alf_max(4)

    INTEGER :: istatus, i, j, k, n_ice

!-----------------------------------------------------------------------
!  Allocate arrays
!-----------------------------------------------------------------------

    n_lam = nlam
    n_alpha = nalpha 

    ALLOCATE(lam_list      (nlam,4),stat=istatus)
    ALLOCATE(alf_list (nalpha,4),stat=istatus)
    ALLOCATE(rlss   (ns,nlam,nalpha),stat=istatus) 
    ALLOCATE(slss(ns,nlam,nalpha,nfw),stat=istatus)

    ALLOCATE(scatth_r(nlam,nalpha),stat=istatus)
    ALLOCATE(scattv_r(nlam,nalpha),stat=istatus)
    ALLOCATE(scatthv_r(nlam,nalpha),stat=istatus)
    ALLOCATE(scattk_r(nlam,nalpha),stat=istatus)

    ALLOCATE(scatth_s(nlam,nalpha,nfw),stat=istatus)
    ALLOCATE(scattv_s(nlam,nalpha,nfw),stat=istatus)
    ALLOCATE(scatthv_s(nlam,nalpha,nfw),stat=istatus)
    ALLOCATE(scattk_s(nlam,nalpha,nfw),stat=istatus)

    IF(on_hail == 1 .AND. on_grpl == 1) THEN
      n_ice = 2
      ALLOCATE(hglss(2,ns,nlam,nalpha,nfw),stat = istatus)
    ELSE
      n_ice = 1
      ALLOCATE(hglss(1,ns,nlam,nalpha,nfw),stat = istatus)
    END IF

    IF(on_hail == 1) THEN
      ALLOCATE(scatth_ha (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scattv_ha (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scatthv_ha (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scattk_ha (nlam,nalpha,nfw),stat = istatus)
    END IF 

    IF(on_grpl == 1) THEN
      ALLOCATE(scatth_g (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scattv_g (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scatthv_g (nlam,nalpha,nfw),stat = istatus)
      ALLOCATE(scattk_g (nlam,nalpha,nfw),stat = istatus)
    END IF 
    
    intvlrl = (l_max(1)-l_min(1))/nlam
    lam_list(1,1) = l_min(1)
    DO i = 2,nlam
       lam_list(i,1) = lam_list(i-1,1)+intvlrl
    END DO 

    IF(nalpha > 1) THEN 
      intvlra = (alf_max(1)-alf_min(1))/nalpha
      alf_list(1,1) = alf_min(1)
      DO i = 2,nalpha
         alf_list(i,1) = alf_list(i-1,1)+intvlra 
      END DO 
    ELSE
      intvlra = 1.0
      alf_list(1,1) = alf_min(1)
    END IF 
        
    intvlsl = (l_max(2)-l_min(2))/nlam
    lam_list(1,2) = l_min(1)
    DO i = 2,nlam
       lam_list(i,2) = lam_list(i-1,2)+intvlsl
    END DO

    IF(nalpha > 1) THEN
      intvlsa = (alf_max(2)-alf_min(2))/nalpha
      alf_list(1,2) = alf_min(2)
      DO i = 2,nalpha
         alf_list(i,2) = alf_list(i-1,2)+intvlsa
      END DO
    ELSE
      intvlsa = 1.0
      alf_list(1,2) = alf_min(2) 
    END IF  

    IF(on_hail == 1) THEN
      intvlhl = (l_max(3)-l_min(3))/nlam
      lam_list(1,3) = l_min(3)
      DO i = 2,nlam
        lam_list(i,3) = lam_list(i-1,3)+intvlhl
      END DO

      IF(nalpha > 1) THEN
        intvlha = (alf_max(3)-alf_min(3))/nalpha
        alf_list(1,3) = alf_min(3)
        DO i = 2,nalpha
          alf_list(i,3) = alf_list(i-1,3)+intvlha
        END DO
      ELSE
        intvlha = 1.0
        alf_list(1,3) = alf_min(1)
      END IF  
    END IF 

    IF(on_grpl == 1) THEN
      intvlgl = (l_max(4)-l_min(4))/nlam
      lam_list(1,4) = l_min(4)
      DO i = 2,nlam
        lam_list(i,4) = lam_list(i-1,4)+intvlgl
      END DO

      IF(nalpha > 1) THEN
        intvlga = (alf_max(4)-alf_min(4))/nalpha
        alf_list(1,4) = alf_min(4)
        DO i = 2,nalpha
          alf_list(i,4) = alf_list(i-1,4)+intvlga
        END DO
      ELSE
        intvlga = 1.0 
        alf_list(1,4) = alf_min(4)
      END IF  
    END IF

    CALL calc_scatt_sum(nlam,nalpha,on_hail,on_grpl)
  
    scatth_r = rlss(1,:,:)
    scattv_r = rlss(2,:,:)
    scatthv_r = rlss(3,:,:)
    scattk_r = rlss(4,:,:)

    scatth_s = slss(1,:,:,:)
    scattv_s = slss(2,:,:,:)
    scatthv_s = slss(3,:,:,:)
    scattk_s = slss(4,:,:,:) 

    IF(on_hail == 1 .AND. on_grpl == 1) THEN
      scatth_ha = hglss(1,1,:,:,:)
      scattv_ha = hglss(1,2,:,:,:)
      scatthv_ha = hglss(1,3,:,:,:)
      scattk_ha = hglss(1,4,:,:,:) 
   
      scatth_g = hglss(2,1,:,:,:)
      scattv_g = hglss(2,2,:,:,:)
      scatthv_g = hglss(2,3,:,:,:)
      scattk_g = hglss(2,4,:,:,:)
    ELSE IF(on_hail == 1) THEN
      scatth_ha = hglss(1,1,:,:,:)
      scattv_ha = hglss(1,2,:,:,:)
      scatthv_ha = hglss(1,3,:,:,:)
      scattk_ha = hglss(1,4,:,:,:)
    ELSE IF(on_grpl == 1) THEN
      scatth_g = hglss(1,1,:,:,:)
      scattv_g = hglss(1,2,:,:,:)
      scatthv_g = hglss(1,3,:,:,:)
      scattk_g = hglss(1,4,:,:,:)
    END IF  


    DEALLOCATE(rlss,slss,hglss)

  END SUBROUTINE calc_scatt_table

  SUBROUTINE calc_scatt_sum(nlam,nalpha,on_hail,on_grpl)


!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Pre-calculate an array of a partial portion of reflectivity and dual-pol
! variables using pre-determined slope and shape parameter values to
! speed up data assimiliation process.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 03/01/2015
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nlam,nalpha
    INTEGER, INTENT(IN) :: on_hail,on_grpl

    INTEGER :: i,j,k,m,hg

    REAL :: dintv

    REAL,PARAMETER,DIMENSION(21) :: fw_list = (/0.0,0.05,0.1,0.15,0.2,0.25, &
                                      0.3, 0.35,0.4,0.45,0.5,0.55,0.6,0.65, &
                                         0.7,0.75,0.8,0.85,0.9,0.95,1.0/)

    REAL, PARAMETER :: pi = 3.141592

    REAL,PARAMETER :: As = 0.8140
    REAL,PARAMETER :: Bs = 0.0303
    REAL,PARAMETER :: Cs = 0.0778
    REAL,PARAMETER :: Ds = 0.4221
    REAL,PARAMETER :: Cks = 0.7837

    REAL,PARAMETER :: sf = 0.8 

    REAL,PARAMETER :: sigmahd = 1.0472
    REAL,PARAMETER :: Ahd = 0.4308
    REAL,PARAMETER :: Bhd = 0.3192
    REAL,PARAMETER :: Chd = 0.1250
    REAL,PARAMETER :: Dhd = 0.3750
    REAL,PARAMETER :: Ckhd = 0.1116

    REAL :: sigmah, Ah, Bh, Ch, Dh, Ckh

    REAL,PARAMETER :: sigmagd = 1.0472
    REAL,PARAMETER :: Agd = 0.4308
    REAL,PARAMETER :: Bgd = 0.3192
    REAL,PARAMETER :: Cgd = 0.1250
    REAL,PARAMETER :: Dgd = 0.3750
    REAL,PARAMETER :: Ckgd = 0.1116

    REAL :: sigmag, Ag, Bg, Cg, Dg, Ckg

    rlss = 0.0
    slss = 0.0
    hglss = 0.0

    dintv = (dsr(2) - dsr(1))*1.e-3
    DO j=1,nalpha
      DO i=1,nlam
        DO k = 1,nd
           rlss(1,i,j) = rlss(1,i,j) + (ABS(far_b(k))**2)* &
                         sngl((dsr(k)*1.e-3)**dble(alf_list(j,1))*   &
                         exp(dble(-lam_list(i,1))*dsr(k)*1.e-3)*dintv)

           rlss(2,i,j) = rlss(2,i,j) + (ABS(fbr_b(k))**2)*&
                         sngl((dsr(k)*1.e-3)**dble(alf_list(j,1))*  &
                         exp(dble(-lam_list(i,1))*dsr(k)*1.e-3)*dintv)

           rlss(3,i,j) = rlss(3,i,j) + far_b(k)*CONJG(fbr_b(k))*&
                         sngl((dsr(k)*1.e-3)**dble(alf_list(j,1))* &
                         exp(dble(-lam_list(i,1))*dsr(k)*1.e-3)*dintv)

           rlss(4,i,j) = rlss(4,i,j) + REAL(far_f(k) - fbr_f(k))*&
                         sngl((dsr(k)*1.e-3)**dble(alf_list(j,1))*  &
                         exp(dble(-lam_list(i,1))*dsr(k)*1.e-3)*dintv)
        END DO
      END DO
    END DO

    dintv = (dss(2) - dss(1))*1.e-3

    DO m = 1,nfw
      DO j=1,nalpha
        DO i=1,nlam
          DO k = 1,nd

            slss(1,i,j,m) = slss(1,i,j,m) + (As*ABS(fas_b(k,m))**2 + Bs*ABS(fbs_b(k,m))**2 +  &
                         2*Cs*Real(fas_b(k,m)*CONJG(fbs_b(k,m)))) * &
   sngl((dss(k)*1.e-3)**dble(alf_list(j,2))*exp(dble(-lam_list(i,2))*dss(k)*1.e-3)*dintv)

            slss(2,i,j,m) = slss(2,i,j,m) + (Bs*ABS(fas_b(k,m))**2 + As*ABS(fbs_b(k,m))**2 + &
                         2*Cs*REAL(fas_b(k,m)*CONJG(fbs_b(k,m))))* &
  sngl((dss(k)*1.e-3)**dble(alf_list(j,2))*exp(dble(-lam_list(i,2))*dss(k)*1.e-3)*dintv)

            slss(3,i,j,m) = slss(3,i,j,m) + (Cs*(ABS(fas_b(k,m))**2 + ABS(fbs_b(k,m))**2) + &
                   As*fas_b(k,m)*CONJG(fbs_b(k,m)) + Bs*fbs_b(k,m)*CONJG(fas_b(k,m)))* &
  sngl((dss(k)*1.e-3)**dble(alf_list(j,2))*exp(dble(-lam_list(i,2))*dss(k)*1.e-3)*dintv)

            slss(4,i,j,m) = slss(4,i,j,m) + Cks*REAL(fas_f(k,m) - fbs_f(k,m))*&
  sngl((dss(k)*1.e-3)**dble(alf_list(j,2))*exp(dble(-lam_list(i,2))*dss(k)*1.e-3)*dintv)
          END DO
        END DO
      END DO
    END DO

    IF(on_hail == 1) THEN
      dintv = (dsh(2) - dsh(1))*1.e-3
      DO m=1,nfw
        IF(m == 1) THEN
          Ah = Ahd
          Bh = Bhd
          Ch = Chd
          Dh = Dhd
          Ckh = Ckhd
        ELSE 
          sigmah=pi/3*(1-sf*fw_list(m))
          Ah=.125*(3+4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
          Bh=.125*(3-4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
          Ch=.125*(1-exp(-8*sigmah**2))
          Dh=.125*(3+exp(-8*sigmah**2))
          Ckh=exp(-2*sigmah**2)
        END IF 

        DO j=1,nalpha
          DO i = 1,nlam   
            DO k = 1,nd

              hglss(1,1,i,j,m) = hglss(1,1,i,j,m) + (Ah*ABS(fah_b(k,m))**2 + Bh*ABS(fbh_b(k,m))**2 +  &
                         2*Ch*Real(fah_b(k,m)*CONJG(fbh_b(k,m)))) * &
                         sngl((dsh(k)*1.e-3)**dble(alf_list(j,3))*exp(-dble(lam_list(i,3))*dsh(k)*1.e-3)*dintv)

              hglss(1,2,i,j,m) = hglss(1,2,i,j,m) + (Bh*ABS(fah_b(k,m))**2 + Ah*ABS(fbh_b(k,m))**2 + &
                         2*Ch*REAL(fah_b(k,m)*CONJG(fbh_b(k,m))))* &
                         sngl((dsh(k)*1.e-3)**dble(alf_list(j,3))*exp(dble(-lam_list(i,3))*dsh(k)*1.e-3)*dintv)

              hglss(1,3,i,j,m) = hglss(1,3,i,j,m) + (Ch*(ABS(fah_b(k,m))**2 + ABS(fbh_b(k,m))**2) + &
                 Ah*fah_b(k,m)*CONJG(fbh_b(k,m)) + Bh*fbh_b(k,m)*CONJG(fah_b(k,m)))* &
                          sngl((dsh(k)*1.e-3)**dble(alf_list(j,3))*exp(dble(-lam_list(i,3))*dsh(k)*1.e-3)*dintv)


              hglss(1,4,i,j,m) = hglss(1,4,i,j,m) + Ckh*REAL(fah_f(k,m) - fbh_f(k,m))*  &
                         sngl((dsh(k)*1.e-3)**dble(alf_list(j,3))*exp(dble(-lam_list(i,3))*dsh(k)*1.e-3)*dintv)

            END DO
          END DO
        END DO
      END DO
    ENDIF

    IF(on_grpl == 1) THEN

     IF(on_hail == 1) THEN
       hg = 2
     ELSE
       hg = 1 
     END IF

     dintv = (dsg(2) - dsg(1))*1.e-3
     DO m=1,nfw
       IF(m == 1) THEN
         Ag = Agd
         Bg = Bgd
         Cg = Cgd
         Dg = Dgd 
         Ckg = Ckgd 
       ELSE     
         sigmag=pi/3*(1-sf*fw_list(m))
         Ag=.125*(3+4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
         Bg=.125*(3-4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
         Cg=.125*(1-exp(-8*sigmag**2))
         Dg=.125*(3+exp(-8*sigmag**2))
         Ckg=exp(-2*sigmag**2)
       END IF 

       DO j=1,nalpha
         DO i=1,nlam
           DO k=1,nd
             hglss(hg,1,i,j,m) = hglss(hg,1,i,j,m) + (Ag*ABS(fag_b(k,m))**2 + Bg*ABS(fbg_b(k,m))**2 +  &
                         2*Cg*Real(fag_b(k,m)*CONJG(fbg_b(k,m)))) * &
                         sngl((dsg(k)*1.e-3)**dble(alf_list(j,4))*exp(dble(-lam_list(i,4))*dsg(k)*1.e-3)*dintv)

             hglss(hg,2,i,j,m) = hglss(hg,2,i,j,m) + (Bg*ABS(fag_b(k,m))**2 + Ag*ABS(fbg_b(k,m))**2 + &
                         2*Cg*REAL(fag_b(k,m)*CONJG(fbg_b(k,m))))* &
                         sngl((dsg(k)*1.e-3)**dble(alf_list(j,4))*exp(dble(-lam_list(i,4))*dsg(k)*1.e-3)*dintv)

             hglss(hg,3,i,j,m) = hglss(hg,3,i,j,m) + (Cg*(ABS(fag_b(k,m))**2 + ABS(fbg_b(k,m))**2) + &
                 Ag*fag_b(k,m)*CONJG(fbg_b(k,m)) + Bg*fbg_b(k,m)*CONJG(fag_b(k,m)))* &
                          sngl((dsg(k)*1.e-3)**dble(alf_list(j,4))* exp(dble(-lam_list(i,4))*dsg(k)*1.e-3)*dintv)

             hglss(hg,4,i,j,m) = hglss(hg,4,i,j,m) + Ckg*REAL(fag_f(k,m) - fbg_f(k,m))* &
                        sngl((dsg(k)*1.e-3)**dble(alf_list(j,4))*exp(dble(-lam_list(i,4))*dsg(k)*1.e-3)*dintv)
            END DO
          END DO
        END DO
      END DO
    ENDIF

  END SUBROUTINE calc_scatt_sum

END MODULE scatt_table

