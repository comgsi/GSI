!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module dualpara                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE DUALPARA

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declare some constants used for calculaion of dual polarization
! parameters such as Zhh, Zdr, and Kdp. (It can be expanded...)
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
!-----------------------------------------------------------------------
! Declare parameters.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

  REAL, PARAMETER :: pi = 3.141592   ! pi

  REAL :: lambda           ! wavelength of radar (mm)

  REAL,PARAMETER :: Kw2 = 0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: alphaa = 4.28e-4   ! backscattering amplitude constant
                                       ! along major axis for rain
  REAL,PARAMETER :: beta_ra = 3.04
  REAL,PARAMETER :: alphab = 4.28e-4   ! backscattering amplitude constant
                                       ! along minor axis for rain
  REAL,PARAMETER :: beta_rb = 2.77
  REAL,PARAMETER :: alphak = 3.88e-4   ! differential forward scattering
                                       ! amplitude for rain
  REAL,PARAMETER :: alphask = 8.53e-7   ! differential forward scattering
                                        ! amplitude for snow
  REAL,PARAMETER :: alphaa_ds = 1.94e-5 ! for dry snow at horz plane
  REAL,PARAMETER :: alphab_ds = 1.91e-5 ! for dry snow at vert plane
  REAL,PARAMETER :: alphaa_tom_ds = 2.8e-5 !for dry snow at horz plane for Thomposon scheme
  REAL,PARAMETER :: alphab_tom_ds = 2.6e-5 !for dry snow at vert plane for the Thompson scheme 

  REAL, PARAMETER :: beta_sa = 3.0
  REAL, PARAMETER :: beta_sb = 3.0
  REAL, PARAMETER :: beta_tom_dsa = 1.95 ! Special expon for Thompson scheme
  REAL, PARAMETER :: beta_tom_dsb = 1.965! Special expon for Thompson scheme

  REAL,PARAMETER :: alphaa_dh = 1.91e-4 ! for dry hail at horz plane
  REAL,PARAMETER :: alphab_dh = 1.65e-4 ! for dry hail at vert plane

  REAL,PARAMETER :: beta_ha = 3.0
  REAL,PARAMETER :: beta_hb = 3.0

  REAL,PARAMETER :: alphaa_dg = 0.81e-4 ! for dry graupel at horz plane
  REAL,PARAMETER :: alphab_dg = 0.76e-4 ! for dry graupel at vert plane

  REAL,PARAMETER :: beta_ga = 3.0
  REAL,PARAMETER :: beta_gb = 3.0

  REAL,PARAMETER :: alphak_ds = 0.03e-5 ! alphaa_ds - alphab_ds
  REAL,PARAMETER :: alphak_tom_ds = 1.05e-6 !alphaa_ds - alphab_ds for Thompson scheme
  REAL,PARAMETER :: alphak_dh = 0.26e-4 ! alphaa_dh - alphab_dh
  REAL,PARAMETER :: alphak_dg = 0.05e-4 ! alphaa_dh - alphab_dh
  REAL,PARAMETER :: betak_s = 3.0
  REAL,PARAMETER :: betak_tom_ds = 2.04 !For Thompson Scheme 
  REAL,PARAMETER :: betak_h = 3.0
  REAL,PARAMETER :: betak_g = 3.0

  REAL,PARAMETER :: rho_0r = 1.0      ! rho_0 for rain
  REAL,PARAMETER :: rho_0s = 1.0      ! rho_0 for snow
  REAL,PARAMETER :: rho_0h = 0.97     ! rho_0 for hail
  REAL,PARAMETER :: rho_0g = 0.95     ! rho_0 for hail
  REAL,PARAMETER :: rho_0rsi = 0.82   ! lower limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rsf = 0.95   ! upper limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rhi = 0.85   ! lower limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rhf = 0.95   ! upper limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rgi = 0.82   ! lower limit of rho_0rg (rain-graupel mixture)
  REAL,PARAMETER :: rho_0rgf = 0.95   ! upper limit of rho_0rg (rain-graupel mixture)

  REAL,PARAMETER :: degKtoC=273.15 ! Conversion factor from degrees K to
                                   !   degrees C

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)

  REAL,PARAMETER :: mm3todBZ=1.0E+9 ! Conversion factor from mm**3 to
                                    !   mm**6 m**-3.
 
  REAL,PARAMETER :: thom_lam0 = 20.78
  REAL,PARAMETER :: thom_lam1 = 3.29
  REAL,PARAMETER :: thom_k0 = 490.6
  REAL,PARAMETER :: thom_k1 = 17.46 

  REAL,PARAMETER :: unit_factor = 1.e-2  ! Unit conversion factor not addressed
                                         ! in the T-matrix scattering amplitude (size D is in cm in T-matrix)

  REAL :: kdpCoefIce 

  REAL :: c_x(6)  !(PI/6)*rho_qx

  REAL :: ta = 273.16  

  REAL,PARAMETER :: missing = -9999.0
  REAL :: grpl_miss
  REAL :: hl_miss 
  
  LOGICAL :: firstcall = .true.
  INTEGER :: grpl_ON
  INTEGER :: hl_ON 
  INTEGER :: qgh_opt

  INTEGER :: attn_ON

  INTEGER :: dualpol_opt

!-----------------------------------------------------------------------
! Precalculated complete gamma function values
!-----------------------------------------------------------------------
  REAL,PARAMETER :: gamma7_08 = 836.7818
  REAL,PARAMETER :: gamma6_81 = 505.8403
  REAL,PARAMETER :: gamma6_54 = 309.3308
  REAL,PARAMETER :: gamma5_63 = 64.6460
  REAL,PARAMETER :: gamma4_16 = 7.3619
  REAL,PARAMETER :: gamma3_97 = 5.7788

!-----------------------------------------------------------------------
! Variables to can be changed by parameter retrieval
!-----------------------------------------------------------------------
  REAL :: N0r        ! Intercept parameter in 1/(m^4) for rain
  REAL :: N0s        ! Intercept parameter in 1/(m^4) for snow
  REAL :: N0h        ! Intercept parameter in 1/(m^4) for hail
  REAL :: N0g        ! Intercept parameter in 1/(m^4) for hail
  REAL :: N0s2       ! Second intercept parameter in 1/(m^4) for snow

  REAL :: N0ms       !Intercept parameter for melting species 
  REAL :: N0ms2 
  REAL :: N0mh
  REAL :: N0mg 

  REAL :: rhor=1000. ! Density of rain (kg m**-3)
  REAL :: rhoh       ! Density of hail (kg m**-3)
  REAL :: rhos       ! Density of snow (kg m**-3)
  REAL :: rhog       ! Density of graupel (kg m**-3)

  REAL :: alphar     !Shape parameter for rain
  REAL :: alphas     !Shape parameter for snow
  REAL :: alphah     !SHape parameter for hail
  REAL :: alphag     !SHape parameter for graupel
  REAL :: alphas2

  REAL :: lamdar     !slope parameter for rain (1/m)
  REAL :: lamdas     
  REAL :: lamdas2  
  REAL :: lamdams
  REAL :: lamdams2  
  REAL :: lamdag
  REAL :: lamdamg
  REAL :: lamdah
  REAL :: lamdamh 

!-----------------------------------------------------------------------
! Variables to can be changed for meling ice
!-----------------------------------------------------------------------
  REAL :: fos        ! Maximum fraction of rain-snow mixture
  REAL :: foh        ! Maximum fraction of rain-hail mixture
  REAL :: fog        ! Maximum fraction of rain-hail mixture

!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------

 
  REAL :: radar_const    !(4*lambda**4)/(pi*kw2)
  
  REAL :: constKdpr

!-----------------------------------------------------------------------
! Scattering matrix coefficient for snow
!
! phi=0.       (Mean orientation)
! sigmas=pi/9
! As=1/8*(3+4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Bs=1/8*(3-4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Cs=1/8*(1-cos(4*phi)*exp(-8*sigmas**2))
! Ds=1/8*(3+cos(4*phi)*exp(-8*sigmas**2))
! Cks=cos(2*phi)*exp(-2*sigmas**2)
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmas = 0.3491
  REAL,PARAMETER :: As = 0.8140
  REAL,PARAMETER :: Bs = 0.0303
  REAL,PARAMETER :: Cs = 0.0778
  REAL,PARAMETER :: Ds = 0.4221
  REAL,PARAMETER :: Cks = 0.7837

!-----------------------------------------------------------------------
! Scattering matrix coefficient for hail
!
! phi=0.     (Mean orientation)
! sigmah=pi/3*(1-sf*fw)
! Ah=1/8*(3+4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Bh=1/8*(3-4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Ch=1/8*(1-cos(4*phi)*exp(-8*sigmah**2))
! Dh=1/8*(3+cos(4*phi)*exp(-8*sigmah**2))
! Ckh=cos(2*phi)*exp(-2*sigmah**2)
!
! corresponding coefficient for dry hail: Ahd, Bhd, Chd, Dhd, Ckhd
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmahd = 1.0472
  REAL,PARAMETER :: Ahd = 0.4308
  REAL,PARAMETER :: Bhd = 0.3192
  REAL,PARAMETER :: Chd = 0.1250
  REAL,PARAMETER :: Dhd = 0.3750
  REAL,PARAMETER :: Ckhd = 0.1116

  REAL,PARAMETER :: q_threshold = 2.e-4
  REAL :: sf
  REAL :: sigmah, Ah, Bh, Ch, Dh, Ckh

!-----------------------------------------------------------------------
! Scattering matrix coefficient for graupel
! 
! phi=0.     (Mean orientation)
! sigmag=pi/3*(1-sf*fw)
! Ag=1/8*(3+4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Bg=1/8*(3-4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Cg=1/8*(1-cos(4*phi)*exp(-8*sigmag**2))
! Dg=1/8*(3+cos(4*phi)*exp(-8*sigmag**2))
! Ckg=cos(2*phi)*exp(-2*sigmag**2)
! 
! corresponding coefficient for dry graupel: Agd, Bgd, Cgd, Dgd, Ckgd
!-----------------------------------------------------------------------
  
  REAL,PARAMETER :: sigmagd = 1.0472
  REAL,PARAMETER :: Agd = 0.4308
  REAL,PARAMETER :: Bgd = 0.3192
  REAL,PARAMETER :: Cgd = 0.1250
  REAL,PARAMETER :: Dgd = 0.3750
  REAL,PARAMETER :: Ckgd = 0.1116
  
  REAL :: sigmag, Ag, Bg, Cg, Dg, Ckg

!-----------------------------------------------------------------------
!  Declare new observation type
!-----------------------------------------------------------------------

  TYPE T_obs_dual
       REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
       REAL :: T_log_zdr, T_sum_ref_hv, T_kdp
       REAL :: T_Ahh,     T_Avv
       REAL :: T_ref_r_h, T_ref_s_h, T_ref_h_h,T_ref_g_h
       REAL :: T_ref_rs_h,T_ref_rh_h,T_ref_rg_h
       REAL :: T_ref_r_v, T_ref_s_v, T_ref_h_v, T_ref_g_v
       REAL :: T_ref_rs_v, T_ref_rh_v, T_ref_rg_v
  END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

  TYPE T_para_dsd
    REAL :: T_qr, T_qs, T_qh, T_qg
    REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
    REAL :: T_alfr,T_alfs,T_alfh,T_alfg
  END TYPE T_para_dsd

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SUBROUTINES AND FUNCTIONS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  CONTAINS


  SUBROUTINE calcMDR()
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates mass-diameter relation based on MP scheme. 
!
!-----------------------------------------------------------------------
!
! Author:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE radaremul_cst, only: mphyopt

  IMPLICIT NONE

  c_x(1) = (pi/6.)*rhor
  c_x(2) = (pi/6.)*rhor
  c_x(3) = 440.0

  SELECT CASE (mphyopt)
  CASE(1:12,106,109:110,116)
    c_x(4) = (pi/6.)*rhos
  CASE(108)
    c_x(4) = .069 
  END SELECT 

  c_x(5) = (pi/6.)*rhog
  c_x(6) = (pi/6.)*rhoh

  END SUBROUTINE calcMDR

  SUBROUTINE calcConstants()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Precalculate commonly unsed constants to save computations.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/28/2005
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Constant in front of dual pol calculations (4*lambda**4)/(pi*kw2)
!-----------------------------------------------------------------------

   radar_const = (4. * lambda**4.)/(pi**4 * Kw2)

!-----------------------------------------------------------------------
! For Kdp constants
!-----------------------------------------------------------------------

    constKdpr = 180. * lambda  * alphak * 1.0e6 / pi !rain
    kdpCoefIce = (180*lambda*1.e6)/pi !ice 

  END SUBROUTINE calcConstants

  SUBROUTINE init_fox()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can vary depend on whether graupel/hail exists. 
!-----------------------------------------------------------------------
  fos = 0.3             ! Maximum fraction of rain-snow mixture
  foh = 0.2              ! Maximum fraction of rain-hail mixture
  fog = 0.25             ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox

  SUBROUTINE init_fox_no_grpl()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice 
! when graupel is suppressed.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  fos = 0.5              ! Maximum fraction of rain-snow mixture
  foh = 0.3              ! Maximum fraction of rain-hail mixture
  fog = 0.0              ! Maximum fraction of rain-hail mixture
      
  END SUBROUTINE init_fox_no_grpl


  SUBROUTINE init_fox_no_hail() 

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
!  Setup default maximum fraction of water in the melting ice 
!  when hail is suprressed. 
!
!-----------------------------------------------------------------------
!
! AUTHOR: Bryan Putnam, 12/14/10
!
!-----------------------------------------------------------------------
! Force explicit declarations. 
!-----------------------------------------------------------------------

  IMPLICIT NONE 

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval 
!-----------------------------------------------------------------------

  fos = 0.5             ! Maximum fraction of rain-snow mixture
  foh = 0.0             ! Maximum fraction of rain-hail mixture
  fog = 0.3             ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox_no_hail

  SUBROUTINE init_dsd()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default dsd values or reinialize default dsd values
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  N0r=8.0E+06 ! Intercept parameter in 1/(m^4) for rain.
  N0h=4.0E+04 ! Intercept parameter in 1/(m^4) for hail.
  N0s=3.0E+06 ! Intercept parameter in 1/(m^4) for snow.
  N0g=4.0E+05 ! Intercept parameter in 1/(m^4) for graupel.

  rhoh=913.  ! Density of hail (kg m**-3)
  rhos=100.  ! Density of snow (kg m**-3)
  rhog=400.  ! Density of graupel (kg m**-3)

  END SUBROUTINE init_dsd

  SUBROUTINE model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl, &
             alpharain,alphasnow,alphagrpl,alphahail)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set dsd values to those used in the arps forecasts
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------


  REAL :: n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl,   &
          alpharain,alphasnow,alphagrpl,alphahail

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r=n0rain
  N0s=n0snow
  N0h=n0hail
  N0g=n0grpl

  rhos=rhosnow
  rhoh=rhohail
  rhog=rhogrpl

  alphar = alpharain
  alphas = alphasnow
  alphag = alphagrpl
  alphah = alphahail

  END SUBROUTINE model_dsd

  SUBROUTINE coeff_hail(fw,qml)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for hail
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: fw, qml

! DTD: test 02/16/2012, do not adjust coefficient downward for small hail
!  IF(qml < q_threshold) THEN
!     sf = 4*qml*1.e3
!  ELSE
     sf = 0.8
!  ENDIF

  sigmah=pi/3*(1-sf*fw)
  Ah=.125*(3+4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Bh=.125*(3-4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Ch=.125*(1-exp(-8*sigmah**2))
  Dh=.125*(3+exp(-8*sigmah**2))
  Ckh=exp(-2*sigmah**2)

  END SUBROUTINE coeff_hail

SUBROUTINE coeff_grpl(fw,qml)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for graupel
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
! MODIFIED: Dan Dawson, 02/16/2012
!           Made separate version for graupel.
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: fw, qml

! DTD: test 02/16/2012, do not adjust coefficient downward for small hail
!  IF(qml < q_threshold) THEN
!     sf = 4*qml*1.e3
!  ELSE
     sf = 0.8
!  ENDIF

  sigmag=pi/3*(1-sf*fw)
  Ag=.125*(3+4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
  Bg=.125*(3-4*exp(-2*sigmag**2)+exp(-8*sigmag**2))
  Cg=.125*(1-exp(-8*sigmag**2))
  Dg=.125*(3+exp(-8*sigmag**2))
  Ckg=exp(-2*sigmag**2)

  END SUBROUTINE coeff_grpl

SUBROUTINE index_search(nvals,intv,val_list,param_val,ibeg,iend)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Subroutine to get the location of a value in an ordered list with a
! consistent interval between items.  
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 03/01/2015
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
 
  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nvals
  REAL, INTENT(IN) :: intv 
  REAL, INTENT(IN) :: val_list(nvals)
  REAL, INTENT(IN) :: param_val
  INTEGER, INTENT(OUT) :: ibeg,iend

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  IF(nvals == 1) THEN
    ibeg = 1
    iend = 1
  ELSE IF(param_val < val_list(1)) THEN
    ibeg = 1
    iend = 1
  ELSE IF(param_val >= val_list(nvals)) THEN
    ibeg = nvals
    iend = nvals
  ELSE
    ibeg = floor((param_val - val_list(1))/intv) + 1 
    iend = ibeg + 1  
  END IF 
 
  END SUBROUTINE index_search  

SUBROUTINE scatt_interp(x,x1,x2,y,y1,y2,f11,f12,f21,f22,res_val) 

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Perform a quadratic interpolation of pre-calculated scattering values
! based on a range of  slope and shape parameter values. 
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 03/01/2015
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: x,x1,x2,y,y1,y2,f11,f12,f21,f22

  REAL, INTENT(OUT) :: res_val 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(x1 == x2 .AND. y1 == y2) THEN
    res_val = f11
  ELSE IF(x1 == x2) THEN
    res_val = f11 + (f12-f11)*((y-y1)/(y2-y1))
  ELSE IF(y1 == y2) THEN
    res_val = f11 + (f21-f11)*((x-x1)/(x2-x1))
  ELSE
  
  res_val = (1/((x2-x1)*(y2-y1)))*(f11*(x2-x)*(y2-y) + &
            f21*(x-x1)*(y2-y) + f12*(x2-x)*(y-y1) + f22*(x-x1)*(y-y1))

  END IF 

  END SUBROUTINE scatt_interp


  TYPE(T_obs_dual) FUNCTION assign_Refl(var1,var2,var3,var4)
       REAL :: var1,var2,var3,var4,var5
       assign_Refl%T_sum_ref_h = var1
       assign_Refl%T_sum_ref_v = var2
       assign_Refl%T_log_zdr = var3
       assign_Refl%T_log_ref = var4
  END FUNCTION assign_Refl

  TYPE(T_obs_dual) FUNCTION init_Refl()
       init_Refl%T_sum_ref_h = 0.
       init_Refl%T_sum_ref_v = 0.
       init_Refl%T_log_zdr = missing
       init_Refl%T_log_ref = 0.
       init_Refl%T_sum_ref_hv = 0.
       init_Refl%T_kdp = 0.
       init_Refl%T_Ahh = 0.
       init_Refl%T_Avv = 0.
       init_Refl%T_ref_r_h = 0.
       init_Refl%T_ref_s_h = 0.
       init_Refl%T_ref_h_h = 0.
       init_Refl%T_ref_g_h = 0.
       init_Refl%T_ref_rs_h = 0.
       init_Refl%T_ref_rh_h = 0.
       init_Refl%T_ref_rg_h = 0.
       init_Refl%T_ref_r_v = 0.
       init_Refl%T_ref_s_v = 0.
       init_Refl%T_ref_h_v = 0.
       init_Refl%T_ref_g_v = 0.
       init_Refl%T_ref_rs_v = 0.
       init_Refl%T_ref_rh_v = 0.
       init_Refl%T_ref_rg_v = 0.
  END FUNCTION init_Refl

  TYPE(T_para_dsd) FUNCTION init_para_dsd()
    init_para_dsd%T_qr = 0.0
    init_para_dsd%T_qs = 0.0
    init_para_dsd%T_qh = 0.0
    init_para_dsd%T_qg = 0.0
    init_para_dsd%T_Ntr = 0.0
    init_para_dsd%T_Nts = 0.0
    init_para_dsd%T_Nth = 0.0
    init_para_dsd%T_Ntg = 0.0
    init_para_dsd%T_alfr = 0.0
    init_para_dsd%T_alfs = 0.0
    init_para_dsd%T_alfh = 0.0
    init_para_dsd%T_alfg = 0.0
  END FUNCTION init_para_dsd

  TYPE(T_para_dsd) FUNCTION assign_para_dsd_SM(var1,var2,var3,var4)
    REAL :: var1,var2,var3,var4

    assign_para_dsd_SM%T_qr = var1
    assign_para_dsd_SM%T_qs = var2
    assign_para_dsd_SM%T_qh = var3
    assign_para_dsd_SM%T_qg = var4
    assign_para_dsd_SM%T_Ntr = 0.0 
    assign_para_dsd_SM%T_Nts = 0.0 
    assign_para_dsd_SM%T_Nth = 0.0 
    assign_para_dsd_SM%T_Ntg = 0.0 
    assign_para_dsd_SM%T_alfr = 0.0
    assign_para_dsd_SM%T_alfs = 0.0
    assign_para_dsd_SM%T_alfh = 0.0
    assign_para_dsd_SM%T_alfg = 0.0
  END FUNCTION assign_para_dsd_SM

  TYPE(T_para_dsd) FUNCTION assign_para_dsd_TM(var1,var2,var3,var4, &
                            var5,var6,var7,var8,var9,var10,var11,var12)
    REAL :: var1,var2,var3,var4,var5,var6,var7,var8
    REAL :: var9,var10,var11,var12

    assign_para_dsd_TM%T_qr = var1
    assign_para_dsd_TM%T_qs = var2
    assign_para_dsd_TM%T_qh = var3
    assign_para_dsd_TM%T_qg = var4
    assign_para_dsd_TM%T_Ntr = var5
    assign_para_dsd_TM%T_Nts = var6
    assign_para_dsd_TM%T_Nth = var7
    assign_para_dsd_TM%T_Ntg = var8
    assign_para_dsd_TM%T_alfr = var9
    assign_para_dsd_TM%T_alfs = var10
    assign_para_dsd_TM%T_alfh = var11
    assign_para_dsd_TM%T_alfg = var12
  END FUNCTION assign_para_dsd_TM


TYPE(T_obs_dual) FUNCTION rainIceRefl(var_dsd,rho,flg)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the partial reflectivity factor
! of melting(wet) snow/hail at horizontal polarization
! and compute total reflectivity as a sum of those.
! The same formula used in shfactor is used with different
! alpha and beta coefficients that contain the effect of the fraction
! of water in the melting snow to take the melting layer into account.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  use radaremul_cst, only: mphyopt, MFflg

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: snow_alpha_a, hail_alpha_a, grpl_alpha_a
  REAL, EXTERNAL :: snow_alpha_b, hail_alpha_b, grpl_alpha_b
  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
 
  TYPE(T_para_dsd) :: var_dsd
  REAL :: qr,qs,qh,qg,rho,ntr,nts,nth,ntg 
  REAL :: rainIceRefl_hh,rainIceRefl_vv,rainIceRefl_hv,zdr
  REAL :: fracqrs,fracqrh,fracqrg
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg,fws,fwh,fwg,rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf
  REAL :: alphaa_ws,alphab_ws,alphaa_wh,alphab_wh,alphaa_wg,alphab_wg
  REAL :: alphak_ws,alphak_wh,alphak_wg
  REAL :: rainReflH,ZdrysnowH,ZwetsnowH
  REAL :: rainReflV,ZdrysnowV,ZwetsnowV
  REAL :: ZdryhailH,ZwethailH,ZdrygrplH,ZwetgrplH
  REAL :: ZdryhailV,ZwethailV,ZdrygrplV,ZwetgrplV
  REAL :: rainReflHV,ZdrysnowHV,ZwetsnowHV
  REAL :: ZdryhailHV,ZwethailHV,ZdrygrplHV,ZwetgrplHV
  REAL :: log_ref
  REAL :: rho_0rs,rho_0rh,rho_0rg,temp
  REAL :: temH,temV,temHV

  INTEGER :: flg

  REAL :: tair_C
  REAL :: z_snow_thom
  REAL*8 :: gamma,exp_term,gam_term

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  fracqs = 0.; fracqh = 0.; fracqg = 0.
  fracqrs = 0.; fracqrh = 0.; fracqrg = 0.

  fms = 0.; fmh = 0.; fmg = 0.
  fws = 0.; fwh = 0.; fwg = 0.
  rhoms = 100.; rhomh = 913.; rhomg = 400.

  rainReflH = 0.
  rainReflV = 0.
  rainReflHV = 0.
  ZdrysnowH = 0.
  ZdrysnowV = 0.
  ZdrysnowHV = 0.
  ZwetsnowH = 0.
  ZwetsnowV = 0.
  ZwetsnowHV = 0.
  ZdryhailH = 0.
  ZdryhailV = 0.
  ZdryhailHV = 0.
  ZwethailH = 0.
  ZwethailV = 0.
  ZwethailHV = 0.
  ZdrygrplH = 0.
  ZdrygrplV = 0.
  ZdrygrplHV = 0.
  ZwetgrplH = 0.
  ZwetgrplV = 0.
  ZwetgrplHV = 0.

  temH = 0.
  temV = 0.
  temHV = 0. 

  rainIceRefl_hh = 0.
  rainIceRefl_vv = 0.
  rainIceRefl_hv = 0.
  zdr = missing
  log_ref = 0.

  rho_0rs = rho_0rsf
  rho_0rh = rho_0rhf
  rho_0rg = rho_0rgf

  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg
  ntr = var_dsd%T_Ntr
  nts = var_dsd%T_Nts
  nth = var_dsd%T_Nth
  ntg = var_dsd%T_Ntg
  
  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
!   qrf  pure rain water mixing ratio
!   qsf  dry snow mixing ratio
!   qhf  dry hail mixing ratio
!   qgf  dry graupel mixing ratio
!   fms  wet snow mixing ratio
!   fmh  wet hail mixing ratio
!   fmg  wet graupel mixing ratio
!   rhoms  density of wet snow (kg/m-3)
!   rhomh  density of wet hail (kg/m-3)
!   rhomg  density of wet graupel (kg/m-3)
!-----------------------------------------------------------------------

  IF (MFflg == 0) THEN

    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hl_ON == 1)  &
      CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
    IF(grpl_ON == 1) &
      CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)

    qrf = qr - fracqrs - fracqrh - fracqrg
    if(qrf < 0.0) qrf = 0.0

    qsf = qs - fracqs
    if(qsf < 0.0) qsf = 0.0
    qhf = qh - fracqh
    if(qhf < 0.0) qhf = 0.0
    qgf = qg - fracqg
    if(qgf < 0.0) qgf = 0.0

  ELSE IF (MFflg == 2) THEN

    qrf = qr
    qsf = qs
    IF(hl_ON == 1)   qhf = qh
    IF(grpl_ON == 1) qgf = qg

  ELSE IF (MFflg == 3) THEN    ! Temperature-based melting.

    tair_C = ta - degKtoC

    CALL fractionWater_temperature_snow(qr,qs,rhos,fms,fws,rhoms,tair_C)
    IF(hl_ON == 1)  &
    CALL fractionWater_temperature_hail(qr,qh,rhoh,fmh,fwh,rhomh,tair_C)
    IF(grpl_ON == 1)  &
    CALL fractionWater_temperature_hail(qr,qg,rhog,fmg,fwg,rhomg,tair_C)

    qrf = qr
    qsf = qs-fms
    qhf = qh-fmh
    qgf = qg-fmg

  END IF

  qrf = qr - fracqrs - fracqrh - fracqrg
  if(qrf < 0.0) qrf = 0.0
  qsf = qs - fracqs
  if(qsf < 0.0) qsf = 0.0
  qhf = qh - fracqh
  if(qhf < 0.0) qhf = 0.0
  qgf = qg - fracqg
  if(qgf < 0.0) qgf = 0.0

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  IF(hl_ON == 1)   CALL coeff_hail(fwh,fmh)
  IF(grpl_ON == 1) CALL coeff_grpl(fwg,fmg)

!-----------------------------------------------------------------------
! Calculate alpha values
!-----------------------------------------------------------------------
  IF(fms > 0.) THEN
    alphaa_ws = snow_alpha_a(fws)
    alphab_ws = snow_alpha_b(fws)
    alphak_ws = alphaa_ws - alphab_ws
  ENDIF

  IF(hl_ON == 1 .and. fmh > 0.) THEN
    alphaa_wh = hail_alpha_a(fwh)
    alphab_wh = hail_alpha_b(fwh)
    alphak_wh = alphaa_wh - alphab_wh
  ENDIF

  IF(grpl_ON == 1 .and. fmg > 0.) THEN
    alphaa_wg = grpl_alpha_a(fwg)
    alphab_wg = grpl_alpha_b(fwg)
    alphak_wg = alphaa_wg - alphab_wg
  ENDIF

!-----------------------------------------------------------------------
! Calculate rho_0rs, rho_0rh, and rho_0rg
!-----------------------------------------------------------------------
  IF(flg > 2 .and. fms > 0.) THEN
    temp=rho*fms*1.e3
    if(temp > 1.) then
      rho_0rs = rho_0rsi
    else if (1.e-2 > temp .and. temp <= 1.) then
      rho_0rs = rho_0rsi - .5*log10(temp)*(rho_0rsf-rho_0rsi)
    endif
  ENDIF

  IF(hl_ON == 1 .and. flg > 2 .and. fmh > 0.) THEN
    temp=rho*fmh*1.e3
    if(temp > 1.) then
      rho_0rh = rho_0rhi
    else if (1.e-2 > temp .and. temp <= 1.) then
      rho_0rh = rho_0rhi - .5*log10(temp)*(rho_0rhf-rho_0rhi)
    endif
  ENDIF

  IF(grpl_ON == 1 .and. flg > 2 .and. fmg > 0.) THEN
    temp=rho*fmg*1.e3
    if(temp > 1.) then
      rho_0rg = rho_0rgi
    else if (1.e-2 > temp .and. temp <= 1.) then
      rho_0rg = rho_0rgi - .5*log10(temp)*(rho_0rgf-rho_0rgi)
    endif
  ENDIF

!-----------------------------------------------------------------------
! Calculate reflectivity (Zhh and Zvv (and Zhv, if necessary))
!-----------------------------------------------------------------------

  CALL calc_N0x_mp(rho,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,qsf,  &
                    fms,qhf,fmh,qgf,fmg)

  CALL calc_lamda_mp(rho,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,  &
                    qrf,qsf,fms,qhf,fmh,qgf,fmg) 


  SELECT CASE (mphyopt)
  CASE(1:12,106,109:110,116)
    IF(lamdas > 0.) THEN
      CALL partialRefIce(N0s,alphas,As,Bs,Cs,alphaa_ds,       &
                         alphab_ds,beta_sa, beta_sb,lamdas,   & 
                         ZdrysnowH,ZdrysnowV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0s,alphas,Cs,Ds,alphaa_ds,        &
                         alphab_ds,beta_sa,beta_sb,rho_0s,    &
                         lamdas,ZdrysnowHV)
      ENDIF
    ENDIF
    IF(lamdams > 0.) THEN
      CALL partialRefIce(N0ms,alphas,As,Bs,Cs,alphaa_ws,       &
                         alphab_ws,beta_sa,beta_sb,lamdams,    &
                         ZwetsnowH,ZwetsnowV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0ms,alphas,Cs,Ds,alphaa_ws,        &
                         alphab_ws,beta_sa,beta_sb,rho_0rs,    &
                         lamdams,ZwetsnowHV)
      ENDIF
    ENDIF
  CASE(108)
    IF(lamdas > 0. .and. qsf > 0.) THEN
      CALL partialRefIce(N0s,alphas,As,Bs,Cs,alphaa_tom_ds,     &
                        alphab_tom_ds,beta_tom_dsa,beta_tom_dsb,  &
                        lamdas,ZdrysnowH,ZdrysnowV)

      CALL partialRefIce(N0s2,alphas2,As,Bs,Cs,alphaa_tom_ds,     &
                        alphab_tom_ds,beta_tom_dsa,beta_tom_dsb,  &
                        lamdas2,temH,temV) 

      ZdrysnowH = ZdrysnowH + temH
      ZdrysnowV = ZdrysnowV + temV

      IF(flg > 2) THEN
        CALL partialRhoIce(N0s,alphas,Cs,Ds,alphaa_ds,        &
                        alphab_ds,beta_tom_dsa,beta_tom_dsb,    &
                        rho_0s,lamdas,ZdrysnowHV)
        CALL partialRhoIce(N0s2,alphas2,Cs,Ds,alphaa_ds,        &
                        alphab_ds,beta_tom_dsa,beta_tom_dsb,    &
                        rho_0s,lamdas2,temHV)

        ZdrysnowHV = ZdrysnowHV + temHV
      ENDIF
    END IF 
    IF(lamdams > 0. .and. fms > 0.) THEN
       CALL partialRefIce(N0ms,alphas,As,Bs,Cs,alphaa_ws,      &
                         alphab_ws,beta_sa,beta_sb,lamdams,    &
                         ZwetsnowH,ZwetsnowV) 
       IF(flg > 2) THEN
       CALL partialRhoIce(N0ms,alphas,Cs,Ds,alphaa_ws,           &
                         alphab_ws,beta_sa,beta_sb,rho_0rs,      &
                         lamdams,ZwetsnowHV)
       END IF
    ENDIF 
  END SELECT 


  IF(hl_ON == 1) THEN
    IF(lamdah > 0. .and. qhf > 0.)THEN
      CALL partialRefIce(N0h,alphah,Ahd,Bhd,Chd,alphaa_dh,    &
                         alphab_dh,beta_ha,beta_hb,lamdah,    &
                         ZdryhailH,ZdryhailV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0h,alphah,Chd,Dhd,alphaa_dh,      &
                         alphab_dh,beta_ha,beta_hb,rho_0h,    &
                         lamdah,ZdryhailHV)
      ENDIF
    ENDIF
    IF(lamdamh > 0. .and. fmh > 0.) THEN
      CALL partialRefIce(N0mh,alphah,Ah,Bh,Ch,alphaa_wh,       &
                         alphab_wh,beta_ha,beta_hb,lamdamh,    &
                         ZwethailH,ZwethailV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0mh,alphah,Ch,Dh,alphaa_wh,        &
                         alphab_wh,beta_ha,beta_hb,rho_0rh,    &
                         lamdamh,ZwethailHV)
      ENDIF
    ENDIF
  ENDIF 

  IF(grpl_ON == 1) THEN
    IF(lamdag > 0. .and. qgf > 0.)THEN
      CALL partialRefIce(N0g,alphag,Agd,Bgd,Cgd,alphaa_dg,    &
                         alphab_dg,beta_ga, beta_gb,lamdag,   &
                         ZdrygrplH,ZdrygrplV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0g,alphag,Cgd,Dgd,alphaa_dg,      &
                         alphab_dg,beta_ga,beta_gb,rho_0g,    &
                         lamdag,ZdrygrplHV)
      ENDIF
    ENDIF
     IF(lamdamg > 0. .and. fmg > 0.) THEN 
      CALL partialRefIce(N0mg,alphag,Ag,Bg,Cg,alphaa_wg,       &
                         alphab_wg,beta_ga,beta_gb,lamdamg,    &
                         ZwetgrplH,ZwetgrplV)
      IF(flg > 2) THEN
        CALL partialRhoIce(N0mg,alphag,Cg,Dg,alphaa_wg,        &
                         alphab_wg,beta_ga,beta_gb,rho_0rg,    &
                         lamdamg,ZwetgrplHV)
      ENDIF
    ENDIF
  ENDIF

  IF(lamdar > 0.) THEN
    CALL partialRefRain(N0r,alphar,alphaa,alphab,beta_ra,beta_rb,  &
                       lamdar,rainReflH,rainReflV)
    rainReflV = MIN(rainReflV,rainReflH)
    IF(flg > 2) THEN

    CALL partialRhoRain(N0r,alphar,alphaa,alphab,beta_ra,beta_rb,  &
                        lamdar,rainReflHV)
    ENDIF
  ENDIF

  rainIceRefl_hh=rainReflH+ZdrysnowH+ZwetsnowH+ZdryhailH+ZwethailH &
                 +ZdrygrplH+ZwetgrplH
 
  log_ref = 10.*LOG10(MAX(1.0,rainIceRefl_hh))

  IF(flg == 1) THEN
    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)

  ELSE IF(flg > 1) THEN
!-----------------------------------------------------------------------
! Calculate differential reflectivity (Zdr)
!-----------------------------------------------------------------------
    rainIceRefl_vv=rainReflV+ZdrysnowV+ZwetsnowV+ZdryhailV+ZwethailV &
                  +ZdrygrplV+ZwetgrplV

    if(rainIceRefl_vv > 0.) then
      zdr = 10.*LOG10(MAX(1.0,rainIceRefl_hh/rainIceRefl_vv))
    endif

    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)

    IF(flg > 2) THEN

      rainIceRefl_hv=rainReflHV+ZdrysnowHV+ZwetsnowHV                  &
                   +ZdryhailHV+ZwethailHV+ZdrygrplHV+ZwetgrplHV

!-----------------------------------------------------------------------
! Safety block to ensure r_hv <= 1.
!-----------------------------------------------------------------------
      IF(rainIceRefl_hv > SQRT(rainIceRefl_hh*rainIceRefl_vv)) &
         rainIceRefl_hv = SQRT(rainIceRefl_hh*rainIceRefl_vv)
!-----------------------------------------------------------------------

      rainIceRefl%T_sum_ref_hv = rainIceRefl_hv

    ENDIF
  ENDIF

END FUNCTION rainIceRefl


REAL FUNCTION rainIceKdp(var_dsd,rho)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates specific differential phase.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  use radaremul_cst, only: mphyopt

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: snow_alpha_a, hail_alpha_a, grpl_alpha_a
  REAL, EXTERNAL :: snow_alpha_b, hail_alpha_b, grpl_alpha_b
  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  TYPE(T_para_dsd) :: var_dsd

  REAL :: qr,qs,qh,qg,rho,kdp,ntr,nts,nth,ntg
  REAL :: fracqrs,fracqrh,fracqrg
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg,fws,fwh,fwg,rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf
  REAL :: alphaa_ws,alphab_ws,alphaa_wh,alphab_wh,alphaa_wg,alphab_wg
  REAL :: alphak_ws,alphak_wh,alphak_wg
  REAL :: rainKdp,drysnowKdp,wetsnowKdp
  REAL :: dryhailKdp,wethailKdp,drygrplKdp,wetgrplKdp
  REAL :: temKdp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  fracqs = 0.; fracqh = 0.; fracqg = 0.
  fracqrs = 0.; fracqrh = 0.; fracqrg = 0.

  fms = 0.; fmh = 0.; fmg = 0.
  fws = 0.; fwh = 0.; fwg = 0.
  rhoms = 100.; rhomh = 913.; rhomg = 400.

  drysnowKdp = 0.
  wetsnowKdp = 0.
  dryhailKdp = 0.
  wethailKdp = 0.
  drygrplKdp = 0.
  wetgrplKdp = 0.
  rainKdp = 0.
  rainIceKdp = 0.

  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg
  ntr = var_dsd%T_Ntr
  nts = var_dsd%T_Nts
  nth = var_dsd%T_Nth
  ntg = var_dsd%T_Ntg

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! For variable names, see FUNCTION rainIceRefl
!-----------------------------------------------------------------------
  CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
  IF(hl_ON == 1) &
    CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
  IF(grpl_ON == 1) &
    CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)

  qrf = qr - fracqrs - fracqrh - fracqrg
  if(qrf < 0.0) qrf = 0.0
  qsf = qs - fracqs
  if(qsf < 0.0) qsf = 0.0
  qhf = qh - fracqh
  if(qhf < 0.0) qhf = 0.0
  qgf = qg - fracqg
  if(qgf < 0.0) qgf = 0.0

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  IF(hl_ON == 1)   CALL coeff_hail(fwh,fmh)
  IF(grpl_ON == 1) CALL coeff_grpl(fwg,fmg)

!-----------------------------------------------------------------------
! Calculate alpha values
!-----------------------------------------------------------------------
  IF(fms > 0.) THEN
    alphaa_ws = snow_alpha_a(fws)
    alphab_ws = snow_alpha_b(fws)
    alphak_ws = alphaa_ws - alphab_ws
  ENDIF

  IF(hl_ON == 1 .and. fmh > 0.) THEN
    alphaa_wh = hail_alpha_a(fwh)
    alphab_wh = hail_alpha_b(fwh)
    alphak_wh = alphaa_wh - alphab_wh
  ENDIF

  IF(grpl_ON == 1 .and. fmg > 0.) THEN
    alphaa_wg = grpl_alpha_a(fwg)
    alphab_wg = grpl_alpha_b(fwg)
    alphak_wg = alphaa_wg - alphab_wg
  ENDIF

!-----------------------------------------------------------------------
! Calculate specific differential phase (Kdp)
!-----------------------------------------------------------------------
  CALL calc_N0x_mp(rho,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,qsf,      &
                   fms,qhf,fmh,qgf,fmg)
   
  CALL calc_lamda_mp(rho,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,          &
                     qrf,qsf,fms,qhf,fmh,qgf,fmg)


  SELECT CASE (mphyopt)
  CASE(1:12,106,109:110,116)
  IF(lamdas > 0.) THEN
    CALL partialKdpIce(kdpCoefIce,Cks,alphak_ds,betak_s,N0s,alphas,    &
                       lamdas,drysnowKdp)
  ENDIF

  IF(lamdams > 0.) THEN
    CALL partialKdpIce(kdpCoefIce,Cks,alphak_ws,betak_s,N0ms,alphas,   &
                       lamdams,wetsnowKdp)
  ENDIF
  CASE(108)
  IF(lamdas > 0.) THEN

     CALL partialKdpIce(kdpCoefIce,Cks,alphak_ds,betak_tom_ds,N0s,alphas,   &
                        lamdas,drysnowKdp)
     CALL partialKdpIce(kdpCoefIce,Cks,alphak_ds,betak_tom_ds,N0s2,alphas2,  &
                        lamdas2,temKdp)

     drysnowKdp = drysnowKdp + temKdp
   
  END IF
  IF(lamdams > 0.) THEN

     CALL partialKdpIce(kdpCoefIce,Cks,alphak_ws,betak_s,N0ms,alphas,lamdams,wetsnowKdp)
  
  END IF
  END SELECT 

  IF(hl_ON == 1) THEN
    IF(lamdah > 0.) THEN
      CALL partialKdpIce(kdpCoefIce,Ckhd,alphak_dh,betak_h,N0h,alphah,    &
                         lamdah,dryhailKdp)
    ENDIF

    IF(lamdamh > 0.) THEN
      CALL partialKdpIce(kdpCoefIce,Ckh,alphak_wh,betak_h,N0mh,alphah,    &
                         lamdamh,wethailKdp)
    ENDIF
  ENDIF 

  IF(grpl_ON == 1) THEN
    IF(lamdag > 0.) THEN
      CALL partialKdpIce(kdpCoefIce,Ckgd,alphak_dg,betak_g,N0g,alphag,    &
                         lamdag,drygrplKdp)
    ENDIF

    IF(lamdamg > 0.) THEN
      CALL partialKdpIce(kdpCoefIce,Ckg,alphak_wg,betak_g,N0mg,alphag,    &
                         lamdamg,wetgrplKdp)
    ENDIF
  ENDIF

  IF(lamdar > 0.) THEN
   
    CALL partialKdpRain(constKdpr,N0r,alphar,lamdar,rainKdp)
    IF(rainKdp < 0.0) rainKdp = 0.0
  ENDIF

  rainIceKdp=rainKdp+drysnowKdp+wetsnowKdp+dryhailKdp+wethailKdp &
            +drygrplKdp+wetgrplKdp


END FUNCTION rainIceKdp

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION calculate_obs                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

 TYPE(T_obs_dual) FUNCTION calculate_obs(rho,var_dsd,flg)

!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
! flg == (1: Zh, 2: Zdr, 3: rho_hv)
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho ! Air density (kg m**-3)

  TYPE(T_para_dsd) :: var_dsd

  INTEGER, INTENT(IN) :: flg   ! flag for ref(1) and zdr(2)

  REAL :: qr
  REAL :: qs
  REAL :: qh
  REAL :: qg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  calculate_obs = init_Refl()

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------
 
  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg


  IF (rho > 0.0 .and. (qr > 0. .or. qs > 0. .or. qh > 0. &
      .or. qg > 0.)) THEN 
    calculate_obs = rainIceRefl(var_dsd,rho,flg)

  END IF

END FUNCTION  calculate_obs

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION refl_rsa                       #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

TYPE(T_obs_dual) FUNCTION refl_rsa (rhoa,var_dsd)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Compute radar observations by integrating radar scattering amplitudes
! over the drop size range. Radar scattering amplitues were calculated
! using T-matrix method and stored in the table.
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  USE scatt_table
  USE radaremul_cst, only: mphyopt, MFflg

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, PARAMETER :: pi4 = 97.409           ! pi^4
  REAL, DIMENSION (nd) :: Ndr, Nds, Ndh, Ndg, Ndrs, Ndrh, Ndrg

  INTEGER, EXTERNAL :: get_qgh_opt

  REAL :: rhoa
  REAL :: qr,qs,qh,qg
  REAL :: ntr,nts,nth,ntg

  REAL*8 :: alfr,alfs,alfh,alfg,alfs2

  REAL*8 :: db_N0,db_N02
  REAL :: Ntw,Ntd,temNtd,temNtw
 
  REAL*8 :: gamma

  REAL :: fracqrs,fracqrh,fracqrg
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg
  REAL :: fws,fwh,fwg
  REAL :: rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf
  REAL*8 :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL*8 :: lams2,lamrs2
  REAL :: tsar_h,tsas_h,tsah_h,tsag_h,tsars_h,tsarh_h,tsarg_h
  REAL :: tsar_v,tsas_v,tsah_v,tsag_v,tsars_v,tsarh_v,tsarg_v
  COMPLEX :: tsar_hv,tsas_hv,tsah_hv,tsag_hv,tsars_hv,tsarh_hv,tsarg_hv
  REAL :: tfsar,tfsas,tfsah,tfsag,tfsars,tfsarh,tfsarg
  REAL :: D,intv,far,lambda4
  REAL :: fa2,fb2,temph,tempv,temphv,tempk,temp
  REAL :: scatth_val,scattv_val,scatthv_val,scattk_val
  COMPLEX :: fab,fba,fconj
  INTEGER :: i,j,k,m,idx
  REAL :: tempAhh,tempAvv
  REAL :: Ar_h,As_h,Ag_h,Ah_h,Ars_h,Arg_h,Arh_h
  REAL :: Ar_v,As_v,Ag_v,Ah_v,Ars_v,Arg_v,Arh_v

  TYPE(T_para_dsd) :: var_dsd

  REAL*8 :: term_exp, term_gam

  REAL :: tair_C

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  lambda4 = lambda**4.

  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg
  ntr = var_dsd%T_Ntr
  nts = var_dsd%T_Nts
  nth = var_dsd%T_Nth
  ntg = var_dsd%T_Ntg
  alfr = dble(var_dsd%T_alfr)
  alfs = dble(var_dsd%T_alfs)
  alfh = dble(var_dsd%T_alfh)
  alfg = dble(var_dsd%T_alfg)
  alfs2 = dble(alphas2)

  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  Ndr = 0.; Nds = 0.; Ndh = 0.; Ndg = 0.
  Ndrs = 0.; Ndrh = 0.; Ndrg = 0.
  temph = 0.; tempv = 0.; temphv = 0.; temp = 0.; tempk = 0.
  tempAhh = 0.; tempAvv = 0.
  fracqrs = 0.; fracqs = 0.; fms = 0.; fws = 0.; rhoms = 100.
  fracqrh = 0.; fracqh = 0.; fmh = 0.; fwh = 0.; rhomh = 913.
  fracqrg = 0.; fracqg = 0.; fmg = 0.; fwg = 0.; rhomg = 400.
  refl_rsa = init_Refl()

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
! For the variable definition, see "FUNCTION rainIceRefl".
!-----------------------------------------------------------------------

  IF (MFflg == 0) THEN

    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hl_ON == 1)  &
      CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
    IF(grpl_ON == 1) &
      CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)

    qrf = qr - fracqrs - fracqrh - fracqrg
    if(qrf < 0.0) qrf = 0.0

    qsf = qs - fracqs
    if(qsf < 0.0) qsf = 0.0
    qhf = qh - fracqh
    if(qhf < 0.0) qhf = 0.0
    qgf = qg - fracqg
    if(qgf < 0.0) qgf = 0.0

  ELSE IF (MFflg == 2) THEN

    qrf = qr 
    qsf = qs; fms = 0.0; 
    IF(hl_ON == 1)   qhf = qh; fmh = 0.0;
    IF(grpl_ON == 1) qgf = qg; fmg = 0.0;

  ELSE IF (MFflg == 3) THEN    ! Temperature-based melting.

    tair_C = ta - degKtoC

    CALL fractionWater_temperature_snow(qr,qs,rhos,fms,fws,rhoms,tair_C)
    IF(hl_ON == 1)  &
    CALL fractionWater_temperature_hail(qr,qh,rhoh,fmh,fwh,rhomh,tair_C)
    IF(grpl_ON == 1)  &
    CALL fractionWater_temperature_hail(qr,qg,rhog,fmg,fwg,rhomg,tair_C)
          
    qrf = qr 
    qsf = qs-fms
    qhf = qh-fmh
    qgf = qg-fmg
        
  END IF

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  IF(hl_ON == 1)   CALL coeff_hail(fwh,fmh)
  IF(grpl_ON == 1) CALL coeff_grpl(fwg,fmg)

!-----------------------------------------------------------------------
! Calculate slopes of the invers exponential curves of N(d)
!-----------------------------------------------------------------------

  CALL calc_N0x_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,qsf,   &
                     fms,qhf,fmh,qgf,fmg)


  CALL calc_lamda_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,qsf, &
                     fms,qhf,fmh,qgf,fmg)

  lamr = dble(lamdar)
  lams = dble(lamdas)
  lamh = dble (lamdah)
  lamrs = dble(lamdams)
  lamrh = dble(lamdamh)
  lams2 = dble(lamdas2)
  lamrs2 = dble(lamdams2) 
  lamg = dble(lamdag)
  lamrg = dble(lamdamg)

  db_N0 = 0.d0
  db_N02 = 0.d0 

  IF(dualpol_opt == 2) THEN
    if(qrf > 0.) then
      intv = (dsr(2) - dsr(1))*1.e-3
      Ntw = 0.
      if(ntr > 0.0) then
        Ntw = ntr
      else
        db_N0 = dble(N0r)
        CALL cal_Nt(rhoa,qrf,db_N0,c_x(2),alfr,Ntw)
      endif

       db_N0 = dble(N0r)
     

      Do i = 1,nd
        Ndr(i) = sngl(db_N0*(dsr(i)*1.e-3)**alfr*exp(-lamr*(dsr(i)*1.e-3))*intv)
        Ntw = Ntw - Ndr(i)
        if(Ntw <= 0) Ndr(i) = 0.
      ENDDO
    endif

    db_N0 = 0.d0

    SELECT CASE (mphyopt)
    CASE(1:12,106,109:110,116)
      if(qsf > 0.) then
        intv = (dss(2) - dss(1))*1.e-3 
        Ntd = 0.
      if(nts > 0.0) then
        Ntd = nts
      else
        db_N0 = dble(N0s)
        CALL cal_Nt(rhoa,qsf,db_N0,c_x(4),alfs,Ntd)
      endif
       db_N0 = dble(N0s)

       Do i = 1,nd
         Nds(i) = sngl(db_N0*(dss(i)*1.e-3)**alfs*exp(-lams*(dss(i)*1.e-3))*intv)
         Ntd = Ntd - Nds(i)
         if(Ntd <= 0) Nds(i) = 0.
       ENDDO
     endif

      db_N0 = 0.d0

     if(fms > 0.) then
       intv = (dss(2) - dss(1))*1.e-3
        Ntw = 0.
      if(nts > 0.0) then
        Ntw = nts
      else
        db_N0 = dble(N0ms)
        CALL cal_Nt(rhoa,fms,db_N0,c_x(4),alfs,Ntw)
      endif
       
       db_N0 = dble(N0ms)

       Do i = 1,nd
         Ndrs(i) = sngl(db_N0*(dss(i)*1.e-3)**alfs*exp(-lamrs*(dss(i)*1.e-3))*intv)
         Ntw = Ntw - Ndrs(i)
         if(Ntw <= 0) Ndrs(i) = 0.
       ENDDO
     endif
    CASE(108)
      if(qsf > 0.) then
        intv = (dss(2) - dss(1))*1.e-3
     
        Ntd = 0. 
        db_N0 = dble(N0s)
        db_N02 = dble(N0s2)

        temNtd = sngl((db_N0/lams))
        Ntd = sngl(dble(temNtd) + (db_N02/(lams2**(1.d0+alfs2)))*gamma(1.d0+alfs2))

        db_N0 = dble(N0s)
        db_N02 = dble(N0s2)

       DO i = 1,nd
        term_exp = db_N0*exp(-lams*(dss(i)*1.e-3))*intv*unit_factor
        term_gam = db_N02*(dss(i)*1.e-3)**alfs2*exp(-lams2*(dss(i)*1.e-3))*intv*unit_factor
        Nds(i) = sngl(term_exp + term_gam)
        Ntd = Ntd - Nds(i)
        if(Ntd <= 0) Nds(i) = 0. 
       enddo
      endif 

      if(fms > 0.) then
         intv = (dss(2) - dss(1))*1.e-3
   
        Ntw = 0.
         if(nts > 0.0) then
            Ntw = nts
         else
           db_N0 = dble(N0ms)
           CALL cal_Nt(rhoa,fms,db_N0,c_x(4),alfs,Ntw)
         end if             

         db_N0 = dble(N0ms)

       DO i = 1,nd
          Ndrs(i) = sngl(db_N0*(dss(i)*1.e-3)**alfs*exp(-lamrs*(dss(i)*1.e-3))* &
                    intv)
         Ntw = Ntw - Ndrs(i)
         if(Ntw <= 0) Ndrs(i) = 0.
        enddo
      endif
    END SELECT 
    

    IF(hl_ON == 1) THEN  
      db_N0 = 0.d0

      if(qhf > 0.) then
        intv = (dsh(2) - dsh(1))*1.e-3
        Ntd = 0.
        if(nth > 0.0) then
          Ntd = nth 
        else
          db_N0 = dble(N0h)
          CALL cal_Nt(rhoa,qhf,db_N0,c_x(6),alfh,Ntd)
        endif

         db_N0 = dble(N0h)

        Do i = 1,nd
          Ndh(i) = sngl(db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamh*(dsh(i)*1.e-3))*intv)
          Ntd = Ntd - Ndh(i)
          if(Ntd <= 0) Ndh(i) = 0.
        ENDDO
      endif

      db_N0 = 0.d0

      if(fmh > 0.) then
        intv = (dsh(2) - dsh(1))*1.e-3
        Ntw = 0.
        if(nth > 0.0) then
          Ntw = nth
        else
          db_N0 = dble(N0mh)
          CALL cal_Nt(rhoa,fmh,db_N0,c_x(6),alfh,Ntw)
        endif

         db_N0 = dble(N0mh)
      
        Do i = 1,nd
          Ndrh(i) = sngl(db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamrh*(dsh(i)*1.e-3))*intv)
          Ntw = Ntw - Ndrh(i)
          if(Ntw <= 0) Ndrh(i) = 0.
        ENDDO
      endif
    ENDIF 

    IF(grpl_ON == 1) THEN
      db_N0 = 0.d0

      if(qgf > 0.) then
        intv = (dsg(2) - dsg(1))*1.e-3
        Ntd = 0.
        if(ntg > 0.0) then
          Ntd = ntg
        else
          db_N0 = dble(N0g)
          CALL cal_Nt(rhoa,qgf,db_N0,c_x(5),alfg,Ntd)
        endif

          db_N0 = dble(N0g)

        Do i = 1,nd
          Ndg(i) = sngl(db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamg*(dsg(i)*1.e-3))*intv)
          Ntd = Ntd - Ndg(i)
          if(Ntd <= 0) Ndg(i) = 0.
        ENDDO
      endif

      db_N0 = 0.d0

      if(fmg > 0.) then
        intv = (dsg(2) - dsg(1))*1.e-3
        Ntw = 0.
        if(ntg > 0.0) then
          Ntw = ntg
        else
          db_N0 = dble(N0mg)
          CALL cal_Nt(rhoa,fmg,db_N0,c_x(5),alfg,Ntw)
        endif

          db_N0 = dble(N0mg) 

        Do i = 1,nd
          Ndrg(i) = sngl(db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3))*intv)
          Ntw = Ntw - Ndrg(i)
          if(Ntw <= 0) Ndrg(i) = 0.
        ENDDO
      endif
    ENDIF
  END IF

!-----------------------------------------------------------------------
! Calculate radar observations.
!-----------------------------------------------------------------------
  tsar_h=0.; tsas_h=0.; tsah_h=0.; tsag_h=0.
  tsars_h=0.; tsarh_h=0.; tsarg_h=0.
  tsar_v=0.; tsas_v=0.; tsah_v=0.; tsag_v=0.
  tsars_v=0.; tsarh_v=0.; tsarg_v=0.
  tsar_hv=0.; tsas_hv=0.; tsah_hv=0.; tsag_hv = 0.
  tsars_hv=0.; tsarh_hv=0.; tsarg_hv=0.
  tfsar=0.; tfsas=0.; tfsah=0.; tfsag=0.
  tfsars=0.; tfsarh=0.; tfsarg=0.
  fa2=0.; fb2=0.; fab=0.; far=0.
  Ar_h=0.; As_h=0.; Ag_h=0.; Ah_h=0.
  Ars_h=0.; Arg_h=0.; Arh_h=0.
  Ar_v=0.; As_v=0.; Ag_v=0.; Ah_v=0.
  Ars_v=0.; Arg_v=0.; Arh_v=0.

  IF(dualpol_opt == 2) THEN
    if(qrf > 0.) then
      do i=1,nd
        fa2 = ABS(far_b(i))**2
        fb2 = ABS(fbr_b(i))**2
        fab = far_b(i)*CONJG(fbr_b(i))
        fba = fbr_b(i)*CONJG(far_b(i))
        tsar_h = tsar_h + fa2*Ndr(i)
        tsar_v = tsar_v + fb2*Ndr(i)
        tsar_hv = tsar_hv + fab*Ndr(i)
        far=REAL(far_f(i) - fbr_f(i))
        tfsar = tfsar + far*Ndr(i)
        IF(attn_ON == 2) THEN
          Ar_h = Ar_h + AIMAG(far_f(i)*Ndr(i))
          Ar_v = Ar_v + AIMAG(fbr_f(i)*Ndr(i))
        ENDIF
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.

    if(qsf > 0.) then
      do i=1,nd
        fa2 = ABS(fas_b(i,1))**2
        fb2 = ABS(fbs_b(i,1))**2
        fab = fas_b(i,1)*CONJG(fbs_b(i,1))
        fba = fbs_b(i,1)*CONJG(fas_b(i,1))
        far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
        tsas_h = tsas_h + far*Nds(i)
        far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
        tsas_v = tsas_v + far*Nds(i)
        fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
        tsas_hv = tsas_hv + fconj*Nds(i)
        far=Cks*REAL(fas_f(i,1) - fbs_f(i,1))
        tfsas = tfsas + far*Nds(i)
        IF(attn_ON == 2) THEN
          As_h = As_h + AIMAG(fas_f(i,1)*Nds(i))
          As_v = As_v + AIMAG(fbs_f(i,1)*Nds(i))
        ENDIF
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.

    if(fms > 0.) then
      idx = INT(fws * 20 + 0.5) + 1
      do i=1,nd
        fa2 = ABS(fas_b(i,idx))**2
        fb2 = ABS(fbs_b(i,idx))**2
        fab = fas_b(i,idx)*CONJG(fbs_b(i,idx))
        fba = fbs_b(i,idx)*CONJG(fas_b(i,idx))
        far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
        tsars_h = tsars_h + far*Ndrs(i)
        far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
        tsars_v = tsars_v + far*Ndrs(i)
        fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
        tsars_hv = tsars_hv + fconj*Ndrs(i)
        far=Cks*REAL(fas_f(i,idx)-fbs_f(i,idx))
        tfsars = tfsars + far*Ndrs(i)
        IF(attn_ON == 2) THEN
          Ars_h = Ars_h + AIMAG(fas_f(i,idx)*Ndrs(i))
          Ars_v = Ars_v + AIMAG(fbs_f(i,idx)*Ndrs(i))
        ENDIF
      enddo
    endif

    IF(hl_ON == 1) THEN
      fa2=0.; fb2=0.; fab=0.; far=0.

      if(qhf > 0.) then
        do i=1,nd
          fa2 = ABS(fah_b(i,1))**2
          fb2 = ABS(fbh_b(i,1))**2
          fab = fah_b(i,1)*CONJG(fbh_b(i,1))
          fba = fbh_b(i,1)*CONJG(fah_b(i,1))
          far=(Ahd*fa2 + Bhd*fb2 + 2*Chd*REAL(fab))
          tsah_h = tsah_h + far*Ndh(i)
          far=(Bhd*fa2 + Ahd*fb2 + 2*Chd*REAL(fab))
          tsah_v = tsah_v + far*Ndh(i)
          fconj=(Chd*(fa2+fb2)+Ahd*fab+Bhd*fba)
          tsah_hv = tsah_hv + fconj*Ndh(i)
          far=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
          tfsah = tfsah + far*Ndh(i)
          IF(attn_ON == 2) THEN
            Ah_h = Ah_h + AIMAG(fah_f(i,1)*Ndh(i))
            Ah_v = Ah_v + AIMAG(fbh_f(i,1)*Ndh(i))
          ENDIF
        enddo
      endif

      fa2=0.; fb2=0.; fab=0.; far=0.

      if(fmh > 0.) then
        idx = INT(fwh * 20 + 0.5) + 1
        do i=1,nd
          fa2 = ABS(fah_b(i,idx))**2
          fb2 = ABS(fbh_b(i,idx))**2
          fab = fah_b(i,idx)*CONJG(fbh_b(i,idx))
          fba = fbh_b(i,idx)*CONJG(fah_b(i,idx))
          far=(Ah*fa2 + Bh*fb2 + 2*Ch*REAL(fab))
          tsarh_h = tsarh_h + far*Ndrh(i)
          far=(Bh*fa2 + Ah*fb2 + 2*Ch*REAL(fab))
          tsarh_v = tsarh_v + far*Ndrh(i)
          fconj=(Ch*(fa2+fb2)+Ah*fab+Bh*fba)
          tsarh_hv = tsarh_hv + fconj*Ndrh(i)
          far=Ckh*REAL(fah_f(i,idx)-fbh_f(i,idx))
          tfsarh = tfsarh + far*Ndrh(i)
          IF(attn_ON == 2) THEN
            Arh_h = Arh_h + AIMAG(fah_f(i,idx)*Ndrh(i))
            Arh_v = Arh_v + AIMAG(fbh_f(i,idx)*Ndrh(i))
          ENDIF
        enddo
      endif
    ENDIF 

    IF(grpl_ON == 1) THEN
      fa2=0.; fb2=0.; fab=0.; far=0.

      if(qgf > 0.) then
        do i=1,nd
          fa2 = ABS(fag_b(i,1))**2
          fb2 = ABS(fbg_b(i,1))**2
          fab = fag_b(i,1)*CONJG(fbg_b(i,1))
          fba = fbg_b(i,1)*CONJG(fag_b(i,1))
          far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*REAL(fab))
          tsag_h = tsag_h + far*Ndg(i)
          far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*REAL(fab))
          tsag_v = tsag_v + far*Ndg(i)
          fconj=(Cgd*(fa2+fb2)+Agd*fab+Bgd*fba)
          tsag_hv = tsag_hv + fconj*Ndg(i)
          far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
          tfsag = tfsag + far*Ndg(i)
          IF(attn_ON == 2) THEN
            Ag_h = Ag_h + AIMAG(fag_f(i,1)*Ndg(i))
            Ag_v = Ag_v + AIMAG(fbg_f(i,1)*Ndg(i))
          ENDIF
        enddo
      endif

      fa2=0.; fb2=0.; fab=0.; far=0.

      if(fmg > 0.) then
        idx = INT(fwg * 20 + 0.5) + 1
        do i=1,nd
          fa2 = ABS(fag_b(i,idx))**2
          fb2 = ABS(fbg_b(i,idx))**2
          fab = fag_b(i,idx)*CONJG(fbg_b(i,idx))
          fba = fbg_b(i,idx)*CONJG(fag_b(i,idx))
          far=(Ag*fa2 + Bg*fb2 + 2*Cg*REAL(fab))
          tsarg_h = tsarg_h + far*Ndrg(i)
          far=(Bg*fa2 + Ag*fb2 + 2*Cg*REAL(fab))
          tsarg_v = tsarg_v + far*Ndrg(i)
          fconj=(Cg*(fa2+fb2)+Ag*fab+Bg*fba)
          tsarg_hv = tsarg_hv + fconj*Ndrg(i)
          far=Ckg*REAL(fag_f(i,idx)-fbg_f(i,idx))
          tfsarg = tfsarg + far*Ndrg(i)
          IF(attn_ON == 2) THEN
            Arg_h = Arg_h + AIMAG(fag_f(i,idx)*Ndrg(i))
            Arg_v = Arg_v + AIMAG(fbg_f(i,idx)*Ndrg(i))
          ENDIF
        enddo
      endif
    ENDIF
  endif !dualpol opt 2 

  if(dualpol_opt == 3) then

    i = 0; j = 0; k = 0; m = 0

    if(qrf > 0.) then   
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlrl,lam_list(:,1),sngl(lamr),i,j)  

      CALL index_search(n_alpha,intvlra,alf_list(:,1),sngl(alfr),k,m)

      CALL scatt_interp(sngl(lamr),lam_list(i,1),lam_list(j,1),sngl(alfr),alf_list(k,1),alf_list(m,1),&
                        scatth_r(i,k),scatth_r(i,m),scatth_r(j,k), &
                        scatth_r(j,m),scatth_val)
 
      tsar_h = N0r*scatth_val 

      CALL scatt_interp(sngl(lamr),lam_list(i,1),lam_list(j,1),sngl(alfr),alf_list(k,1),alf_list(m,1), &
                        scattv_r(i,k),scattv_r(i,m),scattv_r(j,k), &
                        scattv_r(j,m),scattv_val)

      tsar_v = N0r*scattv_val

      CALL scatt_interp(sngl(lamr),lam_list(i,1),lam_list(j,1),sngl(alfr),alf_list(k,1),alf_list(m,1), &
                        scatthv_r(i,k),scatthv_r(i,m),scatthv_r(j,k), &
                        scatthv_r(j,m),scatthv_val)

      tsar_hv = N0r*scatthv_val

      CALL scatt_interp(sngl(lamr),lam_list(i,1),lam_list(j,1),sngl(alfr),alf_list(k,1),alf_list(m,1), &
                        scattk_r(i,k),scattk_r(i,m),scattk_r(j,k), &
                        scattk_r(j,m),scattk_val)

      tfsar = N0r*scattk_val 
    end if !end qrf 

    if(qsf > 0.) then
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0. 

      CALL index_search(n_lam,intvlsl,lam_list(:,2),sngl(lams),i,j)
      CALL index_search(n_alpha,intvlsa,alf_list(:,2),sngl(alfs),k,m)
      CALL scatt_interp(sngl(lams),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2), &
                        scatth_s(i,k,1),scatth_s(i,m,1),scatth_s(j,k,1), &
                        scatth_s(j,m,1),scatth_val)

      tsas_h = N0s*scatth_val

      CALL scatt_interp(sngl(lams),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scattv_s(i,k,1),scattv_s(i,m,1),scattv_s(j,k,1), &
                        scattv_s(j,m,1),scattv_val)

      tsas_v = N0s*scattv_val

      CALL scatt_interp(sngl(lams),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scatthv_s(i,k,1),scatthv_s(i,m,1),scatthv_s(j,k,1), &
                        scatthv_s(j,m,1),scatthv_val)

      tsas_hv = N0s*scatthv_val

      CALL scatt_interp(sngl(lams),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scattk_s(i,k,1),scattk_s(i,m,1),scattk_s(j,k,1), &
                        scattk_s(j,m,1),scattk_val)

      tfsas = N0s*scattk_val
      
    end if !end qsf 

    if(fms > 0.0) then
      idx = INT(fws * 20 + 0.5) + 1 
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlsl,lam_list(:,2),sngl(lamrs),i,j)
      CALL index_search(n_alpha,intvlsa,alf_list(:,2),sngl(alfs),k,m)
      CALL scatt_interp(sngl(lamrs),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scatth_s(i,k,idx),scatth_s(i,m,idx),scatth_s(j,k,idx), &
                        scatth_s(j,m,idx),scatth_val)

      tsars_h = N0ms*scatth_val

      CALL scatt_interp(sngl(lamrs),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scattv_s(i,k,idx),scattv_s(i,m,idx),scattv_s(j,k,idx), &
                        scattv_s(j,m,idx),scattv_val)

      tsars_v = N0ms*scattv_val

      CALL scatt_interp(sngl(lamrs),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                    scatthv_s(i,k,idx),scatthv_s(i,m,idx),scatthv_s(j,k,idx), &
                        scatthv_s(j,m,idx),scatthv_val)

      tsars_hv = N0ms*scatthv_val

      CALL scatt_interp(sngl(lamrs),lam_list(i,2),lam_list(j,2),sngl(alfs),alf_list(k,2),alf_list(m,2),&
                        scattk_s(i,k,idx),scattk_s(i,m,idx),scattk_s(j,k,idx), &
                        scattk_s(j,m,idx),scattk_val)

      tfsars = N0ms*scattk_val
     
    end if !end fms

    if(qhf > 0.0) then
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlhl,lam_list(:,3),sngl(lamh),i,j)
      CALL index_search(n_alpha,intvlha,alf_list(:,3),sngl(alfh),k,m)
      CALL scatt_interp(sngl(lamh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3), &
                        scatth_ha(i,k,1),scatth_ha(i,m,1),scatth_ha(j,k,1), &
                        scatth_ha(j,m,1),scatth_val)

      tsah_h = N0h*scatth_val

      CALL scatt_interp(sngl(lamh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
                        scattv_ha(i,k,1),scattv_ha(i,m,1),scattv_ha(j,k,1), &
                        scattv_ha(j,m,1),scattv_val)

      tsah_v = N0h*scattv_val

      CALL scatt_interp(sngl(lamh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
                        scatthv_ha(i,k,1),scatthv_ha(i,m,1),scatthv_ha(j,k,1), &
                        scatthv_ha(j,m,1),scatthv_val)

      tsah_hv = N0h*scatthv_val

      CALL scatt_interp(sngl(lamh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
                        scattk_ha(i,k,1),scattk_ha(i,m,1),scattk_ha(j,k,1), &
                        scattk_ha(j,m,1),scattk_val)

      tfsah = N0h*scattk_val

    end if !end qhf

    if(fmh > 0.0) then
      idx = INT(fwh * 20 + 0.5) + 1
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlhl,lam_list(:,3),sngl(lamrh),i,j)
      CALL index_search(n_alpha,intvlha,alf_list(:,3),sngl(alfh),k,m)
      CALL scatt_interp(sngl(lamrh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
                     scatth_ha(i,k,idx),scatth_ha(i,m,idx),scatth_ha(j,k,idx), &
                        scatth_ha(j,m,idx),scatth_val)

      tsarh_h = N0mh*scatth_val

      CALL scatt_interp(sngl(lamrh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
                    scattv_ha(i,k,idx),scattv_ha(i,m,idx),scattv_ha(j,k,idx), &
                        scattv_ha(j,m,idx),scattv_val)

      tsarh_v = N0mh*scattv_val

      CALL scatt_interp(sngl(lamrh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
            scatthv_ha(i,k,idx),scatthv_ha(i,m,idx),scatthv_ha(j,k,idx), &
                        scatthv_ha(j,m,idx),scatthv_val)

      tsarh_hv = N0mh*scatthv_val

      CALL scatt_interp(sngl(lamrh),lam_list(i,3),lam_list(j,3),sngl(alfh),alf_list(k,3),alf_list(m,3),&
              scattk_ha(i,k,idx),scattk_ha(i,m,idx),scattk_ha(j,k,idx), &
                        scattk_ha(j,m,idx),scattk_val)

      tfsarh = N0mh*scattk_val

    end if !end fms

    if(qgf > 0.0) then
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlgl,lam_list(:,4),sngl(lamg),i,j)
      CALL index_search(n_alpha,intvlga,alf_list(:,4),sngl(alfg),k,m)
      CALL scatt_interp(sngl(lamg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4), &
                        scatth_g(i,k,1),scatth_g(i,m,1),scatth_g(j,k,1), &
                        scatth_g(j,m,1),scatth_val)

      tsag_h = N0g*scatth_val

      CALL scatt_interp(sngl(lamg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scattv_g(i,k,1),scattv_g(i,m,1),scattv_g(j,k,1), &
                        scattv_g(j,m,1),scattv_val)

      tsag_v = N0g*scattv_val

      CALL scatt_interp(sngl(lamg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scatthv_g(i,k,1),scatthv_g(i,m,1),scatthv_g(j,k,1), &
                        scatthv_g(j,m,1),scatthv_val)

      tsag_hv = N0g*scatthv_val

      CALL scatt_interp(sngl(lamg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scattk_g(i,k,1),scattk_g(i,m,1),scattk_g(j,k,1), &
                        scattk_g(j,m,1),scattk_val)

      tfsag = N0g*scattk_val

    end if !end qgf

    if(fmg > 0.0) then
      idx = INT(fwg * 20 + 0.5) + 1   
      scatth_val = 0.
      scattv_val = 0.
      scatthv_val = 0.
      scattk_val = 0.

      CALL index_search(n_lam,intvlgl,lam_list(:,4),sngl(lamrg),i,j)
      CALL index_search(n_alpha,intvlga,alf_list(:,4),sngl(alfg),k,m) 
      CALL scatt_interp(sngl(lamrg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scatth_g(i,k,idx),scatth_g(i,m,idx),scatth_g(j,k,idx), &
                        scatth_g(j,m,idx),scatth_val)

      tsarg_h = N0mg*scatth_val

      CALL scatt_interp(sngl(lamrg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scattv_g(i,k,idx),scattv_g(i,m,idx),scattv_g(j,k,idx), &
                        scattv_g(j,m,idx),scattv_val)

      tsarg_v = N0mg*scattv_val

      CALL scatt_interp(sngl(lamrg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                    scatthv_g(i,k,idx),scatthv_g(i,m,idx),scatthv_g(j,k,idx), &
                        scatthv_g(j,m,idx),scatthv_val)

      tsarg_hv = N0mg*scatthv_val

      CALL scatt_interp(sngl(lamrg),lam_list(i,4),lam_list(j,4),sngl(alfg),alf_list(k,4),alf_list(m,4),&
                        scattk_g(i,k,idx),scattk_g(i,m,idx),scattk_g(j,k,idx), &
                        scattk_g(j,m,idx),scattk_val)

      tfsarg = N0mg*scattk_val

    end if !end fmg 
  end if !dualpol opt 3 

  temph = 4*lambda4/(pi4*Kw2)*(tsar_h+tsas_h+tsah_h+tsag_h+tsars_h+tsarh_h+tsarg_h)
  tempv = 4*lambda4/(pi4*Kw2)*(tsar_v+tsas_v+tsah_v+tsag_v+tsars_v+tsarh_v+tsarg_v)
  refl_rsa%T_sum_ref_h = temph
  refl_rsa%T_sum_ref_v = tempv

  IF(attn_ON == 2) THEN
    tempAhh = 2*8.686*lambda*(Ar_h+As_h+Ag_h+Ah_h+Ars_h+Arg_h+Arh_h)*1.e-3  ! 2-way
    tempAvv = 2*8.686*lambda*(Ar_v+As_v+Ag_v+Ah_v+Ars_v+Arg_v+Arh_v)*1.e-3  ! 2-way
    refl_rsa%T_Ahh = tempAhh
    refl_rsa%T_Avv = tempAvv
  ENDIF

  temphv = 4*lambda4/(pi4*Kw2)*ABS(tsar_hv+tsas_hv+tsah_hv+tsag_hv+tsars_hv+tsarh_hv+tsarg_hv) 
  refl_rsa%T_sum_ref_hv = temphv

  if(temph > 0.) refl_rsa%T_log_ref = 10*log10(MAX(1.0,temph))

  if(tempv > 0.) then 
    refl_rsa%T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
  endif

!JYS  if(tempk < 0.) tempk = 0.0
  tempk = 180.*lambda/pi*(tfsar+tfsas+tfsah+tfsag+tfsars+tfsarh+tfsarg)*1.e-3
  refl_rsa%T_kdp = tempk

! JYS check this part!!! Activate this part after verification!!!
!  IF(attn_ON == 2) THEN
!    refl_rsa%T_Adp = tempAhh - tempAvv
!  ENDIF

END FUNCTION refl_rsa

END MODULE DUALPARA
