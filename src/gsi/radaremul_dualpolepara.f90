! Adopted from RADREMUL library of ARPS system
! Used for direct reflectity DA capability
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module dualpolepara                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE DUALPOLEPARA

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

END MODULE DUALPOLEPARA
