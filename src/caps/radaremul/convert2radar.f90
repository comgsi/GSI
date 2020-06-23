!########################################################################
!#########                                                      #########
!#########               SUBROUTINE ZhhFromDualPolCG            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE ZhhFromDualPolCG(nx,ny,nz,rho,t,qscalar,rff,ref_h, &
                            dualpol,rdrwave,alpha)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine is same as ZhhFromDualPol except for
! no truncation of DSD for raindrop at 8 mm.
! CG stands for Complete Gamma function.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2005
!
! MODIFICATION HISTORY:
!
! 12/7/2010 Bryan Putnam
!  Added code to comply with new option to turn off hail.
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
!External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) !air tempreature (K)
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: alpha(nx,ny,nz,6)  ! shape parameter

  REAL, INTENT(OUT) :: rff(nx,ny,nz) ! Reflectivity (dBZ)
  REAL, INTENT(OUT) :: ref_h(nx,ny,nz) ! Reflectivity (mm**6/m**3)

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: rdrwave

  INTEGER :: i,j,k,dualpol

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lambda = rdrwave

  SELECT CASE (mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116)                      
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()

  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF

  rff = 0.0
  ref_h = 0.0

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------


  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

          ta = t(i,j,k)

          IF(mphyopt == 11) THEN
             alphar = alpha(i,j,k,P_QR)
             alphas = alpha(i,j,k,P_QS)
             alphah = alpha(i,j,k,P_QH)
             alphag = alpha(i,j,k,P_QG)
          END IF

          obs_dual = init_Refl()
          var_dsd = init_para_dsd()

          CALL rdr_obs(rho(i,j,k),qscalar(i,j,k,:), &
                           0.,obs_dual,var_dsd,1,dualpol)

          rff(i,j,k) = obs_dual%T_log_ref
          ref_h(i,j,k) = obs_dual%T_sum_ref_h


      END DO ! DO i
    END DO ! DO j
  END DO ! DO  k

  RETURN
END SUBROUTINE ZhhFromDualPolCG


!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE ZdpFromDualPolCG            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE ZdpFromDualPolCG(nx,ny,nz,rho,t,qscalar,Zdp,dualpol, &
                            rdrwave,alpha)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates differential reflectivity Zdp,
! which is obtained from the difference btw Zhh and Zvv.
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/6/2004
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) !air temperature (K)
  REAL :: alpha(nx,ny,nz,6)  ! shape parameter
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL, INTENT(OUT) :: Zdp(nx,ny,nz) ! Differential reflectivity

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: rdrwave
  INTEGER :: i,j,k,dualpol

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  REAL :: hcomp, vcomp

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lambda = rdrwave

  SELECT CASE (mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116)
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()

  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF

  Zdp = 0.

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------
  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1
          hcomp = 0.
          vcomp = 0.

          ta = t(i,j,k)

          IF(mphyopt == 11) THEN
            alphar = alpha(i,j,k,P_QR)
            alphas = alpha(i,j,k,P_QS)
            alphah = alpha(i,j,k,P_QH)
            alphag = alpha(i,j,k,P_QG)
          ENDIF

          obs_dual = init_Refl()
          var_dsd = init_para_dsd()

          CALL rdr_obs(rho(i,j,k),qscalar(i,j,k,:),  &
                          0.,obs_dual,var_dsd,2,dualpol)

          hcomp = obs_dual%T_sum_ref_h
          vcomp = obs_dual%T_sum_ref_v

          Zdp(i,j,k) = (hcomp - vcomp)**0.2

      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE ZdpFromDualPolCG

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE ZdrFromDualPolCG            #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE ZdrFromDualPolCG(nx,ny,nz,rho,t,qscalar,Zdr,dualpol, &
                            rdrwave,alpha)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates differential reflectivity Zdr,
! which is obtained from the ratio of Zhh and Zvv.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/6/2004
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) !air temperature (K)
  REAL :: alpha(nx,ny,nz,6)  ! shape parameter
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL, INTENT(OUT) :: Zdr(nx,ny,nz) ! Differential reflectivity

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: rdrwave

  INTEGER :: i,j,k,dualpol

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lambda = rdrwave

  SELECT CASE(mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116) 
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()


  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF


  Zdr = missing

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

          ta = t(i,j,k)

          IF(mphyopt == 11) THEN
             alphar = alpha(i,j,k,P_QR)
             alphas = alpha(i,j,k,P_QS)
             alphah = alpha(i,j,k,P_QH)
             alphag = alpha(i,j,k,P_QG)
          END IF

          obs_dual = init_Refl()
          var_dsd = init_para_dsd()

          CALL rdr_obs (rho(i,j,k),qscalar(i,j,k,:), &
                           0.,obs_dual,var_dsd,2,dualpol)

          Zdr(i,j,k) = obs_dual%T_log_zdr


      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE ZdrFromDualPolCG

!#######################################################################
!#######################################################################
!#########                                                      ########
!#########               SUBROUTINE rhvFromDualPolCG            ########
!#########                                                      ########
!#########                     Developed by                     ########
!#########     Center for Analysis and Prediction of Storms     ########
!#########                University of Oklahoma                ########
!#########                                                      ########
!#######################################################################
!#######################################################################

SUBROUTINE rhvFromDualPolCG(nx,ny,nz,rho,t,qscalar,rhv,dualpol, &
                            rdrwave,alpha)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates cross-correlation coefficient, rho_hv.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/6/2004
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) !air temperature (K)
  REAL :: alpha(nx,ny,nz,6)  ! shape parameter
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL, INTENT(OUT) :: rhv(nx,ny,nz) ! cross-correlation coefficient

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: rdrwave

  INTEGER :: i,j,k,dualpol


  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lambda = rdrwave

  SELECT CASE (mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116)
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()

  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF

  rhv = missing

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

          ta = t(i,j,k)

          IF(mphyopt == 11) THEN
              alphar = alpha(i,j,k,P_QR)
              alphas = alpha(i,j,k,P_QS)
              alphah = alpha(i,j,k,P_QH)
              alphag = alpha(i,j,k,P_QG)
          END IF

          obs_dual = init_Refl()
          var_dsd = init_para_dsd()

          CALL rdr_obs(rho(i,j,k),qscalar(i,j,k,:), &
                           0.,obs_dual,var_dsd,3,dualpol)

          IF(obs_dual%T_sum_ref_h*obs_dual%T_sum_ref_v > 0.) THEN
            rhv(i,j,k) = obs_dual%T_sum_ref_hv/                      &
                       SQRT(obs_dual%T_sum_ref_h*obs_dual%T_sum_ref_v)
          ELSE
            rhv(i,j,k) = missing
          ENDIF

      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE rhvFromDualPolCG

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for snow                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION snow_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting snow.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_a = (0.194 + 7.094*fw + 2.135*fw**2. - 5.225*fw**3.)*10.**(-4)

END FUNCTION snow_alpha_a

REAL FUNCTION snow_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting snow.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_b = (0.191 + 6.916*fw - 2.841*fw**2. - 1.160*fw**3.)*10.**(-4)

END FUNCTION snow_alpha_b


!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for hail                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION hail_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_a = (0.191 + 2.39*fw - 12.57*fw**2. + 38.71*fw**3.   &
                  - 65.53*fw**4. + 56.16*fw**5. - 18.98*fw**6.)*10.**(-3)

END FUNCTION hail_alpha_a

REAL FUNCTION hail_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_b = (0.165 + 1.72*fw - 9.92*fw**2. + 32.15*fw**3.        &
                  - 56.0*fw**4. + 48.83*fw**5. - 16.69*fw**6.)*10.**(-3)

END FUNCTION hail_alpha_b
!
!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for graupel               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
!

REAL FUNCTION grpl_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting graupel.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 3/10/2010
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  grpl_alpha_a = (0.081 + 2.04*fw - 7.39*fw**2. + 18.14*fw**3.   &
                  - 26.02*fw**4. + 19.37*fw**5. - 5.75*fw**6.)*10.**(-3)

END FUNCTION grpl_alpha_a

REAL FUNCTION grpl_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting graupel.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 3/10/2010
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  grpl_alpha_b = (0.076 + 1.74*fw - 7.52*fw**2. + 20.22*fw**3.        &
                  - 30.42*fw**4. + 23.31*fw**5. - 7.06*fw**6.)*10.**(-3)

END FUNCTION grpl_alpha_b

SUBROUTINE partialRefRain(N0,alpha,alp_a,alp_b,beta_a,beta_b,lamda,    &
                           refRainHH,refRainVV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial reflectivity for rain species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: N0,alpha,alp_a,alp_b
  REAL,INTENT(IN) :: lamda,beta_a,beta_b
  REAL,INTENT(OUT) :: refRainHH,refRainVV

  !local variables
  REAL :: expon_h
  REAL :: gamma_h
  REAL :: expon_v
  REAL :: gamma_v
  REAL :: N0_units

  REAL*8 :: gamma

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0_units = (1.e-3)**(4.0+alpha)

  gamma_h = sngl(gamma(dble(alpha)+2.d0*dble(beta_a)+1.d0))
  expon_h = -(alpha+2*beta_a+1)
  gamma_v = sngl(gamma(dble(alpha)+2.d0*dble(beta_b)+1.d0))
  expon_v = -(alpha+2*beta_b+1)

   refRainHH = mm3todBZ*radar_const*alp_a**2*(N0*N0_units)*gamma_h* &
               (lamda*1.e-3)**expon_h

   refRainVV = mm3todBZ*radar_const*alp_b**2*(N0*N0_units)*gamma_v* &
               (lamda*1.e-3)**expon_v


END SUBROUTINE partialRefRain


SUBROUTINE partialRhoRain(N0,alpha,alp_a,alp_b,beta_a,beta_b,         &
                          lamda,refRainHV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the cross components, Z_hv, for rain species
! for rho_hv calculation.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: N0,alpha,alp_a,alp_b,beta_a,beta_b,lamda

  REAL,INTENT(OUT) :: refRainHV

  !local variables
  REAL :: expon_hv
  REAL :: gamma_hv
  REAL :: N0_units

  REAL*8 :: gamma

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0_units = (1.e-3)**(4.0+alpha)

  gamma_hv = sngl(gamma(dble(beta_a)+dble(beta_b)+dble(alpha)+1.d0))
  expon_hv = -(alpha+beta_a+beta_b+1)

  refRainHV = mm3todBZ*radar_const*alp_a*alp_b*(N0*N0_units)*           &
                gamma_hv*(lamda*1.e-3)**expon_hv

END SUBROUTINE partialRhoRain


SUBROUTINE partialRefIce(N0,alpha,Ai,Bi,Ci,alp_a,alp_b,beta_a,beta_b,   &
                         lamda,refIceHH,refIceVV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial reflectivity for each species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  USE DUALPARA

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: N0,alpha,Ai,Bi,Ci,alp_a,alp_b,beta_a,beta_b
  REAL,INTENT(IN) :: lamda

  REAL,INTENT(OUT) :: refIceHH,refIceVV

  !local variables
  REAL :: gamma_h, gamma_v, expon_h, expon_v
  REAL :: N0_units
  REAL*8 :: gamma

  gamma_h = sngl(gamma(dble(alpha) + 2.d0*dble(beta_a)+1.d0))
  expon_h = -(alpha+2*beta_a+1)
  gamma_v = sngl(gamma(dble(alpha) + 2.d0*dble(beta_b)+1.d0))
  expon_v = -(alpha + 2*beta_b+1)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0_units = (1.e-3)**(4.0+alpha)

  refIceHH = mm3toDBZ*radar_Const*gamma_h*(N0*N0_units)*                     &
              (Ai*alp_a**2+Bi*alp_b**2+2*Ci*alp_a*alp_b)*                 &
             (lamda*1.e-3)**expon_h

  refIceVV = mm3toDBZ*radar_Const*gamma_v*(N0*N0_units)*                    &
             (Bi*alp_a**2+Ai*alp_b**2+2*Ci*alp_a*alp_b)*                  &
             (lamda*1.e-3)**expon_v


END SUBROUTINE partialRefIce

SUBROUTINE partialRhoIce(N0,alpha,Ci,Di,alp_a,alp_b,beta_a,beta_b,rho_0,lamda,refIceHV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the cross components, Z_hv, for each species
! for rho_hv calculation.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/16/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE DUALPARA

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: N0,alpha,Ci,Di,alp_a,alp_b,beta_a,beta_b
  REAL,INTENT(IN) :: rho_0,lamda

  REAL,INTENT(OUT) :: refIceHV

  !local variables
   REAL :: gamma_hv, expon
   REAL :: N0_units
   REAL*8 :: gamma

   gamma_hv = sngl(gamma(dble(beta_a)+dble(beta_b)+dble(alpha) + 1.d0))
   expon = -(alpha + beta_a + beta_b + 1)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   N0_units = (1.e-3)**(4.0+alpha)

   refIceHV = mm3todBZ*radar_Const*gamma_hv*(N0*N0_units)*             &
               (Ci*alp_a**2+Ci*alp_b**2+2*Di*alp_a*alp_b*rho_0)*    &
               (lamda*1.e-3)**expon

END SUBROUTINE partialRhoIce


SUBROUTINE partialKdpRain(coef,N0,alpha,lamda,KdpRain)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial Kdp for rain species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

 USE DUALPARA

 IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: coef,N0,alpha,lamda
  REAL, INTENT(OUT) :: KdpRain

  !local variables
  REAL :: expon_at,expon_bt,expon_a,expon_b
  REAL :: gamma3_16,gamma2_97
  REAL :: N0_units
  REAL*8 :: gamma

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  expon_at = alpha+3.16+1
  expon_a = -expon_at
  expon_bt = alpha+2.97+1
  expon_b = -expon_bt
  gamma3_16 = sngl(gamma(dble(expon_at)))
  gamma2_97 = sngl(gamma(dble(expon_bt)))

  N0_units = (1.e-3)**(4.0+alpha)

  KdpRain = coef*(N0*N0_units)*((((lamda*1.e-3)**expon_a)*gamma3_16) - &
                 (((lamda*1.e-3)**expon_b)*gamma2_97))

END SUBROUTINE

SUBROUTINE partialKdpIce(iceCof,Ck,alp_k,beta_k,N0,alpha,lamda,  &
                         KdpIce)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial Kdp for each species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: iceCof,Ck,alp_k,beta_k,N0,alpha,lamda
  REAL,INTENT(OUT) :: KdpIce


  !local variables
  REAL :: expon
  REAL :: N0_units

  REAL*8 :: gamma
  REAL :: gamma_k

  expon = -(alpha+beta_k+1)

  gamma_k = sngl(gamma(dble(alpha)+beta_k+1.d0))

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   N0_units = (1.e-3)**(4.0 + alpha)

   KdpIce = iceCof*Ck*alp_k*gamma_k*(N0*N0_units)*(lamda*1.e-3)**expon

END SUBROUTINE partialKdpIce

SUBROUTINE fractionWater(qr,qi,fo,density_ice,fracqr,fracqi,fm,fw,rhom)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, fo, density_ice
  REAL,INTENT(OUT) :: fracqr, fracqi, fm, fw, rhom

  REAL :: fr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fr = 0.
  fw = 0.
  fracqr = 0.
  fracqi = 0.
  fm = 0.
  rhom = 0.

!-----------------------------------------------------------------------
! Calculate the fraction of mleting ice (fr) based on the ratio between
! qr and qi. fo is the maximum allowable fraction of melting snow.
!-----------------------------------------------------------------------
  IF (qr > 0. .AND. qi > 0.) THEN
    fr = fo*(MIN(qi/qr,qr/qi))**.3
  ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------
  fracqr = fr * qr
  fracqi = fr * qi
  fm = fracqr + fracqi

  IF (fm .EQ. 0. .AND. qr > 0.) THEN
    fw = 1.
  ELSE IF (fm > 0.) THEN
    fw = fracqr/fm
  ENDIF

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater

SUBROUTINE fractionWater_temperature_snow (qr,qi,density_ice,fm,fw,rhom,tair_C)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture based on the air temperature.
! It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 7/25/2014
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, density_ice, tair_C
  REAL,INTENT(OUT) :: fm, fw, rhom
  REAL :: layer_tmax, layer_tmin

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fw = 0.
  fm = 0.
  rhom = 0.

!-----------------------------------------------------------------------
! Calculate the faction of water.
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  fm = qi

! Compute the degree of wet in percentage based on air temperature
  layer_tmax = 2.5
  layer_tmin = -2.5
  if(tair_C >= layer_tmin .and. tair_C < layer_tmax) then
    fw = (tair_C - layer_tmin)/(layer_tmax-layer_tmin)
  else if(tair_C >= layer_tmax) then
    fw = 1.
  else
    fm = 0.
    fw = 0.
  endif

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater_temperature_snow

SUBROUTINE fractionWater_temperature_hail(qr,qi,density_ice,fm,fw,rhom,tair_C)
  
!-----------------------------------------------------------------------
! 
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture based on the air temperature.
! It also calculate the density of mixture.
! 
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 7/25/2014
! 
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, density_ice, tair_C
  REAL,INTENT(OUT) :: fm, fw, rhom
  REAL :: layer_tmax, layer_tmin
  REAL :: maxfrac

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fw = 0.
  fm = 0.
  rhom = 0.
  maxfrac = 0.6

!-----------------------------------------------------------------------
! Calculate the faction of water.
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------

  fm = qi

! Compute the degree of wet in percentage based on air temperature
  layer_tmax = 5.0
  layer_tmin = 0.0
  if(tair_C >= layer_tmin .and. tair_C < layer_tmax) then
    fw = (tair_C - layer_tmin)/(layer_tmax-layer_tmin) * maxfrac
  else if(tair_C >= layer_tmax) then
    fw = maxfrac
  else
    fm = 0.
    fw = 0.
  endif

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater_temperature_hail

!########################################################################
!########################################################################
!#########                                                      #########
!#########                    SUBROUTINE Kdp                    #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE Kdp(nx,ny,nz,rho,t,qscalar,kdph,dualpol,rdrwave,alpha)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates specific differential phase.
!
! The equations are as follows:
!
! Kdp = Kdpr + Kdps (contributions from rain and snow).
!
!            180 * lambda  * alphak  * No
!   Kdpr = ------------------------------- * gamma(betak + 1)
!                           (betak + 1)
!                   pi * gam
!
!           where, alphak = 1.33e-5,   betak = 4.61
!
!
!                                 / pi * No * rho  \1/4
!           slop parameter gam = | --------------- |
!                                 \      w        /
!              No, rho, and w depend on the state of particle
!              w = particle contents ( = air density * mixing ratio)
!
!           gamma : complete gamma function
!
!-----------------------------------------------------------------------
!
! REFERENCES:
!
! Zhang, G., J. Vivekanandan, and E. Brandes, 2001: A method for
!   estimating rain rate and drop size distribution from polarnx-2tric
!   radar measurements. IEEE Trans. Geosci. Remote Sens., 39,
!   830-841
!
! Ryzhkov, A. V., D. S. Zrnic, and B. A. Gordon, 1998: Polarnx-2tric
!   method for ice water content determination. J. Appl. Meteor., 37,
!   125-134
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/25/2005
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE global_paraest
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: calculate_kdp
  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid


  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) !Air Temperature
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: alpha(nx,ny,nz,6)  ! shape parameter
  REAL, INTENT(OUT) :: kdph(nx,ny,nz) ! Specific differential phase

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: rdrwave

  INTEGER :: i,j,k,dualpol

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (dsdparaopt)
  CASE (0)
    CALL init_dsd()
  CASE (1)

    CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl,  &
                   alpharain,alphasnow,alphagrpl,alphahail)
  END SELECT

  lambda = rdrwave

  SELECT CASE (mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116)
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()

  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF

!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

  kdph = 0.0

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------

  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

        ta = t(i,j,k)

        IF(mphyopt == 11) THEN
          alphar = alpha(i,j,k,P_QR)
          alphas = alpha(i,j,k,P_QS)
          alphah = alpha(i,j,k,P_QH)
          alphag = alpha(i,j,k,P_QG)
        END IF

        obs_dual = init_Refl()
        var_dsd = init_para_dsd()

        CALL rdr_obs(rho(i,j,k),qscalar(i,j,k,:),  &
                         kdph(i,j,k),obs_dual,var_dsd,4,dualpol)

        if(dualpol == 2) then
           kdph(i,j,k) = obs_dual%T_kdp
        endif


      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE Kdp


!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE ZdrFromDualPol              #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE ZdrFromDualPol(nx,ny,nz,rho,t,qscalar,rff,sumhcomp,&
                          sumvcomp,sumhvcomp,dualpol,alpha)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates differential reflectivity Zdr,
! which is obtained from the ratio of Zhh and Zvv.
! Equation for Zvv is same as Zhh with different alpha and beta.
! (alpha = 4.76e-4,   beta = 2.69)
!
! The equations are as follows:
!
!                 / Zeh_rain + Zeh_snow + Zeh_hail \
! Zdr = 10*LOG10 | -------------------------------- |
!                 \ Zev_rain + Zeh_snow + Zeh_hail /
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/6/2004
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE global_paraest
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nz ! Dimensions of grid

  REAL, INTENT(IN) :: rho(nx,ny,nz) ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t(nx,ny,nz) ! Air temperature (K)
  REAL :: qscalar(nx,ny,nz,nscalar)
  REAL :: alpha(nx,ny,nz,6)

  REAL, INTENT(OUT) :: rff(nx,ny,nz)  ! horizontal reflectivity
  REAL, INTENT(OUT) :: sumhcomp(nx,ny,nz)
  REAL, INTENT(OUT) :: sumvcomp(nx,ny,nz)
  REAL, INTENT(OUT) :: sumhvcomp(nx,ny,nz)

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  INTEGER :: i,j,k,dualpol
  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (dsdparaopt)
  CASE (0)
    CALL init_dsd()
  CASE (1)
    CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl,   &
                  alpharain,alphasnow,alphagrpl,alphahail)
  END SELECT

  lambda = wavelen

  SELECT CASE (mphyopt)
  CASE(1:4)
    graupel_ON = 0
    hail_ON    = 1
  CASE(5:7,106,108,110,116)
    graupel_ON = 1
    hail_ON    = 0
  END SELECT

  qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  CALL set_dsd_para()

  CALL calcMDR()

  IF(firstcalled) THEN
    CALL calcConstants()
    firstcalled = .false.
  ENDIF

!-----------------------------------------------------------------------
! Initilization of arrays
!-----------------------------------------------------------------------
  sumhcomp = 0.0
  sumvcomp = 0.0
  rff = 0.0

!-----------------------------------------------------------------------
! Now loop through the scalar grid points.
!-----------------------------------------------------------------------


  DO k = 1,nz-1
    DO j = 1,ny-1
      DO i = 1,nx-1

        ta = t(i,j,k)

        IF(mphyopt == 11) THEN
          alphar = alpha(i,j,k,P_QR)
          alphas = alpha(i,j,k,P_QS)
          alphah = alpha(i,j,k,P_QH)
          alphag = alpha(i,j,k,P_QG)
       END IF


        obs_dual = init_Refl()
        var_dsd = init_para_dsd()

        CALL rdr_obs (rho(i,j,k),qscalar(i,j,k,:), &
                           0.,obs_dual,var_dsd,3,dualpol)

        sumhcomp(i,j,k) = obs_dual%T_sum_ref_h
        sumvcomp(i,j,k) = obs_dual%T_sum_ref_v
        sumhvcomp(i,j,k) = obs_dual%T_sum_ref_hv
        rff(i,j,k) = obs_dual%T_log_ref

      END DO ! DO i
    END DO ! DO j
  END DO ! DO k

  RETURN
END SUBROUTINE ZdrFromDualPol

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION calculate_kdp                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################
 REAL FUNCTION calculate_kdp(rho,var_dsd)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates specific differential phase.
! For details, see "Kdp" in convertZ.f90.
!
!
! AUTHOR:  Youngsun Jung, 2/28/2007
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho ! Air density (kg m**-3)
  TYPE(T_para_dsd) :: var_dsd

  REAL :: qr
  REAL :: qs
  REAL :: qh
  REAL :: qg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  calculate_kdp = 0.0

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------

  qr = var_dsd%T_qr
  qs = var_dsd%T_qs
  qh = var_dsd%T_qh
  qg = var_dsd%T_qg

  IF (rho > 0.0 .and. (qs > 0.0 .or. qr > 0.0 .or. qh > 0.0 .or.    &
      qg > 0.0)) THEN

     calculate_kdp = rainIceKdp(var_dsd,rho)

  END IF                             ! outer if ends

  RETURN
END FUNCTION calculate_kdp

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE q2Zdr1                      #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE q2Zdr1(rho,t,qscalar,rff,sumhcomp,sumvcomp,sumhvcomp,l)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates differential reflectivity Zdr,
! which is obtained from the ratio of Zhh and Zvv.
! See SUBROUTINE ZdrFromDualPol
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 3/17/2005
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE global_paraest
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho       ! Air density (kg m**-3)
  REAL, INTENT (IN) :: t     !air temperature (K)
  REAL :: qscalar(nscalar)

  REAL, INTENT(OUT) :: sumhcomp      ! Refl. at horizontal polarization
  REAL, INTENT(OUT) :: sumvcomp      ! Refl. at vertical polarizationa
  REAL, INTENT(OUT) :: sumhvcomp     ! Zhv
  REAL, INTENT(OUT) :: rff           ! Refl. in dBZ

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd
  INTEGER :: l,nq
  REAL    :: Ntr,Nts,Nth,Ntg
  REAL    :: qrp,qsp,qhp,qgp

  REAL*8  :: solveAlpha

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if(dsdparaopt == 0)then
    CALL init_dsd()
  else if(dsdparaopt == 1)then


    CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl, &
                   alpharain,alphasnow,alphagrpl,alphahail)
  else
    WRITE(6,*) 'ERROR: dsdparaopt is not specified correctly'
    rff = -999.99
    RETURN
    !STOP
  endif

  IF(firstcalled) THEN
    lambda = wavelen
    SELECT CASE (mphyopt)
    CASE(1:4)
      graupel_ON = 0
      hail_ON    = 1
    CASE(5:7,106,108,110,116)
      graupel_ON = 1
      hail_ON    = 0
    END SELECT

    qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

    grpl_ON = graupel_ON
    hl_ON = hail_ON

    CALL set_dsd_para()

    CALL calcConstants()

    CALL calcMDR()

    firstcalled = .false.
  ENDIF

  rff = 0.0
  sumhcomp = 0.0
  sumvcomp = 0.0
  sumhvcomp = 0.0
  obs_dual = init_Refl()
  var_dsd = init_para_dsd()
  Ntr = 0.0; Nts = 0.0; Nth = 0.0; Ntg = 0.0
  qrp = 0.0; qsp = 0.0; qhp = 0.0; qgp = 0.0
  alpha_dsd = 0.0

!-----------------------------------------------------------------------
! Call observation calculation subroutine
!-----------------------------------------------------------------------

  IF(rfopt == 1 .and. mphyopt >= 9) THEN

!JYS    CALL rev_pow( Ntr, MAX(0.0,qscalar(P_NR)) )
!JYS    CALL rev_pow( Nts, MAX(0.0,qscalar(P_NS)) )
!JYS    CALL rev_pow( Nth, MAX(0.0,qscalar(P_NH)) )
    Ntr = MAX( 0.0,qscalar(P_NR) )
    Nts = MAX( 0.0,qscalar(P_NS) )
    Nth = MAX( 0.0,qscalar(P_NH) )
    Ntg = MAX( 0.0,qscalar(P_NG) )

    qrp = MAX( 0.0,qscalar(P_QR) )
    qsp = MAX( 0.0,qscalar(P_QS) )
    qhp = MAX( 0.0,qscalar(P_QH) )
    qgp = MAX( 0.0,qscalar(P_QG) )


    qscalar(P_QR) = qrp
    qscalar(P_QS) = qsp
    qscalar(P_QH) = qhp
    qscalar(P_QG) = qgp
    qscalar(P_NR) = Ntr
    qscalar(P_NS) = Nts
    qscalar(P_NH) = Nth
    qscalar(P_NG) = Ntg


  ENDIF

  IF(rfopt == 2 .OR. rfopt == 3) THEN
    SELECT CASE (mphyopt)
    CASE (2:4)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphahail
      alpha_dsd(6) = alphahail
    CASE (5:10)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphagrpl
      alpha_dsd(6) = alphahail
    CASE (11)
        alpha_dsd(2) = solveAlpha(qscalar(2),qscalar(8),      &
                               qscalar(13),c_x(2),rho)
        alpha_dsd(4) = solveAlpha(qscalar(4),qscalar(10),     &
                               qscalar(15),c_x(4),rho)
        IF(graupel_ON == 1) THEN
          alpha_dsd(5) = solveAlpha(qscalar(5),qscalar(11),   &
                                 qscalar(16),c_x(5),rho)
        ENDIF
        IF(hail_ON == 1) THEN
          alpha_dsd(6) = solveAlpha(qscalar(6),qscalar(12),   &
                                 qscalar(17),c_x(6),rho)
        ENDIF
    END SELECT
  ENDIF

  ta = t

  alphar = alpha_dsd(2)
  alphas = alpha_dsd(4)
  alphag = alpha_dsd(5)
  alphah = alpha_dsd(6)

  CALL rdr_obs(rho,qscalar,0.,obs_dual,var_dsd,3,rfopt)

  sumhcomp  = obs_dual%T_sum_ref_h
  sumvcomp  = obs_dual%T_sum_ref_v
  sumhvcomp = obs_dual%T_sum_ref_hv

  rff = obs_dual%T_log_ref

  RETURN
END SUBROUTINE q2Zdr1

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE q2Z1                        #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE q2Z1(rho,t,qscalar,rff,l)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine is same as calculate Z at some point in the space.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2005
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE global_paraest
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------
  REAL, INTENT(IN) :: rho          ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t            !air temperature (K)
  REAL :: qscalar(nscalar)

  REAL, INTENT(OUT) :: rff         ! Reflectivity (dBZ)

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd
  INTEGER :: l,nq
  REAL    :: Ntr,Nts,Nth,Ntg
  REAL    :: qrp,qsp,qhp,qgp   ! positive definied qr, qs, qh etc.

  REAL*8  :: solveAlpha

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if(dsdparaopt == 0)then
    CALL init_dsd()
  else if(dsdparaopt == 1)then

    CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl, &
                   alpharain, alphasnow,alphagrpl,alphahail)
  else
    WRITE(6,'(1x,a)') 'dsdparaopt is not specified correctly'
    rff = -999.9
    RETURN
!    STOP
  endif

  IF(firstcalled) THEN
    lambda = wavelen

    SELECT CASE (mphyopt)
    CASE(1:4)
      graupel_ON = 0
      hail_ON    = 1
    CASE(5:7,106,108,110,116)
      graupel_ON = 1
      hail_ON    = 0
    END SELECT

    qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

    grpl_ON = graupel_ON
    hl_ON = hail_ON

    CALL calcConstants()

    CALL set_dsd_para()

    CALL calcMDR()
  
    firstcalled = .false.
  ENDIF

  rff = 0.0
  obs_dual = init_Refl()
  var_dsd = init_para_dsd()
  Ntr = 0.0; Nts = 0.0; Nth = 0.0; Ntg = 0.0
  qrp = 0.0; qsp = 0.0; qhp = 0.0; qgp = 0.0
  alpha_dsd = 0.0

!-----------------------------------------------------------------------
! Calculate reflectivity
!-----------------------------------------------------------------------

  IF(rfopt == 1 .and. mphyopt >= 9) THEN

!JYS    CALL rev_pow( Ntr, MAX(0.0,qscalar(P_NR)) )
!JYS    CALL rev_pow( Nts, MAX(0.0,qscalar(P_NS)) )
!JYS    CALL rev_pow( Nth, MAX(0.0,qscalar(P_NH)) )
    IF(P_NR > 0) Ntr = MAX( 0.0,qscalar(P_NR) )
    IF(P_NS > 0) Nts = MAX( 0.0,qscalar(P_NS) )
    IF(P_NH > 0) Nth = MAX( 0.0,qscalar(P_NH) )
    IF(P_NG > 0) Ntg = MAX( 0.0,qscalar(P_NG) )

    IF(P_QR > 0) qrp = MAX( 0.0,qscalar(P_QR) )
    IF(P_QS > 0) qsp = MAX( 0.0,qscalar(P_QS) )
    IF(P_QH > 0) qhp = MAX( 0.0,qscalar(P_QH) )
    IF(P_QG > 0) qgp = MAX( 0.0,qscalar(P_QG) )


    IF(P_QR > 0) qscalar(P_QR) = qrp
    IF(P_QS > 0) qscalar(P_QS) = qsp
    IF(P_QH > 0) qscalar(P_QH) = qhp
    IF(P_QG > 0) qscalar(P_QG) = qgp
    IF(P_NR > 0) qscalar(P_NR) = Ntr
    IF(P_NS > 0) qscalar(P_NS) = Nts
    IF(P_NH > 0) qscalar(P_NH) = Nth
    IF(P_NG > 0) qscalar(P_NG) = Ntg

  END IF


  IF(rfopt == 2 .OR. rfopt == 3) THEN
    SELECT CASE (mphyopt)
    CASE (2:4)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphahail
      alpha_dsd(6) = alphahail
    CASE (5:10)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphagrpl
      alpha_dsd(6) = alphahail
    CASE (11)
        alpha_dsd(2) = solveAlpha(qscalar(2),qscalar(8),      &
                               qscalar(13),c_x(2),rho)
        alpha_dsd(4) = solveAlpha(qscalar(4),qscalar(10),     &
                               qscalar(15),c_x(4),rho)
        IF(graupel_ON == 1) THEN
          alpha_dsd(5) = solveAlpha(qscalar(5),qscalar(11),   &
                                 qscalar(16),c_x(5),rho)
        ENDIF
        IF(hail_ON == 1) THEN
          alpha_dsd(6) = solveAlpha(qscalar(6),qscalar(12),   &
                                 qscalar(17),c_x(6),rho)
        ENDIF
    END SELECT
  ENDIF

  ta = t

  alphar = alpha_dsd(2)
  alphas = alpha_dsd(4)
  alphag = alpha_dsd(5)
  alphah = alpha_dsd(6)

  CALL rdr_obs(rho,qscalar,0.,obs_dual,var_dsd,1,rfopt)

  rff = obs_dual%T_log_ref

  RETURN
END SUBROUTINE q2Z1

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE q2kdp1                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE q2kdp1(rho,t,qscalar,kdp,l)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates specific differential phase.
! For details, see "Kdp" in convertZ.f90.
!
!
! AUTHOR:  Youngsun Jung, 5/25/2005
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE global_paraest
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: calculate_kdp
  INTEGER, EXTERNAL :: get_qgh_opt

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho ! Air density (kg m**-3)
  REAL, INTENT(IN) :: t   ! Air temperature (K)
  REAL :: qscalar(nscalar)

  REAL, INTENT(OUT) :: kdp ! Specific differential phase

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd

  REAL    :: rcomp,scomp,hcomp,sumcomp
  INTEGER :: l,nq
  REAL    :: Ntr,Nts,Nth,Ntg
  REAL    :: qrp,qsp,qhp,qgp   ! positive definied qr, qs, qh etc.

  REAL*8  :: solveAlpha

  LOGICAL :: firstcalled
  SAVE firstcalled
  DATA firstcalled/.true./

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  if(dsdparaopt == 0)then
    CALL init_dsd()
  else if(dsdparaopt == 1)then

    CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl, &
                    alpharain,alphasnow,alphagrpl,alphahail)
  else
    write(6,*) 'ERROR: dsdparaopt is not specified correctly'
    kdp = -999.99
    RETURN
    !STOP
  endif

  IF(firstcalled) THEN
    lambda = wavelen

    SELECT CASE (mphyopt)
    CASE(1:4)
      graupel_ON = 0
      hail_ON    = 1
    CASE(5:7,106,108,110,116)
      graupel_ON = 1
      hail_ON    = 0
    END SELECT

    qgh_opt = get_qgh_opt(graupel_ON,hail_ON)

    grpl_ON = graupel_ON
    hl_ON = hail_ON

    CALL set_dsd_para()

    CALL calcConstants()

    CALL calcMDR()

    firstcalled = .false.
  ENDIF

  kdp = 0.0
  obs_dual = init_Refl()
  var_dsd = init_para_dsd()
  Ntr = 0.0; Nts = 0.0; Nth = 0.0; Ntg = 0.0
  qrp = 0.0; qsp = 0.0; qhp = 0.0; qgp = 0.0
  alpha_dsd = 0.0

  IF(rfopt == 1 .and. mphyopt >= 9) THEN
!JYS    CALL rev_pow( Ntr,qscalar(P_NR) )
!JYS    CALL rev_pow( Nts,qscalar(P_NS) )
!JYS    CALL rev_pow( Nth,qscalar(P_NH) )
    Ntr = MAX( 0.0,qscalar(P_NR) )
    Nts = MAX( 0.0,qscalar(P_NS) )
    Nth = MAX( 0.0,qscalar(P_NH) )
    Ntg = MAX( 0.0,qscalar(P_NG) )

    qrp = MAX( 0.0,qscalar(P_QR) )
    qsp = MAX( 0.0,qscalar(P_QS) )
    qhp = MAX( 0.0,qscalar(P_QH) )
    qgp = MAX( 0.0,qscalar(P_QG) )

    qscalar(P_QR) = qrp
    qscalar(P_QS) = qsp
    qscalar(P_QH) = qhp
    qscalar(P_QG) = qgp
    qscalar(P_NR) = Ntr
    qscalar(P_NS) = Nts
    qscalar(P_NH) = Nth
    qscalar(P_NG) = Ntg

  ENDIF

  IF(rfopt == 2 .OR. rfopt == 3) THEN
    SELECT CASE (mphyopt)
    CASE (2:4)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphahail
      alpha_dsd(6) = alphahail
    CASE (5:10)
      alpha_dsd(2) = alpharain
      alpha_dsd(3) = alphaice
      alpha_dsd(4) = alphasnow
      alpha_dsd(5) = alphagrpl
      alpha_dsd(6) = alphahail
    CASE (11)
        alpha_dsd(2) = solveAlpha(qscalar(2),qscalar(8),      &
                               qscalar(13),c_x(2),rho)
        alpha_dsd(4) = solveAlpha(qscalar(4),qscalar(10),     &
                               qscalar(15),c_x(4),rho)
        IF(graupel_ON == 1) THEN
          alpha_dsd(5) = solveAlpha(qscalar(5),qscalar(11),   &
                                 qscalar(16),c_x(5),rho)
        ENDIF
        IF(hail_ON == 1) THEN
          alpha_dsd(6) = solveAlpha(qscalar(6),qscalar(12),   &
                                 qscalar(17),c_x(6),rho)
        ENDIF
    END SELECT
  ENDIF

  ta = t

  alphar = alpha_dsd(2)
  alphas = alpha_dsd(4)
  alphag = alpha_dsd(5)
  alphah = alpha_dsd(6)

  CALL rdr_obs(rho,qscalar,kdp,obs_dual,var_dsd,4,rfopt)

  if(rfopt == 2 .OR. rfopt == 3) then
    kdp = obs_dual%T_kdp
  endif

  RETURN
END SUBROUTINE q2kdp1


SUBROUTINE calc_moment(power,lamda,N0,alpha,moment)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate a given psd moment
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!


  INTEGER, INTENT(IN) :: power

  REAL, INTENT(IN) :: lamda,N0,alpha

  REAL, INTENT(OUT) :: moment

  !local variables
  REAL :: gamma_moment,expon
  REAL*8 :: gamma

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  gamma_moment = sngl(gamma(dble(alpha)+dble(power)+1.d0))
  expon = (alpha+power+1)

  moment = N0*(lamda**expon)*gamma_moment


END SUBROUTINE



SUBROUTINE power_mom(power,cx,t,rhoa,q,moment)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculates moments of the PSD based on the Field et al. 2005 power law
! relations. Used for Thompson scheme.
!
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst, only: mphyopt

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rhoa
  INTEGER, INTENT(IN) :: power
  REAL, INTENT(IN) :: t,q,cx
  REAL, INTENT(OUT) :: moment

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: a,b
  REAL :: rpower  
  REAL*8 :: log_a
  REAL :: second_moment,test
  REAL :: T_c

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  T_c = t-273.16

  SELECT CASE (mphyopt)
  CASE(108)

  second_moment = rhoa * (q/cx)

  IF(power == 2) THEN
    moment = second_moment
  ELSE
 
     rpower = REAL(power)

     log_a = dble(5.065339-.062659*T_c - 3.032362*rpower +                 &
                   0.029469*T_c*rpower -  &
     0.000285*(T_c**2.) + 0.312550*(rpower**2.) + 0.000204*(T_c**2.)*rpower + &
     0.003199*T_c*(rpower**2.) + 0.000000*(T_c**3.) - 0.015952*(rpower**3.))

     a = sngl(10.d0**log_a)

     b = 0.476221 - 0.015896*T_c + 0.165977*rpower + 0.007468*T_c*rpower -   &
      0.000141*(T_c**2.) + 0.060366*(rpower**2.) + 0.000079*(T_c**2.)*rpower + &
      0.000594*T_c*(rpower**2.) + 0.000000*(T_c**3.) - 0.003577*(rpower**3.)


    moment = a*(second_moment)**b
  END IF

  END SELECT

END SUBROUTINE


SUBROUTINE calc_N0x_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,  &
                       qsf,fms,qhf,fmh,qgf,fmg)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates intercep parameter based on MP scheme.
!
!
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst, only: mphyopt

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL :: rhoa,rhoms,rhomh,rhomg
  REAL :: ntr,nts,nth,ntg
  REAL :: qrf,qsf,qhf,qgf
  REAL :: fms,fmh,fmg

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------
  REAL :: moma,momb
  REAL :: no_value = missing

  REAL, PARAMETER :: D0r = 50.e-5
  REAL, PARAMETER :: R1 = 1.e-12
  REAL, PARAMETER :: R2 = 1.e-6
  REAL, PARAMETER :: gonv_min = 1.e4
  REAL, PARAMETER :: gonv_max = 3.e6 
  REAL, PARAMETER :: bm_g = 3.0 

  LOGICAL :: L_qr
  REAL :: mvd_r  
  REAL*8 :: dble_alfr
  REAL*8 :: lamr 
  REAL*8 :: gamma  
  REAL :: xslwq,ygra1,zans1
  REAL :: N0_exp,N0_min
  REAL :: rg,am_g,oge1,cgg_1,cgg_2,cgg_3,ogg1,ogg2,ogmg,cge_1
  REAL :: lam_exp,lamg,ilamg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   SELECT CASE (mphyopt)
   CASE(9:12,109)
     CALL calc_N0x_melt(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,   &
                         qsf,fms,qhf,fmh,qgf,fmg)
   CASE(106)

    N0r = 8.0E06 

    N0g = 4.0E06
    N0mg = N0g

    N0s = 2.0E06*exp((.12*(273.16-ta)))
    N0ms = N0s 
   CASE(108) 

     CALL calc_N0x_melt(rhoa,no_value,no_value,no_value,ntr,no_value, &
                        no_value,no_value,qrf,no_value,no_value,      &
                        no_value,no_value,no_value,no_value)

     
     IF(qrf > R1) THEN
       L_qr = .true.  
       dble_alfr = dble(alphar)
       lamr = 0.0 
       CALL cal_lamda(rhoa,qrf,ntr,rhor,dble_alfr,lamr)
       mvd_r = (3.0 + alphar + 0.672)/sngl(lamr) 
         IF(mvd_r > 2.5e-3) THEN
           mvd_r = 2.5e-3
         ELSE IF(mvd_r < ((D0r)*(0.75))) THEN
           mvd_r =  D0r*0.75
         END IF 
     ELSE
       L_qr = .false. 
       qrf = 0.0
     END IF

     IF(qgf > R1) THEN
       rg = qgf * rhoa
     ELSE
       rg = R1
     END IF 

     IF((ta < 270.65) .and. L_qr .and. (mvd_r > 100.0e-6)) THEN
        xslwq = 4.01 + log10(mvd_r)
     ELSE
        xslwq = 0.01
     END IF

     N0_min = gonv_max 
     ygra1 = 4.31 + log10(max(5.e-5,rg))
     zans1 = 3.1 + (100.0/(300.0*xslwq*ygra1/(10.0/xslwq+1.0+0.25* &
             ygra1)+30.0+10.0*ygra1))         
     N0_exp = 10.0**zans1
     N0_exp = MAX(gonv_min,MIN(N0_exp,gonv_max))
     N0_min = MIN(N0_exp,N0_min)
     N0_exp = N0_min        
     am_g = c_x(5)
     oge1 = 1./(bm_g + 1.)
     cgg_1 = sngl(gamma(dble(bm_g) + 1.d0)) 
     cgg_2 = sngl(gamma(dble(alphag) + 1.d0))  
     cgg_3 = sngl(gamma(dble(bm_g) + dble(alphag) + 1.d0))
     ogg1 = 1./cgg_1
     ogg2 = 1./cgg_2
     ogmg = 1./bm_g 
     cge_1 = alphag + 1.0
     lam_exp = (N0_exp*am_g*cgg_1/rg)**oge1
     lamg = lam_exp*(cgg_3*ogg2*ogg1)**ogmg
     N0g = N0_exp/(cgg_2*lam_exp)*lamg**cge_1

     IF(fmg > R1) THEN
       rg = fmg * rhoa
     ELSE
       rg = R1
     END IF

     N0_min = gonv_max
     ygra1 = 4.31 + log10(max(5.e-5,rg))
     zans1 = 3.1 + (100.0/(300.0*xslwq*ygra1/(10.0/xslwq+1.0+0.25* &
             ygra1)+30.0+10.0*ygra1))
     N0_exp = 10.0**zans1
     N0_exp = MAX(gonv_min,MIN(N0_exp,gonv_max))
     N0_min = MIN(N0_exp,N0_min)
     N0_exp = N0_min
     am_g = c_x(5)
     oge1 = 1./(bm_g + 1.)
     cgg_1 = sngl(gamma(dble(bm_g) + 1.d0))
     cgg_2 = sngl(gamma(dble(alphag) + 1.d0))
     cgg_3 = sngl(gamma(dble(bm_g) + dble(alphag) + 1.d0))
     ogg1 = 1./cgg_1
     ogg2 = 1./cgg_2
     ogmg = 1./bm_g
     cge_1 = alphag + 1.0
     lam_exp = (N0_exp*am_g*cgg_1/rg)**oge1
     lamg = lam_exp*(cgg_3*ogg2*ogg1)**ogmg
     N0mg = N0_exp/(cgg_2*lam_exp)*lamg**cge_1

     IF(qsf >= 1.e-14) THEN  

       CALL  power_mom(2,c_x(4),ta,rhoa,qsf,moma) 
       CALL  power_mom(3,c_x(4),ta,rhoa,qsf,momb)

       N0s = sngl(((dble(moma)**4.d0)/(dble(momb)**3.d0))*dble(thom_k0))
       N0s2 = sngl(((dble(moma)**4.d0)/(dble(momb)**3.d0))*dble(thom_k1)*      &
                  ((dble(moma)/dble(momb))**dble(alphas2)))

     ELSE
       N0s = 3.0E06
       N0s2 = 3.0E06 
     END IF 

   CASE(110)
     CALL calc_N0x_melt(rhoa,rhoms,no_value,rhomg,ntr,nts,no_value,ntg,   &
                        qrf,qsf,fms,no_value,no_value,qgf,fmg)
   CASE(116)
     CALL calc_N0x_melt(rhoa,no_value,no_value,no_value,ntr,no_value,       &
                         no_value,no_value,qrf,no_value,no_value,no_value,   &
                         no_value,no_value,no_value)


    N0g = 4.0E06
    N0mg = N0g 
 
    N0s = 2.0E06*exp((.12*(273.16-ta)))
    N0ms = N0s 

   END SELECT 

END SUBROUTINE

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CALC_N0X                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE calc_N0x(rhoa,ntr,nts,nth,ntg,qr,qs,qh,qg,alfr,alfs,alfh,alfg)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    03/25/2008.
!
!-----------------------------------------------------------------------
!
  USE DUALPARA

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL   :: rhoa
  REAL   :: ntr,nts,nth,ntg,qr,qs,qh,qg
  REAL*8 :: alfr,alfs,alfh,alfg
  REAL*8 :: db_N0r, db_N0s, db_N0h, db_N0g
  REAL*8 :: db_alfr,db_alfs,db_alfh,db_alfg
  REAL   :: pow1,pow2
  REAL, PARAMETER :: epsQ  = 1.e-14
  REAL, PARAMETER :: epsN  = 1.e-3
  REAL, PARAMETER :: maxN0 = 4.e+37

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  db_alfr = dble(alfr); db_alfs = dble(alfs); db_alfh = dble(alfh); db_alfg = dble(alfg)


  IF(qr >= epsQ .AND. ntr >= epsN) THEN
     CALL cal_N0(rhoa,qr,ntr,rhor,db_alfr,db_N0r)
     N0r = MIN(maxN0,sngl(db_N0r))
  ELSE
     qr = 0.0; ntr = 0.0; N0r = 8.0E+06
  ENDIF

  IF(qs >= epsQ .AND. nts >= epsN) THEN
     CALL cal_N0(rhoa,qs,nts,rhos,db_alfs,db_N0s)
     N0s = MIN(maxN0,sngl(db_N0s))
  ELSE
     qs = 0.0; nts = 0.0; N0s = 3.0E+06
  ENDIF

  IF(qh >= epsQ .AND. nth >= epsN) THEN
     CALL cal_N0(rhoa,qh,nth,rhoh,db_alfh,db_N0h)
     N0h = MIN(maxN0,sngl(db_N0h))
  ELSE
     qh = 0.0; nth = 0.0; N0h = 4.0E+04
  ENDIF

  IF(qg >= epsQ .AND. ntg >= epsN) THEN
     CALL cal_N0(rhoa,qg,ntg,rhog,db_alfg,db_N0g)
     N0g = MIN(maxN0,sngl(db_N0g))
  ELSE
     qg = 0.0; ntg = 0.0; N0g = 4.0E+05
  ENDIF

END SUBROUTINE calc_N0x

SUBROUTINE calc_N0x_melt(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,qrf,   &
                         qsf,fms,qhf,fmh,qgf,fmg)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate intercept parameter including melting species
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Bryan Putnam
!    04/16/2013.
!
!-----------------------------------------------------------------------

  USE DUALPARA

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------

  REAL   :: rhoa,rhoms,rhomh,rhomg
  REAL   :: ntr,nts,nth,ntg,qrf,qsf,qhf,qgf,fms,fmh,fmg
  REAL*8   :: db_N0r, db_N0s, db_N0h, db_N0g
  REAL*8   :: db_alfr,db_alfs,db_alfh,db_alfg
  REAL   :: pow1,pow2
  REAL, PARAMETER :: epsQ  = 1.e-14
  REAL, PARAMETER :: epsN  = 1.e-3
  REAL, PARAMETER :: maxN0 = 4.e+37


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 db_alfr = dble(alphar); db_alfs = dble(alphas); db_alfh = dble(alphah);
 db_alfg = dble(alphag)


  IF(qrf >= epsQ .AND. ntr >= epsN) THEN
     CALL cal_N0(rhoa,qrf,ntr,rhor,db_alfr,db_N0r)
     N0r = MIN(maxN0,sngl(db_N0r))
  ELSE
     qrf = 0.0
  ENDIF

  IF(qsf >= epsQ .AND. nts >= epsN) THEN
     CALL cal_N0(rhoa,qsf,nts,rhos,db_alfs,db_N0s)
     N0s = MIN(maxN0,sngl(db_N0s))
  ELSE
     qsf = 0.0
  ENDIF

  IF(fms >= epsQ .AND. nts >= epsN) THEN
     CALL cal_N0(rhoa,fms,nts,rhoms,db_alfs,db_N0s)
     N0ms = MIN(maxN0,sngl(db_N0s))
  ELSE
     fms = 0.0
  ENDIF

  IF(qhf >= epsQ .AND. nth >= epsN) THEN
     CALL cal_N0(rhoa,qhf,nth,rhoh,db_alfh,db_N0h)
     N0h = MIN(maxN0,sngl(db_N0h))
  ELSE
     qhf = 0.0
  ENDIF

  IF(fmh >= epsQ .AND. nth >= epsN) THEN
     CALL cal_N0(rhoa,fmh,nth,rhomh,db_alfh,db_N0h)
     N0mh = MIN(maxN0,sngl(db_N0h))
  ELSE
     fmh = 0.0
  ENDIF

  IF(qgf >= epsQ .AND. ntg >= epsN) THEN
     CALL cal_N0(rhoa,qgf,ntg,rhog,db_alfg,db_N0g)
     N0g = MIN(maxN0,sngl(db_N0g))
  ELSE
     qgf = 0.0
  ENDIF

  IF(fmg >= epsQ .AND. ntg >= epsN) THEN
     CALL cal_N0(rhoa,fmg,ntg,rhomg,db_alfg,db_N0g)
     N0mg = MIN(maxN0,sngl(db_N0g))
  ELSE
     fmg = 0.0
  ENDIF

END SUBROUTINE calc_N0x_melt


!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE CALC_Dnx                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE calc_Dnx(rhoa,ntr,nts,nth,qr,qs,qh,Dnr,Dns,Dnh)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    03/25/2008.
!
!-----------------------------------------------------------------------
!
  USE DUALPARA

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL :: rhoa
  REAL :: ntr,nts,nth,qr,qs,qh,Dnr,Dns,Dnh
  REAL  , parameter :: epsQ  = 1.e-14
  real  , parameter :: epsN  = 1.e-3

  IF(qr >= epsQ .AND. ntr >= epsN) THEN
     N0r=ntr**(4./3.)*(pi*rhor/(qr*rhoa))**(1./3.)
     Dnr = ntr/N0r
  ELSE
     qr = 0.0; ntr = 0.0
  ENDIF
  IF(qs >= epsQ .AND. nts >= epsN) THEN
     N0s=nts**(4./3.)*(pi*rhos/(qs*rhoa))**(1./3.)
     Dns = nts/N0s
  ELSE
     qs = 0.0; nts = 0.0
  ENDIF
  IF(qh >= epsQ .AND. nth >= epsN) THEN
     N0h=nth**(4./3.)*(pi*rhoh/(qh*rhoa))**(1./3.)
     Dnh = nth/N0h
  ELSE
     qh = 0.0; nth = 0.0
  ENDIF

END SUBROUTINE calc_Dnx

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE TAKE_REV                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE take_rev(nx,ny,nz,xq,xen)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Bring the total number concentrations back to original domain
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    05/20/2008
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,i,j,k
  REAL :: xq(nx,ny,nz),xen(nx,ny,nz)

  do i=1,nx
    do j=1,ny
      do k=1,nz
        if(xq(i,j,k) >= 1.e-14 .and. xen(i,j,k) >= 0.0631) then
          CALL rev_pow(xen(i,j,k),xen(i,j,k))
        else
          xen(i,j,k)=0.
          xq(i,j,k)=0.
        endif
      enddo
    enddo
  enddo

END SUBROUTINE take_rev

SUBROUTINE rev_pow(Ntpow,Nt)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take some power of variable.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    06/05/2008
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL :: Ntpow, Nt

  Ntpow=Nt**(1./0.4)

END SUBROUTINE rev_pow

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE TAKE_POW                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                 Unxversity of Oklahoma               ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE take_pow(nx,ny,nz,xen)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take power of the total number concentrations for retrieval
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Youngsun Jung
!    05/20/2008
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz,i,j,k
  REAL :: xen(nx,ny,nz)

  do i=1,nx
    do j=1,ny
      do k=1,nz
        if(xen(i,j,k) > 1.e-3) then
          xen(i,j,k)=xen(i,j,k)**0.4
        else
          xen(i,j,k)=0.0
        endif
      enddo
    enddo
  enddo

END SUBROUTINE take_pow

FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gamma
  INTEGER  :: j
  DOUBLE PRECISION  :: ser,stp,tmp,x,y,cof(6)


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gamma=tmp+log(stp*ser/x)
  gamma= exp(gamma)

END FUNCTION gamma

SUBROUTINE cal_N0(rhoa,q,Ntx,rhox,alpha,N0)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  (03/26/2008)
!  Recast N0 as a double precision variable, and used double precision for
!  all intermediate calculations.  The calling subroutine should
!  also define it as double precision.  For situations with large alpha,
!  N0 can become very large, and loss of precision can result.
!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
!  With Jason Milbrandt's calculation of N0 just before evaporation in
!  the multi-moment code.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q,Ntx
  REAL*8 :: alpha,N0
  REAL :: rhox
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma

  DOUBLE PRECISION :: lamda

  gamma1 = gamma(1.d0+dble(alpha))
  gamma4 = gamma(4.d0+dble(alpha))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
        dble(q)))**(1.d0/3.d0)
  ELSE
    lamda = 0.d0
  END IF

  N0 = dble(Ntx)*lamda**(0.5d0*(1.d0+dble(alpha)))*                         &
              (1.d0/gamma1)*lamda**(0.5d0*(1.d0+dble(alpha)))

END SUBROUTINE cal_N0


SUBROUTINE calc_lamda_mp(rhoa,rhoms,rhomh,rhomg,ntr,nts,nth,ntg,  &
                             qrf,qsf,fms,qhf,fmh,qgf,fmg)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Calculate slope parameter for PSD based on MP scheme.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Bryan Putnam, 4/16/2013
!
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst, only: mphyopt

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------


  REAL :: rhoa,rhoms,rhomh,rhomg
  REAL :: ntr,nts,nth,ntg
  REAL :: qrf,qsf,fms,qhf,fmh,qgf,fmg

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------


  REAL*8 :: db_N0,dble_alfr,dble_alfs,dble_alfg,dble_alfh
  REAL*8 :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL :: Ntw,Ntd

  REAL :: tem1,tem2

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  dble_alfr = dble(alphar)
  dble_alfs = dble(alphas)
  dble_alfg = dble(alphag)
  dble_alfh = dble(alphah)

  if(qrf > 0.0) then

   Ntw = 0.
   if(ntr > 0.0) then
    Ntw = ntr
    CALL cal_lamda(rhoa,qrf,Ntw,rhor,dble_alfr,lamr)
     lamdar = sngl(lamr)
   else
    db_N0 = dble(N0r)
    CALL cal_Nt(rhoa,qrf,db_N0,c_x(2),dble_alfr,Ntw)
    CALL cal_lamda(rhoa,qrf,Ntw,rhor,dble_alfr,lamr)
    lamdar = sngl(lamr)
   end if
  else
   lamdar = 0.0
  end if

  SELECT CASE (mphyopt)
  CASE(1:11,106,109,110,116)
   if(qsf > 0.0) then
    Ntd = 0.
    if (nts > 0.0) then
     Ntd = nts
     CALL cal_lamda(rhoa,qsf,Ntd,rhos,dble_alfs,lams)
     lamdas = sngl(lams)
    else
     db_N0 = dble(N0s)
     CALL cal_Nt(rhoa,qsf,db_N0,c_x(4),dble_alfs,Ntd)
     CALL cal_lamda(rhoa,qsf,Ntd,rhos,dble_alfs,lams)
     lamdas = sngl(lams)
    end if
   else
    lamdas = 0.0
   end if

    if(fms > 0.0) then
     Ntw = 0.
     if(nts > 0.0) then
      Ntw = nts
      CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs)
      lamdams = sngl(lamrs)
    else
     db_N0 = dble(N0s)
     CALL cal_Nt(rhoa,fms,db_N0,c_x(4),dble_alfs,ntw)
     CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs)
     lamdams = sngl(lamrs)
    end if
   else
    lamdams = 0.0
   end if

  CASE(108)
   if(qsf > 0.0) then

    CALL power_mom(2,c_x(4),ta,rhoa,qsf,tem1)
    CALL power_mom(3,c_x(4),ta,rhoa,qsf,tem2)
    lamdas = (tem1/tem2)*thom_lam0
    lamdas2  = (tem1/tem2)*thom_lam1
   else
    lamdas = 0.0
    lamdas2 = 0.0
   end if

   if(fms > 0.0) then
     Ntw = 0.
     if(nts > 0.0) then
      Ntw = nts
      CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs)
      lamdams = sngl(lamrs)
    else
     db_N0 = dble(N0ms)
     CALL cal_Nt(rhoa,fms,db_N0,c_x(4),dble_alfs,ntw)
     CALL cal_lamda(rhoa,fms,Ntw,rhoms,dble_alfs,lamrs)
     lamdams = sngl(lamrs)
    end if
   end if

  END SELECT

 if(hl_ON == 1) then
   if(qhf > 0.) then
    Ntd = 0.
    if(nth > 0.0) then
     Ntd = nth
     CALL cal_lamda(rhoa,qhf,Ntd,rhoh,dble_alfh,lamh)
     lamdah = sngl(lamh)
    else
     db_N0 = dble(N0h)
     CALL cal_Nt(rhoa,qhf,db_N0,c_x(6),dble_alfh,Ntd)
     CALL cal_lamda(rhoa,qhf,Ntd,rhoh,dble_alfh,lamh)
     lamdah = sngl(lamh)
    end if
   else
    lamdah = 0.0
   end if

   if(fmh > 0.) then
    Ntw = 0.
    if(nth > 0.0) then
     Ntw = nth
     CALL cal_lamda(rhoa,fmh,Ntw,rhomh,dble_alfh,lamrh)
     lamdamh = sngl(lamrh)
    else
     db_N0 = dble(N0mh)
     CALL cal_Nt(rhoa,fmh,db_N0,c_x(6),dble_alfh,Ntw)
     CALL cal_lamda(rhoa,fmh,Ntw,rhomh,dble_alfh,lamrh)
     lamdamh = sngl(lamrh)
    end if
   else
    lamdamh = 0.0
   end if
 end if

 if(grpl_ON == 1) then

   if(qgf > 0.) then
    Ntd = 0.
    if(ntg > 0.0) then
     Ntd = ntg
     CALL cal_lamda(rhoa,qgf,Ntd,rhog,dble_alfg,lamg)
     lamdag = sngl(lamg)
    else
     db_N0 = dble(N0g)
     CALL cal_Nt(rhoa,qgf,db_N0,c_x(5),dble_alfg,Ntd)
     CALL cal_lamda(rhoa,qgf,Ntd,rhog,dble_alfg,lamg)
     lamdag = sngl(lamg)
    end if
  else
   lamdag = 0.0
  end if

   if(fmg > 0.) then
    Ntw = 0.
    if(ntg > 0.0) then
     Ntw = ntg
     CALL cal_lamda(rhoa,fmg,Ntw,rhomg,dble_alfg,lamrg)
     lamdamg = sngl(lamrg)
    else
     db_N0 = dble(N0mg)
     CALL cal_Nt(rhoa,fmg,db_N0,c_x(5),dble_alfg,Ntw)
     CALL cal_lamda(rhoa,fmg,Ntw,rhomg,dble_alfg,lamrg)
     lamdamg = sngl(lamrg)
    end if
   else
    lamdamg = 0.0
   end if
  end if

END SUBROUTINE

SUBROUTINE cal_Nt(rhoa,q,N0,cx,alpha,Ntx)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration at scalar points
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  03/31/08 - converted intermediate calculations to double precision
!             as well as a few of the input arguments.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL :: rhoa,q
  REAL*8 :: alpha,N0
  REAL :: cx
  REAL :: Ntx
  REAL*8 :: gamma1,gamma4

  REAL*8 :: gamma

  gamma1 = gamma(1.d0+dble(alpha))
  gamma4 = gamma(4.d0+dble(alpha))

   Ntx = sngl((dble(N0)*gamma1)**(3.d0/(4.d0+dble(alpha)))*   &
             ((gamma1/gamma4)*dble(rhoa)* &
             dble(q)/dble(cx))**((1.d0+dble(alpha))/(4.d0+dble(alpha))))

END SUBROUTINE cal_Nt

SUBROUTINE cal_lamda(rhoa,q,Ntx,rhox,alpha,lamda)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter lamda
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Converted intermediate calculations and arrays alpha and lamda to
!  double precision.
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q
  REAL*8 :: alpha,lamda
  REAL :: rhox
  REAL :: Ntx
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma

  gamma1 = gamma(1.d0+dble(alpha))
  gamma4 = gamma(4.d0+dble(alpha))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = sngl(((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
          dble(q)))**(1.d0/3.d0))

  ELSE
    lamda = 0.d0
  END IF

END SUBROUTINE cal_lamda

SUBROUTINE solve_alpha(nx,ny,nz,rhoa,cx,q,Ntx,Z,alpha)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Changed alpha array to double precision
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  INTEGER :: nx,ny,nz
  REAL :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL :: alpha(nx,ny,nz)


  !Local Variables
  REAL*8 :: solveAlpha
  REAL*8 :: dsA

  REAL :: cx

  INTEGER i,j,k

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN

          alpha(i,j,k) = sngl(solveAlpha(q(i,j,k),Ntx(i,j,k),Z(i,j,k),cx,rhoa(i,j,k)))

        ELSE
          alpha(i,j,k) = 0.0
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE solve_alpha

FUNCTION solveAlpha(Q,N,Z,Cx,rho)

 IMPLICIT NONE

! PASSING PARAMETERS:
  real, INTENT(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real*8 :: solveAlpha
  real   :: a,g,a1,g1,g2,tmp1
  integer :: i
  real, parameter :: alphaMax= 40.
  real, parameter :: epsQ    = 1.e-14
  real, parameter :: epsN    = 1.e-3
  real, parameter :: epsZ    = 1.e-32

!  Q         mass mixing ratio
!  N         total concentration
!  Z         reflectivity
!  Cx        (pi/6)*RHOx
!  rho       air density
!  a         alpha (returned as solveAlpha)
!  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
!              where g = (Cx/(rho*Q))**2.*(Z*N)


!JYS  if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
!JYS  ! For testing/debugging only; this module should never be called
!JYS  ! if the above condition is true.
!JYS    print*,'*** STOPPED in MODULE ### solveAlpha *** '
!JYS    print*,'*** : ',Q,N,Z,Cx*1.9099,rho
!JYS    stop
!JYS  endif

  IF (Q>epsQ .and. N>epsN .and. Z>epsZ ) THEN

     tmp1= Cx/(rho*Q)
     g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2

 !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force; for testing only)
!      a= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            a = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of g(a) to solve for a:
     if (g>=20.) then
       a= 0.
     else
       g2= g*g
       if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
       if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
       if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
       if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
       if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
       if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
       if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
       if (g<1.230) a= alphaMax
     endif

     solveAlpha= max(0.,min(a,alphaMax))

  ELSE

     solveAlpha= 0.

  ENDIF

END FUNCTION solveAlpha

!
!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE RSET_DSD_PARA               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE set_dsd_para()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine sets intercept parameters for rain/snow/hail and
! densities for snow/hail based on values in history dump.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, Spring 2010
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE rsa_table
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Include files.
!-----------------------------------------------------------------------

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  IF( mphyopt < 2 .OR. mphyopt >= 200) THEN
!   TAS: Complain here
!   IF (myproc == 0) WRITE(6,'(/a,i4,/a,/a/)')                         &
!         ' WARNING: mphyopt (or rfopt) = ',mphyopt,                   &
!         ' is not valid for this work.',                              &
!         '          Reset mphyopt!!!'
!   CALL arpsstop(' Program stopped in set_dsd_para.',1)
!   STOP
  ENDIF

  CALL model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl,  &
       alpharain,alphasnow,alphagrpl,alphahail)

  IF (rhos <= 0.0) THEN
    rhos = 100.
  END IF

  IF (rhoh <= 0.0) THEN
    rhoh = 913.
  END IF

  IF (rhog <= 0.0) THEN

    SELECT CASE (mphyopt)
    CASE(1:12,108:110)
    rhog = 400.
    CASE(106,116)
    rhog = 500.

    END SELECT

  END IF

  IF (N0r <= 0.0) THEN
    N0r = 8.0E+06
  END IF

  IF (N0s <= 0.0) THEN
    N0s = 3.0E+06
  SELECT CASE (mphyopt)
  CASE(1:12,106,108:110,116)
    N0s2 = 0.0
  END SELECT
  END IF

  IF (N0h <= 0.0) THEN
    N0h = 4.0E+04
  END IF

  IF (N0g <= 0.0) THEN
    SELECT CASE (mphyopt)
    CASE(1:12,108:110)
    N0g = 4.0E+05
    CASE(106,116)
    N0g = 4.0E+06
    END SELECT
  END IF

   N0ms = N0s
   N0ms2 = N0s2
   N0mh = N0h
   N0mg = N0g

  IF (alphar <= 0.0) THEN
    SELECT CASE (mphyopt)
    CASE(1:12,106,108:110)
      alphar = 0.0
    CASE(116)
      alphar = 1.0
    END SELECT
   END IF

   IF (alphas <= 0.0) THEN
       alphas = 0.0
   END IF

   SELECT CASE (mphyopt)
   CASE(1:12,106,109,110,116)
     alphas2 = 0.0
   CASE(108)
     alphas2 = 0.6357 
   END SELECT

   IF (alphah <= 0.0) THEN
     alphah = 0.0
   END IF

   IF (alphag <= 0.0) THEN
      alphag = 0.0
   END IF

   lamdar = 0.0
   lamdas = 0.0
   lamdas2 = 0.0
   lamdams = 0.0
   lamdams2 = 0.0
   lamdag = 0.0
   lamdamg = 0.0
   lamdah = 0.0
   lamdamh = 0.0


  RETURN
END SUBROUTINE set_dsd_para


!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE rdr_obs                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE rdr_obs (rho,qscalar,kdph,obs_dual,var_dsd,    &
                       var_idx,dualpol)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! A shell subroutine to assign DSD parameters for the simulated
! radar parameters using parameterized formula based on Jung et al.(2008a).
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/14/2010
!
! MODIFICATION HISTORY:
!
!  Bryan Putnam 4/16/2013: Added in information for all radar parameters and
!  all operators, replaces rdr_obs_SM.
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE


 !added external fucntion calculate_kdp for additional kdp calculation
 !now passed through this subroutine
 REAL, EXTERNAL :: calculate_kdp
!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------
  REAL :: qscalar(nscalar)
  REAL :: rho

  INTEGER :: var_idx,dualpol

  TYPE(T_obs_dual) :: obs_dual
  TYPE(T_para_dsd) :: var_dsd
  REAL :: kdph



 !local variables
 REAL :: no_value = missing

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (mphyopt)
  CASE(2:8,106)  ! single moment schemes 
    SELECT CASE (qgh_opt)
      CASE (1)                       ! graupel off, hail off
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),no_value, &
              no_value, no_value, no_value, no_value, no_value, alphar,    &
              alphas,no_value,no_value)
      CASE (2)                       ! graupel off, hail on
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   qscalar(P_QH),no_value,no_value,no_value,no_value,  &
                   no_value,alphar,alphas,alphah,no_value)
      CASE (3)                       ! graupel on, hail off
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),no_value, &
                   qscalar(P_QG),no_value,no_value,no_value,no_value,    &
                   alphar,alphas,no_value,alphag)
      CASE (4)                       ! graupel on, hail on
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                  qscalar(P_QH),qscalar(P_QG),no_value,no_value,        &
                  no_value,no_value,alphar,alphas,alphah,alphag)
    END SELECT
  CASE(9:12,109:110) !double moment schemes
    SELECT CASE (qgh_opt)
      CASE (1)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   no_value,no_value,qscalar(P_NR),qscalar(P_NS),      &
                   no_value,no_value,alphar,alphas,          &
                   no_value,no_value)
      CASE (2)
         var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),       &
                   qscalar(P_QH),no_value,qscalar(P_NR),              &
                   qscalar(P_NS),qscalar(P_NH),no_value,alphar,       &
                   alphas,alphah,no_value)
      CASE (3)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   no_value,qscalar(P_QG),qscalar(P_NR),              &
                   qscalar(P_NS),no_value,qscalar(P_NG),alphar,       &
                   alphas,no_value,alphag)
      CASE (4)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   qscalar(P_QH),qscalar(P_QG),qscalar(P_NR),        &
                   qscalar(P_NS),qscalar(P_NH),qscalar(P_NG),        &
                   alphar,alphas,alphah,alphag)
     END SELECT
  CASE(108,116) ! double moment for rain only
    SELECT CASE (qgh_opt)
      CASE(1)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   no_value,no_value,qscalar(P_NR),no_value,no_value,    &
                   no_value,alphar,alphas,no_value,no_value)
      CASE(2)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                    qscalar(P_QH),no_value,qscalar(P_NR),no_value,     &
                    no_value,no_value,alphar,alphas,         &
                    alphah,no_value)
      CASE(3)
        var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),        &
                   no_value,qscalar(P_QG),qscalar(P_NR),no_value,      &
                   no_value,no_value,alphar,alphas,no_value,  &
                   alphag)
      CASE(4)
       var_dsd = assign_para_dsd_TM(qscalar(P_QR),qscalar(P_QS),         &
                   qscalar(P_QH),qscalar(P_QG),qscalar(P_NR),        &
                   no_value,no_value,no_value,alphar,alphas,  &
                   alphah,alphag)
    END SELECT
  END SELECT

  dualpol_opt = dualpol
  IF(dualpol == 1) THEN
     IF(var_idx <= 3) THEN
        obs_dual = calculate_obs(rho,var_dsd,var_idx)
     ELSE
        !kdph = calculate_kdp(rho,var_dsd,var_idx)
        kdph = calculate_kdp(rho,var_dsd)      ! Youngsun, Please check this, var_idx
                                               ! is not in the definition
     END IF
   ELSE !dualpol 2 or 3
      obs_dual = refl_rsa(rho,var_dsd)
   END IF

  RETURN

END SUBROUTINE rdr_obs


!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE para_dsd_SM                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE para_dsd_SM (qscalar,var_dsd)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! A shell subroutine to assign DSD parameters for the simulated
! radar parameters using T-matrix method.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/14/2010
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------
  REAL :: qscalar(nscalar)

  TYPE(T_para_dsd) :: var_dsd

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  SELECT CASE (qgh_opt)
    CASE (1)               ! graupel off, hail off
      var_dsd = assign_para_dsd_SM(qscalar(P_QR),qscalar(P_QS),     &
                missing,missing)
    CASE (2)               ! graupel off, hail on
      var_dsd = assign_para_dsd_SM(qscalar(P_QR),qscalar(P_QS),     &
                qscalar(P_QH),missing)
    CASE (3)               ! graupel on, hail off
      var_dsd = assign_para_dsd_SM(qscalar(P_QR),qscalar(P_QS),     &
                missing,qscalar(P_QG))
    CASE (4)               ! graupel on, hail on
      var_dsd = assign_para_dsd_SM(qscalar(P_QR),qscalar(P_QS),     &
                qscalar(P_QH),qscalar(P_QG))
  END SELECT

  RETURN

END SUBROUTINE para_dsd_SM

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE para_dsd_MM                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE para_dsd_MM (rhoa,qscalar,alpha,var_dsd)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! A shell subroutine to assign DSD parameters for the simulated
! radar parameters using T-matrix method.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/14/2010
!
! MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA
  USE radaremul_cst

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------
  REAL   :: qscalar(nscalar)
  REAL*8 :: alpha(6)
  REAL   :: rhoa

  TYPE(T_para_dsd) :: var_dsd
  REAL   :: qr,  qs,  qg,  qh
  REAL   :: Ntr, Nts, Ntg, Nth
  REAL   :: alfr,alfs,alfg,alfh

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  call calcMDR()

  SELECT CASE (mphyopt)
  CASE(2:4)
    qr = qscalar(P_QR)
    qs = qscalar(P_QS)
    qg = missing
    qh = qscalar(P_QH)

    alfr = alpha(P_QR)
    alfs = alpha(P_QS)
    alfg = 0.0
    alfh = alpha(P_QH)

    Ntr = 0.0; Nts = 0.0; Ntg = 0.0; Nth = 0.0

  CASE(5:7,106)
    qr = qscalar(P_QR)
    qs = qscalar(P_QS)
    qg = qscalar(P_Qg)
    qh = missing

    alfr = alpha(P_QR)
    alfs = alpha(P_QS)
    alfg = alpha(P_QG)
    alfh = 0.0

    Ntr = 0.0; Nts = 0.0; Ntg = 0.0; Nth = 0.0

  CASE(108)
    qr = qscalar(P_QR)
    qs = qscalar(P_QS)
    qg = qscalar(P_QG)
    qh = qscalar(P_QH)

    alfr = alpha(P_QR)
    alfs = alpha(P_QS)
    alfg = alpha(P_QG)
    alfh = alpha(P_QH)

    Ntr = 0.0; Nts = 0.0; Ntg = 0.0; Nth = 0.0

  CASE(9:11,109)
    qr = qscalar(P_QR)
    qs = qscalar(P_QS)
    qg = qscalar(P_QG)
    qh = qscalar(P_QH)

    alfr = alpha(P_QR)
    alfs = alpha(P_QS)
    alfg = alpha(P_QG)
    alfh = alpha(P_QH)

    Ntr = qscalar(P_NR)
    Nts = qscalar(P_NS)
    Ntg = qscalar(P_NG)
    Nth = qscalar(P_NH)

  CASE(110,116)
    qr = qscalar(P_QR)
    qs = qscalar(P_QS)
    qg = qscalar(P_QG)
    qh = missing

    alfr = alpha(P_QR)
    alfs = alpha(P_QS)
    alfg = alpha(P_QG)
    alfh = 0.0

    Ntr = qscalar(P_NR)
    Nts = qscalar(P_NS)
    Ntg = qscalar(P_NG)
    Nth = 0.0

  END SELECT


  SELECT CASE (qgh_opt)
    CASE (1)                    ! graupel off, hail off
      qg=missing; qh=missing; Ntg=0.0; Nth=0.0; alfg=0.0; alfh=0.0
    CASE (2)                    ! graupel off, hail on
      qg=missing; Ntg=0.0; alfg=0.0
    CASE (3)                    ! graupel on, hail off
      qh=missing; Nth=0.0; alfh=0.0
  END SELECT

  var_dsd = assign_para_dsd_TM(qr,qs,qh,qg,Ntr,Nts,Nth,Ntg,     &
                               alfr,alfs,alfh,alfg)

  RETURN

END SUBROUTINE para_dsd_MM

INTEGER FUNCTION get_qgh_opt(graupel_ON, hail_ON)

  INTEGER :: graupel_ON,hail_ON

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF(graupel_ON == 0 .and. hail_ON == 0) THEN
    get_qgh_opt = 1
  ELSE IF(graupel_ON == 0 .and. hail_ON == 1) THEN
    get_qgh_opt = 2
  ELSE IF(graupel_ON == 1 .and. hail_ON == 0) THEN
    get_qgh_opt = 3
  ELSE IF(graupel_ON == 1 .and. hail_ON == 1) THEN
    get_qgh_opt = 4
  ENDIF

END FUNCTION get_qgh_opt
