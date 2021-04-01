! Adopted from RADREMUL library of ARPS system
! Used for direct reflectity DA capability
!########################################################################
!########################################################################
!#########                                                      #########
!#########                     convert2radar                    #########
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

  USE DUALPOLEPARA

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

  USE DUALPOLEPARA

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

  USE DUALPOLEPARA

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
  USE DUALPOLEPARA

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

  USE DUALPOLEPARA
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

  USE DUALPOLEPARA
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

  USE DUALPOLEPARA

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

  USE DUALPOLEPARA
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

  USE DUALPOLEPARA
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

  USE DUALPOLEPARA
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
!        kdph = calculate_kdp(rho,var_dsd)      ! Youngsun, Please check this, var_idx
                                               ! is not in the definition
     END IF
!   ELSE !dualpol 2 or 3
!      obs_dual = refl_rsa(rho,var_dsd)
   END IF

  RETURN

END SUBROUTINE rdr_obs

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
