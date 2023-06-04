#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! subroutine DUSTRESB

! Purpose:
!   To calculate the surface layer resistance for mineral dust

! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards

! Documentation: "Modelling the atmospheric lifecycle..."
!                 Woodward, JGR106, D16, pp18155-18166
!---------------------------------------------------------------------

 SUBROUTINE dustresb(                                             &
  row_length,rows,                                                &
  pstar,tstar,rhostar,aresist,vshr,cd_std_dust,                   &
  r_b_dust                                                        &
  )

USE c_sulchm
USE dust_param, ONLY: ndiv, drep, rhop
USE c_r_cp
USE c_0_dg_c
USE c_pi
USE c_g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
        !IN
 row_length                                                       &
            !IN
,rows       !IN

REAL                                                              &
     !IN
 pstar(row_length,rows)                                           &
                              !IN surface pressure
,tstar(row_length,rows)                                           &
                              !IN surface temperature
,rhostar(row_length,rows)                                         &
                              !IN surface air density
,aresist(row_length,rows)                                         &
                              !IN aerodynamic resistance
,vshr(row_length,rows)                                            &
                              !IN surface to lowest lev windspeed
!                                   !   difference
,cd_std_dust(row_length,rows) !IN surface transfer coeffient for
!                                   !   momentum, excluding orographic
!                                   !   form drag

REAL                                                              &
     !OUT
 r_b_dust(row_length,rows,ndiv) !OUT surface layer resistance for
!                                     !    mineral dust

!     local variables

INTEGER                                                           &
 idiv                                                             &
      !loop counter, dust divisions
,i                                                                &
      !loop counter
,j                                                                &
      !loop counter
,lev1 !number of levels for vstokes calculation

REAL                                                              &
 nu(row_length,rows)                                              &
                           !kinematic viscosity
,etaa(row_length,rows)                                            &
                           !dynamic viscosity of air
,lamdaa(row_length,rows)                                          &
                           !mean free path of air molecules
,vstokes1(row_length,rows)                                        &
                           !gravitational settling velocity, lev1
,nstokes(row_length,rows)                                         &
                           !stokes number = VstokesVshrVshr/nu g
,nschmidt(row_length,rows)                                        &
                           !schmidt number = nu/diffusivit
,tc(row_length,rows)                                              &
                           !temperature in deg C
,alphaccf(row_length,rows)                                        &
                           !alpha in cunningham correction factor
,ccf(row_length,rows)                                             &
                           !Cunningham correction factor
,fvsq(row_length,rows)                                            &
                           !friction velocity squared
,work(row_length,rows)                                            &
                           !workspace
,stokes_exp                                                       &
                           !stokes term in R_B_DUST equation
,smallp                    !small +ve number, negligible compared to 1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

 EXTERNAL vgrav

IF (lhook) CALL dr_hook('DUSTRESB',zhook_in,zhook_handle)

!... epsilon() is defined as almost negligible, so eps/100 is negligible
smallp = EPSILON(1.0) / 100.0

!...calc stokes number, schmidt number and finally resistance

lev1=1

DO idiv=1,ndiv
! DEPENDS ON: vgrav
  CALL vgrav(                                                     &
  row_length,rows,lev1,drep(idiv),rhop,pstar,tstar,               &
  vstokes1,ccf,etaa                                               &
  )

!CDIR NOVECTOR
  DO j = 1,rows
    DO i= 1,row_length
      nschmidt(i,j)=3.*pi*etaa(i,j)*etaa(i,j)*drep(idiv)/         &
       (rhostar(i,j)*boltzmann*tstar(i,j)*ccf(i,j))
      nstokes(i,j)=vstokes1(i,j)*cd_std_dust(i,j)*rhostar(i,j)*   &
       vshr(i,j)*vshr(i,j)/(etaa(i,j)*g)
      ! Avoid underflow in Stokes term by setting to zero if
      ! negligible compared to Schmidt term, i.e., if NSTOKES
      ! is too small.
      IF ( 3.0 / nstokes(i,j) <                                   &
           - LOG10( smallp *nschmidt(i,j)**(-2./3.) ) ) THEN
         stokes_exp = 10.**(-3./nstokes(i,j))
      ELSE
         stokes_exp = 0.0
      END IF
      r_b_dust(i,j,idiv)=1./( SQRT(cd_std_dust(i,j)) *            &
       (nschmidt(i,j)**(-2./3.)+stokes_exp) )
    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NDIV

IF (lhook) CALL dr_hook('DUSTRESB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dustresb
#endif
