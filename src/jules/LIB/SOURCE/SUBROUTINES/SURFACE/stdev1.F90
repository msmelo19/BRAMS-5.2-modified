#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE STDEV1 ------------------------------------------------
!
!  Purpose: Calculate the standard deviations of layer 1 turbulent
!           fluctuations of temperature and humidity using approximate
!           formulae from first order closure.
!
!  -------------------------------------------------------------------

!  Subroutine interface
SUBROUTINE stdev1 (                                               &
 row_length,rows,points,tile_pts,pts_index,tile_index,fld_sea,    &
 bq_1,bt_1,fqw_1,ftl_1,rhokm_1,rhostar,vshr,z0m,z1_tq,tile_frac,  &
 q1_sd,t1_sd,ltimer                                               &
 )

USE c_g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER                                                           &
 row_length                                                       &
                       ! IN Number of X points?
,rows                                                             &
                       ! IN Number of Y points?
,points                                                           &
                       ! IN Total number of points.
,tile_pts                                                         &
                       ! IN Number of tile points.
,pts_index(points)                                                &
                       ! IN Index of points.
,tile_index(points)  ! IN Index of tile points.

LOGICAL                                                           &
 ltimer                ! IN logical for TIMER

REAL                                                              &
 fld_sea(row_length,rows)                                         &
                       ! IN Fraction of land or sea
,bq_1(row_length,rows)                                            &
                       ! IN Buoyancy parameter.
,bt_1(row_length,rows)                                            &
                       ! IN Buoyancy parameter.
,fqw_1(points)                                                    &
                       ! IN Surface flux of QW.
,ftl_1(points)                                                    &
                       ! IN Surface flux of TL.
,rhokm_1(points)                                                  &
                       ! IN Surface momentum exchange coefficient.
,rhostar(row_length,rows)                                         &
                         ! IN Surface air density.
,vshr(row_length,rows)                                            &
                       ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
,z0m(points)                                                      &
                       ! IN Roughness length for momentum.
,z1_tq(row_length,rows)                                           &
                       ! IN Height of lowest tq level.
,tile_frac(points)   ! IN Tile fraction.

REAL                                                              &
 q1_sd(row_length,rows)                                           &
                       ! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       specific humidity (kg/kg).
,t1_sd(row_length,rows)! INOUT Standard deviation of turbulent
!                            !       fluctuations of surface layer
!                            !       temperature (K).


!  External routines called :-
EXTERNAL timer

!  Workspace --------------------------------------------------------
INTEGER                                                           &
 i,j                                                              &
                       ! Horizontal field index.
,k                                                                &
                       ! Tile index.
,l                     ! Points field index.
REAL                                                              &
 vs                                                               &
                       ! Surface layer friction velocity
,vsf1_cubed                                                       &
                       ! Cube of surface layer free convective
!                            ! scaling velocity
,ws1                   ! Turbulent velocity scale for surface
!                            ! layer

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('STDEV1',zhook_in,zhook_handle)

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('STDEV1  ',3)
END IF

DO k=1,tile_pts
  l = tile_index(k)
  j=(pts_index(l)-1)/row_length + 1
  i = pts_index(l) - (j-1)*row_length

  vs = SQRT ( rhokm_1(l)/rhostar(i,j) * vshr(i,j) )
  vsf1_cubed = 1.25*g*(z1_tq(i,j) + z0m(l)) *                     &
             ( bt_1(i,j)*ftl_1(l) + bq_1(i,j)*fqw_1(l) ) /        &
                 rhostar(i,j)
  IF ( vsf1_cubed  >   0.0 ) THEN
    ws1 = ( vsf1_cubed + vs*vs*vs ) ** (1.0/3.0)
    t1_sd(i,j) = t1_sd(i,j) + MAX ( 0.0 ,                         &
                 fld_sea(i,j)*tile_frac(l)*1.93*ftl_1(l) /        &
                                        (rhostar(i,j)*ws1) )
    q1_sd(i,j) = q1_sd(i,j) + MAX ( 0.0 ,                         &
                 fld_sea(i,j)*tile_frac(l)*1.93*fqw_1(l) /        &
                                        (rhostar(i,j)*ws1) )
  END IF

END DO

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('STDEV1  ',4)
END IF

IF (lhook) CALL dr_hook('STDEV1',zhook_out,zhook_handle)
RETURN
END SUBROUTINE stdev1
#endif
