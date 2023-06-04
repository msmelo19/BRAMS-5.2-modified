#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  SUBROUTINE SF_RIB -----------------------------------------------
!
!  Purpose: Calculate bulk Richardson number for surface layer
!
! ------------------------------------------------------------------

!!   Subroutine interface
SUBROUTINE sf_rib (                                               &
 row_length,rows,points,tile_pts,                                 &
 pts_index,tile_index,                                            &
 bq_1,bt_1,qstar,q_elev,resft,t_elev,tstar,vshr,z0h,z0m,          &
 z1_tq,z1_uv,rib,db,ltimer                                        &
 )

USE c_g
USE c_r_cp

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
,tile_index(points) ! IN Index of tile points.

LOGICAL                                                           &
 ltimer              ! IN logical for TIMER

REAL                                                              &
 bq_1(row_length,rows)                                            &
                       ! IN A buoyancy parameter for lowest atm
!                            !    level. ("beta-q twiddle").
,bt_1(row_length,rows)                                            &
                       ! IN A buoyancy parameter for lowest atm
!                            !    level. ("beta-T twiddle").
,qstar(points)                                                    &
                       ! IN Surface saturated sp humidity.
,q_elev(points)                                                   &
                       ! IN Total water content of lowest
!                            !    atmospheric layer (kg per kg air).
,resft(points)                                                    &
                       ! IN Total resistance factor.
,t_elev(points)                                                   &
                       ! IN Liquid/frozen water temperature for
!                            !    lowest atmospheric layer (K).
,tstar(points)                                                    &
                       ! IN Surface temperature (K).
,vshr(row_length,rows)                                            &
                       ! IN Magnitude of surface-to-lowest-level
!                            !    wind shear.
,z0h(points)                                                      &
                       ! IN Roughness length for heat and
!                            !    moisture m
,z0m(points)                                                      &
                       ! IN Effective roughness length for
!                            !    momentum
,z1_tq(row_length,rows)                                           &
                       ! IN Height of lowest TQ level (m).
,z1_uv(row_length,rows)! IN Height of lowest UV level (m).

REAL                                                              &
 rib(points)                                                      &
                     ! OUT Bulk Richardson number for lowest layer
,db(points)          ! OUT Buoyancy difference between surface
!                          !     and lowest atmospheric level.


!  External routines called :-
EXTERNAL timer


!  Workspace --------------------------------------------------------
INTEGER                                                           &
 i,j                                                              &
                     ! Horizontal field index.
,k                                                                &
                     ! Tile field index.
,l                   ! Points field index.

REAL                                                              &
 dq(points)                                                       &
                     ! Sp humidity difference between surface
!                          ! and lowest atmospheric level (Q1 - Q*).
,dtemp(points)       ! Modified temperature difference between
!                            surface and lowest atmospheric level.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_RIB',zhook_in,zhook_handle)

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('SF_RIB  ',3)
END IF

!-----------------------------------------------------------------------
!!  1 Calculate temperature (strictly, liquid/ice static energy) and
!!    humidity jumps across the surface layer.
!-----------------------------------------------------------------------
DO k=1,tile_pts
  l = tile_index(k)
  j=(pts_index(l)-1)/row_length + 1
  i = pts_index(l) - (j-1)*row_length
  dtemp(l) = t_elev(l) - tstar(l) +                               &
                  (g/cp)*(z1_tq(i,j)+z0m(l)-z0h(l))
!                                                             ! P243.118
  dq(l) = q_elev(l) - qstar(l)                          ! P243.119
END DO

!-----------------------------------------------------------------------
!!  2 Calculate bulk Richardson numbers for the surface layer.
!-----------------------------------------------------------------------
DO k=1,tile_pts
  l = tile_index(k)
  j=(pts_index(l)-1)/row_length + 1
  i = pts_index(l) - (j-1)*row_length
  db(l) = g*(bt_1(i,j)*dtemp(l) + bq_1(i,j)*resft(l)*dq(l))
  rib(l) = z1_uv(i,j)*db(l) / ( vshr(i,j)*vshr(i,j) )
END DO

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('SF_RIB  ',4)
END IF

IF (lhook) CALL dr_hook('SF_RIB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_rib
#endif
