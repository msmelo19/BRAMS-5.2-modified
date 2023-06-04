#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE IM_SF_PT ----------------------------------------------

!  Purpose: Calculate implicit increments to surface variables

!  Documentation: UM Documentation Paper No 24.

!---------------------------------------------------------------------
!  Arguments :-
SUBROUTINE im_sf_pt (                                             &
 off_x,off_y,row_length,rows,n_rows,land_pts                      &
,land_index,ntiles,tile_index,tile_pts                            &
,flandg,tile_frac,snow_tile,ice_fract                             &
,gamma_in,alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile         &
,resft,dtstar_tile,dtstar                                         &
,rhokm_u_1,rhokm_v_1,rhokh_1,rhokh1_sice                          &
,ct_ctq_1,dqw_1,dtl_1,cq_cm_u_1,cq_cm_v_1,du_1,dv_1               &
,flandg_u,flandg_v                                                &
,fqw_gb,ftl_gb                                                    &
,taux_1,taux_land,taux_ssi,tauy_1,tauy_land,tauy_ssi              &
,fqw_tile,epot_tile,ftl_tile,fqw_ice,ftl_ice,e_sea,h_sea          &
,l_flux_bc,ltimer                                                 &
)

USE c_r_cp
USE c_lheat
USE surf_param, ONLY : ls
USE switches  , ONLY : l_epot_corr

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

LOGICAL ltimer                                                    &
,l_flux_bc                    ! Logical for prescribed
!                                   ! surface fluxes (SCM)

INTEGER                                                           &
 row_length                                                       &
                             ! IN Number of X points?
,rows                                                             &
                             ! IN Number of Y points?
,n_rows                                                           &
                             ! Local number of rows in a v field
,off_x                                                            &
                             ! Size of small halo in i
,off_y                                                            &
                             ! Size of small halo in j.
,land_pts                                                         &
                             ! IN Total number of land points.
,land_index(land_pts)                                             &
                             ! IN Index of land points.
,ntiles                                                           &
                             ! IN Number of land surface tiles.
,tile_index(land_pts,ntiles)                                      &
                             ! IN Index of tile points.
,tile_pts(ntiles)            ! IN Number of tiles.


REAL                                                              &
 flandg(row_length,rows)                                          &
                             ! IN Land fraction
,tile_frac(land_pts,ntiles)                                       &
                             ! IN Tile fraction
,snow_tile(land_pts,ntiles)                                       &
                             ! IN Lying snow on land tiles (kg/m2)
,ice_fract(row_length,rows)                                       &
                             ! IN Fraction of grid-box which is
!                                  !    sea-ice (decimal fraction).
,gamma_in                                                         &
                             ! IN Implicit weighting.
,alpha1(land_pts,ntiles)                                          &
                             ! IN Gradient of saturated specific
!                                  !    humidity with respect to
!                                  !    temperature between the bottom
!                                  !    model layer and the surface.
,alpha1_sice(row_length,rows)                                     &
                             ! IN ALPHA1 for sea-ice
,ashtf_prime(row_length,rows)                                     &
                             ! IN Adjusted SEB coefficient for
!                                  !    sea ice

,ashtf_prime_tile(land_pts,ntiles)                                &
                             ! IN Adjusted SEB coefficient for
!                                  !    land tiles
,dtstar_tile(land_pts,ntiles)                                     &
                             ! IN Change in TSTAR over timestep
!                                  !    for land tiles
,dtstar(row_length,rows)                                          &
                             ! IN Change in TSTAR over timestep
!                                  !    for sea-ice
,resft(land_pts,ntiles)                                           &
                             ! IN Total resistance factor
,epot_tile(land_pts,ntiles)                                       &
                             ! INOUT surface tile potential
,e_epot_tile(land_pts,ntiles)                                     &
                             ! Work ratio of explicit
!                                  !      EPOT/E
!                                  !    evaporation

,rhokm_u_1(row_length,rows)                                       &
                             ! IN Level 1 exchange coefficient for
!                                  !    momentum
,rhokm_v_1(row_length,n_rows)                                     &
                             ! IN Level 1 exchange coefficient for
!                                  !    momentum
,rhokh_1(land_pts,ntiles)                                         &
                             ! IN Surface exchange coeffs for FTL

,rhokh1_sice(row_length,rows)                                     &
                             ! IN Sea and sea-ice surface exchange
,ct_ctq_1(row_length,rows)                                        &
                             ! IN Coefficient in T and q
!                                  !     tri-diagonal implicit matrix
,cq_cm_u_1(row_length,rows)                                       &
                             ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
,cq_cm_v_1(row_length,n_rows)                                     &
                             ! IN Coefficient in U and V
!                                  !     tri-diagonal implicit matrix
,dqw_1(row_length,rows)                                           &
                             ! IN Level 1 increment to q field
,dtl_1(row_length,rows)                                           &
                             ! IN Level 1 increment to T field
,du_1(1-off_x:row_length+off_x,                                   &
      1-off_y:rows+off_y)                                         &
                             ! IN Level 1 increment to u wind
!                                  !    field
,dv_1(1-off_x:row_length+off_x,                                   &
      1-off_y:n_rows+off_y)                                       &
                             ! IN Level 1 increment to v wind
!                                  !    field
,flandg_u(row_length,rows)                                        &
                             ! IN Land fraction on U grid.
,flandg_v(row_length,n_rows) ! IN Land fraction on V grid.


REAL                                                              &
 fqw_gb(row_length,rows)                                          &
                             ! INOUT Grid-box value of QW flux at
!                                  !       Kg/sq m/s
,ftl_gb(row_length,rows)                                          &
                             ! INOUT Grid-box value of TL flux at
!                                  !       i.e. H/Cp where H is sensible
!                                  !       in W per sq m).
,taux_1(row_length,rows)                                          &
                             ! OUT   x-component of turbulent
!                                  !       stress at surface.
,taux_land(row_length,rows)                                       &
                             ! INOUT x-component of turbulent
!                                  !       stress at land surface.
,taux_ssi(row_length,rows)                                        &
                             ! INOUT x-component of turbulent
!                                  !       stress at sea surface.
,tauy_1(row_length,n_rows)                                        &
                             ! OUT   y-component of turbulent
!                                  !       stress at surface.
,tauy_land(row_length,n_rows)                                     &
                             ! INOUT y-component of turbulent
!                                  !       stress at land surface.
,tauy_ssi(row_length,n_rows)                                      &
                             ! INOUT y-component of turbulent
!                                  !       stress at sea surface.
,fqw_tile(land_pts,ntiles)                                        &
                             ! INOUT Tile flux of QW. Kg/sq m/s
,ftl_tile(land_pts,ntiles)                                        &
                             ! INOUT Tile flux of TL
,e_sea(row_length,rows)                                           &
                             ! INOUT Evaporation from sea times
!                                  !       leads fraction (kg/m2/s).
!                                  !       Zero over land.
,h_sea(row_length,rows)      ! INOUT Surface sensible heat flux ov
!                                  !       sea times leads fraction (W/m
!                                  !       Zero over land.



!  External references :-
EXTERNAL timer


!  Local and other symbolic constants :-

! Workspace :-
REAL                                                              &
 fqw_ice(row_length,rows)                                         &
                             ! "Explicit" surface flux of QW for
!                                  !  sea-ice fraction of gridsquare.
,ftl_ice(row_length,rows)                                         &
                             ! "Explicit" surface flux of TL for
!                                  !  sea-ice fraction of gridsquare.
,rhokpm(land_pts,ntiles)                                          &
                             !  Surface exchange coeff for tiles
,rhokpm_sice(row_length,rows)                                     &
                             !  Sea-ice surface exchange coeff
,lat_ht                                                           &
          ! Latent heat of evaporation for snow-free land
!               ! or sublimation for snow-covered land and ice.
,apart(row_length,rows,2)                                         &
                             ! Tempary array
,bpart(row_length,rows,2)                                         &
                             ! Tempary array
,recip(row_length,rows)                                           &
                             ! Tempary array
,ftl_land(row_length,rows)                                        &
                                ! Tempary array
,fqw_land(row_length,rows)      ! Tempary array

LOGICAL                                                           &
 epot_calc(land_pts,ntiles)  ! flag to connect first and second
                             ! epot calculations

!  Local scalars :-
INTEGER                                                           &
 i,j                                                              &
          ! Loop counter (horizontal field index).
,k                                                                &
          ! Loop counter (tile index).
,l                                                                &
          ! Loop counter (horizontal land index).
,n        ! Loop counter (tile counter).

REAL                                                              &
 ftl_old                                                          &
          ! Used to hold current value of FTL_GB before updating
,gamma_1     ! local implicit weight

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!!  0.  Check that the scalars input to define the grid are consistent.
!       See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('IM_SF_PT',zhook_in,zhook_handle)

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('IM_SF_PT ',3)
END IF

gamma_1 = gamma_in
IF (l_flux_bc) THEN
  gamma_1 = 0.0       ! use GAMMA=0 for scalars (specified fluxes)
END IF
! Initialise APART and BPART to zero
DO j=1,rows
  DO i=1,row_length
    apart(i,j,1)=0.0
    apart(i,j,2)=0.0
    bpart(i,j,1)=0.0
    bpart(i,j,2)=0.0
    ftl_land(i,j)=0.0
    fqw_land(i,j)=0.0
  END DO
END DO

!-------------------------------------------------------------------------
! initialise epot_calc
!-------------------------------------------------------------------------
epot_calc(:,:) = .FALSE.


! Land tiles
DO n=1,ntiles
!CDIR NODEP
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    lat_ht = lc
    IF (snow_tile(l,n) >  0.) lat_ht = ls

      rhokpm(l,n) = rhokh_1(l,n) / ( ashtf_prime_tile(l,n) +        &
                 rhokh_1(l,n)*(lat_ht*alpha1(l,n)*resft(l,n) + cp) )

      apart(i,j,1)=apart(i,j,1) - tile_frac(l,n) *                  &
                 gamma_1 * rhokpm(l,n) *                            &
              ( lat_ht*resft(l,n)*rhokh_1(l,n)*alpha1(l,n) +        &
                           ashtf_prime_tile(l,n) )
      apart(i,j,2)=apart(i,j,2) + tile_frac(l,n) *                  &
                 gamma_1 * rhokpm(l,n) *                            &
                 lat_ht*resft(l,n)*rhokh_1(l,n)
      bpart(i,j,1)=bpart(i,j,1) + tile_frac(l,n) *                  &
                 gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
                 cp*rhokh_1(l,n)*alpha1(l,n)
      bpart(i,j,2)=bpart(i,j,2) - tile_frac(l,n) *                  &
                 gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
                 ( cp*rhokh_1(l,n) + ashtf_prime_tile(l,n) )

  END DO
END DO



! Sea points
DO j=1,rows
  DO i=1,row_length

    IF(flandg(i,j) <  1.0 .AND. ice_fract(i,j) >  0.0) THEN
! Sea ice point
      rhokpm_sice(i,j) = rhokh1_sice(i,j) / ( ashtf_prime(i,j) +    &
                      rhokh1_sice(i,j)*(ls*alpha1_sice(i,j) + cp) )

      apart(i,j,1)=flandg(i,j)*apart(i,j,1)                         &
         - gamma_1 * (1.0-flandg(i,j)) * ice_fract(i,j)             &
        * rhokpm_sice(i,j) *                                        &
        ( ls*rhokh1_sice(i,j)*alpha1_sice(i,j) + ashtf_prime(i,j) ) &
         - gamma_1 * (1.0-flandg(i,j)) * ( 1.0 - ice_fract(i,j) )   &
        * rhokh1_sice(i,j)

      apart(i,j,2)=flandg(i,j)*apart(i,j,2)                         &
         + gamma_1 * (1.0-flandg(i,j)) *ice_fract(i,j)              &
         * rhokpm_sice(i,j) * ls*rhokh1_sice(i,j)

      bpart(i,j,1)=flandg(i,j)*bpart(i,j,1)                         &
         + gamma_1 * ice_fract(i,j) * ( 1.0 - flandg(i,j) )         &
         * rhokpm_sice(i,j) *cp*rhokh1_sice(i,j)*alpha1_sice(i,j)

      bpart(i,j,2)=flandg(i,j)*bpart(i,j,2)                         &
         - gamma_1 * ice_fract(i,j) * ( 1.0 - flandg(i,j) )         &
         * rhokpm_sice(i,j)                                         &
         * ( cp*rhokh1_sice(i,j) + ashtf_prime(i,j) )               &
         - gamma_1 * ( 1.0 - ice_fract(i,j) )                       &
         * ( 1.0 - flandg(i,j) ) * rhokh1_sice(i,j)

    ELSE IF(flandg(i,j) <  1.0 .AND..NOT.ice_fract(i,j) >  0.0) THEN
! Ordinary sea point
      apart(i,j,1)= flandg(i,j)*apart(i,j,1)                        &
         - gamma_1 * ( 1.0 - flandg(i,j) ) * rhokh1_sice(i,j)
      apart(i,j,2)= flandg(i,j)*apart(i,j,2)

      bpart(i,j,1)= flandg(i,j)*bpart(i,j,1)
      bpart(i,j,2)= flandg(i,j)*bpart(i,j,2)                        &
         - gamma_1 * ( 1.0 - flandg(i,j) ) * rhokh1_sice(i,j)

    END IF
  END DO
END DO

! Land tiles
DO n=1,ntiles
!CDIR NODEP
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    e_epot_tile(l,n)=1.0
    IF(      (epot_tile(l,n) >  0.0)                              &
        .AND.(fqw_tile(l,n)  > SQRT(TINY(0.0))) )THEN
      e_epot_tile(l,n)=epot_tile(l,n)/fqw_tile(l,n)
      epot_calc(l,n)=.TRUE.
    ENDIF
  END DO
END DO



! Calculate grid-box fluxes of heat and moisture
DO j=1,rows
  DO i=1,row_length
    recip(i,j)=( 1.0 + ct_ctq_1(i,j)*apart(i,j,1) ) *               &
             ( 1.0 + ct_ctq_1(i,j)*bpart(i,j,2) ) -                 &
             ct_ctq_1(i,j)*apart(i,j,2)*ct_ctq_1(i,j)*bpart(i,j,1)

    ftl_old=ftl_gb(i,j)

    ftl_gb(i,j) = ( ( 1.0 + ct_ctq_1(i,j)*bpart(i,j,2) ) *          &
                 ( ftl_old + apart(i,j,1)*dtl_1(i,j) +              &
                   apart(i,j,2)*dqw_1(i,j)) -                       &
                   ct_ctq_1(i,j)*apart(i,j,2) * ( fqw_gb(i,j) +     &
                   bpart(i,j,1)*dtl_1(i,j) +                        &
                   bpart(i,j,2)*dqw_1(i,j)) ) / recip(i,j)

    fqw_gb(i,j) = ( ( 1.0 + ct_ctq_1(i,j)*apart(i,j,1) ) *          &
                  ( fqw_gb(i,j) + bpart(i,j,1)*dtl_1(i,j) +         &
                    bpart(i,j,2)*dqw_1(i,j)) -                      &
                    ct_ctq_1(i,j)*bpart(i,j,1) * ( ftl_old +        &
                    apart(i,j,1)*dtl_1(i,j) +                       &
                    apart(i,j,2)*dqw_1(i,j)) ) / recip(i,j)

  END DO
END DO


! Make implicit correction to tile fluxes

! Land tiles
DO n=1,ntiles
!CDIR NODEP
  DO k=1,tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    lat_ht = lc
    IF (snow_tile(l,n) >  0.) lat_ht = ls

      ftl_tile(l,n) = ftl_tile(l,n) -                               &
                 gamma_1 * rhokpm(l,n) *                            &
              ( lat_ht*resft(l,n)*rhokh_1(l,n)*alpha1(l,n) +        &
                         ashtf_prime_tile(l,n) ) *                  &
           ( dtl_1(i,j) - ct_ctq_1(i,j)*ftl_gb(i,j) ) +             &
                 gamma_1 * rhokpm(l,n) *                            &
                 lat_ht*resft(l,n)*rhokh_1(l,n) *                   &
           ( dqw_1(i,j) - ct_ctq_1(i,j)*fqw_gb(i,j) )

      fqw_tile(l,n) = fqw_tile(l,n) +                               &
                 gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
                 cp*rhokh_1(l,n)*alpha1(l,n) *                      &
           ( dtl_1(i,j) - ct_ctq_1(i,j)*ftl_gb(i,j) ) -             &
                 gamma_1 * resft(l,n)*rhokpm(l,n) *                 &
                 ( cp*rhokh_1(l,n) + ashtf_prime_tile(l,n) ) *      &
           ( dqw_1(i,j) - ct_ctq_1(i,j)*fqw_gb(i,j) )


      IF (l_epot_corr) THEN
        IF (epot_calc(l,n)) THEN
          epot_tile(l,n)=fqw_tile(l,n)*e_epot_tile(l,n)
        END IF
      ELSE
        epot_tile(l,n)=fqw_tile(l,n)*e_epot_tile(l,n)
      END IF

      fqw_land(i,j)=fqw_land(i,j)+fqw_tile(l,n)*tile_frac(l,n)
      ftl_land(i,j)=ftl_land(i,j)+ftl_tile(l,n)*tile_frac(l,n)

      dtstar_tile(l,n) = dtstar_tile(l,n) + gamma_1 *               &
               ( cp * rhokh_1(l,n) *                                &
                 ( dtl_1(i,j) - ct_ctq_1(i,j) * ftl_gb(i,j) ) +     &
                 lat_ht * resft(l,n) * rhokh_1(l,n) *               &
                 ( dqw_1(i,j) - ct_ctq_1(i,j) * fqw_gb(i,j) ) ) /   &
               ( cp * rhokh_1(l,n) *                                &
                 ( cp + lat_ht * resft(l,n) * alpha1(l,n) ) +       &
                 ashtf_prime_tile(l,n) )

  END DO
END DO



! Sea points
DO j=1,rows
  DO i=1,row_length

    IF(flandg(i,j) <  1.0 .AND. ice_fract(i,j) >  0.0) THEN
! Sea ice point
      h_sea(i,j) = h_sea(i,j) - gamma_1 * rhokh1_sice(i,j) *        &
                       ( dtl_1(i,j) - ct_ctq_1(i,j) * ftl_gb(i,j) )

      e_sea(i,j) = e_sea(i,j) - gamma_1 * rhokh1_sice(i,j) *        &
                       ( dqw_1(i,j) - ct_ctq_1(i,j) * fqw_gb(i,j) )

      ftl_ice(i,j) = ( ( ftl_gb(i,j)                                &
         - ftl_land(i,j) * flandg(i,j)) / ( 1. - flandg(i,j) )      &
         - h_sea(i,j) * (1. - ice_fract(i,j)) ) /                   &
           ice_fract(i,j)

      fqw_ice(i,j) = ( ( fqw_gb(i,j)                                &
          - fqw_land(i,j) * flandg(i,j) ) / (1. - flandg(i,j) )     &
          - e_sea(i,j) * ( 1. - ice_fract(i,j) ) ) /                &
            ice_fract(i,j)

    ELSE IF(flandg(i,j) <  1.0 .AND. .NOT.ice_fract(i,j) >  0.0) THEN
! Ordinary sea point
      h_sea(i,j) = (ftl_gb(i,j)                                     &
         - ftl_land(i,j) * flandg(i,j)) / ( 1. - flandg(i,j) )

      e_sea(i,j) = ( fqw_gb(i,j)                                    &
         - fqw_land(i,j) * flandg(i,j) ) / ( 1. - flandg(i,j) )

      ftl_ice(i,j)=0.0
      fqw_ice(i,j)=0.0

    END IF
  END DO
END DO


gamma_1=gamma_in       ! use input GAMMA for winds

DO j=1,rows
  DO i=1,row_length

  IF(flandg_u(i,j) >  0.0)THEN
    taux_land(i,j) = ( taux_land(i,j) +                           &
                 gamma_1*rhokm_u_1(i,j)*du_1(i,j) ) /             &
                ( 1.0 + gamma_1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
  ELSE
    taux_land(i,j) = 0.0
  END IF

  IF(flandg_u(i,j) <  1.0)THEN
    taux_ssi(i,j) = ( taux_ssi(i,j) +                             &
                 gamma_1*rhokm_u_1(i,j)*du_1(i,j) ) /             &
                ( 1.0 + gamma_1*rhokm_u_1(i,j)*cq_cm_u_1(i,j) )
  ELSE
    taux_ssi(i,j) = 0.0
  END IF

  taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                      &
                + ( 1.0-flandg_u(i,j))*taux_ssi(i,j)

  END DO
END DO


DO j=1,n_rows
  DO i=1,row_length

  IF(flandg_v(i,j) >  0.0)THEN
    tauy_land(i,j) = ( tauy_land(i,j) +                           &
                 gamma_1*rhokm_v_1(i,j)*dv_1(i,j) ) /             &
                ( 1.0 + gamma_1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
  ELSE
    tauy_land(i,j) = 0.0
  END IF

  IF(flandg_v(i,j) <  1.0)THEN
    tauy_ssi(i,j) = ( tauy_ssi(i,j) +                             &
                 gamma_1*rhokm_v_1(i,j)*dv_1(i,j) ) /             &
                ( 1.0 + gamma_1*rhokm_v_1(i,j)*cq_cm_v_1(i,j) )
  ELSE
    tauy_ssi(i,j) = 0.0
  END IF

  tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                      &
                + ( 1.0-flandg_v(i,j))*tauy_ssi(i,j)

  END DO
END DO

IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('IM_SF_PT ',4)
END IF

IF (lhook) CALL dr_hook('IM_SF_PT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE im_sf_pt
#endif
