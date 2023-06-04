#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE SF_EXPL------------------------------------------------
!
!  Purpose: Calculate explicit surface fluxes of heat, moisture and
!           momentum. Also calculates surface exchange coefficients
!           required for implicit update of surface fluxes and surface
!           information required by the explicit boundary layer
!           routine
!
!
!  Documentation: UMDP 24.
!
!---------------------------------------------------------------------

!    Arguments :-
SUBROUTINE sf_expl (                                              &
! IN values defining field dimensions and subset to be processed :
 halo_i,halo_j,off_x,off_y,row_length,rows,n_rows,                &
 land_pts,land_pts_trif,npft_trif,                                &
 dim_cs1, dim_cs2,                                                &
! IN  parameters for iterative SISL scheme
 numcycles, cycleno,                                              &
! IN parameters required from boundary-layer scheme :
 bq_1,bt_1,z1_uv,z1_uv_top,z1_tq,z1_tq_top,qw_1,tl_1,             &
! IN soil/vegetation/land surface data :
 land_index,land_mask,formdrag,fd_stab_dep,orog_drag_param,       &
 ntiles,sm_levels,                                                &
 canopy,catch,catch_snow,hcon,ho2r2_orog,                         &
 fland,flandg,                                                    &
 snow_tile,sil_orog_land,smvccl,smvcst,smvcwt,                    &
 sthf,sthu,z0_tile,                                               &
! IN sea/sea-ice data :
 ice_fract,u_0,v_0,u_0_p,v_0_p,charnock,seasalinityfactor,        &
! IN everything not covered so far :
 pstar,lw_down,rad_sice,sw_tile,timestep,zh,ddmfx,                &
 co2_mmr,co2_3d,co2_dim_len,co2_dim_row,l_co2_interactive,        &
 l_phenol,l_triffid,l_q10,asteps_since_triffid,                   &
 cs,frac,canht_ft,photosynth_act_rad,lai_ft,lq_mix_bl,            &
 t_soil,ti,tstar,                                                 &
 tstar_land,tstar_sea,tstar_sice,tstar_ssi,                       &
 tstar_tile,z_land,l_ctile,cor_ust,                               &
 albsoil,cos_zenith_angle, ilayers,                               &
 u_1,v_1,u_1_p,v_1_p,                                             &
 l_dust,l_dust_diag,anthrop_heat,soil_clay,o3,                    &
! IN STASH flags :-
 sfme,sq1p5,st1p5,su10,sv10,sz0heff,                              &
! INOUT data :
 z0msea,l_spec_z0,z0m_scm,z0h_scm,gs,                             &
 g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                  &
! OUT Diagnostic not requiring STASH flags :
 cd,ch,recip_l_mo_sea,e_sea,fqw_1,                                &
 ftl_1,ftl_tile,le_tile,h_sea,radnet_sice,radnet_tile,            &
 rhokm_1,rhokm_u_1,rhokm_v_1,rib,rib_tile,taux_1,tauy_1,          &
 taux_land,taux_ssi,tauy_land,tauy_ssi,                           &
! OUT diagnostic requiring STASH flags :
 fme,                                                             &
! OUT diagnostics required for soil moisture nudging scheme :
 wt_ext,ra,                                                       &
! OUT data required for tracer mixing :
 rho_aresist,aresist,resist_b,                                    &
 rho_aresist_tile,aresist_tile,resist_b_tile,                     &
!OUT data required for mineral dust scheme
 r_b_dust,cd_std_dust,u_s_std_tile,                               &
! OUT data required for 4D-VAR :
 rho_cd_modv1,                                                    &
! OUT data required elsewhere in UM system :
 fb_surf,u_s,t1_sd,q1_sd,                                         &
! OUT data required elsewhere in boundary layer or surface code
 alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile,fqw_tile,        &
 epot_tile,fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,             &
 rhokh,rhokh_tile,rhokh_sice,rhokpm,rhokpm_pot,rhokpm_sice,       &
 rhokh_mix,dtstar_tile,dtstar,                                    &
 h_blend_orog,z0hssi,z0h_tile,z0h_eff,z0m_gb,z0mssi,z0m_tile,     &
 z0m_eff,cdr10m_u,cdr10m_v,chr1p5m,chr1p5m_sice,smc,hcons,        &
 vshr,vshr_land,vshr_ssi,                                         &
 gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,                             &
 resp_p_ft,resp_s,resp_s_tot,resp_w_ft,                           &
 gc,canhc_tile,wt_ext_tile,flake,                                 &
 tile_index,tile_pts,tile_frac,fsmc,                              &
 flandg_u,flandg_v,                                               &
 emis_tile,emis_soil,                                             &
! OUT data for ozone
 flux_o3_ft,fo3_ft,                                               &
! LOGICAL LTIMER
 ltimer,                                                          &
 l_ukca                                                           &
 )

USE dust_param, ONLY: ndiv
USE c_0_dg_c
USE c_r_cp
USE c_g
USE csigma
USE c_mdi

! #include "soil_thick.h"
USE soil_param, ONLY : dzsoil

USE snow_param, ONLY : ds

USE nstypes, ONLY :                                               &
!      imported scalars with intent(in)
   npft,ntype

USE fldtype
USE veg_param, ONLY : secs_per_360days

USE switches, ONLY: l_aggregate, can_model, can_rad_mod, &
                    buddy_sea

USE ancil_info, ONLY: nsmax,ssi_pts,sea_pts,sice_pts,             &
    ssi_index,sea_index,sice_index,fssi,sea_frac,sice_frac
USE prognostics, ONLY: nsnow,sice,sliq,snowdepth,tsnow
USE c_elevate, ONLY: surf_hgt

USE blopt8a, ONLY : on
USE solinc_data, ONLY: sky, l_skyview

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.
!    Checked for consistency.

INTEGER                                                           &
 row_length                                                       &
            ! Local number of points on a row
,rows                                                             &
            ! Local number of rows in a theta field
,n_rows                                                           &
            ! Local number of rows in a v field
,off_x                                                            &
            ! Size of small halo in i.
,off_y                                                            &
            ! Size of small halo in j.
,halo_i                                                           &
            ! Size of halo in i direction.
,halo_j                                                           &
            ! Size of halo in j direction.
,land_pts                                                         &
            ! IN No of land points being processed.
,land_pts_trif                                                    &
!                 ! IN For dimensioning land fields
,npft_trif  ! IN For dimensioning PFT fields available only
!                 !    with TRIFFID. Set to NPFT when TRIFFID on,
!                 !    set to 1 when TRIFFID off.

INTEGER                                                           &
  ilayers                                                         &
!                                  !No of layers in canopy radiation model
,numcycles                                                        &
            ! Number of cycles (iterations) for iterative SISL.
,cycleno                                                          &
            ! Iteration no
,dim_cs1, dim_cs2      ! soil carbon dimensions


! Defining vertical grid of model atmosphere.
REAL                                                              &
 bq_1(row_length,rows)                                            &
                             ! IN A buoyancy parameter
!                                  !    (beta q tilde).
,bt_1(row_length,rows)                                            &
                             ! IN A buoyancy parameter
!                                  !    (beta T tilde).
,z1_uv(row_length,rows)                                           &
                             ! IN Height of lowest uv level (m).
,z1_tq(row_length,rows)                                           &
                             ! IN Height of lowest tq level (m).
!                                  !    Note, if the grid used is
!                                  !    staggered in the vertical,
!                                  !    Z1_UV and Z1_TQ can be
!                                  !    different.
,qw_1(row_length,rows)                                            &
                             ! IN Total water content
,tl_1(row_length,rows)       ! IN Ice/liquid water temperature

REAL, INTENT(IN) :: z1_uv_top(row_length, rows)
                             ! Height of top of lowest uv-layer
REAL, INTENT(IN) :: z1_tq_top(row_length, rows)
                             ! Height of top of lowest Tq-layer


! (c) Soil/vegetation/land surface parameters (mostly constant).

LOGICAL                                                           &
 land_mask(row_length,rows)                                       &
                             ! IN T if land, F elsewhere.
,l_co2_interactive                                                &
                             ! IN Switch for 3D CO2 field
,l_phenol                                                         &
                             ! IN Indicates whether phenology
!                                  !    in use
,l_triffid                                                        &
                             ! IN Indicates whether TRIFFID
!                                  !    in use.
,l_ctile                                                          &
                             ! IN True if coastal tiling
,l_spec_z0                                                        &
                             ! IN T if using prescribed
!                                  !    sea surface roughness lengths
,l_q10                       ! IN True if using Q10 for soil resp

INTEGER                                                           &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                  !    point in ROW_LENGTH,ROWS is the
!                                  !    land point.

INTEGER                                                           &
 sm_levels                                                        &
                             ! IN No. of soil moisture levels
,ntiles                                                           &
                             ! IN No. of land-surface tiles
,co2_dim_len                                                      &
                             ! IN Length of a CO2 field row.
,co2_dim_row                                                      &
                             ! IN Number of CO2 field rows.
,asteps_since_triffid                                             &
                             ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
,formdrag                                                         &
                             ! IN Switch for orographic drag
,fd_stab_dep                 ! IN Switch to implement stability
!                                  !    dependence of orog form drag

REAL                                                              &
 canopy(land_pts,ntiles)                                          &
                             ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
,catch(land_pts,ntiles)                                           &
                             ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
,catch_snow(land_pts,ntiles)                                      &
                             ! IN Snow interception capacity of
!                                  !    tiles (kg/m2).
,hcon(land_pts)                                                   &
                             ! IN Soil thermal conductivity
!                                  !    (W/m/K).
,snow_tile(land_pts,ntiles)                                       &
                             ! IN Lying snow on tiles (kg/m2)
,smvccl(land_pts,sm_levels)                                       &
                             ! IN Critical volumetric SMC
!                                  !    (cubic m per cubic m of soil).
,smvcst(land_pts,sm_levels)                                       &
                             ! IN Volumetric saturation point
!                                  !    (m3/m3 of soil).
,smvcwt(land_pts,sm_levels)                                       &
                             ! IN Volumetric wilting point
!                                  !    (cubic m per cubic m of soil).
,sthf(land_pts,sm_levels)                                         &
                             ! IN Frozen soil moisture content of
!                                  !    each layer as a fraction of
!                                  !    saturation.
,sthu(land_pts,sm_levels)                                         &
                             ! IN Unfrozen soil moisture content
!                                  !    of each layer as a fraction of
!                                  !    saturation.
,z0_tile(land_pts,ntiles)                                         &
                             ! IN Tile roughness lengths (m).
,sil_orog_land(land_pts)                                          &
                             ! IN Silhouette area of unresolved
!                                  !    orography per unit horizontal
!                                  !    area on land points only.
,ho2r2_orog(land_pts)                                             &
                             ! IN Standard Deviation of orography.
!                                  !    equivilent to peak to trough
!                                  !    height of unresolved orography
,orog_drag_param                                                  &
!                                  ! IN Orographic form drag coefficient
,fland(land_pts)                                                  &
                             ! IN Land fraction on land tiles.
,flandg(row_length,rows)
!                                  ! IN Land fraction on all tiles.
!                                  !    divided by 2SQRT(2) on land
!                                  !    points only (m)

! (d) Sea/sea-ice data.

REAL                                                              &
 ice_fract(row_length,rows)                                       &
                             ! IN Fraction of gridbox covered by
!                                  !    sea-ice (decimal fraction).
,u_0(row_length,rows)                                             &
                             ! IN W'ly component of surface
!                                  !    current (m/s).
,v_0(row_length,n_rows)                                           &
                             ! IN S'ly component of surface
!                                  !    current (m/s).
,u_0_p(row_length,rows)                                           &
                             ! IN W'ly component of surface
!                                       current (m/s). P grid
,v_0_p(row_length,rows)                                           &
                             ! IN S'ly component of surface
!                                       current (m/s). P grid
 ,charnock   ! Charnock parameter for sea surface


! (f) Atmospheric + any other data not covered so far, incl control.

REAL                                                              &
 pstar(row_length,rows)                                           &
                             ! IN Surface pressure (Pascals).
,lw_down(row_length,rows)                                         &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
,rad_sice(row_length,rows)                                        &
                             ! IN Surface net shortwave and
!                                  !    downward LWradiation for
!                                  !    sea-ice (W/sq m).
,sw_tile(land_pts,ntiles)                                         &
                             ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
,timestep                                                         &
                             ! IN Timestep (seconds).
,zh(row_length,rows)                                              &
                             ! IN Height above surface of top of
!                            !    boundary layer (metres).
,ddmfx(row_length,rows)                                           &  
!                            ! IN Convective downdraught  
!                            !    mass-flux at cloud base
,co2_mmr                                                          &
                             ! IN CO2 Mass Mixing Ratio
,co2_3d(co2_dim_len,co2_dim_row)                                  &
!                                  ! IN 3D CO2 field if required.
,cs(land_pts,dim_cs1)                                             &
                        ! IN Soil carbon (kg C/m2).
,frac(land_pts,ntype)                                             &
                             ! IN Fractions of surface types.
,canht_ft(land_pts,npft)                                          &
                             ! IN Canopy height (m)
,photosynth_act_rad(row_length,rows)                              &
!                                  ! IN Net downward shortwave radiation
!                                  !    in band 1 (w/m2).
,lai_ft(land_pts,npft)                                            &
                             ! IN Leaf area index
,t_soil(land_pts,sm_levels)                                       &
                             ! IN Soil temperatures (K).
,ti(row_length,rows)                                              &
                             ! IN Sea-ice surface layer
,tstar_land(row_length,rows)                                      &
                             ! IN Land mean surface temperature (K
,tstar_sea(row_length,rows)                                       &
                             ! IN Open sea surface temperature (K)
,tstar_sice(row_length,rows)                                      &
                             ! IN Sea-ice surface temperature (K).
,tstar_ssi(row_length,rows)                                       &
                             ! IN mean sea surface temperature (K)
!                                  !    temperature (K).
,tstar(row_length,rows)                                           &
                             ! IN GBM surface temperature (K).
,tstar_tile(land_pts,ntiles)                                      &
                             ! IN Surface tile temperatures
,z_land(row_length,rows)                                          &
                             ! IN Land height (m).
,albsoil(land_pts)                                                &
!                                  ! Soil albedo.
, cos_zenith_angle(row_length, rows)                              &
!                                  ! Cosine of the zenith angle
,u_1(1-off_x:row_length+off_x,1-off_y:rows+off_y)                 &
!                                  ! IN W'ly wind component (m/s)
,v_1(1-off_x:row_length+off_x,1-off_y:n_rows+off_y)               &
!                                  ! IN S'ly wind component (m/s)
,u_1_p(row_length,rows)                                           &
                             ! IN U_1 on P-grid.
,v_1_p(row_length,rows)                                           &
                             ! IN V_1 on P-grid.
,anthrop_heat(land_pts,ntiles)                                    &
!                                  ! IN Additional heat source on tiles
!                                  !    used for anthropgenic urban
!                                  !    heat source (W/m2)

,soil_clay ( row_length, rows )                                   &
                             ! IN Soil clay fraction
,o3(land_pts)                ! IN Surface ozone concentration (ppb).


LOGICAL                                                           &
 ltimer                                                           &
                             ! IN Logical switch for TIMER diags
,l_dust                                                           &
                             ! IN switch for mineral dust
,l_dust_diag                 ! IN Switch for diagnostic mineral dust
                             !    lifting

REAL, INTENT(IN) :: seasalinityfactor
!                                  ! Factor allowing for the effect
!                                  ! of the salinity of sea water
!                                  ! on the evaporative flux.
INTEGER, INTENT(IN) :: cor_ust
!                        ! Switch to use correct friction velocity
!                        ! NOT USED IN JULES
!                        ! Remains here to maintain a consistent
!                        ! subroutine interface for the UM

LOGICAL                                                           &
 lq_mix_bl
LOGICAL :: l_ukca   ! switch for UKCA scheme
!  STASH flags :-

LOGICAL                                                           &
 sfme                                                             &
         ! IN Flag for FME (q.v.).
,sz0heff                                                          &
         ! IN Flag for Z0H_EFF
,sq1p5                                                            &
         ! IN Flag for Q1P5M (q.v.)
,st1p5                                                            &
         ! IN Flag for T1P5M (q.v.)
,su10                                                             &
         ! IN Flag for U10M (q.v.)
,sv10    ! IN Flag for V10M (q.v.)

!  In/outs :-

REAL                                                              &
 z0msea(row_length,rows)                                          &
                             ! INOUT Sea-surface roughness
!                                  !       length for momentum (m).
,z0m_scm(row_length,rows)                                         &
                             ! IN Fixed Sea-surface roughness
!                                  !    length for momentum (m).(SCM)
,z0h_scm(row_length,rows)                                         &
                             ! IN Fixed Sea-surface roughness
!                                  !    length for heat (m). (SCM)
,gs(land_pts)                                                     &
                             ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
,g_leaf_acc(land_pts,npft)                                        &
                             ! INOUT Accumulated G_LEAF
,npp_ft_acc(land_pts_trif,npft_trif)                              &
!                                  ! INOUT Accumulated NPP_FT
,resp_w_ft_acc(land_pts_trif,npft_trif)                           &
!                                  ! INOUT Accum RESP_W_FT
,resp_s_acc(land_pts_trif,dim_cs1) ! INOUT Accumulated RESP_S

!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

INTEGER                                                           &
 tile_index(land_pts,ntype)                                       &
                             ! OUT Index of tile points
,tile_pts(ntype)             ! OUT Number of tile points


!  (a) Calculated anyway (use STASH space from higher level) :-

REAL                                                              &
 cd(row_length,rows)                                              &
                             ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     momentum.
,ch(row_length,rows)                                              &
                             ! OUT Turbulent surface exchange
!                                  !     (bulk transfer) coefficient for
!                                  !     heat and/or moisture.
,recip_l_mo_sea(row_length,rows)                                  &
!                                  ! OUT Reciprocal of the surface
!                                  !     Obukhov  length at sea
!                                  !     points. (m-1).
,e_sea(row_length,rows)                                           &
                             ! OUT Evaporation from sea times
!                                  !     leads fraction. Zero over land.
!                                  !     (kg per square metre per sec).
,fqw_1(row_length,rows)                                           &
                             ! OUT Moisture flux between layers
!                                  !     (kg per square metre per sec).
!                                  !     FQW(,1) is total water flux
!                                  !     from surface, 'E'.
,ftl_1(row_length,rows)                                           &
                             ! OUT FTL(,K) contains net turbulent
!                                  !     sensible heat flux into layer K
!                                  !     from below; so FTL(,1) is the
!                                  !     surface sensible heat, H.(W/m2)
,ftl_tile(land_pts,ntiles)                                        &
                             ! OUT Surface FTL for land tiles
,le_tile(land_pts,ntiles)                                         &
                             ! OUT Surface latent heat flux for
!                                  !     land tiles
,h_sea(row_length,rows)                                           &
                             ! OUT Surface sensible heat flux over
!                                  !     sea times leads fraction (W/m2)
,radnet_sice(row_length,rows)                                     &
                             ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
,radnet_tile(land_pts,ntiles)                                     &
                             ! OUT Surface net radiation on
,rhokm_land(1-off_x:row_length+off_x,1-off_y:rows+off_y)          &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
,rhokm_ssi(1-off_x:row_length+off_x,1-off_y:rows+off_y)           &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
,flandg_u(row_length,rows)                                        &
                             ! OUT Land frac (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
!                                  !     land tiles (W/m2)
,rhokm_1(1-off_x:row_length+off_x,1-off_y:rows+off_y)             &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum on P-grid
,rhokm_u_land(row_length,rows)                                    &
                                ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rhokm_u_ssi(row_length,rows)                                     &
                               ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,flandg_v(row_length,n_rows)                                      &
                             ! OUT Land frac (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rhokm_u_1(row_length,rows)                                       &
                             ! OUT Exchange coefficients for
!                                  !     momentum (on U-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rhokm_v_land(row_length,n_rows)                                  &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rhokm_v_ssi(row_length,n_rows)                                   &
!                                  ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rhokm_v_1(row_length,n_rows)                                     &
                             ! OUT Exchange coefficients for
!                                  !     momentum (on V-grid, with 1st
!                                  !     and last rows undefined or, at
!                                  !     present, set to "missing data")
,rib(row_length,rows)                                             &
                             ! OUT Mean bulk Richardson number for
!                                  !     lowest layer.
,rib_tile(land_pts,ntiles)                                        &
                             ! OUT RIB for land tiles.
,taux_1(row_length,rows)                                          &
                             ! OUT W'ly component of surface wind
,taux_land(row_length,rows)                                       &
                             ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
,taux_ssi(row_length,rows)                                        &
                             ! OUT W'ly component of sfc wind
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
!                                  !     set to missing data
!                                  !     stress (N/sq m). (On U-grid
!                                  !     with first and last rows
!                                  !     undefined or, at present,
,tauy_land(row_length,n_rows)                                     &
                             ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
,tauy_ssi(row_length,n_rows)                                      &
                             ! OUT S'ly component of sfc wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
!                                  !     set to missing data
,tauy_1(row_length,n_rows)                                        &
                             ! OUT S'ly component of surface wind
!                                  !     stress (N/sq m).  On V-grid;
!                                  !     comments as per TAUX.
,rho_cd_modv1(row_length,rows)                                    &
!                                  ! OUT Surface air density * drag coef
!                                  !     *mod(v1 - v0) before interp
,rho_aresist(row_length,rows)                                     &
                             ! OUT RHOSTAR*CD_STD*VSHR for Sulphur
!                                  !     cycle
,aresist(row_length,rows)                                         &
                             ! OUT 1/(CD_STD*VSHR) for Sulphur
!                                  !     cycle
,resist_b(row_length,rows)                                        &
                             ! OUT (1/CH-1/(CD_STD)/VSHR for
!                                  !     Sulphur cycle
,rho_aresist_tile(land_pts,ntiles)                                &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles
,aresist_tile(land_pts,ntiles)                                    &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
,resist_b_tile(land_pts,ntiles)                                   &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles
, r_b_dust(row_length,rows,ndiv)                                  &
                             ! OUT surf layer res for dust
, cd_std_dust(row_length,rows)                                    &
                             ! OUT Bulk transfer coef. for
!                                  ! momentum, excluding orographic effects
, u_s_std_tile(land_pts,ntiles)                                   &
                             ! OUT Surface friction velocity
!                                  !     (standard value)
,emis_tile(land_pts,ntiles)                                       &
                             ! OUT Emissivity for land tiles
,emis_soil(land_pts)                                              &
                             ! OUT Emissivity of underlying soil
,wt_ext(land_pts,sm_levels)                                       &
                             ! OUT cumulative fraction of transp'n
,ra(land_pts)                ! OUT Aerodynamic resistance (s/m).


!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

REAL                                                              &
 fme(row_length,rows)        ! OUT Wind mixing "power" (W/m2).

!-2 Genuinely output, needed by other atmospheric routines :-

REAL                                                              &
 fb_surf(row_length,rows)                                         &
                             ! OUT Surface flux buoyancy over
!                                  !     density (m^2/s^3)
,u_s(row_length,rows)                                             &
                             ! OUT Surface friction velocity (m/s)
,t1_sd(row_length,rows)                                           &
                             ! OUT Standard deviation of turbulent
!                                  !     fluctuations of layer 1 temp;
!                                  !     used in initiating convection.
,q1_sd(row_length,rows)      ! OUT Standard deviation of turbulent
!                                  !     flucs of layer 1 humidity;
!                                  !     used in initiating convection.

REAL                                                              &
 alpha1(land_pts,ntiles)                                          &
                             ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
,alpha1_sice(row_length,rows)                                     &
                             ! OUT ALPHA1 for sea-ice.
,ashtf_prime(row_length,rows)                                     &
                             ! OUT Coefficient to calculate
!                                  !     surface heat flux into soil or
!                                  !     sea-ice.
,ashtf_prime_tile(land_pts,ntiles)                                &
                             ! OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
,fqw_tile(land_pts,ntiles)                                        &
                             ! OUT Surface FQW for land tiles
,epot_tile(land_pts,ntiles)                                       &
                             ! OUT Local EPOT for land tiles.
,rhokpm_pot(land_pts,ntiles)                                      &
!                                  ! OUT Potential evaporation
!                                  !     exchange coeff.
!                                  !     (Note used with JULES)
,fqw_ice(row_length,rows)                                         &
                             ! OUT Surface FQW for sea-ice
,ftl_ice(row_length,rows)                                         &
                             ! OUT Surface FTL for sea-ice
,fraca(land_pts,ntiles)                                           &
                             ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
,rhostar(row_length,rows)                                         &
                             ! OUT Surface air density
,resfs(land_pts,ntiles)                                           &
                             ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
,resft(land_pts,ntiles)                                           &
                             ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
,rhokh(row_length,rows)                                           &
                             ! OUT Grid-box surface exchange
!                                  !     coefficients
,rhokh_tile(land_pts,ntiles)                                      &
                             ! OUT Surface exchange coefficients
!                                  !     for land tiles
,rhokh_sice(row_length,rows)                                      &
                             ! OUT Surface exchange coefficients
!                                  !     for sea and sea-ice
,rhokpm(land_pts,ntiles)                                          &
                             ! OUT Land surface exchange coeff.
!                                  !     (Note used with JULES)
,rhokpm_sice(row_length,rows)                                     &
                             ! OUT Sea-ice surface exchange coeff.
!                                  !     (Note used with JULES)
,rhokh_mix(row_length,rows)                                       &
                             ! OUT Exchange coeffs for moisture.
,dtstar_tile(land_pts,ntiles)                                     &
                             ! OUT Change in TSTAR over timestep
!                                  !     for land tiles
,dtstar(row_length,rows)                                          &
                             ! OUT Change is TSTAR over timestep
!                                  !     for sea-ice
,h_blend_orog(row_length,rows)                                    &
!                                  ! OUT Blending height used as part of
!                                  !     effective roughness scheme
,z0hssi(row_length,rows)                                          &
                             ! OUT Roughness length for heat and
!                                  !     moisture over sea (m).
,z0mssi(row_length,rows)                                          &
                             ! OUT Roughness length for momentum
!                                  !     over sea (m).
,z0h_tile(land_pts,ntiles)                                        &
                             ! OUT Tile roughness lengths for heat
!                                  !     and moisture (m).
,z0h_eff(row_length,rows)                                         &
                             ! OUT Effective grid-box roughness
!                                  !     length for heat, moisture (m)
,z0m_gb(row_length,rows)                                          &
                             ! OUT Gridbox mean roughness length
!                                  !     for momentum (m).
,z0m_tile(land_pts,ntiles)                                        &
                             ! OUT Tile roughness lengths for
!                                  !     momentum.
,z0m_eff(row_length,rows)                                         &
                             ! OUT Effective grid-box roughness
!                                  !     length for momentum
,cdr10m_u(row_length,rows)                                        &
                             ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     U-grid; comments as per RHOKM.
,cdr10m_v(row_length,n_rows)                                      &
                             ! OUT Ratio of CD's reqd for
!                                  !     calculation of 10 m wind. On
!                                  !     V-grid; comments as per RHOKM.
,chr1p5m(land_pts,ntiles)                                         &
                             ! OUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
,chr1p5m_sice(row_length,rows)                                    &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
,smc(land_pts)                                                    &
                             ! OUT Available moisture in the
!                                  !     soil profile (mm).
,hcons(land_pts)                                                  &
                             ! OUT Soil thermal conductivity
!                                  !     including water and ice
,vshr(row_length,rows)                                            &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,vshr_land(row_length,rows)                                       &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,vshr_ssi(row_length,rows)                                        &
                             ! OUT Magnitude of surface-to-lowest
!                                  !     atm level wind shear (m per s).
,gpp(land_pts)                                                    &
                             ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
,npp(land_pts)                                                    &
                             ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
,resp_p(land_pts)                                                 &
                             ! OUT Plant respiration (kg C/m2/s).
,g_leaf(land_pts,npft)                                            &
                             ! OUT Leaf turnover rate (/360days).
,gpp_ft(land_pts,npft)                                            &
                             ! OUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
,npp_ft(land_pts,npft)                                            &
                             ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
,resp_p_ft(land_pts,npft)                                         &
                             ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
,resp_s(land_pts,dim_cs1)                                         &
                          ! OUT Soil respiration (kg C/m2/s).
,resp_s_tot(dim_cs2)                                              &
                            ! OUT Total soil respiration
                            ! (kg C/m2/s).
,resp_w_ft(land_pts,npft)                                         &
                             ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
,gc(land_pts,ntiles)                                              &
                             ! OUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
,canhc_tile(land_pts,ntiles)                                      &
                             ! IN Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
,wt_ext_tile(land_pts,sm_levels,ntiles)                           &
!                                  ! IN Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
,flake(land_pts,ntiles)                                           &
                             ! IN Lake fraction.
,tile_frac(land_pts,ntiles)                                       &
                             ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
,fsmc(land_pts,npft)                                              &
                            ! OUT Moisture availability factor.
,flux_o3_ft(land_pts,npft)                                        &
                            ! OUT Flux of O3 to stomata (nmol O3/m2/s).
,fo3_ft(land_pts,npft)      ! OUT Ozone exposure factor.

!---------------------------------------------------------------------
!  External routines called :-

EXTERNAL tilepts,physiol,heat_con,sf_exch,snowtherm
#if defined(UM_RUN)
EXTERNAL swap_bounds,fill_external_halos,p_to_u,p_to_v
#endif
EXTERNAL timer


!-----------------------------------------------------------------------

!  Workspace :-

REAL :: work_clay         ! working variable

REAL                                                              &
 flandg_x(1-off_x:row_length+off_x,1-off_y:rows+off_y)            &
!                               ! Land fraction on P-grid
!                               ! the effects of water and ice (W/m2)
,cdr10m(1-off_x:row_length+off_x,1-off_y:rows+off_y)              &
!                               ! Ratio of CD's reqd for calculation
!                               ! of 10 m wind. On P-grid
,vfrac_tile(land_pts,ntiles)                                      &
!                               ! Fractional canopy coverage for
!                               ! land tiles.
,csnow(land_pts,nsmax)                                            &
                          ! Areal heat capacity of snow (J/K/m2)
,ksnow(land_pts,nsmax)                                            &
                          ! Thermal conductivity of snow (W/m/K)
,hcons_snow(land_pts,ntiles)                                      &
!                               ! Snow thermal conductivity
,resp_frac(dim_cs2)                                               &
                          ! respired fraction of RESP_S
,clay_land(dim_cs2)       ! Clay fraction on land points


!  Expanded wind arrays - only needed for Buddy_sea coastal tiling
REAL                                                              &
 u_1_px(1-off_x:row_length+off_x,1-off_y:rows+off_y)              &
,v_1_px(1-off_x:row_length+off_x,1-off_y:rows+off_y)              &
,u_0_px(1-off_x:row_length+off_x,1-off_y:rows+off_y)              &
,v_0_px(1-off_x:row_length+off_x,1-off_y:rows+off_y)

!     Factors used to weight windspeeds over coastal points
REAL                                                              &
  flandfac(row_length,rows)                                       &
 ,flandfac_u(row_length,rows)                                     &
 ,flandfac_v(row_length,n_rows)                                   &
 ,flandfac_x(1-off_x:row_length+off_x,1-off_y:rows+off_y)         &
 ,fseafac(row_length,rows)                                        &
 ,fseafac_u(row_length,rows)                                      &
 ,fseafac_v(row_length,n_rows)                                    &
 ,fseafac_x(1-off_x:row_length+off_x,1-off_y:rows+off_y)


!  Local scalars :-

INTEGER                                                           &
 i,j,k,l,n                                                        &
            ! LOCAL Loop counter (horizontal field index).
,is,js                                                            &
            ! Loop counter for coastal point stencil
,COUNT      ! Counter for average wind speed

REAL                                                              &
 ushear                                                           &
              ! U-component of surface-to-lowest-level wind shear.
,vshear                                                           &
              ! V-component of surface-to-lowest-level wind shear.
,vshr2        ! Square of magnitude of surface-to-lowest-level
!                   ! wind shear.

REAL seawind  ! average wind speed adjacent to coast
REAL fseamax  ! Maximum factor to apply to coast wind speed

     ! Minimum factor allowed to convert coastal wind speed to land part
REAL flandmin
PARAMETER(flandmin=0.2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SF_EXPL_L',zhook_in,zhook_handle)
IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('SF_EXPL ',3)
END IF

!-----------------------------------------------------------------------
! Call TILEPTS to calculate TILE_PTS and TILE_INDEX for surface types
!-----------------------------------------------------------------------
! DEPENDS ON: tilepts
CALL tilepts(land_pts,frac,tile_pts,tile_index)

!-----------------------------------------------------------------------
! Calculate wind shear between level 1 and the surface
!-----------------------------------------------------------------------

DO j=1,rows
  DO i=1,row_length
  IF(flandg(i,j) <  1.0)THEN
    ushear = u_1_p(i,j) - u_0_p(i,j)
    vshear = v_1_p(i,j) - v_0_p(i,j)
    vshr2 = MAX (1.0e-6 , ushear*ushear + vshear*vshear)
    vshr_ssi(i,j) = SQRT(vshr2)
  ELSE
    vshr_ssi(i,j) = 0.0
  END IF

  IF(flandg(i,j) >  0.0)THEN
  vshr2 = MAX (1.0e-6 , u_1_p(i,j)*u_1_p(i,j)                     &
    + v_1_p(i,j)*v_1_p(i,j))
  vshr_land(i,j) = SQRT(vshr2)
  ELSE
    vshr_land(i,j) = 0.0
  END IF

  vshr(i,j)= flandg(i,j)*vshr_land(i,j)                           &
    + (1.0 - flandg(i,j))*vshr_ssi(i,j)
  END DO
END DO

#if !defined(SCMA)

DO j= 1, rows
  DO i= 1, row_length
    flandg_x(i,j)  = flandg(i,j)
  END DO
END DO

! DEPENDS ON: swap_bounds
CALL swap_bounds(flandg_x, row_length, rows, 1,                   &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: fill_external_halos
CALL fill_external_halos(flandg_x,row_length,rows,1,              &
                         off_x,off_y)

IF (l_ctile .AND. buddy_sea == on) THEN

  DO j= 1, rows
    DO i= 1, row_length
      u_1_px(i,j) = u_1_p(i,j)
      v_1_px(i,j) = v_1_p(i,j)
      u_0_px(i,j) = u_0_p(i,j)
      v_0_px(i,j) = v_0_p(i,j)
    END DO
  END DO

! DEPENDS ON: swap_bounds
  CALL swap_bounds(u_1_px, row_length, rows, 1,                   &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: swap_bounds
  CALL swap_bounds(v_1_px, row_length, rows, 1,                   &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: swap_bounds
  CALL swap_bounds(u_0_px, row_length, rows, 1,                   &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: swap_bounds
  CALL swap_bounds(v_0_px, row_length, rows, 1,                   &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(u_1_px,row_length,rows,1,              &
                         off_x,off_y)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(v_1_px,row_length,rows,1,              &
                         off_x,off_y)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(u_0_px,row_length,rows,1,              &
                         off_x,off_y)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(v_0_px,row_length,rows,1,              &
                         off_x,off_y)

  DO j=1,rows
   DO i=1,row_length
    fseafac(i,j)  = 1.0
    flandfac(i,j) = 1.0

    IF ( flandg(i,j) > 0.01 .AND. flandg(i,j) < 0.99 ) THEN
!           !-----------------------------------------------------
!           ! Calculate average windspeed over adjacent sea points
!           !-----------------------------------------------------
      seawind=0.0
      COUNT = 0
      DO is=i-1,i+1
      DO js=j-1,j+1
        IF( flandg_x(is,js) < 0.001 )THEN
!               ! ie. this is basically a sea point
          ushear = u_1_px(is,js) - u_0_px(is,js)
          vshear = v_1_px(is,js) - v_0_px(is,js)
          vshr2 = MAX (1.0e-10 , ushear*ushear + vshear*vshear)
          seawind = seawind + SQRT( vshr2 )
          COUNT = COUNT + 1
        END IF
      END DO
      END DO
!           !-----------------------------------------------------
!           ! Calculate multiplicative factor, FSEAFAC, to convert
!           ! from the GBM VSHR to an appropriate marine VSHR
!           !-----------------------------------------------------
      IF (COUNT > 0) THEN
        seawind = seawind/FLOAT(COUNT)
!             ! Restrict FSEAFAC so FLANDFAC>FLANDMIN
        fseamax = MIN( 1./flandmin,                               &
                      (1.-flandmin*flandg(i,j))/(1.-flandg(i,j)) )
!             ! First limit is to keep fseamax sensible as FLANDG -> 1
!             ! Second limit is to keep fland > flandmin, remembering
!             !   that the we want FLANDG-weighted sum of factors =1
!             !   to preserve the gridbox mean VSHR
        fseafac(i,j) = MAX(1.0,                                   &
                       MIN( fseamax, seawind/vshr(i,j) ))
      END IF

      vshr_ssi(i,j) = vshr(i,j) * fseafac(i,j)

      flandfac(i,j) = ( 1.0 - fseafac(i,j)*(1.0-flandg(i,j)) )    &
                    / flandg(i,j)
      vshr_land(i,j) = vshr(i,j) * flandfac(i,j)


      vshr(i,j)= flandg(i,j)*vshr_land(i,j)                       &
        + (1.0 - flandg(i,j))*vshr_ssi(i,j)
      END IF
  END DO
  END DO

END IF  ! test on buddy_sea switch
#endif


! This IF-test is a fix, but is it the right fix?
IF (l_triffid) THEN
   DO l = 1,land_pts
     j = (land_index(l)-1)/row_length + 1
     i = land_index(l) - (j-1)*row_length
     clay_land(l) = soil_clay(i,j)
! Set default value for missing data points
     IF (clay_land(l) == rmdi) clay_land(l)=0.23
   END DO !LAND_PTS
END IF



!-----------------------------------------------------------------------
! Call MOSES II physiology routine to calculate surface conductances
! and carbon fluxes.
!-----------------------------------------------------------------------

! DEPENDS ON: physiol
CALL physiol (                                                    &
 row_length,rows,land_pts,land_index,                             &
 sm_levels,ntiles,tile_pts,tile_index,l_aggregate,                &
 dim_cs1, dim_cs2,                                                &
 co2_mmr,co2_3d,co2_dim_len,co2_dim_row,l_co2_interactive,        &
 l_triffid, l_q10,                                                &
 can_model,cs,frac,canht_ft,photosynth_act_rad,                   &
 lai_ft,pstar,qw_1,sthu,t_soil(:,1),tstar_tile,                   &
 smvccl,smvcst,smvcwt,vshr,z0_tile,z1_uv,o3,                      &
 canhc_tile,vfrac_tile,emis_tile,emis_soil,flake,                 &
 g_leaf,gs,gc,gpp,gpp_ft,npp,npp_ft,                              &
 resp_p,resp_p_ft,resp_s,resp_w_ft,smc,wt_ext_tile,fsmc,          &
 wt_ext,ra,albsoil,cos_zenith_angle                               &
,can_rad_mod,ilayers,flux_o3_ft,fo3_ft)


!----------------------------------------------------------------------
! If TRIFFID is being used apply any correction to the land-atmosphere
! fluxes on the first timestep after the last TRIFFID call. Such a
! correction will typically be associated with a total depletion of
! carbon or with maintanence of the seed fraction. The corrections
! are stored in the accumulation variables after the call to TRIFFID.
! The correction is added to the instantaneous land-atmosphere fluxes
! (so that the atmospheric carbon budget is corrected) but is not
! included in the accumulation variables which drive TRIFFID, since
! this has already been dealt with during the last TRIFFID call.
!----------------------------------------------------------------------
IF (l_triffid .AND.(asteps_since_triffid==1)                      &
    .AND. cycleno==numcycles) THEN
  DO n=1,npft
    DO l=1,land_pts
      npp_ft(l,n)=npp_ft(l,n)+npp_ft_acc(l,n)/timestep
      resp_p_ft(l,n)=resp_p_ft(l,n)-npp_ft_acc(l,n)/timestep
      npp_ft_acc(l,n)=-npp_ft_acc(l,n)
    END DO
  END DO
  DO n=1,dim_cs1
    DO l=1,land_pts
      resp_s(l,n)=resp_s(l,n)+resp_s_acc(l,n)/timestep
      resp_s_acc(l,n)=-resp_s_acc(l,n)
    END DO
  END DO
END IF

!----------------------------------------------------------------------
! Increment accumulation of leaf turnover rate.
! This is required for leaf phenology and/or TRIFFID, either of
! which can be enabled independently of the other.
!----------------------------------------------------------------------
IF ( cycleno == numcycles ) THEN
IF (l_phenol.OR.l_triffid) THEN
  DO n=1,npft
    DO l=1,land_pts
      g_leaf_acc(l,n) = g_leaf_acc(l,n) +                         &
      g_leaf(l,n)*(timestep/secs_per_360days)
    END DO
  END DO
END IF

!----------------------------------------------------------------------
! Increment accumulation prognostics for TRIFFID
!----------------------------------------------------------------------
IF (l_triffid) THEN
  DO n=1,npft
    DO l=1,land_pts
      npp_ft_acc(l,n) = npp_ft_acc(l,n) + npp_ft(l,n)*timestep
      resp_w_ft_acc(l,n) = resp_w_ft_acc(l,n)                     &
                                      + resp_w_ft(l,n)*timestep
    END DO
  END DO
  DO n=1,dim_cs1
    DO l=1,land_pts
      resp_s_acc(l,n)=resp_s_acc(l,n)+resp_s(l,n)*timestep
    END DO
  END DO
END IF
END IF ! CycleNo == NumCycles

!-----------------------------------------------------------------------
! calculate CO2:(BIO+HUM) ratio, dependent on soil clay content, and
! sum soil respiration components
! (RESP_FRAC here then contains the fraction of soil respiration which
! is respired to the atmos. the rest is re-partitioned into BIO+HUM)

! RESP_S_ACC contains the full amount, and this is carried forward to
! VEG_CTL for use in updating soil carbon pools. RESP_S_TOT calculated
! here is passed to BL_TRMIX as the fraction which is respired as CO2
! to the atmosphere. RESP_S_TOT, and RESP_S are also passed out for
! storage in diagnostics 3293, and 3467-470.

!-----------------------------------------------------------------------
IF (l_triffid) THEN
  DO i=1,land_pts
    work_clay = EXP(-0.0786 * 100.0*clay_land(i))
    resp_frac(i) = (3.0895 + 2.672*work_clay) /                &
                   (4.0895 + 2.672*work_clay)
    resp_s(i,1)  = resp_s(i,1) * resp_frac(i)
    resp_s(i,2)  = resp_s(i,2) * resp_frac(i)
    resp_s(i,3)  = resp_s(i,3) * resp_frac(i)
    resp_s(i,4)  = resp_s(i,4) * resp_frac(i)
    resp_s_tot(i) = resp_s(i,1) + resp_s(i,2) +                   &
                    resp_s(i,3) + resp_s(i,4)
  END DO
END IF
!-----------------------------------------------------------------------
! Reset TILE_PTS and TILE_INDEX and set tile fractions to 1 if aggregate
! tiles are used (L_AGGREGATE=.T.).
! Otherwise, set tile fractions to surface type fractions.
!-----------------------------------------------------------------------
IF (l_aggregate) THEN
  tile_pts(1) = land_pts
  DO l=1,land_pts
    tile_frac(l,1) = 1.
    tile_index(l,1) = l
  END DO
ELSE
  DO n=1,ntype
    DO l = 1, land_pts
      tile_frac(l,n) = frac(l,n)
    END DO
  END DO
END IF

IF (land_pts >  0) THEN    ! Omit if no land points

!-----------------------------------------------------------------------
! Calculate the thermal conductivity of the top soil layer.
!-----------------------------------------------------------------------
! DEPENDS ON: heat_con
  CALL heat_con (land_pts,hcon,sthu(:,1),sthf(:,1),               &
                 smvcst(:,1),hcons,ltimer)

! Thermal conductvity of top snow layer if nsmax > 0
  IF(nsmax > 0) THEN
    DO n=1,ntiles
! DEPENDS ON: snowtherm
      CALL snowtherm(land_pts,tile_pts(n),nsnow(:,n),             &
                     tile_index(:,n),ds(:,n,:),sice(:,n,:),       &
                     sliq(:,n,:),csnow,ksnow)
      DO l=1,land_pts
        hcons_snow(l,n) = ksnow(l,1)
      END DO
    END DO
  END IF

END IF                     ! End test on land points

!-----------------------------------------------------------------------
!! Calculate net radiation on land tiles and sea-ice
!-----------------------------------------------------------------------

    radnet_tile(:,:) = 0.
    le_tile(:,:) = 0.

IF (l_skyview) THEN
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/row_length + 1
      i = land_index(l) - (j-1)*row_length
      radnet_tile(l,n) = sw_tile(l,n) + emis_tile(l,n)*           &
        sky(i,j)*( lw_down(i,j) - sbcon*tstar_tile(l,n)**4 )
    END DO
  END DO
ELSE
  DO n=1,ntiles
    DO k=1,tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/row_length + 1
      i = land_index(l) - (j-1)*row_length
      radnet_tile(l,n) = sw_tile(l,n) + emis_tile(l,n)*           &
                 ( lw_down(i,j) - sbcon*tstar_tile(l,n)**4 )
    END DO
  END DO
END IF

DO j=1,rows
 DO i=1,row_length
  radnet_sice(i,j) = 0.
  IF (flandg(i,j) <  1.0 .AND. ice_fract(i,j) >  0.)              &
    radnet_sice(i,j) = rad_sice(i,j) -                            &
                       sbcon*tstar_sice(i,j)**4
 END DO
END DO

!-----------------------------------------------------------------------
!! 4.  Surface turbulent exchange coefficients and "explicit" fluxes
!!     (P243a, routine SF_EXCH).
!!     Wind mixing "power" and some values required for other, later,
!!     diagnostic calculations, are also evaluated if requested.
!-----------------------------------------------------------------------

! DEPENDS ON: sf_exch
CALL sf_exch (                                                    &
 row_length,rows,off_x,off_y,halo_i,halo_j,                       &
 land_pts,ntiles,land_index,                                      &
 tile_index,tile_pts,fland,flandg,                                &
 ssi_pts,sea_pts,sice_pts,                                        &
 ssi_index,sea_index,sice_index,fssi,sea_frac,sice_frac,          &
 nsmax,nsnow,ds,hcons_snow,                                       &
 bq_1,bt_1,canhc_tile,canopy,catch,dzsoil(1),flake,gc,hcons,      &
 can_model,catch_snow,lq_mix_bl,                                  &
 ho2r2_orog,ice_fract,snowdepth,snow_tile,pstar,qw_1,radnet_sice, &
 radnet_tile,sil_orog_land,smvcst(:,1),tile_frac,timestep,        &
 surf_hgt,emis_tile,emis_soil,                                    &
 tl_1,ti,t_soil(:,1),                                             &
 tsnow,                                                           &
 tstar_tile,tstar_land,tstar_sea,tstar_sice,tstar_ssi,z_land,     &
 l_ctile,seasalinityfactor,                                       &
 tstar,l_aggregate,l_spec_z0,z0m_scm,z0h_scm,l_dust,l_dust_diag,  &
 vfrac_tile,vshr_land,vshr_ssi,zh,ddmfx,                          &
 z0_tile,z1_uv,z1_uv_top,z1_tq,z1_tq_top,land_mask,               &
 su10,sv10,sq1p5,st1p5,sfme,sz0heff,ltimer,formdrag,fd_stab_dep,  &
 orog_drag_param,z0msea,                                          &
 alpha1,alpha1_sice,ashtf_prime,ashtf_prime_tile,cd,ch,           &
 recip_l_mo_sea,cdr10m,chr1p5m,                                   &
 chr1p5m_sice,e_sea,fme,fqw_1,fqw_tile,epot_tile,fqw_ice,         &
 ftl_1,ftl_tile,ftl_ice,fraca,h_blend_orog,h_sea,charnock,        &
 rhostar,resfs,resft,rib,rib_tile,                                &
 fb_surf,u_s,q1_sd,t1_sd,z0hssi,z0h_tile,z0h_eff,                 &
 z0m_gb,z0mssi,z0m_tile,z0m_eff,rho_aresist,aresist,resist_b,     &
 rho_aresist_tile,aresist_tile,resist_b_tile,                     &
 r_b_dust,cd_std_dust,u_s_std_tile,                               &
 rho_cd_modv1,rhokh_tile,rhokh_sice,rhokm_1,rhokm_land,rhokm_ssi, &
 dtstar_tile,dtstar,rhokh,rhokh_mix,anthrop_heat                  &
 )




#if !defined(SCMA)

! DEPENDS ON: swap_bounds
CALL swap_bounds(rhokm_1, row_length, rows,                       &
                 1 , off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: swap_bounds
CALL swap_bounds(rhokm_land, row_length, rows,                    &
                 1 , off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: swap_bounds
CALL swap_bounds(rhokm_ssi, row_length, rows,                     &
                 1 , off_x, off_y, fld_type_p, .FALSE.)

! DEPENDS ON: fill_external_halos
CALL fill_external_halos(rhokm_1,row_length,rows,1,               &
                         off_x,off_y)
! DEPENDS ON: fill_external_halos
CALL fill_external_halos(rhokm_land,row_length,rows,1,            &
                         off_x,off_y)
! DEPENDS ON: fill_external_halos
CALL fill_external_halos(rhokm_ssi,row_length,rows,1,             &
                         off_x,off_y)

! DEPENDS ON: p_to_u
CALL p_to_u(rhokm_1,row_length,rows,1,                            &
            off_x, off_y, rhokm_u_1)
! DEPENDS ON: p_to_u_land
CALL p_to_u_land(rhokm_land,flandg_x,row_length,rows,1,           &
            off_x, off_y, rhokm_u_land)
! DEPENDS ON: p_to_u_sea
CALL p_to_u_sea(rhokm_ssi,flandg_x,row_length,rows,1,             &
            off_x, off_y, rhokm_u_ssi)

! DEPENDS ON: p_to_v
CALL p_to_v(rhokm_1,row_length, rows, n_rows,                     &
            1, off_x, off_y, rhokm_v_1)
! DEPENDS ON: p_to_v_land
CALL p_to_v_land(rhokm_land,flandg_x,row_length, rows, n_rows,    &
            1, off_x, off_y, rhokm_v_land)
! DEPENDS ON: p_to_v_sea
CALL p_to_v_sea(rhokm_ssi,flandg_x,row_length, rows, n_rows,      &
            1, off_x, off_y, rhokm_v_ssi)

! DEPENDS ON: p_to_u
CALL p_to_u(flandg_x,row_length,rows,1,                           &
            off_x, off_y, flandg_u)
! DEPENDS ON: p_to_v
CALL p_to_v(flandg_x,row_length, rows, n_rows,                    &
            1, off_x, off_y, flandg_v)

IF (l_ctile .AND. buddy_sea == on) THEN
!       ! Interpolate wind speed factors to u and v columns

  DO j= 1, rows
  DO i= 1, row_length
    flandfac_x(i,j)= flandfac(i,j)
    fseafac_x(i,j) = fseafac(i,j)
  END DO
  END DO

! DEPENDS ON: swap_bounds
  CALL swap_bounds(flandfac_x, row_length, rows, 1,               &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(flandfac_x,row_length,rows,1,          &
                         off_x,off_y)
! DEPENDS ON: p_to_u
  CALL p_to_u(flandfac_x,row_length,rows,1,                       &
                         off_x, off_y, flandfac_u)
! DEPENDS ON: p_to_v
  CALL p_to_v(flandfac_x,row_length, rows, n_rows, 1,             &
                         off_x, off_y, flandfac_v)
! DEPENDS ON: swap_bounds
  CALL swap_bounds(fseafac_x, row_length, rows, 1,                &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: fill_external_halos
  CALL fill_external_halos(fseafac_x,row_length,rows,1,           &
                         off_x,off_y)
! DEPENDS ON: p_to_u
  CALL p_to_u(fseafac_x,row_length,rows,1,                        &
                         off_x, off_y, fseafac_u)
! DEPENDS ON: p_to_v
  CALL p_to_v(fseafac_x,row_length, rows, n_rows, 1,              &
                         off_x, off_y, fseafac_v)

     END IF  ! test on Buddy_sea

#else
DO j= 1, rows
 DO i= 1, row_length
   rhokm_u_1(i,j) = rhokm_1(i,j)
   rhokm_u_land(i,j) = rhokm_land(i,j)
   rhokm_v_land(i,j) = rhokm_land(i,j)
   rhokm_u_ssi(i,j) = rhokm_ssi(i,j)
   rhokm_v_ssi(i,j) = rhokm_ssi(i,j)
   rhokm_v_1(i,j) = rhokm_1(i,j)
   flandg_u(i,j) = flandg(i,j)
   flandg_v(i,j) = flandg(i,j)
   flandfac_u(i,j) = 1.0
   flandfac_v(i,j) = 1.0
   fseafac_u(i,j)  = 1.0
   fseafac_v(i,j)  = 1.0
 END DO
END DO
#endif

IF (l_ctile .AND. buddy_sea == on) THEN

  DO j=1,rows
  DO i=1,row_length
    taux_land(i,j) = rhokm_u_land(i,j)* u_1(i,j) * flandfac_u(i,j)
    taux_ssi(i,j)  = rhokm_u_ssi(i,j) * ( u_1(i,j) - u_0(i,j) )   &
                                      * fseafac_u(i,j)
    taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                    &
                  + (1.-flandg_u(i,j))*taux_ssi(i,j)
  END DO
  END DO

  DO j=1,n_rows
  DO i=1,row_length
    tauy_land(i,j) = rhokm_v_land(i,j)* v_1(i,j) * flandfac_v(i,j)
    tauy_ssi(i,j)  = rhokm_v_ssi(i,j) * ( v_1(i,j) - v_0(i,j) )   &
                                      * fseafac_v(i,j)
    tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                    &
                  + (1.-flandg_v(i,j))*tauy_ssi(i,j)
  END DO
  END DO

ELSE   ! Standard code

  DO j=1,rows
  DO i=1,row_length
    taux_land(i,j) = rhokm_u_land(i,j) * u_1(i,j)
    taux_ssi(i,j) = rhokm_u_ssi(i,j) * ( u_1(i,j) - u_0(i,j) )

    taux_1(i,j) = flandg_u(i,j)*taux_land(i,j)                    &
                  + (1.-flandg_u(i,j))*taux_ssi(i,j)
  END DO
  END DO

  DO j=1,n_rows
  DO i=1,row_length
   tauy_land(i,j) = rhokm_v_land(i,j) * v_1(i,j)
   tauy_ssi(i,j) = rhokm_v_ssi(i,j) * ( v_1(i,j) - v_0(i,j) )

   tauy_1(i,j) = flandg_v(i,j)*tauy_land(i,j)                     &
                 + (1.-flandg_v(i,j))*tauy_ssi(i,j)
  END DO
  END DO

END IF


#if !defined(SCMA)
IF (su10 .OR. sv10) THEN
! DEPENDS ON: swap_bounds
  CALL swap_bounds(                                               &
                 cdr10m, row_length, rows,                        &
                 1, off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: fill_external_halos
CALL fill_external_halos(cdr10m,row_length,rows,1,off_x,off_y)
END IF

IF (su10)THEN
! DEPENDS ON: p_to_u
  CALL p_to_u(cdr10m,row_length,rows,1,                           &
              off_x,off_y,                                        &
              cdr10m_u)

END IF

IF (sv10)THEN
! DEPENDS ON: p_to_v
  CALL p_to_v(cdr10m,row_length,rows,n_rows,1,                    &
            off_x,off_y, cdr10m_v)
END IF
#else
 DO j= 1, rows
  DO i= 1, row_length
    cdr10m_u(i,j) = cdr10m(i,j)
    cdr10m_v(i,j) = cdr10m(i,j)
  END DO
 END DO
#endif


IF (ltimer) THEN
! DEPENDS ON: timer
  CALL timer('SF_EXPL ',4)
END IF

IF (lhook) CALL dr_hook('SF_EXPL_L',zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_expl
#endif
