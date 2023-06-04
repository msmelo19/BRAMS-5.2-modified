! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Module containing switches imported directly by JULES subroutines
! (as opposed to receiving the switches as arguments)
!

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE switches

IMPLICIT NONE

  LOGICAL ::                                                      &
   l_point_data  = .FALSE.                                        &
                     ! Switch for using point rainfall data
  ,l_snow_albedo = .FALSE.                                        &
                     ! Switch for prognostic snow albedo
  ,l_phenol      = .FALSE.                                        &
                     ! Switch for phenology
  ,l_triffid     = .FALSE.                                        &
                     ! Switch for TRIFFID
  ,l_vg_soil     = .FALSE.                                        &
                     ! Switch for using Van Genuchten soil scheme
  ,l_aggregate   = .FALSE.                                        &
                     ! Switch for setting an aggregate surface
                     ! scheme
  ,l_360         = .FALSE.                                        &
                     ! Switch for setting a 360 day year
  ,l_um_jules    = .FALSE.                                        &
                     ! Switch for using JULES in the UM
  ,l_o3_damage   = .FALSE.                                        &
                     ! Switch for ozone damage
  ,l_veg_compete = .FALSE.                                        &
                     ! Switch for competing vegetation
                     ! Setting l_triffid = .TRUE. and leaving this
                     ! as .FALSE. means that the carbon pools evolve
                     ! but the PFT distribution does not change
  ,l_epot_corr   = .FALSE.                                        &
                     ! Switch for using a correction to the
                     ! calculation of potential evaporation
  ,l_snowdep_surf= .FALSE.
                     ! use equivalent canopy snow depth for surface
                     ! calculations on tiles with a snow canopy

  INTEGER ::                                                      &
   can_model                                                      &
!                            switch for thermal vegetation
  ,soilhc_method       = 1                                        &
!                            switch for the calculation method
!                            of soil thermal conductivity
!        SOILHC_METHOD=1: Method of Cox et al (1999).
!        SOILHC_METHOD=2: Simplified Johansen (1975).
  ,i_modiscopt         = 0                                        &
!                            Method of discretization 
!                            in the surface layer 
  ,frac_snow_subl_melt = 0                                        &
!                            switch for use of snow-cover
!                            fraction in the calculation of
!                            sublimation and melting
!                            0 = off
!                            1 = on
  ,all_tiles           = 0                                        &
!                            switch for doing calculations
!                            of tile properties on all tiles
!                            for all gridpoints even when the
!                            tile fraction is zero
!                            (except for land ice).
  ,can_rad_mod         = 4                                        &
!                            Canopy radiation model
  ,cor_mo_iter         = 1                                        &
!                            Switch for MO iteration correction
  ,iseaz0t             = 0                                        &
!                            Switch for the definition of
!                            the thermal roughness length over the sea.
  ,buddy_sea           = 0                                        &
!                            Switch to use the wind speed from
!                            adjacent sea points for the
!                            sea part of coastal grid points
  ,iscrntdiag          = 0
!                            Method of diagnosing the screen temperature

  REAL ::                                                         &
    dz_pdm = 1.0                                                  &             
!                            Soil layer thickness for PDM (m):
   ,b_pdm = 1.0
!                            Shape factor for PDM:

!----------------------------------------------------------------
! Switch for IMOGEN (never changed from default in the UM)
!----------------------------------------------------------------
  LOGICAL ::                                                     &
   l_imogen = .FALSE.


#if !defined(UM_RUN)
!-----------------------------------------------------------------------
! START OF NON-UM SECTION
!
! If we are not in UM, define everything else that is needed
!-----------------------------------------------------------------------

  INTEGER ::                                                      &
   ilayers      ! Number of layers for canopy radiation model

  LOGICAL ::                                                      &
   route                                                          &
                !  Switch for runoff routing.
  ,routeonly                                                      &
                !  Switch to only do runoff routing, nothing else
  ,l_spec_albedo                                                  &
                ! .TRUE. for spectral albedos
                ! and prognostic snow albedo
  ,l_ssice_albedo                                                 &
                ! Switch for including the effect of
                ! snow on the sea-ice albedo
  ,l_trif_eq                                                      &
                ! Switch for vegetation equilibrium
  ,l_top                                                          &
                ! Switch for TOPMODEL-based hydrology
  ,l_pdm                                                          &
                ! Switch for PDM hydrology
  ,l_co2_interactive                                              &
                ! Switch for 3D CO2 field
  ,l_dust                                                         &
                ! Switch for mineral dust
  ,l_z0_orog                                                      &
                ! T to use orog.roughness in surface calcs
  ,l_neg_tstar                                                    &
                ! Switch for -ve TSTAR error check
  ,ltimer                                                         &
                ! Switch for timing information.
  ,l_cosz       ! Switch for turning on calculations of cosz

!-----------------------------------------------------------------------
! Switches introduced during reconciliation with UM that can be set
! via the control file but are passed down in the UM
!-----------------------------------------------------------------------
  LOGICAL ::                                                      &
   l_soil_sat_down     =.FALSE.                                   &
                ! Direction of super_saturated soil moisture
!                 TRUE excess water is pushed down
!                 FALSE excess water is pushed up (as in JULES2.0)
  ,l_q10               = .TRUE.                                   &
                ! True if using Q10 for soil resp
  ,l_anthrop_heat_src  = .FALSE.
                ! Switch for anthropogenic heat source on urban
                ! tile
!-----------------------------------------------------------------------
! Switches introduced during reconciliation with UM that cannot be
! changed via the control file
!-----------------------------------------------------------------------
  LOGICAL ::                                                      &
   l_mod_barker_albedo = .FALSE.                                  &
                ! Use modified Barker albedo (open sea).
  ,l_sice_meltponds    = .FALSE.                                  &
                ! switch for seaice albedo meltponds
  ,l_sice_scattering   = .FALSE.                                  &
                ! switch for seaice albedo internal scatter
  ,lq_mix_bl           = .FALSE.                                  &
                ! TRUE if mixing ratios used in
                ! boundary layer code
  ,l_ctile             = .FALSE.                                  &
                ! True if coastal tiling
  ,l_spec_z0           = .FALSE.                                  &
                ! T if using prescribed sea surface
                ! roughness lengths
  ,l_sice_heatflux     = .FALSE.                                  &
                ! T: semi-implicit update of TI
  ,l_inland            = .FALSE.
                ! True if re-routing inland basin flow
                ! to soil moisture

  INTEGER ::                                                      &
   formdrag      = 0                                              &
                   ! Switch for orographic drag
  ,fd_stab_dep   = 0
                   ! Switch to implement stability
                   ! dependence of orog form drag
                   
!-----------------------------------------------------------------------
! END OF NON-UM SECTION
!-----------------------------------------------------------------------
#endif

!-----------------------------------------------------------------------
! Set up a namelist to allow switches to be set
! ONLY USED IN UM (for now!!)
!-----------------------------------------------------------------------
  NAMELIST /jules_switches/ l_point_data,l_snow_albedo,l_phenol,  &
                            l_triffid,l_vg_soil,l_aggregate,      &
                            l_epot_corr,                          &
                            can_model,soilhc_method,              &
                            i_modiscopt,                          & 
                            frac_snow_subl_melt,all_tiles,        &
                            can_rad_mod, cor_mo_iter,             &
                            iseaz0t, buddy_sea,                   &
                            iscrntdiag, b_pdm, dz_pdm
                            

END MODULE switches
