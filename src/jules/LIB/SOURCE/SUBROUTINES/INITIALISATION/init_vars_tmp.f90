  SUBROUTINE init_vars_tmp(sst,nx,ny)

  USE Ancil_info
  USE Switches
  USE Top_pdm, ONLY : inlandout_atm
  USE Forcing, ONLY : u_0,v_0
  USE Prognostics
  USE Aero
  USE Orog
  USE Trifctl
  USE Coastal
  USE C_Densty
  USE Sea_ice
  USE C_Rough
  USE inout, ONLY : echo
  USE p_s_parms, ONLY : satcon,soil_clay
  USE c_elevate, ONLY : z_land

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER I,J,L,N     ! Loop counters

  
  !{DSM
   REAL    :: sst(ROW_LENGTH,ROWS)
   INTEGER :: nx,ny
   
   IF (nx /= ROW_LENGTH .or. ny /= ROWS) THEN
      PRINT*, 'nx/=ROW_LENGTH ou ny/=ROWS - nx,ROW_LENGTH,ny,ROWS =',nx,ROW_LENGTH,ny,ROWS
      STOP
   ENDIF
   
   !print*,sst
   !pause
  !DSM}
    
!-----------------------------------------------------------------------
! If only doing routing, nothing need be set in here.
!-----------------------------------------------------------------------
  IF ( routeOnly ) RETURN

!-----------------------------------------------------------------------
! Set number of timesteps since TRIFFID to be zero
!-----------------------------------------------------------------------
  ASTEPS_SINCE_TRIFFID=0

!-----------------------------------------------------------------------
! Initialise accumulated fluxes for TRIFFID and phenology.
! This is not necessary if these are read from a restart file - but at
! present they're not.
!-----------------------------------------------------------------------
  g_leaf_acc(:,:) = 0.0
  g_leaf_phen_acc(:,:) = 0.0
  npp_ft_acc(:,:) = 0.0
  resp_s_acc(:,:) = 0.0
  resp_w_ft_acc(:,:) = 0.0

!-----------------------------------------------------------------------
! Set saturated hydraulic conductivity to zero at land ice points.
!-----------------------------------------------------------------------
  IF ( lice_pts > 0 ) THEN
    IF ( echo ) WRITE(*,*) 'Setting satcon to zero at land ice points.'
    DO i=1,lice_pts
      l = lice_index(i)
      satcon(l,:) = 0.0
    ENDDO
  ENDIF

!-----------------------------------------------------------------------
! Set surface velocity to be zero
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      U_0(I,J)=0.0
      V_0(I,J)=0.0
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set CO2 variables
!-----------------------------------------------------------------------
  DO I=1,co2_dim_len
    DO J=1,co2_dim_row
      CO2_3D(I,J)=0.0
    ENDDO
  ENDDO


!-----------------------------------------------------------------------
! Set coastal tiling variables
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      !DSM TSTAR_SEA(I,J)=280.0
      TSTAR_SEA(I,J)=sst(I,J)  !DSM
      TSTAR_SICE(I,J)=270.0
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set orographic roughness variables
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      H_BLEND_OROG(I,J)=0.0
    ENDDO
  ENDDO
  DO L=1,LAND_PTS
    SIL_OROG_LAND(L)=0.0
    HO2R2_OROG(L)=0.0
  ENDDO

!-----------------------------------------------------------------------
! Set up prognostics which are not currently in dump
!-----------------------------------------------------------------------
  DO I=1,ROW_LENGTH
    DO J=1,ROWS
      TI(I,J)=270.0
      Z0MSEA(I,J)=Z0HSEA
      SNOW_MASS(I,J)=0.0
      SNOW_MASS_SEA(I,J)=0.0
      DO L=1,NICE
        ICE_FRACT_NCAT(I,J,L)=0.0
        DI_NCAT(I,J,L)=0.1
      ENDDO
    ENDDO
  ENDDO

!-----------------------------------------------------------------------
! Set logical switch for snow affected sea-ice albedo to false
!-----------------------------------------------------------------------
  L_SSICE_ALBEDO= .FALSE.

!-----------------------------------------------------------------------
! Set up sea-ice parameter variables
!-----------------------------------------------------------------------
  IF(L_SSICE_ALBEDO) THEN
    SW_alpham=0.65
    SW_alphac=0.80
    SW_alphab=0.57
    SW_dtice=2.00
  ELSE
    SW_alpham=0.50
    SW_alphac=0.80
    SW_alphab=-1.00
    SW_dtice=10.00
  ENDIF

  DO I=1,ROW_LENGTH*ROWS
    SOOT(I)=0.0
  ENDDO

!------------------------------------------------------------------
! Temporary initialisation of variables added during UM integration
!------------------------------------------------------------------
  z_land(:,:)     = 0.0
  soil_clay(:,:)  = 0.0
  
  INLANDOUT_ATM(:) = 0.0

  RETURN
  END SUBROUTINE init_vars_tmp


