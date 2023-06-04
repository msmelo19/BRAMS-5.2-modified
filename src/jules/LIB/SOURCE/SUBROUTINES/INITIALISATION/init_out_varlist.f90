!###############################################################################
!###############################################################################
! subroutine init_out_varlist
! Creates a list of available variables.
! Every variable listed here should also have code in whatever routine(s) actually
! do the output (e.g. subroutine loadout).
! This routine could be replaced by a simple list - but this code makes it easier to
! add new variables (can simply be entered anywhere in list) and ensures that the
! variables are dimensioned correctly by counting beforehand.
! At present, available variable types are: LA, PF, RG, RP, SC, SI, SO, TI, TY

  SUBROUTINE init_out_varlist()

  USE ancil_info, ONLY : &
!  imported scalars with intent(in)
     dim_cs1

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE inout, ONLY :  &
!   imported scalars with intent(in)
     echo  &
!   imported arrays with intent(out)
    ,nvar,varDescList,varNameList,varTypeList,varUnitsList

  USE soil_param, ONLY : zsmc

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  SCALARS
    i              &!  counter
   ,iloop          &!  counter
   ,j               !  work
!-------------------------------------------------------------------------------

! On first loop simply count the number of variables.
! On second loop, allocate space and store values.
  i = 0  !  initialise counter
  DO iloop=1,2

    IF ( iloop==2 ) THEN
      nvar = i
      CALL allocate_arrays( 'init_out_varlist' )
      i = 0  !  reset counter for this (2nd) loop
    ENDIF

!###############################################################################
!###############################################################################
!   Variables that have a single value at LAND gridpoints only.
!   varTypeList = 'LA'
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'albedoLand'
      varDescList(i) = 'Gridbox albedo (as used for net shortwave calcualtion)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'canopy'
      varDescList(i) = 'Gridbox canopy water content'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'cs'
      varDescList(i) = 'Gridbox soil carbon (total)'
      varUnitsList(i) = 'kgC m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'cv'
      varDescList(i) = 'Gridbox mean vegetation carbon'
      varUnitsList(i) = 'kgC m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'depthFrozen'
      varDescList(i) = 'Gridbox depth of frozen ground at surface'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'depthUnfrozen'
      varDescList(i) = 'Gridbox depth of unfrozen ground at surface'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'drain'
      varDescList(i) = 'Drainage from bottom (nshyd) soil layer'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'elake'
      varDescList(i) = 'Gridbox mean evaporation from lakes'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'emis'
      varDescList(i) = 'Gridbox emissivity'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fch4_wetl'
      varDescList(i) = 'Scaled methane flux from wetland fraction'
      varUnitsList(i) = '10^-9 kgC m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fsat'
      varDescList(i) = 'Surface saturated fraction'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fsmc'
      varDescList(i) = 'Gridbox soil moisture availability factor (beta)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fwetl'
      varDescList(i) = 'Wetland fraction'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gpp'
      varDescList(i) = 'Gridbox gross primary productivity'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gs'
      varDescList(i) = 'Gridbox surface conductance to evaporation'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'hfSnowMelt'
      varDescList(i) = 'Gridbox snowmelt heat flux'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'landIndex'
      varDescList(i) = 'Index (gridbox number) of land points'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'liceIndex'
      varDescList(i) = 'Index (gridbox number) of land ice points'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'litCMn'
      varDescList(i) = 'Gridbox mean carbon litter'
      varUnitsList(i) = 'kgC  m-2 per 360days'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'LWnet'
      varDescList(i) = 'Gridbox net downward longwave radiation at surface'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'LWup'
      varDescList(i) = 'Gridbox surface upward LW radiation of land points'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'npp'
      varDescList(i) = 'Gridbox net primary productivity'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'qbase'
      varDescList(i) = 'Baseflow (lateral subsurface runoff)'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'qbase_zw'
      varDescList(i) = 'Baseflow from deep layer'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'radnet'
      varDescList(i) = 'Surface net radiation of land points'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respP'
      varDescList(i) = 'Gridbox plant respiration'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respS'
      varDescList(i) = 'Gridbox soil respiration (total)'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respSDrOut'
      varDescList(i) = 'Gridbox mean soil respiration for driving TRIFFID'
      varUnitsList(i) = 'kgC  m-2 per 360days'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'runoff'
      varDescList(i) = 'Gridbox runoff rate'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'sat_excess_roff'
      varDescList(i) = 'Saturation excess surface ("Dunne") runoff'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'smcAvailTop'
      varDescList(i) = 'Gridbox available moisture in top'
!     Include a depth in the variable description. Assume varDescList is sufficiently long.
      j = LEN_TRIM( varDescList(i) ) + 1
      WRITE(varDescList(i)(j+1:j+3),"(f3.1)") zsmc
      varDescList(i) = TRIM(varDescList(i)) // 'm  of soil'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'smcAvailTot'
      varDescList(i) = 'Gridbox available moisture in soil column'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'smcTot'
      varDescList(i) = 'Gridbox total soil moisture in column'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snomltSubHtf'
      varDescList(i) = 'Gridbox sub-canopy snowmelt heat flux'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowCan'
      varDescList(i) = 'Gridbox snow on canopy'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowDepth'
      varDescList(i) = 'Snow depth (on ground)'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowFrac'
      varDescList(i) = 'Gridbox snow-covered fraction of land points'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowFracAlb'
      varDescList(i) = 'Gridbox average weight given to snow in albedo'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowGrCan'
      varDescList(i) = 'Gridbox snow below canopy (snow_grnd)'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowIceTot'
      varDescList(i) = 'Gridbox ice content of snow on ground'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowLiqTot'
      varDescList(i) = 'Gridbox liquid content of snow on ground'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowMelt'
      varDescList(i) = 'Gridbox rate of snowmelt'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'soilIndex'
      varDescList(i) = 'Index (gridbox number) of soil points'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'sthzw'
      varDescList(i) = 'Soil wetness in deep (water table/TOPMODEL) layer'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'subSurfRoff'
      varDescList(i) = 'Gridbox sub-surface runoff'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'surfRoff'
      varDescList(i) = 'Gridbox surface runoff'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'surfRoffInf'
      varDescList(i) = 'Gridbox infiltration excess surface runoff'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'swetLiqTot'
      varDescList(i) = 'Gridbox unfrozen soil moisture as fraction of saturation'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'swetTot'
      varDescList(i) = 'Gridbox soil moisture as fraction of saturation'
      varUnitsList(i) = '-'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'SWnet'
      varDescList(i) = 'Gridbox net downward shortwave radiation at surface'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tfall'
      varDescList(i) = 'Gridbox throughfall'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'trad'
      varDescList(i) = 'Gridbox effective radiative temperature (assuming emissivity=1)'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'wFluxSfc'
      varDescList(i) = 'Gridbox downwards moisture flux at soil surface'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'LA'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'zw'
      varDescList(i) = 'Gridbox mean depth to water table'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'LA'
    ENDIF
!###############################################################################
!###############################################################################
! Vegetation (PFT) diagnostics (on land grid).
! varTypeList = 'PF'
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'cVegP'
      varDescList(i) = 'PFT total carbon content of the vegetation'
      varUnitsList(i) = 'kgC m-2'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'canhtP'
      varDescList(i) = 'PFT canopy height'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ciP'
      varDescList(i) = 'PFT internal CO2 pressure '
      varUnitsList(i) = 'Pa'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fluxO3Stom'
      varDescList(i) = 'Flux of O3 to stomata'
      varUnitsList(i) = 'mol m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fsmcP'
      varDescList(i) = 'PFT soil moisture availability factor (beta)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gLeafP'
      varDescList(i) = 'PFT leaf turnover rate'
      varUnitsList(i) = 'per 360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gLeafDayP'
      varDescList(i) = 'PFT mean leaf turnover rate for input to PHENOL'
      varUnitsList(i) = 'per 360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gLeafDrOutP'
      varDescList(i) = 'PFT mean leaf turnover rate for driving TRIFFID'
      varUnitsList(i) = 'per 360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gLeafPhenP'
      varDescList(i) = 'PFT mean leaf turnover rate over phenology period'
      varUnitsList(i) = 'per360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gstomP'
      varDescList(i) = 'PFT bulk (canopy) stomatal conductance for water vapour'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gppP'
      varDescList(i) = 'PFT gross primary productivity'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'laiP'
      varDescList(i) = 'PFT leaf area index'
      varUnitsList(i) = '-'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'laiPhenP'
      varDescList(i) = 'PFT leaf area index after phenology.'
      varUnitsList(i) = '-'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'litCP'
      varDescList(i) = 'PFT carbon Litter'
      varUnitsList(i) = 'kgC  m-2 per 360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'nppDrOutP'
      varDescList(i) = 'PFT mean NPP for driving TRIFFID'
      varUnitsList(i) = 'kgC  m-2 per 360days'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'nppP'
      varDescList(i) = 'PFT net primary productivity'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'o3ExpFac'
      varDescList(i) = 'Ozone exposure factor'
      varUnitsList(i) = '-'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rdcP'
      varDescList(i) = 'Canopy dark respiration, without soil water dependence'
      varUnitsList(i) = 'molCO2 m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respPP'
      varDescList(i) = 'PFT plant respiration'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respWDrOutP'
      varDescList(i) = 'PFT mean wood respiration for driving TRIFFID (kg C/m2/360days)  '
      varUnitsList(i) = '-'
      varTypeList(i) = 'PF'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respWP'
      varDescList(i) = 'PFT wood maintenance respiration'
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'PF'
    ENDIF
!###############################################################################
!###############################################################################
!###############################################################################
! Variables on the channel routing grid.
! varTypeList = 'RG'
!-------------------------------------------------------------------------------
!    i = i + 1
!    if ( iloop == 2 ) then
!      varNameList(i) = 'bflow'
!      varDescList(i) = 'baseflow in river gridcells'
!      varUnitsList(i) = 'm3 s-1'
!      varTypeList(i) = 'RG'
!    endif
!-------------------------------------------------------------------------------
!    i = i + 1
!    if ( iloop == 2 ) then
!      varNameList(i) = 'bflowin'
!      varDescList(i) = 'subsurface inflow of G2G model'
!      varUnitsList(i) = 'mm'
!      varTypeList(i) = 'RG'
!    endif
!-------------------------------------------------------------------------------
!    i = i + 1
!    if ( iloop == 2 ) then
!      varNameList(i) = 'flowin'
!      varDescList(i) = 'surface inflow of G2G model'
!      varUnitsList(i) = 'mm'
!      varTypeList(i) = 'RG'
!    endif
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'flowdir'
      varDescList(i) = 'coded flow direction on routing grid (1-8, clockwise from N)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rflow'
      varDescList(i) = 'river flow rate'
      varUnitsList(i) = 'kg s-1'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'routeIndex'
      varDescList(i) = 'routeIndex (location on routing grid)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rInflow'
      varDescList(i) = 'river inflow rate (in channel)'
      varUnitsList(i) = 'kg s-1'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'routeNext'
      varDescList(i) = 'routeNext (location on routing grid of downstream point)'
      varUnitsList(i) = '-'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'runoffR'
      varDescList(i) = 'runoff on routing grid'
      varUnitsList(i) = 'kg s-1'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rstore'
      varDescList(i) = 'river channel storage'
      varUnitsList(i) = 'kg'
      varTypeList(i) = 'RG'
    ENDIF
!-------------------------------------------------------------------------------
!    i = i + 1
!    if ( iloop == 2 ) then
!      varNameList(i) = 'substore'
!      varDescList(i) = 'subsurface store of G2G model'
!      varUnitsList(i) = 'mm'
!      varTypeList(i) = 'RG'
!    endif
!-------------------------------------------------------------------------------
!    i = i + 1
!    if ( iloop == 2 ) then
!      varNameList(i) = 'surfstore'
!      varDescList(i) = 'surface store of G2G model'
!      varUnitsList(i) = 'mm'
!      varTypeList(i) = 'RG'
!    endif
!###############################################################################
!###############################################################################
! Variables at selected points on the channel routing grid.
! varTypeList = 'RP'
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rflowP'
      varDescList(i) = 'riverflow rate at point'
      varUnitsList(i) = 'kg s-1'
      varTypeList(i) = 'RP'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rstoreP'
      varDescList(i) = 'river channel storage at a point'
      varUnitsList(i) = 'kg'
      varTypeList(i) = 'RP'
    ENDIF
!###############################################################################
!   Variables that have dim_cs1 values at every land gridpoint.
!   These are soil carbon and related respiration variables.
!   dim_cs1=4 if TRIFFID (and RothC) are used, else =1.
!   varTypeList = 'SC'
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'csPool'
      IF ( dim_cs1 == 4 ) THEN
        varDescList(i) = 'Gridbox soil carbon in each pool (DPM,RPM,bio,hum)'
      ELSEIF ( dim_cs1 == 1 ) THEN
        varDescList(i) = 'Gridbox soil carbon (single pool)'
      ELSE
        WRITE(*,*)'ERROR: init_out_varlist: only have code for dim_cs1=1 or 4.'
        WRITE(*,*)'This is not too important, and you could proceed past this'
        WRITE(*,*)'stage using existing code. But with the present code, you''ll'
        WRITE(*,*)'fall over later in the code!'
        STOP
      ENDIF
      varUnitsList(i) = 'kgC m-2'
      varTypeList(i) = 'SC'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'respSPool'
      IF ( dim_cs1 == 4 ) THEN
        varDescList(i) = 'Gridbox soil respiration from each pool (DPM,RPM,bio,hum)'
      ELSE
        varDescList(i) = 'Gridbox soil respiration from each pool'
      ENDIF
      varUnitsList(i) = 'kgC m-2 s-1'
      varTypeList(i) = 'SC'
    ENDIF
!###############################################################################
!   Variables that have a single value at every gridpoint over land or sea.
!   varTypeList = 'SI'
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'conRain'
      varDescList(i) = 'Gridbox convective rainfall'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'conSnow'
      varDescList(i) = 'Gridbox convective snowfall'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'cosz'
      varDescList(i) = 'Cosine of the zenith angle'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'diffFrac'
      varDescList(i) = 'Gridbox fraction of radiation that is diffuse'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
!   Make the "extra" driving variables available.
!    do ivar=1,ndriveExtra
!      i = i + 1
!      IF ( iloop == 2 ) THEN
!        WRITE(varnameList(i),"(a10,i2.2)") 'extraDrive',ivar
!        WRITE(varDescList(i),"(a24,i2.2)") 'Extra driving variable #',ivar
!        varUnitsList(i) = 'unknown'
!      varTypeList(i) = 'SI'
!      ENDIF
!    enddo
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ecan'
      varDescList(i) = 'Gridbox mean evaporation from canopy/surface store'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ei'
      varDescList(i) = 'Gridbox sublimation from lying snow or sea-ice'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'esoil'
      varDescList(i) = 'Gridbox surface evapotranspiration from soil moisture store'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fqw'
      varDescList(i) = 'Gridbox moisture flux from surface'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ftl'
      varDescList(i) = 'Gridbox surface sensible heat flux'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'landAlbedo1'
      varDescList(i) = 'Gridbox albedo for waveband 1'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'landAlbedo2'
      varDescList(i) = 'Gridbox albedo for waveband 2.'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'landAlbedo3'
      varDescList(i) = 'Gridbox albedo for waveband 3.'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'landAlbedo4'
      varDescList(i) = 'Gridbox albedo for waveband 4.'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'latentHeat'
      varDescList(i) = 'Gridbox surface latent heat flux'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'latitude'
      varDescList(i) = 'Gridbox latitude'
      varUnitsList(i) = 'degrees N'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'longitude'
      varDescList(i) = 'Gridbox longitude'
      varUnitsList(i) = 'degrees E'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'lsRain'
      varDescList(i) = 'Gridbox large-scale rainfall'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'lsSnow'
      varDescList(i) = 'Gridbox large-scale snowfall'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'LWdown'
      varDescList(i) = 'Gridbox surface downward LW radiation'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'precip'
      varDescList(i) = 'Gridbox precipitation rate'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'pstar'
      varDescList(i) = 'Gridbox surface pressure'
      varUnitsList(i) = 'Pa'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'q1p5m'
      varDescList(i) = 'Gridbox specific humidity at 1.5m height'
      varUnitsList(i) = 'kg kg-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'qw1'
      varDescList(i) = 'Gridbox specific humidity (total water content)'
      varUnitsList(i) = 'kg kg-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rainfall'
      varDescList(i) = 'Gridbox rainfall rate'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snomltSurfHtf'
      varDescList(i) = 'Gridbox heat flux used for surface melting of snow'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowfall'
      varDescList(i) = 'Gridbox snowfall rate'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowMass'
      varDescList(i) = 'Gridbox snowmass (on land this is lying_snow=snow_tile+snow_grnd)'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'surfHtFlux'
      varDescList(i) = 'Gridbox net downward heat flux at surface over land and sea-ice fraction of gridbox'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'SWdown'
      varDescList(i) = 'Gridbox surface downward SW radiation'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 't1p5m'
      varDescList(i) = 'Gridbox temperature at 1.5m height'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'taux1'
      varDescList(i) = 'Gridbox westerly component of surface wind stress'
      varUnitsList(i) = 'N m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tauy1'
      varDescList(i) = 'Gridbox southerly component of surface wind stress'
      varUnitsList(i) = 'N m-2'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tl1'
      varDescList(i) = 'Gridbox ice/liquid water temperature'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tstar'
      varDescList(i) = 'Gridbox surface temperature'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'u1'
      varDescList(i) = 'Gridbox westerly wind component'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'u10m'
      varDescList(i) = 'Gridbox westerly wind component at 10 m height'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'v1'
      varDescList(i) = 'Gridbox southerly wind component'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'v10m'
      varDescList(i) = 'Gridbox southerly wind component at 10m height'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'SI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'wind'
      varDescList(i) = 'Gridbox wind speed'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'SI'
    ENDIF
!###############################################################################
!###############################################################################
! Snow layer diagnostics (on tiles).
! varTypeList = 'SN'
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rgrainL'
      varDescList(i) = 'Grain size in snow layers'
      varUnitsList(i) = 'microns'
      varTypeList(i) = 'SN'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowDs'
      varDescList(i) = 'Depth of snow layers'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'SN'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowIce'
      varDescList(i) = 'Ice mass in snow layers'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'SN'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowLiq'
      varDescList(i) = 'Liquid mass in snow layers'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'SN'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tsnow'
      varDescList(i) = 'Temperature of snow layers'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'SN'
    ENDIF
!###############################################################################
!###############################################################################
! Soil layer diagnostics.
! varTypeList = 'SO'
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'bSoil'
      varDescList(i) = 'Gridbox Brooks-Corey exponent for each soil layer'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ext'
      varDescList(i) = 'Gridbox extraction of water from each soil layer'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'hCapSoil'
      varDescList(i) = 'Gridbox soil heat capacity for each soil layer'
      varUnitsList(i) = 'J K-1 m-3'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'hConSoil'
      varDescList(i) = 'Gridbox soil thermal conductivity for each soil layer'
      varUnitsList(i) = 'W m-1 K-1'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'satCon'
      varDescList(i) = 'Gridbox saturated hydraulic conductivity for each soil layer'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'sathh'
      varDescList(i) = 'Gridbox saturated soil water pressure for each soil layer'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'smcl'
      varDescList(i) = 'Gridbox moisture content of each soil layer'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'soilWet'
      varDescList(i) = 'Gridbox total moisture content of each soil layer, as fraction of saturation'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'sthf'
      varDescList(i) = 'Gridbox frozen moisture content of each soil layer as a fraction of saturation'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'sthu'
      varDescList(i) = 'Gridbox unfrozen moisture content of each soil layer as a fraction of saturation'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tSoil'
      varDescList(i) = 'Gridbox sub-surface temperature of each layer'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'vsmcCrit'
      varDescList(i) = 'Gridbox volumetric moisture content at critical point for each soil layer'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'vsmcSat'
      varDescList(i) = 'Gridbox volumetric moisture content at saturation for each soil layer'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'vsmcWilt'
      varDescList(i) = 'Gridbox  volumetric moisture content at wilting point for each soil layer'
      varUnitsList(i) = '-'
      varTypeList(i) = 'SO'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'wFlux'
      varDescList(i) = 'Downwards moisture flux at bottom of each soil layer'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'SO'
    ENDIF
!###############################################################################
!###############################################################################
! Tile diagnostics (on land grid).
! varTypeList = 'TI'
!-------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'alb1T'
      varDescList(i) = 'Tile land albedo, waveband 1'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'alb2T'
      varDescList(i) = 'Tile land albedo, waveband 2'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'alb3T'
      varDescList(i) = 'Tile land albedo, waveband 3'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'alb4T'
      varDescList(i) = 'Tile land albedo, waveband 4'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'anthropHtFluxT'
      varDescList(i) = 'Anthropogenic heat flux for each tile'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'canopyT'
      varDescList(i) = 'Tile surface/canopy water for snow-free land tiles'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'catchT'
      varDescList(i) = 'Tile surface/canopy water capacity of snow-free land tiles'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ecanT'
      varDescList(i) = 'Tile evaporation from canopy/surface store for snow-free land tiles'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'eiT'
      varDescList(i) = 'Tile sublimation from lying snow for land tiles'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'emisT'
      varDescList(i) = 'Tile emissivity'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'esoilT'
      varDescList(i) = 'Tile surface evapotranspiration from soil moisture store for snow-free land tiles'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'fqwT'
      varDescList(i) = 'Tile surface moisture flux for land tiles'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'ftlT'
      varDescList(i) = 'Tile surface sensible heat flux for land tiles'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'gcT'
      varDescList(i) = 'Tile surface conductance to evaporation for land tiles'
      varUnitsList(i) = 'm s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'leT'
      varDescList(i) = 'Tile surface latent heat flux for land tiles'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'nsnow'
      varDescList(i) = 'Number of snow layers'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'q1p5mT'
      varDescList(i) = 'Tile specific humidity at 1.5m over land tiles'
      varUnitsList(i) = 'kg kg-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'radnetT'
      varDescList(i) = 'Tile surface net radiation'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'rgrainT'
      varDescList(i) = 'Tile snow surface grain size'
      varUnitsList(i) = 'microns'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowCanMeltT'
      varDescList(i) = 'Tile melt of canopy snow'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowCanT'
      varDescList(i) = 'Tile canopy snow'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowDepthT'
      varDescList(i) = 'Tile snow depth (on ground)'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowGrCanMeltT'
      varDescList(i) = 'Tile melt of snow under canopy'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowGroundRhoT'
      varDescList(i) = 'Tile density of snow on ground'
      varUnitsList(i) = 'kg m-3'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowGrCanT'
      varDescList(i) = 'Tile snow on ground below canopy (snow_grnd)'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowGroundT'
      varDescList(i) = 'Tile snow on ground (snow_tile or snow_grnd)'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowIceT'
      varDescList(i) = 'Tile total ice content in snow on ground'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowLiqT'
      varDescList(i) = 'Tile total liquid content in snow on ground'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowMassT'
      varDescList(i) = 'Tile lying snow (snow_tile+snow_grnd)'
      varUnitsList(i) = 'kg m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'snowMeltT'
      varDescList(i) = 'Tile snow melt (melt_tile)'
      varUnitsList(i) = 'kg m-2 s-1'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'surfHtFluxT'
      varDescList(i) = 'Downward heat flux for each tile'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'surfHtStoreT'
      varDescList(i) = 'C*(dT/dt) for each tile'
      varUnitsList(i) = 'W m-2'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 't1p5mT'
      varDescList(i) = 'Tile temperature at 1.5m over land tiles'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tstarT'
      varDescList(i) = 'Tile surface temperature'
      varUnitsList(i) = 'K'
      varTypeList(i) = 'TI'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'z0T'
      varDescList(i) = 'Tile surface roughness'
      varUnitsList(i) = 'm'
      varTypeList(i) = 'TI'
    ENDIF
!###############################################################################
!###############################################################################
! Surface type diagnostics (on land grid).
! varTypeList = 'TY'
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'frac'
      varDescList(i) = 'Fractional cover of each surface type.  '
      varUnitsList(i) = '-'
      varTypeList(i) = 'TY'
    ENDIF
!--------------------------------------------------------------------------------
    i = i + 1
    IF ( iloop == 2 ) THEN
      varNameList(i) = 'tileIndex'
      varDescList(i) = 'Index (gridbox number) of land points with each surface type'
      varUnitsList(i) = '-'
      varTypeList(i) = 'TY'
    ENDIF
!###############################################################################
!###############################################################################
!###############################################################################
  ENDDO  !  iloop
!###############################################################################
!###############################################################################

! Some checks - guarding against coding errors.
  DO i=1,nvar
!   Check that all fields are present.
    IF ( varNameList(i)=='' .OR. varDescList(i)=='' .OR. varTypeList(i)=='' ) THEN
      WRITE(*,*)'ERROR: incorrect specification of a variable.'
      WRITE(*,*) i,' name,desc,type=',TRIM(varNameList(i)),TRIM(varDescList(i)),TRIM(varTypeList(i))
      WRITE(*,*)'Stopping in init_out_varlist'
      STOP
    ENDIF
!   Check that type of variable is recognised.
!   Note this has to be updated if new variable types are added!
    SELECT CASE ( varTypeList(i) )
      CASE( 'LA','PF','RG','RP','SC','SI','SN','SO','TI','TY' )
!       OK, nothing to be done.
      CASE default
        WRITE(*,*) 'ERROR: do not recognise variable type ',varTypeList(i)
        WRITE(*,*)'Stopping in init_out_varlist'
        STOP
    END SELECT
!    if ( .NOT. any(varTypeListList(:)==varTypeList(i)) ) then
!      write(*,*) 'ERROR: do not recognise variable type ',varTypeList(i)
!      write(*,*)'Stopping in init_out_varlist'
!      stop
!    endif

!   Check that variable name is not repeated.
    DO j=i+1,nvar
      IF ( varNameList(i) == varNameList(j) ) THEN
        WRITE(*,*) 'ERROR: repeated variable: ',TRIM(varNameList(i))
        WRITE(*,*) 'NB This means that the code uses the same name for two variables!'
        WRITE(*,*)'Stopping in init_out_varlist'
        STOP
      ENDIF
    ENDDO
  ENDDO
!--------------------------------------------------------------------------------

! Output list to screen.
  IF ( echo ) THEN
    WRITE(*,"(/,a)") '--------Available variables (see init_out_varlist):----------------'
    WRITE(*,*)'varNameList,varTypeList,varDescList'
    DO i=1,nvar
      WRITE(*,"(i3,3(tr1,a))") i,varNameList(i),varTypeList(i),TRIM(varDescList(i))
    ENDDO
    WRITE(*,"(70('-'))")
  ENDIF

  END SUBROUTINE init_out_varlist
!###############################################################################
!###############################################################################
