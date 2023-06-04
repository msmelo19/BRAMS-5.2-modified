!#######################################################################
!#######################################################################
! subroutine init_opts
! Subroutine to read model options and misc other values.

!-----------------------------------------------------------------------
  SUBROUTINE init_opts

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY :  &
!   imported scalars with intent(out)
     ntiles,sm_levels,nsmax

  USE file_utils, ONLY:  &
!   imported procedures
     findTag

  USE inout, ONLY : &
!   imported scalars with intent(out)
     echo,nxIn,nyIn,print_step,jinUnit,yrevIn

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(out)
    ncType

  USE misc_utils, ONLY :  &
!  imported procedures
     repeatVal

  USE nstypes, ONLY :  &
!   imported scalars with intent(out)
     nnvg,npft,ntype

  USE nvegparm, ONLY :  &
!   imported arrays with intent(out)
     nvgName

  USE pftparm, ONLY :  &
!   imported arrays with intent(out)
     pftName

  USE route_mod, ONLY : &
!   imported scalars with intent(out)
     nxRoute,nyRoute

  USE switches, ONLY :  &
!   imported scalars with intent(out)
     can_model,can_rad_mod,ilayers,l_co2_interactive,l_cosz  &
    ,l_dust,l_pdm,l_phenol,l_snow_albedo,l_spec_albedo,l_top &
    ,l_triffid,l_trif_eq,l_z0_orog,route,routeOnly,l_aggregate &
    ,l_anthrop_heat_src,ISCRNTDIAG,l_o3_damage,l_veg_compete &
    ,l_imogen,l_epot_corr,l_snowdep_surf

  USE switches_urban, ONLY :  &
!   imported scalars with intent(out)
     l_moruses

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    ierr    !  error value

  CHARACTER(len=LEN(pftName)), ALLOCATABLE ::  &!  local arrays
    typeName(:)    !   names of surface types
!------------------------------------------------------------------------------

! Prevent the use of certain options.
  L_CO2_INTERACTIVE=.FALSE.
  L_DUST=.FALSE.
  L_Z0_OROG=.FALSE.

!-------------------------------------------------------------------------------
! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init','>INIT_OPTS' )

  READ(jinUnit,*) npft,nnvg
  READ(jinUnit,*) l_aggregate

! Allocate space.
  CALL allocate_arrays( 'init_opts' )

  READ(jinUnit,*) pftName(:)
  READ(jinUnit,*) nvgName(:)

  READ(jinUnit,*) nxIn,nyIn
  READ(jinUnit,*) sm_levels
  READ(jinUnit,*) NSMAX
  READ(jinUnit,*) CAN_MODEL
  READ(jinUnit,*) CAN_RAD_MOD,ILAYERS
  READ(jinUnit,*) L_COSZ,L_SPEC_ALBEDO
  READ(jinUnit,*) L_PHENOL,L_TRIFFID,l_veg_compete,L_TRIF_EQ
  READ(jinUnit,*) l_top,l_pdm
  READ(jinUnit,*) route,routeOnly
  READ(jinUnit,*) l_anthrop_heat_src,l_moruses
  READ(jinUnit,*) l_o3_damage
  READ(jinUnit,*) l_imogen
  READ(jinUnit,*) l_epot_corr,l_snowdep_surf
  
  READ(jinUnit,*) ISCRNTDIAG
  READ(jinUnit,*) yrevIn
  READ(jinUnit,*) ncType
  READ(jinUnit,*) echo
  READ(jinUnit,*) PRINT_STEP

! Have finished reading this section.
!-------------------------------------------------------------------------------

! If no routing, make sure routeOnly is FALSE.
  IF ( .NOT. route ) routeOnly=.FALSE.

!-------------------------------------------------------------------------------
! Set and check values.
!-------------------------------------------------------------------------------
  IF ( .NOT. routeOnly ) THEN

!   Calculate number of surface types.
    ntype = npft + nnvg
    IF(l_aggregate) THEN
      ntiles=1
    ELSE
      ntiles=ntype
    ENDIF

    IF ( l_veg_compete .AND. ( .NOT. l_triffid ) ) THEN
      WRITE(*,*)'Cannot use competing vegetation with l_triffid = false'
      WRITE(*,*)'Stopping in INIT_OPTS'
      STOP
    END IF

    IF ( l_veg_compete .AND. npft /= 5 ) THEN
      WRITE(*,*)'Competition in TRIFFID is hardwired to expect certain PFTs.'
      WRITE(*,*)'npft must be 5 (and they must be the "expected" PFTs).'
      WRITE(*,*)'Stopping in INIT_OPTS'
      STOP
    ENDIF

    IF ( l_aggregate .AND. (l_phenol .OR. l_triffid .OR. l_trif_eq) ) THEN
      WRITE(*,*)'Phenology or TRIFFID can not used with'
      WRITE(*,*)'the aggregated surface scheme (i.e. l_aggregate = true).'
      WRITE(*,*)'Stopping in INIT_OPTS'
      STOP
    ENDIF

    IF ( can_model<1 .OR. can_model>4 ) THEN
      WRITE(*,*)'ERROR: can_model should be in range 1 to 4.'
      WRITE(*,*)'Stopping in INIT_OPTS'
      STOP
    ENDIF

    IF ( can_model==4 .AND. l_aggregate ) THEN
      WRITE(*,*)'ERROR: init_opts: can_model=4 cannot be used with'
      WRITE(*,*)'the aggregated surface scheme (i.e. l_aggregate = true).'
      STOP
    ENDIF

    IF ( can_rad_mod<1 .OR. can_rad_mod>5 ) THEN
      WRITE(*,*)'ERROR: can_rad_mod should be in range 1 to 5.'
      WRITE(*,*)'Stopping in INIT_OPTS'
      STOP
    ENDIF

!   Check that names of surface types are unique.
    ALLOCATE( typeName(npft+nnvg), stat=ierr )
    IF ( ierr /= 0 ) THEN
      WRITE(*,*)'ERROR: init_opts: could not allocate.'
      STOP
    ENDIF
    typeName(1:npft) = pftname(:)
    typeName(npft+1:) = nvgName(:)
    IF ( repeatVal( typeName(:) ) ) THEN
      WRITE(*,*)'ERROR: init_opts: repeated name of surface type.'
      WRITE(*,*)'All chosen surface types must have unique names.'
      STOP
    ENDIF
    DEALLOCATE( typeName, stat=ierr )
    IF ( ierr /= 0 ) WRITE(*,*)'WARNING: init_opts: error on deallocate.'

!   If spectral albedo selected, use spectral snow albedo.
    l_snow_albedo = .FALSE.
    IF ( l_spec_albedo ) l_snow_albedo = .TRUE.

!   Check that at most one of TOPMODEL and PDM have been selected.
    IF ( l_pdm .AND. l_top ) THEN
      WRITE(*,*)'ERROR: init_opts: l_pdm and l_top'
      WRITE(*,*)'At most one of PDM and TOPMODEL can be selected.'
      STOP
    ENDIF

!   Reset some sizes that were read from file, so that they do not cause work
!   space to be unnecessarily large.
    nxRoute = 1
    nyRoute = 1

  ELSE   !  routeOnly

!   Ensure switches are off (to prevent trying to initialise etc).
    l_pdm = .FALSE.
    l_top = .FALSE.

!   Reset some sizes that were read from file, so that they do not cause work
!   space to be unnecessarily large.
    npft = 1
    nnvg = 1
    ntype = 1
    ntiles = 1
    sm_levels = 1
    ilayers = 1
    nsmax = 0

  ENDIF  !  routeOnly

!-------------------------------------------------------------------------------
! Messages to screen.
!-------------------------------------------------------------------------------
  IF ( yrevIn ) WRITE(*,*)'WARNING: yrevIn=T: input data are in "north to south" order'

  IF ( routeOnly ) THEN
    WRITE(*,*)'########## routeOnly ############'
    WRITE(*,*)'This means that only runoff routing is performed, all other'
    WRITE(*,*)'sections of model are inactive.'
  ENDIF

  WRITE(*,*)'Any netCDF files to be read will be of type: ',TRIM(ncType)

  END SUBROUTINE init_opts
!###############################################################################
!###############################################################################
