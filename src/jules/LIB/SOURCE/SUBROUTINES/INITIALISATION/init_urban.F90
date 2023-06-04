SUBROUTINE init_urban()
! Description:
!   Routine to initialise URBAN parameters.
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! -----------------------------------------------------------------------------
!!!  6.1   10/01/07   First written. Peter Clark and Aurore Porson
!******************************************************************************
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Land

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE misc_utils, ONLY :  &
!  imported procedures
     checkVarPos,checkVars,read_list,repeatVal,varInfo,varList,varValue

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2dcomp

  USE file_utils, ONLY:  &
!  imported procedures
     closeFile, fileUnit, findTag, openFile  &
!  imported arrays with intent(inout)
    ,irecPrev

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatBin,formatLen,formatNc,formatPP,jinUnit,tagAscBin,tagNc  &
!  imported scalars with intent(in)
    ,echo,nxIn,nyIn  &
!  imported arrays with intent(in)
    ,mapInLand

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     routeOnly,         &
!  intent(inout)
     l_cosz

  USE nvegparm, ONLY :  &
     albsnc_nvg, albsnf_nvg, emis_nvg, z0_nvg, ch_nvg, vf_nvg

  USE c_z0h_z0m, ONLY : z0h_z0m

  USE urban_param, ONLY :                                                     &
     hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr,                    &
     a, cdz, kappa2, z0m_mat, anthrop_heat_scale

  USE switches_urban, ONLY :                                                  &
     l_urban2t, l_moruses, l_moruses_albedo, l_moruses_emissivity,            &
     l_moruses_rough, l_moruses_storage, l_moruses_storage_thin,              &
     l_moruses_macdonald, l_urban_empirical

  USE ancil_info, ONLY : land_pts, frac

  USE nstypes, ONLY : urban, urban_canyon, urban_roof, ntype, npft, ice

  IMPLICIT NONE

!------------------------------------------------------------------------------
! Local scalar parameters.
!------------------------------------------------------------------------------
  INTEGER, PARAMETER ::  &
     nvarMax = 9     !  Maximum number of MORUSES parameters

!------------------------------------------------------------------------------
! Local scalar variables.
!------------------------------------------------------------------------------
  INTEGER ::      &
     i,l          &  !  Loop counter
    ,ierr         &  !  Error flag
    ,ierrSum      &  !  Error flag
    ,inUnit       &  !  Unit used to connect to input file
    ,iorder       &  !  Work
!------------------------------------------------------------------------------
!  Every expected variable (nvarMax) must have an "iposX" variables declared
    ,iposwrr      &  !  Position in list of all expected variables of wrr
    ,iposhwr      &  !  Position in list of all expected variables of hwr
    ,iposhgt      &  !  Position in list of all expected variables of hgt
    ,iposztm      &  !  Position in list of all expected variables of ztm
    ,iposdisp     &  !  Position in list of all expected variables of disp
    ,iposalbwl    &  !  Position in list of all expected variables of albwl
    ,iposalbrd    &  !  Position in list of all expected variables of albrd
    ,iposemisw    &  !  Position in list of all expected variables of emisw
    ,iposemisr    &  !  Position in list of all expected variables of emisr
!------------------------------------------------------------------------------
    ,ivar,jvar    &  !  Work
    ,nfieldFile   &  !  Number of fields per time in a file
    ,nheaderField &  !  Number of header records before each field in file
    ,nheaderFile  &  !  Number of header records at start of file
    ,nheaderT     &  !  Number of headers at start of each time
    ,nlineField   &  !  Work
    ,nvar         &  !  Number of variables required for chosen configuration
    ,nvarFound    &  !  Number of variables found in list
    ,nvarIn       &  !  Number of variables that are read in
!                       Often this = nvar, but may be < nvar if we can derive
!                       a field from other fields, or set by some other
!                       assumption.
    ,readT        &  !  Time level to be read from file
    ,useIndex     &  !  Index in irecPrev
    ,dumZ            !  Dummy z level

  REAL ::         &
     lambdaf, lambdap         !  Frontal and planar area index

  LOGICAL ::                &
     readFile               & !  Flag indicating if another file is to be read
    ,sdfFile                & !  TRUE if dealing with an SDF
!                                  (self-describing file)
    ,errFound               & !  Flag indicating an error
    ,checkNames,checkPos    & !  Work
    ,summary                  !  Work

  CHARACTER(len = formatLen) ::  &
     fileFormat                    !  Format of file

  CHARACTER(len = 150) ::        &
     fileName                      !  The name of a file

!------------------------------------------------------------------------------
! Local array variables
!------------------------------------------------------------------------------
  INTEGER ::                    &
     ival(nvarMax)              & !  Work
    ,varStashCode(nvarMax)      & !  STASH code for each variable
    ,varUse(nvarMax)            & !  Flag indicating if a variable is required
    ,varFlag(nvarMax)           & !  Flags indicating source of data
    ,varFlagTmp(nvarMax)        & !  Work
    ,varInSort(nvarMax)         & !  A list of which variables in master list
!                                    are to be read in, sorted into ascending
!                                    order of location in file. Values 1:varIn
!                                    give a list of variables that are to be
!                                    read from file (in terms of location in
!                                    masterlist of expected variables).
    ,tile_pts(ntype)            & !  Number of land points which
                                  !      include the nth surface type.
    ,tile_index(land_pts,ntype)   !  Indices of land points which
                                  !      include the nth surface type.

  REAL ::                    &
     varConst(nvarMax)       & !  Constant value for setting urban parameters
    ,varConstTmp(nvarMax)    & !  Work
    ,tmpval(land_pts)        & !  Workspace to read in MORUSES data
    ,sc_hwr(land_pts)        & !  Working variable
    ,d_h(land_pts)             !  Working variable


  LOGICAL ::                 &
     foundVar(nvarMax)       & ! Flag indicating if a variable is listed in
!                                  run control file
    ,done(nvarMax)             ! Work


  CHARACTER(len = 150) ::    &
     varName(nvarMax)        & ! List of expected variables names
    ,varDesc(nvarMax)        & ! Description of each variable
    ,varNameTmp(nvarMax)     & ! Variable names as read in
    ,varNameSDF(nvarMax)     & ! Names of variables in input file
!                                 (SDF=self-describing file)
    ,varNameSDFTmp(nvarMax)    ! Work


!------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

! If the two-tile urban schemes are not used, make sure the appropriate
! switches are set to .FALSE. to avoid unnecessary calculations and leave.
! The urban arrays are not allocated yet so there is no need to initialise
! them to zero
  IF ( .NOT. l_urban2t ) THEN
    l_urban_empirical      = .FALSE.
    l_moruses_albedo       = .FALSE.
    l_moruses_emissivity   = .FALSE.
    l_moruses_rough        = .FALSE.
    l_moruses_storage      = .FALSE.
    l_moruses_storage_thin = .FALSE.
    l_moruses_macdonald    = .FALSE.
    RETURN
  END IF

  PRINT "(50('-'),/,a)", 'init_urban'

! First check that the run actually has urban and that urban schemes are not
! run in error

  CALL tilepts(land_pts,frac,tile_pts,tile_index)

! urban_canyon tile number is set in init_nonveg so should be present if
! two-tiles schemes are used
  IF ( l_urban2t ) THEN
    IF ( tile_pts(urban_canyon) == 0 ) THEN
      PRINT *, 'WARNING: init_urban: Number of urban land points = 0.'
      PRINT *, '  Either URBAN-2T or MORUSES has been selected in init_opts,'
      PRINT *, '  but there are no urban land points. Extra calculations may'
      PRINT *, '  be being performed that will not impact on the results.'
      PRINT *
    END IF
  END IF

! init_top followed as a template
! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_urban','>INIT_URBAN' )

  READ(jinUnit,*) l_urban_empirical,l_moruses_macdonald
  READ(jinUnit,*) l_moruses_albedo,l_moruses_emissivity,l_moruses_rough
  READ(jinUnit,*) l_moruses_storage,l_moruses_storage_thin
  READ(jinUnit,*) anthrop_heat_scale

!!   Check urban switch logic is sensible
  IF ( l_moruses ) THEN
!   Check that MORUSES has some parametrisations turned on
    IF ( .NOT. l_moruses_albedo .AND. .NOT. l_moruses_emissivity &
       .AND. .NOT. l_moruses_rough .AND. .NOT. l_moruses_storage ) THEN
      PRINT *, 'ERROR: init_urban: MORUSES has no parametrisations turned on'
      STOP
    END IF
  ELSE
!   If MORUSES is not used then ALL MORUSES switches MUST be .FALSE.
    PRINT *,'WARNING: l_moruses = .FALSE.'
    PRINT *,'Making sure individual MORUSES switches are off'
    l_moruses_albedo       = .FALSE.
    l_moruses_emissivity   = .FALSE.
    l_moruses_rough        = .FALSE.
    l_moruses_storage      = .FALSE.
    l_moruses_storage_thin = .FALSE.
    l_moruses_macdonald    = .FALSE.
  END IF

!   Check MORUSES switch logic
  IF ( l_moruses .AND. l_urban_empirical .AND.  &
     .NOT. l_moruses_macdonald ) THEN
    PRINT *, 'WARNING: If l_urban_empirical with MORUSES then the'
    PRINT *, 'MacDonald (1998) formulation, for roughness length for'
    PRINT *, 'momentum and displacement height, must be used for'
    PRINT *, 'consistency.'
    PRINT *, '  Setting l_moruses_macdonald = T'
    l_moruses_macdonald = .TRUE.
  END IF

  IF ( l_moruses_albedo .AND. .NOT. l_cosz ) THEN
    PRINT *, 'WARNING: l_moruses_albedo = .TRUE. with l_cosz = .FALSE.'
    PRINT *, 'Setting l_cosz = .TRUE.'
    l_cosz = .TRUE.
  END IF

  IF ( echo ) THEN
    PRINT '(/,a)', 'Urban switches used'
    PRINT *, 'l_urban2T             ', l_urban2T
    PRINT *, 'l_urban_empirical     ', l_urban_empirical
    PRINT *, 'l_moruses             ', l_moruses
    IF ( l_moruses ) THEN
      PRINT *, 'l_moruses_albedo      ', l_moruses_albedo
      PRINT *, 'l_moruses_emissivity  ', l_moruses_emissivity
      PRINT *, 'l_moruses_rough       ', l_moruses_rough
      PRINT *, 'l_moruses_storage     ', l_moruses_storage
      PRINT *, 'l_moruses_storage_thin', l_moruses_storage_thin
      PRINT *, 'l_moruses_macdonald   ', l_moruses_macdonald
    END IF
  END IF

! Issue warnings on logic

  IF ( l_moruses_storage_thin .AND. .NOT. l_moruses_storage ) THEN
    PRINT '(/,a)', 'WARNING: MORUSES storage parametrisation not used.'
    PRINT *, 'l_moruses_storage      = .FALSE. when'
    PRINT *, 'l_moruses_storage_thin = .TRUE.'
  END IF

! Some URBAN-2T parameters will be used by MORUSES depending on MORUSES
! switches selected. Issue warnings if this is the case.
  IF ( l_moruses ) THEN
    PRINT '(/,a)', 'WARNING: init_urban: if the list below is populated then'
    PRINT *, 'some MORUSES switches are turned off'
    PRINT *, 'MORUSES is using the following URBAN-2T parameters instead ', &
       'for canyon & roof:'
    IF ( .NOT. l_moruses_albedo ) THEN
      PRINT *,'albsnc_nvg =', &
         albsnc_nvg(urban_canyon-npft),albsnc_nvg(urban_roof-npft)
      PRINT *,'albsnf_nvg =', &
         albsnf_nvg(urban_canyon-npft),albsnc_nvg(urban_roof-npft)
    END IF
    IF ( .NOT. l_moruses_emissivity ) THEN
      PRINT *,'emis_nvg   =',   &
         emis_nvg(urban_canyon-npft),emis_nvg(urban_roof-npft)
    END IF
    IF ( .NOT. l_moruses_rough ) THEN
      PRINT *,'z0_nvg     =',z0_nvg(urban_canyon-npft),z0_nvg(urban_roof-npft)
      PRINT *,'z0h_z0m    =',z0h_z0m(urban_canyon),z0h_z0m(urban_roof)
    END IF
    IF ( .NOT. l_moruses_storage ) THEN
      PRINT *,'ch_nvg     =',ch_nvg(urban_canyon-npft),ch_nvg(urban_roof-npft)
      PRINT *,'vf_nvg     =',vf_nvg(urban_canyon-npft),vf_nvg(urban_roof-npft)
    END IF
  END IF

! Allocate & initialise arrays
  CALL allocate_arrays( 'init_urban' )

!------------------------------------------------------------------------------
! Set parameters for urban morphology & materials (mostly MORUSES)
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Set up list of expected variables.
!------------------------------------------------------------------------------

  PRINT '(/,a)', &
     'Setting urban geometry and material characteristics where appropriate'

  i = 0
  varName(:) = ''

  i = i + 1; iposwrr = i; varName(i) = 'wrr'
  varStashcode(i) = 496
  varDesc(i)      = 'Urban canyon width to repeating width ratio'

  i = i + 1; iposhwr = i; varName(i) = 'hwr'
  varStashcode(i) = 495
  varDesc(i)      = 'Urban height-to-width ratio'

  i = i + 1; iposhgt  = i; varName(i) = 'hgt'
  varStashcode(i) = 494
  varDesc(i)      = 'Urban building height'

  i = i + 1; iposztm = i; varName(i) = 'ztm'
  varStashcode(i) = 498
  varDesc(i)      = 'Urban effective roughness length'

  i = i + 1; iposdisp = i; varName(i) = 'disp'
  varStashcode(i) = 497
  varDesc(i)      = 'Urban displacement height'

  i = i + 1; iposalbwl = i; varName(i) = 'albwl'
  varStashcode(i) = 499
  varDesc(i)      = 'Urban wall albedo'

  i = i + 1; iposalbrd = i; varName(i) = 'albrd'
  varStashcode(i) = 500
  varDesc(i)      = 'Urban road albedo'

  i = i + 1; iposemisw = i; varName(i) = 'emisw'
  varStashcode(i) = 501
  varDesc(i)      = 'Urban wall emissivity'

  i = i + 1; iposemisr = i; varName(i) = 'emisr'
  varStashcode(i) = 502
  varDesc(i)      = 'Urban road emissivity'

  IF ( l_moruses ) THEN
!   All data variables are required by MORUSES
    varUse(:)      = 1
  ELSE
!   Only wrr is required by URBAN-2T
    varUse(:)       = 0
    varUse(iposwrr) = 1
  END IF

! Count number of variables required for chosen configuration.
  nvar = COUNT( varUse(:) > 0 )

! Establish where data are to be read from.
  READ(jinUnit,*) readFile

  IF ( readFile ) THEN
!------------------------------------------------------------------------------
!   An external file will be read.
!------------------------------------------------------------------------------
    READ(jinUnit,*) fileFormat
    READ(jinUnit,*) fileName

    PRINT *, fileFormat
    PRINT *, fileName

  ELSE
    fileFormat = formatAsc
  END IF

  SELECT CASE ( fileFormat )

    CASE ( formatAsc,formatBin,formatPP )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_urban',tagAscBin,preInit=.TRUE. )
      READ(jinUnit,*) nheaderFile,nheaderField

!     Read variable names, flags and constVals from a blank-delimited list.
      CALL read_list( jinUnit,3,nvarMax,'>VARS','>ENDVARS',' ','init_urban'  &
        ,nvarFound,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2     &
        ,rvar=varConstTmp,rvarPos=3 )

      sdfFile = .FALSE.

    CASE ( formatNc )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_urban',tagNc,preInit=.TRUE. )

!     Read variable names, flags, constVals and names in file, from
!     blank-delimited list.

      CALL read_list( jinUnit,4,nvarMax,'>VARS','>ENDVARS',' ','init_urban'  &
         ,nvarFound,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2    &
         ,rvar=varConstTmp,rvarPos=3,cvar2=varNameSDFTmp,cvar2Pos=4 )
      sdfFile = .TRUE.

    CASE default
      WRITE(*,*)'ERROR: init_urban: no code for fileFormat=',TRIM(fileFormat)
      STOP

  END SELECT

! Check that there are no repeated names in input list.
  IF ( repeatVal( varNameTmp(1:nvarFound) ) ) THEN
    PRINT *, 'ERROR: init_urban: repeated variable names.'
    STOP
  END IF

! Work out what variables are indicated.
  CALL varList( nvarFound,sdfFile,varNameTmp,varName,errFound,foundVar       &
     ,varFlagIn=varFlagTmp,varConstIn=varConstTmp,varNameSDFin=varNameSDFtmp &
     ,varFlag=varFlag,varConst=varConst,varNameSDF=varNameSDF )
  IF ( errFound ) THEN
    WRITE(*,*)'ERROR: init_urban: error returned by varList'
    STOP
  END IF

! Check that all necessary variables are provided.
  CALL checkVars( nvarMax,.TRUE.,.FALSE.,varUse,foundVar,varDesc,varName  &
                 ,nvarFound,varFlag,errFound )
! Stop if an error was detected.
  IF ( errFound ) THEN
    WRITE(*,*)'ERROR: init_urban: an error was found by checkVars.'
    STOP
  END IF

! Calculate number of variables to be read in (not including spatially
! constant).
  nvarIn = COUNT( varFlag(:) > 0 .AND. varUse(:) > 0 )

!------------------------------------------------------------------------------

! If no fields are to be read from external file, reset flags.
  IF ( nvarIn == 0 .AND. readFile ) THEN
    readFile   = .FALSE.
    fileFormat = formatAsc
    WRITE(*,*)'WARNING: init_urban: No fields to be read from external file,'
    WRITE(*,*)'but readFile = .TRUE.'
    WRITE(*,*)'All values will be set via flag < 0 (e.g. spatially constant).'
    WRITE(*,*)'Setting readFile to FALSE.'
  END IF

! Get a sorted list, giving input variables according to position in input
! file(s). This is needed for reading from run control file (and makes ASCII
! reading easier) and is always done.
  done(:)      = .TRUE.
  varInSort(:) = 0
! Only consider variables that are to be read in ( flag > 0 ). SDF variables
! have flag set > 0.
  WHERE ( varFlag(:) > 0 ) done(:) = .FALSE.
  DO ivar = 1, nvarIn
    ival(1:1) = MINLOC( varFlag(:), .NOT. done(:) )
    i         = ival(1)
    done(i)   = .TRUE.
    varInSort(ivar) = i
  END DO

!------------------------------------------------------------------------------
! Some output to screen.
  IF ( echo ) THEN
    DO ivar = 1, nvarMax
      IF ( varUse(ivar) > 0 ) &
         CALL varInfo( varName(ivar),varFlag(ivar),varStashcode(ivar)  &
         ,varNameSDF(ivar),fileFormat,constVal=varConst(ivar)          &
         ,varDesc=varDesc(ivar) )
    END DO
  END IF
!------------------------------------------------------------------------------

! Check we have an acceptable file format, and set flags indicating what other
! variables to check.
  checkNames = .FALSE.
  checkPos   = .FALSE.
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )
!     ASCII. Need positions of variables.
      checkPos = .TRUE.
    CASE ( formatBin )
!     GrADS. For now, need positions of variables.
      checkPos = .TRUE.
    CASE ( formatNc )
!     netCDF. Need names of variables.
      checkNames = .TRUE.
    CASE ( formatPP )
!     PP format. For now, need positions of variables.
      checkPos = .TRUE.
    CASE default
      WRITE(*,*)'ERROR: init_urban: no code for fileFormat=',TRIM(fileFormat)
      STOP
  END SELECT

  IF ( checkPos ) THEN
!   Check positions of variables within file.
!   Check for repeated locations.
!
!   First argument to checkVarPos indicates how important the order of the
!   variables is. 0 means order is not important. If reading from run control
!   file, insist that variables are 1,2,3,.... - indicated by iorder=2. Note
!   that for run control file, we interpret varFlag as variable number, not
!   field number, since this is required for current code to pass checkVarPos
!   with iorder=2 (want vars to appear adjacent, which wouldn't be if
!   had multi-level data and gave non-consecutive field numbers. Admittedly,
!   this is becoming a bit of a faff, and needs to be rewritten, one day.).
!   In all cases, values < 1 are ignored.
    iorder = 0
    IF ( .NOT. readFile ) iorder = 2

!   We need to load file locations in sorted order.
    ival(:) = 0
    DO ivar = 1, nvarIn
      ival(ivar) = varFlag( varInSort(ivar) )
    END DO

    ierr = checkVarPos( iorder,ival(:),' init_urban: ival' )
    IF ( ierr < 0 ) THEN
      WRITE(*,*)'ERROR: init_urban: error from checkvarPos.'
      WRITE(*,*)'If error was repeated use of same varPos, but you do in fact want to reuse the'
      WRITE(*,*)'same data for more than one variable, comment out this stop!'
      STOP
    END IF
  END IF

  IF ( checknames ) THEN
!   Check names only for variables that are to be read in ( varFlag > 0 ).
    DO ivar = 1, nvarMax
      DO jvar = ivar + 1, nvarMax
        IF ( ( varFlag(ivar) > 0 .AND. varFlag(jvar) > 0 ) .AND.  &
             ( varNameSDF(ivar) == varNameSDF(jvar) ) ) THEN
          WRITE(*,*)'ERROR: init_urban: repeated varNameSDF: ',TRIM(varNameSDF(ivar))
          WRITE(*,*)'If you really do want to use the same variable from file to set values'
          WRITE(*,*)'of more than one FORTRAN variable, comment out this check!'
          STOP
        END IF
      END DO
    END DO
  END IF
!------------------------------------------------------------------------------

! Open file.
  IF ( readFile ) THEN
    inUnit = fileUnit( fileFormat )  !  Get unit
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old'         &
      ,'init_urban',ncType )
  ELSE IF ( nvarIn > 0 ) THEN
    WRITE(*,*) 'Reading urban morphology and building material characteristics'
    WRITE(*,*) 'from the run control file.'
!    WRITE(*,*) 'Data must be in fields 1,2,3,.... (i.e. sequential from 1).'
    inUnit = jinUnit
    CALL findTag( inUnit,'init_urban','>DATA',preInit=.TRUE. )
  ELSE
    WRITE(*,*) 'WARNING: init_urban: No fields to be read from ',            &
       'run control file'
    WRITE(*,*) 'All values set via flag < 0 (e.g. spatially constant).'
  END IF

!------------------------------------------------------------------------------

! Read data.
! All levels of a variable are read into temporary space, then loaded into
! final variables. Only data for the required points are returned. Only need
! to do anything for variables that are to be read in, since currently all
! (required) variables must be read in (no option to set as constant).
! For non-SDF files, variables are read in the order in which they appear in
! file - this is done so that values can be read from the run control file
! (which is possibly open on stdIn and therefore can't be backspaced).

  DO ivar = 1, nvarIn

!   Get location in master list of next variable to be read.
    jvar = varInSort(ivar)

    IF ( varUse(jvar) > 0 .AND. varFlag(jvar) > 0 ) THEN

      IF ( inUnit == jinUnit .AND. nxIn*nyIn == 1 .AND. &
           fileFormat == formatAsc ) THEN
!       If the input grid is a single point and reading from run control file,
!       expect no new line between fields (eg all on one line). Calling
!       readVar means we could cope with headers in the run control file.
!       But there's no need since annotation is already simple in this case.
        READ(jinUnit,*) tmpval(:)

      ELSE

!------------------------------------------------------------------------------
! Simplifying assumptions regarding input file. Otherwise have to read in.
!------------------------------------------------------------------------------
         readT      = 1         !  Time level to read from file.
         dumZ       = 1         !  'z' level to read from file.
         nfieldFile = nvar      !  # of fields in file
         nheaderT   = 0         !  No headers at top of each time.
         nlineField = 0         !  Will not attempt to read ASCII line-by-line.

!        Set index to use for irecPrev with netCDF files - this irecPrev isn't
!        changed, but need to keep index within bounds.
         useIndex = inUnit
         IF ( fileFormat == formatNC ) useIndex = 1

!        Read data.
         CALL readVar2dComp( readT,dumZ,varFlag(jvar),varStashcode(jvar)  &
            ,irecPrev(useIndex),nfieldFile,nheaderFile,nheaderT           &
            ,nheaderField,nlineField,nxIn,nyIn,inUnit,varNameSDF(jvar)    &
            ,mapInLand(:),(/(i,i=1,land_pts)/),fileFormat,tmpval(:)      &
            ,'init_urban','init_urban',ncType )

      END IF

!     Copy data into final variable.
      IF ( jvar == iposwrr ) THEN
        wrr(:) = tmpval(:)
      ELSE IF ( jvar == iposhwr ) THEN
        hwr(:) = tmpval(:)
      ELSE IF ( jvar == iposhgt ) THEN
        hgt(:) = tmpval(:)
      ELSE IF ( jvar == iposztm ) THEN
        ztm(:) = tmpval(:)
      ELSE IF ( jvar == iposdisp ) THEN
        disp(:) = tmpval(:)
      ELSE IF ( jvar == iposalbwl ) THEN
        albwl(:) = tmpval(:)
      ELSE IF ( jvar == iposalbrd ) THEN
        albrd(:) = tmpval(:)
      ELSE IF ( jvar == iposemisw ) THEN
        emisw(:) = tmpval(:)
      ELSE IF ( jvar == iposemisr ) THEN
        emisr(:) = tmpval(:)
      ELSE
        WRITE(*,*)'ERROR: init_urban: no code to load jvar=',jvar
        STOP
      END IF

    END IF   !  Read variable

  END DO  !  Variables

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!------------------------------------------------------------------------------
! Now deal with variables that are to be set as spatially constant.
! This involves a list of variables that is similar to that used above - it
! would be much better to avoid this replication by dealing with all variables
! inside one loop, regard less of whether a constant value is to be used.
! But for now it's easier to do this....
  IF ( COUNT( varFlag(:) == -1 ) > 0 ) THEN
    DO ivar = 1, nvarMax
      IF ( varFlag(ivar) == -1 ) THEN
        IF ( ivar == iposwrr ) THEN
          wrr(:) = varConst(ivar)
        ELSE IF ( ivar == iposhwr ) THEN
          hwr(:) = varConst(ivar)
        ELSE IF ( ivar == iposhgt ) THEN
          hgt(:) = varConst(ivar)
        ELSE IF ( ivar == iposztm ) THEN
          ztm(:) = varConst(ivar)
        ELSE IF ( ivar == iposdisp ) THEN
          disp(:) = varConst(ivar)
        ELSE IF ( ivar == iposalbwl ) THEN
          albwl(:) = varConst(ivar)
        ELSE IF ( ivar == iposalbrd ) THEN
          albrd(:) = varConst(ivar)
        ELSE IF ( ivar == iposemisw ) THEN
          emisw(:) = varConst(ivar)
        ELSE IF ( ivar == iposemisr ) THEN
          emisr(:) = varConst(ivar)
        ELSE
          WRITE(*,*)'ERROR: init_urban: no code to load ivar=',ivar
          STOP
        END IF
      END IF
    END DO
  END IF

! Empirical relationships derived from correlating CEH urban fraction and
! LUCID urban geometry data for London. Obtained from collaboration with the
! University of Reading. See:
!     Bohnenstengel, S.I., Evans, S., Clark, P., Belcher, S.E. (2010);
!     Simulations of the London urban heat island, Q.J.R.Meteorol. Soc., to
!     be submitted.
! for more information

! Check for ice has been left in to be consistent with UM, but is not actually
! required here
  IF ( l_urban_empirical ) THEN
    PRINT *, 'Using empirical relationships for urban geometry: wrr'
    IF ( l_moruses ) PRINT *, 'Using empirical relationships for urban geometry: hwr'
    DO l = 1, land_pts
      IF ( frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
        lambdap = 22.878*frac(l,urban)**6 - 59.473*frac(l,urban)**5      &
           + 57.749*frac(l,urban)**4 - 25.108*frac(l,urban)**3           &
           + 4.3337*frac(l,urban)**2 + 0.1926*frac(l,urban)              &
           + 0.036
        wrr(l) = 1.0 - lambdap
        IF ( l_moruses ) THEN
          lambdaf = 16.412*frac(l,urban)**6 - 41.855*frac(l,urban)**5    &
             + 40.387*frac(l,urban)**4 - 17.759*frac(l,urban)**3         &
             + 3.2399*frac(l,urban)**2 + 0.0626*frac(l,urban)            &
             + 0.0271
          hwr(l) = 4.0 * atan(1.0)/2.0 * lambdaf / ( 1.0 - lambdap )
        END IF
      END IF
    END DO
  END IF

  IF ( l_moruses ) THEN
    PRINT *, "Setting MORUSES parameters"

    IF ( l_urban_empirical ) THEN
      PRINT *, 'Using empirical relationships for urban geometry: hgt'
      DO l = 1, land_pts
        IF (frac(l,urban) > 0.0 .AND. frac(l,ice) == 0.0 ) THEN
          hgt(l) =                                                       &
               167.409  * frac(l,urban)**5 - 337.853  * frac(l,urban)**4 &
             + 247.813  * frac(l,urban)**3 -  76.3678 * frac(l,urban)**2 &
             +  11.4832 * frac(l,urban)    +   4.48226
        END IF
      END DO
    END IF

    IF ( l_moruses_macdonald ) THEN
      !       Macdonald Formulation
      PRINT *, 'Using MacDonald formulation'
      sc_hwr(:) = 0.5 * ( hwr(:) / (2.0 * ATAN(1.0)) )
      d_h(:)    = 1.0 - wrr(:) * ( a**(wrr(:) - 1.0) )
      disp(:)   = d_h(:) * hgt(:)
      DO l = 1, land_pts
        IF ( wrr(l) > 0.0 ) THEN
          ztm(l)    = (cdz * (1.0 - d_h(l)) *                               &
             sc_hwr(l) * wrr(l) / kappa2)**(-0.5)
          ztm(l)    = (1.0 - d_h(l))*EXP(-ztm(l))
          ztm(l)    = ztm(l) * hgt(l)
          ztm(l)    = MAX(ztm(l),z0m_mat)
        END IF
      END DO
    END IF

  END IF ! l_moruses

!------------------------------------------------------------------------------
! Optional writing of fields to screen.
!------------------------------------------------------------------------------
  IF ( echo ) THEN

    DO i = 1, 2
!     i = 1 writes full field
!     i = 2 writes summary statistics for whole field
!     Don't bother summarizing single point, single level field (i.e. 1 value).
      IF ( nxIn*nyIn == 1 .AND. i == 2 ) EXIT
      IF ( i == 1 ) THEN
        summary = .FALSE.
      ELSE
        summary = .TRUE.
        WRITE(*,*)'### NB The ranges below include any ice points. ###'
      END IF

      DO ivar = 1, nvarMax

        IF ( varUse(ivar) > 0 ) THEN

          IF ( ivar == iposwrr ) THEN
            CALL varValue( summary,wrr,varFormat='f5.3'                      &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposhwr ) THEN
            CALL varValue( summary,hwr,varFormat='f6.3'                      &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposhgt ) THEN
            CALL varValue( summary,hgt,varFormat='f7.3'                      &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposztm ) THEN
            CALL varValue( summary,ztm,varFormat='f6.3'                      &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposdisp ) THEN
            CALL varValue( summary,disp,varFormat='f6.3'                     &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposalbwl ) THEN
            CALL varValue( summary,albwl,varFormat='f5.3'                    &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposalbrd ) THEN
            CALL varValue( summary,albrd,varFormat='f5.3'                    &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposemisw ) THEN
            CALL varValue( summary,emisw,varFormat='f5.3'                    &
               ,varName=varName(ivar) )

          ELSE IF ( ivar == iposemisr ) THEN
            CALL varValue( summary,emisr,varFormat='f5.3'                    &
               ,varName=varName(ivar) )

          ELSE
            WRITE(*,*)'varName=',TRIM( varName(ivar) )
            WRITE(*,*)'INIT_URBAN: no code to echo this variable to screen.'
            WRITE(*,*)'Not important - this is just diagnostic output.'

          END IF

        END IF  !  varUse

      END DO  !  ivar
    END DO  !  i

  END IF  !  echo

!------------------------------------------------------------------------------
! Expand urban tile to two tiles based on WRR: has been moved to init_ic so
! that it is compatible with triffid
!------------------------------------------------------------------------------

  RETURN

END SUBROUTINE init_urban
