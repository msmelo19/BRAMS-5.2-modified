!###############################################################################
!###############################################################################
!###############################################################################
! subroutine init_route
! Read details of any channel routing.
! The input routing grid may differ from the "main" input grid that is used
! elsewhere in the model.
!###############################################################################

  SUBROUTINE init_route

!-------------------------------------------------------------------------------
  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     nx=>row_length,ny=>rows

  USE file_utils, ONLY : &
!  imported procedures
    closeFile,findTag,fileUnit,openFile  &
!  imported arrays with intent(out)
   ,irecPrev

  USE grid_utils, ONLY :  &
!  imported procedures
     getXYPos

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatBin,formatLen,formatNc,formatPP,optLen,jinUnit,tagAscBin,tagNc   &
!  imported scalars with intent(in)
    ,echo

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error,init_count

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2d

  USE route_mod, ONLY : &
!  imported scalar parameters
     routeTypeTRIP  &
!  imported scalars with intent(out)
    ,npRoute,nxRoute,nxInRoute,nyRoute,nyInRoute,routeCount,routeDlat  &
    ,routeDlon,routeLat1,routeLon1,routeMeander  &
    ,routeRegLatLon,routeRegrid,routeSpeed,routeStep,routeTimestep,routeType  &
    ,yrevInRoute  &
!  imported arrays with intent(out)
    ,mapInRoute,routeIndex,routeMask,routeNext,routeOrder

  USE switches, ONLY : &
!  imported scalars with intent(in)
     route

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     regDlat,regDlon,regLatLon,subRegLatLon,timeStep  &
!  imported arrays with intent(in)
    ,dateRun,latitude,longitude,timeRun

  USE timeConst, ONLY : &
     iSecInDay

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar parameters.

  INTEGER, PARAMETER ::  &
   stashCodeDir = 0   &!  STASH code for flow direction - not known
  ,stashCodeSeq = 0    !  STASH code for flow sequence - not known

  LOGICAL, PARAMETER :: byteSwap = .FALSE.   !  FALSE means that byte order of
!                    data in a binary file is not reversed after reading

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------

  INTEGER ::  &
    fieldDir     &!  location (field number) of flow direction variable in input
   ,fieldSeq     &!  location (field number) of sequence variable in input
   ,i            &!  work
   ,ierr         &!  error value
   ,ierrSum      &!  work
   ,inUnit       &!  unit connected to file
   ,ip           &!  work
   ,iseq         &!  work
   ,ix           &!  loop counter
   ,ixr          &!  work
   ,iy           &!  loop counter
   ,iyr          &!  work
   ,iyy          &!  work
   ,maxRouteInY  &!  work
   ,minRouteInY  &!  work
   ,nfieldFile   &!  number of fields per time in a file
   ,nheaderField &!  number of header records before each field in file
   ,nheaderFile  &!  number of header records at start of file
   ,nheaderT     &!  number of headers at start of each time
   ,nlineField   &!  work
   ,nsect        &!  number of sections in which grid is read
   ,readT        &!  time level to be read from file
   ,readZ        &!  'z' level to be read from file
   ,routeSeqMax  &!  maximum value of input sequence field
   ,useIndex      !  index in irecPrev

  REAL ::  &
    dx,dy            &!  work
   ,latMax           &!  work
   ,latMin           &!  work
   ,latVal           &!  work
   ,lonMaxVal        &!  work
   ,lonMinVal        &!  work
   ,lonVal,lonVal2   &!  work
   ,routeLat1In      &!  latitude of point (1,1) in input routing grid
   ,routeLon1In       !  longitude of point (1,1) in input routing grid

  LOGICAL ::   &
    globalLon   &!  T means routing grid covers 360 deg of longitude
   ,globalLonIn &!  T means input routing grid covers 360 deg of longitude
   ,readFile    &!  flag indicating if another file is to be read
   ,route360    &!  T means longitudes of routing grid are given in range 0 to 360 deg
!                   F means in range -180 to 180
   ,useFullRoute !  T means the full input grid is to be used as the routing grid
!                   F means the smallest possible part of the input grid (consistent
!                    with covering full model domain) is used as routing grid

  CHARACTER(len=formatLen) ::  &
    fileFormat   !  format of file

  CHARACTER(len=optLen) ::  &
    flowDirType   !   indicates coding scheme used for flow direction data

  CHARACTER(len=150) ::  &
    fileName    &!  the name of a file
   ,varNameDir  &!  name of variable used for input direction field
   ,varNameSeq   !  name of variable used for input sequence field
!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER :: maxRouteInX(2)  !  work
  INTEGER :: minRouteInX(2)  !  work

  INTEGER, ALLOCATABLE ::  &
    gridFlowDir(:,:)      &!  flow direction on routing grid
   ,gridIndexG(:,:)       &!  point number (i.e. location) in routing grid of each "active"
!                               land point (i.e. each point where routing required).
!                               This is zero if not land.
   ,gridNextPointG(:,:)   &!  location in routing grid of the next downstream point
   ,gridSeq(:,:)           !  routing sequence on grid. Point with gridSeq=1 are routed first,
!                               then points with gridSeq=2, and so on. This ensures that the
!                               inflow to a gridbox is known before the outflow is calculated.

  REAL :: lonMax(2)    !  work
  REAL :: lonMin(2)    !  work

  REAL, ALLOCATABLE ::  &
    gridWork(:,:)       !  work space

  LOGICAL, ALLOCATABLE ::  &
    gridMask(:,:)  !  T at land points on routing grid

!------------------------------------------------------------------------------------------

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_route'

! If no routing selected, nothing to do here.
  IF ( .NOT. route ) THEN
    IF ( echo ) WRITE(*,*)'No runoff routing has been selected.'
    RETURN
  ENDIF

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_route','>INIT_ROUTE' )

! Read choice of routing model.
  READ(jinUnit,*) routeType

! Check that routing type is recognised.
  SELECT CASE ( routeType )
    CASE ( routeTypeTRIP )
!     OK
    CASE default
      WRITE(*,*)'ERROR: init_route: do not recognise routing type (routeType)=',TRIM(routeType)
      STOP
  END SELECT

! For now, insist that the routing grid is a regular lat/lon grid.
! This can be relaxed in future versions of the code.
  routeRegLatLon = .TRUE.

! Establish what the input routing grid is.
  READ(jinUnit,*) nxInRoute,routeLon1In,routeDlon
  READ(jinUnit,*) nyInRoute,routeLat1In,routeDlat
  READ(jinUnit,*) useFullRoute

! Read timestep and parameters for routing.
  READ(jinUnit,*) routeTimestep
  READ(jinUnit,*) routeSpeed,routeMeander

! Establish coding used for flow direction.
  READ(jinUnit,*) flowDirType

  READ(jinUnit,*) yrevInRoute

! Establish how data will be input.
  READ(jinUnit,*) readFile

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.

  IF ( readFile ) THEN

!   An external file will be read.
    READ(jinUnit,*) fileFormat
    READ(jinUnit,*) fileName

    SELECT CASE ( fileFormat )

      CASE ( formatAsc,formatBin,formatPP )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_route',tagAscBin,preInit=.TRUE. )
        READ(jinUnit,*) nheaderFile,nheaderField
        READ(jinUnit,*) fieldDir,fieldSeq

      CASE ( formatNc )
!         Locate the information in run control file.
          CALL findTag( jinUnit,'init_route',tagNc,preInit=.TRUE. )
          READ(jinUnit,*) varNameDir,varNameSeq

      CASE default
        WRITE(*,*)'ERROR: init_route: no code for fileFormat=',TRIM(fileFormat)
        STOP

    END SELECT

  ELSE   !  NOT readFile

!   Data will be read from run control file.
!   The first fields will be read, and no headers are expected. Field numbers
!   are redundant for stdIn, but are used to set nfieldFile.
    fileFormat = formatAsc
    nheaderFile = 0
    nheaderField = 0
    fieldDir = 1
    fieldSeq = 2

  ENDIF

!-------------------------------------------------------------------------------
! Check that the routing timestep satisfies various conditions, including that it
! is a multiple of main timestep. Also convert to number of "main" timesteps.
!-------------------------------------------------------------------------------

  IF ( routeTimeStep >= 0 ) THEN

    IF ( routeTimeStep == 0 ) THEN
!     Use main timestep.
      routeTimeStep = NINT( timeStep )
    ELSE
!     Impose a limit of 30 days - makes it easier to initialise counter!
      IF ( routeTimeStep > 30*iSecInDay ) THEN
        WRITE(*,*)'ERROR: init_route: constant timestep must be <= 30 days.'
        WRITE(*,*)'Timestep=',REAL(routeTimeStep)/REAL(iSecInDay),' days.'
        STOP
      ENDIF
      IF ( MOD( routeTimeStep,NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: init_route: routing timestep  must be a multiple of main timestep'
        STOP
      ENDIF
    ENDIF

!   Insist that periods can be synchronised (kept in phase with) with days.
!   Periods <= 1 day must fit into 1 day, periods > 1 day must be multiples of a day.
!   This is not essential to the rest of the code (at present!!), but seems sensible.
    IF ( routeTimeStep <= iSecInDay ) THEN
      IF ( MOD( iSecInDay,routeTimeStep ) /= 0 ) THEN
        WRITE(*,*)'ERROR: one day must be a multiple of period for routing model (for period <= 1day)'
        WRITE(*,*)'init_route'
        STOP
      ENDIF
    ELSE
      IF ( MOD( routeTimeStep,iSecInDay ) /= 0 ) THEN
        WRITE(*,*)'ERROR: periods longer than one day must be a multiple of days.'
        WRITE(*,*) 'init_route'
        STOP
      ENDIF
    ENDIF

!   Convert to number of timesteps.
    routeTimeStep = routeTimeStep / NINT(timeStep)

  ELSE  !  routeTimeStep<0

    WRITE(*,*) 'ERROR: init_route: routeTimeStep must be >= 0.'
    STOP

  ENDIF

!-------------------------------------------------------------------------------
! Check we have positive speed and meander.
!-------------------------------------------------------------------------------
  IF ( routeSpeed<=0.0 .OR. routeMeander<0.0 ) THEN
    WRITE(*,*)'ERROR: init_route: speed and meander ratio must be >0.'
    WRITE(*,*)'routeSpeed=',routeSpeed,' routeMeander=',routeMeander
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Check the specification of the input grid.
!-------------------------------------------------------------------------------

! Establish what range of longitude is used for routing longitude: -180 to 180, or 0 to 360.
  route360 = .TRUE.
  lonVal = 0.0
  lonVal2 = 360.0
  IF ( routeLon1In < 0.0 ) THEN
    route360=.FALSE.
    lonVal = -180.0
    lonVal2 = 180.0
  ENDIF

! Check lat and lon values are within range and grid is not unfeasibly large.
! Note that this will not allow a grid with rows at +/-90.
  IF ( routeLat1In-0.5*routeDlat < -90.0-EPSILON(routeDlat) .OR.  &
       routeLat1In+(REAL(nyInRoute-1)+0.5)*routeDlat > 90.0+EPSILON(routeDlat) ) THEN
    WRITE(*,*)'ERROR: init_route: edge of routing grid is at latitude<=-90 or >=90 degN'
    WRITE(*,*)'Have edges at ',routeLat1In-0.5*routeDlat  &
               ,routeLat1In+(REAL(nyInRoute-1)+0.5)*routeDlat
    STOP
  ENDIF
  IF ( REAL(nyInRoute)*routeDlat > 180.0+EPSILON(routeDlat) ) THEN
    WRITE(*,*)'ERROR: init_route: latitudinal extent of routing grid is > 180 degrees'
    WRITE(*,*)'Have extent=',REAL(nyInRoute)*routeDlat
    STOP
  ENDIF
  IF ( routeLon1In-0.5*routeDlon < lonVal-EPSILON(lonVal) .OR.  &
       routeLon1In+(REAL(nxInRoute-1)+0.5)*routeDlon > lonVal2+EPSILON(lonVal2) ) THEN
    WRITE(*,*)'ERROR: init_route: longitude of edge of routing grid must be in'
    WRITE(*,*)'range ',lonVal,' to ',lonVal2,'degE.'
    WRITE(*,*)'Have edges at ',routeLon1In-0.5*routeDlon  &
                              ,routeLon1In+(REAL(nxInRoute-1)+0.5)*routeDlon
    STOP
  ENDIF
  IF ( REAL(nxInRoute)*routeDlon > 360.0+EPSILON(routeDlon) ) THEN
    WRITE(*,*)'ERROR: init_route: longitudinal extent of routing grid is > 360 degrees'
    WRITE(*,*)'Have extent=',REAL(nxInRoute)*routeDlon
    STOP
  ENDIF

! Establish if the input routing grid is global, in longitudinal direction only.
  globalLonIn = .FALSE.
  IF ( ABS(REAL(nxInRoute)*routeDlon-360.0) < EPSILON(routeDlon) ) globalLonIn=.TRUE.

!###############################################################################

!-------------------------------------------------------------------------------
! Set up the routing grid.
!-------------------------------------------------------------------------------

! For now, only allow "regular lat/lon" grids. The model grid can be a subset
! of a regular grid - e.g. land points.
  IF ( .NOT. subRegLatLon .OR. .NOT. routeRegLatLon ) THEN
    WRITE(*,*)'ERROR: init_route: code only exists for regular lat/lon grids.'
    WRITE(*,*)'Sorry. Maybe you''d like to write some new code...?'
    STOP
  ENDIF

! Decide if model and routing grids are identical.
  routeRegrid = .FALSE.
! Are both grids "regular"? Note that a regular grid uses the standard points order.
  IF ( .NOT.regLatLon .OR. .NOT.routeRegLatLon ) routeRegrid = .TRUE.
! Is the gridbox size the same in both grids?
  IF ( routeDlon /=regDlon .OR. routeDlat==regDlat  ) routeRegrid = .TRUE.
! Are origins "consistent"? Check if (1,1) on model grid is an integer number of
! gridboxes away from (1,1) on routing grid.
  dy = ( latitude(1,1) - routeLat1In ) / routeDlat   !  number of gridboxes separating points
  dx = ( longitude(1,1) - routeLon1In ) / routeDlon
  IF ( ABS(dx-REAL(NINT(dx)))>EPSILON(dx) .OR.  &
       ABS(dy-REAL(NINT(dy)))>EPSILON(dy) ) routeRegrid = .TRUE.

!-------------------------------------------------------------------------------
! Work out what part of the input routing grid is covered by the model grid.
!-------------------------------------------------------------------------------
  nsect = 1
  IF ( useFullRoute ) THEN

!   Whole routing grid is to be used - take as a single section of input grid.
    nxRoute = nxInRoute
    nyRoute = nyInRoute
!   For now, assume that grid will be used as presented, i.e. same origin.
    routeLat1 = routeLat1In
    routeLon1 = routeLon1In
    minRouteInX(1) = 1
    maxRouteInX(1) = nxInRoute
    minRouteInY = 1

  ELSE

!   Work out what part of the input routing grid is covered by the model grid.

!   Get bounds of model domain (outside edges of gridboxes).
!   Latitude - no cyclic boundary condition allowed.
    latMin = MINVAL( latitude(:,:) ) - 0.5*regDlat
    latMax = MAXVAL( latitude(:,:) ) + 0.5*regDlat

!   Check that required area all lies within the routing grid.
    IF ( latMin < routeLat1In-0.5*routeDlat-EPSILON(routeLat1In) .OR.  &
         latMax > routeLat1In+(REAL(nyInRoute-1)+0.5)*routeDlat+EPSILON(routeLat1In) ) THEN
      WRITE(*,*)'ERROR: init_route: model domain extends beyond edge of routing grid.'
      WRITE(*,*)'N edge of model grid=',latMax
      WRITE(*,*)'N edge of input routing grid=',routeLat1In+(REAL(nyInRoute-1)+0.5)*routeDlat
      WRITE(*,*)'S edge of model grid=',latMin
      WRITE(*,*)'S edge of input routing grid=',routeLat1In-0.5*routeDlat
      STOP
    ENDIF
!   Locate the routing gridbox that covers each extreme and work out size
!   of routing grid.
!   Locate S edge - first gridbox with N edge > latMin.
    minRouteInY = 0
    latval = routeLat1In + 0.5*routeDlat   !  N edge of 1st gridbox
    DO iy=1,nyInRoute
      IF ( latval > latMin ) THEN
        minRouteInY = iy
        EXIT
      ENDIF
      latval = latval + routeDlat
    ENDDO
!   Locate N edge - first gridbox with N edge >= latMax.
    maxRouteInY = 0
    latval = routeLat1In + 0.5*routeDlat   !  N edge of 1st gridbox
    DO iy=1,nyInRoute
      IF ( latval >= latMax ) THEN
        maxRouteInY = iy
        EXIT
      ENDIF
      latval = latval + routeDlat
    ENDDO
    nyRoute = maxRouteInY - minRouteInY + 1
    routeLat1 = routeLat1In + REAL(minRouteInY-1)*routeDlat

!   Longitude - cyclic boundary condition may be invoked.
!   Express longitude in same range as for routing grid.
    lonMin(1) = MINVAL( longitude(:,:) ) - 0.5*regDlon
    lonMinVal = lonMin(1)
    lonMax(1) = MAXVAL( longitude(:,:) ) + 0.5*regDlon
    IF ( route360 ) THEN
      IF ( lonMin(1) < 0.0 ) THEN
        lonMinVal = lonMin(1)
        lonMin(1)  = 360.0 + lonMin(1)
      ENDIF
      IF ( lonMax(1)  < 0.0 ) lonMax(1)  = 360.0 + lonMax(1)
    ENDIF
    lonMaxVal = lonMax(1)

!   Note: If input routing grid is cyclic (covers all longitudes), the
!   model grid will cross at most one edge of the routing grid, and points
!   that are "cut off" by crossing an edge will definitely "reappear" at the
!   other edge of the grid.

    IF ( lonMinVal < routeLon1In-0.5*routeDlon-EPSILON(routeDlon) .OR.  &
         lonMax(1) > routeLon1In+(REAL(nxInRoute-1)+0.5)*routeDlon+EPSILON(routeDlon) ) THEN
      IF ( .NOT. globalLonIn ) THEN
        WRITE(*,*)'ERROR: init_route: model domain extends beyond edge of routing grid.'
        WRITE(*,*)'LHS of model grid=',lonMin(1)
        WRITE(*,*)'LHS of input routing grid=',routeLon1In-0.5*routeDlon
        WRITE(*,*)'RHS of model grid=',lonMax(1)
        WRITE(*,*)'RHS of input routing grid=',routeLon1In+(REAL(nxInRoute-1)+0.5)*routeDlon
        STOP
      ENDIF

!     Use cyclic boundary condition.
      nsect = 2

!     Routing grid crosses edge of input grid, but input is cyclic.
!     LHS of routing grid comes from RHS of input grid.
      lonMax(1) = routeLon1In + (REAL(nxInRoute-1)+0.5)*routeDlon
!     RHS of routing grid comes from LHS of input grid.
      lonMin(2) = routeLon1In - 0.5*routeDlon
      lonMax(2) = lonMaxVal
      nxRoute = 0
!      routeLon1 = routeLon1In + real(minRouteInX(1)-1)*routeDlon

    ENDIF   !  crosses x edge

!   Locate the routing gridbox that covers each extreme and work out size
!   of routing grid.
    DO i=1,nsect
!     Locate W edge - first gridbox with E edge > lonMin.
      minRouteInX(i) = 0
      lonval = routeLon1In + 0.5*routeDlon   !  E edge of 1st gridbox
      DO ix=1,nxInRoute
        IF ( lonval > lonMin(i) ) THEN
          minRouteInX(i) = ix
          EXIT
        ENDIF
        lonval = lonval + routeDlon
      ENDDO
!     Locate E edge - first gridbox with E edge >= lonMax.
      maxRouteInX(i) = 0
      lonval = routeLon1In + 0.5*routeDlon   !  E edge of 1st gridbox
      DO ix=1,nxInRoute
        IF ( lonval >= lonMax(i) ) THEN
          maxRouteInX(i) = ix
          EXIT
        ENDIF
        lonval = lonval + routeDlon
      ENDDO
      nxRoute = nxRoute + maxRouteInX(i) - minRouteInX(i) + 1
    ENDDO

    routeLon1 = routeLon1in + REAL( minRouteInX(1) - 1 ) * routeDlon

!   If input routing grid is cyclic in longitude, check if a smaller grid can
!   be used.
!    if ( globalLon ) call init_route_smallgrid()

  ENDIF  !  useFullRoute

! If present status indicates that regridding will not be used, check that
! grid sizes are equal.
  IF ( .NOT.routeRegrid .AND. (nx/=nxRoute .OR. ny/=nyRoute ) ) routeRegrid = .TRUE.

! Establish if routing grid can use cyclic BC in longitude.
  globalLon = .FALSE.
  IF ( ABS(REAL(nxRoute)*routeDlon-360.0) < EPSILON(routeDlon) ) globalLon=.TRUE.
  PRINT*,'globalLon=',globalLon,' nxRoute=',nxRoute
  PRINT*,'routeRegrid=',routeRegrid

! Grids can be effectively identical but with different "base" longitude,
! meaning regridding is required. See if the routing grid can be shifted
! (columns reordered) so that regridding is not required. Only do this if
! the routing grid is currently global (the most likely remaining case, I think).
!  endShift = 0
  IF ( routeRegrid .AND. globalLon .AND. nxRoute==nx .AND. nyRoute==ny .AND.  &
       routeRegLatLon .AND. regLatLon .AND. routeDlat==regDlat .AND. &
       routeDlon==regDlon .AND. routeLat1==latitude(1,1) .AND.  &
       routeLon1/=longitude(1,1) ) THEN
!   Calculate how many columns are to be removed from the R edge of the
!   routing grid and moved to L edge.
!    endShift = nint( (routeLon1-longitude(1,1))/routeDlon )
!   Set new base longitude.
    routeLon1 = longitude(1,1)
    PRINT*,'Reset routeLon1=',routeLon1
  ENDIF

! Check values are in bounds.
  IF ( nxRoute<1 .OR. nxRoute>nxInRoute .OR. nyRoute<1 .OR. nyRoute>nyInRoute ) THEN
    WRITE(*,*)'ERROR: init_route: error in limits of routing grid.'
    WRITE(*,*)'nxRoute,nyRoute=',nxRoute,nyRoute
    WRITE(*,*)'Input grid has size=',nxInRoute,'*',nyInRoute
    STOP
  ENDIF

!###############################################################################

! Allocate space for mappings for routing grid.
  CALL allocate_arrays( 'init_route 1' )

! Get mapping from input routing grid to routing grid, i.e. for each point on
! routing grid get point to be used from input routing grid.
  mapInRoute(:) = 0
  WRITE(*,*)'minRouteInX/Y=',minRouteInX(:),minRouteInY
  WRITE(*,*)'maxRoouteInx=',maxRouteInX(:)
  WRITE(*,*)'nxyRoute=',nxRoute,nyRoute,' size mapInRoute=',SIZE(mapInRoute)
  ip = 0
  DO iy=1,nyRoute
    IF ( .NOT. yrevInRoute ) THEN
      iyy = minRouteInY + iy - 2
    ELSE
      iyy = nyInRoute - ( minRouteInY + iy -2 )
    ENDIF
!   We use nsect sections from the input grid.
    DO i=1,nsect
      DO ix=minRouteInX(i),maxRouteInX(i)
        ip = ip + 1
        mapInRoute(ip) = iyy*nxInRoute + ix
      ENDDO
    ENDDO
  ENDDO

  PRINT*,'Range of mapinRoute=',MINVAL(mapInRoute),MAXVAL(mapInRoute)

!###############################################################################

! Allocate space for local workspace.
  ierrSum = 0
  ALLOCATE( gridIndexG(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( gridNextPointG(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( gridFlowDir(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( gridSeq(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( gridWork(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( gridMask(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum/=0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_route: work 2' )
    STOP
  ENDIF

!###############################################################################
! Read the input.

! Check that each input field uses different data from file. Currently trivial.
  IF ( readFile ) THEN
    SELECT CASE ( fileFormat )
      CASE ( formatAsc,formatBin,formatPP )
        IF ( fieldDir == fieldSeq ) THEN
          WRITE(*,*)'ERROR: init_route: repeated use of data'
          WRITE(*,*)'Each variable must use different data from file.'
          STOP
        ENDIF
      CASE ( formatNc )
!       Check that variable names are different.
        IF ( varNameDir == varNameSeq ) THEN
          WRITE(*,*)'ERROR: init_grid_latlon: repeated variable name.'
          WRITE(*,*)'Each variable must use different data from file.'
          STOP
        ENDIF
    END SELECT
  ENDIF

! Open file.
  IF ( readFile ) THEN
!   Get unit
    inUnit = fileUnit( fileFormat )
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_route',ncType )
  ELSE
    WRITE(*,*)'Reading data from the run control file.'
    WRITE(*,*)'Fields must be in the correct order, and be the first fields in file.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_grid latlon','>DATA' )
  ENDIF

! Simplifying assumptions regarding input file. Otherwise have to read these in.
  readT = 1      !   time level to read from file
  readZ = 1      !   'z' level to read from file
  nfieldFile = MAX(fieldDir,fieldSeq)  !  # of fields in file. Use max of fields needed - OK while readT=1
  nheaderT = 0    !  no headers at top of each time
  nlineField = 0  !  will not attempt to read ASCII line-by-line

! Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
! but need to keep index within bounds.
  useIndex = inUnit
  IF ( fileFormat == formatNC ) useIndex = 1

! Read data.
! Read flow direction.
  CALL readVar2d( readT,readZ,fieldDir,stashCodeDir,irecPrev(useIndex)  &
                 ,nfieldFile,nheaderFile,nheaderT,nheaderField  &
                 ,nlineField,nxInRoute,nyInRoute,inUnit,varNameDir  &
                 ,mapInRoute(:),(/ (ip,ip=1,nxRoute*nyRoute) /),fileFormat  &
                 ,gridWork(:,:),byteSwap,'init_route','init_route',ncType )
! Convert to integer. Note that this may give an arithmetic exception if there
! are large numbers which are beyond range of integers (e.g. large missing data value!).
  gridFlowDir(:,:) = NINT( gridWork(:,:) )

! Read sequence field.
  CALL readVar2d( readT,readZ,fieldSeq,stashCodeSeq,irecPrev(useIndex)  &
                 ,nfieldFile,nheaderFile,nheaderT,nheaderField  &
                 ,nlineField,nxInRoute,nyInRoute,inUnit,varNameSeq  &
                 ,mapInRoute(:),(/ (ip,ip=1,nxRoute*nyRoute) /),fileFormat  &
                 ,gridWork(:,:),byteSwap,'init_route','init_route',ncType )

! Convert to integer. Note that this may give an arithmetic exception if there
! are large numbers which are beyond range of integers (e.g. large missing data value!).
  gridSeq(:,:) = NINT( gridWork(:,:) )

  WRITE(*,*)'Range of gridFlowDir=',MINVAL(gridFlowDir(:,:)),MAXVAL(gridFlowDir(:,:))
  WRITE(*,*)'Range of gridSeq=',MINVAL(gridSeq(:,:)),MAXVAL(gridSeq(:,:))

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!###############################################################################

! Create temporary (gridded) copies of index and downstream point index.
  CALL init_route_map( globalLon,flowDirType,gridFlowDir  &
                      ,gridIndexG,gridMask,gridNextPointG )

  IF ( npRoute < 1 ) THEN
    WRITE(*,*)'ERROR: init_route: npRoute=0. There are no locations at which routing will be done.'
    STOP
  ENDIF

! Allocate space.
  CALL allocate_arrays( 'init_route 2' )

!-------------------------------------------------------------------------------
! Load variables from grid into vector variables.
! Order of vector is found by scanning across grid.
  ip = 0
  DO iy=1,nyRoute
    DO ix=1,nxRoute
      IF ( gridMask(ix,iy) ) THEN
        ip = ip + 1
        routeIndex(ip) = (iy-1)*nxRoute + ix
        routeNext(ip) = gridNextPointG(ix,iy)
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------

! Get maximum value of sequence number.
  routeSeqMax = MAXVAL( gridSeq(:,:) )

! Get order in which to take points.
  routeOrder(:) = 0
  i = 0
! Loop over sequence
  DO iseq=1,routeSeqMax
!   Loop over routing points.
    DO ip=1,npRoute
!     Get location on routing grid
      CALL getXYPos( routeIndex(ip),nxRoute,nyRoute,ix,iy )
      IF ( gridSeq(ix,iy) == iseq ) THEN
!       Increment counter.
        i = i + 1
        routeOrder(i) = ip
      ENDIF
    ENDDO
  ENDDO

!###############################################################################

! Some checks that we've done everything correctly, and that sequence field
! is correct.

  DO i=1,npRoute

!   Get location in index vector.
    ip = routeOrder(i)

!   Get location on routing grid
    CALL getXYPos( routeIndex(ip),nxRoute,nyRoute,ix,iy )

!   Check the prescription of any downsteam point.
    IF ( routeNext(ip) > 0 ) THEN

!     Get location on routing grid of downstream point.
      CALL getXYPos( routeNext(ip),nxRoute,nyRoute,ixr,iyr )

!     Check that the next downstream point is not of lower sequence.
!     We check using the work variable gridSeq, although this is not used later in the model.
      IF ( gridSeq(ixr,iyr) < gridSeq(ix,iy) ) THEN
        WRITE(*,*)'ERROR: init_route: routeSeq appears to be incorrect.'
        WRITE(*,*)'Next downstream point must not have lower sequence (order) number.'
        WRITE(*,*)'ip=',ip,' ix,iy=',ix,iy,' lat/lon=',routeLat1+REAL(iy-1)*routeDlat  &
             ,routeLon1+REAL(ix-1)*routeDlon,' gridSeq=',gridSeq(ix,iy)
        WRITE(*,*)'Next point=',routeNext(ip),' ix,iy=',ixr,iyr,' lat/lon='  &
                   ,routeLat1+REAL(iyr-1)*routeDlat,routeLon1+REAL(ixr-1)*routeDlon  &
                  ,' gridSeq=',gridSeq(ixr,iyr)
        WRITE(*,*)    !    blank
        WRITE(*,*)'Although this is indeed an error, you could, feasibly, continue!'
        WRITE(*,*)'I only suggest this because the sequence field is incorrect in'
        WRITE(*,*)'some UM ancillary files (as of 2006). The routing pathways (river'
        WRITE(*,*)'directions are not affected, but the points will not be taken in'
        WRITE(*,*)'downstream order as expected. I don''t know how large an error this'
        WRITE(*,*)'introduces, but I suspect it''s not "first order"!'
        WRITE(*,*)'It''s your decision....'
        STOP
      ENDIF

    ENDIF

  ENDDO  !  i  (points)

!-------------------------------------------------------------------------------
! Initialise counters.
  CALL init_count( routeTimeStep,dateRun(1),timeRun(1),routeStep )
  routeCount = 0
  PRINT*,'INIT routeStep=',routeStep

!-------------------------------------------------------------------------------

! Copy gridMask to saved space.
  routeMask(:,:) = gridMask(:,:)

! Deallocate work space.
  ierrSum = 0
  DEALLOCATE( gridIndexG,gridNextPointG,stat=ierr ); ierrSum=ierrSum+ierr
  DEALLOCATE( gridFlowDir,gridSeq,gridWork,stat=ierr ); ierrSum=ierrSum+ierr
  DEALLOCATE( gridMask,stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) CALL allocate_error( 'dealloc',ierrSum,'init_route' )

  IF ( echo ) THEN
    WRITE(*,*)  !   blank
    WRITE(*,*)'routeTimestep=',routeTimestep*NINT(timeStep),' seconds'
    WRITE(*,*)'Routing grid has size nx,ny=',nxRoute,nyRoute
    WRITE(*,*)'and gridbox SIZE (delta long,lat)=',routeDlon,routeDlat
    WRITE(*,*)'Gridpoint in "SW" corner of routing grid is at '  &
        ,routeLat1,routeLon1,' (degN,degE)'
    IF ( routeRegrid ) THEN
      WRITE(*,*)'Model grid will be regridded to routing grid.'
    ELSE
      WRITE(*,*)'No regridding required: model and routing grids are indentical.'
    ENDIF
    WRITE(*,*)'Number of land points on routing grid=',npRoute
    WRITE(*,*)'routeSeqMax=',routeSeqMax
    WRITE(*,*)'mapInRoute=',mapInRoute(:)
  ENDIF

  END SUBROUTINE init_route
!###############################################################################
!###############################################################################
!###############################################################################

! subroutine init_route_map
! Interpret the coded flow direction to identify next downstream point.
! Identify "active" or land points on routing grid.

  SUBROUTINE init_route_map( globalLon,flowDirType,gridFlowDir  &
                            ,gridIndexG,gridMask,gridNextPointG )

  USE inout, ONLY :  &
!  imported scalar parameters
     optLen

  USE grid_utils, ONLY :  &
!  imported procedures
     getXYPos

  USE route_mod, ONLY :  &
!  imported array parameters
     flowDirDelta   &
!  imported scalars with intent(in)
    ,nxRoute,nyRoute,routeDlat,routeDlon,routeLat1,routeLon1  &
!  imported scalars with intent(out)
    ,npRoute  &
!  imported arrays with intent(out)
    ,flowDirSet

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar parameters.

!  integer, parameter :: noDefDir = 9   ! index of "no defined flow direction"
!            in flowDirSet. Expect 9 ( ubound( flowDirSet(:) ).

  CHARACTER(len=optLen), PARAMETER ::  &
    flowDirTypeTRIP = 'trip'   !  value of flowDirType indicating that TRIP
!                                   coding sequence is used (see below)

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
  LOGICAL, INTENT(in) ::  &
    globalLon   !  T means routing grid covers 360 deg of longitude

  CHARACTER(len=*), INTENT(in) ::   &
   flowDirType    !  indicates scheme used to code flow direction

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
  INTEGER, INTENT(in) ::  &
    gridFlowDir(nxRoute,nyRoute)  !  flow direction (coded)

!-------------------------------------------------------------------------------
! Array arguments with intent(out)
  INTEGER, INTENT(out) ::  &
    gridIndexG(nxRoute,nyRoute)      &!  location (in routing grid) of each "active"
!                                      point (i.e. each point where routing required).
!                                      This is zero if not an active point.
   ,gridNextPointG(nxRoute,nyRoute)   !  location in routing grid of the next
!                                         downstream point.
!                                         If there is no data or undefined flow direction,
!                                         this is zero.

  LOGICAL, INTENT(out) ::  &
    gridMask(nxRoute,nyRoute)        !  T at land points on routing grid
!-------------------------------------------------------------------------------
! Local scalar variables.
  INTEGER ::  &
    dirErr  &!  used to detect error in direction field
   ,idir    &!  loop counter
   ,ip      &!  work
   ,ix      &!  loop counter
   ,ixNext  &!  x-location of next downstream point
   ,iy      &!  loop counter
   ,iyNext   !  y-location of next downstream point

!-------------------------------------------------------------------------------
! Initialise.
!-------------------------------------------------------------------------------
  ip = 0
! Set indices to zero. This value is retained at points that have the "no data"
! (sea) value.
  gridIndexG(:,:) = 0
  gridNextPointG(:,:) = 0
  gridMask(:,:)=.FALSE.

!-------------------------------------------------------------------------------
! Get the values of flow direction for the coding scheme indicated.
!-------------------------------------------------------------------------------
  SELECT CASE ( flowDirType )

    CASE ( flowDirTypeTRIP )
!     TRIP: 1-8, clockwise from N.
!     No data=0. Undefined direction=9 (as for original 1deg data; For 0.5deg
!     data, this value was 12 - in which case change 9 to 12 below!).
      flowDirSet(:) = (/ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 /)

!    case ( flowDirNim )
!      Nimrod: 1-8, clockwise from E

!    case ( flowDirArc )
!      Arc: 1-128, clockwise from E

    CASE default
      WRITE(*,*)'ERROR: init_route_map: DO not recognise routeType=',TRIM(flowDirType)
      STOP

  END SELECT

!-------------------------------------------------------------------------------
! Check that coding scheme does not have values < 0.
!-------------------------------------------------------------------------------
  IF ( MINVAL( flowDirSet(:) ) < 0 ) THEN
    WRITE(*,*)'ERROR: init_route_map: minimum of flowDirSet < 0'
    WRITE(*,*)'Flow direction codes must be >= 0 for current code.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Check that "real" flow directions (see definition above) are coded with values >0.
!-------------------------------------------------------------------------------
  IF ( MINVAL( flowDirSet(1:8) ) < 1 ) THEN
    WRITE(*,*)'ERROR: init_route_map: minimum of flowDirSet < 1'
    WRITE(*,*)'"Real" flow direction codes must be > 0 for current code.'
    STOP
  ENDIF

! Set error value to be out of range of all other values used.
  dirErr = -10

!-------------------------------------------------------------------------------
! Find any downstream point.
!-------------------------------------------------------------------------------

  DO iy=1,nyRoute
    DO ix=1,nxRoute

      IF ( gridFlowDir(ix,iy) /= flowDirSet(0) ) THEN

!-------------------------------------------------------------------------------
!       This is a land point on routing grid.
!       Note we have assumed that the flow direction field is set (not missing) at
!       every location. The only "no data" value is that of the flow direction
!       coding scheme, we do not have any gaps in this field.
!-------------------------------------------------------------------------------
        gridMask(ix,iy)=.TRUE.

!       Set index.
        gridIndexG(ix,iy) = (iy-1)*nxRoute + ix

!-------------------------------------------------------------------------------
!       Find the next downstream point.
!-------------------------------------------------------------------------------
!       Initialise to error value.
        ixNext = dirErr
        iyNext = dirErr
        IF ( gridFlowDir(ix,iy) == flowDirSet(9) ) THEN
!         No defined flow direction (e.g. a depression or flow into sea).
!         Set "next" values to -9 (i.e. 1 * no flow direction ).
          ixNext = -9
          iyNext = ixNext
        ELSE
!         Loop over 8 possible flow directions.
          DO idir=1,8
            IF ( gridFlowDir(ix,iy) == flowDirSet(idir) ) THEN
               ixNext = ix + flowDirDelta(idir,1)
               iyNext = iy + flowDirDelta(idir,2)
              EXIT
            ENDIF
          ENDDO
        ENDIF

!-------------------------------------------------------------------------------
!       Check that we have a valid flow direction.
!       This would not generally be the case if the file has missing data values
!       with a value other than the coding scheme's missing data value.
!-------------------------------------------------------------------------------
        IF ( ixNext==dirErr .OR. iyNext==dirErr ) THEN
           WRITE(*,*)'ERROR: init_route_map: unexpected value of coded flow direction.'
          WRITE(*,*)'ix,iy=',ix,iy,' gridFlowDir=',gridFlowDir(ix,iy)
          WRITE(*,*)'flowDirType=',TRIM(flowDirType)
          WRITE(*,*)'Acceptable values for flow direction=flowDirSet=',flowDirSet(:)
          WRITE(*,*)'A possible cause is if the missing data value in the file'
          WRITE(*,*)'is not the missing data value used by the direction encoding.'
          WRITE(*,*)'For flowDirType=',TRIM(flowDirType),' missing data value=',flowDirSet(0)
          STOP
        ENDIF

!-------------------------------------------------------------------------------
!       Look for routing past the edges of the grid.
!-------------------------------------------------------------------------------
        IF ( ixNext==0 .OR. ixNext==nxRoute+1 ) THEN

!         Action depends upon type of grid.

          IF ( globalLon ) THEN
!           This is a global grid, so impose cyclic boundary conditions.
            IF ( ixNext == 0 ) THEN
              ixNext = nxRoute
            ELSEIF ( ixNext == nxRoute+1 ) THEN
              ixNext = 1
            ENDIF
          ELSE
!-------------------------------------------------------------------------------
!           If the next downstream point lies off the grid, set indices to
!           -1 * flow direction index. We will not try to access this next point
!           since it is off the edge of the grid, but we would still like to know
!           where it is so that we can calculate the length of the river channel
!           exiting the current point - so that, as far as possible, regional
!           and global implementations will give same answers.
!           NB This works for outflow from grid, but if a regional run is
!           missing an inflow from outside the grid, answers will
!           generally be different in global and regional implementations.
!-------------------------------------------------------------------------------
            DO idir=1,8
              IF ( gridFlowDir(ix,iy) == flowDirSet(idir) ) THEN
                ixNext = -1 * idir
                iyNext = ixNext
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        IF ( iyNext<1 .OR. iyNext>nyRoute ) THEN

!-------------------------------------------------------------------------------
!           If the next downstream point lies off the grid, set indices to
!           -1 * flow direction index. We will not try to access this next point
!           since it is off the edge of the grid, but we would still like to know
!           where it is so that we can calculate the length of the river channel
!           exiting the current point - so that, as far as possible,regional and
!           global implementations will give same answers.
!           NB This works for outflow from grid, but if
!           a regional run is missing an inflow from outside the grid, answers will
!           generally be different in global and regional implementations.
!           For a global grid, this condition means flow over the pole of the grid,
!           which shouldn't be allowed, but we let pass here as we haven't tested
!           if a grid is fully global (only tested in longitude).
!-------------------------------------------------------------------------------
            DO idir=1,8
              IF ( gridFlowDir(ix,iy) == flowDirSet(idir) ) THEN
                iyNext = -1 * idir
                ixNext = iyNext
                EXIT
              ENDIF
            ENDDO

        ENDIF
!-------------------------------------------------------------------------------
!       ixNext and iyNext are now:
!         > 0 at all points with defined flow direction, where the immediately
!             downstream point is also on the grid.
!         -8 <= ixNext (or iyNext) <= -1 at points with defined flow direction,
!             where the immediately downstream point lies off the edge of the
!             grid
!         -9 at points with no defined flow direction

!-------------------------------------------------------------------------------
!       Get index (location on grid) of the next downstream point.
!       Values <0 are used for flow across the edge of the grid, or no defined
!       flow direction.
!-------------------------------------------------------------------------------

        IF ( ixNext>0 .AND. iyNext>0 ) THEN
          gridNextPointG(ix,iy) = (iyNext-1)*nxRoute + ixNext
        ELSE
!         ixNext will be <0.
          gridNextPointG(ix,iy) = ixNext
        ENDIF
!-------------------------------------------------------------------------------

      ENDIF  !  flowDir (land points)
    ENDDO  !  ix
  ENDDO    !  iy
!-------------------------------------------------------------------------------

! Count number of routing points.
  npRoute=COUNT( gridMask(:,:) )

!-------------------------------------------------------------------------------
! Reset gridNextPointG to 0 (no flow direction, or no data) at all points
! other than routing points.
!-------------------------------------------------------------------------------
!  WHERE ( .NOT. gridMask(:,:) ) gridNextPointG(:,:) = 0

!-------------------------------------------------------------------------------
! Check that the downstream point is land. Some input datasets will already
! have enforced this (e.g. TRIP, if correct), but it's worth checking!
! If downstream point is sea, indicate that the current point is an outlet by
! setting next location to indicate no defined flow direction.
!-------------------------------------------------------------------------------
  DO iy=1,nyRoute
    DO ix=1,nxRoute
      IF ( gridMask(ix,iy) ) THEN
!       This is a land point.
!       Get location on grid of next downstream point.
        IF ( gridNextPointG(ix,iy) > 0 ) THEN
!         There is a downstream point, and it is in the routing grid.
          CALL getXYPos( gridNextPointG(ix,iy),nxRoute,nyRoute,ixNext,iyNext )
!         If downstream point is sea, reset flow direction to "no direction".
          IF ( .NOT. gridMask(ixNext,iyNext) )   &
               gridNextPointG(ix,iy) = flowDirSet(9)
        ENDIF
      ENDIF  !  gridMask
    ENDDO
  ENDDO

  END SUBROUTINE init_route_map
!###############################################################################
!###############################################################################
!###############################################################################
