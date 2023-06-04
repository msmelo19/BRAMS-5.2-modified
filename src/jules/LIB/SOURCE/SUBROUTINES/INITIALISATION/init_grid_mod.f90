! module init_grid_mod
! Contains subroutines used to set up the model grid.
!
!###############################################################################
!###############################################################################
  MODULE init_grid_mod

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType

  IMPLICIT NONE

! Local scalar parameters.
  LOGICAL, PARAMETER :: byteSwap = .FALSE.   !  FALSE means that byte order of
!                    data in a binary file are not reversed after reading

  CONTAINS

!###############################################################################
!###############################################################################

! subroutine init_grid
! Internal procedure in module init_grid_mod.
! Read in details of and set up the model grid.
!
! Note that the land fraction (flandg) is always read, even in cases when it is
! not required - e.g. if specifying points via a list of point numbers, all
! points must be land (and 100% land at present), yet land fraction is still
! read in. This could be changed in future (i.e. set flandg to 1 if readPoints).

!###############################################################################

SUBROUTINE Init_Grid(nia,nja,npatch,patch_area,glon,glat)

!-------------------------------------------------------------------------------
  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY : co2_dim_len,co2_dim_row,land_index,land_mask,land_pts  &
                        ,land_pts_trif,n_rows  &
                        ,nice,npft_trif,nx=>row_length,ny=>rows &
                        ,dim_cs1,dim_cs2

  USE coastal, ONLY : &
!  imported arrays with intent(out)
    fland,flandg

  USE file_utils, ONLY :  &
!  imported procedures
    closeFile,fileUnit,findTag,openFile

  USE grid_utils, ONLY :  &
!  imported procedures
     getGridPosLL,getXYpos

  USE inout, ONLY : echo,formatAsc,formatBin,formatLen,formatNc,formatPP  &
          ,mapIn,mapInLand,npoints  &
          ,nxIn,nyIn,jinUnit,tagAscBin,tagNc,yrevIn

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     npft

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2d

  USE switches, ONLY : l_co2_interactive,l_point_data,l_triffid

  USE time_loc, ONLY :  &
!  imported scalars with intent(out)
     nxGrid,nyGrid,regDlat,regDlon,regLat1,regLat1in,regLatLon  &
    ,regLon1,regLon1in,subRegLatLon  &
!  imported arrays with intent(out)
    ,latitude,longitude
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Local scalar parameters.
  CHARACTER(len=formatLen), PARAMETER :: fileFormatPoints = formatAsc  !  format of file for list of points

  INTEGER ::  &!  local SCALARS
    fieldLand    &!  field number in file of land fraction
   ,fieldLat     &!  field number in file of latitude
   ,fieldLon     &!  field number in file of longitude
   ,i            &!  loop counter
   ,ierr         &!  error value
   ,inUnit       &!  unit used to connect to file
   ,ix,iy        &!  loop counters/work
   ,l            &!  loop counter/work
   ,maxDx,maxDy  &!  work
   ,nheaderFieldLand &!  number of headers before each field in flandg file
   ,nheaderFieldLL   &!  number of headers before each field in lat/lon file
   ,nheaderFileLand  &!  number of headers at start of flandg file
   ,nheaderFileLL    &!  number of headers at start of lat/lon file
   ,npointsList  &!  number of points specified via lists
   ,npointsOrig   !  initial value of npoints

  INTEGER ::  &!  local arrays
    xIn(2)  &!  first and last selected columns in input grid
   ,yIn(2)   !  first and last selected rows in the input grid

  REAL ::  &!  local scalars
    dx,dy          &!  work
   ,latMin,lonMin  &!  minimum latitude and longitude values in a grid
   ,rdx,rdy         !  work

  REAL, ALLOCATABLE :: tmpCoord(:,:)  !  Coordinate values

  REAL ::  &!  local arrays
    latIn(2)  &!  min and max latitude of points to be modelled
   ,lonIn(2)  &!  min and max longitude of points to be modelled
   ,xcoord(2) &!  x co-ordinates used to select a sub-area
   ,ycoord(2)  !  y co-ordinates used to select a sub-area

  LOGICAL ::  &!  local SCALARS
    coord           &!  T means points are specified as coordinate pairs
!                       F means points are specified as point number (index)
   ,landOnly        &!  T means that only land points will be modelled
!                       F means all (chosen) points will be modelled
   ,coordLL          &!  T means coordinates are latitude and longitude
!                       F means coordinates are (x,y) indices
   ,llFlag          &!  work
   ,readFileLand    &!  flag indicating if another file is to be read for
!                         land fraction or list of points
   ,readFileLL      &!  flag indicating if another file is to be read for lat/lon data
   ,readFilePoints  &!  flag indicating if another file is to be read for list of points
!                  !  F means a number of points is specified (and all are assumed to be land points)
   ,regLatLonIn   &!  T means grid is specified by coords of one point and interval
!                       Points are filled in standard JULES W-E S-N order.
!                     F means grid is specified by reading lat/lon data for each gridpoint
!                     Note that a regular lat/lon grid can be specified by either method, but an
!                     irregular grid must be specified via regLatLonIn=FALSE.
   ,readPoints    &!  T means points to be modelled are indicated via a list of point numbers
!                     F means maps are read to determine points
   ,subArea       &!  T means that a subarea of the full input grid is used
   ,subAreaLatLon  !  Only used if subArea=TRUE.
!                       T means subArea is specified via lat/lon of rectangle,
!                       F means subArea is specified by x/y indices of rectangle.

  LOGICAL, ALLOCATABLE ::  &!  local allocatable arrays
    usePoint(:,:)    !  used to select points from grid

  CHARACTER(LEN=FORMATLEN) ::  &! local SCALARS
    fileFormatLand  &!  format of file for land fraction
   ,fileFormatLL     !  format of file for lat/lon data

  CHARACTER(LEN=150) ::  &! local SCALARS
    fileNameLand     &!  the name of file for land fraction
   ,fileNameLL       &!  format of file for lat/lon data
   ,fileNamePoints   &!  the name of file with list of points
   ,varNameLand      &!  name of land fraction variable in input file
   ,varNameLat       &!  name of latitude variable in input file
   ,varNameLon        !  name of longitude variable in input file

   INTEGER, INTENT(IN)  :: nia,nja,npatch
   REAL,    INTENT(IN)  :: patch_area(nia,nja,npatch),glon(nia,nja),glat(nia,nja)

!------------------------------------------------------------------------------------------

  if ( echo ) WRITE(*,"(50('-'),/,a)") 'init_grid'

!-------------------------------------------------------------------------------
! Locate the start of this section in input file, establish whether land
! fraction or a list of points is to be read, and establish if a sub-area is to
! be selected.
!-------------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_grid','>INIT_GRID' )
  READ(jinUnit,*) readPoints,coord,coordLL
  READ(jinUnit,*) landOnly
  READ(jinUnit,*) subArea,subAreaLatLon
  READ(jinUnit,*) xcoord(1:2),ycoord(1:2)
  READ(jinUnit,*) npoints

  READ(jinUnit,*) readFilePoints
  READ(jinUnit,*) fileNamePoints

  IF ( readPoints .AND. landOnly ) THEN
    WRITE(*,*)'ERROR: init_grid: readPoints .AND. landOnly'
    WRITE(*,*)'If points are specified by a list, all must be modelled - '
    WRITE(*,*)'cannot then select only land points.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Establish details of land fraction.
!-------------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_grid','>INIT_LAND' )
  READ(jinUnit,*) readFileLand

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.
!-------------------------------------------------------------------------------
  IF ( readFileLand ) THEN

!   An external file will be read.
    READ(jinUnit,*) fileFormatLand
    READ(jinUnit,*) fileNameLand

    SELECT CASE ( fileFormatLand )

      CASE ( formatAsc,formatBin,formatPP )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_grid land',tagAscBin,preInit=.TRUE. )
        READ(jinUnit,*) nheaderFileLand,nheaderFieldLand
        READ(jinUnit,*) fieldLand

      CASE ( formatNc )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_grid land',tagNc,preInit=.TRUE. )
        READ(jinUnit,*) varNameLand

      CASE default
        WRITE(*,*)'ERROR: init_grid: no code for fileFormat=',TRIM(fileFormatLand)
        STOP
    END SELECT

  ELSE   !  NOT readFileLand
!-------------------------------------------------------------------------------
!   Data will be read from run control file.
!   The first field will be read, no headers expected. Field numbers are
!   redundant for stdIn, but will be used to set nfieldFile.
!-------------------------------------------------------------------------------
    fileFormatLand = formatAsc
    fieldLand = 1
    nheaderFieldLand = 0
    nheaderFileLand = 0

  ENDIF   !   readFileLand

!-------------------------------------------------------------------------------
! Get information about how latitude and longitude will be set.
!-------------------------------------------------------------------------------
! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_grid','>INIT_LATLON' )

! Establish how grid is to be specified.
  READ(jinUnit,*) regLatLonIn

! Read details of regular input grid.
  READ(jinUnit,*) regLat1in,regLon1in
  READ(jinUnit,*) regDlat,regDlon

! Ensure we have intervals > 0.
  regDlat = ABS( regDlat )
  regDlon = ABS( regDlon )
  IF ( regDlat<EPSILON(regDlat) .OR. regDlon<EPSILON(regDlon) ) THEN
    WRITE(*,*)'ERROR: init_grid: regDlat or regDlon is too small (effectively zero).'
    STOP
!   Note: regDlat/lon are used for river routing, and also GrADS output.
!   If you really want to run on a very fine grid, the code is probably almost ready
!   and just needs a few details such as these sorted out!
  ENDIF

  IF ( .NOT. regLatLonIn ) THEN
!-------------------------------------------------------------------------------
!   Establish where data will be read from.
!-------------------------------------------------------------------------------
    READ(jinUnit,*) readFileLL

    IF ( readFileLL ) THEN
!-------------------------------------------------------------------------------
!     An external file will be read.
!     Only read parameters for the file format indicated.
!-------------------------------------------------------------------------------

      READ(jinUnit,*) fileFormatLL
      READ(jinUnit,*) fileNameLL

      SELECT CASE ( fileFormatLL )

        CASE ( formatAsc,formatBin,formatPP )
!         Locate the information in run control file.
          CALL findTag( jinUnit,'init_grid latlon',tagAscBin,preInit=.TRUE. )
          READ(jinUnit,*) nheaderFileLL,nheaderFieldLL
          READ(jinUnit,*) fieldLat,fieldLon

        CASE ( formatNc )
!         Locate the information in run control file.
          CALL findTag( jinUnit,'init_grid latlon',tagNc,preInit=.TRUE. )
          READ(jinUnit,*) varNameLat,varNameLon

        CASE default
          WRITE(*,*)'ERROR: init_grid: no code for fileFormat=',TRIM(fileFormatLL)
          STOP
      END SELECT

    ELSE   !  NOT readFileLL

!-------------------------------------------------------------------------------
!     Data will be read from run control file.  The first fields will be
!     read, no headers expected. Field numbers are redundant for stdIn, but will
!     be used to set nfieldFile
!-------------------------------------------------------------------------------
      fileFormatLL = formatAsc
      nheaderFileLL = 0
      nheaderFieldLL = 0
      fieldLat = 1
      fieldLon = 2

    ENDIF   !   readFileLL
  ENDIF   !  regLatLonIn

!###############################################################################
!###############################################################################

!-------------------------------------------------------------------------------
! We now know how all necessary fields are to be read.
! Order of next steps depends upon what information is to be used.
!-------------------------------------------------------------------------------

! Initialise flag.
  llFlag = .FALSE.

! Allocate space for coordinates.
! If coords are not used, allocate minimal space so we can still pass as argument.
  npointsList = npoints
  IF ( .NOT. readPoints ) npointsList = 1
  ALLOCATE( tmpCoord(2,npointsList), stat=ierr )
  IF ( ierr /= 0 ) CALL allocate_error( 'alloc',ierr,'init_grid tmpCoord' )

  IF ( readPoints ) THEN
!-------------------------------------------------------------------------------
!   Set up grid via a list of points.
!-------------------------------------------------------------------------------
    IF ( .NOT. coord ) THEN
      WRITE(*,*)'A list of ',npoints,' point numbers will be read.'
    ELSE
      IF ( coordLL ) THEN
        WRITE(*,*)'A list of ',npoints,' pairs of lat/lon coords will be read.'
      ELSE
        WRITE(*,*)'A list of ',npoints,' pairs of x/y coords wil be read.'
      ENDIF
    ENDIF
    WRITE(*,*)'Setting ny to 1.'

    subArea = .FALSE.
    subAreaLatLon = .FALSE.
    npointsList = npoints
    npointsOrig = npoints
!   Create a vector of points.
    xIn(1) = 1
    nx = npoints
    yIn(1) = 1
    ny = 1

!   Set flag to identify a particular(ly awkward) case.
    IF ( readPoints .AND. coord .AND. coordLL .AND. .NOT.regLatLonIn ) llFlag = .TRUE.

!   Allocate space.
!   If flag set, don't allocate yet. This space is not used until after later
!   allocation.
    IF ( .NOT. llFlag ) CALL allocate_arrays( 'init_grid 1' )

!-------------------------------------------------------------------------------
!   Open file containing list of points.
!-------------------------------------------------------------------------------
    IF ( readFilePoints ) THEN
!     File will be ASCII.
      inUnit = fileUnit( formatAsc )   !  get unit
!     Use first arg to openFile to set recl (for unformatted file) to be
!     enough for a single value.
      CALL openFile( 1,.FALSE.,inUnit,'read',formatAsc,fileNamePoints,'old' )
    ELSE
      WRITE(*,*)'Reading list of points/coords from the run control file.'
      WRITE(*,*)'Data must be the first field(s) encountered.'
      inUnit = jinUnit
!     Locate the start of the data in the run control file.
      CALL findTag( inUnit,'init_grid','>DATA_POINTS',preInit=.TRUE. )
    ENDIF

!-------------------------------------------------------------------------------
!   Read a list of points from ASCII file or run control file.
!-------------------------------------------------------------------------------
    IF ( .NOT. coord ) THEN
!     Read a list of points.
      READ(inUnit,*) mapIn(:)
    ELSE
!     Read coordinate pairs.
!     Note that this will read tmpCoord(1,1), (2,1), (1,2), (2,2), etc.
!     lat/lon pairs should be provided as lat(1),lon(1), lat(2),lon(2),...
!     xy pairs should be provided as x(1),y(1), x(2),y(2),....
      READ(inUnit,*) tmpCoord(:,:)
    ENDIF

!------------------------------------------------------------------------------
! Only close the file if it is not the JULES in file
!------------------------------------------------------------------------------
    IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormatPoints )

!------------------------------------------------------------------------------
!   If coordinates were input, calculate point numbers.
!------------------------------------------------------------------------------
    IF ( coord ) THEN

      IF ( coordLL ) THEN
!-------------------------------------------------------------------------------
!       If grid is regular lat/lon, we get point numbers here. For other grids
!       we need to wait until we have read in latitude and longitude fields
!       later, so we set a flag here.
!       Arguably it would be tidier to ignore the possibility of doing here
!       and always process later.
!-------------------------------------------------------------------------------
        IF ( regLatLonIn ) THEN

!         Locate each point in the grid.
          mapIn(:) = getGridPosLL( nxIn,nyIn,tmpCoord(1,:),tmpCoord(2,:)  &
                           ,lat1=regLat1in,lon1=regLon1in  &
                           ,dlat=regDlat,dlon=regDlon )

        ELSE

!         This case will now be dealt with via subArea.
          subArea = .TRUE.
          subAreaLatLon = .TRUE.

        ENDIF

      ELSE
!-------------------------------------------------------------------------------
!       .NOT. coordLL
!       Calculate index from (x,y) coordinates.
!-------------------------------------------------------------------------------
        DO i=1,nx*ny
          ix = NINT( tmpCoord(1,i) )
          iy = NINT( tmpCoord(2,i) )
          mapIn(i) = (iy-1)*nxIn + ix
        ENDDO

      ENDIF   !  coordLL
    ENDIF  !  coord

!------------------------------------------------------------------------------
!   Next section not done if haven't yet read lat/lon fields.
!------------------------------------------------------------------------------

    IF ( .NOT. llFlag ) THEN

!-------------------------------------------------------------------------------
!     Check for points not found.
!-------------------------------------------------------------------------------
      IF ( coord ) THEN
        DO i=1,nx*ny
          IF ( mapIn(i) < 1 ) THEN
            WRITE(*,*)'ERROR: init_grid: point not found.'
            IF ( coordLL ) THEN
              WRITE(*,*)'Input lat and lon=',tmpCoord(:,i)
              WRITE(*,*) 'Check these correspond exactly to a gridpoint.'
            ELSE
              WRITE(*,*)'Input x and y=',tmpCoord(:,i)
              WRITE(*,*) 'Size of input grid=',nxIn,nyIn
            ENDIF
            STOP
          ENDIF
        ENDDO
      ENDIF  !  coord

!------------------------------------------------------------------------------
!     Recreate list if order of rows is to be reversed.
!------------------------------------------------------------------------------
      IF ( yrevIn ) THEN
        DO i=1,nx*ny
          CALL getXYPos( mapIn(i),nxIn,nyIn,ix,iy )
          mapIn(i) = (nyIn-iy)*nxIn + ix
        ENDDO
      ENDIF

!-------------------------------------------------------------------------------
!     Check that mapping appears reasonable.
!-------------------------------------------------------------------------------
      IF ( MAXVAL(mapIn)>nxIn*nyIn .OR. MINVAL(mapIn)<0 ) THEN
        WRITE(*,*)'Input mapping (mapIn) includes a value that is'
        WRITE(*,*)'out of range. Must be in range 1 to nxIn*nyIn,'
        WRITE(*,*)'Range of mapIn=',MINVAL(mapIn),' to ',MAXVAL(mapIn),' nxIn*nyIn=',nxIn*nyIn
        WRITE(*,*)'init_grid: error in mapIn'
        STOP
      ENDIF
!-------------------------------------------------------------------------------
!     Read and process land fraction.
      CALL init_grid_land( nia,nja,npatch,patch_area, fieldLand,nheaderFieldLand,nheaderFileLand  &
               ,nx,ny,landOnly,readFileLand  &
               ,fileFormatLand,fileNameLand,varNameLand )

!-------------------------------------------------------------------------------
!     Get latitude and longitude of each gridpoint.
      CALL init_grid_latlon( nia,nja,glon,glat, fieldLat,fieldLon,nheaderFieldLL,nheaderFileLL  &
               ,npointsList,nx,ny,latIn,lonIn,llFlag  &
               ,readFileLL,regLatLonIn,subAreaLatLon  &
               ,fileFormatLL,fileNameLL,varNameLat,varNameLon )
!-------------------------------------------------------------------------------

    ENDIF   !  llFlag
  ENDIF  !  readPoints

!###############################################################################

  IF ( .NOT. readPoints .OR. llFlag ) THEN

!###############################################################################

!   Under some circumstances (NOT readPoints AND NOT subAreaLatlon AND NOT
!   landOnly) we already know how many points are to be modelled and which
!   points they are. Otherwise, we need to read entire fields and then inspect
!   these to work out which points to model.

    IF ( .NOT. subArea ) THEN
      subAreaLatLon = .FALSE.
!     Assume we will read and use the full input grid.
      xIn(1)=1; xIn(2)=nxIn; yIn(1)=1; yIn(2)=nyIn
    ELSE
      IF ( .NOT. llFlag ) THEN
!       Check range was specified in increasing order.
!       Don't just swap values, since swapping longitudes gives a differnt domain!
        IF ( xcoord(2)<xcoord(1) .OR. ycoord(2)<ycoord(1) ) THEN
          WRITE(*,*)'ERROR: init_grid: range not in increasing order.'
          WRITE(*,*)'To select a subarea, give coords in increasing order.'
          WRITE(*,*)'e.g. lon=-30 to 30, NOT 30 to -30.'
          STOP
        ENDIF
        IF ( .NOT. subAreaLatLon ) THEN
          xIn(1) = NINT( xcoord(1) )
          xIn(2) = NINT( xcoord(2) )
          yIn(1) = NINT( ycoord(1) )
          yIn(2) = NINT( ycoord(2) )
        ELSE
          latIn(:) = ycoord(:)
          lonIn(:) = xcoord(:)
!         Check range of longitude is <= 360.0.
          IF ( lonIn(2)-lonIn(1)-360.0 > EPSILON(lonIn) ) THEN
            WRITE(*,*)'ERROR: init_grid: requested range of longitude > 360.0.'
            WRITE(*,*)'Requested lonIn=',lonIn(:)
            STOP
          ENDIF
        ENDIF   !  subAreaLatLon
      ENDIF   !  llFlag

      IF ( llFlag .OR. subAreaLatLon ) THEN
!       Assume we will read and use the full input grid.
        xIn(1)=1
        xIn(2)=nxIn
        yIn(1)=1
        yIn(2)=nyIn
      END IF

    ENDIF    !  subArea

    nx = xIn(2) - xIn(1) + 1
    ny = yIn(2) - yIn(1) + 1
    npoints = nx * ny
    npointsOrig = npoints

!   Allocate space.
    CALL allocate_arrays( 'init_grid 1' )

!   Allocate work space.
    ALLOCATE( usePoint(nx,ny), STAT=ierr )
    IF ( ierr /= 0 ) THEN
      WRITE(*,*)'ERROR: init_grid: error allocating usePoint. ierr=',ierr
      STOP
    ENDIF

!-------------------------------------------------------------------------------
!   Create input mapping to select the desired points from input grid.
!-------------------------------------------------------------------------------
    mapIn(:) = 0
    usePoint(:,:) = .TRUE.
    i = 0
    DO iy=yIn(1),yIn(1)+ny-1
      DO ix=xIn(1),xIn(1)+nx-1
        i = i + 1
        IF ( .NOT. yrevIn ) THEN
          mapIn(i) = (iy-1)*nxIn + ix
        ELSE
          mapIn(i) = (nyIn-iy)*nxIn + ix
        ENDIF
      ENDDO
    ENDDO
!-------------------------------------------------------------------------------
!   Read and process land fraction.
    CALL init_grid_land( nia,nja,npatch,patch_area, fieldLand,nheaderFieldLand,nheaderFileLand  &
             ,nx,ny,landOnly,readFileLand  &
             ,fileFormatLand,fileNameLand,varNameLand,usePoint )

!-------------------------------------------------------------------------------
!   Get latitude and longitude of each gridpoint.
    CALL init_grid_latlon( nia,nja,glon,glat ,fieldLat,fieldLon,nheaderFieldLL,nheaderFileLL  &
             ,npointsList,nx,ny,latIn,lonIn,llFlag  &
             ,readFileLL,regLatLonIn,subAreaLatLon  &
             ,fileFormatLL,fileNameLL,varNameLat,varNameLon  &
             ,tmpCoord(1,:),tmpCoord(2,:),usePoint )

!-------------------------------------------------------------------------------
!   Establish if we have excluded points and therefore need to alter the model
!   grid.
!-------------------------------------------------------------------------------
!   Recalculate grid size.
    npoints = COUNT( usePoint(:,:) )

    IF ( npoints < npointsOrig ) CALL init_grid_realloc( npointsOrig,usePoint )

!   Deallocate work space.
    DEALLOCATE( usePoint, STAT=ierr )
    IF ( ierr /= 0 ) WRITE(*,*) &
      'WARNING: init_grid: error deallocating usePoint. stat=',ierr

  ENDIF   !  NOT readPoints OR llFlag
!###############################################################################
!###############################################################################
!###############################################################################

!-------------------------------------------------------------------------------
! Check that values are reasonable.
! These are intended as checks for coding errors, now that code is rather complicated!
!-------------------------------------------------------------------------------
  IF ( npoints < 1 ) THEN
    WRITE(*,*)'ERROR: init_grid: npoints < 1. npoints=',npoints
    STOP
  ENDIF
  IF ( land_pts > npoints ) THEN
    WRITE(*,*)'ERROR: init_grid: land_pts > npoints'
    WRITE(*,*)'land_pts=',land_pts,' npoints=',npoints
    STOP
  ENDIF
  IF ( land_pts == 0 ) THEN
    WRITE(*,*)'ERROR: init_grid: land_pts == 0'
    WRITE(*,*)'This means there are no land points in the model grid.'
    WRITE(*,*)'JULES possibly should be capable of running with land_pts=0,'
    WRITE(*,*)'but at present cannot do so.'
    !DSM <vide Joe> STOP  
  ENDIF

! Check for doublers.
  DO i=1,npoints
    IF ( ANY(mapIn(i+1:npoints) == mapIn(i) ) ) THEN
      WRITE(*,*)'mapIn should not have any repeated values.'
      WRITE(*,*)'Repeated value=',mapIn(i)
      WRITE(*,*)'init_grid: repeated value'
      STOP
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! This chunk of code is largely a repeat of something above, only now done
! for llFlag.
!-------------------------------------------------------------------------------
  IF ( llFlag ) THEN
!   Check for points not found.
    DO i=1,nx*ny
      IF ( mapIn(i) < 1 ) THEN
        WRITE(*,*)'ERROR: init_grid: point not found.'
        WRITE(*,*)'Input lat and lon=',tmpCoord(:,i)
        WRITE(*,*) 'Check these correspond exactly to a gridpoint.'
        STOP
      ENDIF
    ENDDO

!   Recreate list if order of rows is to be reversed.
    IF ( yrevIn ) THEN
      DO i=1,nx*ny
        CALL getXYPos( mapIn(i),nxIn,nyIn,ix,iy )
        mapIn(i) = (nyIn-iy)*nxIn + ix
      ENDDO
    ENDIF
  ENDIF  !  llFlag

!-------------------------------------------------------------------------------
! Check lat/lon values are within range! Currently allowing lon=-180 to 360.
!-------------------------------------------------------------------------------
  IF ( ANY(latitude(:,:)<-90.0) .OR. ANY(latitude(:,:)>90.0) .OR.  &
       ANY(longitude(:,:)<-180.0) .OR. ANY(longitude(:,:)>360.0) ) THEN
    WRITE(*,*)'ERROR: init_grid: latitude or longitude out of range'
    WRITE(*,*)'lat must be in range -90 to 90, lon must be in range -180 to 360.'
    WRITE(*,*)'lat range=',MINVAL(latitude),MAXVAL(latitude)
    WRITE(*,*)'lon range=',MINVAL(longitude),MAXVAL(longitude)
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Establish if the grid is a regular lat/lon grid or a subset of such a grid.
! This can be useful for certain applications (e.g. routing).
! Note that a regular grid can (also) be set up by giving a list of coords in the
! standard order (i.e. with regLatLonIn=F).
! If the input grid is regular (specified via regLatLonIn), the model grid is
! a subset of a regular grid, but we still have to determine whether the model
! grid is itself a regular grid.
! If grid is a subset of a regular grid, work out size of the grid.
!-------------------------------------------------------------------------------
! Initially assume a regular grid.
  regLatLon = .TRUE.
  subRegLatLon = .TRUE.
  nxGrid = nx
  nyGrid = ny
  regLat1 = regLat1in
  regLon1 = regLon1in

! Test if this is a "regular" lat/lon grid - this has to be done in all cases
! except if the input grid was regular and the model grid equals this input grid.
  IF ( .NOT. ( regLatLonIn .AND. nx==nxIn .AND. ny==nyIn) ) THEN
    maxDx = 0
    maxDy = 0
!   Find minima.
    latMin = MINVAL( latitude(:,:) )
    lonMin = MINVAL( longitude(:,:) )
!   To set subRegLatLon, loop over grid, checking that each gridpoint is located
!   on a grid that is consistent with the mimimum values and the grid interval.
!   To set regLatLon, further check that points are in the correct order.
    DO iy=1,ny
      DO ix=1,nx
        dy = ( latitude(ix,iy) - latMin ) / regDlat   !  number of dlat separating points
        rdy = dy - INT(dy)                            !  non-integer part of dy
        dx = ( longitude(ix,iy) - lonMin ) / regDlon
        rdx = dx - INT(dx)
        IF ( ABS(rdy)>EPSILON(rdy) .OR. ABS(rdx)>EPSILON(rdx) ) THEN
!         Not on regular grid.
          regLatLon = .FALSE.
          subRegLatLon = .FALSE.
        ELSEIF ( NINT(dx)/=ix .OR. NINT(dy)/=iy ) THEN
!         Point is on an otherwise regular grid, but not in the the order
!         required for a regular grid.
          regLatLon = .FALSE.
          maxDx = MAX( NINT(dx),maxDx )
          maxDy = MAX( NINT(dy),maxDy )
        ENDIF
      ENDDO
    ENDDO

!   If grid is a subset of a larger grid, get size and location of larger grid.
    IF ( .NOT.regLatLon .AND. subRegLatLon ) THEN
      nxGrid = maxDx + 1
      nyGrid = maxDy + 1
      regLat1 = latMin
      regLon1 = lonMin
    ENDIF

  ENDIF  !  regLatLonIn etc

! For clarity, reset size of regular grid if no such grid has been identified.
  IF ( .NOT. subRegLatLon ) THEN
    nxGrid = -1
    nyGrid = -1
  ENDIF

!-------------------------------------------------------------------------------
! Set up grid information
!-------------------------------------------------------------------------------
  n_rows = ny
  nice = 1

!-------------------------------------------------------------------------------
! Set up user dependent grid information
!-------------------------------------------------------------------------------
  IF ( l_triffid ) THEN
    land_pts_trif = land_pts
    npft_trif = npft
!   Soil carbon dimensions
!   Set dim_cs1=4 for the 4 pools of RothC.
    dim_cs1 = 4
    dim_cs2 = land_pts
  ELSE
    land_pts_trif = 1
    npft_trif = 1
!   Soil carbon dimensions
!   Set dim_cs1=1 to use a single soil C pool.
!   Set dim_cs2=1 to save space (variables with this dimension are not used).
    dim_cs1 = 1
    dim_cs2 = 1
  ENDIF

  IF( l_co2_interactive ) THEN
    co2_dim_len = nx
    co2_dim_row = ny
  ELSE
    co2_dim_len = 1
    co2_dim_row = 1
  ENDIF

!-------------------------------------------------------------------------------
! Allocate space.
!-------------------------------------------------------------------------------
  CALL allocate_arrays( 'init_grid 2' )

!-------------------------------------------------------------------------------
! Set up land index and fland.
!-------------------------------------------------------------------------------
  L = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( land_mask(ix,iy) ) THEN
        l = l + 1
        land_index(l) = (iy-1)*nx + ix
        fland(l) = flandg(ix,iy)
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Create input mapping for land points - use mapIn for this point number.
!-------------------------------------------------------------------------------
  DO l=1,land_pts
    mapInLand(l) = mapIn(land_index(l))
  ENDDO

!-------------------------------------------------------------------------------
! Deallocate memory.
!-------------------------------------------------------------------------------
  IF ( ALLOCATED(tmpCoord) ) THEN
    DEALLOCATE( tmpCoord,stat=ierr )
    IF ( ierr/=0 ) WRITE(*,*)'WARNING: init_grid: error deallocating tmpCoord. ierr=',ierr
  ENDIF

!-------------------------------------------------------------------------------
! If requested, write (to stdout) grid information for each gridbox.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN
    WRITE(*,"(/,a,i6,a3,i6,a3,i6,a)")'Size of input grid=',nxIn,' x ',nyIn  &
              ,' = ',nxIn*nyIn,' points.'
    WRITE(*,"(a,i6,a3,i6)")'Size of model grid=',nx,' x',ny
    WRITE(*,"(a,i6)")'Number of points in model grid=',npoints
    IF ( landOnly ) WRITE(*,"(a)")'Only land points are modelled.'
    WRITE(*,"(a,i6)")'Number of land points=',land_pts
    !DSM <vide Joe> WRITE(*,"(a,(10i7))") 'land_index=',land_index(:)
    WRITE(*,"(a,(10i7))") 'mapIn(:)=',mapIn(:)
    WRITE(*,"(a,(10i7))") 'mapInLand(:)=',mapInLand(:)
    WRITE(*,"(a,(10f8.2))") 'latitude=',latitude(:,:)
    WRITE(*,"(a,(10f8.2))") 'longitude=',longitude(:,:)
    IF ( npoints > 1 ) THEN
      WRITE(*,"(2(a,f8.2))") 'Range of latitude='  &
                ,MINVAL(latitude(:,:)),' to ',MAXVAL(latitude(:,:))
      WRITE(*,"(2(a,f8.2))") 'Range of longitude='  &
                ,MINVAL(longitude(:,:)),' to ',MAXVAL(longitude(:,:))
    ENDIF
    IF ( regLatLon ) THEN
      WRITE(*,*)'Grid is "regular" in lat and lon.'
    ELSEIF ( subRegLatLon ) THEN
      WRITE(*,*)'Grid is a subset of a "regular" lat/lon grid of shape '  &
                ,nxGrid,nyGrid
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! A gentle warning in the case of...
!-------------------------------------------------------------------------------
  IF ( npoints>1 .AND. l_point_data ) THEN
    WRITE(*,"(50('-'))")
    WRITE(*,*)'init_grid: WARNING: l_point_data is selected but there is more than one point.'
    WRITE(*,*)'This may be OK! i.e. if several POINTS (rather than gridboxes) are being modelled by one run.'
    WRITE(*,"(50('-'))")
  ENDIF

END SUBROUTINE Init_Grid

!###############################################################################
!###############################################################################
! subroutine init_grid_land
! Internal procedure in module init_grid_mod.
! Read land fraction and set up land points.

  SUBROUTINE init_grid_land( nia,nja,npatch,patch_area, fieldLand,nheaderFieldLand,nheaderFileLand  &
             ,nx,ny,landOnly,readFileLand  &
             ,fileFormatLand,fileNameLand,varNameLand,usePoint )

  USE ancil_info, ONLY :  &
!   imported scalars with intent(out)
     land_pts  &
!   imported arrays with intent(inout)
    ,land_mask

  USE coastal, ONLY : &
!  imported arrays with intent(out)
    flandg

  USE file_utils, ONLY :  &
!   imported procedures
     closeFile,fileUnit,findTag,openFile  &
!  imported arrays with intent(out)
    ,irecPrev

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatNc,jinUnit  &
!  imported scalars with intent(in)
    ,echo,npoints,nxIn,nyIn  &
!  imported arrays with intent(inout)
    ,mapIn

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2d
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  scalar parameters
    stashLand = 205   !  STASH code for land fraction (correct?). PP code=36.

  INTEGER, INTENT(in) ::  &!  in scalars
    fieldLand   &!
   ,nheaderFieldLand  &!
   ,nheaderFileLand   &!
   ,nx          &!
   ,ny           !

  INTEGER ::  &!  local scalars
    i         &!  work
   ,inUnit    &!  unit used to connect to file
   ,ix        &!  work
   ,iy        &!  work
   ,nfieldFile,nheaderT,nlineField  &!  work
   ,np        &!  work
   ,readT,readZ   &!  work
   ,useIndex       !  index in irecPrev

  LOGICAL, INTENT(in) ::  &!  in scalars
    landOnly      &!  T means
   ,readFileLand   !  T means

  LOGICAL, OPTIONAL, INTENT(inout) ::  &!  optional inout arrays
    usePoint(:,:)  !  T at points that are to be modelled
!                         Will be present if landOnly=T (and in other circumstances).

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    fileFormatLand   &!  file format
   ,fileNameLand     &!  name of file
   ,varnameLand       !  name of netCDF variable for land fraction

   INTEGER, INTENT(IN)  :: nia,nja,npatch
   REAL,    INTENT(IN)  :: patch_area(nia,nja,npatch)

!-------------------------------------------------------------------------------
  if ( echo ) WRITE(*,"(a)") 'init_grid_land'
!-------------------------------------------------------------------------------
! Open file containing grid information.
!-------------------------------------------------------------------------------
  IF ( readFileLand ) THEN
    inUnit = fileUnit( fileFormatLand )   !  get unit
!   Use first arg to openFile to set recl (for unformatted file) to be enough for a single value.
    !DSM CALL openFile( 1,.FALSE.,inUnit,'read',fileFormatLand,fileNameLand,'old','init_grid',ncType )
  ELSE
    WRITE(*,*)'Reading land fraction from the run control file.'
    WRITE(*,*)'Data must be the first field encountered.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_grid_land','>DATA_LAND',preInit=.TRUE. )
  ENDIF

!-------------------------------------------------------------------------------
! Simplifying assumptions regarding input file. Otherwise have to read these in.
!-------------------------------------------------------------------------------
  readT = 1       !   time level to read from file
  readZ = 1       !   'z' level to read from file
  nfieldFile = fieldLand  !  # of fields in file. use required field - OK while readT=1
  nheaderT = 0    !  no headers at top of each time
  nlineField = 0  !  will not attempt to read ASCII line-by-line

! Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
! but need to keep index within bounds.
  useIndex = inUnit
  IF ( fileFormatLand == formatNc ) useIndex = 1

!-------------------------------------------------------------------------------
! Read land fraction.
!-------------------------------------------------------------------------------
  !DSM CALL readVar2d( readT,readZ,fieldLand,stashLand,irecPrev(useIndex),nfieldFile  &
  !DSM                ,nheaderFileLand,nheaderT,nheaderFieldLand  &
  !DSM                ,nlineField,nxIn,nyIn,inUnit,varNameLand  &
  !DSM                ,mapIn(:),(/ (i,i=1,npoints) /),fileFormatLand  &
  !DSM                ,flandg,byteSwap,'init_grid_land','init_grid',ncType )

flandg=1-patch_area(:,:,1)  !DSM

!-------------------------------------------------------------------------------
! Close file, as long as it is not the JULES in file
!-------------------------------------------------------------------------------
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormatLand )

!-------------------------------------------------------------------------------
! Count land points. If only allowing land points, remove any sea points.
! We indicate that points are to be removed via usePoint.
!
! Note: We could remake mapIn here, to avoid reading lat/lon for unwanted sea
! points, but we're not doing that - mapIn is needed in init_grid_latlon, and
! needs to have npoints values as there is currently no code to read a small
! number of values and scatter across grid).
!-------------------------------------------------------------------------------
  LAND_PTS = 0
  LAND_MASK(:,:) = .FALSE.
  np = 0
  DO iy=1,ny
    DO ix=1,nx
      !DSM <vide Joe> IF ( flandg(ix,iy) > 0.0 ) THEN
!       Set land fraction to 1. This is necessary while JULES only deals with
!       land points.
        flandg(ix,iy) = 1.0
        LAND_PTS = LAND_PTS + 1
        land_mask(ix,iy) = .TRUE.
      !DSM <vide Joe> ELSEIF ( landOnly ) THEN
!       Indicate that this sea point will not be modelled.
        !DSM <vide Joe> usePoint(ix,iy) = .FALSE.
        !DSM <vide Joe> np = np + 1
      !DSM <vide Joe> ENDIF
    ENDDO
  ENDDO

  IF ( echo ) THEN
    IF ( landOnly .AND. np>0 ) WRITE(*,*) np,' sea points will be removed.'
    WRITE(*,*)'Number of land points=',land_pts
  ENDIF

  END SUBROUTINE init_grid_land

!###############################################################################
!###############################################################################
!###############################################################################

! subroutine init_grid_latlon
! Internal procedure in module init_grid_mod.
! Read latitude and longitude of each gridpoint.
! If necessary, reject points outside given lat/lon range.

  SUBROUTINE init_grid_latlon( nia,nja,glon,glat ,fieldLat,fieldLon,nheaderFieldLL,nheaderFileLL  &
             ,npointsList,nx,ny,latIn,lonIn,llFlag  &
             ,readFileLL,regLatLonIn,subAreaLatLon  &
             ,fileFormatLL,fileNameLL,varNameLat,varNameLon  &
             ,latValIn,lonValIn,usePoint )

  USE file_utils, ONLY :  &
!   imported procedures
     closeFile,fileUnit,findTag,openFile  &
!  imported arrays with intent(out)
    ,irecPrev

  USE grid_utils, ONLY :  &
!  imported procedures
     getGridPosLL,getXYpos

  USE inout, ONLY : echo,formatAsc,formatBin,formatNc,formatPP,mapIn  &
              ,npoints,nxIn,nyIn,jinUnit,yrevIn

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2d

  USE time_loc, ONLY :  &
!  imported scalars with intent(out)
     regDlat,regDlon,regLat1in,regLon1in  &
!  imported arrays with intent(out)
    ,latitude,longitude
!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar parameters.
  INTEGER, PARAMETER ::  stashLat = -1  !  STASH code for latitude - not known. PP code = 176.
  INTEGER, PARAMETER ::  stashLon = -1  !  STASH code for longitude - not known. PP code = 177.

! Scalar arguments with intent(in).
  INTEGER, INTENT(in) ::  &
    fieldLat        &!
   ,fieldlon        &!
   ,nheaderFieldLL  &!
   ,nheaderFileLL   &!
   ,npointsList     &!  number of points in lists
   ,nx              &!
   ,ny               !

  LOGICAL, INTENT(in) ::  &
    llFlag        &!  TRUE if grid was specified via a list of lat/lon values AND
!                        the input grid is not "regular" in lat/lon
   ,readFileLL    &!
   ,regLatLonIn     &!
   ,subAreaLatLon  !

  CHARACTER(len=*), INTENT(in) ::  &
    fileFormatLL   &!  file format
   ,fileNameLL     &!  name of file
   ,varNameLat,varNamelon      !  name of netCDF variables

! Array arguments with intent(in)
  REAL, INTENT(in) :: latIn(2),lonIn(2)

! Optional array arguments with intent(inout)
  REAL, OPTIONAL, INTENT(in) :: latValIn(:)  !  latitudes of the requested points
  REAL, OPTIONAL, INTENT(in) :: lonValIn(:)  !  longitudes of the requested points
  LOGICAL, OPTIONAL, INTENT(inout) ::  usePoint(:,:)
!                    Mask indicating points that are to be modelled.
!                    Will be present if subAreaLatLon=TRUE (and
!                    in other circumstances).

! Scalar variables.
  INTEGER ::  &
    i         &!  work
   ,inUnit    &!  unit used to connect to file
   ,ip        &!  work
   ,ix        &!  work
   ,ixx       &!  work
   ,iy        &!  work
   ,iyy       &!  work
   ,nfieldFile,nheaderT,nlineField  &!  work
   ,npRemove  &!  counter of number of points excluded
   ,readT,readZ   &!  work
   ,useIndex       !  index in irecPrev

  LOGICAL :: lon360   !  work

! Array variables.
  INTEGER :: gIndex(npointsList)   !   index in grid
  REAL :: lonLimit(2)  !  local version of lonIn
  REAL :: lonVal(npointsList)  !  local version of lonValIn
  LOGICAL :: useP(nx,ny)   !  mask for selected points

  INTEGER, INTENT(IN)  :: nia,nja
  REAL,    INTENT(IN)  :: glon(nia,nja),glat(nia,nja)
!-------------------------------------------------------------------------------

  if ( echo ) WRITE(*,"(a)") 'init_grid_latlon'

  IF ( regLatLonIn ) THEN
!-------------------------------------------------------------------------------
!   Calculate coords for all points on grid.
!-------------------------------------------------------------------------------
    DO iy=1,ny
      DO ix=1,nx
!       Get point number in model grid.
        ip = (iy-1)*nx + ix
!       Get location on input grid.
        CALL getXYPos( mapIn(ip),nxIn,nyIn,ixx,iyy )
!       If input is not in default S-to_N order, change iyy to this default.
        IF ( yrevIn ) iyy = nyIn - iyy + 1
!       Get lat and lon
        latitude(ix,iy) = regLat1in + REAL(iyy-1)*regDlat
        longitude(ix,iy) = regLon1in + REAL(ixx-1)*regDlon
      ENDDO
    ENDDO

  ELSE  !  .NOT. regLatLonIn

!   Read coordinates of each point.

!-------------------------------------------------------------------------------
!   Open file.
!-------------------------------------------------------------------------------
    IF ( readFileLL ) THEN

      SELECT CASE ( fileFormatLL )
        CASE ( formatAsc,formatBin,formatPP )
!         Check that each input field uses different data from file.
          IF ( fieldLat == fieldLon ) THEN
            WRITE(*,*)'ERROR: init_grid_latlon: repeated use of data'
            WRITE(*,*)'Each variable must use different data from file.'
            STOP
          ENDIF
        CASE ( formatNc )
!         Check that variable names are different.
          IF ( varNameLat == varNameLon ) THEN
            WRITE(*,*)'ERROR: init_grid_latlon: repeated variable name.'
            WRITE(*,*)'Each variable must use different data from file.'
            STOP
          ENDIF
      END SELECT

!     Get unit
      inUnit = fileUnit( fileFormatLL )
!     Use first arg to openFile to set recl (for unformatted file) to be enough for a single value.
      !DSM CALL openFile( 1,.FALSE.,inUnit,'read',fileFormatLL,fileNameLL,'old','init_grid',ncType )

    ELSE

      WRITE(*,*)'Reading lat/lon from the run control file.'
      WRITE(*,*)'These must be the first two fields, in order lat then lon.'
      inUnit = jinUnit
!     Locate the start of the data in the run control file.
      CALL findTag( inUnit,'init_grid_latlon','>DATA_LATLON',preInit=.TRUE. )

    ENDIF  !  readFileLL

!-------------------------------------------------------------------------------
!   Simplifying assumptions regarding input file. Otherwise have to read these in.
!-------------------------------------------------------------------------------
    readT      = 1          !   time level to read from file
    readZ      = 1          !   'z' level to read from file
    nfieldFile = MAX(fieldLat,fieldLon)  !  # of fields in file. Use max
!                                           required field - OK while readT=1
    nheaderT   = 0          !  no headers at top of each time
    nlineField = 0          !  will not attempt to read ASCII line-by-line

!   Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!   but need to keep index within bounds.
    useIndex = inUnit
    IF ( fileFormatLL == formatNc ) useIndex = 1

!-------------------------------------------------------------------------------
!   Read data.
!-------------------------------------------------------------------------------
    !DSM CALL readVar2d( readT,readZ,fieldLat,stashLat,irecPrev(UseIndex),nfieldFile  &
    !DSM                ,nheaderFileLL,nheaderT,nheaderFieldLL  &
    !DSM                ,nlineField,nxIn,nyIn,inUnit,varNameLat  &
    !DSM                ,mapIn(:),(/(i,i=1,npoints)/),fileFormatLL  &
    !DSM                ,latitude,byteSwap,'init_grid_latlon','init_grid',ncType )

    latitude=glat  !DSM

    !DSM CALL readVar2d( readT,readZ,fieldLon,stashLon,irecPrev(useIndex),nfieldFile  &
    !DSM                ,nheaderFileLL,nheaderT,nheaderFieldLL  &
    !DSM                ,nlineField,nxIn,nyIn,inUnit,varNameLon  &
    !DSM                ,mapIn(:),(/(i,i=1,npoints)/),fileFormatLL  &
    !DSM                ,longitude,byteSwap,'init_grid_latlon','init_grid',ncType )

    longitude=glon  !DSM

!-------------------------------------------------------------------------------
!   Close file unless it is the JULES in file
!-------------------------------------------------------------------------------
    IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormatLL )

  ENDIF  !  regLatLonIn
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! If requested, subset points according to latitude and longitude.
!-------------------------------------------------------------------------------
  IF ( subAreaLatLon ) THEN

!   Establish whether grid longitudes are in range 0 to 360, or -180 to 180.
!   Note that other, more unusual ranges (eg 360-720), will not be dealt with
!   correctly.
    lon360 = .TRUE.
    IF ( MINVAL(longitude(:,:)) < 0.0 ) lon360 = .FALSE.

    IF ( .NOT. llFlag ) THEN

!     Express requested longitudes in same range as grid.
      DO ix=1,2
        lonLimit(ix) = lonIn(ix)
        IF ( lonLimit(ix)<0.0 .AND. lon360 ) THEN
            lonLimit(ix) = lonLimit(ix) + 360.0
        ELSEIF ( lonLimit(ix)>180.0  .AND. .NOT. lon360 ) THEN
          lonLimit(ix) = lonLimit(ix) - 360.0
        ENDIF
      ENDDO

      npRemove = 0
      DO iy=1,ny
        DO ix=1,nx
          IF ( usePoint(ix,iy) ) THEN
            IF ( ( latitude(ix,iy) < latIn(1) ) .OR.  &
                 ( latitude(ix,iy) > latIn(2) ) .OR.  &
                 ( longitude(ix,iy) < lonLimit(1) ) .OR.  &
                 ( longitude(ix,iy) > lonLimit(2) ) ) THEN
!             This point will not be modelled.
              usePoint(ix,iy) = .FALSE.
              npRemove = npRemove + 1
            ENDIF
          ENDIF  !  usePoint
        ENDDO
      ENDDO
!     Check that we still have some points!
      IF ( npRemove == npoints ) THEN
        WRITE(*,*)'ERROR: init_grid_latlon: no points left!'
        WRITE(*,*)'All points lie outside the selected lat/lon range.'
        WRITE(*,*)'Range of latitude=',latIn(1),' to ',latIn(2)
        WRITE(*,*)'Range of longitude=',lonIn(1),' to ',lonIn(2)
        STOP
      ENDIF

    ELSE

!-------------------------------------------------------------------------------
!     llFlag (i.e. we have a list of requested locations)
!-------------------------------------------------------------------------------

!     Express requested longitudes in same range as grid.
      DO i=1,npointsList
        lonVal(i) = lonValIn(i)
        IF ( lonVal(i)<0.0 .AND. lon360 ) THEN
            lonVal(i) = lonVal(i) + 360.0
        ELSEIF ( lonVal(i)>180.0  .AND. .NOT. lon360 ) THEN
          lonVal(i) = lonVal(i) - 360.0
        ENDIF
      ENDDO

!     Get index of each point.
      gIndex(:) = getGridPosLL( nx,ny,latValIn(:),lonVal(:)  &
                               ,latitude=latitude,longitude=longitude )

!     Initialise mask.
      useP(:,:) = .FALSE.

!     Check for points not found, and reset mask.
      mapIn(:) = 0
      DO i=1,npointsList
        IF ( gIndex(i) < 1 ) THEN
!         The requested location was not a gridpoint.
          WRITE(*,*)'ERROR: init_grid_latlon: selected location is not a gridpoint.'
          WRITE(*,*)'i.e. requested coords do not match any point in lat/lon fields.'
          WRITE(*,*)'Requested latitude=',latValIn(i),' longitude=',lonValIn(i)
          STOP
        ELSE
!         Get location in grid.
          CALL getXYPos( gIndex(i),nx,ny,ix,iy )
!         Check if this point had previously been excluded (for example,
!         if it is sea and only land points requested).
          IF ( .NOT. usePoint(ix,iy) ) THEN
            WRITE(*,*)'ERROR: init_grid_latlon: selected location is not available.'
            WRITE(*,*)'It is a valid gridpoint, but has been excluded previously,'
            WRITE(*,*)'probably because it is a sea point.'
            WRITE(*,*)'Requested latitude=',latValIn(i),' longitude=',lonValIn(i)
            STOP
          ENDIF
!         This location is OK. Set mask.
          useP(ix,iy) = .TRUE.
          mapIn(gIndex(i)) = gIndex(i)
        ENDIF
      ENDDO

!     Copy mask.
      usePoint(:,:) = useP(:,:)

    ENDIF  !  llFlag

  ENDIF  !  subAreaLatLon

!-------------------------------------------------------------------------------
! Optional printing to screen.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN
    IF ( subAreaLatLon .AND. .NOT.llFlag ) THEN
      IF ( npRemove > 0 ) THEN
        WRITE(*,*) npRemove,' points fall outside requested lat/lon range and will be removed.'
      ELSE
        WRITE(*,*)'All points fall within requested lat/lon range.'
      ENDIF
    ENDIF
  ENDIF

  END SUBROUTINE init_grid_latlon
!###############################################################################
!###############################################################################

! subroutine init_grid_realloc
! Internal procedure in module init_grid_mod.
! Change the size of the model grid, reallocating memory as required.

  SUBROUTINE init_grid_realloc( npointsOrig,usePoint )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(inout)
     land_pts,nx=>row_length,ny=>rows  &
!  imported arrays with intent(inout)
    ,land_mask

  USE coastal, ONLY : &
!  imported arrays with intent(inout)
     flandg

  USE inout, ONLY :  &
!  imported scalars with intent(in)
     echo  &
!  imported scalars with intent(inout)
    ,npoints  &
!  imported arrays with intent(inout)
    ,mapIn

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error

  USE time_loc, ONLY :  &
!  imported arrays with intent(inout)
     latitude,longitude
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::  &!  in scalars
    npointsOrig   !  value of npoints at time of allocation

  INTEGER ::  &!  local scalars
    i,ierr           &!  work
   ,ix,ix1,ix2,ixx   &!  work
   ,iy,iy1,iy2,iyy   &!  work
   ,np               &!  work
   ,nxNew            &!  x size of new grid (number of columns)
   ,nyNew             !  y size of new grid (number of rows)

  INTEGER ::  &!  local arrays
    iCopy(npointsOrig)  !  work space

  REAL ::  &!  local arrays
    rCopy(nx,ny)  !  work space

  LOGICAL, INTENT(IN) ::  &!  in arrays
    usePoint(:,:)

  LOGICAL ::  &!  local scalars
    vector   !  T means new grid is a vector (nx>1,ny=1)

  LOGICAL ::  &!  local arrays
    lCopy(nx,ny)  !  work space
!-------------------------------------------------------------------------------

  IF ( echo ) WRITE(*,*)'Points have been excluded. Reshaping grid.'

! Initialise assuming new grid will be a vector.
  vector = .TRUE.

  IF ( ny == 1 ) THEN
!   Grid was already a vector. Now make smaller.
    nxNew = npoints
    nyNew = 1
    IF ( echo ) WRITE(*,*)'New grid will be a smaller vector.'
  ELSE
!   Grid was a rectangle.
!   We either move to a vector, or we can try to use a smaller rectangle.
!   The latter is not really required, but is coded here in case someone
!   creates an application that is finds it easier to use a grid since
!   the spatial relationship of points is easier to determine (although
!   this can still be done on a vector, it can become awkward).
!   Establish bounds for rectangle.
    ix1=nx+1; ix2=0; iy1=ny+1; iy2=0
    DO ix=1,nx
      IF ( ANY( usePoint(ix,:) ) ) THEN
        ix1 = ix
        EXIT
      ENDIF
    ENDDO
    DO ix=nx,ix1,-1
      IF ( ANY( usePoint(ix,:) ) ) THEN
        ix2 = ix
        EXIT
      ENDIF
    ENDDO
    DO iy=1,ny
      IF ( ANY( usePoint(:,iy) ) ) THEN
        iy1 = iy
        EXIT
      ENDIF
    ENDDO
    DO iy=ny,iy1,-1
      IF ( ANY( usePoint(:,iy) ) ) THEN
        iy2 = iy
        EXIT
      ENDIF
    ENDDO

!   Establish if all points within the rectangle are to be modelled.
    nxNew = ix2 - ix1 + 1
    nyNew = iy2 - iy1 + 1
    np = nxNew * nyNew
    IF ( np == npoints ) THEN
!     New grid will be a smaller rectangle.
      vector = .FALSE.
      IF ( echo ) WRITE(*,*)'New grid will be a smaller rectangle.'
    ELSE
!     New grid will be a vector.
      nxNew = npoints
      nyNew = 1
      IF ( echo ) WRITE(*,*)'New grid will be a vector (not rectangle).'
    ENDIF

  ENDIF  !  ny

!-------------------------------------------------------------------------------
! Move to new grid.
!-------------------------------------------------------------------------------
! Take a copy of existing fields, deallocate existing space, allocate new space,
! copy existing data to new space.
! Fields to be processed are: mapIn,land_mask, flandg, latitude, longitude.
!
!-------------------------------------------------------------------------------
! mapIn
  iCopy(:) = mapIn(:)
  DEALLOCATE( mapIn, STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'dealloc',ierr,'init_grid_realloc: mapIn' )
    STOP
  ENDIF
  ALLOCATE( mapIn(npoints), STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_grid_realloc: mapIn' )
    STOP
  ENDIF
  ixx = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( usePoint(ix,iy) ) THEN
        ixx = ixx + 1
        i = (iy-1)*nx + ix
        mapIn(ixx) = iCopy(i)
      ENDIF
    ENDDO
  ENDDO
!-------------------------------------------------------------------------------
! land_mask
  lCopy(:,:) = land_mask(:,:)
  DEALLOCATE( land_mask, STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'dealloc',ierr,'init_grid_realloc: land_mask' )
    STOP
  ENDIF
  ALLOCATE( land_mask(nxNew,nyNew), STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_grid_realloc: land_mask' )
    STOP
  ENDIF
  ixx = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( usePoint(ix,iy) ) THEN
        IF ( vector ) THEN
          ixx = ixx + 1
          land_mask(ixx,1) = lCopy(ix,iy)
        ELSE
          ixx = ix - ix1 + 1
          iyy = iy - iy1 + 1
          land_mask(ixx,iyy) = lCopy(ix,iy)
        ENDIF
      ENDIF  !  usePoint
    ENDDO
  ENDDO
!-------------------------------------------------------------------------------
! flandg
  rCopy(:,:) = flandg(:,:)
  DEALLOCATE( flandg, STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'dealloc',ierr,'init_grid_realloc: flandg' )
    STOP
  ENDIF
  ALLOCATE( flandg(nxNew,nyNew), STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_grid_realloc: flandg' )
    STOP
  ENDIF
  ixx = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( usePoint(ix,iy) ) THEN
        IF ( vector ) THEN
          ixx = ixx + 1
          flandg(ixx,1) = rCopy(ix,iy)
        ELSE
          ixx = ix - ix1 + 1
          iyy = iy - iy1 + 1
          flandg(ixx,iyy) = rCopy(ix,iy)
        ENDIF
      ENDIF  !  usePoint
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! latitude
  rCopy(:,:) = latitude(:,:)
  DEALLOCATE( latitude, STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'dealloc',ierr,'init_grid_realloc: latitude' )
    STOP
  ENDIF
  ALLOCATE( latitude(nxNew,nyNew), STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_grid_realloc: latitude' )
    STOP
  ENDIF
  ixx = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( usePoint(ix,iy) ) THEN
        IF ( vector ) THEN
          ixx = ixx + 1
          latitude(ixx,1) = rCopy(ix,iy)
        ELSE
          ixx = ix - ix1 + 1
          iyy = iy - iy1 + 1
          latitude(ixx,iyy) = rCopy(ix,iy)
        ENDIF
      ENDIF  !  usePoint
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! longitude
  rCopy(:,:) = longitude(:,:)
  DEALLOCATE( longitude, STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'dealloc',ierr,'init_grid_realloc: longitude' )
    STOP
  ENDIF
  ALLOCATE( longitude(nxNew,nyNew), STAT=ierr )
  IF ( ierr /= 0 ) THEN
    CALL allocate_error( 'alloc',ierr,'init_grid_realloc: longitude' )
    STOP
  ENDIF
  ixx = 0
  DO iy=1,ny
    DO ix=1,nx
      IF ( usePoint(ix,iy) ) THEN
        IF ( vector ) THEN
          ixx = ixx + 1
          longitude(ixx,1) = rCopy(ix,iy)
        ELSE
          ixx = ix - ix1 + 1
          iyy = iy - iy1 + 1
          longitude(ixx,iyy) = rCopy(ix,iy)
        ENDIF
      ENDIF  !  usePoint
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Update grid size.
!-------------------------------------------------------------------------------
  nx = nxNew
  ny = nyNew
  land_pts = COUNT( land_mask(:,:) )

  END SUBROUTINE init_grid_realloc

!###############################################################################
!###############################################################################

  END MODULE init_grid_mod

!###############################################################################
!###############################################################################
