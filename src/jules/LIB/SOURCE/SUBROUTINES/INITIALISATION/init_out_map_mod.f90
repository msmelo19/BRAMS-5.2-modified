!###############################################################################
!###############################################################################
! module init_out_map_mod
!
! Module containing routines used to initialise output mappings.
!
! Originally gathered in a module so as to get explicit interface...but no
! longer needed!
!###############################################################################
!###############################################################################

  MODULE init_out_map_mod

  CONTAINS

!###############################################################################
!###############################################################################

! subroutine init_out_map
! Internal procedure in module init_out_map_mod.
! Read details of output mappings and set up mappings.

!-------------------------------------------------------------------------------

  SUBROUTINE init_out_map( iout,callNumber,tmpLoc,tmpSuffix )

!-------------------------------------------------------------------------------
  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     nx=>row_length,ny=>rows

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile,fileUnit,findTag,openFile

  USE inout, ONLY : coord,coordList,coordLL,echo,formatAsc,formatBin  &
                   ,mapIn,mapOut,mapOutCompress  &
                   ,nxIn,outNpWrite,nyIn,outAreaLL,outCompress  &
                   ,outName,outGridDxy,outGridNxy,outGridXY,outLLorder,outName  &
                   ,outRangeX,outRangeY,pointsFlag,pointsOut,pointsOutMax  &
                   ,rgProfile,rpProfile,jinUnit,yrevIn,yrevOut

  USE grid_utils, ONLY :   &
!  imported procedures
     getGridPosLL,getXYpos

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     nxInRoute,nxRoute,nyInRoute,nyRoute,routeDlat,routeDlon  &
    ,routeLat1,routeLon1,routeRegLatLon,yrevInRoute  &
!  imported arrays with intent(in)
    ,mapInRoute

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     regDlat,regDlon,regLat1in,regLon1in,subRegLatLon  &
!  imported arrays with intent(in)
    ,latitude,longitude

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar parameters.
  LOGICAL :: doCheck = .FALSE.  !  Switch controlling whether maps are checked.
!               Use TRUE while developing/debugging.

! Scalar arguments with intent(in):
  INTEGER, INTENT(in) ::  &
    callNumber    &!  identifies where this procedure has been called from
   ,iout           !  number of output profile

! Array arguments with intent(in)
  CHARACTER(len=*), INTENT(in) :: tmpLoc(:,:)  !  Coordinates of chosen points (for RP variables)

  CHARACTER(len=*), INTENT(in) :: tmpSuffix(:) !  name for location (for RP variables)
!-------------------------------------------------------------------------------
! Local scalars:

  INTEGER ::  &
    i,ierr,ierrSum,ip,ix,ixout,ixx,iy,iyout,iyy    &!  work
   ,jp                       &!  loop counter
   ,mapVal                   &!  work
   ,np                       &!  work
   ,npGrid                   &!  number of points in model grid
   ,nxGrid,nyGrid            &!  shape of model grid
   ,nxGridIn,nyGridIn        &!  shape of input grid
   ,npMapIn                  &!  size of input mapping for this grid
   ,unit                     &!  work
   ,x1,x2,y1,y2               !  work

  REAL ::  &
    lat1,latMin,latMax,lon1,lonMax,lonMin,lonVal   !  work

  LOGICAL ::  &
    found          &!  work
   ,lon360         &!  work
   ,readFile       &!  flag indicating if another file is to be read
   ,regLatLonGrid  &!  flag indicating if model grid is 'regular' in lat/lon
!                        (or a subset of such a grid). The points needs not be
!                        in the JULES default W-E S-N order.
   ,yrevInput       !  flag indicating if input grid did not follow JULES default S-N order

  CHARACTER(len=150) ::  &
    fileName   !  the name of a file
!-------------------------------------------------------------------------------
! Local arrays.

  INTEGER, ALLOCATABLE :: mapInUse(:)  !  input mapping for this grid
  REAL, ALLOCATABLE :: gridLat(:,:)  !  latitude
  REAL, ALLOCATABLE :: gridLon(:,:)  !  longitude

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! This routine is called twice for each output profile, to do different tasks.

  IF ( iout==1 .AND. callNumber==2 ) WRITE(*,"(50('-'),/,a)") 'init_out_map call2'

  SELECT CASE ( callNumber )

!-------------------------------------------------------------------------------
    CASE ( 1 )

!     Read as much as possible from run control file.

!     Establish how output is to be specified.
      READ(jinUnit,*) pointsFlag(iout,1:2)
      READ(jinUnit,*) outAreaLL(iout)
      READ(jinUnit,*) outRangeX(iout,1:2),outRangeY(iout,1:2)
      READ(jinUnit,*) outCompress(iout),outLLorder(iout)

!-------------------------------------------------------------------------------
!     Check values are acceptable.
!-------------------------------------------------------------------------------
      IF ( pointsFlag(iout,1)<0 .OR. pointsFlag(iout,1)>2 ) THEN
        WRITE(*,*)'ERROR: init_out_map: pointsFlag(1) out of range.'
        WRITE(*,*)'Acceptable values are 0,1 and 2.'
        WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
        STOP
      ENDIF
      IF ( pointsFlag(iout,2)<0 .OR. pointsFlag(iout,1)>5 ) THEN
        WRITE(*,*)'ERROR: init_out_map: pointsFlag(2) out of range.'
        WRITE(*,*)'Acceptable values are 0 to 5.'
        WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
        STOP
      ENDIF
      IF ( pointsFlag(iout,1) == 1 ) THEN
!       Check range was specified in increasing order.
!       Don't just swap values, since swapping longitudes gives a different domain!
        IF ( outRangeX(iout,2) < outRangeX(iout,1) .OR.  &
             outRangeX(iout,2) < outRangeX(iout,1) ) THEN
          WRITE(*,*)'ERROR: init_out_map: range not in increasing order.'
          WRITE(*,*)'To select a subarea, give coords in increasing order.'
          WRITE(*,*)'e.g. lon=-30 to 30, NOT 30 to -30.'
          STOP
        ENDIF
      ENDIF

!-------------------------------------------------------------------------------
!     Stop if certain options are indicated.
!-------------------------------------------------------------------------------
!     The mapping can only be read if we know how many points are to be output,
!     which can only be guaranteed if pointsFlag(1)=2 (otherwise, we need to know
!     what type of variables are chosen, so that we know grid size).
      IF ( pointsFlag(iout,2)==2 .AND. pointsFlag(iout,1)/=2 ) THEN
        WRITE(*,*)'ERROR: init_out_map: pointsFlag(2)==2 .AND. pointsFlag(1)/=2'
        WRITE(*,*)'pointsFlag(2)=2 can only be used with pointsFlag(1)=2.'
        WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
        STOP
      ENDIF

!     pointsFlag(2)=1 can only be used with pointsFlag(1)=1.
!     i.e. output grid can only be a given subarea if subarea was specified!
      IF ( pointsFlag(iout,2)==1 .AND. pointsFlag(iout,1)/=1 ) THEN
        WRITE(*,*)'ERROR: init_out_map: pointsFlag(2)==1 .AND. pointsFlag(1)/=1'
        WRITE(*,*)'pointsFlag(2)==1 can only be used with pointsFlag(1)=1.'
        WRITE(*,*)'pointsFlag(1)=1 sets up the chosen subarea.'
        WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
        STOP
      ENDIF

!     Only allow a simple specification of output vector If a list of points was
!     read in.
      IF ( pointsFlag(iout,2)==5 .AND. pointsFlag(iout,1)/=2 ) THEN
        WRITE(*,*)'ERROR: init_out_map: pointsFlag(2)==5 .AND. pointsFlag(1)/=2'
        WRITE(*,*)'pointsFlag(2)=5 can only be used with pointsFlag(1)=2.'
        WRITE(*,*)'i.e. a list of points must have been read in.'
        WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
        STOP
      ENDIF

!-------------------------------------------------------------------------------
!     If reading a list of points, establish where to get data from and read.
!-------------------------------------------------------------------------------
      IF ( pointsFlag(iout,1) == 2 ) THEN

        READ(jinUnit,*) readFile
        READ(jinUnit,*) fileName

!       Open an ASCII file, if requested.
        IF ( readFile ) THEN
          unit = fileUnit( formatAsc )  !  get unit
          CALL openFile( 1,.FALSE.,unit,'read',formatAsc,fileName,'old' )
        ELSE
          IF ( echo ) WRITE(*,*)'Reading output mapping info from the run&
                      & control file for profile ',TRIM(outName(iout))
          unit = jinUnit
        ENDIF

!       Note that pointsOut and the maps are read via successive reads, so that all
!       reading of data from external file is  kept together in the code (for clarity -
!       hopefully!).

!       Read number of points to be output.
        READ(unit,*) pointsOut(iout)

!       Read flags describing coordinates.
        READ(unit,*) coord(iout),coordLL(iout)

!       If needed, allocate more space for map.
        IF ( pointsOut(iout) > pointsOutMax ) THEN
!         Allocate more space and copy maps.
          CALL allocate_map( iout )
          pointsOutMax = pointsOut(iout)
        ENDIF

!       Read output mapping, or points to be output.

        IF ( coord(iout) ) THEN
!         Read a list of coordinates of points to be output.
          READ(unit,*) coordList(:,iout,1:pointsOut(iout))
        ELSE
!         Read a list of model points to be output.
          READ(unit,*) mapOut(iout,1:pointsOut(iout),1)
        ENDIF

!       Read list of destinations, if required.
        IF ( pointsFlag(iout,2) == 2 )  &
                READ(unit,*) mapOut(iout,1:pointsOut(iout),2)
!-------------------------------------------------------------------------------
!       Close file if it is not the run control file
!-------------------------------------------------------------------------------
        IF ( unit /= jinUnit ) CALL closeFile( unit,formatAsc )

      ENDIF  !  pointsFlag

!-------------------------------------------------------------------------------
!     Find tag indicating next section of input.
!-------------------------------------------------------------------------------
      CALL findTag( jinUnit,'init_out_map >GRID','>GRID',.TRUE. )

!     Establish shape of (uncompressed) output grid. This is only used
!     with certain options.
      READ(jinUnit,*) outGridNxy(iout,1),outGridNxy(iout,2)

!###############################################################################
!###############################################################################

    CASE ( 2 )

!-------------------------------------------------------------------------------
!     Now that we know what variables are in each profile, we can finish mappings.
!-------------------------------------------------------------------------------

!     Note that the input and model grids may be the same.
!     And for a routing grid variable, the relevant grid is the routing grid (or
!     input routing grid.

!-------------------------------------------------------------------------------
!     Get size and number of points in the model and input grids.
!     Set gridbox interval and location of point (1,1) of output grid to equal
!     values given under regDlat and regDlon for input.
!     These may later be changed.
!     Also load lat/lon variables.
!-------------------------------------------------------------------------------
      nxGrid = nx
      nyGrid = ny
      outGridXY(iout,1) = regLon1in
      outGridXY(iout,2) = regLat1in
      outGridDxy(iout,1) = regDlon
      outGridDxy(iout,2) = regDlat
      nxGridIn = nxIn
      nyGridIn = nyIn
      yrevInput = yrevIn
      npMapIn = SIZE( mapIn )
!     Set flag indicating a regular lat/lon grid.
      regLatLonGrid = subRegLatLon

      IF ( rgProfile(iout) .OR. rpProfile(iout) ) THEN
!       Use the full routing grid, even although routing might be "active" at
!       only some of these points.
        nxGrid = nxRoute
        nyGrid = nyRoute
        outGridXY(iout,1) = routeLon1
        outGridXY(iout,2) = routeLat1
        outGridDxy(iout,1) = routeDlon
        outGridDxy(iout,2) = routeDlat
        nxGridIn = nxInRoute
        nyGridIn = nyInRoute
        yrevInput = yrevInRoute
        npMapIn = SIZE( mapInRoute )
        regLatLonGrid = routeRegLatLon
      ENDIF

      npGrid = nxGrid * nyGrid

!     Allocate space for local lat/lon variables and input map.
      ierrSum = 0
      ALLOCATE( gridLat(nxGrid,nyGrid), gridLon(nxGrid,nyGrid), stat=ierr )
      ierrSum = ierrSum + ierr
      ALLOCATE( mapInUse(npMapIn), stat=ierr )
      ierrSum = ierrSum + ierr
      IF ( ierrSum /= 0 ) THEN
        WRITE(*,*)'ERROR: init_out_map: could not allocate. ierrSum=',ierrSum
        STOP
      ENDIF

!     Load input mapping.
      IF ( rgProfile(iout) .OR. rpProfile(iout) ) THEN
        mapInUse(:) = mapInRoute(:)
      ELSE
        mapInUse(:) = mapIn(:)
      ENDIF

!     Load lat/lon values.
!XX   We don't always need these. Would be better to only allocate and set if required.
      IF ( rgProfile(iout) .OR. rpProfile(iout) ) THEN
        IF ( regLatLonGrid ) THEN
          DO iy=1,nyGrid
            gridLat(:,iy) = routeLat1 + REAL(iy-1)*routeDlat
          ENDDO
          DO ix=1,nxGrid
            gridLon(ix,:) = routeLon1 + REAL(ix-1)*routeDlon
          ENDDO
        ELSE
          WRITE(*,*)'ERROR: init_out_map: NOT regLatLonGrid'
          WRITE(*,*)'We need some code to calculate lat and lon on grid!'
          STOP
        ENDIF
      ELSE
        gridLat(:,:) = latitude(:,:)
        gridLon(:,:) = longitude(:,:)
      ENDIF

!-------------------------------------------------------------------------------
!     Routing variables at selected points (type=RP) are treated quite
!     differently from other types of variables - so deal with separately.
!-------------------------------------------------------------------------------
      IF ( rpProfile(iout) ) THEN

        CALL init_out_map_rp( iout,tmpLoc,tmpSuffix )

      ELSE

!-------------------------------------------------------------------------------
!       Check that any attempts to order output using lat/lon order have a
!       suitable model grid, i.e. regular lat/lon.
!-------------------------------------------------------------------------------
        SELECT CASE ( pointsFlag(iout,2) )
          CASE ( 1, 3 )
            IF ( outLLorder(iout) .AND. .NOT.regLatLonGrid ) THEN
              WRITE(*,*)'ERROR: init_out_map: for output to use the lat/lon order'
              WRITE(*,*)'of points, we must have a regular lat/lon model grid.'
              STOP
            ENDIF
        END SELECT

!-------------------------------------------------------------------------------
!       Check that any attempts to use given lat/lon range to define output grid
!       have a suitable model grid, i.e. regular lat/lon or subset thereof).
!-------------------------------------------------------------------------------
        IF ( pointsFlag(iout,2)==1  .AND. outAreaLL(iout) .AND.  &
             .NOT.regLatLonGrid ) THEN
          WRITE(*,*)'ERROR: init_out_map: in order to specify output grid using'
          WRITE(*,*)'a lat/lon range, we must have a regular lat/lon model grid.'
          STOP
        ENDIF

!-------------------------------------------------------------------------------
!       Set the number of model points to be output, and set the part of mapping
!       that indicates which points these are (mapOut(:,:,1)).
!-------------------------------------------------------------------------------

        SELECT CASE ( pointsFlag(iout,1) )

!-------------------------------------------------------------------------------
!         All points in the grid are to be output.
!-------------------------------------------------------------------------------
          CASE ( 0 )

            pointsOut(iout) = npGrid
!           If needed, allocate more space for map.
            IF ( pointsOut(iout) > pointsOutMax ) THEN
!             Allocate more space and copy maps.
              CALL allocate_map( iout )
              pointsOutMax = pointsOut(iout)
            ENDIF
!           Map is simply a list of all grid point numbers.
            mapOut(iout,1:npGrid,1) = (/ (ip,ip=1,npGrid) /)

!-------------------------------------------------------------------------------
!         All points in a given area are to be output.
!-------------------------------------------------------------------------------
          CASE ( 1 )

!           Identify these points and set map.
            CALL init_out_map_subarea( iout,nxGrid,nxGridIn,nyGrid,nyGridIn  &
                    ,yrevInput  &
                    ,gridLat,gridLon,mapInUse )

!-------------------------------------------------------------------------------
!         Points to be output were listed.
!-------------------------------------------------------------------------------
          CASE ( 2 )

!-------------------------------------------------------------------------------
!           If coordinate pairs were input, convert these to grid point numbers.
!-------------------------------------------------------------------------------
            IF ( coord(iout) ) THEN

              IF ( coordLL(iout) ) THEN

!               Coords are lat,lon.
!               Make coords consistent with the range of the longitude variable.
                lon360 = .TRUE.
                IF ( MINVAL(gridLon(:,:) ) < 0.0 ) lon360 = .FALSE.
                DO ip=1,pointsOut(iout)
                  IF ( coordList(2,iout,ip)<0.0 .AND. lon360 ) THEN
                    coordList(2,iout,ip) = coordList(2,iout,ip) + 360.0
                  ELSEIF ( coordList(2,iout,ip)>180.0  .AND. .NOT. lon360 ) THEN
                    coordList(2,iout,ip) = coordList(2,iout,ip) - 360.0
                  ENDIF
                ENDDO

!               Get grid index.
!               NB Always do by providing lat/lon grids, even if grid is regular,
!               since here regLatLonGrid=TRUE even if the grid is actually just
!               a subset of a regular grid (in which case providing origin and
!               gridbox size is no use).
                mapOut(iout,:,1) = getGridPosLL( nxGrid,nyGrid  &
                         ,coordList(1,iout,1:pointsOut(iout))  &
                         ,coordList(2,iout,1:pointsOut(iout))  &
                         ,latitude=gridLat,longitude=gridLon )
!               Check for points not found.
                DO ip=1,pointsOut(iout)
                  IF ( mapOut(iout,ip,1) < 1 ) THEN
                    WRITE(*,*)'ERROR: init_out_map: requested output point not found on grid.'
                    WRITE(*,*)'Lat=',coordList(1,iout,ip),' lon=',coordList(2,iout,ip)
                    WRITE(*,*)'Error for profile #',iout
                    STOP
                  ENDIF
!                 Deal with yrev input grid.
                  IF ( yrevInput ) THEN
                    CALL getXYPos( mapOut(iout,ip,1),nxGridIn,nyGridIn,ix,iy )
                    mapOut(iout,ip,1) = (nyGridIn-iy)*nxGridIn + ix
                  ENDIF
                ENDDO

              ELSE

!               Coords are x,y. Convert to grid index.
!               Note there's no need to worry about yrevIn, since the input (x,y)
!               already refer to the input grid.
                DO ip=1,pointsOut(iout)
                  mapOut(iout,ip,1) = ( NINT(coordList(2,iout,ip)) - 1 ) * nxGridIn + &
                                        NINT(coordList(1,iout,ip))
                ENDDO

              ENDIF  !  coordLL

            ENDIF  !  coord

!-------------------------------------------------------------------------------
!           The lists give point numbers in the input grid. Convert these to point
!           numbers in model grid, and check that all points requested for output
!           are in the model! Not needed if coords were lat/lon.
!-------------------------------------------------------------------------------
            IF ( .NOT. ( coord(iout) .AND. coordLL(iout) ) ) THEN
              DO ip=1,pointsOut(iout)

!               Check if the output point is in the model grid.
                found = .FALSE.

                DO i=1,npGrid
!                 mapIn gives order as in input file, whereas mapOut uses default
!                 (S to N) order. Make consistent.
                  mapVal = mapInUse(i)
                  IF ( yrevInPut ) THEN
                    CALL getXYPos( mapVal,nxGridIn,nyGridIn,ix,iy )
                    mapVal = (nyGridIn-iy)*nxGridIn + ix
                  ENDIF
                  IF ( mapOut(iout,ip,1) == mapVal ) THEN
!                   Output point is in model grid.
                    found = .TRUE.
!                   Convert point number to number in model grid.
                    mapOut(iout,ip,1) = i
                    EXIT
                  ENDIF
                ENDDO

                IF ( .NOT. found ) THEN
                  WRITE(*,*)'ERROR: init_out_map: requested output point is not in'
                  WRITE(*,*)'in the model grid.'
                  WRITE(*,*)'ip=',ip,' mapOut=',mapOut(iout,ip,1)
                  IF ( mapOut(iout,ip,1) > nxGridIn*nyGridIn ) THEN
                    WRITE(*,*)'In fact, the requested point was not even in the'
                    WRITE(*,*)'input grid, which was of size ',nxGridIn*nyGridIn
                  ENDIF
                  WRITE(*,*)'Problem with output profile=',TRIM(outName(iout))
                  STOP
                ENDIF

              ENDDO  !  ip
            ENDIF      !   .NOT. ( coord(iout) .AND. coordLL(iout) )
!-------------------------------------------------------------------------------
          CASE default

            WRITE(*,*)'ERROR: init_out_map: no code for pointsFlag(1)=',pointsFlag(iout,1)
            WRITE(*,*)'Error for iout=',iout
            STOP

!-------------------------------------------------------------------------------

        END SELECT    !  pointsFlag(1)

!-------------------------------------------------------------------------------
!       In all cases, now check that number of output points is reasonable.
!-------------------------------------------------------------------------------
        IF ( pointsOut(iout)<1 .OR. pointsOut(iout)>npGrid ) THEN
          WRITE(*,*)'ERROR: init_out_map: pointsOut is out of range.'
          WRITE(*,*)'pointsOut=',pointsOut(iout)
          WRITE(*,*)'Acceptable values are 1 to number of points in grid=',npGrid
          WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
          !DSM <vide Joe> STOP
        ENDIF

!-------------------------------------------------------------------------------
!       Don't allow the use of compression (e.g. GrADS pdef) if a single point is
!       to be output and/or the input grid is a single point. Really, what I'd
!       like to avoid is if the output grid is a single point - attempting to
!       compress is then an unnecessary complication - but that is not known until
!       later, so we test for some similar conditions!
!       If you really want to compress in this case, just override here!
!-------------------------------------------------------------------------------
        IF ( outCompress(iout) .AND.  &
            ( pointsOut(iout)==1 .OR. nxGridIn*nyGridIn==1 ) ) THEN
          WRITE(*,*)'WARNING: init_out_map: do not indicate a compressed grid'
          WRITE(*,*)'(e.g. GrADS pdef) for a single point or single point input,'
          WRITE(*,*)'as there''s little advantage and added complications!'
          WRITE(*,*)'WARNING: changing outCompress to FALSE for output profile ' &
                   ,TRIM(outName(iout))
          outCompress(iout) = .FALSE.
        ENDIF

!-------------------------------------------------------------------------------
!      Don't attempt to compress a "simple vector" (pointsFlag(2)=5), as no
!      compression is possible.
!-------------------------------------------------------------------------------
       IF ( outCompress(iout) .AND. pointsFlag(iout,2)==5 )  &
              outCompress(iout) = .FALSE.

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       Set 2nd part of mapping (the destination of each point) and the shape of
!       the output grid. Also set location of output grid.
!-------------------------------------------------------------------------------

        SELECT CASE ( pointsFlag(iout,2) )

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
          CASE ( 0 )

!-------------------------------------------------------------------------------
!           The output grid is the model grid.
!-------------------------------------------------------------------------------
            IF ( echo ) WRITE(*,*)'Profile ',TRIM(outName(iout))  &
                                 ,' output grid is the model grid.'

!-------------------------------------------------------------------------------
!           Set shape of (full, uncompressed) output grid to be that of model grid.
!-------------------------------------------------------------------------------
            outGridNxy(iout,1) = nxGrid
            outGridNxy(iout,2) = nyGrid

!-------------------------------------------------------------------------------
!           Set location of output grid by using lat/lon of a point in model grid.
!-------------------------------------------------------------------------------
!           Note that this may not make much sense (the model grid may not be
!           regular in lat/lon, and the point (1,1) might not lie in the SW corner
!           if the points are not stored in lat/lon order) - but that's OK, as
!           the location is not that important in this case! We proceed by getting
!           the correct answer for a regular lat/lon grid with (1,1) being the
!           point in the SW corner.
            outGridXY(iout,1) = gridLon(1,1)
            outGridXY(iout,2) = gridLat(1,1)

!-------------------------------------------------------------------------------
!           Set size of output variable.
!-------------------------------------------------------------------------------
            IF ( .NOT. outCompress(iout) ) THEN
!             Set number of points written to be number in model grid.
              outNpWrite(iout) = npGrid
            ELSE
!             Set number of points written to be number requested.
              outNpWrite(iout) = pointsOut(iout)
            ENDIF

!-------------------------------------------------------------------------------
!           Set mapping.
!-------------------------------------------------------------------------------
            IF ( .NOT. outCompress(iout) ) THEN
!             Set order for storage in output variable. This is the order of
!             the model grid.
              mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)
            ELSE
!             First save current mapping.
              mapOutCompress(iout,:) = mapOut(iout,:,2)
!             Set order for storage in output variable.
              mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)
            ENDIF

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
          CASE ( 2 )

!-------------------------------------------------------------------------------
!           Output mapping and shape of output grid were read in.
!           Mapping can be changed by compression.
!-------------------------------------------------------------------------------

!           Calculate location of point (1,1) of output grid.
!           We take any point in the output (use the 1st), work out where it lies in the
!           output grid, then knowing the lat/lon of that point we can work out the
!           lat/lon of point (1,1) in the output grid.
!           Work out where the point lies in the model grid. We then know the
!           lat/lon of the point.
            CALL getXYPos( mapOut(iout,1,1),nxGrid,nyGrid,ix,iy )
!           Work out where this point lies in the output grid.
            CALL getXYPos( mapOut(iout,1,2),outGridNxy(iout,1),outGridNxy(iout,2),ixOut,iyOut )
            lat1 = gridLat(ix,iy)-REAL(iyOut-1)*outGridDxy(iout,2)
            lon1 = gridLon(ix,iy)-REAL(ixOut-1)*outGridDxy(iout,1)

            IF ( .NOT. outCompress(iout) ) THEN

!             Not compressing, so mapping is as read in.
!             Set number of points written, including any padding.
              outNpWrite(iout) = outGridNxy(iout,1) *  outGridNxy(iout,2)

            ELSE

!             Compressing, so recalculate mapping.
!             We have a list of arbitrary points in the model grid and locations
!             in output grid. Compression means we only create output for those
!             points, not for any "empty" locations in the output grid.

!             Set number of points written to be number requested.
              outNpWrite(iout) = pointsOut(iout)

!             First save current mapping.
              mapOutCompress(iout,:) = mapOut(iout,:,2)
!             Set order for storage in output variable.
              mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)

            ENDIF    !  outCompress

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

          CASE ( 1, 3 )

!           1 = The output grid is given by the selected subarea.
!           3 = The output grid is the smallest rectangle that contains all requested
!                 points.
!           The specification of the grid location is done separately for these
!           two cases, but the mapping is calculated by common code.

            IF ( pointsFlag(iout,2) == 1 ) THEN
!-------------------------------------------------------------------------------
!             The output grid is the area chosen for output.
!             This was defined by lat/lon, or row/column coords.
!-------------------------------------------------------------------------------
              IF ( outAreaLL(iout) ) THEN
!               Chosen area was defined by lat/lon range.
!               Define grid using location of gridpoints.
!               The model grid is known to be regular in lat/lon (or a subset
!               of such a grid).
!
!               Note that the use of given lat/lon values should be consistent
!               with code in init_out_map_subarea which actually selects the points.
!               Get coords of gridpoints that just lie inside the output grid.
!               We already know we have a regular lat/lon grid (but could be a
!               subset of regular grid).

!               Identify a lat/lon on the output grid - this will be the origin.
!               Use first output point.
!               Get row and column that this point represented in the model grid,
!               so we can get lat/lon.
                CALL getXYPos( mapOut(iout,1,1),nxGrid,nyGrid,ix,iy )
                lon1 = gridLon(ix,iy)
                lat1 = gridLat(ix,iy)

!               Locate the gridpoints in the corners of the output grid.
                ix = FLOOR( ( lon1 - outRangeX(iout,1) ) / outGridDxy(iout,1) )
                lonMin = lon1 - REAL(ix)*outGridDxy(iout,1)
                ix = FLOOR( ( outRangeX(iout,2) - lon1 ) / outGridDxy(iout,1) )
                lonMax = lon1 + REAL(ix)*outGridDxy(iout,1)
                iy = FLOOR( ( lat1 - outRangeY(iout,1) ) / outGridDxy(iout,2) )
                latMin = lat1 - REAL(iy)*outGridDxy(iout,2)
                iy = FLOOR( ( outRangeY(iout,2) - lat1 ) / outGridDxy(iout,2) )
                latMax = lat1 + REAL(iy)*outGridDxy(iout,2)

                outGridXY(iout,1) = lonMin
                outGridXY(iout,2) = latMin
                outGridNxy(iout,1) = NINT( ( lonMax - lonMin ) / outGridDxy(iout,1) ) + 1
                outGridNxy(iout,2) = NINT( ( latMax - latMin ) / outGridDxy(iout,2) ) + 1

              ELSE

!               Chosen area was defined by row/column range in input.
!               Lat/lon description of grid is not critical. Use lat/lon of first output point.
                CALL getXYPos( mapOut(iout,1,1),nxGrid,nyGrid,ix,iy )
                outGridXY(iout,1) = gridLon(ix,iy)
                outGridXY(iout,2) = gridLat(ix,iy)
                outGridNxy(iout,1) = NINT(outRangeX(iout,2)) - NINT(outRangeX(iout,1)) + 1
                outGridNxy(iout,2) = NINT(outRangeY(iout,2)) - NINT(outRangeY(iout,1)) + 1
!               Load minima frm specified range into local variables for use below.
                x1 = NINT( outRangeX(iout,1) )
                y1 = NINT( outRangeY(iout,1) )

              ENDIF  !  outAreaLL

!-------------------------------------------------------------------------------
            ELSE
!-------------------------------------------------------------------------------
!             pointsFlag(2) = 3
!             Find the smallest rectangle that can be drawn around the points.
!             This may be done using lat/lon, or row/column order.
!-------------------------------------------------------------------------------
              IF ( outLLorder(iout) ) THEN
!-------------------------------------------------------------------------------
!               Use lat/lon order.
!-------------------------------------------------------------------------------
!               The "best" (=smallest) output grid might be one that crosses the
!               eastern edge of the domain - i.e. that acknowledges that the Earth
!               is round.  e.g. if the input data starts at 0degE, but we want to
!               model Africa only, rather than the output grid being
!               an equatorial strip around the world, it would be better to take
!               an output grid with edges at the likes of 20W and 60E.

!               Initially, assume a "simple" output grid (that doesn't cross grid
!               edge). Find this by locating min and max lat and lon for points
!               that are to be output.
                latMin = 1.0e6
                latMax = -1.0e6
                lonMin = 1.0e6
                lonMax = -1.0e6
                DO ip=1,pointsOut(iout)
!               Get row and column that this point represented in the model grid.
                  CALL getXYPos( mapOut(iout,ip,1),nxGrid,nyGrid,ix,iy )
                  latMin = MIN( latMin, Gridlat(ix,iy) )
                  latMax = MAX( latMax, Gridlat(ix,iy) )
                  lonMin = MIN( lonMin, GridLon(ix,iy) )
                  lonMax = MAX( lonMax, GridLon(ix,iy) )
                ENDDO

                outGridXY(iout,1) = lonMin
                outGridXY(iout,2) = latMin
                outGridNxy(iout,1) = NINT( ( lonMax - lonMin )  &
                             / outGridDxy(iout,1) ) + 1
                outGridNxy(iout,2) = NINT( ( latMax - latMin )   &
                             / outGridDxy(iout,2) ) + 1

!               If the longitudinal extent of this output grid is >180deg, it may be
!               better to consider an alternative grid. If grid is full, nothing
!               smaller will do!
                IF ( REAL(outGridNxy(iout,1))*regDlon > 180.0        .AND.  &
                   pointsOut(iout) < outGridNxy(iout,1)*outGridNxy(iout,2) ) &
                  CALL init_out_map_grid( iout,nxGrid,nyGrid,gridLon )
!-------------------------------------------------------------------------------
              ELSE
!-------------------------------------------------------------------------------
!               Use x/y order from input.
!-------------------------------------------------------------------------------
!               Locate output points with "smallest" and "largest" x and y locations
!               in input grid.
!               First, initialise.
                x1 = nxGridIn
                x2 = 1
                y1 = nyGridIn
                y2 = 1
                DO ip=1,pointsOut(iout)
!                 Get model grid point number.
                  i = mapOut(iout,ip,1)
!                 Get row and column that this point represented in the input grid.
                  CALL getXYPos( mapInUse(i),nxGridIn,nyGridIn,ix,iy )
!                 If input grid was presented in N to S order (yrevInput), adjust
!                 locations to give location in input grid under default
!                 (S to N ) order.
                  IF ( yrevInput ) iy = nyGridIn - iy + 1
                  x1 = MIN( x1, ix )
                  x2 = MAX( x2, ix )
                  y1 = MIN( y1, iy )
                  y2 = MAX( y2, iy )
                ENDDO
!               Set location of grid using lat/lon of first output point. This may
!               not be particularly sensible, but grid has to be located somewhere
!               (so that GrADS can plot it). Lat/lon of this grid are presumably
!               not essential - at any rate order is being set using x,y from input grid.
                CALL getXYPos( mapOut(iout,1,1),nxGrid,nyGrid,ix,iy )
                outGridXY(iout,1) = gridLon(ix,iy)
                outGridXY(iout,2) = gridLat(ix,iy)
                outGridNxy(iout,1) = x2 - x1 + 1
                outGridNxy(iout,2) = y2 - y1 + 1
              ENDIF

            ENDIF   !  pointsFlag

!-------------------------------------------------------------------------------
!           Set mapping from model to the chosen rectangle.
!-------------------------------------------------------------------------------
!           First assume no compression, later altered if compressing.

            IF ( outLLorder(iout) ) THEN
!-------------------------------------------------------------------------------
!             Use lat/lon to define order.
!-------------------------------------------------------------------------------

!             Make sure the longitude is expressed in same range as requested range.
              lon360 = .TRUE.
              IF ( outGridXY(iout,1) < 0.0 ) lon360 = .FALSE.  !  expect -180 to 180
              DO ip=1,pointsOut(iout)
!               Get location in model grid, so as to get lat/lon.
                CALL getXYPos( mapOut(iout,ip,1),nxGrid,nyGrid,ix,iy )
!               Get longitude in same "convention" as used to define output grid.
                lonVal = gridlon(ix,iy)
                IF ( lonVal<0.0 .AND. lon360 ) THEN
                  lonVal = lonVal + 360.0
                ELSEIF ( lonVal>180.0  .AND. .NOT. lon360 ) THEN
                  lonVal = lonVal - 360.0
                ENDIF
!               Get location in output lat/lon grid.
                ixx = NINT( ( lonVal - outGridXY(iout,1) )  &
                          / outGridDxy(iout,1) ) + 1
                iyy = NINT( ( gridLat(ix,iy) - outGridXY(iout,2) )   &
                          / outGridDxy(iout,2) ) + 1
!               Set mapping.
                mapOut(iout,ip,2) = (iyy-1)*outGridNxy(iout,1) + ixx
              ENDDO

            ELSE
!-------------------------------------------------------------------------------
!             Use position in input grid (minus any offset) to define order.
!-------------------------------------------------------------------------------
              DO ip=1,pointsOut(iout)
!               Get model grid point number.
                i = mapOut(iout,ip,1)
!               Get row and column that this point represented in the input grid.
                CALL getXYPos( mapInUse(i),nxGridIn,nyGridIn,ix,iy )
!               If input grid was presented in N to S order (yrevInput), adjust
!               locations to give location in input grid under default
!               (S to N ) order.
                IF ( yrevInput ) iy = nyGridIn - iy + 1
!               Remove offset.
                ixx = ix - x1 + 1
                iyy = iy - y1 + 1
                mapOut(iout,ip,2) = (iyy-1)*outGridNxy(iout,1) + ixx
              ENDDO

            ENDIF   !  outLLorder

!-------------------------------------------------------------------------------
!           Set number of points written and alter mapping for compression
!-------------------------------------------------------------------------------
            IF ( .NOT. outCompress(iout) ) THEN

!             Set number of points written, including any padding.
              outNpWrite(iout) = outGridNxy(iout,1) *  outGridNxy(iout,2)

            ELSE

!             Set number of points written to be number requested.
              outNpWrite(iout) = pointsOut(iout)

!             First save current mapping.
              mapOutCompress(iout,:) = mapOut(iout,:,2)
!             Set order for storage in output variable.
              mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)

            ENDIF   !  outCompress
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

          CASE ( 4 )

!           The output grid is the input grid.
            IF ( echo ) WRITE(*,*)'Profile ',TRIM(outName(iout)),' output grid is the input grid.'

!-------------------------------------------------------------------------------
!           Set shape of (full, uncompressed) output grid to be that of input grid.
!-------------------------------------------------------------------------------
            outGridNxy(iout,1) = nxGridIn
            outGridNxy(iout,2) = nyGridIn

!-------------------------------------------------------------------------------
!           Set lat/lon location of output grid.
!-------------------------------------------------------------------------------
!           In fact, use values set above.
!           If input grid is regular in lat/lon, this is trivial. If grid was not
!           regular, the values set here are for convenience only (e.g. may want
!           to pretend it is a regular grid for plotting).

!-------------------------------------------------------------------------------
!           Set size of output variable.
!-------------------------------------------------------------------------------
            IF ( .NOT. outCompress(iout) ) THEN
!             Set number of points written to be number in input grid.
              outNpWrite(iout) = nxGridIn * nyGridIn
            ELSE
!             Set number of points written to be number requested.
              outNpWrite(iout) = pointsOut(iout)
            ENDIF   !  outCompress

!-------------------------------------------------------------------------------
!           Set output mapping. This is later revised if outCompress=TRUE.
!-------------------------------------------------------------------------------
            IF ( .NOT. yrevInput ) THEN
!             Reuse input mapping.
!             Note that mapOut(1) is the model gridpoint number of each output point.
!             mapInUse( mapOut(1) ) is thus the input point for this output point.
              DO i=1,pointsOut(iout)
                mapOut(iout,i,2) = mapInUse( mapOut(iout,i,1) )
              ENDDO
            ELSE
!             If input grid used yrevIn, need to alter mapping to get a S-N grid
!             (which is the default output).
!             If yrevIn and yrevOut, there's nothing to do, but since yrevOut reverses
!             the grid below, we actually have to reverse it here so that later this
!             is undone! Sorry about this...
!             Note that yrevIn, used with pointsFlag(2)=4, does NOT result in a N-S order
!             of gridpoints, but in the (default) S-N (GrADS) order. The output
!             is in N-S order ONLY if yrevOut=T.
              DO i=1,pointsOut(iout)
               CALL getXYpos( mapInUse( mapOut(iout,i,1) ),nxGridIn,nyGridIn,ix,iy )
               mapOut(iout,i,2) = (nyIn-iy)*nxGridIn + ix
              ENDDO
            ENDIF

!-------------------------------------------------------------------------------
!           Set mapping if compressing output.
!-------------------------------------------------------------------------------
            IF ( outCompress(iout) ) THEN
!             First save current mapping.
              mapOutCompress(iout,:) = mapOut(iout,:,2)
!             Set order for storage in output variable.
              mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)
            ENDIF

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

          CASE ( 5 )

!           The input points were listed, and will simply be written to a vector,
!           in same order.
            outGridNxy(iout,1) = pointsOut(iout)
            outGridNxy(iout,2) = 1
!           Set location using an arbitrary point.
            outGridXY(iout,1) = gridLat(1,1)
            outGridXY(iout,2) = gridLon(1,1)
            outNpWrite(iout) = pointsOut(iout)
            mapOut(iout,1:pointsOut(iout),2) = (/ (ip,ip=1,pointsOut(iout)) /)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
          CASE default

            WRITE(*,*)'ERROR: init_out_map: no code for pointsFlag(2)=',pointsFlag(iout,2)
            WRITE(*,*)'Error for iout=',iout
            STOP
        END SELECT   !  pointsFlag(2)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!cxyz use some of this checkig for RP too?
!       Check output variable is large enough for suggested mapping.
!-------------------------------------------------------------------------------
        IF ( pointsOut(iout) > outNpWrite(iout) ) THEN
          WRITE(*,*)'ERROR: init_out_map: Output variable too small for profile',TRIM(outName(iout))
          WRITE(*,*)'pointsOut=',pointsOut(iout),'outNpWrite=',outNpWrite(iout)
          STOP
        ENDIF

!-------------------------------------------------------------------------------
!       Redo mapping if output order is to be N to S.
!       Not done if output is compressed, since compression mapping just gives
!       location in output vector for each location in input grid.
!-------------------------------------------------------------------------------
        IF ( yrevOut ) THEN
          IF ( .NOT. outCompress(iout) ) THEN
            DO i=1,pointsOut(iout)
              CALL getXYpos( mapOut(iout,i,2),outGridNxy(iout,1),outGridNxy(iout,2),ix,iy )
              mapOut(iout,i,2) = (outGridNxy(iout,2)-iy)*outGridNxy(iout,1) + ix
            ENDDO
          ENDIF
!         Older code - retained until I'm quite sure it's not needed!
!          ELSE
!!         When compressing, change mapOutCompress.
!            DO i=1,pointsOut(iout)
!              CALL getXYpos( mapOutCompress(iout,i),outGridNxy(iout,1),outGridNxy(iout,2),ix,iy )
!              mapOutCompress(iout,i) = (outGridNxy(iout,2)-iy)*outGridNxy(iout,1) + ix
!            ENDDO
!          ENDIF
        ENDIF

!-------------------------------------------------------------------------------
!       Check that mappings appear reasonable.
!       Values must be in range and not repeated.
!       Check even if maps were created here, rather than read.
!-------------------------------------------------------------------------------

!       Check mapOut(1).
        IF ( ( MINVAL(mapOut(iout,1:pointsOut(iout),1)) < 1 ) .OR.  &
             ( MAXVAL(mapOut(iout,1:pointsOut(iout),1)) > npGrid ) ) THEN
          WRITE(*,*)'ERROR: init_out_map: mapping error for output profile ',TRIM(outName(iout))
          WRITE(*,*)'Error for mapOut(iout,:,1)'
          WRITE(*,*)'Min and max values:',MINVAL(mapOut(iout,1:pointsOut(iout),1))  &
                                         ,MAXVAL(mapOut(iout,1:pointsOut(iout),1))
          WRITE(*,*)'Values must be in range 1 to number of points in model grid=',npGrid
          STOP
        ENDIF

        IF ( ( MINVAL(mapOut(iout,1:pointsOut(iout),2)) < 1 ) .OR.  &
             ( MAXVAL(mapOut(iout,1:pointsOut(iout),2)) > outNpWrite(iout) ) ) THEN
          WRITE(*,*)'ERROR: init_out_map: mapping error for output profile ',TRIM(outName(iout))
          WRITE(*,*)'Error for mapOut(iout,:,2)'
          WRITE(*,*)'Min and max values:',MINVAL(mapOut(iout,1:pointsOut(iout),2))  &
                                         ,MAXVAL(mapOut(iout,1:pointsOut(iout),2))
          WRITE(*,*)'Values must be in range 1 to size of output variable=',outNpWrite(iout)
          STOP
        ENDIF

        IF ( outCompress(iout) ) THEN
          np = outGridNxy(iout,1) * outGridNxy(iout,2)
          IF ( ( MINVAL(mapOutCompress(iout,1:pointsOut(iout))) < 1 ) .OR.  &
               ( MAXVAL(mapOutCompress(iout,1:pointsOut(iout))) > np ) ) THEN
            WRITE(*,*)'ERROR: init_out_map: mapping error for output profile ',TRIM(outName(iout))
            WRITE(*,*)'Error for mapOutCompress.'
            WRITE(*,*)'Min and max values:',MINVAL(mapOutCompress(iout,1:pointsOut(iout)))  &
                                           ,MAXVAL(mapOutCompress(iout,1:pointsOut(iout)))
            WRITE(*,*)'Values must be in range 1 to size of full, uncompressed output grid=',np
            STOP
          ENDIF
        ENDIF

!       Check for repeats.
!       This can be slow for a large grid and unoptimised executable, and is arguably only
!       useful when the mapping might be suspect - i.e when it is read from a file rather
!       than calculated by the code. Also useful to check code during development!

        IF ( doCheck ) THEN

          IF ( echo ) WRITE(*,*)'Checking output mapping for profile #',iout
          DO i=1,2
            DO ip=1,pointsOut(iout)
              DO jp=ip+1,pointsOut(iout)
                IF ( mapOut(iout,jp,i) == mapOut(iout,ip,i) ) THEN
                  WRITE(*,*)'ERROR: init_out_map: repeat in mapOut.'
                  WRITE(*,*)'This map should not have repeated values.'
                  WRITE(*,*)'Output profile ',outName(iout),' i=',i
                  WRITE(*,*)'Repeated value=', mapOut(iout,:,i)
                  WRITE(*,*)'Found at points #',ip,' and ',jp
                  STOP
                ENDIF
              ENDDO  !  jp
!             Older code. Avoid using ANY - it can be slow with large grids!
!             IF ( ANY(mapOut(iout,ip+1:pointsOut(iout),i) == mapOut(iout,ip,i) ) ) THEN
            ENDDO  ! ip
          ENDDO    ! i

!         Check mapOutCompress for repeats.
          IF ( outCompress(iout) ) THEN
            IF ( echo ) WRITE(*,*)'Checking mapOutCompress for profile #',iout
            DO ip=1,pointsOut(iout)
              DO jp=ip+1,pointsOut(iout)
                IF ( mapOutCompress(iout,jp) == mapOutCompress(iout,ip) ) THEN
                  WRITE(*,*)'ERROR: init_out_map: repeat in mapOutCompress.'
                  WRITE(*,*)'This map should not have repeated values.'
                  WRITE(*,*)'Output profile ',outName(iout),' repeated value=',mapOutCompress(iout,ip)
                  WRITE(*,*)'Found at points #',ip,' and ',jp
                  STOP
                ENDIF
              ENDDO  !  jp
!             Older code. Avoid using ANY - it can be slow with large grids!
!              IF ( ANY(mapOutCompress(iout,ip+1:pointsOut(iout)) == mapOutCompress(iout,ip) ) ) THEN
            ENDDO  ! ip
          ENDIF

        ENDIF  !  doCheck

!-------------------------------------------------------------------------------
     ENDIF   !  rpProfile

!-------------------------------------------------------------------------------

!    Deallocate space.
     DEALLOCATE( gridLat,gridLon,mapInUse, stat=ierr )
     IF (ierr /= 0) THEN
       WRITE(*,*) 'ERROR: init_out_map: error on deallocate.'
       STOP
     ENDIF

!-------------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'ERROR: init_out_map: no code for callNumber=',callNumber
      STOP
!-------------------------------------------------------------------------------
  END SELECT   !  callNumber

  END SUBROUTINE init_out_map

!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################

! subroutine allocate_map
! Internal procedure in module init_out_map_mod.
! Used to allocate space for mapping variables. Subsequent calls are used to
! allocate more space if it is found to be needed.

  SUBROUTINE allocate_map( iout )

!-------------------------------------------------------------------------------
  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE inout, ONLY :  &
!   imported scalars with intent(in)
      echo,nout,pointsOutMax  &
!   imported arrays with intent(in)
     ,pointsOut  &
!   imported arrays with intent(inout)
     ,coordList,mapOut,mapOutCompress
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout   !  the number of the current output profile

  INTEGER ::  &!  local SCALARS
    ierr  !  work

  INTEGER :: &!  local ARRAYS
    iwork(nout,pointsOutMax,2)   !  workspace

  REAL :: &!  local ARRAYS
    rwork(2,nout,pointsOutMax)   !  workspace

  LOGICAL ::  &!  local SCALARS
    exists    !  T if space already allocated

!-------------------------------------------------------------------------------

! The operations are repeated for each variable.

!-------------------------------------------------------------------------------

! mapOut

! Find out if space is already allocated (for an earlier profile).
  exists = ALLOCATED( mapOut )

  IF ( exists ) THEN
!   Copy existing space to work space.
    iwork(:,:,:) = mapOut(:,:,:)
    DEALLOCATE( mapOut, stat=ierr )
    IF ( ierr/=0 ) THEN
      WRITE(*,*) 'ERROR: allocate_map: could not deallocate mapOut'
      STOP
    ENDIF
  ENDIF

! Allocate larger space.
  IF ( echo ) WRITE(*,*)'allocate_map: allocating mapOut',pointsOut(iout)
  CALL allocate_arrays( 'init_out_map allocate_mapOut',pointsOut(iout) )

  mapOut(:,:,:) = 0

! Copy any existing data back.
  IF ( exists ) mapOut(:,1:pointsOutMax,:) = iwork(:,1:pointsOutMax,:)
!-------------------------------------------------------------------------------

! mapOutCompress

! Find out if space is already allocated (for an earlier profile).
  exists = ALLOCATED( mapOutCompress )

  IF ( exists ) THEN
!   Copy existing space to work space.
    iwork(:,:,1) = mapOutCompress(:,:)
    DEALLOCATE( mapOutCompress, stat=ierr )
    IF ( ierr/=0 ) THEN
      WRITE(*,*) 'ERROR: allocate_map: could not deallocate mapOutCompress'
      STOP
    ENDIF
  ENDIF

! Allocate larger space.
  IF ( echo ) WRITE(*,*)'allocate_map: allocating mapOutCompress ',pointsOut(iout)
  CALL allocate_arrays( 'init_out_map allocate_mapOutCompress',pointsOut(iout) )

  mapOutCompress(:,:) = 0

! Copy any existing data back.
  IF ( exists ) mapOutCompress(:,1:pointsOutMax) = iwork(:,1:pointsOutMax,1)
!-------------------------------------------------------------------------------

! coordList

! Find out if space is already allocated (for an earlier profile).
  exists = ALLOCATED( coordList )

  IF ( exists ) THEN
!   Copy existing space to work space.
    rwork(:,:,:) = coordList(:,:,:)
    DEALLOCATE( coordList, stat=ierr )
    IF ( ierr/=0 ) THEN
      WRITE(*,*) 'ERROR: allocate_map: could not deallocate coordList'
      STOP
    ENDIF
  ENDIF

! Allocate larger space.
  IF ( echo ) WRITE(*,*)'allocate_map: allocating coordList',pointsOut(iout)
  CALL allocate_arrays( 'init_out_map allocate_coordList',pointsOut(iout) )

  coordList(:,:,:) = 0.0

! Copy any existing data back.
  IF ( exists ) coordList(:,:,1:pointsOutMax) = rwork(:,:,1:pointsOutMax)
!-------------------------------------------------------------------------------

  END SUBROUTINE allocate_map

!###############################################################################
!###############################################################################

! subroutine init_out_map_land
! Internal procedure in module init_out_map_mod.
! Create output mappings for land points.
! This process has been separated from init_out_map so that we can more easily
! allocate space for mapOutLand before we get here.

  SUBROUTINE init_out_map_land( nvarMax,tmpVarType )

  USE ancil_info, ONLY : land_index,land_pts
  USE inout, ONLY : echo,mapOut,mapOutLand,nout,nvarOut  &
                   ,pointsOut,pointsOutLand

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalars with intent(in)
  INTEGER, INTENT(in) ::  &
    nvarMax    !  size of array

! Arrays with intent(in)
  CHARACTER(len=*), INTENT(in) ::  &
    tmpVarType(nvarMax)  !  type of each variable
!-------------------------------------------------------------------------------
! Local scalars.

  INTEGER ::  &
    i,iout,ip,ivar,jvar,n  !  work

  LOGICAL ::  &
    needMap   !  TRUE if mapOutLand needed for a profile

!-------------------------------------------------------------------------------

  WRITE(*,"(50('-'),/,a)") 'init_out_map_land'

  ivar = 0

  DO iout=1,nout

!   The nested loops over pointsOut and land_pts below can be slow if either is
!   large. Avoid unnecessary work by only setting mapOutLand if it is needed for
!   this profile - i.e. if there are variables only on land points.
    needMap = .FALSE.

!   Find out what type of variables have been selected.
!   Information on the chosen variables is still held in "temporary" space.
    DO jvar=1,nvarOut(iout)
      ivar = ivar + 1
!     Reset needMap for any variable type that needs mapOutLand.
!     NB This code must be consistent with subroutines loadOut and writeOut,
!     where mapOutLand is used or not used for each varType.
      SELECT CASE ( tmpVarType(ivar) )
        CASE ( 'LA','PF','SC','SN', 'SO','TI','TY' )
          needMap = .TRUE.
      END SELECT
    ENDDO

!   If mapOutLand is not needed for this profile, move to next profile.
    IF ( .NOT. needMap ) CYCLE

    n = 0
!   This is a slow, nested loop if the grid is large.
!   Note we can't assume anything about the order of the
!   points in mapOut, which cuts off one potential simplification - e.g. they
!   may or may not run in order of increasing land_index.
    IF ( echo ) WRITE(*,*)'Creating mapOutLand for profile #',iout
    DO ip=1,pointsOut(iout)
      DO i=1,land_pts
        IF ( mapOut(iout,ip,1) == land_index(i) ) THEN
!           This is a land point that is to be output.
            n = n + 1
            pointsOutLand(iout) =  pointsOutLand(iout) + 1
            mapOutLand(iout,n,1) = i
            mapOutLand(iout,n,2) = mapOut(iout,ip,2)
!           Note there is no need to take account of yrevOut here - that is
!           already included in mapOut.
          EXIT
        ENDIF
      ENDDO   !  land points
    ENDDO    !   points
  ENDDO    !   iout

  END SUBROUTINE init_out_map_land

!##########################################################################################
!##########################################################################################

! subroutine init_out_map_subarea
! Internal procedure in module init_out_map_mod.
! Identifies points that lie within given area and sets (part of) output mapping.

!-------------------------------------------------------------------------------
  SUBROUTINE init_out_map_subarea( iout,nxGrid,nxGridIn,nyGrid,nyGridIn  &
                            ,yrevInput  &
                            ,latitude,longitude,mapInUse )

!-------------------------------------------------------------------------------
  USE grid_utils, ONLY :   &
!  imported procedures
     getXYpos

  USE inout, ONLY :  &
!  imported arrays with intent(in)
     outAreaLL,outName,outRangeX,outRangeY  &
!  imported arrays with intent(out)
    ,pointsOut  &
!  imported scalars with intent(inout)
    ,pointsOutMax  &
!  imported arrays with intent(inout)
    ,mapOut

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) ::  &
    iout       &!  output profile number
   ,nxGrid     &!  x size of grid
   ,nxGridIn   &!  x size of input grid (as input to model)
   ,nyGrid     &!  y size of input grid
   ,nyGridIn    !  y size of input grid (as input to model)

  LOGICAL, INTENT(in) :: yrevInput   !  T indicates input grid was in N-S order

! Array arguments with intent(in)
  INTEGER, INTENT(in) :: mapInUse(:)  !  input mapping for this grid
  REAL, INTENT(in) ::  &
    latitude(nxGrid,nyGrid)   &!  latitude of each point
   ,longitude(nxGrid,nyGrid)   !  longitude of each point

! Local scalars.
  INTEGER ::  &
    ip   &!  loop counter
   ,ix   &!  loop counter
   ,ixx  &!  work
   ,iy   &!  loop counter
   ,iyy  &!  work
   ,np    !  counter

  REAL ::  &
    lonVal   !  work (a longitude)

  LOGICAL ::  &
    lon360   !  work

! Local arrays.
  LOGICAL ::  &
    outPointMask(nxGrid,nyGrid)  !  mask indicating which points are to be output
!                             TRUE = point is to be output
!-------------------------------------------------------------------------------
  IF ( outAreaLL(iout) ) THEN
!-------------------------------------------------------------------------------
!   Range is given in terms of lat and lon.
!-------------------------------------------------------------------------------

!   Check requested range of longitude is <= 360.0.
    IF ( outRangeX(iout,2)-outRangeX(iout,1)-360.0 > EPSILON(outRangeX) ) THEN
      WRITE(*,*)'ERROR: init_out_map_subarea: requested range of longitude > 360.0.'
      WRITE(*,*)'Requested ',outRangeX(iout,:)
      STOP
    ENDIF

!   Establish whether requested longitudes are in range 0 to 360, or -180 to 180.
    lon360 = .TRUE.
    IF ( MINVAL(outRangeX(iout,:) ) < 0.0 ) lon360 = .FALSE.

!   Start by including all points.
    outPointMask(:,:) = .TRUE.

!   Exclude points outside range.
    DO iy=1,nyGrid
      DO ix=1,nxGrid
!       Make sure the longitude is expressed in same range as requested range.
        lonVal = longitude(ix,iy)
        IF ( lonVal<0.0 .AND. lon360 ) THEN
          lonVal = lonVal + 360.0
        ELSEIF ( lonVal>180.0  .AND. .NOT. lon360 ) THEN
          lonVal = lonVal - 360.0
        ENDIF
        IF ( ( latitude(ix,iy) < outRangeY(iout,1) ) .OR.  &
             ( latitude(ix,iy) > outRangeY(iout,2) ) .OR.  &
             ( lonVal < outRangeX(iout,1) )          .OR.  &
             ( lonVal > outRangeX(iout,2) ) )              &
!            This point will not be modelled.
             outPointMask(ix,iy) = .FALSE.
      ENDDO
    ENDDO

  ELSE

!-------------------------------------------------------------------------------
!   Range is given in terms of row and column numbers in input grid.
!-------------------------------------------------------------------------------
!   Initialise.
    outPointMask(:,:) = .TRUE.
!   Reset mask at points outside range.
!   Loop over points to be output, identifying which relate to input points
!   within the chosen range of x/y coords.
    DO iy=1,nyGrid
      DO ix=1,nxGrid
!       Get point number.
        ip = (iy-1)*nxGrid + ix
!       Get position in input grid of this point.
        CALL getXYPos( mapInUse(ip),nxGridIn,nyGridIn,ixx,iyy )
!       If input grid was presented in N to S order (yrevInput), adjust locations to
!       give location in input grid under default (S to N ) order.
        IF ( yrevInput ) iyy = nyGridIn - iyy + 1
        IF ( ( iyy < outRangeY(iout,1) ) .OR.  &
             ( iyy > outRangeY(iout,2) ) .OR.  &
             ( ixx < outRangeX(iout,1) ) .OR.  &
             ( ixx > outRangeX(iout,2) ) )     &
!            This point will not be modelled.
             outPointMask(ix,iy) = .FALSE.
      ENDDO
    ENDDO

  ENDIF   !  outAreaLL

!-------------------------------------------------------------------------------
! Count points.
!-------------------------------------------------------------------------------
  pointsOut(iout) = COUNT( outPointMask(:,:) )

  IF ( pointsOut(iout) == 0 ) THEN
    WRITE(*,*)'ERROR: init_out_map_subarea: pointsOut=0'
    WRITE(*,*)'There are no points inside the area selected for output.'
    IF ( outAreaLL(iout) ) THEN
      WRITE(*,*)'Selected area lat=',outRangeY(iout,:),' lon=',outRangeX(iout,:)
    ELSE
      WRITE(*,*)'Selected area x=',outRangeX(iout,:),' y=',outRangeY(iout,:)
      IF ( yrevInput ) WRITE(*,*)'NB yrevInput=T, but outRangeY uses model S to N order.'
    ENDIF
    WRITE(*,*)'Error for output profile ',TRIM(outName(iout))
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! If needed, allocate more space for map.
!-------------------------------------------------------------------------------
  IF ( pointsOut(iout) > pointsOutMax ) THEN
!   Allocate more space and copy maps.
    CALL allocate_map( iout )
    pointsOutMax = pointsOut(iout)
  ENDIF

!-------------------------------------------------------------------------------
! Set map. We only set the part of the mapping that identifies the points to
! be output (mapOut(:,:,1), not their destinations.
!-------------------------------------------------------------------------------
  np = 0
  DO iy=1,nyGrid
    DO ix=1,nxGrid
      IF ( outPointMask(ix,iy) ) THEN
        np = np + 1
        mapOut(iout,np,1) = (iy-1)*nxGrid + ix
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE init_out_map_subarea

!###############################################################################
!###############################################################################

! subroutine init_out_map_compress
! Internal procedure in module init_out_map_mod.
! Write files containing the mappings used to "scatter" compressed output
! across a larger grid (e.g. GrADS pdef files). The output is designed to be
! used with GrADS pdef and may well be less useful for other packages.

!-------------------------------------------------------------------------------
  SUBROUTINE init_out_map_compress

!  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
!     nx=>row_length,ny=>rows

  USE file_utils, ONLY :  &
!   imported procedures
     closeFile

  USE inout, ONLY :  &
!  imported scalar parameters
     formatBin,formatNc  &
!  imported scalars with intent(in)
    ,echo,nout,outDir,outFormat,outGrADS,outName,runID,undefOut  &
!  imported arrays with intent(in)
    ,mapOut,mapOutCompress,outgridNxy,pointsOut,outCompress,rgProfile  &
!  imported arrays with intent(out)
    ,compressGridFile,useCompressGrid

  USE grid_utils, ONLY :  &
!   imported procedures
     getXYpos

  USE output_mod, ONLY :  &
!   imported procedures
     newOutFile

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error

  USE readwrite_mod, ONLY :  &
!   imported procedures
     writeVar

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
    dateNext,timeNext

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    ierr,ierrSum,iout,ip,ix,iy,jout,nxPdef,nyPdef,compUnit,varID

  CHARACTER(len=LEN(formatBin)) :: fileFormat  !  work

  INTEGER :: ncCount(1),ncStart(1)

  INTEGER, ALLOCATABLE ::  &!  local arrays
    ival(:,:)  !  the location on the model vector for each point in the grid

  REAL, ALLOCATABLE ::  &!  local arrays
    wt(:,:)   !  weights. =1 at model points, elsewhere =0.

  LOGICAL ::  &!  local scalars
    logWork,same     !  work
!-------------------------------------------------------------------------------

  WRITE(*,"(50('-'),/,a)") 'init_out_map_compress'

! Initialise.
  DO iout=1,nout
    IF ( outCompress(iout) ) useCompressGrid(iout) = iout
  ENDDO

!-------------------------------------------------------------------------------
! If using compressed output, look for identical output grids.
!-------------------------------------------------------------------------------
  DO iout=1,nout-1
    IF ( outCompress(iout) ) THEN
      DO jout=iout+1,nout

!-------------------------------------------------------------------------------
!       Check grid compression, number of points, and type of output are the same.
!-------------------------------------------------------------------------------
        IF ( outCompress(jout) .AND.  pointsOut(iout)==pointsOut(jout) .AND. &
             rgProfile(iout).EQV.rgProfile(jout) ) THEN

!-------------------------------------------------------------------------------
!         Check that mappings are identical.
!         Do not look for mappings that may have the same effect, but are not
!         identical (i.e. different order).
!-------------------------------------------------------------------------------
          same = .TRUE.
          DO ip=1,pointsOut(iout)
            IF ( mapOut(iout,ip,1) /=  mapOut(jout,ip,1) ) same=.FALSE.
            IF ( mapOut(iout,ip,2) /=  mapOut(jout,ip,2) ) same=.FALSE.
            IF ( mapOutCompress(iout,ip) /=  mapOutCompress(jout,ip) ) same=.FALSE.
          ENDDO

          IF ( same ) THEN
!-------------------------------------------------------------------------------
!           We have identical output grids, both using compression.
!           Don't write a compression mapping data file for the latter grid, reuse
!           the earlier file.
!-------------------------------------------------------------------------------
!           Test that we haven't already changed this profile.
            IF ( useCompressGrid(jout) == jout ) THEN
              IF ( echo ) THEN
                WRITE(*,*)'Output grid for profile ',TRIM(outname(jout))  &
                     ,' is identical to that for ',TRIM(outName(iout))
                WRITE(*,*)'Will reuse the compression mapping data written for '  &
                    ,TRIM(outName(iout))
              ENDIF
              useCompressGrid(jout) = iout
            ENDIF
          ENDIF  !  same

        ENDIF  !  outCompress etc
      ENDDO  !  jout
    ENDIF  !  outCompress(iout)
  ENDDO  !  iout

!-------------------------------------------------------------------------------
! Work out how much space is required - enough for the largest grid (the grids being
! the "full" grids, rather than the output vectors).
!-------------------------------------------------------------------------------
  nxPdef = 0
  nyPdef = 0
  DO iout=1,nout
    IF ( outCompress(iout) ) THEN
      nxPdef = MAX( nxPdef, outgridNxy(iout,1) )
      nyPdef = MAX( nyPdef, outgridNxy(iout,2) )
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Allocate space.
!-------------------------------------------------------------------------------
  IF ( ANY( outCompress(:) ) ) THEN
    ALLOCATE( ival(nxPdef,nyPdef),wt(nxPdef,nyPdef), stat=ierr )
    IF ( ierr /= 0 ) THEN
      WRITE(*,*) 'ERROR: init_out_map_compress: could not allocate for pdef variables.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! For each profile, get name of the supplementary file to be used and write file.
!-------------------------------------------------------------------------------
  DO iout=1,nout

    IF ( outCompress(iout) ) THEN

!-------------------------------------------------------------------------------
!     Get name of supplementary file. This may be that from an earlier profile.
!-------------------------------------------------------------------------------
      ip = useCompressGrid(iout)
      compressGridFile(iout) = TRIM(outDir)//'/'//TRIM(runID) // '.' //   &
                              TRIM(outName(ip)) // '.pdefData'
!     Add extension. For now, name binary files with 'gra' (rather than 'bin').
!     Note that even when netCDF output is requested, the supplementary (pdef)
!     data file is a flat binary file (because the test is on outGrADS),
!     so that GrADS can read the netCDF files.
!     If you want the supplementary file to also be netCDF, you'll have to
!     write more code here!
      IF ( outGrADS .OR. outFormat==formatNc ) THEN
        fileFormat = formatBin
        compressGridFile(iout) = TRIM(compressGridFile(iout)) // '.gra'
      ELSE
        fileFormat = outFormat
        compressGridFile(iout) = TRIM(compressGridFile(iout)) //  '.' // TRIM(outFormat)
      ENDIF

      IF ( useCompressGrid(iout)==iout ) THEN

!-------------------------------------------------------------------------------
!       Write a compression mapping data file for this output profile.
!       Also write a corresponding GrADS ctl file.
!-------------------------------------------------------------------------------

        nxPdef = outgridNxy(iout,1)
        nyPdef = outgridNxy(iout,2)

!-------------------------------------------------------------------------------
!       For each point in the output vector, set values at corresponding point in grid.
!-------------------------------------------------------------------------------
!       Initialise.
        ival(:,:) = 0
        wt(:,:) = 0.0

        DO ip=1,pointsOut(iout)

!         Work out where this point lies in the model grid.
          CALL getXYpos( mapOutCompress(iout,ip),nxPdef,nyPdef,ix,iy )

!         Save location in vector.
          ival(ix,iy) = ip
!         Set weight to 1.
          wt(ix,iy) = 1.0

        ENDDO  !  ip

!-------------------------------------------------------------------------------
!       Open data file, and  write a GrADS ctl file.
!       A ctl file is not needed for pdef to work, but it does make it easier to
!       check that the pdef data file is correct.
!-------------------------------------------------------------------------------
        CALL newOutFile( iout,dateNext,timeNext,.FALSE.,.TRUE.,logWork,compUnit )

!-------------------------------------------------------------------------------
!       Set netCDF values - not currently used (except to pass as unused args)
!       as the mapping file is not netCDF.
!-------------------------------------------------------------------------------
        varID = 1
        ncCount(:) = 1
        ncStart(:) = 1

!-------------------------------------------------------------------------------
!       Write fields to data file.
!       If file type is not compatible (losely defined!) with GrADS, only write
!       the first field as later fields are specific to GrADS pdef.
!       Note that the first call to writeVar (for ival) is for an integer!
!-------------------------------------------------------------------------------
        CALL writeVar( 1,compUnit  &
                      ,PACK(ival(1:nxPdef,1:nyPdef),.TRUE.)  &
                      ,fileFormat,'init_out_map_compress' )

        IF ( outGrADS .OR. outFormat==formatNc ) THEN
          CALL writeVar( 2,varID,compUnit,ncCount,ncStart  &
                        ,PACK(wt(1:nxPdef,1:nyPdef),.TRUE.)    &
                        ,fileFormat,'init_out_map_compress' )
!         Create a field of missing data to represent wind rotation values.
          wt(:,:) = undefOut
          CALL writeVar( 3,varID,compUnit,ncCount,ncStart  &
                        ,PACK(wt(1:nxPdef,1:nyPdef),.TRUE.)  &
                        ,fileFormat,'init_out_map_compress' )
        ENDIF

!-------------------------------------------------------------------------------
!       Close data file.
!-------------------------------------------------------------------------------
        CALL closeFile( compUnit,fileFormat )

      ENDIF !  useCompressGrid
    ENDIF   !  outCompress
  ENDDO     !  iout

!-------------------------------------------------------------------------------
! Deallocate space.
!-------------------------------------------------------------------------------
  ierrSum = 0
  IF ( ALLOCATED( ival ) ) THEN
    DEALLOCATE( ival,stat=ierr )
    ierrSum=ierrSum+ierr
  ENDIF
  IF ( ALLOCATED( wt ) ) THEN
    DEALLOCATE( wt,stat=ierr )
    ierrSum=ierrSum+ierr
  ENDIF

  IF ( ierrSum /= 0 ) CALL allocate_error( 'dealloc',ierrSum,'init_out_map_compress' )

  END SUBROUTINE init_out_map_compress

!###############################################################################
!###############################################################################

! subroutine init_out_map_grid
! Internal procedure in module init_out_map_mod.
! Detect if a smaller output grid is possible because of cyclicity in
! longitudinal direction.

  SUBROUTINE init_out_map_grid( iout,nxGrid,nyGrid,lon )

  USE grid_utils, ONLY :   &
!  imported procedures
     getXYpos

  USE inout, ONLY :  &
!  imported scalars with intent(in)
     pointsOut  &
!  imported arrays with intent(in)
    ,mapOut,outGridDxy  &
!  imported arrays with intent(inout)
    ,outGridNxy,outGridXY

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) ::  &
    iout     &!  output profile number
   ,nxGrid   &!  number of columns in model grid
   ,nyGrid    !  number of rows in model grid

! Array arguments with intent(in)
  REAL, INTENT(in) ::  &
    lon(nxGrid,nyGrid)    !  longitude at each point on model grid

!-------------------------------------------------------------------------------
! Local scalars.
  INTEGER ::  &
    ip,ix,iy  &!  work
   ,nrun      &!  length of current run of empty columns
   ,nrunMax   &!  length of longest run of empty columns
   ,nxNew     &!  width (# of gridboxes) of proposed new grid
!   ,x1        &!  location of start of current run of empty columns
!   ,x1Max     &!  start of longest run of empty columns
   ,x2        &!  location of end of current run of empty columns
   ,x2Max      !  end of longest run of empty columns

  LOGICAL ::  &
    inRun    !  flag indicating if currently in a run of empty columns

! Local arrays.
  LOGICAL ::  &
    colMask(outGridNxy(iout,1))  !  mask showing which columns in original
!                                     output grid contain occupied points

!-------------------------------------------------------------------------------
! Create a mask showing which columns in the currently proposed output grid have
! occupied points.
!-------------------------------------------------------------------------------
  colMask(:) = .FALSE.

  DO ip=1,pointsOut(iout)
!   Get location in model grid, so as to get longitude.
    CALL getXYPos( mapOut(iout,ip,1),nxGrid,nyGrid,ix,iy )
!   Get location in output lat/lon grid.
    ix = NINT( ( lon(ix,iy) - outGridXY(iout,1) ) / outGridDxy(iout,1) ) + 1
    colMask(ix) = .TRUE.
  ENDDO

!-------------------------------------------------------------------------------
! If there are no empty columns, no smaller grid is possible, so nothing to do.
!-------------------------------------------------------------------------------
  IF ( ALL(colMask(:)) ) RETURN

!-------------------------------------------------------------------------------
! Look for the largest run of empty columns.
!-------------------------------------------------------------------------------
  nrun = 0
  nrunMax = 0
  inRun = .FALSE.
  DO ix=1,outGridNxy(iout,1)
    IF ( .NOT. colMask(ix) ) THEN
!     Increment length of empty run.
      nrun = nrun + 1
!     If this is start of empty run, save location.
!      IF ( .NOT. inRun ) x1 = ix
!     Set flag showing we are in empty run.
      inRun = .TRUE.
!     Update current estimate of end of run.
      x2 = ix
    ELSE
      IF ( inRun ) THEN
!      We were in a run, but it has just ended.
       inRun = .FALSE.
!      Compare size of run with previous maximum.
       IF ( nrun > nrunMax ) THEN
         nrunMax = nrun
!         x1Max = x1
         x2Max = x2
       ENDIF
!      Reset nrun, ready for start of possible next run.
       nrun = 0
      ENDIF  !  inRun
    ENDIF  !  colMask
  ENDDO

! Deal with end point.
! If end was still in a run, indicate this as last point in run.
  IF ( inRun ) THEN
!   Compare size of run with previous maximum, and update max.
    IF ( nrun > nrunMax ) THEN
      nrunMax = nrun
!      x1Max = x1
      x2Max = x2
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Calculate length of proposed new output grid.
! This is 360deg minus length of longest empty run.
!-------------------------------------------------------------------------------
  nxNew = NINT(360.0/outGridDxy(iout,1)) - nrunMax

!-------------------------------------------------------------------------------
! If proposed grid is smaller than current output grid, adopt proposed grid.
! Don't adopt new grid if it is only slightly smaller - such a change can make
! analysis of output data much harder, for little saving in storage.
!-------------------------------------------------------------------------------
  IF ( nxNew < NINT(0.95*REAL(outGridNxy(iout,1))) ) THEN

!   Set longitude of first point in new grid.
!   First point in new grid is one point after end of longest empty run.
    outGridXY(iout,1) = outGridXY(iout,1) + REAL(x2Max)*outGridDxy(iout,1)

!   Set width of new grid.
    outGridNxy(iout,1) = nxNew

!   Revise representation of start longitude to keep within "usual" range.
    IF ( outGridXY(iout,1) + REAL(outGridNxy(iout,1)-1)*outGridDxy(iout,1)  &
         > 360.0 ) outGridXY(iout,1) = outGridXY(iout,1) - 360.0

  ENDIF

  END SUBROUTINE init_out_map_grid

!###############################################################################
!###############################################################################

! subroutine init_out_map_rp
! Internal procedure in module init_out_map_mod.
! Set output mapping for routing variables at a point (type=RP).

  SUBROUTINE init_out_map_rp( iout,tmpLoc,tmpSuffix )

  USE grid_utils, ONLY :   &
!  imported procedures
     getGridPosLL,getXYpos

  USE inout, ONLY :  &
!  imported arrays with intent(in)
     nvarOut,pointsFlag,varPos  &
!  imported scalars with intent(inout)
    ,pointsOutMax  &
!  imported arrays with intent(inout)
    ,mapOut,outCompress,outGridNxy,outNpWrite,pointsOut,varDesc

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     nxRoute,nyRoute,routeDlat,routeDlon,routeLat1,routeLon1,routeRegLatLon  &
!  imported arrays with intent(in)
    ,routeMask

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar parameters.
  LOGICAL :: llCoord = .TRUE. !  T indicates coords are lat and lon

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) :: iout   !  output profile number

! Array arguments with intent(in)
  CHARACTER(len=*), INTENT(in) :: tmpLoc(:,:)  !  Coordinates of chosen points.

  CHARACTER(len=*), INTENT(in) :: tmpSuffix(:) !  name for location

!-------------------------------------------------------------------------------
! Scalar variables.
  INTEGER :: i,i1,i2,ILEN   !  work
  INTEGER :: ivar,jvar    !  work
  INTEGER :: rx,ry        !  work

  LOGICAL :: activePoint  !  work
  LOGICAL :: onGrid
!  LOGICAL :: route360     !  T means longitudes of routing grid are given in range 0 to 360 deg
!                            F means in range -180 to 180

  CHARACTER(len=7) :: cformat     !  format for a read

! Array variables.
  INTEGER :: gindex(SIZE(tmpSuffix))  !  index on routing grid of each chosen point
  REAL :: tmpLat(SIZE(tmpSuffix)),tmpLon(SIZE(tmpSuffix))  !  Coordinates of chosen points.

!-------------------------------------------------------------------------------

! Each variable is given space in mapOut (unlike the case for other types
! of variables, for which all variables in a profile use all values of mapOut).
! mapOutLand and mapOutCompress are not used.

!-------------------------------------------------------------------------------
! Insist that this profile has pointsFlag(1)=0 - just to ensure that
! mappings were not read.
!-------------------------------------------------------------------------------
  IF ( pointsFlag(iout,1) /= 0 ) THEN
    WRITE(*,*)'ERROR: init_out_map_rp: rpProfile but pointsFlag(1)/=0'
    WRITE(*,*)'Routing diagnostics at individual points can only'
    WRITE(*,*)'be specified if pointsFlag(1)=0.'
    WRITE(*,*)'So...change pointsFlag(1) to 0 for profile #',iout
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Override values read from file.
!-------------------------------------------------------------------------------
  outCompress(iout) = .FALSE.  !  no compression of output (it's only 1 point)
  outGridNxy(iout,:) = 1       !  output grid = 1x1
  outNpWrite(iout) = 1         !  number of points to be written = 1

!-------------------------------------------------------------------------------
! Set "number of points" in this profile to be the number of variables,
! so that we can use mapOut.
!-------------------------------------------------------------------------------
  pointsOut(iout) = nvarOut(iout)

!-------------------------------------------------------------------------------
! If needed, allocate more space for map.
!-------------------------------------------------------------------------------
  IF ( pointsOut(iout) > pointsOutMax ) THEN
!   Allocate more space and copy maps.
    CALL allocate_map( iout )
    pointsOutMax = pointsOut(iout)
  ENDIF

!-------------------------------------------------------------------------------
! Establish what range of longitude is used for routing grid: -180 to 180, or 0 to 360.
!-------------------------------------------------------------------------------
!  route360 = .TRUE.
!  IF ( routeLon1 < 0.0 ) route360=.FALSE.

!-------------------------------------------------------------------------------
! Move location from character to real variables.
!-------------------------------------------------------------------------------
! Loop over all variables in this profile (all are of type RP).
  DO jvar=1,nvarOut(iout)
    ivar = varPos(iout,jvar)
    IF ( llcoord ) THEN
!     Character holds decimal lat/lon, followed by "N" or "E".
      DO i=1,2
        ILEN = LEN_TRIM( tmpLoc(ivar,i) ) - 1   !   omit last character (N or E)
        WRITE(cformat,"( '(f' , i2, '.0)' )") ILEN
        IF ( i == 1 ) THEN
          READ(tmpLoc(ivar,i),cformat) tmpLat(jvar)
        ELSE
          READ(tmpLoc(ivar,i),cformat) tmpLon(jvar)
        ENDIF
      ENDDO
    ELSE
      WRITE(*,*) 'init_out_map_rp: No code for NOT llcoord.'
      STOP
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Convert lat/lon to grid index.
!-------------------------------------------------------------------------------
  IF ( llcoord ) THEN
    IF ( routeRegLatLon ) THEN
      gindex(:) = getGridPosLL( nxRoute,nyRoute,tmpLat,tmpLon  &
                       ,lat1=routeLat1,lon1=routeLon1  &
                       ,dlat=routeDlat,dlon=routeDlon )
    ELSE
      WRITE(*,*) 'init_out_map_rp: No code for NOT lrouteRegLatLon'
      STOP
    ENDIF
  ELSE
    WRITE(*,*) 'init_out_map_rp: No code for NOT llcoord.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Finish processing.
!-------------------------------------------------------------------------------
  DO jvar=1,nvarOut(iout)

    ivar = varPos(iout,jvar)

!-------------------------------------------------------------------------------
!   Get location of point on routing grid.
!-------------------------------------------------------------------------------
    activePoint = .FALSE.
    onGrid = .FALSE.
    IF ( gindex(jvar) > 0 ) THEN
      onGrid = .TRUE.
!     Convert index to (x,y) pos.
      CALL getXYPos( gindex(jvar),nxRoute,nyRoute,rx,ry )
!     Check that this is an "active" routing point, not just on grid.
      IF ( routeMask(rx,ry) ) activePoint = .TRUE.
    ENDIF

!-------------------------------------------------------------------------------
!   Report any error.
!-------------------------------------------------------------------------------
    IF ( .NOT. onGrid .OR. .NOT. activePoint ) THEN
      WRITE(*,*)'ERROR: init_out_map_rp: chosen location is not a valid location on routing grid.'
      WRITE(*,*)'Profile #',iout,' variable #',ivar,' location='  &
             ,TRIM(tmpLoc(ivar,1)),' ',TRIM(tmpLoc(ivar,2))
      IF ( .NOT. onGrid ) THEN
        WRITE(*,*)'Location is not a gridpoint.'
      ELSE
        WRITE(*,*)'Location is a gridpoint, but it is not one that is "active"'
        WRITE(*,*)'e.g. it is sea.'
      ENDIF
      STOP
    ENDIF   !  not valid location

!-------------------------------------------------------------------------------
!   Set maps.
!-------------------------------------------------------------------------------

!   Set mapOut(1) to identify this location.
    mapOut(iout,jvar,1) = (ry-1)*nxRoute + rx

!   Set mapOut(2) = 1, indicating that first location in output grid is to be used.
    mapOut(iout,jvar,2) = 1

!-------------------------------------------------------------------------------
!   Add extra information to the description of this variable - the location and
!   a name/descriptiopn.
!-------------------------------------------------------------------------------
    ILEN = LEN_TRIM( varDesc(ivar) )
!   Make sure we get the extra info in - even at the expense of the description.
!   Work out how much space we need.
    i1 = 1 + LEN_TRIM( tmpLoc(ivar,1) )      !  with 1 preceeding space
    i1 = i1 + 1 + LEN_TRIM( tmpLoc(ivar,2) )
    i1 = i1 + 1 + LEN_TRIM( tmpSuffix(ivar) )
!   Work out where to start adding the extra info.
    i2 = ILEN + 1
!   Adjust if off end of variable.
    IF ( i2+i1 > LEN(varDesc) ) i2=LEN(varDesc)-i1+1
!   Error if not enough space - don't allow into first 10 characters (arbitrary!).
    IF ( i2 < 11 ) THEN
      WRITE(*,*)'ERROR: init_out_map_rp: description too long.'
      WRITE(*,*)'Character variable not long enough.'
      WRITE(*,*)'Reduce amount to be written, or increase length of varDesc.'
      WRITE(*,*)'Error raised for tmpLoc=',tmpLoc(ivar,:),' tmpSuffix=',TRIM(tmpSuffix(ivar))
      STOP
    ENDIF
    varDesc(ivar) = TRIM(varDesc(ivar)(1:i2-1)) // ' ' // TRIM(tmpLoc(ivar,1))  &
                   // ' ' // TRIM(tmpLoc(ivar,2)) // ' ' // TRIM(tmpSuffix(ivar))

  ENDDO   !  jvar

  END SUBROUTINE init_out_map_rp

!###############################################################################
!###############################################################################

 END MODULE init_out_map_mod
!###############################################################################
!###############################################################################
