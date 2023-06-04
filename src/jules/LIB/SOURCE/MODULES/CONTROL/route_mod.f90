! module route_mod
! Contains variables used by channel routing schemes.
! At present holds both "physics" and control variables - which will likely have
! to be separated before the UM arrives...

  MODULE route_mod

!-------------------------------------------------------------------------------
  USE inout, ONLY :  &
!  imported scalar parameters
     optLen
!-------------------------------------------------------------------------------

! Scalar parameters.

  CHARACTER(len=optLen), PARAMETER ::  &
    routeTypeTRIP = 'trip'    !  value of routeType indicating that the TRIP
!             (aka STRIP) linear model is to be used.
!              TRIP really only refers to the river network data, but is used
!              here to refer to a linear model that has often been used with
!              the TRIP data.

!-------------------------------------------------------------------------------
! Array parameters.

  INTEGER,PARAMETER :: flowDirDelta(0:9,2) = RESHAPE( (/  &
!      missing  N    NE    E    SE    S    SW     W   NW   no_direction
         -9,    0,   1,    1,    1,   0,   -1,   -1,  -1,    0,    &!  x offsets
         -9,    1,   1,    0,   -1,  -1,   -1,    0,   1,    0    &!  y offfsets
                 /), (/ 10,2/) )
!      Components of displacement (number of gridboxes) in x and
!                   y directions to the immediately downstream gridbox. These
!                   correspond to the directions in flowDirSet. The elements in
!                   the 2nd dimension are x and y respectively.
!             e.g. flowDirDelta(1,1:2)=(/0,1/)
!             flowDirSet(1) refers to a displacement to the N
!             flowDirDelta(1,2)=1 means 1 gridbox in the N (y) direction
!             flowDirDelta(1,1)=0 means no displacement in the W-E (x) direction
!     The zeroth column (flowDirDelta(0,:)) is not currently used, but is
!     included so as to make size(1) of flowDirDelta and flowDirSet equal - 
!     in an attempt for clarity.
!     The "no defined direction" must be the last column.
!-------------------------------------------------------------------------------
! Scalar variables.

  INTEGER ::  &
    npRoute        &!  number of points in the routing grid at which routing is done
!                        This is the number of "active" routing points, e.g. not sea
   ,nxInRoute      &!  row length of input grid for routing
   ,nxRoute        &!  row length for routing grid
   ,nyInRoute      &!  column length of input grid for routing
   ,nyRoute        &!  column length for routing grid
   ,routeCount     &!  counter of timesteps done in current routing timestep
!                        Generally equals routeStep, but is < routeStep if part
!                        of a routing timestep was missed (e.g. at start of run).
   ,routeStep      &!  counter of timesteps since routing last done (actually a
!                        counter that triggers call to routing. May be > steps
!                        actually done - see routeCount).
   ,routeTimestep   !  timestep for runoff routing (number of model timesteps)

  REAL ::  &
    routeDlat      &!  size of gridbox of (regular) routing grid in latitude (degrees)
   ,routeDlon      &!  size of gridbox of (regular) routing grid in longitude (degrees)
   ,routeLat1      &!  latitude of southernmost row of gridpoints on a regular
!                        routing grid (degrees)
   ,routeLon1      &!  longitude of westernmost (first) column of gridpoints on a
!                        regular routing grid (degrees)
   ,routeMeander   &!  meander ratio for routing. This is the ratio of actual
!                         river length to the calculated river length.
   ,routeSpeed      !  flow speed for routing (m s-1)

  LOGICAL ::    &
   routeRegLatLon  &!  flag indicating if routing grid is regular in latitude and
!                         longitude See above for definition of a regular grid.
  ,routeRegrid     &!  flag indicating if model and routing grids are identical
!                        TRUE grids are identical
!                        FALSE grids differ and regridding required
   ,yrevInRoute  !  T means routing input data are arranged N to S,
!                   F means routing input data are arranged S to N

  CHARACTER(len=optLen) ::  &!  scalars
     routeType       !  choice of routing model

!-------------------------------------------------------------------------------
! Array variables.

  INTEGER :: flowDirSet(0:9)   !  values of flow direction code that represent the
!       displacement given by the corresponding positions in flowDirDelta.
!       The following are the meanings of each element of the array:
!        0: no data value (e.g. over sea)
!        1-8: N,NE,E,SE,S,SW,W,NW  i.e. clockwise from N
!             Although referred to via these compass directions, they are used as "grid-relative"
!             directions, i.e. it is assumed that columns run S-N, and rows W-E, so "N" means
!             "same column, one row up". If the grid is actually rotated (so that columns do not
!             run S-N on the earth), the point that is "same column, one row up" in fact does not
!             lie immediately N.
!        9: undefined flow direction (i.e. no outflow, or outflow to sea)
!           For some encoding schemes, there is not a single value that represents undefined flow
!           (e.g. I think ARC combines values of all directions with equal slope - so "undefined" flow
!           direction depends upon slope to neighbouring points).
!           This is not a problem here, since none of these schemes is currently encoded, BUT may have
!           to be accounted for in future - code could be added below where currently we stop with
!           "unexpected value".
!        Note that the values of flowDirSet that correspond to "real" flow
!        directions (as opposed to missing data or no defined flow direction)
!        must be >0. Further, ALL values of flowDirSet (incl sea and no flow)
!        must be >=0. In that case, values <0 can be used to indicate points that have
!        a flow direction that points off the edge of the grid (for non-global
!        applications).  These restrictions are necessary to make the current
!        algorithm work. In
!        particular, a run with a smaller grid (e.g. regional) should give the
!        same routing pathways as a run with a larger grid (at the points that
!        are on both grids). This was not the case with code that set the flow
!        direction to "no defined direction" when flow was off the edge of the
!        smaller grid.

  INTEGER, ALLOCATABLE ::  &
    mapInRoute(:)  &!  mapping from routing input grid to routing grid
!                        i.e. for each point in routing grid, gives location in
!                        input routing grid
!                        This is only used during initialsation, but as it is
!                        potentially needed by several routines, we store it.
   ,routeIndex(:)  &!  Index of routing points. This gives location on routing
!                        grid of each active routing point.
   ,routeNext(:)   &!  Index of the next downstream point.
   ,routeOrder(:)   !  Gives order in which points in vector are taken for routing.
!                        e.g. routeOrder(1)=n indicates that the first point to be
!                        routed is the nth in routeIndex. The location of this point
!                        in the routing grid is then given by routeIndex(n).
!                        For (S)TRIP, the points are taken in order of increasing
!                        "sequence" number.

! NOTE on routing grid
! At present, the routing grid must be regular in latitude and longitude.
! A grid that is regular in latitude and longitude (not rotated) can be defined by the coords
! of a single point (the origin), the gridbox size in lat and lon, and the number of rows and
! columns on the grid. Columns in a regular grid run N-S, rows run W-E, i.e. not rotated.

  REAL, ALLOCATABLE ::  &
    roffAccumLand(:)    !  average runoff (production) rate between calls to
!                           routing (kg m-2 s-1). This is held on land points.

  LOGICAL, ALLOCATABLE ::  &
    routeMask(:,:)  !  T at points on routing grid where routing is "active"
!                        (i.e. are in routing vector)
!                      F at all other points (e.g. sea)
!                        Only used during initialisation.


  END MODULE route_mod
