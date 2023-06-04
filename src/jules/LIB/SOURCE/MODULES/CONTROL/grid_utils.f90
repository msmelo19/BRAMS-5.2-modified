!#######################################################################
!#######################################################################
!#######################################################################

! module grid_utils
! Contains procedures for manipulating points - e.g. mapping  between
! grids (e.g. input grid to model grid), calculating point numbers and the likes.
!
!###############################################################################

  MODULE grid_utils
  IMPLICIT NONE

! Interfaces to switch between integer or real variables.

  INTERFACE mapAtoB
  MODULE PROCEDURE mapAtoBInt,mapAtoBReal
  END INTERFACE

  INTERFACE reverseCols
  MODULE PROCEDURE reverseColsInt,reverseColsReal
  END INTERFACE

  CONTAINS

!#######################################################################
!#######################################################################
!###############################################################################
!###############################################################################
! function getLandValues
! Internal procedure in module grid_utils.
! Given a field on the (x,y) model grid, extract values at land points into a vector.

  FUNCTION getLandValues( inval ) RESULT ( landval )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts,nx=>row_length,ny=>rows  &!
!  imported arrays with intent(in)
    ,land_index

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL ::  &!  array function result
    landval(land_pts)

  INTEGER ::  &!  local scalars
    ix,iy,l

  REAL, INTENT(in) ::  &!  in arrays
    inval(nx,ny)
!-------------------------------------------------------------------------------

  DO l=1,land_pts
    CALL getXYPos( land_index(l),nx,ny,ix,iy )
    landval(l) = inval(ix,iy)
  ENDDO

  END FUNCTION getLandValues

!###############################################################################
!###############################################################################
! function mapAtoBInt
! Internal procedure in module grid_utils
! Module procedure with generic name mapAtoB.
!
! Map input array onto output array (both arrays integer).
! From an input array (length=npIn), select npMap values (not
! necessarily the first npMap values) as described by mapIn, and copy
! them to the output array (length=npOut) in positions given by mapOut.
!
! Points that are not given a new value by the mapping are given the
! default value undefOut.
! It is assumed that the mappings supplied are correct - errors could
! lead to attempts to access non-existent array elements.
!
! Examples:
! 1) getting input data for model
!    Get data for 4 model points from an input of 10 points.
!    npIn=10,npOut=4,npMap=4
!    giving inval(10),outval(4) and maps(4)
!    mapIn=1,5,7,9  selects points 1,5,7,9 in input grid.
!    mapOut=1,2,3,4 writes to points 1 to 4 in model grid
!
! 2) writing output (all of model grid)
!    Model grid of 4 points. Write 4 points in output grid of 20 points.
!    npIn=4,npOut=20,npMap=4
!    giving inval(4),outval(20) and maps(4)
!    mapIn=1,2,3,4      selects points 1-4 in input grid
!    mapOut=5,10,11,20  writes to points 5, 10,11 and 20 in output grid
!
! 3) writing output (part of model grid)
!    Model grid of 4 points. Write 3 points in output grid of 20 points.
!    npIn=4,npOut=20,npMap=3
!    giving inval(4),outval(20) and maps(3)
!    mapIn=1,2,4     selects points 1,2 and 4 in input grid
!    mapOut=5,10,20  writes to points 5, 10 and 20 in output grid

  FUNCTION mapAtoBInt( npIn,npMap,npOut,mapIn,mapOut,undefOut,inval )  &
               RESULT( outval )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  IN scalars
    npIn   &!   number of points (size) in input grid
   ,npMap  &!   number of points to map
   ,npOut   !   number of points (size) in output grid

  INTEGER ::    &!   array function result
    outval(npOut) !  data on output grid

  INTEGER, INTENT(in) ::  &!  IN arrays
    mapIn(npMap)   &!  points to select from input grid
   ,mapOut(npMap)   !  points to write to in output grid

  INTEGER, INTENT(in) :: &!  IN scalars
    undefOut   !  value written to all output points that are not
!        included in mapping e.g. if there are 10 output points and only 4 are
!        assigned values from input, the other 6 will have value undefOut

  INTEGER, INTENT(in) ::  &!   IN arrays
    inval(npIn)    !  data on input grid

  INTEGER :: i  !  local loop counter
!-------------------------------------------------------------------------------
! Check that don't have more points to map than are in output.
  IF ( npOut < npMap ) THEN
    WRITE(*,*)'mapAtoBInt cannot map npMap=',npMap,' points to npOut=',npOut,' points.'
    WRITE(*,*)'Output grid must be large enough for all mapped points.'
    WRITE(*,*)'(Mapping is 1-to-1 at all (npMap) points in map.)'
    WRITE(*,*)'Check actual arguments are correct.'
    WRITE(*,*)'mapAtoBInt: npOut < npMap'
    WRITE(*,*)'Stopping in mapAtoBInt'
    STOP
  ENDIF

! Initialise the output to all have the undefined value.
  outval(:) = undefOut

  DO i=1,npOut
    outval( mapOut(i) ) = inval( mapIn(i) )
  ENDDO

  END FUNCTION mapAtoBInt
!#######################################################################
!#######################################################################
!#######################################################################
! function mapAtoBReal
! Internal procedure in module grid_utils
! Module procedure with generic name mapAtoB.
!
! Map input array onto output array (both arrays real).
! From an input array (length=npIn), select npMap values (not
! necessarily the first npMap values) as described by mapIn, and copy
! them to the output array (length=npOut) in positions given by mapOut.
!
! Points that are not given a new value by the mapping are given the
! default value undefOut.
! It is assumed that the mappings supplied are correct - errors could
! lead to attempts to access non-existent array elements.
!
! Examples:
! 1) getting input data for model
!    Get data for 4 model points from an input of 10 points.
!    npIn=10,npOut=4,npMap=4
!    giving inval(10),outval(4) and maps(4)
!    mapIn=1,5,7,9  selects points 1,5,7,9 in input grid.
!    mapOut=1,2,3,4 writes to points 1 to 4 in model grid
!
! 2) writing output (all of model grid)
!    Model grid of 4 points. Write 4 points in output grid of 20 points.
!    npIn=4,npOut=20,npMap=4
!    giving inval(4),outval(20) and maps(4)
!    mapIn=1,2,3,4      selects points 1-4 in input grid
!    mapOut=5,10,11,20  writes to points 5, 10,11 and 20 in output grid
!
! 3) writing output (part of model grid)
!    Model grid of 4 points. Write 3 points in output grid of 20 points.
!    npIn=4,npOut=20,npMap=3
!    giving inval(4),outval(20) and maps(3)
!    mapIn=1,2,4     selects points 1,2 and 4 in input grid
!    mapOut=5,10,20  writes to points 5, 10 and 20 in output grid

  FUNCTION mapAtoBReal( npIn,npMap,npOut,mapIn,mapOut,undefOut,inval )  &
                RESULT( outval )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  IN scalars
    npIn   &!   number of points (size) in input grid
   ,npMap  &!   number of points to map
   ,npOut   !   number of points (size) in output grid

  REAL ::    &!  array function result
    outval(npOut) !  data on output grid

  INTEGER, INTENT(in) ::  &!  IN arrays
    mapIn(npMap)  &!  points to select from input grid
   ,mapOut(npMap)  !  points to write to in output grid

  REAL, INTENT(in) ::  &!  IN scalars
    undefOut      !  value written to all output points that are not
!     included in mapping e.g. if there are 10 output points and only 4 are
!     assigned values from input, the other 6 will have value undefOut

  REAL, INTENT(in) ::  &!   IN arrays
    inval(npIn)    !  data on input grid

  INTEGER :: i  !  local loop counter
!-------------------------------------------------------------------------------
!  print*,'mapAtoB: npIn,npMap,npOut=',npIn,npMap,npOut
!  print*,'mapIn=',mapIn
!  print*,'mapOut=',mapOut
!  print*,'inval=',inval

! Check that don't have more points to map than are in output.
  IF ( npOut < npMap ) THEN
    WRITE(*,*)'mapAtoBReal cannot map npMap=',npMap,' points to npOut=',npOut,' points.'
    WRITE(*,*)'Output grid must be large enough for all mapped points.'
    WRITE(*,*)'(Mapping is 1-to-1 at all (npMap) points in map.)'
    WRITE(*,*)'Check actual arguments are correct.'
    WRITE(*,*)'mapAtoBReal: npOut < npMap'
    WRITE(*,*)'Stopping in mapAtoBReal'
    STOP
  ENDIF

! Initialise the output to all have the undefined value.
  outval(:) = undefOut

  DO i=1,npMap
    outval( mapOut(i) ) = inval( mapIn(i) )
  ENDDO

  END FUNCTION mapAtoBReal
!#######################################################################
!#######################################################################
!#######################################################################
! function reverseColsInt
! Module procedure with generic name reverseCols
! Take a 2-D (x,y) field and reverse the order of the rows, that is reverse each y column.
! i.e. var(nx,ny), swap var(:,1) and var(:,nx), swap var(:,2) and var(:,nx-1) etc.

  FUNCTION reverseColsInt(inval) RESULT( outval )

  IMPLICIT NONE

  INTEGER  ::   &!  local SCALARS
    iy   &!  loop counter
   ,ny    !  extent of 2nd dimension of arrays

  INTEGER, INTENT(in) ::  &!  in arrays
    inval(:,:)   !  input data

  INTEGER ::   &!   array function result
    outval(SIZE(inval,1),SIZE(inval,2))   !  output data
!-------------------------------------------------------------------------------

  ny = SIZE(inval,2)

  DO iy=1,ny
    outval(:,iy) = inval(:,ny-iy+1)
  ENDDO

  END FUNCTION reverseColsInt
!#######################################################################
!#######################################################################
!#######################################################################

! function reverseColsReal
! Module procedure with generic name reverseCols
! Take a 2-D (x,y) field and reverse the order of the rows, that is reverse each y column.
! i.e. var(nx,ny), swap var(:,1) and var(:,nx), swap var(:,2) and var(:,nx-1) etc.

  FUNCTION reverseColsReal(inval) RESULT( outval )

  IMPLICIT NONE

  INTEGER  ::   &!  local SCALARS
    iy   &!  loop counter
   ,ny    !  extent of 2nd dimension of arrays

  REAL, INTENT(in) ::   &!  in arrays
    inval(:,:)   !  input data

  REAL ::    &!   array function result
    outval(SIZE(inval,1),SIZE(inval,2))   !  output data
!-------------------------------------------------------------------------------
  ny = SIZE(inval,2)

  DO iy=1,ny
    outval(:,iy) = inval(:,ny-iy+1)
  ENDDO

  END FUNCTION reverseColsReal
!###############################################################################
!###############################################################################
! function getGridPosLL
! Internal procedure in module grid_utils.
! Given the latitude and longitude of gridpoints, calculate location in
! a grid.
! NB The locations must be gridpoints, rather than just on the grid.

  FUNCTION getGridPosLL( nx,ny,latVal,lonVal,lat1,lon1,dlat,dlon  &
                        ,latitude,longitude)  &
           RESULT( gridIndex )

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar arguments with intent(in)

  INTEGER, INTENT(in) :: nx   !  x size of grid
  INTEGER, INTENT(in) :: ny   !  y size of grid

! Array arguments with intent(in).

  REAL, INTENT(in) :: latVal(:)   !  latitude of points
  REAL, INTENT(in) :: lonVal(:)   !  longitude of points

!-------------------------------------------------------------------------------
! Array function result.

  INTEGER :: gridIndex(SIZE(latVal))   !  index (point number) in the grid

!-------------------------------------------------------------------------------
! Optional scalar arguments with intent(in).

  REAL, INTENT(in), OPTIONAL :: lat1   !  latitude of point (1,1) (gridpoint in
!                SW corner) in grid
  REAL, INTENT(in), OPTIONAL :: lon1   !  longitude of point (1,1) (gridpoint in
!                SW corner) in grid
  REAL, INTENT(in), OPTIONAL :: dlat   !  gridbox size in latitude (degrees)
  REAL, INTENT(in), OPTIONAL :: dlon   !  gridbox size in longitude (degrees)

! Optional array arguments with intent(in).

  REAL, INTENT(in), OPTIONAL :: latitude(:,:)   !  latitudes of gridpoints
  REAL, INTENT(in), OPTIONAL :: longitude(:,:)  !  longitudes of gridpoints

!-------------------------------------------------------------------------------
! Scalar variables.

  INTEGER :: ip,ix,iy  !  worl/loop counters
  REAL :: dx,dy,lonWork   !  work

  LOGICAL :: lon360    !  flag indicating range of values for longitude
!      TRUE means longitude is in range 0 to 360 degrees
!      FALSE means longitude is in range -180 to 180 degrees

  LOGICAL :: regGrid   !  flag indicating if grid is "regular". In this context,
!       a regular grid means an equal angle grid which can be specified by the
!       size of the grid, the gridbox size and the location of a point on the
!       grid. Rotated grids are NOT included.
!          
!-------------------------------------------------------------------------------

! Initialise result to show points are not in grid. 
  gridIndex(:) = 0

! Determine if grid is "regular".
! Assume that the presence of one optional argument indicates that all the correct
! optional arguments have been provided!
  regGrid = .FALSE.
  IF ( PRESENT(lat1) ) regGrid = .TRUE.

!-------------------------------------------------------------------------------
! Decide what "convention" (range) is used for longitude values.
  lon360 = .TRUE.  
  IF ( regGrid ) THEN
    lonWork = lon1
  ELSE
!   Get minimum longitude on grid.
    lonWork = MINVAL( longitude(:,:) )
  ENDIF
  IF ( lonWork < 0.0 ) lon360 = .FALSE.

  DO ip=1,SIZE( latVal(:) )

!   Get longitude in same "convention" as is used for grid.
    lonWork = lonVal(ip)
    IF ( lonWork < 0.0 ) THEN
      IF ( lon360 ) lonWork=360.0+lonWork
    ELSEIF ( lonWork > 180.0 ) THEN
      IF ( .NOT. lon360 ) lonWork=lonWork-360.0
    ENDIF

!   Check if this location is a grid point.
    IF ( regGrid ) THEN
!     Compare coords with those of origin of grid.
      dy = ( latVal(ip) - lat1 ) / dlat
      dx = ( lonWork - lon1 ) / dlon
      IF ( ABS(dx-REAL(NINT(dx)))<=EPSILON(dx) .AND. ABS(dy-REAL(NINT(dy)))<=EPSILON(dy) ) THEN
!       This location is a gridpoint. Locate in grid.
        ix = NINT( dx ) + 1
        iy = NINT( dy ) + 1
        IF ( ix<1 .OR. ix>nx .OR. iy<1 .OR. iy>ny ) THEN
!         Location is off edge of grid.
          gridIndex(ip) = 0
        ELSE
!         Calculate index.
          gridIndex(ip) = (iy-1)*nx + ix
        ENDIF
      ENDIF
!-------------------------------------------------------------------------------
    ELSE

!     Not a regular grid. Loop through coordinates until find a match.
      DO iy=1,ny
        DO ix=1,nx
          IF ( ( ABS(latVal(ip)-latitude(ix,iy)) < EPSILON(latVal) ).AND.  &
               ( ABS(lonWork-longitude(ix,iy)) < EPSILON(lonVal) ) ) THEN
            gridIndex(ip) = (iy-1)*nx + ix
            CYCLE
          ENDIF
        ENDDO
      ENDDO
     
    ENDIF   !  routeRegLatLon

  ENDDO   !  points

  END FUNCTION getGridPosLL

!###############################################################################
!###############################################################################
! subroutine getXYpos
! Internal procedure in module grid_utils
! Given a point number and the extents (size) of a 2-D (xy) grid, returns the
! x and y indices (coords) of the point. All coords are relative to (1,1)
! at bottom left of grid, with numbering running left to right, bottom to top.

  SUBROUTINE getXYpos( i,nx,ny,ix,iy )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    i   &!  the point number
   ,nx  &!  the extent (size) of the grid in the x direction
   ,ny   !  the extent (size) of the grid in the y direction
!             NB ny is only required for (cursory) error checking.
    
  INTEGER, INTENT(out) ::  &!  out SCALARS
    ix,iy   !  the x and y coordinates of the point
!-------------------------------------------------------------------------------
! Assume that i>=1!

!  iy = ceiling( real(i) / real(nx) )
! Or avoiding function call....
  iy = ( i - 1 ) / nx + 1
  
  ix = i - (iy-1)*nx

! If locations are out of range (grid is too small for this point number), set to -1.
  IF ( ix>nx .OR. iy>ny ) THEN
    ix = -1
    iy = -1
  ENDIF

  END SUBROUTINE getXYpos
!#######################################################################
!###############################################################################
!###############################################################################
!###############################################################################


  END MODULE grid_utils
