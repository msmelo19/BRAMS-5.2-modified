! module earth_utils
! Contains declarations and procedures used to calculate distances on earth etc.

! Contains:
!   variable declarations
!   function earthArea
!   function giveEarthRadius
!   function giveLength
!   function giveLatLength
!   function giveLonLength

!###############################################################################
  MODULE earth_utils

  USE c_pi, ONLY :  &!
!  imported scalar parameters
     pi,pi_over_180

  IMPLICIT NONE

  REAL, PARAMETER ::  &!  scalar parameters
   earthRadius = 6378.136  &!  radius of Earth (km), from Appendix A of
!                                Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.
  ,eccenSq = 0.00669447     !  square of eccentricity of Earth spheroid
!                                (spheroid - lines of constant latitude are circular,
!                                meridional cross-section is an ellipse). From Appendix A of
!                                Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

  CONTAINS

!###############################################################################
!###############################################################################
!###############################################################################
! function giveLength
! Internal procedure in module earth_utils.
! Calculates distance on surface of Earth between two locations, assuming
! Earth is a spheroid.
! Based on function giveLen by Taikan Oki (26/August/1996).

  FUNCTION giveLength( lat1,lon1,lat2,lon2 ) result ( length )

  REAL ::  &!  scalar function result
    length   !  distance (km)

  REAL, INTENT(in) ::  &!  in scalars
    lat1,lat2,lon1,lon2   !  latitude and longitude of the two points (degrees)

  REAL ::  &!  local scalars
    dlat   &!  difference in latitude (degrees)
   ,dlon   &!  difference in longitude (degrees)
   ,dx     &!  work
   ,dy     &!  work
   ,lat    &!  work: latitude (degrees)
!   ,lat1d,lat2d,lon1d,lon2d   &!  latitude and longitude of the two points (radians)
   ,radius  !  equivalent radius of Earth at given latitude (km)

!-------------------------------------------------------------------------------
! Convert degrees to radians.
!  lat1 = lat1d * pi_over+180
!  lat2 = lat2d * pi_over+180
!  lon1 = lon1d * pi_over+180
!  lon2 = lon2d * pi_over+180

  dlon = ABS( lon2 - lon1 )
  IF ( dlon>=180.0 ) dlon = 360.0 - dlon
  dlat = ABS( lat2 - lat1 )

  IF ( dlon < EPSILON( dlon ) ) THEN
!   Constant longitude.
    lat = ( lat1 + lat2 ) * 0.5
    length = giveLatLength(lat) * dlat
  ELSEIF ( dlat < EPSILON(dlat) ) THEN
!   Constant latitude.
    length = giveLonLength( lat1 ) * dlon
  ELSE
!   Both lat and lon change.
!   Use equation A8 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1
    lat = ( lat1 + lat2 ) * 0.5
    radius = giveEarthRadius( lat )
    dx = giveLonLength( lat ) * dlon / radius
    dy = giveLatLength( lat ) * dlat / radius
    length = ACOS(COS(dx)*COS(dy)) * radius
  ENDIF

  END FUNCTION giveLength

!###############################################################################  
!###############################################################################  
! function giveLatLength
! Internal procedure in module earth_utils.
! Calculates the distance (km) along the surface of the Earth between two points
! separated by 1 degree of latitude and at the same longitude.
! Based on function givelat by Taikan Oki (23/April/1996).
! See EqnA2 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

  FUNCTION giveLatLength( lat )  result (latLen)

  IMPLICIT NONE

  REAL ::  &!  scalar function result
    latLen  !  distance between points (km) (actually same units as earthRadius)

  REAL, INTENT(in) ::  &!  in scalars
    lat   !  effective (e.g. average) latitude (degrees)

  latLen = pi_over_180 * earthRadius * ( 1 - eccenSq)  &
           / (SQRT(1.0 - eccenSq*SIN(lat*pi_over_180)*SIN(lat*pi_over_180)))**3

  END FUNCTION giveLatLength
!###############################################################################
!###############################################################################
! function giveLonLength
! Internal procedure in module earth_utils.
! Calculates the distance (km) along the surface of the Earth between two points
! separated by 1 degree of longitude and at the same latitude.
! Based on function givelat by Taikan Oki (23/April/1996).
! See EqnA2 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

  FUNCTION giveLonLength( latDeg )  result (lonLen)

  IMPLICIT NONE

  REAL ::  &!  scalar function result
    lonLen    !  length (km) (actually same units as earthRadius)

  REAL, INTENT(in) ::  &!  in scalars
    latDeg   !  latitude (degrees)

  REAL ::  &!  local scalars
    lat  !  latitude (radians)
!-------------------------------------------------------------------------------

  lat = latDeg * pi_over_180

  lonLen = pi_over_180 * earthRadius * COS(lat)  &
           / SQRT(1.0 - eccenSq*SIN(lat)*SIN(lat))

  END FUNCTION giveLonLength
!###############################################################################
!###############################################################################
!###############################################################################
! function giveEarthRadius
! Internal procedure in module earth_utils.
! Calculates equivalent radius of Earth at a given latitude, assuming the Earth
! is a spheroid.
! From equations A10 and A11 in Appendix A of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.
! Based on function giverade by Taikan Oki (26/August/1996).

  FUNCTION giveEarthRadius( latDeg ) result( radius )

  IMPLICIT NONE

  REAL ::  &!  scalar function result
    radius   !  equivalent radius (km) (actually same units as earthRadius)

  REAL, INTENT(in) ::  &!  in scalars
    latDeg   !  latitude (degrees)

  REAL ::  &!  local scalars
    lat &!  latitude (radians)
   ,rn    !  work
!-------------------------------------------------------------------------------

  lat = latDeg * pi_over_180
  rn = earthRadius / SQRT(1.0 - eccenSq*SIN(lat)*SIN(lat))
  radius = rn * SQRT(1.0 - 2.0*eccenSq*SIN(lat) + eccenSq*eccenSq*SIN(lat)*SIN(lat))

  END FUNCTION giveEarthRadius
!###############################################################################
!###############################################################################
!###############################################################################
! function earthArea
! Internal procedure in module earth_utils.
! Calculates area of Earth surface between two lines of latitude and two of longitude,
! assuming that the Earth is a spheriod.
! Uses equation A6 of Oki and Sud, 1998, Earth Interactions, Vol.2, Paper 1.

  FUNCTION earthArea( lat1d,lat2d,lon1d,lon2d ) result( area )

  IMPLICIT NONE

  REAL ::  &!  scalar function result
    area   !  area (km2) (actually in units of earthRadius**2)

  REAL, INTENT(in) ::  &!  in scalars
    lat1d    &!  latitude of southern edge of strip (degrees)
   ,lat2d    &!  latitude of northern edge of strip (degrees)
   ,lon1d    &!  longitude of western edge of strip (degrees)
   ,lon2d     !  longitude of eastern edge of strip (degrees)

  REAL ::  &!  local scalars
    eccen        &!  eccentricity of Earth spheroid
   ,lat1,lat2    &!  lat1 and lat2 in radians
   ,val1,val2,val3     !  work

!-------------------------------------------------------------------------------
  eccen = SQRT( eccenSq )

! Convert latitudes to radians.
  lat1 = lat1d * pi_over_180
  lat2 = lat2d * pi_over_180

! Evaluate terms at each of lat1 and lat2.
  val1 =  0.5*eccen*SIN(lat1) / ( 1.0-eccenSq*SIN(lat1)*SIN(lat1) )
  val1 = val1 + 0.25 * LOG( ABS( (1.0+eccen*SIN(lat1) ) / (1.0-eccen*SIN(lat1) ) ) )

  val2 =  0.5*eccen*SIN(lat2) / ( 1.0-eccenSq*SIN(lat2)*SIN(lat2) )
  val2 = val2 + 0.25 * LOG( ABS( (1.0+eccen*SIN(lat2) ) / (1.0-eccen*SIN(lat2) ) ) )

  val3 = ABS(lon2d-lon1d) * pi_over_180 * earthRadius * earthRadius * (1.0-eccenSq) / eccen

  area = val3 * ( val2 - val1 )

  END FUNCTION earthArea

!###############################################################################
!###############################################################################
!###############################################################################
  END MODULE earth_utils
