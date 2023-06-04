!###############################################################################

! Module containing variables relating to time or location of grid.

  MODULE time_loc

  IMPLICIT NONE

!###############################################################################

! Time variables.

!-------------------------------------------------------------------------------
! Scalars.
!-------------------------------------------------------------------------------
  INTEGER ::   &
    date       &!  date (yyyymmdd)
   ,dateNext   &!  date on next timestep (assuming end of spin up does not intervene)
   ,datePrev   &!  date on previous timestep
   ,IJULIAN    &!  Julian day
   ,stepFlag   &!  used to indicate "special" timesteps (e.g. first in a cycle of spin up)
!                    Used in output, but bit of a faff.
!                    0  default
!                    1  first timestep in a new cycle of spin up
!                    3  first timestep in the "main" run (after any spin up)
   ,time       &!  time of day (UTC, s)
   ,timeNext   &!  time on next timestep (assuming end of spin up does not intervene)
   ,timePrev    !  time on previous timestep

  REAL ::          &
    timestep       &!  timestep length (seconds)
   ,utc_time_secs   !  UTC time at end of timestep (seconds)

  LOGICAL ::    &
    endMonth    &!  T on the last timestep in a month (1st timestep that ends
!                     at t>=0H 1st next month)
   ,endSec      &!  T on the last timestep in any 'section' of the run
!                       The "main" run and each cycle of spin up are considered separate sections.
!                       Used for output - a new section often triggers e.g. a new output file.
   ,endYear     &!  T on the last timestep in a year (1st timestep that ends
!                      at t>=0H 1st Jan next year)
   ,newMonth    &!  T on the first timestep in a month (generally 1st timestep after endMonth)
   ,newSec      &!  T on the first timestep in any 'section' of the run
   ,newYear      !  T on the first timestep in a year (generally 1st timestep after endYear)
!    Note that both endMonth/Year and newMonth/Year can be true at once - e.g. if run
!    starts at end of month, newMonth is also set because it is start of run.

!-------------------------------------------------------------------------------
! Arrays.
!-------------------------------------------------------------------------------
  INTEGER ::        &
    dateMainRun(2)  &!  start and end dates for the 'main' section of the run (yyyymmdd)
!                         (excluding any spin up)
   ,dateRun(2)      &!  start and end dates for the run (yyyymmdd) (including any spin up)
   ,dateSpin(2)     &!  start and end dates for the spin up. Note that spin up is over a complete
!                    !    number of days, i.e. timeRun(1) on dateSpin(1) to timeRun(1) on dateSpin(2).
   ,timeRun(2)       !  UTC time of day at start and end of run (s)


!###############################################################################

! Location / space variables.

! A grid that is regular in latitude and longitude (not rotated) can be defined
! by the coords of a single point (the origin - or more likely the point (1,1)),
! the gridbox size in lat and lon, and the number of rows and
! columns on the grid. Columns in a regular grid run S-N, rows run W-E, and the
! points use the standard JULES order - across rows W to E, starting at S.

!-------------------------------------------------------------------------------
! Scalars.
!-------------------------------------------------------------------------------

! nxGrid and nyGrid are required by offline (uncoupled) applications in which
! the "main" model grid can be a vector of points from a larger, "regular" grid
! - e.g. land points from a larger grid. In this case, they specify the size of
! the larger grid from which the model points have been extracted.

  integer :: nxGrid  !  size of grid in x direction (row length, or number of columns)
  integer :: nyGrid  !  size of grid in y direction (column length, or number or rows)

  REAL ::     &
    regDlat   &!  gridbox size in latitude, for a regular grid (as defined
!                   above - also used for some other grids) (degrees)
   ,regDlon   &!  gridbox size in longitude, for a regular grid (as defined
!                   above - also used for some other grids) (degrees)
   ,regLat1   &!  latitude of southernmost (first) row of gridpoints
!                   If  model grid is a subset of a regular grid, this is the
!                   first row in the regular grid.
   ,regLat1In &!  latitude of southernmost (first) row of gridpoints in the
!                   input grid (which is a regular grid, as defined above)
!                   (degrees)
   ,regLon1   &!  longitude of westernmost (first) column of gridpoints
!                   If  model grid is a subset of a regular grid, this is the
!                   first column in the regular grid.
   ,regLon1In  !  longitude of westernmost (first) column of gridpoints in the
!                   input grid (which is a regular grid, as defined above)
!                   (degrees)

  LOGICAL ::     &
    regLatLon    &!  flag indicating a regular lat/lon grid
   ,subRegLatLon  !  flag indicating if the model points are a subset of
!                      a regular grid (as defined above). They needn't be stored
!                      in the "standard" W-E S-N order demanded of a regular grid.
!                      A subset can equal the full set, i.e. all points on the
!                      regular grid may be present.

!-------------------------------------------------------------------------------
! Arrays.
!-------------------------------------------------------------------------------
  REAL, ALLOCATABLE ::  &!
    latitude(:,:)       &!  Latitude of model points
   ,longitude(:,:)       !  Longitude of model points

  END MODULE time_loc
!###############################################################################
!###############################################################################
