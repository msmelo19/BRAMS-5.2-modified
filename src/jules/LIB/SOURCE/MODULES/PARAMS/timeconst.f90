!########################################################################
! MODULE timeConst.
! This module provides variables for standard number of seconds in
! hour, day, year to avoid these numbers appearing explicitly in JULES
! code.  At present the normal JULES year is 365 days or 366 days in
! leap years.  There is an option in the JULES-in file to run with a
! 360 day year, as per the Unified Model, comprising 12 30-day months.
! ppha 19/02/2007.
!########################################################################
MODULE timeConst
  IMPLICIT NONE

  INTEGER,PARAMETER :: iSecInMin  = 60                       ! Number of seconds in a minute.
  INTEGER,PARAMETER :: iMinInHour = 60                       ! Number of minutes in an hour.
  INTEGER,PARAMETER :: iHourInDay = 24                       ! Number of hours in a day.
  INTEGER,PARAMETER :: iDayInYear = 365                      ! Number of days in a year.
  INTEGER,PARAMETER :: iSecInHour = iSecInMin  * iMinInHour  ! Number of seconds in an hour.
  INTEGER,PARAMETER :: iSecInDay  = iSecInHour * iHourInDay  ! Number of seconds in a day.
  INTEGER,PARAMETER :: iSecInYear = iSecInDay  * iDayInYear  ! Number of seconds in a year.
  INTEGER,PARAMETER :: iHourInYear = iHourInDay * iDayInYear ! Number of hours in a year.

END MODULE timeConst
