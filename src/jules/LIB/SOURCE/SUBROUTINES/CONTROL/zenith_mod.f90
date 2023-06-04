  MODULE ZENITH_MOD

! This module contains subroutine zenith and some supporting time/calendar code.

  CONTAINS

!###############################################################################
!###############################################################################

  SUBROUTINE ZENITH(NPOINTS,COSZ)

! Subroutine to calculate the cosine of the zenith angle

  USE switches, ONLY : l_360
  USE c_pi, ONLY : pi,pi_over_180
  USE ancil_info, ONLY : row_length
  USE time_loc, ONLY : ijulian,latitude,longitude,utc_time_secs,date
  USE time_mod, ONLY : dateToBits
  USE timeConst, ONLY : iSecInDay

  IMPLICIT NONE

  INTEGER                   :: NPOINTS  ! Number of points
  REAL, DIMENSION(NPOINTS)  :: COSZ     ! Cosine of zenith angle

  REAL    ::  HH     ! Hour angle of sun relative to solar noon (radians)
  REAL    ::  SD     ! Solar declination angle (radians)

  INTEGER ::  DAY    ! Current day
  INTEGER ::  MONTH  ! Current month
  INTEGER ::  YEAR   ! Current year

  INTEGER :: I,J,L   ! Loop counters

! Get current year from date
  CALL dateToBits( date,day,month,year,l_360,'Zenith' )

  DO L=1,NPOINTS
    J=(L-1)/ROW_LENGTH + 1
    I= L - (J-1)*ROW_LENGTH

    HH = PI * (2.0*REAL(UTC_time_secs)/REAL(iSecInDay) + LONGITUDE(I,J)/180.0 - 1.0)

    SD=-23.4*PI_OVER_180* COS(2.0*PI*(IJULIAN+10)/DAYS_IN_YEAR(year,L_360))

    COSZ(L)=SIN(PI_OVER_180*LATITUDE(I,J))*SIN(SD)  &
          + COS(PI_OVER_180*LATITUDE(I,J))*COS(SD)*COS(HH)

    IF(COSZ(L) < 0.0) COSZ(L)=0.0

  ENDDO

  RETURN
  END SUBROUTINE ZENITH


!###############################################################################
!###############################################################################

! Function to calculate the number of days in the given year

  INTEGER FUNCTION DAYS_IN_YEAR(IYEAR,L_360)

  IMPLICIT NONE

  INTEGER :: IYEAR   !  Year
  LOGICAL :: L_360   ! True if running with a 360 day year

  DAYS_IN_YEAR=365

! Add a day if it is a leap year
  IF((IYEAR - 4*(IYEAR/4)) == 0) THEN
    IF((IYEAR - 100*(IYEAR/100)) == 0) THEN
      IF((IYEAR - 400*(IYEAR/400)) ==0) DAYS_IN_YEAR=DAYS_IN_YEAR+1
    ELSE
      DAYS_IN_YEAR=DAYS_IN_YEAR+1
    ENDIF
  ENDIF

  IF(L_360) DAYS_IN_YEAR=360

  RETURN
  END FUNCTION DAYS_IN_YEAR


!###############################################################################
!###############################################################################

  INTEGER FUNCTION JULIAN_DAY(IDAY,IMONTH,IYEAR,L_360)

! Function to calculate the julain day given the year, month and day of month

  IMPLICIT NONE

  INTEGER :: IDAY      ! Day of month
  INTEGER :: IMONTH    ! Month
  INTEGER :: IYEAR     ! YEAR

  LOGICAL :: L_360     ! True if running with a 360 day year

  INTEGER, PARAMETER  :: NMONTHS=12  ! Number of months in a year
  INTEGER, PARAMETER :: DAYS_IN_MONTH(NMONTHS) = &
!                                        ! Number of days in each month
                    (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)

  JULIAN_DAY=0

! Add up number of days in previous months
  IF(IMONTH > 1) THEN
    IF(L_360) THEN
      JULIAN_DAY=(IMONTH-1)*30
    ELSE
      JULIAN_DAY=SUM(DAYS_IN_MONTH(1:IMONTH-1))
! Add on leap day if required
      IF(DAYS_IN_YEAR(IYEAR,L_360) == 366 .AND. IMONTH > 2) JULIAN_DAY=JULIAN_DAY+1
    ENDIF
  ENDIF

! Add on days in current month
  JULIAN_DAY=JULIAN_DAY+IDAY

  RETURN
  END FUNCTION JULIAN_DAY

!###############################################################################
!###############################################################################

  END MODULE zenith_mod




