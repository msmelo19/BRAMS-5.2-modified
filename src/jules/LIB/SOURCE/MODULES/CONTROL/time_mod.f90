! module time_mod
! Contains various date and time procedures.
!
! Contains:
!   subroutine timeDate
!   function dInMonth
!   subroutine dateToBits
!   subroutine timeToBits
!   function hourOfYear
!   function SecOfYear
!   function daysInYear
!   function dayOfYear
!   function getHours
!   function monthName3
!   function timeDate_cmp
!   subroutine timeDate_diff
!   subroutine secToSMH
!   function s_to_chhmmss
!   function chhmmss_to_s
!   function get_period

  MODULE time_mod

  IMPLICIT NONE

  CONTAINS


!#######################################################################
! subroutine timeDate
! Internal procedure in module time_mod.
! Given time of day (seconds) and date, calculate time and date a
! given number of seconds or months before or after.

  SUBROUTINE timeDate( timeIn,dateIn,inc,units,l_360,timeOut,dateOut,errMess )

  USE timeConst, ONLY : iSecInDay

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  IN scalars
    dateIn   &!  input date (YYYYMMDD)
   ,inc      &!  increment to apply
   ,timeIn    !  input time of day (s)

  INTEGER, INTENT(out) ::  &!  OUT scalars
    dateOut   &!  date (YYYYMMDD)  
   ,timeOut    !  output time of day (s)

  INTEGER ::  &!  local SCALARS  
    dateYear,dateMonth,dateDay,dday,ddayOld,incMon,incSec,tmpTime

  LOGICAL, INTENT(in) ::  &!  IN scalars
    l_360  !  T if 360 day calendar used

  LOGICAL ::  &!  local scalars
    forward  !  T if inc>0 and units are seconds

  CHARACTER(len=*), INTENT(in) ::  &!  IN scalars
    errMess  &!  message printed on error  
   ,units     !  units of increment
!                  Must be either mon,MON (for months) or sec or SEC (seconds)
!-----------------------------------------------------------------------

  incMon = 0
  incSec = 0

  IF ( units=='mon' .OR. units=='MON' ) THEN
    incMon = inc
  ELSEIF ( units=='sec' .OR. units=='SEC' ) THEN
    incSec = inc
  ELSE
    WRITE(*,*)'ERROR: No units or wrong unit supplied to timeDate.'
    WRITE(*,*)'Units=',units
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

! Break date into components.
  CALL dateToBits(dateIn,dateDay,dateMonth,dateYear,l_360,errMess)

! Default output time is as input.
  timeOut = timeIn

! Proceed now according to what type of increment is specified.
! If inc=0 this IF block is avoided and date bits remain unchanged.

  IF ( incSec /= 0 ) THEN
!   Increment by given number of seconds.
    dday = 0

!   Decide whether time is in the past or ahead of input time.
    forward = .TRUE.
    IF ( incSec < 0 ) forward=.FALSE.

!   Start by incrementing the time of day.
    timeOut = timeIn + incSec

!   Establish whether day has changed.
    IF ( timeOut<0 .OR. timeOut>=iSecInDay ) THEN
      tmpTime = timeOut
      timeOut = MODULO(timeOut,iSecInDay)
      dday = ( tmpTime - timeOut ) / iSecInDay
    ENDIF

!   If day has changed, month and year may have changed.
    IF ( dday /= 0 ) THEN

!   Loop from current month and day, until ddays are exhausted.
    DO
      ddayOld = dday
!     Going forward, reduce dday by days to end of month.
!     Going backward, increase dday by days to start of month.
      IF ( forward ) THEN
         dday = dday - ( dInMonth(dateMonth,dateYear,l_360) - dateDay )
       ELSE
         dday = dday + dateDay - 1
       ENDIF

       IF ( (forward .AND. dday<=0) .OR. (.NOT.forward .AND. dday>=0) ) THEN
!        Have found correct month. Find day and exit.
         IF ( forward ) THEN
           dateDay = dateDay + ddayOld
         ELSE
           dateDay = dateDay + ddayOld  ! + 1
         ENDIF
         EXIT
       ELSE
         IF ( forward ) THEN
           dateMonth = dateMonth + 1
           dateDay = 1
           dday = dday - 1  !  -1 since have already set day to 1
           IF ( dateMonth == 13 ) THEN
             dateMonth = 1
             dateYear = dateYear + 1
           ENDIF
         ELSE  !  not forward
           dateMonth = dateMonth - 1
           IF ( dateMonth == 0 ) THEN
             DateMonth = 12
             dateYear = dateYear - 1
           ENDIF
           dateDay = dInMonth(dateMonth,dateYear,l_360)
           dday = dday + 1 ! +1 since have already set to month end
         ENDIF !  forward
       ENDIF
     ENDDO

     ENDIF  !  dday ne 0

  ELSEIF ( incMon/=0 ) THEN

!   Increment by given number of months.
    dateMonth = dateMonth + incMon
    IF ( dateMonth > 12 ) THEN
      dateYear = dateYear + (dateMonth-1)/12  !  truncates
      dateMonth = MOD(dateMonth,12)
      IF ( dateMonth == 0 ) dateMonth=12
    ELSEIF ( dateMonth < 1 ) THEN
      dateYear = dateYear + dateMonth/12 - 1
      dateMonth = 12 + MOD(dateMonth,12)
    ENDIF
    dateDay = MIN( dateDay, dInMonth(dateMonth,dateYear,l_360) )

  ENDIF  !  increment type

  dateOut = dateYear*10000 + dateMonth*100 + dateDay

  END SUBROUTINE timeDate

!#######################################################################
!#######################################################################

! function dInMonth
! Internal procedure in module time_mod.
! Get the number of days in a given month.

  FUNCTION dInMonth(month,year,l_360) RESULT(day)

  IMPLICIT NONE
  INTEGER :: day  !  function result
  INTEGER, INTENT(in) :: month,year
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used
  INTEGER :: imonth

  INTEGER, PARAMETER ::  &!  13th month is Feb for leap year  
       days(13)=(/31,28,31,30,31,30,31,31,30,31,30,31,29/)
!-------------------------------------------------------------------------------
  IF ( l_360 ) THEN
!   Assume 360 day year is 12*30 day months.
    day = 30
  ELSE
    imonth = month
!   Deal with February in leap years.
    IF ( imonth==2 .AND. leapYear(year,l_360) ) imonth=13
    day = days(imonth)
  ENDIF

  END FUNCTION dInMonth
!#######################################################################

! function monthName3
! Internal procedure in module time_mod.
! Get a character variable for month given month number.
! At present limited to a fixed length result - failed to implement a
! version in which length of result was specified by an argument,
! possibly because it's a Friday afternoon.

  FUNCTION monthName3( imonth ) RESULT( cmonth )

  IMPLICIT NONE

  CHARACTER(len=3) :: cmonth  !  function result

  INTEGER, INTENT(in) ::  &!  IN scalars
     imonth !  month number (1-12)

  CHARACTER(len=9) ::  &!  local arrays
     fullMonth(12)=   (/'  january',' february','    march','    april'   & 
                       ,'      may','     june','     july','   august'   &
                       ,'september','  october',' november',' december'/) &
    ,xmonth  !  workspace
!-------------------------------------------------------------------------------
  xmonth = fullMonth(imonth)
  cmonth = ADJUSTL(xmonth)

  END FUNCTION monthName3
!#######################################################################
! subroutine dateToBits
! Internal procedure in module time_mod.
! Given date (YYYYMMDD), returns YYYY, MM and DD.

  SUBROUTINE DateToBits( date,dateDay,dateMonth,dateYear,l_360,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  IN scalars
    date

  INTEGER, INTENT(out) ::  &! OUT scalars
    dateDay,dateMonth,dateYear

  INTEGER :: ierr  !    (0=no error)

  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used

  CHARACTER(len=*), INTENT(in) :: errMess  !  error message
!-------------------------------------------------------------------------------
  ierr = 0
  dateYear = date / 10000
  dateMonth = ( date - dateYear*10000 ) / 100
  IF ( dateMonth<1 .OR. dateMonth>12 ) THEN
    ierr = 1
    WRITE(*,*)'ERROR: dateToBits: dateMonth=',DateMonth
    WRITE(*,*)'Month is out of range 1-12'
  ENDIF
  dateDay = date - dateYear*10000 - dateMonth*100
  IF ( dateDay < 1 .OR. dateDay>dInMonth(dateMonth,dateYear,l_360) ) THEN
    ierr = 1
    WRITE(*,*)'DateToBits: day=',DateDay,' month=',DateMonth
    WRITE(*,*)'ERROR: day number not valid for this month.'
  ENDIF

  IF ( ierr /= 0 ) THEN
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

  END SUBROUTINE dateToBits
!#######################################################################
!#######################################################################
! subroutine timeToBits
! Internal procedure in module time_mod.
! Given time of day (hhmmss) returns hh,mm,ss.

  SUBROUTINE timeToBits(time,sec,minute,hour,errMess)    

  IMPLICIT NONE

  INTEGER, INTENT(in) :: time
  INTEGER, INTENT(out) :: hour,minute,sec
  INTEGER :: ierr   !  local
  CHARACTER(len=*), INTENT(in) :: errMess  !  error message
!-------------------------------------------------------------------------------

  ierr = 0
  hour = time / 10000
  minute = ( time - hour*10000 ) / 100
  sec = time - hour*10000 - minute*100

  IF ( hour<0 .OR. hour>23 ) THEN
    ierr = 1
    WRITE(*,*) 'timeToBits: hour out of range.'
  ENDIF
  IF ( minute<0 .OR. minute>59 ) THEN
    ierr = 1
    WRITE(*,*) 'timeToBits: minute out of range.'
  ENDIF
  IF ( sec<0 .OR. sec>59 ) THEN
    ierr = 1
    WRITE(*,*) 'timeToBits: second out of range.'
  ENDIF

  IF ( ierr /= 0 ) THEN
    WRITE(*,*)'ERROR: timeToBits'
    WRITE(*,*)'Input time (hhmmss)=',time
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

  END SUBROUTINE timeToBits
!#######################################################################
! function hourOfYear
! Internal procedure in module time_mod.
! Given year,month, day of month and hour of day, returns hour of year.

  FUNCTION HourOfYear( hourOfDay,dateDay,dateMonth,dateYear,l_360 ) &
               RESULT( hour )

  USE timeConst, ONLY : iHourInDay

  IMPLICIT NONE

  REAL :: hour  !  function result
  INTEGER, INTENT(in) :: dateDay,DateMonth,DateYear
  INTEGER :: day  !  work
  REAL, INTENT(in) :: HourOfDay
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used

!-------------------------------------------------------------------------------
! Get day of year.
  day = dayOfYear( dateDay,dateMonth,dateYear,l_360 )
    
  hour = REAL(day-1)*REAL(iHourInDay) + hourOfDay

  END FUNCTION HourOfYear
!#######################################################################
!#######################################################################
!#######################################################################
! function secOfYear
! Internal procedure in module time_mod.
! Given year,month,day of month, and h,m,s of day, returns seconds of year.

  FUNCTION secOfYear( s,m,h,dateDay,dateMonth,dateYear,l_360 )  &
              RESULT( sec )

  USE timeConst, ONLY : iHourInDay,iSecInHour,iSecInMin

  IMPLICIT NONE
  INTEGER :: sec  !  function result
  INTEGER, INTENT(in) :: dateDay,dateMonth,dateYear,h,m,s
  INTEGER :: day  !  work
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used
!-------------------------------------------------------------------------------
! Get day of year.
  day = dayOfYear( dateDay,dateMonth,dateYear,l_360 )

  sec = ((day-1)*iHourInDay + h)*iSecInHour + m*iSecInMin + s

  END FUNCTION secOfYear
!#######################################################################
!#######################################################################
! function dayOfYear
! Internal procedure in module time_mod.
! Given year,month and day of month, returns day of year.

  FUNCTION dayOfYear( dateDay,dateMonth,dateYear,l_360 ) &
              RESULT( day )

  IMPLICIT NONE

  INTEGER :: day  !  function result

  INTEGER, PARAMETER ::  &!  local ARRAYS  
      daysInMonth(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)

  INTEGER, INTENT(in) :: dateDay,dateMonth,dateYear
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used
!-------------------------------------------------------------------------------

  IF ( l_360 ) THEN
    day = (dateMonth-1)*30 + dateDay

  ELSE
!   Count the days in earlier months.
    day = 0
    IF ( dateMonth > 1 ) THEN
      day = SUM( daysInMonth(1:dateMonth-1) )
!     Account for leap years.
      IF ( leapYear(dateYear,l_360) .AND. dateMonth>2 ) day=day+1
    ENDIF
!   Add days in the present month.
    day = day + dateDay
  ENDIF

  END FUNCTION dayOfYear
!#######################################################################
! function daysInYear
! Internal procedure in module time_mod.
! Given the year, returns number of days in year.

  FUNCTION daysInYear(year,l_360) RESULT(days)

  USE timeConst, ONLY : iDayInYear

  IMPLICIT NONE
  INTEGER :: days  !  function result
  INTEGER, INTENT(in) :: year
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used
!-------------------------------------------------------------------------------
  IF ( l_360 ) THEN

    days = 360

  ELSE

    days = iDayInYear
!   Is this a leap year?
    IF ( leapYear(year,l_360) ) days = iDayInYear + 1

  ENDIF

  END FUNCTION daysInYear
!#######################################################################
!#######################################################################
!#######################################################################
! function leapYear
! Internal procedure in module time_mod.
! Given the year, works out if it's a leap year.

  FUNCTION leapYear(year,l_360) RESULT(leap)

  IMPLICIT NONE
  LOGICAL :: leap  !  function result
  INTEGER, INTENT(in) :: year
  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used
!-------------------------------------------------------------------------------  
  leap = .FALSE.
  IF ( .NOT.l_360 .AND. MOD(year,4)==0 ) THEN
    leap = .TRUE.
!   End of century years (eg 1900,2000) are leap only if div by 400.
    IF ( MOD(year,100)==0 .AND. MOD(year,400)/=0 ) leap = .FALSE.
  ENDIF

  END FUNCTION leapYear
!#######################################################################
!#######################################################################
! function getHours
! Internal procedure in module time_mod.
! Calculate the number of hours between two times, each specified as
! seconds of day and date. dates may or may not be in chronological order.

  FUNCTION getHours( time1,date1,time2,date2,l_360,errMess ) &
             RESULT( nhr )

  USE timeConst, ONLY : iHourInDay,iSecInHour

  IMPLICIT NONE

  REAL :: nhr  !  function result

  INTEGER, INTENT(in) :: date1,date2,time1,time2

  INTEGER ::  &!  local SCALARS  
    i1,i2     &!  index of earlier and later date respectively  
   ,iyr

  INTEGER ::  &!  local ARRAYS  
    date(2),day(2),month(2),time(2),year(2)

  REAL :: hod(2)

  LOGICAL, INTENT(in) :: l_360  !  T if 360 day calendar used

  CHARACTER(len=*), INTENT(in) :: errMess  !  error message
!-----------------------------------------------------------------------

! Fill local copies of variables.
  date(1) = date1
  date(2) = date2
  time(1) = time1
  time(2) = time2

! Establish if times are in chronological order.
  i1 = 1  !  chronological order
  i2 = 2
  IF ( timeDate_cmp( time(1),date(1),'>',time(2),date(2),errMess ) ) THEN
    i1 = 2   !   reversed chronological order
    i2 = 1
  ENDIF

! Convert times from s to hours.
  hod(:) = REAL(time(:)) / REAL(iSecInHour)

! Decompose dates.
  CALL dateToBits(date(i1),day(i1),month(i1),year(i1),l_360,errMess)
  CALL dateToBits(date(i2),day(i2),month(i2),year(i2),l_360,errMess)

  IF ( year(i1) == year(i2) ) THEN

!   Dates are in same year, simply get difference.
    nhr = HourOfYear( hod(i2),day(i2),month(i2),year(i2),l_360 ) -  &
          HourOfYear( hod(i1),day(i1),month(i1),year(i1),l_360 )

  ELSE

!   From hod1 on date1 to end of first year.
    nhr = REAL( daysInYear(year(i1),l_360) )*REAL(iHourInDay) -  &
                hourOfYear( hod(i1),day(i1),month(i1),year(i1),l_360 )
!   Complete years other than 1st and last years.
    DO iyr=year(i1)+1,year(i2)-1
      nhr = nhr + REAL( daysInYear(iyr,l_360) ) * REAL(iHourInDay)
    ENDDO

!   Any time in final year.
    nhr  = nhr + hourOfYear( hod(i2),day(i2),month(i2),year(i2),l_360 )

  ENDIF

! If dates were not in chronological order, time difference is < 0.
  IF ( i1 == 2 ) nhr = -1.0 * nhr

  END FUNCTION getHours
!#######################################################################
!#######################################################################
! function timeDate_cmp
! Internal procedure in module time_mod.
! Compare two date/time pairs.

  FUNCTION timeDate_cmp( time1,date1,ineq,time2,date2,errMess ) &
                 RESULT( answer )

  IMPLICIT NONE

  LOGICAL :: answer  !  function result

  INTEGER, INTENT(in) ::  &!  IN scalars
     date1,date2    &!   dates (yyyymmdd)  
    ,time1,time2     !   times of day (s)

  CHARACTER(len=*), INTENT(in) :: &!  IN scalars
     errMess  &!    message printed on error  
    ,ineq      !    the inequality to be tested
!-------------------------------------------------------------------------------
  answer = .FALSE.

  SELECT CASE (ineq)
    CASE ( '<' )
      IF ( date1<date2 .OR. (date1==date2 .AND. time1<time2) ) answer=.TRUE.
    CASE ( '<=' )
      IF ( date1<date2 .OR. (date1==date2 .AND. time1<=time2) ) answer=.TRUE.
    CASE ( '=' )
      IF ( date1==date2 .AND. time1==time2 ) answer=.TRUE.
    CASE ( '>=' )
      IF ( date1>date2 .OR. (date1==date2 .AND. time1>=time2) ) answer=.TRUE.
    CASE ( '>' )
      IF ( date1>date2 .OR. (date1==date2 .AND. time1>time2) ) answer=.TRUE.
    CASE default
      WRITE(*,*)'ERROR: timeDate_cmp: do not recognise inequality.'
      WRITE(*,*) TRIM(errMess)
      STOP
  END SELECT

  END FUNCTION timeDate_cmp
!#######################################################################
!#######################################################################
!#######################################################################
!#######################################################################

! subroutine timeDate_diff
! Internal procedure in module time_mod.
! Calculate the time between two date/time pairs, each specified as
! seconds of day and date. dates may or may not be in chronological order.
! The answer may be given in terms of days/hours/mins/secs (dhms), or years,months,d/h/m/s.
! This is probably more complicated than it needs to be....

  SUBROUTINE timeDate_diff( time1,date1,time2,date2,l_360,errMess & 
                           ,nsec,nmin,nhr,ndy,nmon,nyr,chronoArg )

  USE timeConst, ONLY : iSecInDay,iSecInHour,iSecInMin

  IMPLICIT NONE

  INTEGER, INTENT(in) :: date1,date2,time1,time2

  INTEGER, INTENT(out) ::  &!  out SCALARS  
    nsec,nmin,nhr,ndy      !  number of seconds, minutes, hours, days

  INTEGER, INTENT(out), OPTIONAL ::  &!  out SCALARS  
    nmon,nyr   !  number of months, years

  INTEGER ::   &!  local SCALARS  
    i,i1,i2,tmpHour,tmpMin,tmpSec,tmpTime,yr  !  work

  INTEGER ::  &!  local ARRAYS  
    date(2),dateDone(2),day(2),hour(2),doy(2),mins(2),month(2),sec(2)  &!  work  
   ,time(2),year(2),tmpDate(2),tmpDay(3),tmpMonth(3),tmpYear(3)  !  work

  LOGICAL, INTENT(in) ::  &!  in SCALARS  
    l_360  !  T if 360 day calendar used

  LOGICAL, INTENT(out), OPTIONAL ::   &!  out SCALARS  
    chronoArg !  T means input times/dates were in chronological order (2nd is >= first)

  LOGICAL ::   &!  local SCALARS  
    chrono     &!  local version of chronoArg  
   ,dhms        !  T means answer is given in terms of days/hours/mins/secs
!                  F means answer is given in terms of years/months/days/hours/mins/secs  

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS  
     errMess   !  message printed on error
     
  LOGICAL :: l_neg_increment

!-----------------------------------------------------------------------

! Look for optional arguments for months and years.
  IF ( PRESENT(nmon) .NEQV. PRESENT(nyr) ) THEN
    WRITE(*,*)'timeDate_diff: need all or none of optional args nmon and nyr'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF
  dhms = .TRUE.
  IF ( PRESENT(nmon) ) dhms = .FALSE.

! Load into local variables.
  time(1) = time1
  time(2) = time2
  date(1) = date1
  date(2) = date2

! Establish if times are in chronological order.
  chrono = .TRUE.  !  chronological order
  i1 = 1
  i2 = 2
  IF ( timeDate_cmp( time(1),date(1),'>',time(2),date(2),errMess ) ) THEN
    chrono = .FALSE.
    i1 = 2   !   reversed chronological order
    i2 = 1
  ENDIF

! Decompose dates and get hours,mins,secs from secs in day.
  DO i=1,2
    CALL dateToBits( date(i),day(i),month(i),year(i),l_360,errMess )
    CALL secToSMH( time(i),sec(i),mins(i),hour(i),errMess )
  ENDDO

!-----------------------------------------------------------------------
  IF ( dhms ) THEN

!   Count days, hours, mins and secs.

    IF ( year(i1) == year(i2) ) THEN

!     Dates are in same year, simply get difference.
      nsec = secOfYear( sec(i2),mins(i2),hour(i2),day(i2),month(i2),year(i2),l_360 ) -  &
             secOfYear( sec(i1),mins(i1),hour(i1),day(i1),month(i1),year(i1),l_360 )
      ndy  = nsec / iSecInDay
      nsec = nsec - ndy*iSecInDay
      nhr  = nsec / iSecInHour
      nsec = nsec - nhr*iSecInHour
      nmin = nsec / iSecInMin
      nsec = nsec - nmin*iSecInMin

    ELSE

!     Note that tmpSec must be able to hold the number of s in 2 years,
!     which is ~63 million. This is well within the range of a 4-byte integer.

!     Count to end of first year.
      tmpSec = daysInYear(year(i1),l_360) * iSecInDay -   &
               secOfYear( sec(i1),mins(i1),hour(i1),day(i1),month(i1),year(i1),l_360 )

!     Get days in complete years other than 1st and last years.
      ndy = 0
      DO yr=year(i1)+1,year(i2)-1
        ndy = ndy + daysInYear(yr,l_360)
      ENDDO

!     Count time in final year.
      tmpSec  = tmpSec + SecOfYear( sec(i2),mins(i2),hour(i2),day(i2),month(i2),year(i2),l_360 )

!     Add contributions.
      tmpDay(1) = tmpSec/iSecInDay
      ndy       = ndy + tmpDay(1)
      tmpSec    = tmpSec - tmpDay(1)*iSecInDay
      tmpHour   = tmpSec/iSecInHour
      nhr       = tmpHour
      tmpSec    = tmpSec - tmpHour*iSecInHour
      tmpMin    = tmpSec/iSecInMin
      nmin      = tmpMin
      nsec      = tmpSec - tmpMin*iSecInMin
  
    ENDIF

!-----------------------------------------------------------------------
  ELSE   !  NOT dhms

!   Example: 06:30:00H 10/01/1985 to 04:30:00H 02/01/2003

!   Count complete years.
!   i.e. count from start time/date until same date in last possible complete year <= end.
!   e.g. 06:30:00 10/1/1985 to 06:30:00 10/1/2002
    nyr = year(i2) - year(i1)
!   If final time/date is earlier in year than start, final year was not complete.
!   e.g. 04:30:00H 02/01 is earlier in year than 06:30:00H 10/01.
!   Get date(i1) but in year(i2), eg 10/01/2003. Take care of Feb 29.
    tmpDay(1) = MIN( day(i1), dInMonth(month(i2),year(i2),l_360) )
    dateDone(1) = year(i2)*10000 + month(i1)*100 + tmpDay(1)
!   timeDate_cmp clause takes care of nyr=0 too.
    IF ( timeDate_cmp(time(i2),date(i2),'<',time(i1),dateDone(1),errMess) ) THEN
      dateDone(1) = MAX( dateDone(1)-10000  &!  1 year earlier
                        ,date(i1) )
      nyr = MAX( nyr-1, 0 )
    ENDIF
      

!   Count complete months (in incomplete years)
!   i.e. count up to start time on start day in last possible month <= end.
!   e.g. 06:30:00 10/1/2002 to 06:30:00 10/12/2002
!   Initially assume will count to start time of start day in end month.
!   Use end of month(i2) if month(i1) is longer.
    tmpDay(1) = MIN( day(i1), dInMonth(month(i2),year(i2),l_360) )
    datedone(2) = year(i2)*10000 + month(i2)*100 + tmpDay(1)  !  day(i1) in month(i2), year(i2)
    IF ( timeDate_cmp( time(i1),dateDone(2),'>',time(i2),date(i2),errMess ) ) THEN
!     Incomplete month. Use previous month.
      CALL timeDate( 0,dateDone(2),-1,'mon',l_360,tmpTime,tmpDate(1),errMess )
      dateDone(2) = tmpDate(1)
    ENDIF
    
    IF ( dateDone(2) >  date(i1) ) THEN
!     Get months from dates.
      CALL dateToBits( dateDone(1),tmpDay(1),tmpMonth(1),tmpYear(1),l_360,errMess )
      CALL dateToBits( dateDone(2),tmpDay(2),tmpMonth(2),tmpYear(2),l_360,errMess )
      nmon = tmpMonth(2) - tmpMonth(1)
      IF ( nmon < 0 ) nmon = 12 - tmpMonth(1) + tmpMonth(2)
      dateDone(1) = dateDone(2)
    ELSE
      nmon = 0
!      dateDone(1) = dateDone(2) !  xx what happens here??
    ENDIF

!   Count days, to start time on last day (which may be beyond end time).
!   e.g. 06:30:00 10/12/2002 to 06:30:00 2/1/2003
!   Work out days of year.
    CALL dateToBits( dateDone(1),tmpDay(1),tmpMonth(1),tmpYear(1),l_360,errMess )
    doy(1) = dayOfYear( tmpDay(1),tmpMonth(1),tmpYear(1),l_360 )
    doy(2) = dayOfYear( day(i2),month(i2),year(i2),l_360 )
    ndy = doy(2) - doy(1)
    IF ( ndy < 0 ) ndy = daysInYear( tmpYear(2),l_360 ) - doy(1) + doy(2)

!   Count time on last day.
!   e.g. 06:30:00 2/1/2003 to 04:30:00 2/1/2003
    tmpTime = time(i2) - time(i1)
    IF ( tmpTime < 0 ) THEN
      ndy = ndy - 1
      tmpTime = iSecInDay + tmpTime
    ENDIF
    CALL secToSMH( tmpTime,nsec,nmin,nhr,errMess )

  ENDIF  !  dhms
!-----------------------------------------------------------------------
! Debugging - stop if any increments are less than zero.

! Hack to get round debugging stuff causing the program to execute
! both sides of an .AND. even when left hand side is false
  l_neg_increment = .NOT. dhms
  IF ( l_neg_increment ) l_neg_increment = (nmon<0 .OR. nyr<0)
  
!  IF ( nsec<0 .OR. nmin<0 .OR. nhr<0 .OR. ndy<0 .OR.  &
!       ( .NOT.dhms .AND.(nmon<0 .OR. nyr<0) ) ) THEN
  IF ( nsec<0 .OR. nmin<0 .OR. nhr<0 .OR. ndy<0 .OR. l_neg_increment ) THEN
    WRITE(*,*)'ERROR: timeDate_diff: answer<0!'
    WRITE(*,*) TRIM(errMess)
    WRITE(*,*)'Inputs: time1,date1=',time1,date1,' time2,date2=',time2,date2
    WRITE(*,*)'nsec=',nsec,' nmin=',nmin,' nhr=',nhr,' ndy=',ndy
    IF ( .NOT. dhms ) WRITE(*,*)'nmon=',nmon,' nyr=',nyr
    STOP
  ENDIF
!-----------------------------------------------------------------------

! If dates were not in chronological order, time difference is < 0.
  IF ( .NOT. chrono ) THEN
    nsec = -1 * nsec
    nmin = -1 * nmin
    nhr  = -1 * nhr
    ndy  = -1 * ndy
    IF ( .NOT. dhms ) THEN
      nmon = -1 * nmon
      nyr  = -1 * nyr
    ENDIF
  ENDIF

  IF ( PRESENT(chronoArg) ) chronoArg=chrono

  END SUBROUTINE timeDate_diff
!#######################################################################
!#######################################################################
! subroutine secTosmh
! Given a time (seconds), returns secs,mins and hours.
! One use is to convert time of day (s) to s/m/h.

  SUBROUTINE secToSMH ( time,s,m,h,errMess )

  USE timeConst, ONLY : iSecInDay,iSecInHour,iSecInMin

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS  
      time   !  time of day (s)

  INTEGER, INTENT(out) ::  &!  out SCALARS  
      s,m,h   !  s,mins,hours

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS  
      errMess   !  message printed on error
!-------------------------------------------------------------------------------
  IF ( time<0 .OR. time>iSecInDay ) THEN
    WRITE(*,*)'ERROR: secToSMH time must be in range 0 to 86400s'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

  h = time / iSecInHour
  m = ( time  - h*iSecInHour ) / iSecInMin
  s = time - h*iSecInHour - m*iSecInMin

  END SUBROUTINE secToSMH
!#######################################################################
!#######################################################################
! function w
! Internal procedure in module time_mod.
! Given a time of day (seconds), returns character "hh:mm:ss H", or
! "hhmmss" according to flag.

  FUNCTION s_to_chhmmss ( time,shortArg ) RESULT( hms )

  USE timeConst, ONLY : iSecInDay,iSecInHour,iSecInMin

  IMPLICIT NONE

  CHARACTER(len=10) ::  hms  !  scalar function result of form "hh:mm:ss H"
!    If optional argument is set to TRUE, the result is "hhmmss".
    

  INTEGER, INTENT(in) ::  &!  in scalars
    time   !  time of day (s)

  LOGICAL, INTENT(in), OPTIONAL ::  &
    shortArg  !  flag controlling form of output
!                  If present and TRUE, output is of form "hhmmss".
!                  Else output is of form "hh:mm:ss H".

  INTEGER  ::  &!  local scalars
    s,m,h   &!  s,mins,hours
   ,t        !  time of day (s)

  LOGICAL :: shortForm  !  local version of shortArg (qv)

!-------------------------------------------------------------------------------

! Use optional argument if present, otherwise set flag to FALSE.
  IF ( PRESENT(shortArg) ) THEN
    shortForm = shortArg
  ELSE
    shortForm = .FALSE.
  ENDIF

  t = MOD( time, iSecInDay )

  h = t / iSecInHour
  m = ( t  - h*iSecInHour ) / iSecInMin
  s = t - h*iSecInHour - m*iSecInMin

  IF ( shortForm ) THEN
    WRITE(hms,"(3(i2.2))") h,m,s
  ELSE
    WRITE(hms,"(i2.2,a1,i2.2,a1,i2.2,a2)") h,':',m,':',s,' H'
  ENDIF

  END FUNCTION s_to_chhmmss

!#######################################################################
!#######################################################################
! function chhmmss_to_s
! Given a character time of day ( hh:mm:ss ), returns seconds of day.
! Input should be in range 00:00:00 to 23:59:59.

  FUNCTION chhmmss_to_s( charTime,errMess ) RESULT( time )

  USE timeConst, ONLY : iSecInHour,iSecInMin

  IMPLICIT NONE

  INTEGER ::  &!  SCALAR function result  
     time      !  time of day (s)

  INTEGER ::  &!  local ARRAYS  
     tval(3)   !  work

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS  
     charTime  &!  time (hh:mm:ss)  
    ,errMess    !  message printed on error

!-----------------------------------------------------------------------

  IF ( LEN_TRIM(charTime) < 8 ) THEN
    WRITE(*,*)'ERROR: chhmmss_to_s: input too short.'
    WRITE(*,*) 'Format for time should be hh:mm:ss.'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF
  IF ( charTime(3:3)/=':' .OR. charTime(6:6)/=':' ) THEN
    WRITE(*,*)'ERROR: chhmmss_to_s: format'
    WRITE(*,*)'Format for time should be hh:mm:ss.'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

  READ(charTime(1:2),"(i2)") tval(1)
  READ(charTime(4:5),"(i2)") tval(2)
  READ(charTime(7:8),"(i2)") tval(3)

  IF ( tval(1)>23 .OR. tval(2)>59 .OR. tval(3)>59 ) THEN
    WRITE(*,*)'ERROR: chhmmss_to_s: tval out of range.'
    WRITE(*,*) TRIM(errMess)
    STOP
  ENDIF

  time = tval(1)*iSecInHour + tval(2)*iSecInMin + tval(3)

  END FUNCTION chhmmss_to_s
!#######################################################################
!#######################################################################
! subroutine get_period
! Given a time period (seconds), suggest an appropriate unit to
! describe it and the period in those units, eg 7200s can be described
! in terms of hours.

  SUBROUTINE get_period( period,timeStep,nper,units,longUnits )

  USE inout, ONLY :  &
!  imported scalar parameters
     periodAnn,periodMon

  USE timeConst, ONLY : &
     iSecInDay,iSecInHour,iSecInMin

  IMPLICIT NONE

!--------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!--------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: period    !  time period
!                                     >0 is period as a number of timesteps
!                                    <=0 is a "special" case
  INTEGER, INTENT(in) :: timeStep  !  timestep length (s) for period>0

!--------------------------------------------------------------------------------
! Scalar arguments with intent(out)
!--------------------------------------------------------------------------------
  INTEGER, INTENT(out) :: nper             !  period converted to the units of
!                                             "units"
  CHARACTER(len=2), INTENT(out) :: units   !  units for period

!--------------------------------------------------------------------------------
! Optional scalar arguments with intent(out)
!--------------------------------------------------------------------------------
  CHARACTER(len=7), INTENT(out), OPTIONAL  :: longUnits  !  a longer version of
!                                                           units

!--------------------------------------------------------------------------------

  IF ( period < 0 ) THEN

    IF ( period == periodAnn ) THEN

!     Period is one year.
      units = 'yr'
      nper = 1
      IF ( PRESENT(longUnits) ) longUnits = 'years'

    ELSEIF ( period == periodMon ) THEN

!     Period is one month.
      units = 'mo'
      nper = 1
      IF ( PRESENT(longUnits) ) longUnits = 'months'

    ELSE

     WRITE(*,*)'ERROR: get_period: do not recognise period=',period
     STOP

    ENDIF

  ELSEIF ( MOD( period*timeStep, iSecInDay ) == 0 ) THEN

!   Period is a number of days.
    units = 'dy'
    nper = period * timeStep / iSecInDay
    IF ( PRESENT(longUnits) ) longUnits = 'days'

  ELSEIF ( MOD( period*timeStep, iSecInHour ) == 0 ) THEN
!   Period is a number of hours
    units = 'hr'
    nper = period * timeStep / iSecInHour
    IF ( PRESENT(longUnits) ) longUnits = 'hours'

  ELSEIF ( MOD( period*timeStep, iSecInMin ) == 0 ) THEN

!   Period is a number of minutes
    units = 'mn'
    nper = period * timeStep / iSecInMin
    IF ( PRESENT(longUnits) ) longUnits = 'minutes'

  ELSE

!   Default to using seconds to describe period.
    units = 's'
    nper = period
    IF ( PRESENT(longUnits) ) longUnits = 'seconds'

  ENDIF

  END SUBROUTINE get_period
!###############################################################################
!###############################################################################
!###############################################################################
  END MODULE time_mod
