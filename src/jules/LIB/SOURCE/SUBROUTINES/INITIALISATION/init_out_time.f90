!###############################################################################
!###############################################################################
! subroutine init_out_time
! Read and process the dates and times chosen for output.
!###############################################################################
!###############################################################################

  SUBROUTINE init_out_time( iout )

  USE inout, ONLY :  &
!   imported scalar parameters
      formatBin,periodAnn,periodMon,periodOneFile  &
!   imported scalars with intent(in)
     ,outGrads,jinUnit  &
!   imported arrays with intent(inout)
     ,outDate,outDateFlag,outFilePer,outPer  &
     ,outSamPer,outTime,outTemplate

  USE spin_mod, ONLY :  &
!   imported scalars with intent(in)
     nspin,spinUp

  USE time_loc, ONLY :  &
!   imported scalars with intent(in)
     timeStep  &
!   imported arrays with intent(in)
    ,dateMainRun,dateSpin,timeRun

  USE time_mod, ONLY :  &
!   imported procedures
     chhmmss_to_s,get_period,timeDate_cmp

  USE timeConst, ONLY : &
     iSecInDay,iSecInHour,iSecInMin

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout   !  number of output profile

  INTEGER ::  &!  local SCALARS
    i,iper   !  work

  CHARACTER(len=2) ::  &!  local SCALARS
    periodUnits   !  suggested time unit

  CHARACTER(len=8) ::  &!  local SCALARS
    cOutTime(2)  !   times (hh:mm:ss)


!------------------------------------------------------------------------
! Read output periods, dates etc.
!------------------------------------------------------------------------
  READ(jinUnit,*) outPer(iout),outFilePer(iout)
  READ(jinUnit,*) outSamPer(iout)
  READ(jinUnit,*) outDate(iout,1),cOutTime(1)
  READ(jinUnit,*) outDate(iout,2),cOutTime(2)

!------------------------------------------------------------------------
! Convert input times (hh:mm:ss) to seconds (even if outDate is a special
! case).
!------------------------------------------------------------------------
  DO i=1,2
    outTime(iout,i) = chhmmss_to_s( coutTime(i),'init_out_time' )
  ENDDO

!------------------------------------------------------------------------
! If only the initial state is to be output, set various periods to make
! life easier.
!------------------------------------------------------------------------
  IF ( outDate(iout,1) == -2 ) THEN
    outPer(iout) = NINT( timestep )
    outSamPer(iout) = outPer(iout)
    outFilePer(iout) = periodOneFile
  ENDIF

!------------------------------------------------------------------------
! Check dates and/or interpret special cases. Also set flag outDateFlag.
! First, if there is no spin up, don't use the special case of outDate=-1.
!------------------------------------------------------------------------
  IF ( nspin==0 .AND. outDate(iout,1)==-1 ) outDate(iout,1)=0

  SELECT CASE ( outDate(iout,1) )
!------------------------------------------------------------------------
    CASE ( 1: )
!     outDate>=1 - these are the dates to use.
      outDateFlag(iout) = 1
!------------------------------------------------------------------------
!     If dates/times lie outside those of the main run, set to the start/
!     end of the main run. Probably this should be treated as a fatal
!     error...but we'll allow for now....better too much output than none.
!------------------------------------------------------------------------
      IF ( timeDate_cmp( outTime(iout,1),outDate(iout,1),'<'  &
                        ,timeRun(1),dateMainRun(1),'init_out_time' ) .OR.  &
           timeDate_cmp( outTime(iout,1),outDate(iout,1),'>'  &
                        ,timeRun(2),dateMainRun(2),'init_out_time' ) ) THEN
        WRITE(*,*)'init_out_time: moving start time of output to fall within main run. iout=',iout
        outTime(iout,1) = timeRun(1)
        outDate(iout,1) = dateMainRun(1)
      ENDIF
      IF ( timeDate_cmp( outTime(iout,2),outDate(iout,2),'<'  &
                        ,timeRun(1),dateMainRun(1),'init_out_time' ) .OR.  &
           timeDate_cmp( outTime(iout,2),outDate(iout,2),'>'  &
          ,timeRun(2),dateMainRun(2),'init_out_time' ) ) THEN
        WRITE(*,*)'init_out_time: moving end time of output to fall within main run. iout=',iout
        outTime(iout,2) = timeRun(2)
        outDate(iout,2) = dateMainRun(2)
      ENDIF
!------------------------------------------------------------------------
!     Check for dates/times that are not in chronological order.
!------------------------------------------------------------------------
      IF ( timeDate_cmp(outTime(iout,1),outDate(iout,1),'>',outTime(iout,2),outDate(iout,2),'init_out_time' ) ) THEN
        WRITE(*,*)'ERROR: dates/times for output are incorrect (backwards).'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
!------------------------------------------------------------------------
    CASE ( 0 )
!     Output at all times, including during any spin up.
      outDateFlag(iout) = 0
!     Set dates and times when output starts. End date/time are not needed, but set anyway.
      outDate(iout,1:2) = (/ dateMainRun(1), dateMainRun(2) /)
      outTime(iout,1:2) = (/ timeRun(1), timeRun(2) /)
      IF ( nspin /= 0 ) outDate(iout,1) = dateSpin(1)
!------------------------------------------------------------------------
    CASE ( -1 )
!     Output at all times after spin up.
      outDateFlag(iout) = -1
!     Set dates and times when output starts. End date/time are not needed, but set anyway.
      outDate(iout,1:2) = (/ dateMainRun(1), dateMainRun(2) /)
      outTime(iout,1:2) = (/ timeRun(1), timeRun(2) /)
!------------------------------------------------------------------------
    CASE ( -2 )
!     Output initial state (at start of first timestep).
      outDateFlag(iout) = -2
!     Set dates and times when output starts. End date/time are not needed, but set anyway.
      outDate(iout,1:2) = (/ dateMainRun(1), dateMainRun(1) /)
      outTime(iout,1:2) = (/ timeRun(1), timeRun(1) /)
      IF ( nspin /= 0 ) outDate(iout,1) = dateSpin(1)
!------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'ERROR: init_out_time: do not understand outDate.'
      WRITE(*,*)'iout=',iout,' outDate(iout,1)=',outDate(iout,1)
  END SELECT

!------------------------------------------------------------------------
! Check periods are OK and convert from s to # of timesteps.
!------------------------------------------------------------------------
  IF ( outPer(iout) >= 0 ) THEN

    IF ( outPer(iout) == 0 ) THEN
      outPer(iout) = 1
    ELSE
!     Impose a limit of 30 days - makes it easier to initialise outStep.
      IF ( outPer(iout) > 30*iSecInDay ) THEN
        WRITE(*,*)'ERROR: init_out_time: constant output period must be <= 30 days.'
        WRITE(*,*)'Profile #',iout,' period=',REAL(outper(iout))/REAL(iSecInDay),' days.'
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
      IF ( MOD( outPer(iout),NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: period for output must be a multiple of timestep'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ELSE
        outPer(iout) = outPer(iout) / NINT(timeStep)
      ENDIF
    ENDIF

!------------------------------------------------------------------------
!   Insist that periods can be synchronised (kept in phase with) with days.
!   Periods <= 1 day must fit into 1 day, periods > 1 day must be
!   multiples of a day.  This is not essential to the rest of the code
!   (at present!!), but seems sensible.
!------------------------------------------------------------------------
    IF ( outPer(iout)*NINT(timeStep) <= iSecInDay ) THEN
      IF ( MOD( iSecInDay,outPer(iout)*NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: one day must be a multiple of period for output (for period <= 1day)'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
    ELSE
      IF ( MOD( outPer(iout)*NINT(timeStep),iSecInDay ) /= 0 ) THEN
        WRITE(*,*)'ERROR: periods longer than one day must be a multiple of days.'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
    ENDIF

  ELSE  !  outPer<0

    SELECT CASE( outper(iout) )
      CASE ( periodMon, periodAnn )  !  Allowable values.
      CASE default
        WRITE(*,*)'ERROR: do not recognise value (',outper(iout),') for output period.'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
    END SELECT

  ENDIF  !  outper

!------------------------------------------------------------------------
! Check file periods are OK and convert from s to # of timesteps.
!------------------------------------------------------------------------
  IF ( outFilePer(iout) >= 0 ) THEN

    IF ( outFilePer(iout) == 0 ) THEN
      outFilePer(iout) = 1
    ELSEIF( outFilePer(iout) > 0 ) THEN
!     Impose a limit of 30 days - makes it easier to initialise outFileStep.
      IF ( outFilePer(iout) > 30*iSecInDay ) THEN
        WRITE(*,*)'ERROR: init_out_time: constant output file period must be <= 30 days.'
        WRITE(*,*)'Profile #',iout,' file period=',REAL(outFilePer(iout))/REAL(iSecInDay),' days.'
        STOP
      ENDIF
      IF ( MOD( outFilePer(iout),NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: period for output files must be a multiple of timestep'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
      outFilePer(iout) = outFilePer(iout) / NINT(timeStep)
    ENDIF

!------------------------------------------------------------------------
!   Insist that file periods can be synchronised (kept in phase with)
!   with days.  Periods <= 1 day must fit into 1 day, periods > 1 day
!   must be multiples of a day.  This is not essential to the rest of the
!   code (at present!!), but seems sensible.
!------------------------------------------------------------------------
    IF ( outFilePer(iout)*NINT(timeStep) <= iSecInDay ) THEN
      IF ( MOD( iSecInDay,outFilePer(iout)*NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: one day must be a multiple of period for output files (for period <= 1day)'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
    ELSE
      IF ( MOD( outFilePer(iout)*NINT(timeStep),iSecInDay ) /= 0 ) THEN
        WRITE(*,*)'ERROR: file periods longer than one day must be a multiple of days.'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out_time'
        STOP
      ENDIF
    ENDIF

  ELSE  !  outFilePer<0

    SELECT CASE( outFilePer(iout) )
      CASE ( periodMon, periodAnn, periodOneFile )
!       Special cases, nothing to do.
      CASE ( -7, -8 )
!       Special cases, only for use with spin up. If no spin up, change.
        IF ( .NOT. spinUp ) outFilePer(iout) = periodOneFile
      CASE default
         WRITE(*,*)'ERROR: do not recognise value for output file period.'
         WRITE(*,*)'Error for output profile #',iout
         WRITE(*,*)'Stopping in init_out_time'
         STOP
    END SELECT
  ENDIF

!------------------------------------------------------------------------
! Check that output and file periods are compatible.  These checks are
! not exhaustive. For example, output every 50 days is possible with
! monthly files - but many files will be empty! (NB Haven't actually
! tested this!)  More checking is done later for GrADS, particularly
! template option.
!------------------------------------------------------------------------
  IF ( outFilePer(iout) >= 0 ) THEN

!   NB This test will stop for outPer > outFilePer, as required.
    IF ( MOD( outFilePer(iout),outPer(iout) ) /= 0 ) THEN
      WRITE(*,*)'ERROR: period for output files must be a multiple of output period'
      WRITE(*,*)'Error for output profile #',iout
      WRITE(*,*)'Stopping in init_out_time'
      STOP
    ENDIF

  ELSE  !  outFilePer < 0

!------------------------------------------------------------------------
!   Not very restrictive testing, so plenty of (weird) combinations left.
!   eg monthly files test for annual output, but outper=3456000 (40days)
!   would be allowed.
!------------------------------------------------------------------------
    SELECT CASE( outFilePer(iout) )
      CASE ( periodMon )
!       Monthly output files. Just check we don't have annual output.
        IF ( outPer(iout) == periodAnn ) THEN
          WRITE(*,*)'ERROR: monthly output files, cannot use annual output.'
          WRITE(*,*)'Error for output profile #',iout
          WRITE(*,*)'Stopping in init_out_time'
          STOP
        END IF
    END SELECT

  ENDIF  !  outFilePer

!------------------------------------------------------------------------
! Further checks for output that is to be used with GrADS,
!------------------------------------------------------------------------
  IF ( outGrADS ) THEN

!------------------------------------------------------------------------
!   Check output and file periods are acceptable (i.e. can be represented).
!------------------------------------------------------------------------
    IF ( outPer(iout) > 0 ) THEN
      call get_period( outPer(iout),NINT(timeStep),iper,periodUnits )
      SELECT CASE ( periodUnits )
        CASE ( 'yr','mo','dy','hr','mn' )   !  OK
        CASE default
          WRITE(*,*)'ERROR: output period is not acceptable for GrADS output.'
          WRITE(*,*)'Error for output profile #',iout
          WRITE(*,*)'Stopping in init_out_time'
          STOP
      END SELECT
    ENDIF

    IF ( outFilePer(iout) > 0 ) THEN
      call get_period( outFilePer(iout),NINT(timeStep),iper,periodUnits )
      SELECT CASE ( periodUnits )
        CASE ( 'yr','mo','dy','hr','mn' )   !  OK
        CASE default
          WRITE(*,*)'ERROR: output file period is not acceptable for GrADS output.'
          WRITE(*,*)'Error for output profile #',iout
          WRITE(*,*)'Stopping in init_out_time'
          STOP
      END SELECT
    ENDIF

!------------------------------------------------------------------------
!   Check the use of templating.
!------------------------------------------------------------------------
    IF ( outTemplate(iout) ) THEN

!     No need for templating in certain cases.
      SELECT CASE( outFilePer(iout) )
        CASE ( periodOneFile, -8, -7 )
!         Template output does not make sense for these options.
          outTemplate(iout) = .FALSE.
      END SELECT
!     No need for template if there is (obviously) only one output file.
      IF ( outDateFlag(iout) == -2 ) outTemplate(iout)=.FALSE.

!------------------------------------------------------------------------
!     Templating will not work if output and file periods are incompatible.
!     The following is my present understanding of this, but may be wrong.
!     Almost definitely incomplete - all possible ways of specifying
!     outper and outFilePer probably haven't been tested.  Given that the
!     output file period is described using time units U, and the output
!     file period is n*U, if n=1 any output period <= 1U is allowed,
!     if n>1 the output period must be n*U.  Nothing to test for the
!     "simple" cases of monthly or annual files, since these are 1 x time
!      unit.  We have previously checked that outPer <= outFilePer.
!------------------------------------------------------------------------

      IF ( outFilePer(iout) > 0 ) THEN

!       Get period units again, in case code above is ever changed!
        call get_period( outFilePer(iout),NINT(timeStep),iper,periodUnits )

        IF ( iper>1 .AND. outPer(iout)/=outFilePer(iout) ) THEN
          WRITE(*,*)'ERROR: output period and file period are not suitable for GrADS template option'
          WRITE(*,*)'Error for output profile #',iout
          WRITE(*,*)'Stopping in init_out_time'
          STOP
        ENDIF

      ENDIF  !  outFilePer

    ENDIF  !  outTemplate

  ELSE !  outFormat

!------------------------------------------------------------------------
!   Check output period is acceptable - only needed because have only
!   written code to deal with GrADS cases.
!------------------------------------------------------------------------
    IF ( outPer(iout) > 0 ) THEN
      call get_period( outPer(iout),NINT(timeStep),iper,periodUnits )
      SELECT CASE ( periodUnits )
        CASE ( 'yr','mo','dy','hr','mn' )   !  OK
        CASE default
          WRITE(*,*)'ERROR: output period is not acceptable - no code!'
          WRITE(*,*)'Error for output profile #',iout
          WRITE(*,*)'Stopping in init_out_time'
          STOP
      END SELECT
    ENDIF

  ENDIF  !  outFormat

  END SUBROUTINE init_out_time
!###############################################################################
!###############################################################################
