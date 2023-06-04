! Read details of model timestep and run length.

  SUBROUTINE init_time(dtlong,iyear1,imonth1,idate1,itime1,timmax, time_BRAMS,frqhis_BRAMS,runtype)

  USE file_utils, ONLY: findTag
  USE inout, ONLY : echo,jinUnit

  USE misc_utils, ONLY :  &
!  imported procedures
     read_list,repeatVal,varList

  USE spin_mod, ONLY :  &
!  imported scalar parameters
     nspinVarMax  &
!  imported scalars with intent(out)
    ,iposRouteStore,iposSmcl,iposTsoil,ispin,nspin,nspinFinal,nspinVar,spinEnd,spinFail,spinUp  &
!  imported arrays with intent(out)
    ,spinTol,spinTolPercent,spinVar,spinVarName

  USE switches, ONLY : l_360,route,routeOnly
  USE time_loc, ONLY : date,dateMainRun,dateNext,datePrev,dateRun,dateSpin,endMonth,endSec,endYear  &
                      ,stepFlag,time,timeNext,timePrev,timeRun,timeStep
  USE time_mod, ONLY : chhmmss_to_s,s_to_chhmmss,timeDate,timeDate_cmp
  USE trifctl, ONLY : phenol_period,triffid_period
  USE timeConst, ONLY : iSecInDay

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    i,ivar     &!  loop counters/work
   ,nvarFound   !  number of variables found in list

  LOGICAL ::  endrun   !  work (newTime needs an intent(out) logical variable,
!                         will be FALSE here)
  LOGICAL :: errFound  !  error flag

  CHARACTER(len=8) :: cTimeRun(2)  !   times (hh:mm:ss)

  CHARACTER(len=10) ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
!                                 ! that is required because some compilers can not cope
!                                 ! with recursive write statements

  real,              intent(in)  :: dtlong,timmax, time_BRAMS,frqhis_BRAMS
  integer,           intent(in)  :: iyear1,imonth1,idate1,itime1
  CHARACTER(len=16), intent(in)  :: runtype
  
  integer                           :: timmax_seg, itime1_seg

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  REAL :: spinTolTmp(nspinVarMax)       !  work
  LOGICAL :: spinTolPcTmp(nspinVarMax)  !  work
  CHARACTER(len=LEN(spinVarName)) :: varNameTmp(nspinVarMax)
  CHARACTER(len=50) :: spinVarDesc(nspinVarMax)  !  description of each spin up variable

!------------------------------------------------------------------------
! Print message to screen.
!------------------------------------------------------------------------
  if (echo) WRITE(*,"(50('-'),/,a)") 'init_time'

!------------------------------------------------------------------------
! Locate the start of this section in input file.
!------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_time','>INIT_TIME' )

!------------------------------------------------------------------------
! Initialise.
!------------------------------------------------------------------------
  spinUp = .FALSE.
  spinEnd = .FALSE.
  spinVar(:) = .FALSE.  !  no variables used to determine spin up
  stepFlag = 3  !  first timestep in "main" run

!------------------------------------------------------------------------
! Read data from run control file.
!------------------------------------------------------------------------
  READ(jinUnit,*) timeStep
  READ(jinUnit,*) dateMainRun(1),ctimeRun(1)
  READ(jinUnit,*) dateMainRun(2),cTimeRun(2)

  READ(jinUnit,*) l_360
  READ(jinUnit,*) phenol_period,triffid_period

  READ(jinUnit,*) dateSpin(1:2),nspin

!------------------------------------------------------------------------
! Only read spin-up details if there is to be spin up.
!------------------------------------------------------------------------
  IF ( nspin /= 0 ) THEN
    READ(jinUnit,*) spinFail

!   Establish which variables are to be used to determine spin up.
!   Read variable name, a logical flag and a constant from a blank-delimited list.
    CALL read_list( jinUnit,3,nspinvarMax,'>VARS','>ENDVARS',' ','init_time'  &
                     ,nvarFound,cvar1=varNameTmp,cvar1Pos=1  &
                     ,lvar1=spinTolPcTmp,lvar1Pos=2   &
                     ,rvar=spinTolTmp,rvarPos=3 )
  ENDIF

  !{DSM
      if (nspin/=0) then
         print*
         print*, 'nspin deve ser igual a 0 (zero) nesta versao acoplada!!!'
         stop
      endif
      timeStep=dtlong  !DSM foi colocado em sfclyr_jules.f90 um fator multiplicando o dtlong para nao chamar o JULES a cada timeStep do BRAMS
      dateMainRun(1)=iyear1*10000+imonth1*100+idate1
      write(ctimeRun(1)(1:2),'(i2.2)') itime1/100
      ctimeRun(1)(3:3)=':'
      write(ctimeRun(1)(4:5),'(i2.2)') mod(itime1,100)
      ctimeRun(1)(6:8)=':00'

      timmax_seg=nint(timmax)
      itime1_seg=(itime1/100)*3600 + mod(itime1,100)*60
  
      CALL timeDate( itime1_seg,dateMainRun(1),timmax_seg,'sec',l_360,timmax_seg,dateMainRun(2),'init_time')
           
      write(ctimeRun(2)(1:2),'(i2.2)') timmax_seg/3600
      ctimeRun(1)(3:3)=':'
      write(ctimeRun(2)(4:5),'(i2.2)') mod(timmax_seg,3600)/60
      ctimeRun(2)(6:6)=':'
      write(ctimeRun(2)(7:8),'(i2.2)') mod( mod(timmax_seg,3600) ,60)
  !DSM}

!------------------------------------------------------------------------
! Nothing else read by this subroutine.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! Process time variables. Check that timestep is an integer number of
! seconds.
!------------------------------------------------------------------------
  IF ( timeStep - REAL(INT(timeStep)) >= EPSILON(timeStep) ) THEN
    WRITE(*,*)'ERROR: INIT_TIME: timestep must be an integer number of seconds'
    STOP
  ENDIF

!------------------------------------------------------------------------
! Make sure one day is a multiple of timestep - makes life easier.
!------------------------------------------------------------------------
  IF ( MOD(iSecInDay,NINT(timeStep)) /= 0 ) THEN
    WRITE(*,*)'ERROR: init_time: 24 hours must be a multiple of timestep'
    STOP
  ENDIF

!------------------------------------------------------------------------
! Convert input times (hh:mm:ss) to seconds.
!------------------------------------------------------------------------
  DO i=1,2
    timeRun(i) = chhmmss_to_s( cTimeRun(i),'init_time' )
  ENDDO

!------------------------------------------------------------------------
! At present, any spin up must be over times that either:
!    1. start at the same the "main" run, and end at or before the time
!       when the main run starts.
! or 2. end at the start time of the "main" run.
!
! For #1, the spin up is thus over times that are either the same as
! those for the main run, or are a subset of those starting at same time.
! #1 means that spin up period is < or = to length of main run. For #2
! the spin up can be over a period that is <,=, or > than main run.
!------------------------------------------------------------------------

! Check that run times are in chronological order.
  IF ( timeDate_cmp(timeRun(1),dateMainRun(1),'>=',timeRun(2),dateMainRun(2),'init_time') ) THEN
    WRITE(*,*)'ERROR: init_time: Start time/date for main run is >= end.'
    STOP
  ENDIF

! Set dates for integration, assuming no spin up.
  dateRun(:) = dateMainRun(:)

! Initially assume that all spin-up cycles will be required.
  nspinFinal = nspin

! Check any spinup period is correctly prescribed.
  IF ( nspin /= 0 ) THEN
    spinUp=.TRUE.
    ispin = 0   !  the call to newTime increments this

!------------------------------------------------------------------------
!   Check dates for spin up are reasonable.  Spinup must immediately
!   precede the main run, or have same start date as run (eg if want to
!   spinup over same year). Spinup will be for a complete number of days.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
!   First checking that spinup dates are successive.
!------------------------------------------------------------------------
    IF ( timeDate_cmp( 0,dateSpin(1),'>=',0,dateSpin(2),'init_time' ) ) THEN
      WRITE(*,*)'ERROR: End date for spin up is <= start.'
      WRITE(*,*)'Stopping in init_time'
      STOP
    ENDIF
!   Now check that spinup times agree with run times...
    IF ( dateSpin(2)/=dateMainRun(1) .AND. dateSpin(1)/=dateMainRun(1) ) THEN
        WRITE(*,*)'ERROR: Spinup must immediately precede main run OR start at same date/time as main run.'
      WRITE(*,*)'Stopping in init_time'
      STOP
    ENDIF

!------------------------------------------------------------------------
!   Reset start date for integration.
!------------------------------------------------------------------------
    dateRun(1) = dateSpin(1)

!------------------------------------------------------------------------
!XX At present, output does not work particularly well for the case of a
!XX single cycle of spin up followed by a main run with
!XX dateMainRun(1)=dateSpin(2). Part of the problem is that endSec is set
!XX to T at end of spin up, whereas other bits of code want endSec=F for
!XX this special case.   This needs to be sorted out, possibly by always
!XX setting endSec=T at end of spin up, and always naming spin-up output
!XX with .spin. in filenames. Any solution should probably avoid creating
!XX this "special case", since it has to be correctly accounted for at
!XX all later stages in the code!  At present, if spin up ends mid-month,
!XX monthly output creates a file for the last month of spin up, then
!XX endSec resets outFirstWrite, and then at start of main run a file
!XX with same name is required...and goes wrong.  So, avoid all these
!XX complications for now (although current code is likely OK for cases
!XX such as no output through spin up, or every time to separate file).
!------------------------------------------------------------------------
    IF ( nspin==1 .AND. dateMainRun(1)==dateSpin(2) ) THEN
      WRITE(*,*)'ERROR: init_time: precautionary error!'
      WRITE(*,*)'Output code is not robust/working correctly for nspin==1 .AND. dateMainRun(1)==dateSpin(2).'
      WRITE(*,*)'You can anyway do this integration by setting'
      WRITE(*,*)'dateMainRun(1) to current value of =dateSpin(1), and setting nspin=0'
      WRITE(*,*)'This means the "spin up" part of the run is now considered part of the "main" run.'
      WRITE(*,*)'Sorry about this - hopefully the code will be improved, one day...'
      STOP
    ENDIF

  ENDIF  !  nspin

!------------------------------------------------------------------------
! Set time to one timestep before start, so that first increment (call to
! newTime, below) takes to start time.
!------------------------------------------------------------------------
  IF ( .NOT. spinUp ) THEN
    date = dateRun(1)
  ELSE
    date = dateSpin(1)
  ENDIF
  time = timeRun(1)
  CALL timeDate( time,date,-1*NINT(timeStep),'sec',l_360,timePrev,datePrev,'init_time')
  date = datePrev
  time = timePrev

! Initialise past and next times.
  CALL timeDate( time,date,-1*NINT(timeStep),'sec',l_360,timePrev,datePrev,'init_time' )
  CALL timeDate( time,date,NINT(timeStep),'sec',l_360,timeNext,dateNext,'init_time' )

!------------------------------------------------------------------------
! Establish if we are at the start/end of a month/year. newTime will also
! advance time to the start time.  Set "end flags", so that newTime sets
! "new flags".
!------------------------------------------------------------------------
  endMonth = .TRUE.
  endYear = .TRUE.
  endSec = .TRUE.

!------------------------------------------------------------------------
! 1st argument to newTime (timestep number) is zero, to indicate that
! none of the spin up code is to be activated.
!------------------------------------------------------------------------
  IF (TRIM(runtype)=='INITIAL') CALL newTime ( 0,endRun, time_BRAMS,frqhis_BRAMS)  !DSM - Dah problema no caso do history, portanto nao pode ser chamada na inicializacao

!------------------------------------------------------------------------
! Processing of spin up variables.
! Only required if there is any spin up.
!------------------------------------------------------------------------
  if ( nspin /= 0 ) then

!------------------------------------------------------------------------
!   Set up the "master" list of possible spin up variables.
!------------------------------------------------------------------------
    i = 0

    i=i+1; iposRouteStore=i; spinVarName(i)='routeStore'
    spinVarDesc(i)='routing channel store'

    i=i+1; iposSmcl=i; spinVarName(i)='smcl'
    spinVarDesc(i)='soil moisture content in layer (kg m-2)'

    i=i+1; iposTsoil=i; spinVarName(i)='t_soil'
    spinVarDesc(i)='subsurface (soil) temperature in layer (K)'

!------------------------------------------------------------------------
!   Check for repeated names in list.
!------------------------------------------------------------------------
    IF ( repeatVal(varNameTmp(1:nvarFound)) ) THEN
      WRITE(*,*)'ERROR: init_time: repeated name in list of spin-up variables.'
      STOP
    ENDIF

!------------------------------------------------------------------------
!   Identify which variables have been listed.
!------------------------------------------------------------------------
     CALL varList( nvarFound,.FALSE.,varNameTmp,spinVarName  &
                  ,errFound,spinVar  &
                  ,varConstIn=spinTolTmp,varConst=spinTol  &
                  ,varLog1In=spinTolPcTmp,varLog1=spinTolPercent )
     IF ( errFound ) THEN
       WRITE(*,*)'ERROR: init_time: error found in varList.'
       STOP
     ENDIF

!-------------------------------------------------------------------------------
!   Check that the chosen spin-up variables are available for this model setup.
!-------------------------------------------------------------------------------
    DO ivar=1,nspinVarMax
      IF ( spinVar(ivar) ) THEN
        IF ( ivar==iposSmcl .OR. ivar==iposTsoil ) THEN
!       Not useful if only doing routing.
          IF ( routeOnly ) THEN
            WRITE(*,*) 'init_time: ',TRIM( spinVarName(ivar) )
            WRITE(*,*)'Not useful as routeOnly=TRUE. Switching off.'
            spinVar(ivar) = .FALSE.
          ENDIF
        ELSEIF ( ivar == iposRouteStore ) THEN
!       Not useful if not doing routing.
          IF ( .NOT. route ) THEN
            WRITE(*,*) 'init_time: ',TRIM( spinVarName(ivar) )
            WRITE(*,*)'Not useful as route=FALSE. Switching off.'
            spinVar(ivar) = .FALSE.
          ENDIF
        ELSE
          WRITE(*,*)'ERROR: init_time: no code for name=',TRIM(spinVarName(ivar))
        ENDIF
      ENDIF  !  spinVar

!    If using a % tolerance, ensure it is >=0.
     IF ( spinVar(ivar) .AND. spinTolPercent(ivar) ) spinTol(ivar)=ABS( spinTol(ivar) )

    ENDDO  !  ivar

!-------------------------------------------------------------------------------
!   Count the number of spin up variables.
!-------------------------------------------------------------------------------
    nspinVar = COUNT( spinVar(:) )

    IF ( nspinVar == 0 ) THEN
      WRITE(*,*)'ERROR: init_time: nspinVar=0.'
      WRITE(*,*)'No variables have been selected for use in asessing spin up.'
      STOP
    ENDIF

  endif  !  nspin

!------------------------------------------------------------------------
! Some output to screen.
!------------------------------------------------------------------------
  IF ( echo ) THEN

    IF ( nspin == 0 ) THEN
      WRITE(*,*) 'There is NO spin-up period.'
    ELSE
      time_hms = s_to_chhmmss( timeRun(1) )
      WRITE(*,"(a,i8,a,tr1,a)") 'Spin-up period starts ',dateSpin(1),' time=',time_hms
      time_hms = s_to_chhmmss( timeRun(2) )
      WRITE(*,"(a,i8,a,tr1,a)") 'Spin-up period ends   ',dateSpin(2),' time=',time_hms
      WRITE(*,*) 'Spin up period repeated up to ',ABS(nspin),' times.'
      IF ( spinFail ) THEN
        WRITE(*,*)'If model is not then spun-up, run will stop'
      ELSE
        WRITE(*,*)'If model is not then spun-up, run will CONTINUE'
      ENDIF

      WRITE(*,*)'Variables and tolerances used in spin up:'
      DO ivar=1,nspinVarMax
        IF ( spinVar(ivar) ) THEN
          IF ( spinTolPercent(ivar) ) THEN
            WRITE(*,*) TRIM(spinVarDesc(ivar)),' tolerance=',spinTol(ivar),' percent'
          ELSE
            WRITE(*,*) TRIM(spinVarDesc(ivar)),' tolerance=',spinTol(ivar)
            IF ( spinTol(ivar) < 0.0 ) WRITE(*,*)'WARNING: spinTol<0 is never &
                   &satisfied. The full number of spin up cycles will be used.'
          ENDIF
        ENDIF
      ENDDO

    ENDIF  !  nspin

    time_hms = s_to_chhmmss( timeRun(1) )
    WRITE(*,"(a,i8,a,tr1,a)") 'Main run (after spin up) starts at ',dateMainRun(1),' time=',time_hms
    time_hms = s_to_chhmmss( timeRun(2) )
    WRITE(*,"(a,i8,a,tr1,a)") 'and ends at                        ',dateMainRun(2),' time=',time_hms

  ENDIF  !  echo

  END SUBROUTINE init_time

