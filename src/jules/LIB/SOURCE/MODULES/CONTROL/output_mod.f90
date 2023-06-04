! module output_mod.
! Module containing procedures used for output.
!
!-------------------------------------------------------------------------------
! Some notes, that had to go somewhere, on netCDF output.
! Provision for netCDF output has been bolted on to pre-existing code, which
! was sometimes quite painful. In particular the existing code was set up for
! a "GrADS" model, processing fields "level by level" and with provision for
! at most 1 "vertical" level (which made snow layer variables difficult).
! The i/o code eas in need of an overhaul even before the netCDF code came
! along!
! More:
! Snow layer variables might well be better treated as ntile separate
!   variables. This is done to some extent for GrADS output, but mainly in
!   how it is represented in a ctl file. The actual prcessing considers the
!   variable to have ntiles*nsmax levels. But it would be quite likely be better
!   to represent as ntiles separate variables from init_out onwards.
! Snow layer variables are represented in the code as A(land_pts,ntiles,nsmax),
!   but in netCDF output as A(land_pts,nsmax,ntiles), i.e. the order of the
!   dimensions has been changed. Again this was to fit in with existing code,
!   particularly for GrADS which represents as ntiles variables. But it makes
!   life more awkward than it would have been otherwise!
!-------------------------------------------------------------------------------

!###############################################################################

  MODULE output_mod

  INTEGER, PARAMETER :: nmax = 750  !  a size for various arrays

  CONTAINS

!###############################################################################
!###############################################################################
! subroutine output
! Driver for output routines.

  SUBROUTINE output( a_step,endCall )

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile

  USE inout, ONLY : formatNc,nlevMax,nout,ntCtlNeed,nxyMax,outActivePrev  &
                   ,outDate,outDateFlag,outFormat,outNpWrite,outWarnEarly,outFilePer  &
                   ,periodAnn,periodMon,periodOneFile,outFileStep  &
                   ,outFirstActive,outFirstSection,outFirstWrite,outLenWrite  &
                   ,outPer,outSamPer,outStep,outStepSamPer,outTime  &
                   ,outUnit,outWriteCount,snapProfile,tmeanProfile
  USE spin_mod, ONLY : ispin,spinUp
  USE time_loc, ONLY : date,dateNext,endMonth,endSec,endYear  &
                      ,stepFlag,time,timeNext
  USE time_mod, ONLY : timeDate_cmp
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    a_step   !  timestep number

  INTEGER ::  &!  local SCALARS
    iout    &!  loop counter
   ,outputDate,outputTime  !  date and time when this procedure was called (i.e. the time to
!   which the model has been integrated when this routine was called).

  REAL ::  &!  local ARRAYS
    outvalWrite(outLenWrite)  !  data as written to file.

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    endCall    !  T means this call is at the end of a model timestep
!                 F means this call is at the start of a model timestep - this is only
!                   needed when output is required at the start of the chosen time interval

  LOGICAL ::   &!  local SCALARS
    doAccum    &!  T if values are to be added to accumulation on this timestep
   ,doTmean    &!  T if time averages/accumulations are to be finished and output on this timestep
!                    Most times, doTmean occurs with doAccum (i.e. a final time is added
!                    to the accumulation, and then the average/accumulation is output).
!                    On occasions doTmean is T while doAccum=F - this is likely to be
!                    caused by the "early" end of a section of the run, at a timestep when
!                    nothing was expected to be added to the accumulation. In these cases,
!                    doTmean indicates that output is to be written, but doAccum=F means
!                    there is no further increment to accumulation.
   ,doSnap     &!  T if instantaneous values are to be output on this timestep
   ,doWrite    &!  T if values are to be written to file on this timestep
   ,lateFile   &!  T if a file is opened at a time later than "expected"
!                    (e.g. a daily file opened some time after start of day)
   ,newFile    &!  T is a new output file is to be opened on this timestep
   ,openedFile &!  T means newOutFile opened a new data file
   ,outActive   !  T when a profile is "active" at this time

  CHARACTER(len=7) ::  &!  local scalars
    periodStr   !  work

!-------------------------------------------------------------------------------

! Loop over output profiles.
  DO iout=1,nout

!   Initialise
    outActive = .FALSE.

!   Get time and date of output data - i.e. the time to which the model has been integrated
!   when this procedure was called.
    outputTime = timeNext
    outputDate = dateNext
    IF ( .NOT. endCall ) THEN
      outputTime = time
      outputDate = date
    ENDIF

!   Decide if output profile is currently active.
!   At endCall=T (meaning this call is at the end of a timestep), "active" means that the time is within
!   the range selected for output. The time used is the time at the end of the timestep (i.e.the time to
!   which the model has been integrated when this routine was called).
!   At endCall=F, there is a further condition that the profile must require output at the start of this timestep.
    SELECT CASE ( outDateFlag(iout) )
      CASE ( 1: )
!       The given dates are to be used, only after spinUp.
!       Generally only act at end of timestep (endCall).
        IF ( .NOT. spinUp ) THEN
          IF ( endCall ) THEN
!           Check if time (at the end of the timestep) is within range.
            IF ( timeDate_cmp( outputTime,outputDate,'>=',outTime(iout,1),outDate(iout,1),'output' ) .AND.  &
                 timeDate_cmp( outputTime,outputDate,'<=',outTime(iout,2),outDate(iout,2),'output' )  &
                ) outActive = .TRUE.
          ELSE
!           NOT endCall - this is needed to get output at the start of a time interval (eg initial state).
!           Activate if time (at start of current timestep) is the start date for output.
!           Also use outFirstActive to ensure that this is only used to generate the first output (so as to avoid
!           the case e.g. run starts 00H, output starts at end of first timestep - output is triggered at
!           end of first timestep (correct) and also at start of 2nd timestep (incorrect)).
            IF ( outFirstActive(iout) .AND.  &
                 timeDate_cmp( time,date,'=',outTime(iout,1),outDate(iout,1),'output' ) ) outActive = .TRUE.

          ENDIF
        ENDIF
      CASE ( 0 )
!       Output is always active, but only do anything at the end of timestep, or start of first timestep.
        IF ( endCall .OR. a_step==0 ) outActive = .TRUE.
      CASE ( -1 )
!       Output is active at all times after spin up (and there is spin up). Do something at the end of
!       timestep, or start of first timestep after spin up.
!       Note that stepFlag=3 will only occur with NOT spinUp, but test anyway (in case this changes!).
        IF ( .NOT.spinUp ) THEN
          IF ( endCall .OR. (.NOT.endcall .AND. stepFlag==3) ) outActive = .TRUE.
        ENDIF
      CASE ( -2 )
!       Output the initial state (at start of first timestep). Note that a_step=0 only
!       occurs with .NOT.endCall.
        IF ( a_step == 0 ) outActive = .TRUE.
    END SELECT

!XX For now....generally don't allow a call at start of timestep -  outStep, outStepFile can't always cope.
!XX Only allow a start of timestep call if this is the call at a_step=0 (inititalisation) and
!XX either this is the only output time (outDateFlag=-2),
!XX or if all timesteps are to be output (outDateFlag=0,outper=1) to a single file.
!XX REMOVE THIS RESTRICTION ASAP PLEASE!
    IF ( .NOT. endCall ) THEN
      outActive = .FALSE.
      IF ( a_step == 0 ) THEN
        IF ( outDateFlag(iout)==-2 .OR. ( outDateFlag(iout)==0  .AND.  &
                 outFilePer(iout)==periodOneFile .AND. outPer(iout)==1 ) ) outActive = .TRUE.
      ENDIF
    ENDIF

!    If a profile was active on the previous timestep, but it now inactive, check GrADS ctl file.
!     if ( .NOT.outActive .AND. outActivePrev ) call rewrite_ctl( checkTemplate,iout )

!   If a profile was active on the previous timestep but it now inactive, close
!   output files. This is particularly useful for netCDF files which otherwise
!   remain unfinished until end of run (and for other formats releases the units)
!   Only test this if endCall, as outActive generally FALSE at start of timestep -
!   not sure how robust this endCall bit is!
    IF ( endCall .AND. .NOT.outActive .AND. outActivePrev(iout) ) THEN
       CALL closeFile( outUnit(iout),outFormat )
!      Set outActivePrev to FALSE so the next timestep knows that output
!      on this profile was not active on this timestep.
       outActivePrev(iout) = .FALSE.
    ENDIF

!   If profile is not active at this time, nothing more to do.
    IF ( .NOT. outActive ) CYCLE

!-------------------------------------------------------------------------------
!   This profile is active.

!   Initialise
    doAccum = .FALSE.
    doTmean = .FALSE.
    doSnap = .FALSE.
    doWrite = .FALSE.

!   If this is the first time that this profile has been active in this "section" of the run,
!   reset counters.
    IF ( outFirstSection(iout) ) THEN
      CALL init_out_count( iout,outputTime,outputDate )
!     Reset switch.
      outFirstSection(iout) = .FALSE.
    ENDIF

!   Increment the counters of accumulated times in current output interval.
    outStep(iout) = outStep(iout) + 1
    outFileStep(iout) = outFileStep(iout) + 1

!   Test if it is time to add to accumulations.
!   Don't add to an accumulation if time is exactly at start of output.
!   Accumulation and time-average are defined as for tstart < t <= t_end - WHY??
!   Note that outStep may be <0 at start of output, hence use of mod.
!xx Does this cope with outAlways, now known as outDateFlag=0 and nspin>0?
    IF ( tmeanProfile(iout) .AND. MOD( outStep(iout),outSamPer(iout) )==0 ) doAccum = .TRUE.

!   Test if it is time to calculate and write output.
!   Deal with endSec separately - so we can catch early endMonth etc - done below (I think).
    IF ( ( outPer(iout)>0 .AND. MOD(outStep(iout),outPer(iout))==0 ) .OR.  &
         ( outPer(iout)==periodMon .AND. endMonth ) .OR.  &
         ( outPer(iout)==periodAnn .AND. endYear ) ) THEN
      IF ( tmeanProfile(iout) ) doTmean = .TRUE.
      IF ( snapProfile(iout) ) doSnap = .TRUE.
      doWrite=.TRUE.
!      if ( tstart ) then
!       If time equals the start time for this profile, we can output snapshot values,
!       and time averages will be undefined (because counter=0).
!          doSnap = .TRUE.
!          doWrite = .TRUE.
!          doTmean = .TRUE.
!      endif
    ENDIF

!XX Note that we have missed some cases. e.g. monthly average output, run starts at 00H
!XX on 1st of month and the IC is to be output - endMonth is currently F.
!XX Maybe have to test if this is timestep 0, and endCall=F, and time=00H on !st and then set endMonth=T??

!----------------------------------------------------------------------
!   Deal with incomplete output at the end of a 'section' of the run,
!   eg monthly averages, if run ends mid-month.
!   A monthly or annual variable may be output at the end of the section, even if the full time
!   interval was not simulated, but output on other periods is not written. The monthly or annual average is
!   calculated as long as there are "sufficient" times in the average, since it seems a pity to
!   abandon e.g. the annual average just because it is only midday on 31st December....
!   However, shorter averages (eg 2 hourly) are abandoned and nothing written.
!   Note that the case where data were not available at the start (rather than end) of an averaging period is
!   considered in subroutine loadout, where a threshold for data numbers is used.

!XX Code should probably also deal with incomplete monthly/annual output if run starts at 00H
!XX on 1st of month and the IC is to be output. e.g. monthly ave will only have one value (00H 1st).
!XX This will not be caught by endSec, since it's not the end of a section.
!
!  In future we might want to ensure that time averages or accumulations are always correct, even if
!  if output starts or ends part way through average's sampling period (e.g. daily sampling, run starts midday),
!  or we could make the value undefined. Need to alter weights to give less weight to sampled value from
!  a shorter interval. Possibly: in init_out_count,set a flag if outStep/=0, or if monthly output is not
!  starting at start of month. Or have to reliably calculate how many times are expected in output
!  (maybe harder). May also be hard to extend to "intermittent" variables, such as snow albedo.

    IF ( endSec .AND. .NOT.doWrite ) THEN

      outWarnEarly = .TRUE.

      SELECT CASE ( outper(iout) )

        CASE ( periodAnn, periodMon )
!         Expect endYear/endMonth to be FALSE, since if TRUE, doWrite should be TRUE and we would not be here!
!         Output snapshot at current time. Only output time average if >=90% of the expected data exist.
          doWrite = .TRUE.
          periodStr = 'annual'
          IF ( outper(iout)==periodMon ) periodStr='monthly'
          IF ( snapProfile(iout) ) THEN
             doSnap = .TRUE.
             WRITE(*,"(a,tr1,a,tr1,a,i2,a,i8,i6)")  &
             'WARNING: outWarnEarly:',TRIM(periodStr),'snapshot output is early. iout='  &
                         ,iout,' date,time=',date,time
          ENDIF
          IF ( tmeanProfile(iout) ) THEN
!           Calculate the average, but don't reset the accumulation flag - so only add to accumulation if
!           this time fits into sampling period.
            doTmean = .TRUE.
            WRITE(*,"(a,tr1,a,tr1,a,i2,a,i8,i6)")  &
               'WARNING: outWarnEarly:'  &
               ,TRIM(periodStr),'average calculated before end of interval. iout=',iout  &
               ,' date,time=',date,time
            WRITE(*,*)'All available data (on sampling frequency) to this time are included.'
          ENDIF   !  tmeanProfile

        CASE default
!         Regular output.
          WRITE(*,*) 'WARNING: outWarnEarly: endSec and no output written because of time mismatch.'
          WRITE(*,*) 'Output profile #',iout,' date,time at end of timestep=',dateNext,timeNext
          WRITE(*,*) 'This probably means that the end of a "section" of the run (e.g. spin up)'
          WRITE(*,*) 'does not fall at a timestep when output would be generated - so there''s no output.'

      END SELECT   !  outPer

!     Also indicate if a new file is to be opened for monthly or annual output (at this stage, doWrite=T only for
!     monthly or annual output). Again, nothing done for other output periods (because no output done).
      SELECT CASE ( outFileper(iout) )
        CASE ( periodAnn, periodMon )
          IF ( doWrite ) outWriteCount(iout) = 0
      END SELECT

    ENDIF   !    endSec .AND. .NOT.doWrite

!------------------------------------------------------------------------------
!  If a time-accumulation is to be incremented, also increment the counter.
   IF ( doAccum ) outStepSamPer(iout) = outStepSamPer(iout) + 1
!-------------------------------------------------------------------------------
!   If a file is to be written to, decide if a new file is to be opened.

    IF ( doWrite ) THEN
      newFile = .FALSE.
      lateFile = .FALSE.

      IF (  &
!         Regular file period
           ( outFilePer(iout)>0 .AND. outFilePer(iout)==outFileStep(iout) ) &
!         Monthly files.
!         Open a new file at end Month/endYear (rather than newMonth/newYear), since end flag is set at
!         the timestep that ends at the end of the month - the time is then 00H on 1st of next month,
!         when data should be in next month's file (by GrADS's convention).
          .OR. ( outFilePer(iout)==periodMon .AND. endMonth )  &
!         Annual files.
          .OR. ( outFilePer(iout)==periodAnn .AND. endYear )   &
!         Separate files during spin-up.
          .OR. ( outFilePer(iout)==-7 .AND. ( stepFlag==1 .OR. stepFlag==3 ) )  &
          .OR. ( outFilePer(iout)==-8 .AND. stepFlag==1 .AND. ispin==1 )        &
!         All output to a single file.
          .OR. ( outFilePer(iout)==periodOneFile .AND. outFirstWrite(iout) )    &
         ) newFile = .TRUE.

!      A new file will not have been triggered if this is the first output in a run (or
!      section of a run) and the time is not that usually expected for a new file.
!      e.g. timestep output to daily files, run starting midday - can't wait until start of next day
       IF ( outWriteCount(iout)==0 .AND. .NOT.newFile ) THEN
         newFile = .TRUE.
         lateFile = .TRUE.   !  this file is being opened later than expected
!                               lateFile is used to avoid resetting outFileStep
       ENDIF


    ENDIF  !  doWrite
!-------------------------------------------------------------------------------
!    If it is time to write to file, and a new file is needed, open it.

     IF ( doWrite .AND. newFile ) THEN
!      5th argument (FALSE) is compressGridCall - always FALSE here.
       CALL newOutFile( iout,outputDate,outputTime,endCall,.FALSE.,openedFile )
       IF ( openedFile ) THEN
!        Reset counters.
         outWriteCount(iout) = 0
         IF ( .NOT. lateFile ) outFileStep(iout) = 0    !  only reset if file was opened at "expected" time
         outFirstWrite(iout) = .FALSE.    !  indicates that profile has been written to in this run
       ENDIF
     ENDIF

!-------------------------------------------------------------------------------
!    Call routine to load data into output variables.
!    No need to call unless there is either accumulation, calculation of
!    time-average, or snapshot to do.
!-------------------------------------------------------------------------------
     IF ( doAccum .OR. doTmean .OR. doSnap )  &
            CALL loadOut( doAccum,doTmean,iout,nlevMax(iout),nxyMax(iout)  &
                        ,outputDate,outputTime,outvalWrite )

!-------------------------------------------------------------------------------
!    If any data are to be output, call writing routine.
!-------------------------------------------------------------------------------
     IF ( doWrite ) THEN

!      Increment counter of number of times written to file,
       outWriteCount(iout) = outWriteCount(iout) + 1
       ntCtlNeed(iout) = ntCtlNeed(iout) + 1

       CALL writeOut( a_step,iout,outWriteCount(iout),outNpWrite(iout)  &
                     ,outNpWrite(iout)*nlevMax(iout)  &
                     ,outputDate,outputTime,outvalWrite )

!      Reset output variables and counter.
       CALL reset_outval( iout )

    ENDIF

!-------------------------------------------------------------------------------
!   Reset flag, now that profile has been written to.
!-------------------------------------------------------------------------------
    IF ( outFirstActive(iout) ) outFirstActive(iout) = .FALSE.

!   Set flag so that next timestep knows that output was active on this timestep.
    outActivePrev(iout) = .TRUE.


  ENDDO  !  iout (profiles)

  END SUBROUTINE output

!###############################################################################
!###############################################################################
!###############################################################################

  SUBROUTINE newOutFile( iout,outputDate,outputTime,endCall,compressGridCall &
                        ,openedFile,compUnit )

! Driver to get name of a new output file and open it.

  USE ancil_info, ONLY : dim_cs1,nsmax,nsoil=>sm_levels,ntiles
  USE file_utils, ONLY : fileUnit,openFile
  USE inout, ONLY : echo,formatAsc,formatBin,formatLen,formatNc  &
                   ,gradsNc,haveSCpool,havePFT,haveSnow,haveSoil,haveTile  &
                   ,haveType,irecPrevOut,mapOutCompress,ntCtl,ntCtlNeed  &
                   ,ntOutFilePer,nvarOut,openedFileName  &
                   ,outCompress,outDataFile,outDateFlag,outFilePer,outGrADS  &
                   ,compressGridFile,outGridDxy,outGridNxy  &
                   ,outGridXY,outNpWrite,periodAnn,periodMon,periodOneFile  &
                   ,pointsOut,outFirstWrite  &
                   ,outPer,outStatus,outFormat,outUnit,outVarID,redoTdef  &
                   ,outTemplate,taccumVar,tmeanVar,undefOut  &
                   ,varDesc,varName,varNum,varPos,varType,varUnitsList

  USE jules_netcdf, ONLY :  &
!  imported procedures
     createNcFile

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     npft,ntype

  USE spin_mod, ONLY : ispin,nspin,spinUp
  USE switches, ONLY : l_360
  USE time_loc, ONLY : date,dateMainRun,dateSpin,timeStep,time,timeRun
  USE time_mod, ONLY : timeDate,timeDate_diff,dateToBits,getHours
  USE timeConst, ONLY : iSecInHour

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: iout         !  number of output profile
  INTEGER, INTENT(in) :: outputDate   !  date (yyyymmdd) of output
  INTEGER, INTENT(in) :: outputTime   !  time of day (s) of output

  LOGICAL ,INTENT(in) :: endCall          !  T means this call is at the end of a timestep
  LOGICAL ,INTENT(in) :: compressGridCall !  T means this call is for the
!                    supplementary data file that describes the mapping used
!                    when compressing output (e.g. the supplementary file
!                    referred to in the GrADS pdef entry).
!                    Note this is FALSE when just writing a ctl file that
!                    contains a pdef entry.

!-------------------------------------------------------------------------------
! Scalar arguments with intent(out)
!-------------------------------------------------------------------------------
  LOGICAL, INTENT(out) :: openedFile  !  T means that this subroutine opened
!                                         a new data file - the usual case
!                      F means no file was opened - occasionally happens if the
!                        required file was opened earlier in timestep.

!-------------------------------------------------------------------------------
! Optional scalar arguments with intent(out)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(out), OPTIONAL :: compUnit  !  unit used to connect to file in
!    the case of compressGridCall=TRUE.

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: lenPer = 8  ! length of character variables used to
!                  to represent period. len=8 is sufficient for longest
!                  possible (i4.4'step')

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER :: i1,i2,i3,ivar,jvar,kvar,ndy,nhr,nmin,nmon,nsec,ntMain,ntOutTot  &
   ,ntSpin,nxy,nyr,nz,timeTmp  !  work
  INTEGER :: nvarCtl  !  number of variables, as represented in ctl file
  INTEGER :: useUnit  !  unit used to connect to file

  REAL :: hr    !  number of hours
  LOGICAL :: checkTemplate  !  T means check previous template ctl file has
!                                correct number of times
  LOGICAL :: writeCtl       !  work

  CHARACTER(len=8) :: callType         !  argument to create_annotation
  CHARACTER(len=lenper) :: cper      !  used to hold output period
  CHARACTER(len=1) :: cwork  !  work
  CHARACTER(len=LEN(outDataFile)) :: fileName   !  work
  CHARACTER(len=20) :: dateTimeString  !  date and time as yyyy-mm-dd hh:mm:ss
  CHARACTER(len=200) :: dsetLine       !  part of the line used for DSET in
!                                     GrADS ctl file (hopefully long enough!)
  CHARACTER(len=formatLen) :: fileFormat   !  work
  CHARACTER(len=100) :: optionsLine    !  GrADS options line
  CHARACTER(len=7) :: periodUnitsLong  !  units of dataperiod (e.g. seconds)
  CHARACTER(len=100) :: tdef           !  (most of) the line used in GrADS tdef
  CHARACTER(len=100) :: titleLine      !  GrADS title line
  CHARACTER(len=100) :: zdefLine       !  GrADS zdef line

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER :: dateTmp(2),day(2),month(2),year(2)   !  work
  CHARACTER(len=200) :: annot(nmax)  !  annotation (metadata)
  CHARACTER(len=200) :: tdefAnnot(nmax)!  annotation (metadata)
  CHARACTER(len=200) :: varCtlLine(nmax)   !  descriptions of variables for ctl file
  CHARACTER(len=100) :: varNcAtt(nmax,4)  !  attributes of variables for
!                     netCDf output.
!                     1=axes, 2=units, 3=name, 4=description

!-------------------------------------------------------------------------------
! Initially indicate that a file will be opened and a ctl file written for GrADS.
! Write a ctl file for all netCDF files, for now at least!
!-------------------------------------------------------------------------------
  openedFile = .TRUE.
  IF ( outGrADS .OR. outFormat==formatNc ) THEN
    writeCtl = .TRUE.
  ELSE
    writeCtl = .FALSE.
  ENDIF

!-------------------------------------------------------------------------------
! Get name for new file.
!-------------------------------------------------------------------------------
  CALL newOutFileName( cper,compressGridCall,outputDate,outputTime  &
                      ,dateTimeString,periodUnitsLong,dsetLine,iout,tdef )

!-------------------------------------------------------------------------------
! Decide whether there really is a need for a new file.
!-------------------------------------------------------------------------------
  IF ( .NOT. compressGridCall ) THEN
    IF ( .NOT. endCall ) THEN
!     A file will be opened.
!     Save the name of the file that is to be opened.
      openedFileName(iout) = outDataFile(iout)
    ELSE
!     Decide if a file with this name was opened at the start of this timestep.
      IF ( openedFileName(iout) == outDataFile(iout) ) THEN
        openedFile = .FALSE.  !  this call will not open a new file
        IF ( echo ) THEN
          WRITE(*,*)'newOutFile: not opening file as already open.'
          WRITE(*,*)'file: ',TRIM(outDataFile(iout))
         ENDIF
        RETURN   !  nothing more to do
      ENDIF
    ENDIF  !  endCall

!-------------------------------------------------------------------------------
!    Before a new data file is opened, check that GrADS tdef for the previous
!    ctl file was correct.
!-------------------------------------------------------------------------------
     checkTemplate = .FALSE.   !   don't call rewrite_ctl to check template tdef

!-------------------------------------------------------------------------------
!    Only check a template ctl once all data have been written to the associated
!    data files. At present a template ctl need only be checked at the end of
!    the run (done elsewhere) or at the start of a new "section" of the run -
!    test for the first time data are written in this section, then further make
!    sure that this is not the first section of the run (which won't have an
!    earlier ctl file to check). Note if no spin up, only one section, so
!    nothing to do.
!xx  if there was a call at endCall=F, the correct ctl was done then, and we
!    since have a new ctl.
!-------------------------------------------------------------------------------
     IF ( outTemplate(iout) .AND. outFirstWrite(iout) .AND. nspin/=0 .AND.  &
          outDateFlag(iout)==0 ) THEN
       IF ( .NOT. spinUp .OR. (spinUp .AND. ispin>1) ) THEN
         checkTemplate = .TRUE.
!-------------------------------------------------------------------------------
!        If the file names will carry on a sequence, there is no need to check a
!        template ctl until the end of the run - a single template ctl will
!        serve for the whole run. This is only the case if there is a single
!        cycle of spin up, and the main run start date equals the end date of
!         spin up.
!-------------------------------------------------------------------------------
         IF ( nspin==1 .AND. dateMainRun(1)==dateSpin(2) ) checkTemplate=.FALSE.
       ENDIF
     ENDIF
     IF ( redoTdef .AND. .NOT.compressGridCall )   &
                 CALL rewrite_ctl( checkTemplate,iout )

  ENDIF   !  compressGridCall

!--------------------------------------------------------------------------
! Work out how many times are expected to be written to this file.
! There are umpteen combinations.....
! Some of these are estimated without worrying too much about trying to get
! it right, since if redoTdef=T is used (as it always is now), this number
! is later corrected once we know how many times were actually written.
! For example, interuptions by spin up may not be correctly counted. Also,
! no account is taken of the fact that output may be requested for only
! part of the time - again, this is corrected later.
! It could be argued that if we're not going to make the effort to get
! it exactly right, we'd be as well just using some constant (e.g. 1000)
! and correcting later!
! A template ctl file needs to know the total number of times across all files.
!
! Also set the grid size (to set length of record).
!-------------------------------------------------------------------------------

! Initialise ntSpin.
  ntSpin = 0

  IF ( compressGridCall ) THEN

!   Simple case if this call is to write compression mapping (e.g. GrADS pdef file).
    nxy = outGridNxy(iout,1) * outGridNxy(iout,2)
    ntOutFilePer(iout) = 1
    ntCtl(iout) = ntOutFilePer(iout)
    fileName = compressGridFile(iout)
!   Note that even when netCDF output is requested, the supplementary (pdef)
!   data file is a flat binary file, so that GrADS can read the netCDF files.
!   If you want the supplementary file to also be netCDF, you'll have to
!   write more code here!
    IF ( outGrads .OR. outFormat==formatNc ) THEN
       fileFormat = formatBin
    ELSE
!     This file will be of same type as other requested output (e.g. ASCII).
      fileFormat = outFormat
    ENDIF

  ELSE

!-------------------------------------------------------------------------------
!   Work out if a new ctl file is to be written. May not be needed if template
!   ctl is used. A template ctl is written at first call, but the start of a new
!   section of the run also forces a new template ctl file (if names change - if
!   names don't change newSec=F anyway).
!-------------------------------------------------------------------------------
    IF ( ( outGrADS .OR. outFormat==formatNc ) .AND. outTemplate(iout) ) THEN
      writeCtl = .FALSE.
      IF ( outFirstWrite(iout) ) writeCtl = .TRUE.
    ENDIF

    nxy = outNpWrite(iout)
    filename = outDataFile(iout)
    fileFormat = outFormat

    IF ( .NOT. outTemplate(iout) ) THEN

      SELECT CASE ( outFilePer(iout) )

        CASE ( periodAnn )
!         Assume the full year will be used.
          IF ( outPer(iout) == periodAnn ) THEN
            ntoutFilePer(iout) = 1     !    annual file, annual output
          ELSEIF ( outPer(iout) == periodMon ) THEN
            ntoutFilePer(iout) = 12    !    annual file, monthly output
          ELSEIF ( outPer(iout) > 0 ) THEN
!           Get time from now until just before start of next year (start of
!           next year is in next file).
            CALL dateToBits( date,day(1),month(1),year(1),l_360,'newoutfile' )
            dateTmp(1) = (year(1)+1)*10000 + 101   !  1st Jan next year
            CALL timeDate(0,dateTmp(1),-10,'sec',l_360,timeTmp,dateTmp(2),'newoutfile')
            hr = getHours( time,date,timeTmp,dateTmp(2),l_360,'newoutfile')
            ntOutFilePer(iout) = INT( hr*REAL(iSecInHour)  &
                                     / REAL(outPer(iout)*NINT(timeStep)) )
          ELSE
            WRITE(*,*)'ERROR: newoutfile: 1: no code for outPer=',outPer(iout)
            STOP
          ENDIF

        CASE ( periodMon )
!         Assume the full month will be used.
          IF ( outPer(iout) == periodMon ) THEN
            ntoutFilePer(iout) = 1    !    monthly file, monthly output
          ELSEIF ( outPer(iout) > 0 ) THEN
!           Get time from now until just before start of next month (start of
!           next month is in next file).
            CALL timeDate( time,date,1,'mon',l_360,timeTmp,dateTmp(1),'newoutfile' )
            CALL dateToBits( dateTmp(1),day(1),month(1),year(1),l_360,'newoutfile' )
            dateTmp(1) = year(1)*10000 + month(1)*100 + 1   !  1st of next month
            CALL timeDate( 0,dateTmp(1),-10,'sec',l_360,timeTmp,dateTmp(2)  &
                ,'newoutfile')   !  this gives 10s before next month
            hr = getHours( time,date,timeTmp,dateTmp(2),l_360,'newoutfile' )
            ntOutFilePer(iout) = INT( hr*REAL(iSecInHour)  &
                 / REAL(outPer(iout)*NINT(timeStep)) )
          ELSE
            WRITE(*,*)'ERROR: newoutfile: 1: no code for outPer=',outPer(iout)
            STOP
          ENDIF

        CASE ( 1: )
!         Output with constant period.
          ntOutFilePer(iout) = outFilePer(iout)/outPer(iout)

        CASE ( periodOneFile,-8,-7 )
!         -7: Each spin up cycle to a separate file.
!         -8: All spin up cycles to one file.
!         periodOneFile: All times to one file.
!         Each can have various output periods.
!         First, get number of timesteps in each section.

!         Cases that involve a spin up cycle.
          IF ( spinUp ) THEN
!           We may be here on occasions when we don't need to be - but never mind.
            IF ( outPer(iout) == periodAnn ) THEN
!             Get number of years in a spin up cycle.
              ntSpin = MAX( (dateSpin(2)-dateSpin(1))/10000, 1 )
            ELSEIF ( outPer(iout) == periodMon ) THEN
!             Get number of months in a spin up cycle.
              CALL timeDate_diff( 0,dateSpin(1),0,dateSpin(2),l_360,'newoutfile' &
                           ,nsec,nmin,nhr,ndy,nmon,nyr )
              ntSpin = nyr*12 + nmon
!             Unless interval is complete months, starting 00H on 1st, we need a file
!             for the last month (eg Jan to Jan, often need 13 files).
              CALL dateToBits( dateSpin(1),day(1),month(1),year(1),l_360,'newoutfile' )
              IF ( timeRun(1)>0 .OR. day(1)>1 ) ntSpin=ntSpin+1
            ELSEIF ( outPer(iout) > 0 ) THEN
!             Get number of hours from now until end of spin up cycle.
              hr = getHours( time,date,timeRun(1),dateSpin(2),l_360,'newoutfile' )
              ntSpin = INT( (hr*REAL(iSecInHour)) / REAL(outPer(iout)*NINT(timeStep)) )
            ELSE
              WRITE(*,*)'ERROR: newoutfile: 2: no code for outPer=',outPer(iout)
              STOP
            ENDIF
          ENDIF  !  spinUp

!         Cases that involve the main run.
          IF ( ( (outFilePer(iout)==-7 .OR. outFilePer(iout)==-8) .AND. .NOT.spinUp ) .OR.  &
                outFilePer(iout)==periodOneFile ) THEN
!           We may be here on occasions when we don't need to be - but never mind.
            IF ( outPer(iout) == periodAnn ) THEN
!             Get number of years in main section.
              ntMain = MAX( (dateMainRun(2)-dateMainRun(1))/10000, 1 )
            ELSEIF ( outPer(iout) == periodMon ) THEN
!             Get number of months in main section.
              CALL timeDate_diff( timeRun(1),dateMainRun(1),timeRun(2),dateMainRun(2)  &
                                 ,l_360,'newoutfile',nsec,nmin,nhr,ndy,nmon,nyr )
              ntMain = nyr*12 + nmon
!             Unless interval is complete months, starting 00H on 1st, we need a file
!             for the last month (eg Jan to Jan, often need 13 files).
              CALL dateToBits( dateMainRun(1),day(1),month(1),year(1),l_360,'newoutfile' )
              IF ( timeRun(1)>0 .OR. day(1)>1 ) ntMain=ntMain+1
            ELSEIF ( outPer(iout) > 0 ) THEN
!             Get number of hours in main section, then convert to # of output period.
              hr = getHours( timeRun(1),dateMainRun(1),timeRun(2),dateMainRun(2),l_360,'newoutfile' )
!             Add 1 below in case start of run is also an output time - not testing.
              ntMain = INT( (hr*REAL(iSecInHour)) / REAL(outPer(iout)*NINT(timeStep)) ) + 1
            ELSE
              WRITE(*,*)'ERROR: newoutfile: 3: no code for outPer=',outPer(iout)
              STOP
            ENDIF
          ENDIF

!         Now use ntSpin and ntMain to get final number of times.
          ntOutFilePer(iout) = 0
          IF ( outFilePer(iout) == -7 ) THEN
            IF ( spinUp ) THEN
              ntOutFilePer(iout) = ntSpin
            ELSE
              ntOutFilePer(iout) = ntMain
            ENDIF
          ELSEIF ( outFilePer(iout) == -8 ) THEN
            IF ( spinUp ) THEN
              ntOutFilePer(iout) = ntSpin * nspin
            ELSE
              ntOutFilePer(iout) = ntMain
            ENDIF
          ELSEIF ( outFilePer(iout) == periodOneFile ) THEN
              ntOutFilePer(iout) = ntSpin*nspin + ntMain
          ENDIF

        CASE default
          WRITE(*,*)'ERROR: newoutfile: no code for outFilePer=',outFilePer(iout)
          STOP

      END SELECT  !  outFilePer


!     Save this value.
      ntCtl(iout) = ntOutFilePer(iout)

!  ----------------------------------------------------------------------
    ELSEIF ( outTemplate(iout) .AND. writeCtl ) THEN

!     Template ctl file.
!     Get number of output times over whole run (since only allowing one template ctl per run).
!     Template is not used for the outFilePer=-7,-8,periodOneFile.
!     This code is very similar to that above for cases -7:periodOneFile - could probably combine.

      ntOutTot = 0

!     Cases that involve a spin up cycle - in most cases a separate ctl is written for each spin-up cycle,
!     so just assume that is the case.
      IF ( spinUp ) THEN
        IF ( outPer(iout) == periodAnn ) THEN
!         Get number of years in a spin up cycle.
          ntSpin = MAX( (dateSpin(2)-dateSpin(1))/10000, 1 )
        ELSEIF ( outPer(iout) == periodMon ) THEN
!         Get number of months in a spin up cycle.
          CALL timeDate_diff( 0,dateSpin(1),0,dateSpin(2),l_360,'newoutfile' &
                             ,nsec,nmin,nhr,ndy,nmon,nyr )
          ntSpin = nyr*12 + nmon
!         Unless interval is complete months, starting 00H on 1st, we need a file
!         for the last month (eg Jan to Jan, often need 13 files).
          CALL dateToBits( dateSpin(1),day(1),month(1),year(1),l_360,'newoutfile' )
          IF ( timeRun(1)>0 .OR. day(1)>1 ) ntSpin=ntSpin+1
        ELSEIF ( outPer(iout) > 0 ) THEN
!         Get number of hours from now until end of spin up cycle, then convert to # of output periods.
          hr = getHours( time,date,timeRun(1),dateSpin(2),l_360,'newoutfile' )
          ntSpin = INT( (hr*REAL(iSecInHour)) / REAL(outPer(iout)*NINT(timeStep)) )
        ELSE
          WRITE(*,*)'ERROR: newoutfile: template 1: no code for outPer=',outPer(iout)
          STOP
        ENDIF
        ntOutTot = ntSpin

      ELSE
!       NOT spinUp
!       Main run.
        IF ( outPer(iout) == periodAnn ) THEN
!         Get number of years in main section.
          ntMain = MAX( (dateMainRun(2)-dateMainRun(1))/10000, 1 )
        ELSEIF ( outPer(iout) == periodMon ) THEN
!         Get number of months in main section.
          CALL timeDate_diff( timeRun(1),dateMainRun(1),timeRun(2),dateMainRun(2)  &
                             ,l_360,'newoutfile',nsec,nmin,nhr,ndy,nmon,nyr )
          ntMain = nyr*12 + nmon
!         Unless interval is complete months, starting 00H on 1st, we need a file
!         for the last month (eg Jan to Jan, often need 13 files).
          CALL dateToBits( dateMainRun(1),day(1),month(1),year(1),l_360,'newoutfile' )
          IF ( timeRun(1)>0 .OR. day(1)>1 ) ntMain=ntMain+1
        ELSEIF ( outPer(iout) > 0 ) THEN
!         Get number of hours in main section.
          hr = getHours( timeRun(1),dateMainRun(1),timeRun(2),dateMainRun(2),l_360,'newoutfile' )
          ntMain = INT( (hr*REAL(iSecInHour)) / REAL(outPer(iout)*NINT(timeStep)) )
        ELSE
          WRITE(*,*)'ERROR: template 2: no code for outPer=',outPer(iout)
          STOP
        ENDIF
        ntOutTot = ntMain
      ENDIF !   spinUp

!     Save this value.
      ntCtl(iout) = ntOutTot - 1
!----------------------------------------------------------------------

    ENDIF  !  template
  ENDIF    !  compressGridCall

! Make sure we have at least one time.
  ntCtl(iout) = MAX( ntCtl(iout), 1 )

  IF ( outDateFlag(iout) == -2 ) ntCtl(iout)=1

!--------------------------------------------------------------------------
! Create metadata and annotation that is useful for a GrADS ctl file,
! netCDF attributes, etc. Much of this might not be needed (e.g. if a ctl
! file is not written for this data file, because template ctl used), but
! is always done as it generates some useful stuff.
!--------------------------------------------------------------------------
  IF ( compressGridCall ) THEN
    callType = 'compress'
  ELSE
    callType = 'diag'
  ENDIF
  CALL create_annotation( iout,callType,nvarCtl  &
                  ,optionsLine,titleLine,zdefLine,annot,tdefAnnot  &
                  ,varCtlLine,varNcAtt )

!--------------------------------------------------------------------------
! Get unit number for file.
!--------------------------------------------------------------------------
  IF ( compressGridCall ) THEN
!   Get a new unit to use for this pdef data file. We need to avoid using
!   outUnit in this case because netCDF output does not have a useful unit
!   until the dataset is created.
    useUnit = fileUnit( fileFormat )
  ELSE
!   Use the preset unit for this output profile.
    useUnit = outUnit(iout)
  ENDIF

!--------------------------------------------------------------------------
! Open the output data file.
!--------------------------------------------------------------------------
! 2nd argument=TRUE will close any open netCDF file (needed for netCDF).
  CALL openFile( nxy,.TRUE.,useUnit,'write',fileFormat,fileName,outStatus )
  irecPrevOut(iout) = 0

!--------------------------------------------------------------------------
! Save the unit number.
! This is only required for netCDF files (otherwise has intent(in)), but
! do anyway. For compressGridCall, save to optional argument (which is
! assumed to have been provided!).
!--------------------------------------------------------------------------
  IF ( compressGridCall ) THEN
    compUnit = useUnit
  ELSE
    outUnit(iout) = useUnit
  ENDIF

!--------------------------------------------------------------------------
! Create dimensions etc for a netCDF file.
! This is kept separate from openFile to avoid complicating that interface!
!--------------------------------------------------------------------------
  IF ( fileFormat == formatNc ) THEN
!   Get indices needed so we just pass the bits of arrays that refer to
!   this output profile.
    i1 = varPos(iout,1)
    i2 = i1 + nvarOut(iout) - 1
!   Use nvarCtl to cope with nasty snow layer cases. Note that this requires
!   that the use/calculation of nvarCtl does not change! It's a hack!
    i3 = i1 + nvarCtl - 1
    CALL createNcFile( iout,outUnit(iout),outGridNxy(iout,1),outGridNxy(iout,2)  &
                           ,pointsOut(iout),nvarOut(iout)  &
                           ,outGridDxy(iout,1),outGridDxy(iout,2)  &
                           ,outGridXY(iout,1),outGridXY(iout,2),gradsNc    &
                           ,haveSCpool(iout),havePFT(iout),haveSnow(iout)  &
                           ,haveSoil(iout),haveTile(iout),haveType(iout)   &
! The IMOGEN prognostics are not available for output yet
                           ,.FALSE.,.FALSE.,.FALSE.,.FALSE.  &
                           ,outCompress(iout),callType  &
                           ,compressGridFile(iout)  &
                           ,dateTimeString  &
                           ,fileName,'Diagnostic output'  &
                           ,mapOutCompress(iout,1:pointsOut(iout))  &
                           ,varName(i1:i2),varNcAtt(i1:i3,:)  &
                           ,varType(i1:i2),outVarID(i1:i2) )
  ENDIF

!--------------------------------------------------------------------------
! Reset counter of times needed for GrADS tdef.
!--------------------------------------------------------------------------
  IF ( compressGridCall .OR. .NOT.outTemplate(iout) .OR.  &
      ( outTemplate(iout) .AND. outFirstWrite(iout) ) ) ntCtlNeed(iout) = 0

!--------------------------------------------------------------------------
! Data in compression mapping file are not written via usual code, so
! ntCtlNeed needs to be set here.
! Assuming there is only 1 time of data in mapping file!
!--------------------------------------------------------------------------
  IF ( compressGridCall ) ntCtlNeed(iout) = 1

!--------------------------------------------------------------------------------
! Write a GrADS ctl file, if appropriate.
! This must be done after any new data file is opened, so that redoTdef gets
! correct name of ctl file.
!--------------------------------------------------------------------------------
!  if ( outGrADS .AND.  &
!       ( ( .NOT.outTemplate(iout) ) .OR. ( outWriteCount(iout)==1 .AND. outTemplate(iout) ) .OR.    &
!         ( newSec .AND. outTemplate(iout) ) ) ) call writeGradsCtl (  &
!           iout,outDataFile(iout),dsetLine,cper,tdef,compressGridCall)

  IF ( writeCtl ) CALL writeGradsCtl ( iout,nvarCtl,dsetLine,optionsLine  &
                           ,titleLine,zdefLine  &
                           ,cper,tdef,compressGridCall  &
                           ,annot,tdefAnnot,varCtlLine )

!----------------------------------------------------------------------
! If writing ASCII output, write headers and list of variable names.
! Always write the same number of headers - easier for automatic
! processing of output. Annotation/headers will be particularly useful
! for single point output. For grids would probably also like to know
! grid size, which might vary between variables.
! Also should indicate if zrevOutSoil or zrevOutSnow are set.
!--------------------------------------------------------------------------------

  IF ( outFormat == formatAsc ) THEN
    IF ( .NOT. compressGridCall ) THEN
      WRITE(outUnit(iout),"(i3,a)") nvarOut(iout),' variables in this file.'
      WRITE(outUnit(iout),"(a)") '# timeFlag, number of levels, name, description'
!     NOTE: e.g. monthly output, run ends midmonth, snapshot will be early, so this message is wrong...
      WRITE(outUnit(iout),"(a)") '# S denotes snapshot value at given time.'
      IF ( outPer(iout)>0 ) THEN
        WRITE(outUnit(iout),"(a,i6,a)") '# M denotes backward time average over '  &
           ,outPer(iout)*NINT(timeStep),'s, ending at given time.'
        WRITE(outUnit(iout),"(a,i6,a)") '# A denotes time-accumulation*timestep (s), over '  &
           ,outPer(iout)*NINT(timeStep),'s, ending at given time.'
      ELSE
        WRITE(outUnit(iout),"(a)") '# M denotes backward time average over month/year, ending at given time.'
        WRITE(outUnit(iout),"(a)") '# A denotes time-accumulation*timestep (s), ending at given time.'
      ENDIF
      DO ivar=1,nvarOut(iout)
        jvar = varPos(iout,ivar)
        kvar = varNum(jvar)       !  number in list of available vars.
        cwork = 'S'  ! snapshot variable
        IF ( taccumVar(jvar) ) cwork = 'A'
        IF ( tmeanVar(jvar) ) cwork = 'M'
        SELECT CASE ( varType(jvar) )
          CASE ( 'LA' )
            nz = 1
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a)")  &
                         cwork,nz  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')','defined on land points'
          CASE ( 'PF' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a)") &
                          cwork,npft  &
                         ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                         ,'(',TRIM(varUnitsList(kvar)),')',', for ',npft,' PFTs'
          CASE ( 'RG' )
            nz = 1
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a)")  &
                         cwork,nz  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')',' on river routing grid'
          CASE ( 'RP' )
            nz = 1
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a)")  &
                         cwork,nz  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')'  &
                        ,', for selected point on river routing grid'
          CASE ( 'SC' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a)")   &
                         cwork,dim_cs1  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')',', for ',dim_cs1,' soil pool(s)'
          CASE ('SI' )
            nz = 1
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1)")  &
                         cwork,nz  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')'
          CASE ( 'SO' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a)")   &
                         cwork,nsoil  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')',', for ',nsoil,' soil layers'
          CASE ( 'SN' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a,i3,a)")  &
                         cwork,ntiles*nsmax  &
                         ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                         ,'(',TRIM(varUnitsList(kvar)),')',', for ',nsmax  &
                         ,' snow layers on ',ntiles,' tiles'
          CASE ( 'TI' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a)")  &
                         cwork,ntiles  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')',', for ',ntiles,' tiles'
          CASE ( 'TY' )
            WRITE(outUnit(iout),"(a1,tr1,i3,tr1,a,tr1,a,tr1,a1,a,a1,tr1,a,i3,a)" )  &
                         cwork,ntype  &
                        ,ADJUSTL(varName(jvar)),TRIM(varDesc(jvar))  &
                        ,'(',TRIM(varUnitsList(kvar)),')',', for ',ntype,' surface types'
          CASE default
            WRITE(*,*)'ERROR: newOutFile: varType=',varType(jvar)
            WRITE(*,*)'No code for this varType.'
            WRITE(*,*)'Stopping in newOutFile.'
            STOP
        END SELECT
      ENDDO  !  variables

      WRITE(outUnit(iout),*) undefOut,' missing data value'

    ELSE  !  compressGridCall

!     When writing compression mapping, can use simpler headings.
!     Note that these should be kept consistent with what is actually written to
!     this file - currently done in subroutine init_out_map_compress.

      WRITE(outUnit(iout),"(a)") '1 variable in this file:'
      WRITE(outUnit(iout),"(a)") 'index (position in output vector to use at this location on grid)'
      WRITE(outUnit(iout),"(a,2i5)") 'Grid size (x,y)=',outGridNxy(iout,1),outGridNxy(iout,2)
      WRITE(outUnit(iout),"(a,2f8.2)") 'Lat/lon of point (1,1) (SW corner)=',outGridXY(iout,2),outGridXY(iout,1)

    ENDIF  !  compressGridCall

  ENDIF  !  outFormat ASCII

  END SUBROUTINE newOutFile


!###############################################################################
!###############################################################################
!###############################################################################
! subroutine newOutFileName
! Gets name of new output file and saves further details for use in
! writing a GrADS ctl file (e.g. period).

  SUBROUTINE newOutFileName( cper,compressGridCall,outputDate   &
                            ,outputTime,dateTimeString,periodUnitsLong  &
                            ,dsetLine,iout,tdef )

  USE inout, ONLY : formatBin,numMonth,outDataFile,outDir  &
                   ,outFilePer,compressGridFile,periodAnn,periodMon  &
                   ,outName,outPer,outPerNunits,periodOneFile,outFormat,runID,outTemplate
  USE misc_utils, ONLY : rm_path
  USE spin_mod, ONLY : ispin,nspin,spinUp
  USE switches, ONLY : l_360
  USE time_loc, ONLY : dateMainRun,dateSpin,timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    dateToBits,get_period,monthName3,secToSMH

  USE timeConst, ONLY :  &
!  imported scalar parameters
     iSecInDay,iSecInHour,iSecInMin
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout        &!  number of output profile
   ,outputDate  &!  date (yyyymmdd) of output
   ,outputTime   !  time of day (s) of output

  INTEGER ::   &!  local SCALARS
    day       &!  work
   ,fileDate  &!  date used for file name
   ,fileTime  &!  time used for file name
   ,hh,iper,mm,month,ss,tmpDate,tmpTime,year  ! work

  LOGICAL ,INTENT(in) ::  &!  in SCALARS
    compressGridCall   !  T means this call is for the supplementary data file
!                           that describes the mapping used when compressing output
!                           (e.g. the supplementary file referred to in the GrADS pdef entry).
!                           Note this is FALSE when just writing a ctl file that contains a pdef entry.

  CHARACTER(len=*), INTENT(out) ::  &!  out SCALARS
    cper            &!  used to hold file period for name
   ,dateTimeString  &!  date and time as yyyy-mm-dd hh:mm:ss
   ,dsetLine        &!  part of the line used for DSET in GrADS ctl file
   ,periodUnitsLong &!  longer form of periodUnits (q.v.)
   ,tdef             !  (most of ) the line used for TDEF in GrADS ctl file

  CHARACTER(len=2) ::  &!  local SCALARS
    periodUnits      !  suggested time unit for period

  CHARACTER(len=15) ::  &! local SCALARS
    charWork !  workspace. Longest possible is for a date=ddcccyyyy_hhmmH

  CHARACTER(len=3) ::   &!  local SCALARS
    cmonth3  !  3 character month name

  CHARACTER(len=8) ::  &!  local SCALARS
    tdefPer  !  data interval part of tdef

  CHARACTER(len=LEN(outDataFile)) ::  &!  local SCALARS
    fileName    !  name of file
!-------------------------------------------------------------------------------

! In the following,
!   per is the output period
!     (note this is the period of output, not necessarily that of the output files)
!   ext is the file extension or suffix
!
! All file names (unless compressGridCall=TRUE) begin with:          runID.profile
!              and end with:                                      .ext
! In most cases the output period is then indicated:         runID.profile.per
!   (although per is omitted if the period is not 'obviously'
!    related to one of year, month, day, hour, minutes)
! If all data to one file, name is:                          runID.profile.per.ext
! For other file periods (but see possible simplifications below):
!                                                            runID.profile.per.datetime.ext
!  where datetime is yyyymmdd_hhmm if numMonth=FALSE
!                    ddmonyyyy_hhmm if numMonth=TRUE
!
! If file period is yearly or monthly then "redundant" parts of datetime (shorter time intervals)
! are removed. Similarly, if file period is an integer number of days,hours or minutes, the shorter
! times are not in the name.
! According to outFilePer, spin-up cycles may be placed in files labelled 'spin01' etc
!      eg  runID.profile.per.datetime.spin01.ext

! If compressGridCall=TRUE, the file name is created in subroutine init_out_pdef.
!-------------------------------------------------------------------------------

! Get the 'basic' name.
  fileName = TRIM(outDir)//'/'//TRIM(runID) // '.' // TRIM(outName(iout))

! Deal with compressGridCall first, as most of the later stuff is not needed in that case.
  IF ( compressGridCall ) THEN
    fileName = TRIM(compressGridFile(iout))
    dsetLine = fileName
!   Remove path from dsetline and add"^" because the ctl file will be written
!   in the same directory as the data file.
    dsetLine = '^' // rm_path( dsetLine )
    cper = '-'
    tdef = '01jan2000 1dy'   !   tdef is largely unimportant in this case
    RETURN  !  nothing more to do
  ENDIF

! Get a representation of the output period (for name & GrADS tdef).

! As a general rule, the timestamp on output averages refers to the end of the averaging period,
! or the time of a snapshot.
! e.g. hourly aves, the average from 23H to midnight on last day of month gets timestamp 00H on
! 1st of following month, and would be written to next month's file (if monthly files).
! However, an exception is made for monthly and annual output as it's much easier to follow
! what's going on if the timestamp refers to the actual month/year.
! e.g. monthly averages, by the 'general rule' the Dec average is timestamped 00H 01Jan, and
! might be written to a file with January in name (if monthly files). So the 'general' rule is
! waived and the Dec ave is given a Dec timestamp. This gets rather awkward because a Dec
! snapshot (at 00Z 01Jan) is labelled as Dec, but if time is just specified as Dec, GrADS will
! make it 00Z 01Dec - which is rather confusing. Hence time for monthly data is given as ~23H
! on last day of month, and for annual data is Dec31.

! Set the date and time to be used for output.
  tmpDate = outputDate
  tmpTime = outputTime

  CALL DateToBits( tmpDate,day,month,year,l_360,'newOutFileName' )

  IF ( outPer(iout)==periodMon ) THEN

!   In the old days....I used to think that it was convenient that the monthly average for month M
!   appeared in a file with M in its name. However, a monthly snapshot for month M is output at
!   time 00Z on 1st of month M+1, and plotting this at 00Z on 1st of month M is rather confusing.
!     So....now we label monthly output with time 00Z on 1st of month M+1, with the understanding that
!     averages are backward averages ending at that time (so really average for month before).

!   Generally, monthly data for month M are output with time 00H on 1st of month M+1 (although output
!   is generated on a timestep that ends at that time, but starts 1*timestep earlier).
!   However, if end of run comes mid-month, we need to use 1st of next month for filenames.
    IF ( day/=1 .OR. tmpTime/=0 ) THEN
      month = month + 1
      IF ( month == 13 ) THEN
        month = 1
        year = year + 1
      ENDIF
      tmpDate = year*10000 + month*100 + 1
    ENDIF

  ELSEIF ( outPer(iout)==periodAnn ) THEN

!   Generally annual data for year Y are output at 00H on 1st Jan Y+1.
!   However, if end of run comes mid-year, we need to use 1st Jan of next year for filenames.
    IF ( month/=1 .OR. day/=1 .OR. tmpTime/=0 ) THEN
      year = year + 1
      month = 1
      tmpDate = year*10000 + month*100 + 1
    ENDIF

  ENDIF  !  outPer

  fileDate = tmpDate
  fileTime = tmpTime

  CALL DateToBits( fileDate,day,month,year,l_360,'newOutFileName' )
  CALL secToSMH ( fileTime,ss,mm,hh,'newOutFileName' )
  cmonth3 = monthName3( month )  !  get character month name
  WRITE(dateTimeString,"(i4.4,a1,a3,a1,i2.2,tr1,i2.2,a1,i2.2,a1,i2.2)")  &
       year,'-',cmonth3,'-',day,hh,':',mm,':',ss

  CALL get_period( outPer(iout),NINT(timeStep),outPerNunits(iout),periodUnits,periodUnitsLong )

  IF ( outPer(iout) == periodAnn ) THEN
!   Annual data.
    cper = 'ann'
    tdefPer = '1yr'
!   Include month and day.
    WRITE(tdef,"(i2.2,a3,i4.4,tr1,a)") day,cmonth3,year,TRIM(tdefper)

  ELSEIF ( outPer(iout) == periodMon ) THEN
!   Monthly data.
    cper = 'mon'
    tdefPer = '1mo'
!   Include month and day.
    WRITE(tdef,"(i2.2,a3,i4.4,tr1,a)") day,cmonth3,year,TRIM(tdefper)

  ELSE

!   Regular (not monthly or annual) data.
!   Note that the formats employed below will fail if periods are too large (e.g. if period
!   seems to be 2038 days!).

    IF ( periodUnits == 'dy' ) THEN

!     Period is a number of days.
      WRITE(cper,"(i2.2)") outPerNunits(iout)
      cper=TRIM(cper)//'d'
      WRITE(tdefPer,"(i2.2,'dy')") outPerNunits(iout)
      WRITE(tdef,"(i2.2,a3,i4.4,tr1,a)") day,cmonth3,year,TRIM(tdefper)

    ELSEIF ( periodUnits == 'hr' ) THEN

!     Period is a number of hours
      WRITE(cper,"(i2.2)") outPerNunits(iout)
      cper=TRIM(cper)//'h'
      WRITE(tdefPer,"(i2.2,'hr')") outPerNunits(iout)
      WRITE(tdef,"(i2.2,':',i2.2,'Z',i2.2,a3,i4.4,tr1,a)") hh,mm,day,cmonth3,year,TRIM(tdefper)

    ELSEIF ( periodUnits == 'mn' ) THEN

!     Period is a number of minutes
      WRITE(cper,"(i2.2)") outPerNunits(iout)
      cper=TRIM(cper)//'m'
      WRITE(tdefPer,"(i2.2,'mn')") outPerNunits(iout)
      WRITE(tdef,"(i2.2,':',i2.2,'Z',i2.2,a3,i4.4,tr1,a)") hh,mm,day,cmonth3,year,TRIM(tdefper)

    ELSE

!     Period is not any of the most likely (as far as GrADS is concerned at any rate).
!     If using GrADS, the run should have failed in initialisation, so we shouldn't be here!
!     Label with number of timesteps.
      WRITE(cper,"(i4.4)") outPerNunits(iout)
      cper=TRIM(cper)//'step'
      WRITE(tdefPer,"(i4.4,'step')") outPerNunits(iout)
      WRITE(tdef,"(i2.2,':',i2.2,'Z',i2.2,a3,i4.4,tr1,a)") hh,mm,day,cmonth3,year,TRIM(tdefper)

    ENDIF  !  periodUnits

  ENDIF  !  outPer

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Get (start of) file name.
! The file period affects how much of the date is included in the names.
! In some cases, the data period also alters name - NOTE eg monthly files we could always date
! as 1st of month (e.g. 19820101). This is needed for monthly data (see discussion elsewhere).
! However, for high frequency data (eg daily) we instead choose to omit the day from the file
! name (eg 198201), since this is in fact better if start of month is missing from file.
! It's all a bit of a faff, and mainly of importance if wanting to use GrADS with template ctl files!
!
!XX Could we use replaceTemplate to help with some of this? Would have to create a template name to start with.

  SELECT CASE ( outFilePer(iout) )
    CASE ( periodOneFile,-8,-7 )
!     All data (at least for each section of the run) go into one file. Will not be template type.
      fileName = TRIM(fileName)//'.'//TRIM(cper)
      dsetLine = TRIM(fileName)

!   NOTE: Now that month and day are included for annual and monthly files, code for those can
!         largely be combined. Not doing now...until this has been tested!

    CASE ( periodAnn )
!     Annual files.
      IF ( outPer(iout) < 0 ) THEN
!       Now include month and day.
        IF ( numMonth ) THEN
!         date: yyyymmdd
          WRITE(charWork,"(i4.4,2i2.2)") year,month,day
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.%y40101'
          fileName = TRIM(fileName) // '.' // charWork
        ELSE
!         date: ddcccyyyy
          WRITE(charWork,"(i2.2,a3,i4.4)") day,monthName3( month ),Year
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.01jan%y4'
          fileName = TRIM(fileName) // '.' // charWork
        ENDIF
      ELSE   !  outPer>=0
        WRITE(charWork,"(i4.4)") year
        fileName = TRIM(fileName)//'.'//TRIM(cper)
        IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%y4'
        fileName = TRIM(fileName)// '.' //charWork
      ENDIF

    CASE ( periodMon )
!     Monthly files.
      IF ( outPer(iout) < 0 ) THEN
!       Now include month and day.
        IF ( numMonth ) THEN
!         date: yyyymmdd
          WRITE(charWork,"(i4.4,2i2.2)") year,month,day
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.%y4%m201'
          fileName = TRIM(fileName) // '.' // charWork
        ELSE
!         date: ddcccyyyy
          WRITE(charWork,"(i2.2,a3,i4.4)") day,monthName3( month ),Year
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.01%mc%y4'
          fileName = TRIM(fileName) // '.' // charWork
        ENDIF
      ELSE  !  outPer>=0
        IF ( numMonth ) THEN
!         date: yyyymm
          WRITE(charWork,"(i4.4,i2.2)") year,month
          fileName = TRIM(fileName)//'.'//TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%y4%m2'
          fileName = TRIM(fileName) // '.' //charWork
        ELSE
!         date: cccyyyy
          WRITE(charWork,"(a3,i4.4)") monthName3( month ),year
          fileName = TRIM(fileName)//'.'//TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%mc%y4'
          fileName = TRIM(fileName) // '.' // charWork
        ENDIF
      ENDIF

    CASE default
!     Files at regular number of timesteps. Find out how to represent period.
!     Note that the default formats employed below will fail if periods are too large (e.g. if period
!     seems to be 2038 days!).

      CALL get_period( outFilePer(iout),NINT(timeStep),iper,periodUnits )

      IF ( periodUnits == 'dy' ) THEN

!       Period is a number of days.
        IF ( numMonth ) THEN
!         date: yyyymmdd
          WRITE(charWork,"(i4.4,2i2.2)") year,month,day
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.%y4%m2%d2'
          fileName = TRIM(fileName) // '.' // charWork
        ELSE
!         date: ddcccyyyy
          WRITE(charWork,"(i2.2,a3,i4.4)") day,monthName3( month ),Year
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine =  TRIM(fileName) // '.%d2%mc%y4'
          fileName = TRIM(fileName) // '.' // charWork
        ENDIF

      ELSEIF ( periodUnits == 'hr' ) THEN

!       Period is a number of hours.
        IF ( numMonth ) THEN
!         date: yyyymmdd_hh
          WRITE(charWork,"(i4.4,2i2.2,'_',i2.2,'H')") year,month,day,hh
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%y4%m2%d2_%h2H'
          fileName = TRIM(fileName) // '.' // charWork
        ELSE
!         date: ddcccyyyy_hh
          WRITE(charWork,"(i2.2,a3,i4.4,'_',i2.2,'H')") day,monthName3( month ),year,hh
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%h2H%d2%mc%y4H'
          fileName = TRIM(fileName) // '.' //charWork
        ENDIF

      ELSEIF ( periodUnits == 'mn' ) THEN

!       Period is a number of minutes
        IF ( numMonth ) THEN
!         date: yyyymmdd_hhmm
          WRITE(charWork,"(i2.2,a3,i4.4,2i2.2,'H')") year,month,day,INT(fileTime),hh,mm
          fileName = TRIM(fileName) // '.' // TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%h2%n2H%d2%mc%y4%m2%d2_%h2%n2H'
          fileName = TRIM(fileName) // '.' // charWork
        ELSE
!         date: ddcccyyyy_hhmm
          WRITE(charWork,"(i2.2,a3,i4.4,'_',2i2.2,'H')") day,monthName3( month ),year,hh,mm
          fileName = TRIM(fileName)// '.' //TRIM(cper)
          IF ( outTemplate(iout) ) dsetLine = TRIM(fileName) // '.%d2%mc%y4_%h2%n2H'
          fileName = TRIM(fileName) // '.' // charWork
        ENDIF

      ELSE

!       Period is not any of the most likely (as far as GrADS is concerned at any rate).
!       Don't include in file name.

      ENDIF  !  periodUnits

  END SELECT  !  outFilePer
!-----------------------------------------------------------------------

  IF ( .NOT. outTemplate(iout) ) dsetLine = TRIM(fileName)

! When in spinup, show this in file names. Don't do this if main run dates follow immediately
! from spinup over a single cycle, or if all data are to go in the one file.
  IF ( spinUp .AND. .NOT.(dateMainRun(1)==dateSpin(2) .AND. nspin==1) .AND.  &
               outFilePer(iout)/=periodOneFile ) THEN
    IF ( outFilePer(iout) /= -8 ) THEN
      WRITE(charWork,"(i2.2)") ispin
      fileName = TRIM(fileName)//'.spin'//charWork
      dsetLine = TRIM(dsetLine)//'.spin'//charWork
    ELSE
!     outFilePer=-8. Label is "spin", with no number.
      fileName = TRIM(fileName)//'.spin'
      dsetLine = TRIM(dsetLine)//'.spin'
    ENDIF
  ENDIF

! Add extension. For now, name binary files with 'gra' (rather than 'bin').
  SELECT CASE ( outFormat )
    CASE ( formatBin )
      fileName = TRIM(fileName) // '.gra'
      dsetLine = TRIM(dsetLine) // '.gra'
    CASE default
      fileName = TRIM(fileName) // '.' // TRIM(outFormat)
      dsetLine = TRIM(dsetLine) // '.' // TRIM(outFormat)
  END SELECT

! Save name of data file.
  outDataFile(iout) = fileName

! Remove path from dsetline and add"^" because the ctl file will be written
! in the same directory as the data file..
  dsetLine = '^' // rm_path( dsetLine )

  END SUBROUTINE newOutFileName
!###############################################################################
!###############################################################################
!###############################################################################

! subroutine writeGradsCtl

! Write a GrADS ctl file.
!
! For simplicity, only one template ctl file is allowed per profile (apart from
! dealing with any spinup files).
!-------------------------------------------------------------------------------

  SUBROUTINE writeGradsCtl( iout,nvarCtl,dsetLine,optionsLine,titleLine,zdefLine  &
                           ,cper,tdef,compressGridCall  &
                           ,annot,tdefAnnot,varCtlLine )

  USE file_utils, ONLY : closeFile,fileUnit,openFile
  USE inout, ONLY : formatAsc,formatNc,outCtlFile,outDataFile,outDir,outFormat,outName  &
                   ,outStatus,compressGridFile,ntCtl,outEndian  &
                   ,outGridDxy,outGridXY,outGridNxy,periodOneFile,pointsOut   &
                   ,runID,outTemplate,undefOut,outCompress

  USE misc_utils, ONLY : &
!  imported procedures
     rm_path

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in
     ispin,spinUp

  USE time_loc, ONLY :  &
!  imported arrays with intent(in)
     dateMainRun,dateSpin

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: iout     !  # of output profile
  INTEGER, INTENT(in) :: nvarCtl  !  number of variables, as represented in ctl file

  LOGICAL, INTENT(in) ::  &
    compressGridCall   !  T means this call is for the supplementary data file
!                           that describes the mapping used when compressing output
!                           (e.g. the supplementary file referred to in the GrADS pdef entry).
!                           Note this is FALSE when just writing a ctl file that
!                           contains a pdef entry.  work

  CHARACTER(len=*), INTENT(in) ::   &! in SCALARS
    cper          &!  output period as used in dataFileName
   ,dsetLine      &!  file name as used for DSET
   ,optionsLine   &!  most of the line used for GrADS options line
   ,tdef          &!  (most of) the line used for GrADS tdef
   ,titleLine     &!  the line used for GrADS title
   ,zdefLine       !  most of the line used for GrADS zdef

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  CHARACTER(len=*), INTENT(in) :: annot(:)     !  annotation (metadata)
  CHARACTER(len=*), INTENT(in) :: tdefAnnot(:)   !  annotation (metadata)
  CHARACTER(len=*), INTENT(in) :: varCtlLine(:)     !  descriptions of variables

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::   &
    i,ivar  &!  counters/work
   ,nxVal   &!
   ,nyVal   &!
   ,unit    &!  unit number
   ,xlen     !  used to hold a len_trim value.

  REAL :: &
    dlat,dlon,lat1,lon1   !  work

  CHARACTER(len=LEN(outDataFile)) ::  &
    fileName   !  work

  CHARACTER(len=LEN(compressGridFile)) ::  &
    outGridNxyame   !  work

  CHARACTER(len=200) ::  &
    cwork                                  &!  work
   ,xdefLine,xdefLine2,ydefLine,ydefLine2   !  work

!-------------------------------------------------------------------------------
! Get name for ctl file.
!-------------------------------------------------------------------------------
  IF ( .NOT. outTemplate(iout) .OR. compressGridCall ) THEN
!   This is the same as the data file name except for the extension.
    fileName = outDataFile(iout)
    IF ( compressGridCall ) fileName = compressGridFile(iout)
!   Check!
    xlen = LEN_TRIM( fileName )
    IF ( fileName(xlen-2:xlen) == 'gra' ) THEN
      i = 3
    ELSEIF ( fileName(xlen-1:xlen) == 'nc' ) THEN
      i = 2
    ELSE
      WRITE(*,*)'writeGradsCtl: data file extension is not gra or nc. Unexpected.'
      STOP
    ENDIF
    outCtlFile(iout) = filename(1:xlen-i)//'ctl'
  ELSE
!   Template option.
!   Only one template ctl is generated per profile - except for spin up.
    IF ( spinUp .AND. dateSpin(2)/=dateMainRun(1) ) THEN
      WRITE(cwork,"(i2.2)") ispin
      outCtlFile(iout) = TRIM(outDir) // '/' // TRIM(runID) // '.' // TRIM(outName(iout)) //  &
        '.' // TRIM(cper) // '.spin' // TRIM(cwork) // '.ctl'
    ELSE
      outCtlFile(iout) = TRIM(outDir) // '/' // TRIM(runID) // '.' // TRIM(outName(iout)) //  &
        '.' // TRIM(cper) // '.ctl'
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Get values describing grid, for use in xdef and ydef lines.
!-------------------------------------------------------------------------------
  nxVal = outGridNxy(iout,1)
  nyVal = outGridNxy(iout,2)
  lon1 = outGridXY(iout,1)
  dlon = outGridDxy(iout,1)
  lat1 = outGridXY(iout,2)
  dlat = outGridDxy(iout,2)

!-------------------------------------------------------------------------------
! Get values for xdef and ydef.
!-------------------------------------------------------------------------------
  IF ( compressGridCall .OR. outCompress(iout) ) THEN
!   We are writing a ctl file either for the supplementary file, or a ctl file that refers
!   to a supplementary file (i.e. has a pdef entry).
!   Grid extent is the larger uncompressed grid - set above.
    WRITE(xdefline,"('xdef ',i6,' linear',2(tr1,f9.4))") nxVal,lon1,dlon
    WRITE(ydefline,"('ydef ',i6,' linear',2(tr1,f9.4))") nyVal,lat1,dlat
  ENDIF

  IF ( .NOT.outCompress(iout) .OR. (outCompress(iout).AND..NOT.compressGridCall) ) THEN
!   We are writing a ctl file for an uncompressed output profile.
!   Also used to get values for annotation if is a compressed profile and we are not
!   writing a ctl for the supplementary [pdef] file
!   If values are unreasonable (would not plot well), make up something reasonable(ish).
!   e.g. large vector of 1 degree data with dlon=1 would plot with enormous range of longitude.
    IF ( dlon*nxVal > 360.0 ) dlon = 360.0 / nxVal
    IF ( dlat*nyVal > 180.0 ) dlat = 180.0 / nyVal
    dlon = MAX( dlon, 0.0001 )
    dlat = MAX( dlat, 0.0001 )

!   Create the lines for the ctl file.
    IF ( .NOT. outCompress(iout) ) THEN
      WRITE(xdefline,"('xdef ',i6,' linear',2(tr1,f9.4))") nxVal,lon1,dlon
      WRITE(ydefline,"('ydef ',i6,' linear',2(tr1,f9.4))") nyVal,lat1,dlat
    ELSE
!     Create comment lines.
      WRITE(xdefline2,"('# xdef ',i6,' linear',2(tr1,f9.4))") pointsOut(iout),lon1,dlon
      WRITE(ydefline2,"('# ydef      1 linear',2(tr1,f9.4))") lat1,dlat
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Get unit number and open file.
!-------------------------------------------------------------------------------
  unit = fileUnit( formatAsc )
  CALL openFile(1,.FALSE.,unit,'write',formatAsc,outCtlFile(iout),outStatus)

!-------------------------------------------------------------------------------
! Start writing the ctl file (dset, options, dtype, title, undef).
!-------------------------------------------------------------------------------
  WRITE(unit,"('dset ',a)") TRIM(dsetLine)
  WRITE(unit,"('title ',a)") TRIM(titleLine)
  WRITE(unit,"('undef ',es10.2 )") undefOut

  IF ( LEN_TRIM(optionsLine) > 0 ) WRITE(unit,"(a,tr1,a)") 'options',TRIM(optionsLine)
  IF ( .NOT.compressGridCall .AND. outFormat==formatNc ) WRITE(unit,"(a)") 'dtype netcdf'

!-------------------------------------------------------------------------------
!  Write a pdef line.
!-------------------------------------------------------------------------------
  IF ( outCompress(iout) .AND. .NOT.compressGridCall ) THEN
    cwork = 'binary-big'
    IF ( outEndian == 'little_endian' ) cwork = 'binary-little'
    outGridNxyame = compressGridFile(iout)
!   Remove path from and add"^" because the pdef data file will be written
!   in the same directory as this ctl file.
    outGridNxyame = '^' // rm_path( outGridNxyame )
    WRITE(unit,"('pdef ',i7,' 1 file 1 stream ',a,tr1,a)") pointsOut(iout)  &
                          ,TRIM(cwork),TRIM(outGridNxyame)
!   Write a comment and the xdef/ydef that would be used without pdef.
    WRITE(unit,"(a)") '# The following xdef and ydef would be used WITHOUT the pdef line above.'
    WRITE(unit,"(a)") TRIM(xdefLine2)
    WRITE(unit,"(a)") TRIM(ydefLine2)
  ENDIF

!-------------------------------------------------------------------------------
! Write grid description (xdef, ydef, zdef, tdef).
!-------------------------------------------------------------------------------
  WRITE(unit,"(a)") TRIM(xdefLine)
  WRITE(unit,"(a)") TRIM(ydefLine)
  WRITE(unit,"(a)") TRIM(zdefLine)

! Include any annotation for tdef.
  DO i=1,nmax
    IF ( LEN_TRIM(tdefAnnot(i)) == 0 ) EXIT
    WRITE(unit,"(a)") TRIM( tdefAnnot(i) )
  ENDDO

  WRITE(unit,"('tdef ',i7,' linear ',a)") ntCtl(iout),TRIM(tdef)

!-------------------------------------------------------------------------------
! Write general annotation.
!-------------------------------------------------------------------------------
  DO i=1,nmax
    IF ( LEN_TRIM(annot(i)) == 0 ) EXIT
    WRITE(unit,"(a)") TRIM( annot(i) )
  ENDDO

!-------------------------------------------------------------------------------
! Write details of each variable.
!-------------------------------------------------------------------------------
  WRITE(unit,"('vars ',i3)") nvarCtl

  DO ivar=1,nvarCtl
    WRITE(unit,"(a)") TRIM( varCtlLine(ivar) )
  ENDDO

  WRITE(unit,"('endvars')")

!-------------------------------------------------------------------------------
! Close file.
!-------------------------------------------------------------------------------
  CALL closeFile( unit,formatAsc )

  END SUBROUTINE writeGradsCtl
!################################################################################
!################################################################################
!################################################################################

! subroutine loadOut
! Fill the output variables for a particular output profile.

  SUBROUTINE loadOut( doAccum,doTmean,iout,nlev,nxy,outputDate,outputTime,outvalWrite )

  USE ancil_info, ONLY :  &
!   imported scalars with intent(in)
      dim_cs1,land_pts,lice_pts,nsmax,nsoil=>sm_levels,ntiles  &
     ,nx=>row_length,ny=>rows,soil_pts,tile_pts &
!   imported arrays with intent(in)
     ,frac,land_index,lice_index,soil_index,tile_index

  USE c_0_dg_c, ONLY :  &
!  imported scalar parameters
     tm   !   temperature at which fresh water freezes and ice melts (K)

  USE csigma, ONLY :  &
!  imported scalar parameters
    sbcon

  USE fluxes, ONLY :  &
!   imported arrays with intent(in)
             alb_tile,ecan,ecan_tile,ei,ei_tile,esoil,esoil_tile,ext  &
            ,fqw_1,fqw_tile,fsmc,ftl_1,ftl_tile   &
            ,hf_snow_melt,land_albedo,latent_heat,le_tile,melt_tile  &
            ,radnet_tile,emis_tile,snomlt_sub_htf,snomlt_surf_htf,snow_melt  &
            ,sub_surf_roff,surf_ht_flux,surf_roff  &
            ,taux_1,tauy_1,tot_tfall,tstar,surf_htf_tile  &
            ,surf_ht_store,anthrop_heat

  USE forcing, ONLY :  &
!   imported arrays with intent(in)
              con_rain,con_snow,ls_rain,ls_snow   &
             ,lw_down,pstar,qw_1,sw_down,tl_1,u_1,v_1

  USE grid_utils, ONLY :  &
!  imported procedures
     getXYPos,reverseCols

  USE inout, ONLY :   &
!   imported scalar parameters
     periodMon  &
!   imported scalars with intent(in)
    ,outLenWrite  &
!   imported arrays with intent(in)
    ,mapOut,mapOutLand,ntOutPer,nvarOut,outPer  &
    ,outStep,outStepSamPer,outWarnUnder,pointsOut,pointsOutLand,taccumVar,tmeanVar  &
    ,varNameList,varNum,varPos,varStartPos,varType,zrevOutSnow,zrevOutSoil  &
!   imported arrays with intent(inout)
    ,outval

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     lake,npft,ntype

  USE offline_diag, ONLY :  &
!  imported arrays with intent(in)
     ciDiag,GstomDiag,RdcDiag,RFlowDiag,roffInfDiag,RRunDiag  &
    ,snowGMeltDiag,wfluxDiag,wfluxSfcDIag

  USE p_s_parms, ONLY :   &
!   imported arrays with intent(in)
     b,catch,cosz,hcap,hcon  &
    ,satcon,sathh,smvccl,smvcst,smvcwt,sthf,sthu,z0_tile

  USE prognostics, ONLY :   &
!   imported arrays with intent(in)
     canht_ft,canopy,canopy_gb,cs,gc,gs,lai,nsnow,rgrain,rgrainL,rho_snow_grnd  &
    ,routeStore,sice,sliq,smcl,snowDepth,snow_mass,snow_grnd,snow_tile  &
    ,t_soil,tsnow,tstar_tile

  USE route_mod, ONLY :  &
!  imported array parameters
     flowDirDelta  &
!  imported scalars with intent(in)
    ,npRoute,nxRoute,nyRoute  &
!  imported arrays with intent(in)
    ,flowDirSet,routeIndex,routeNext

  USE screen, ONLY :  &
!   imported arrays with intent(in)
     q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m

  USE snow_param, ONLY :  &
!  imported arrays with intent(in)
     canSnowTile,ds

  USE surf_param, ONLY : diff_frac

  USE soil_param, ONLY : zsmc

  USE switches, ONLY :  &
!   imported scalars with intent(in)
     l_aggregate,routeOnly

  USE time_loc, ONLY :  &
!   imported scalars with intent(in)
     timestep  &
!   imported  arrays with intent(in)
    ,latitude,longitude

  USE top_pdm, ONLY :   &
!  imported arrays with intent(in)
     drain,dun_roff,fch4_wetl,fsat,fwetl,qbase,qbase_zw,sthzw,zw

  USE trifctl, ONLY :   &
!   imported arrays with intent(in)
     c_veg,cv,g_leaf,g_leaf_day,g_leaf_dr_out,g_leaf_phen,gpp,gpp_ft   &
    ,lai_phen,lit_c,lit_c_mn,npp,npp_dr_out,npp_ft,resp_p,resp_p_ft,resp_s  &
    ,resp_s_dr_out,resp_w_dr_out,resp_w_ft

  USE ozone_vars

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) ::   &
    iout   &!  number of output profile to do
   ,nlev   &!  number of levels needed for workspace (enough for variable with largest # of levels)
   ,nxy    &!  number of x,y points needed for workspace (enough for variable with largest nx*ny)
   ,outputDate  &!  date (yyyymmdd) of output
   ,outputTime   !  time (s of day) of output

  LOGICAL, INTENT(in) ::  &
    doAccum    &!  T if values are to be added to accumulation over time
   ,doTmean     !  T if time average is to be calculated


!-------------------------------------------------------------------------------
! Array arguments with intent(out)
!-------------------------------------------------------------------------------
  REAL, INTENT(out) :: outvalWrite(outlenWrite)   !  data as written to file.
!                             Only used on calls when output to file is required.
!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  REAL, PARAMETER :: snowThresh = 0.005  !  threshold depth for snow on ground (m),
!              below which derived snow diagnostics are not calculated - used to avoid
!              small snow depths

  REAL, PARAMETER :: epsilonReal = EPSILON( tstar )   !  epsilon of a REAL variable

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    dx,dy     &!  work
   ,ipos      &!  position in outval of first value of current level of a
!                 variable (i.e. current field)
   ,iposOne   &!  position in outval of 1st value of 1st level of a variable
   ,iposWrite &!  position in outvalWrite of
   ,i,ivar,ixn,ixr,iyn,iyr,iz,j,jvar,kvar,l,n,p  &!  loop counters/work
   ,np        &!  number of points in a variable
   ,npOut     &!  number of points to be output
   ,nz         !  number of levels in a variable

  REAL :: sumFrac   !  holds accumulated frac

  LOGICAL ::  &
    found      &!  T when variable found in list
   ,undefTave   !  T if time-average is to be set to undefined (missing data) value
!                   because of missing data (missing times in average)

  CHARACTER(len=LEN(varNameList)) :: var   !  (short) variable name

  CHARACTER(len=7) :: periodStr   !  work

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    countTile(nxy,ntiles)      &!  Counter of values for time-mean tile variables.
!               Set to outStep where tile_frac>0 - this asumes that all tiles have had a
!               value at all output times, in particular that active tiles have not changed.
!               Note: tile variables are dimensioned e.g. land_pts,ntiles, but only have "useful"
!               values at point where tile_frac>0. countTile is usd to mark those points.
   ,countTileSnap(nxy,ntiles)  &!  Counter of values for snapshot tile variables.
!               Set to 1 where tile_frac>0.
   ,pointsUse(nxy)             &!  a list of the points in the model grid that are to be output
   ,xcount(nxy)                &!  counter
   ,xcountTile(nxy,ntiles)     &!  counter for tile variables
   ,xxcount(nxy,nlev)           !  counter

  REAL ::  &
    tile_frac(land_pts,ntiles) &!  tile fractions
   ,tmpval(nxy,nlev)           &!  workspace
   ,tmpval2(land_pts,nsoil)     !  workspace - enough for a soil layer variable

!--------------------------------------------------------------------------------

! Start by assuming that there are enough data for time averages.
  undefTave = .FALSE.

! Check we have sufficient times of data if calculating time average (or accumulation).
  IF ( doTmean ) THEN
    IF ( outPer(iout) > 0 ) THEN
!     If any times of data are missing, don't calculate averages or accumulations.
      IF ( outStepSamPer(iout) /= ntOutPer(iout) ) THEN
!       Set flags indicating that the time-average is to be set to undefined (missing)
!       value, and error flagged.
        undefTave = .TRUE.
        outWarnUnder = .TRUE.
        WRITE(*,"(a,tr1,i7,tr1,a,i2,a,i8,i6)")  &
            'WARNING: outWarnUnder: insufficient data to calculate average over'  &
           ,ntOutPer(iout)*NINT(timestep),'seconds. iout=',iout  &
           ,' date,time=',outputDate,outputTime
      ENDIF
    ELSE
!     "Special" periods.
!     Check if >=90% of expected data times were added to accumulation.
!     See note in INIT_OUT (where ntOutPer is set) regarding the level at which the
!     threshold for number of data
!     is set (e.g. 90%). In particular, if it is too high, some monthly or annual
!      averages may not be calculated,
!     because of the approximation of the length of the month or year.
      IF ( REAL(outStepSamPer(iout))/REAL(ntOutPer(iout)) < 0.9 ) THEN
        periodStr = 'annual'
        IF ( outper(iout)==periodMon ) periodStr='monthly'
!       Set flags indicating that the time-average is to be set to undefined
!       (missing) value, and error flagged.
        undefTave = .TRUE.
        outWarnUnder = .TRUE.
        WRITE(*,"(a,tr1,a,tr1,a,i2,a,i8,i6)")  &
            'WARNING: outWarnUnder: insufficient data to calculate'  &
           ,TRIM(periodStr),'average. iout=',iout  &
           ,' date,time=',outputDate,outputTime
      ENDIF
    ENDIF  !  outPer
  ENDIF    !  doTmean

!-------------------------------------------------------------------------------
! Calculate tile fractions.

  IF ( .NOT. routeOnly ) THEN
    tile_frac(:,:) = 0.0
    IF ( l_aggregate ) THEN
      tile_frac(:,1) = 1.0
    ELSE
      DO n=1,ntype
        DO j=1,tile_pts(n)
          i = tile_index(j,n)
          tile_frac(i,n) = frac(i,n)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
!-------------------------------------------------------------------------------

! Some/all of these counters could be set within case( varType ) later. But not done yet!

! Counters that are not specific to an output profile
! (e.g. countTileSnap) could be calculated for iout=1 and then saved
! for reuse by iout>1. Currently recalculated for every iout.

  xxcount(:,:) = 0

! Get mask showing tiles with tile_frac>0.
  countTile(:,:) = 0    !    reinit for now, but may later want to accumulate
  countTileSnap(:,:) = 0
  IF ( .NOT. routeOnly ) THEN
    DO n=1,ntiles
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        countTile(l,n) = outStep(iout)  !  assuming always defined (not intermittent)
        countTileSnap(l,n) = 1
      ENDDO
    ENDDO
  ENDIF  !  routeOnly

!-------------------------------------------------------------------------------

! Initialise.
  iposWrite = 0

! Loop over variables, loading each.
  DO ivar=1,nvarOut(iout)

    found = .FALSE.
    jvar = varPos(iout,ivar)  !  position in output.
    kvar = varNum(jvar)       !  number in list of available vars.
    var = varNameList(kvar)
    iposOne = varStartPos(jvar)  !  position of first value (of first level) in output

!   Get counter that is passed to loadOutVarReal - this contains number of valid times in
!   time accumulation for time average, or indicates tile points with frac>0 for tile variables.
!   Counter depends upon whether snapshot/mean/accumulation.
    IF ( tmeanVar(varPos(iout,ivar)) ) THEN
      xcount(:) = outStepSamPer(iout)  !  assuming that all points have had a value at all times
      xcountTile(:,:) = countTile(:,:)
    ELSE  !  snapshot or time-accumulation variable
      xcount(:) = 1   !   assuming that all points have a value - NOT suitable for time-accumulation of intermittent variable
       xcountTile(:,:) = countTileSnap(:,:)
    ENDIF

!   Reinitialise workspace.
    tmpval(:,:) = 0.0

!   NOTE: This procedure extracts a copy of each field requested at all points, and then
!   selects the points that are to be output. This is rather inefficient if a small subset of points
!   is to be output, and the diagnostic needs to be derived, since it is derived at all points,
!   then most are discarded.
!   However, alternative approaches probably also have their drawbacks....!
!   For derived fields, we could still use a "full-size" array, but only calculate values at
!   points that are selected for output - at worst would require a loop rather than a whole-array
!   operation.

!   Get number of points (np) and number of levels (nz) for this variable (over full domain).
!   Also get the number of points to be output (npOut), the points to be output (pointsUse),
!   and appropriate counters and mapping.

!   Load the whole field into workspace.
!   Note that the variable names (var) given here must match those in init_out_varlist,
!   including the case.

    SELECT CASE ( varType(jvar) )

!###############################################################################
!###############################################################################
!###############################################################################
!     Variables that have a single value at LAND gridpoints only.
!     varType is 'LA'

      CASE ( 'LA' )
        np = land_pts
        nz = 1
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(:,1) = xcount(:)

!-------------------------------------------------------------------------------
        IF( var == 'albedoLand' ) THEN
          found = .TRUE.
!         Calculate the albedo as used in subroutine control when calculating
!         the net shortwave on tiles. Obviously this diagnostic is only useful
!         if the form used here and in control are the same!
!         Here we take the average of diffuse albedos in VIS and NIR.
          tmpval(1:np,1) = gbmTileDiag( 0.5*(alb_tile(:,:,2)+alb_tile(:,:,4))  &
                                      , tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'canopy' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = canopy_gb(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'cs' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = SUM( cs(:,:),2 )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'cv' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = cv(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'depthFrozen' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = frozenDepth( t_soil(:,:),'frozen' )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'depthUnfrozen' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = frozenDepth( t_soil(:,:),'unfrozen' )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'drain' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = drain(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'elake' ) THEN
          found = .TRUE.
!         Note this only makes sense if l_aggregate=FALSE.
          IF ( lake > 0 ) tmpval(1:np,1) = fqw_tile(:,lake) * tile_frac(:,lake)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fch4_wetl' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = fch4_wetl(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fsat' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = fsat(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fsmc' ) THEN
          found = .TRUE.
!         Calculate gridbox mean over PFTs.
          DO i=1,land_pts
            sumFrac = 0.0
            DO n=1,npft
              tmpval(i,1) = tmpval(i,1) + fsmc(i,n) * frac(i,n)
              sumFrac = sumFrac + frac(i,n)
            ENDDO
            IF ( sumFrac > epsilonReal ) THEN
              tmpval(i,1) = tmpval(i,1) / sumFrac
            ELSE
!             No veg here. For now....just set to some impossible value.
              tmpval(i,1) = -9.0
            ENDIF
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fwetl' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = fwetl(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gpp' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gpp(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gpp' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gpp(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gs' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gs(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'hfSnowMelt' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = hf_snow_melt(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'landIndex' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = REAL( land_index(:) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'liceIndex' ) THEN
          found = .TRUE.
          DO l=1,lice_pts
            tmpval(lice_index(l),1) = REAL( lice_index(l) )
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'litCMn' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = lit_c_mn(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'emis' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gbmTileDiag( emis_tile, tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'LWnet' ) THEN
          found = .TRUE.
!         First get gridbox upwards longwave.
          tmpval(1:np,1) = sbcon * gbmTileDiag( emis_tile *  &
               tstar_tile*tstar_tile*tstar_tile*tstar_tile,tile_frac )
!         Now get net flux.
          DO l=1,land_pts
            j = ( land_index(l)-1 ) / nx + 1
            i = land_index(l) - (j-1)*nx
            tmpval(l,1) = lw_down(i,j) - tmpval(l,1)
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'LWup' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = sbcon * gbmTileDiag( emis_tile *  &
               tstar_tile*tstar_tile*tstar_tile*tstar_tile,tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'npp' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = npp(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'qbase' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = qbase(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'qbase_zw' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = qbase_zw(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'radnet' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gbmTileDiag( radnet_tile, tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respP' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = resp_p(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respS' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = SUM( resp_s(:,:),2 )
!-------------------------------------------------------------------------------
! HACK: We only output the total respiration for now
        ELSEIF( var == 'respSDrOut' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = resp_s_dr_out(:,5)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'runoff' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = sub_surf_roff(:) + surf_roff(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'sat_excess_roff' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = dun_roff(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'smcAvailTop' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = soilMoistDiag( zsmc,'avail',sthu )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'smcAvailTot' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = soilMoistDiag( -1.0,'avail',sthu )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'smcTot' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = soilMoistDiag( -1.0,'total',sthf+sthu )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snomltSubHtf' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = snomlt_sub_htf(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowCan' ) THEN
          found = .TRUE.
!         Only include tiles where canopy snow model is used.
          tmpval(1:np,1) = gbmTileDiag( snow_tile, tile_frac  &
                                       ,canSnowTile )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowDepth' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = gbmTileDiag( snowDepth, tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowFrac' ) THEN
          found = .TRUE.
!         Sum tile_frac over tiles with snow.
          tmpval(1:np,1) = SUM( tile_frac, 2  &
                               ,snow_tile+snow_grnd>snowThresh )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowFracAlb' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = get_fsnow( tile_frac )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowGrCan' ) THEN
          found = .TRUE.
!         Only include tiles where canopy snow model is used.
          tmpval(1:np,1) = gbmTileDiag( snow_grnd, tile_frac  &
                                       ,canSnowTile )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowIceTot' ) THEN
          found = .TRUE.
          DO n=1,ntiles
            DO j=1,tile_pts(n)
              i = tile_index(j,n)
              tmpval(i,1) = tmpval(i,1) + SUM( sice(i,n,1:nsnow(i,n)) )  &
                               * tile_frac(i,n)
            ENDDO
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowLiqTot' ) THEN
          found = .TRUE.
          DO n=1,ntiles
            DO j=1,tile_pts(n)
              i = tile_index(j,n)
              tmpval(i,1) = tmpval(i,1) + SUM( sliq(i,n,1:nsnow(i,n)) )  &
                               * tile_frac(i,n)
            ENDDO
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowMelt' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = snow_melt(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'soilIndex' ) THEN
          found = .TRUE.
          DO l=1,soil_pts
            tmpval(soil_index(l),1) = REAL( soil_index(l) )
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'sthzw' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = sthzw(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'subSurfRoff' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = sub_surf_roff(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'surfRoff' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = surf_roff(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'surfRoffInf' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = roffInfDiag(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'swetLiqTot' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = soilMoistDiag( -1.0,'total',sthu )   !  column unfrozen smc
!         Divide by saturated moisture content.
          tmpval2(1:np,1:nsoil) = 1.0
          tmpval(1:np,1) = tmpval(1:np,1) / soilMoistDiag( -1.0,'total',tmpval2(1:np,1:nsoil) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'swetTot' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = soilMoistDiag( -1.0,'total',sthf+sthu )   !  column smc
!         Divide by saturated moisture content.
          tmpval2(1:np,1:nsoil) = 1.0
          tmpval(1:np,1) = tmpval(1:np,1) / soilMoistDiag( -1.0,'total',tmpval2(1:np,1:nsoil) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'SWnet' ) THEN
          found = .TRUE.
!         Calculate the albedo as used in subroutine control when calculating
!         the net shortwave on tiles. Obviously this diagnostic is only useful
!         if the form used here and in control are the same!
!         Here we take the average of diffuse albedos in VIS and NIR.
          DO l=1,land_pts
            j = ( land_index(l)-1 ) / nx + 1
            i = land_index(l) - (j-1)*nx
            tmpval(l,1) = ( 1.0 -  &
                  0.5 * ( land_albedo(i,j,2) + land_albedo(i,j,4) )  &
                          ) * sw_down(i,j)
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'tfall' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = tot_tfall(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'trad' ) THEN
          found = .TRUE.
!         Assuming emissivity=1.
          tmpval(1:np,1) = ( gbmTileDiag(  &
               tstar_tile*tstar_tile*tstar_tile*tstar_tile,tile_frac ) ) ** 0.25
!-------------------------------------------------------------------------------
        ELSEIF( var == 'wFluxSfc' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = wfluxSfcDiag(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'zw' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = zw(:)
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################

!     Vegetation (PFT) diagnostics (on land grid).
!     varType is 'PF'

      CASE ( 'PF' )
        np = land_pts
        nz = npft
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
!       Use the tile counters, so that where tile_frac=0 we get undef output.
!       Note that this assumes that the PFTs are tiles 1:npft - as they are at present.
!       For the (generally not well dealt with) case of l_aggregate+ntype>1, these
!       variables are still PFT variables, but we only have a counter variable
!       for tiles - so spread that.
        IF ( ntiles > 1 ) THEN
          xxcount(1:np,1:nz) = xcountTile(1:np,1:nz)
        ELSE
!         Get npft copies of counter.
          xxcount(1:np,1:nz) = SPREAD( xcountTile(1:np,1),2,nz )
        ENDIF

!--------------------------------------------------------------------------------
        IF( var == 'cVegP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = c_veg(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'canhtP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = canht_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'ciP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = ciDiag(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fsmcP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = fsmc(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gLeafP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = g_leaf(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gLeafDayP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = g_leaf_day(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gLeafDrOutP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = g_leaf_dr_out(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gLeafPhenP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = g_leaf_phen(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gppP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = gpp_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'gstomP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = gstomDiag(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'laiP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = lai(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'laiPhenP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = lai_phen(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'litCP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = lit_c(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'nppDrOutP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = npp_dr_out(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'nppP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = npp_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rdcP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = rdcDiag(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respPP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = resp_p_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respWDrOutP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = resp_w_dr_out(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respWP' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = resp_w_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fluxO3Stom' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = flux_o3_ft(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'o3ExpFac' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = fo3_ft(:,:)
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################

!     Variables on the channel routing grid.
!     varType is 'RG'

      CASE ( 'RG' )
        np = nxRoute*nyRoute
        nz = 1
        npOut = pointsOut(iout)
        pointsUse(1:npOut) = mapOut(iout,1:npOut,1)
        xxcount(:,1) = xcount(:)
!-------------------------------------------------------------------------------
        IF( var == 'flowdir' ) THEN
          found = .TRUE.
!         Initialise. This value is retained at sea points.
          tmpval(1:np,1:nz) = flowDirSet(0)
!         Loop over "active" routing points.
          DO n=1,npRoute
!           The location in tmpval of this point is given by routeIndex,
!           since both refer to a space of size nxRoute*nyRoute points.
            IF ( routeNext(n) > 0 ) THEN
!             Back-calculate flow direction from index of next downstream point.
!             Work out where this point lies in the routing grid.
              CALL getXYPos( routeIndex(n),nxRoute,nyRoute,ixr,iyr )
!             Get location in grid of downstream point.
              CALL getXYPos( routeNext(n),nxRoute,nyRoute,ixn,iyn )
              dx = ixn - ixr
              dy = iyn - iyr
!             Deal with cyclic boundary conditions.
              IF ( ABS(dx) > 1 ) dx = SIGN(1,-dx)
              DO p=1,8
                IF ( dx==flowDirDelta(p,1) .AND. dy==flowDirDelta(p,2) ) THEN
                  tmpval(routeIndex(n),1) = REAL( flowDirSet(p) )
                  EXIT
                ENDIF
              ENDDO
            ELSE
!             Next point is off grid, or there is no next point (no defined flow).
!             routeNext is -1*flow direction index.
              ixn = ABS( routeNext(n) )
              tmpval(routeIndex(n),1) = REAL( flowDirSet(ixn) )
            ENDIF
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rflow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( rFlowDiag(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rInflow' ) THEN
          found = .TRUE.
!         Calculate this from the rFlowDiag variable (outflow).
!         Loop over "active" routing points.
          DO n=1,npRoute
!           If there is a downstream gridbox, outflow from the current gridbox
!           contributes to inflow to the downstream box.
            IF ( routeNext(n) > 0 ) THEN
!             Get location of this point in the routing grid.
              CALL getXYPos( routeIndex(n),nxRoute,nyRoute,ixr,iyr )
!             The location in tmpval of the downstream point is given by routeNext,
!             since both refer to a space of size nxRoute*nyRoute points.
              j = routeNext(n)
!             Add to inflow.
              tmpval(j,1) = tmpval(j,1) + rflowDiag(ixr,iyr)
            ENDIF
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'routeIndex' ) THEN
          found = .TRUE.
!         Loop over "active" routing points.
          DO n=1,npRoute
!           The location in tmpval of this point is given by routeIndex,
!           since both refer to a space of size nxRoute*nyRoute points.
            tmpval(routeIndex(n),1) = REAL ( routeIndex(n) )
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'routeNext' ) THEN
          found = .TRUE.
!         Loop over "active" routing points.
          DO n=1,npRoute
!           The location in tmpval of this point is given by routeIndex,
!           since both refer to a space of size nxRoute*nyRoute points.
            tmpval(routeIndex(n),1) = REAL ( routeNext(n) )
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rstore' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) =  RESHAPE( routeStore(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'runoffR' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) =  RESHAPE( rRunDiag(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'bflowin' ) THEN
          found = .TRUE.
          WRITE(*,*)'ERROR: loadout: no code for var=',TRIM(var)
          STOP
!-------------------------------------------------------------------------------
        ELSEIF( var == 'flowin' ) THEN
          found = .TRUE.
          WRITE(*,*)'ERROR: loadout: no code for var=',TRIM(var)
          STOP
!-------------------------------------------------------------------------------
        ELSEIF( var == 'substore' ) THEN
          found = .TRUE.
          WRITE(*,*)'ERROR: loadout: no code for var=',TRIM(var)
          STOP
!-------------------------------------------------------------------------------
        ELSEIF( var == 'surfstore' ) THEN
          found = .TRUE.
          WRITE(*,*)'ERROR: loadout: no code for var=',TRIM(var)
          STOP
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################
!     Variables at selected points on the channel routing grid.
!     varType is 'RP'
      CASE ( 'RP' )
        np =  nxRoute*nyRoute
        nz = 1
        npOut = 1
!       Note that mapOut is used differently for 'RP' variables.
!       Each variable in profile gets one element of mapOut.
        pointsUse(1:npOut) = mapOut(iout,ivar,1)
        xxcount(:,1) = xcount(:)

!       At present, the approach is to load all points on the routing grid,
!       and then extract the one required. The loading is as for variables
!       of type 'RG', but code is just repeated, for now.

!-------------------------------------------------------------------------------
        IF( var == 'rflow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( rFlowDiag(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rInflow' ) THEN
          found = .TRUE.
!         Calculate this from the rFlowDiag variable (outflow).
!         Loop over "active" routing points.
          DO n=1,npRoute
!           If there is a downstream gridbox, outflow from the current gridbox
!           contributes to inflow to the downstream box.
            IF ( routeNext(n) > 0 ) THEN
!             Get location of this point in the routing grid.
              CALL getXYPos( routeIndex(n),nxRoute,nyRoute,ixr,iyr )
!             The location in tmpval of the downstream point is given by routeNext,
!             since both refer to a space of size nxRoute*nyRoute points.
              j = routeNext(n)
!             Add to inflow.
              tmpval(j,1) = tmpval(j,1) + rflowDiag(ixr,iyr)
            ENDIF
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rstore' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) =  RESHAPE( routeStore(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'runoffR' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) =  RESHAPE( rRunDiag(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################
!     Variables that have dim_cs1 values at every land gridpoint.
!     varType is 'SC' (soil carbon)

      CASE ('SC' )
        np = land_pts
        nz = dim_cs1
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(:,1:nz) = SPREAD( xcount(1:np),2,nz )

!-------------------------------------------------------------------------------
        IF( var == 'csPool' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = cs(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'respSPool' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = resp_s(:,:)
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################
!     Variables that have a single value at every gridpoint over land or sea.
!     varType is 'SI'

      CASE ('SI' )
        np = nx*ny
        nz = 1
        npOut = pointsOut(iout)
        pointsUse(1:npOut) = mapOut(iout,1:npOut,1)
        xxcount(:,1) = xcount(:)
!-------------------------------------------------------------------------------
        IF( var == 'conRain' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( con_rain(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'conSnow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( con_snow(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'cosz' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = cosz(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'diffFrac' ) THEN
          found = .TRUE.
          tmpval(1:np,1) = diff_frac(:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'ecan' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ecan(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'ei' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ei(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'esoil' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( esoil(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'fqw' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( fqw_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'ftl' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ftl_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'landAlbedo1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( land_albedo(:,:,1), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'landAlbedo2' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( land_albedo(:,:,2), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'landAlbedo3' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( land_albedo(:,:,3), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'landAlbedo4' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( land_albedo(:,:,4), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'latentHeat' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( latent_heat(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'latitude' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( latitude(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'longitude' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( longitude(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'lsRain' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ls_rain(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'lsSnow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ls_snow(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'LWdown' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( lw_down(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'precip' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ls_rain(:,:)+con_rain(:,:)  &
                                     +ls_snow(:,:)+con_snow(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'pstar' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( pstar(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'q1p5m' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( q1p5m(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'qw1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( qw_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'rainfall' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ls_rain(:,:)+con_rain(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snomltSurfHtf' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( snomlt_surf_htf(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowfall' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( ls_snow(:,:)+con_snow(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowMass' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( snow_mass(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'surfHtFlux' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( surf_ht_flux(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'SWdown' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( sw_down(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 't1p5m' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( t1p5m(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'taux1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( taux_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'tauy1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( tauy_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'tl1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( tl_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'tstar' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( tstar(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'u1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( u_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'u10m' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( u10m(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'v1' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( v_1(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'v10m' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( v10m(:,:), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'wind' ) THEN
          found = .TRUE.
          tmpval(1:np,1:1) = RESHAPE( SQRT( u_1*u_1 + v_1*v_1 ), (/ np,1 /) )
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################

!     Snow layer diagnostics.
!     varType is 'SN'
!     These variables are dimensioned (land_pts,ntiles,nsmax), but are treated
!     as 2-D vars of dimension (land_pts,ntiles*nsmax). The function
!     loadSnowLayers is used to load the values into tmpval such that the vertical
!     dimension (nsmax) is varying faster than the tile dimension, so that the
!     existing 2-D diagnostic code can be used to output a snow layer variable
!     as ntiles 2-D variables. It also accounts for zrevOutSnow and changes
!     values in empty layers.
!     In fact there is no need to call loadSnowLayers here while loading data,
!     something that is potentially done every timestep while accumulating for
!     a time average. A similar job could be done more efficiently by
!     calling a similar routine once, probably in writeOut but anytime before
!     data are written. Could probably also avoid some of the complications
!     for netCDF output too! But it's here for now.
!
!     NB The time average of a layer variable will be incorrect if the layer
!     disappears or is added during the interval! This is because the counter
!     assumes that a value is present at every time, so the average is found
!     simply as accum/nt. To do this properly would require a new counter that
!     counted times when a variable was defined.

      CASE ( 'SN' )
        np = land_pts
        nz = ntiles*nsmax
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(1:np,1:nz) = SPREAD( xcount(1:np),2,nz )

!-------------------------------------------------------------------------------
        IF ( var == 'rgrainL' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = loadSnowLayers( 0.0,zrevoutSnow  &
                                             ,nsnow(:,:),rgrainL(:,:,:) )
!-------------------------------------------------------------------------------
        ELSEIF ( var == 'snowDs' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = loadSnowLayers( 0.0,zrevoutSnow  &
                                             ,nsnow(:,:),ds(:,:,:) )
!-------------------------------------------------------------------------------
        ELSEIF ( var == 'snowIce' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = loadSnowLayers( 0.0,zrevoutSnow  &
                                             ,nsnow(:,:),sice(:,:,:) )
!-------------------------------------------------------------------------------
        ELSEIF ( var == 'snowLiq' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = loadSnowLayers( 0.0,zrevoutSnow  &
                                             ,nsnow(:,:),sliq(:,:,:) )
!-------------------------------------------------------------------------------
        ELSEIF ( var == 'tsnow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = loadSnowLayers( tm,zrevoutSnow  &
                                             ,nsnow(:,:),tsnow(:,:,:) )
        ENDIF
!###############################################################################
!###############################################################################

!     Soil layer diagnostics.
!     varType is 'SO'

      CASE ( 'SO' )
        np = land_pts
        nz = nsoil
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(1:np,1:nz) = SPREAD( xcount(1:np),2,nz )

!-------------------------------------------------------------------------------
        IF( var == 'bSoil' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = b(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'ext' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = ext(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'hCapSoil' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = hcap(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'hConSoil' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = hcon(:,1:nz)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'satCon' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = satcon(:,1:nz)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'sathh' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = sathh(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'smcl' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = smcl(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'soilWet' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = sthf(:,:) + sthu(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'sthf' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = sthf(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'sthu' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = sthu(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'tSoil' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = t_soil(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'vsmcCrit' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = smvccl(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'vsmcSat' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = smvcst(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'vsmcWilt' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = smvcwt(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ELSEIF( var == 'wFlux' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = wfluxDiag(:,:)
          IF ( zrevOutSoil ) tmpval(1:np,1:nz)=reverseCols( tmpval(1:np,1:nz) )
!-------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################

!     Tile diagnostics (on land grid).
!     varType is 'TI'

      CASE ( 'TI' )
        np = land_pts
        nz = ntiles
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(1:np,1:nz) = xcountTile(1:np,1:nz)

!-------------------------------------------------------------------------------
        IF( var == 'alb1T' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = alb_tile(:,:,1)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'alb2T' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = alb_tile(:,:,2)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'alb3T' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = alb_tile(:,:,3)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'alb4T' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = alb_tile(:,:,4)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'canopyT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = canopy(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'catchT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = catch(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'ecanT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = ecan_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'eiT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = ei_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'esoilT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = esoil_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'fqwT' ) THEN
          found = .TRUE.
!         Note that fqw_tile does not (always) do this job!
          DO iz=1,nz
            DO j=1,tile_pts(iz)
              l = tile_index(j,iz)
              tmpval(l,iz) = ecan_tile(l,iz) + ei_tile(l,iz) + esoil_tile(l,iz)
!             Add lake evaporation.
              IF ( iz == lake ) tmpval(l,iz) = tmpval(l,iz) + fqw_tile(l,iz)
            ENDDO
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'ftlT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = ftl_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'gcT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = gc(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'leT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = le_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'nsnow' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = REAL( nsnow(:,:) )
!-------------------------------------------------------------------------------
       ELSEIF( var == 'q1p5mT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = q1p5m_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'radnetT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = radnet_tile(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'emisT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) =  emis_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'rgrainT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = rgrain(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowCanMeltT' ) THEN
          found = .TRUE.
!         Only include tiles where canopy snow model is used.
          DO i=1,ntiles
            IF ( canSnowTile(i) ) tmpval(1:np,i) = melt_tile(:,i)
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowCanT' ) THEN
          found = .TRUE.
!         Only include tiles where canopy snow model is used.
          DO i=1,ntiles
            IF ( canSnowTile(i) ) tmpval(1:np,i) = snow_tile(:,i)
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowDepthT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = snowDepth(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowGrCanMeltT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = snowGMeltDiag(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowGrCanT' ) THEN
          found = .TRUE.
!         Only include tiles where canopy snow model is used.
          DO i=1,ntiles
            IF ( canSnowTile(i) ) tmpval(1:np,i) = snow_grnd(:,i)
          ENDDO
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowGroundRhoT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = rho_snow_grnd(:,:)
!-------------------------------------------------------------------------------
        ELSEIF( var == 'snowGroundT' ) THEN
          found = .TRUE.
          DO i=1,ntiles
            IF ( canSnowTile(i) ) THEN
              tmpval(1:np,i) = snow_grnd(:,i)
            ELSE
              tmpval(1:np,i) = snow_tile(:,i)
            ENDIF
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowIceT' ) THEN
          found = .TRUE.
          DO n=1,ntiles
            DO j=1,tile_pts(n)
              i = tile_index(j,n)
              tmpval(i,n) = SUM( sice(i,n,1:nsnow(i,n)) )
            ENDDO
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowLiqT' ) THEN
          found = .TRUE.
          DO n=1,ntiles
            DO j=1,tile_pts(n)
              i = tile_index(j,n)
              tmpval(i,n) = SUM( sliq(i,n,1:nsnow(i,n)) )
            ENDDO
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowMassT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = snow_tile(:,:)
!         Add snow below canopy.
          DO n=1,ntiles
            IF ( canSnowTile(n) ) tmpval(1:np,n) =  &
             tmpval(1:np,n) + snow_grnd(:,n)
          ENDDO
!-------------------------------------------------------------------------------
       ELSEIF( var == 'snowMeltT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = melt_tile(:,:)
!         Add melting of snow below canopy.
!          IF ( can_model == 4 ) tmpval(1:np,1:nz) =  &
!             tmpval(1:np,1:nz) + snowGMeltDiag(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 't1p5mT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = t1p5m_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'tstarT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = tstar_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'z0T' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = z0_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'surfHtStoreT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = surf_ht_store(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'surfHtFluxT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = surf_htf_tile(:,:)
!-------------------------------------------------------------------------------
       ELSEIF( var == 'anthropHtFluxT' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = anthrop_heat(:,:)
!-------------------------------------------------------------------------------
       ENDIF
!###############################################################################
!###############################################################################

!     Surface type diagnostics (on land grid).
!     varType is 'TY'

      CASE ( 'TY' )
        np = land_pts
        nz = ntype
        npOut = pointsOutLand(iout)
        pointsUse(1:npOut) = mapOutLand(iout,1:npOut,1)
        xxcount(1:np,1:nz) = xcountTile(1:np,1:nz)  !  tmp xx ONLY for ntiles=ntype, I think

!--------------------------------------------------------------------------------
        IF( var == 'frac' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = frac(:,:)
!         Prevent frac=0 from being set to undef.
          xxcount(1:np,1:nz) = MAXVAL(countTile(1:np,1:nz))  !  tmp xx ONLY for ntiles=ntype, I think
!--------------------------------------------------------------------------------
        ELSEIF( var == 'tileIndex' ) THEN
          found = .TRUE.
          tmpval(1:np,1:nz) = REAL( tile_index(:,:) )
!--------------------------------------------------------------------------------
        ENDIF
!###############################################################################
!###############################################################################
      CASE default
        WRITE(*,*)'ERROR: loadout: varType=',varType(jvar)
        WRITE(*,*)'No code for this varType.'
        WRITE(*,*)'Stopping in loadout'
        STOP
    END SELECT
!###############################################################################

    IF ( .NOT. found ) THEN
      WRITE(*,*)'ERROR: loadout: Did not find code for varName=',TRIM(var)
      WRITE(*,*)'Stopping in loadout'
      STOP
!     Note that if code exists for this variable, but under wrong varType, it will not be found.
    ENDIF
!-------------------------------------------------------------------------------

!   If the time average is to be set to missing, and this is a time average variable, reset counter.
!   Note that this does not affect non-time-average variables.
!   (Note: if we later want to extend this so that, e.g. accumulations can be set to missing, might
!    be better to apply this in subroutine writeOut, on writing.)
    IF ( undefTave .AND. tmeanVar(varPos(iout,ivar)) ) xxcount(1:np,1:nz)=0

!   Now extract the desired points, do time averaging etc.
    DO iz=1,nz
      ipos = iposOne + (iz-1)*npOut      !   location of first value in outval
      iposWrite = iposWrite + npOut      !   location of last value in outvalWrite
      CALL loadOutVarReal( iout,np,npOut,xxcount(1:np,iz),pointsUse(1:npOut)  &
                          ,tmpval(1:np,iz),doAccum,doTmean,taccumVar(jvar),tmeanVar(jvar)  &
                          ,outval(ipos:ipos+npOut-1),outvalWrite(iposWrite-npOut+1:iposWrite) )
    ENDDO

  ENDDO   !  ivar (variables)

  END SUBROUTINE loadOut
!######################################################################################
!######################################################################################
!######################################################################################

! subroutine loadOutVarReal
!
! Depending upon input, does one of the following:
! 1) if a variable is being accumulated prior to getting average (or is simply to be accumulated),
!    adds the latest value to the accumulation (in vector outval).
!    At the end of the averaging period, the average is calculated.
!    The final average or accumulation is held in the output vector (outvalWrite)
! 2) if a snapshot variable is to be output this time, loads data into
!    the output vector (outvalWrite).
! Note that only the requested points are processed here - no mapping onto output grid is done.


  SUBROUTINE loadOutVarReal( iout,npIn,npOut,count,mapOut,inval,doAccum,doTmean,taccumVar,tmeanVar  &
                            ,outval,outvalWrite )

  USE inout, ONLY : ntOutPer,outSamPer,undefOut
  USE grid_utils, ONLY : mapAtoB
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout,npIn,npOut
!   ,ntExpect     !  the expected number of times in a time accumulation (for average)

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    COUNT(npIn)    &!  counter of # of valid times at each point
!   ,countSnap(npIn) &!  flag indicating points at which input variable is defined at current time
!    in future may want to use this, so that only add to accumulation where variable is currently defined
!    This would allow for tiles to change within an output period (because frac goes from >0 to 0/frac_min.
!    See snow code in mymoses - that probably does the same. So maybe don't need this var??
   ,mapOut(npOut)   !  mapping: a list of the points in inval that are to be output

  INTEGER ::  &!  local SCALARS
    i  !  loop counter

  INTEGER :: &!  local ARRAYS
    mapOut2(npOut) &!  mapping: a list of the points in outval that are to be filled
   ,xcount(npOut)   !  count mapped onto output grid

  REAL, INTENT(in) ::  &!  in ARRAYS
     inval(npIn)   !   input data

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    doAccum,doTmean,taccumVar,tmeanVar

  REAL, INTENT(inout) ::  &!  inout ARRAYS
    outval(npOut)        &!  used to hold time accumulations (for average or accumulation)
   ,outvalWrite(npOut)    !  values as later written to file
!        The separation of outval and outvalWrite is sometimes unnecessary (e.g. outval alone
!        would often suffice for snapshot variables), but allows e.g.:
!         * time accumulations to be multipied by outSamPer before output, without affecting any
!           continuing accumulation. This alone would probably not merit the extra code...and we
!           could multiply output in writeOut.
!         * the units of variables could be changed before output (not currently coded!).
!            Again could be done anyway in writeOut.
!        So...the case for this extra faff is currently not particularly strong....but we'll do
!        it this way until some more complicated cases (e.g. variables that are only defined
!        when there's snow) are coded.
!
!--------------------------------------------------------------------------------

! Indicate that the 1st npOut positions of output are to be filled.
  mapOut2(:) = (/ (i,i=1,npOut) /)

! Map the counter onto the same grid as the output variable.
  xcount(:) = mapAtoB( npIn,npOut,npOut,mapOut(1:npOut),mapOut2(1:npOut),0,COUNT(1:npIn) )

!-------------------------------------------------------------------------------
  IF ( taccumVar ) THEN
!   Time accumulation variable.

!   Add to the accumulation.
    IF ( doAccum ) outval(1:npOut) = outval(1:npOut) + &
                     mapAtoB( npIn,npOut,npOut,mapOut(1:npOut),mapOut2(1:npOut),0.0,inval(1:npIn) )
!   Get final value.
    IF ( doTmean ) THEN
!     Get final accumulation and store this is the output vector used for writing (outvalWrite).
!     The final value is multiplied by the number of timesteps in sampling period - note that this does not
!     affect any continuing accumulation, since that is held in outval.
!     Note that, as for time averages,  at present no account is taken of the possibility that first or
!     last increment may cover a smaller time interval (e.g. output started midway through averaging period).
      WHERE ( xcount(:) > 0 )
        outvalWrite(:) = outval(:) * REAL(outSamPer(iout))
      ELSEWHERE
        outvalWrite(:) = undefOut
      endwhere
    ENDIF

!-------------------------------------------------------------------------------
  ELSEIF ( tmeanVar )  THEN
!   Time mean variable.

!   Minimise the accumulation by weighting each value by 1/nt.
    IF ( doAccum ) outval(1:npOut) = outval(1:npOut) + &
                     mapAtoB( npIn,npOut,npOut,mapOut(1:npOut),mapOut2(1:npOut),0.0,inval(1:npIn) ) &
                        / ntOutPer(iout)
    IF ( doTmean ) THEN
!     Get final average. If the (estimated!) expected number of times have not been added to the
!     accumulation, correct for this. The final values is stored in the output vector used for
!     writing (outvalWrite).
!     Note that, as for time accumulations, at present no account is taken of the possibility that first or
!     last increment may cover a smaller time interval (e.g. output started midway through averaging period).
      outvalWrite(:) = outval(:)
      WHERE ( xcount(:) > 0 )
!       Only normalise where needed. This gives agreement with earlier release of code.
        WHERE ( xcount(:) /= ntOutPer(iout) ) outvalWrite(:) = outvalWrite(:) * REAL(ntOutPer(iout))/xcount(:)
      ELSEWHERE
        outvalWrite(:) = undefOut
      endwhere
    ENDIF

!-------------------------------------------------------------------------------
  ELSE  !  a snapshot variable (not accum or tmean)

!   Save a snapshot value. This is saved in the output vector (outvalWrite) since a snapshot is only saved
!   when it is due for output.
    outvalWrite(:) = mapAtoB( npIn,npOut,npOut,mapOut(1:npOut),mapOut2(1:npOut),undefOut,inval(1:npIn) )

!   Where count=0 there are no defined values, so make undef.
    WHERE ( xcount(:)==0 ) outvalWrite(:) = undefOut

  ENDIF    ! mean/snapshot
!-------------------------------------------------------------------------------

  END SUBROUTINE loadOutVarReal

!################################################################################
!################################################################################
!################################################################################

! subroutine writeOut
! Manages Writing of data to file.
! May apply a mapping to write output to selected points in a larger grid.

  SUBROUTINE writeOut( a_step,iout,outWriteCount,npOut,nvalMax,outputDate  &
                      ,outputTime,outvalWrite )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    dim_cs1,nsmax,nsoil=>sm_levels,ntiles

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatNc  &
!  imported scalars with intent(in)
    ,gradsNC,outFormat,undefOut  &
!  imported arrays with intent(in)
    ,mapOut,mapOutLand,nvarOut,outCompress,outGridNxy,outNpWrite  &
    ,outPer,outTimeID,outUnit,outVarID,pointsOut  &
    ,pointsOutLand,varNlev,varPos,varType  &
!  imported arrays with intent(inout)
    ,irecPrevOut

  USE grid_utils, ONLY :  &
!  imported procedures
    mapAtoB

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     npft,ntype

  USE readWrite_Mod, ONLY :  &
!  imported procedures
    writeVar

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
    ispin,spinUp

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timestep
!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: ndimMax = 5  !  max possible number of netCDF dimensions
!                                       for a single variable

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) ::  &
    a_step        &!  timestep number
   ,iout          &!  counter of output profile
   ,outWriteCount &!  timestep (time level) in file
   ,npOut         &!  # of points in output grid (written to file). This is >= np.
   ,nvalMax       &!  maximum number of values in any variable (on output grid)
   ,outputDate    &!  date (yyyymmdd) of output data
   ,outputTime     !  time of day (s) of output data

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  REAL, INTENT(in) ::  &
    outvalWrite(:)   !  data values to be output

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    i,ilev,ipos,itile,ivar  &!  work/loop counters
   ,jvar              &!  work
   ,ndim              &!  number of dimensions for netCDF
   ,np                &!  # of points to be output
   ,tpos              &!  work
   ,varID             &!  a netCDF ID
   ,zdim               !  index (1:ndim) of z dimension for certain variables

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER :: ncCount(ndimMax)   !  count argument for netCDF variables
  INTEGER :: ncStart(ndimMax)   !  start argument for netCDF variables
  INTEGER :: tmpMapOut(npOut,2) !  mapping between input points and output locations
!                  !    Only elements 1:np are used.
!                  !    1st dimension - input points.
!                  !    2nd dimension - output locations.

  REAL :: tmpVal(nvalMax)  !  output data (elements of outvalWrite
!                              mapped onto output grid)

!-----------------------------------------------------------------------
! If using ASCII output, write a timestep message.
!-----------------------------------------------------------------------
  IF ( outFormat == formatAsc ) THEN
    IF ( .NOT. spinUp ) THEN
      WRITE(outunit(iout),"(a,i7,a,i8,i7)") 'timestep:',a_step  &
             ,' output date and time(s): ',outputDate,outputTime
    ELSE
      WRITE(outunit(iout),"(a,i7,a,i8,i7)") 'timestep:',a_step  &
            ,' output date and time(s): '  &
            ,outputDate,outputTime,' spin up cycle: ',ispin
    ENDIF
  ENDIF

!-----------------------------------------------------------------------
! Write time variables for netCDF output.
!-----------------------------------------------------------------------
  IF ( outFormat == formatNc ) THEN

    ndim = 1
    ncCount(1) = 1
    ncStart(1) = outWriteCount

!   Write time since start of file.
    tmpval(1) = REAL(outWriteCount-1) * outPer(iout) * timeStep
    varID = outTimeID(iout,1)
    CALL writeVar( irecPrevOut(iout),varID,outUnit(iout)  &
                    ,ncCount(1:ndim),ncStart(1:ndim)  &
                    ,tmpval(1:1),outFormat,'writeOut: time' )

!   Write timesteps since start of file.
!   As usual, I reckon one should use outWriteCount-1, but the convention
!   seems to be otherwise.
    tmpval(1) = REAL(outWriteCount)
    varID = outTimeID(iout,2)
    CALL writeVar( irecPrevOut(iout),varID,outUnit(iout)  &
                    ,ncCount(1:ndim),ncStart(1:ndim)  &
                    ,tmpval(1:1),outFormat,'writeOut: timesteps' )

  ENDIF  !  outFormat

!-----------------------------------------------------------------------
! Loop over variables in this output profile (stream).
!-----------------------------------------------------------------------
  ipos = 0   !  initialise

  DO ivar=1,nvarOut(iout)
    jvar = varPos(iout,ivar)  !  position in overall list of selected variables
!-----------------------------------------------------------------------
!   Get a list of output points.
!-----------------------------------------------------------------------
    SELECT CASE ( varType(jvar) )
      CASE ( 'RG', 'SI' )
        np = pointsOut(iout)
        tmpMapOut(1:np,2) = mapOut(iout,1:np,2)
      CASE ( 'LA', 'PF', 'SC', 'SN', 'SO', 'TI', 'TY' )
        np = pointsOutLand(iout)
        tmpMapOut(1:np,2) = mapOutLand(iout,1:np,2)
      CASE ( 'RP' )
        np = 1
        tmpMapOut(1:np,2) = mapOut(iout,1:np,2)
      CASE default
        WRITE(*,*)'ERROR: writeOut: varType=',varType(jvar)
        WRITE(*,*)'writeOut: No code for this varType.'
        STOP
    END SELECT

!-----------------------------------------------------------------------
!   Indicate that the 1st np positions of each field in outvalWrite are to be used.
!-----------------------------------------------------------------------
    tmpMapOut(1:np,1) = (/ (i,i=1,np) /)

!-----------------------------------------------------------------------
!   For netCDF output, get number of dimensions and set indices.
!   At present there is no possibility to output a subset of vertical
!   levels, and all variables have similar dimensions, namely tyx to tzyz.
!   If a variable has more than one level we use all dims, else omit z.
!-----------------------------------------------------------------------
    IF ( outFormat == formatNc ) THEN

!-----------------------------------------------------------------------
!     Set default starts and counts for all dimensions.
!-----------------------------------------------------------------------
      ncStart(:) = 1
      ncCount(:) = 1

!-----------------------------------------------------------------------
!     Get number of dimensions, assuming not compressed output, and set
!     count for layers.
!-----------------------------------------------------------------------
      SELECT CASE  ( varType(jvar) )
        CASE ( 'LA', 'RG', 'SI' )
          ndim = 3   !  x,y,t
        CASE ( 'PF' )
          ndim = 4   ! x,y,z,t
          ncCount(3) = npft
        CASE ( 'SC' )
          ndim = 4   ! x,y,z,t
          ncCount(3) = dim_cs1
        CASE ( 'SN' )
          IF ( gradsNc ) THEN
            ndim = 4   ! x,y,z,t
            zdim = 3   ! save z dimension
            ncCount(3) = 1  !  we write each layer separately
          ELSE
            ndim = 5   ! x,y,z,tile,t
            ncCount(3) = nsmax
            ncCount(4) = ntiles
          ENDIF
        CASE ( 'SO' )
          ndim = 4   ! x,y,z,t
          ncCount(3) = nsoil
        CASE ( 'TI' )
          ndim = 4   ! x,y,z,t
          ncCount(3) = ntiles
        CASE ( 'TY' )
          ndim = 4   ! x,y,z,t
          ncCount(3) = ntype
        CASE default
          WRITE(*,*)'ERROR: writeOut: ndim: varType=',varType(jvar)
          WRITE(*,*)'writeOut: No code for this varType.'
          STOP
      END SELECT

!-----------------------------------------------------------------------
!     Set counts (sizes) of field for "horizontal" dimensions.
!     Reduce number of dimensions if compressing horizontal dims.
!-----------------------------------------------------------------------
      IF ( outCompress(iout) ) THEN
        ndim = ndim - 1    !   remove a dimension (y)
!       Shift all non-horizontal dimensions.
        ncCount(2:ndim) = ncCount(3:ndim+1)
        ncCount(1) = outNpWrite(iout)          !  # of points in vector
        zdim = zdim - 1
      ELSE
        ncCount(1:2) = outGridNxy(iout,1:2)    !  # of x and y points
      ENDIF

!-----------------------------------------------------------------------
!     Set time level. Time is always last dimension.
!-----------------------------------------------------------------------
      ncCount(ndim) = 1
      ncStart(ndim) = outWriteCount

    ELSE

!-----------------------------------------------------------------------
!     Set defaults for non-netCDF output.
!-----------------------------------------------------------------------
      ndim = 1
      zdim = 1

    ENDIF  !  outFormat

!-----------------------------------------------------------------------
!   Call writing routine.
!   For files that are not netCDF, this is done separately for each
!   level. For netCDF, a whole variable can be output, unless it is a
!   snow layer variables for GrADS-readable output.
!-----------------------------------------------------------------------
    IF ( outFormat==formatNc .AND.  &
         .NOT. (varType(jvar)=='SN' .AND. gradsNc) ) THEN

      tpos = 0
      DO ilev=1,varNlev(jvar)
        tpos = tpos + npOut   !  last location filled in tmpval
        ipos = ipos + np   !  last location used from outValWrite
!       Map data.
        tmpval(tpos-npOut+1:tpos) = mapAtoB( np,np,npOut  &
                                  ,tmpMapOut(1:np,1),tmpMapOut(1:np,2)  &
                                  ,undefOut,outvalWrite(ipos-np+1:ipos) )
      ENDDO  !  ilev

!     Write data.
      CALL writeVar( irecPrevOut(iout),outVarID(jvar),outUnit(iout)  &
                    ,ncCount(1:ndim),ncStart(1:ndim)  &
                    ,tmpval(1:tpos),outFormat,'writeOut: var' )

    ELSE

!-----------------------------------------------------------------------
!     Formats other than netCDF, or a snow layer variable for GrADS+netCDF.
!-----------------------------------------------------------------------
      DO ilev=1,varNlev(jvar)

        ipos = ipos + np   !  last location used from outValWrite
        irecPrevOut(iout) = irecPrevOut(iout) + 1

        varID = outVarID(jvar)

        IF ( outFormat == formatNc ) THEN
!         GrADS-readable netCDF files, snow variables. ntiles variables were
!         created, but only the first given a netCDF ID. Here we assume that
!         the rest have the subsequent netCDF IDs! This is a gross hack!
          itile = CEILING( REAL(ilev) / REAL(nsmax) )
          varID = varID + itile - 1

!         Get start index for z. Note we only have to deal with SN variables.
          IF ( varNlev(jvar) > 1 ) THEN
            ncStart(zdim) = MOD( ilev, nsmax )
            IF ( ncStart(zdim) == 0 ) ncStart(zdim) = nsmax
          ENDIF  !  varNlev
        ENDIF  !  outFormat

!       Map data.
        tmpval(1:npOut) = mapAtoB( np,np,npOut,tmpMapOut(1:np,1),tmpMapOut(1:np,2)  &
                                  ,undefOut,outvalWrite(ipos-np+1:ipos) )

!       Write data.
        CALL writeVar( irecPrevOut(iout),varID,outUnit(iout)  &
                      ,ncCount(1:ndim),ncStart(1:ndim)  &
                      ,tmpval(1:npOut),outFormat,'writeOut: var' )
      ENDDO  !  ilev
    ENDIF  !  outFormat

  ENDDO  !  ivar

  END SUBROUTINE writeOut
!################################################################################
!################################################################################
!################################################################################
! subroutine rewrite_ctl
! Check if a GrADS ctl file had correct number of times in tdef. If not, rewrite.

  SUBROUTINE rewrite_ctl( checkTemplate,iout )

  USE file_utils, ONLY : closeFile,fileUnit,openFile
  USE inout, ONLY : echo,formatAsc,ntCtl,ntCtlNeed,outCtlFile  &
                   ,outDataFile,outWarnCtl,outTemplate
  USE misc_utils, ONLY : getFields

  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  SCALARS
    nlineMax = 200

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout     !  counter of output profile

  INTEGER ::  &!  local SCALARS
    i       &!  loop counter
   ,ierr    &!  error flag
   ,nfield  &!  number of fields (words) found in tdef line
   ,nline   &!  number of lines read from ctl file
   ,unit     !  unit for file

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    checkTemplate   !  T if template ctl file are to be checked on this call
!                   !    This is T at end of either spin up (if separate ctl file) and main run.

  LOGICAL ::  &!  local SCALARS
    abandon   &!  T if try to rewrite a file but encounter an error
   ,errFlag   &!  error flag
   ,found      !  work

  CHARACTER(len=15) ::  &!  local ARRAYS
    field(5)  !  the fields that make up the tdef line
!             !    Length should be enough for longest field - likely the time/date, but
!             !    depends on format used for write in writeGradsCtl.

  CHARACTER(len=250) ::  &!  local arrays
    inLine(nlineMax)    !  lines as read from ctl file (must be long enough!)

!--------------------------------------------------------------------------------

! If the name of the output data file is empty, no output data has yet been written for this profile,
! so nothing to do. (i.e. this call is associated with the first write.)
  IF ( LEN_TRIM(outDataFile(iout)) == 0 ) RETURN

! Initialise.
  abandon = .FALSE.
!-------------------------------------------------------------------------------

  IF ( (checkTemplate .AND. outTemplate(iout)) .OR. .NOT. outTemplate(iout) ) THEN

    IF ( ntCtlNeed(iout) /= ntCtl(iout) ) THEN

!     Number of times as written in tdef is not correct.
      IF ( echo ) THEN
        WRITE(*,"(/,2a)")'rewrite_ctl: altering ',TRIM(outCtlFile(iout))
        WRITE(*,*)'Changing nt from',ntCtl(iout),' to ',ntCtlNeed(iout)
      ENDIF

!     Establish if the previous ctl file still exists.
      IF ( LEN_TRIM(outCtlFile(iout)) == 0 ) THEN
        WRITE(*,*)'ERROR: rewrite_ctl: name of previous ctl file name is empty!'
        WRITE(*,*)'This should not happen - but will continue run anyway!'
        outWarnCtl = .TRUE.
        RETURN
      ENDIF

      INQUIRE( file=outCtlFile(iout), exist=found )

      IF ( found ) THEN

!       Open existing ctl file.
        unit = fileUnit( formatAsc )
        CALL openFile(1,.FALSE.,unit,'read',formatAsc,outCtlFile(iout),'old')
!       Read file.
        nline = 1
        DO
          READ(unit,"(a)",iostat=ierr) inLine(nline)
          IF ( ierr /= 0 ) EXIT
          nline = nline + 1
          IF ( nline > nlineMax ) THEN
            WRITE(*,*)'WARNING: Have not reached end of ctl file! Increase array size!'
!           Could attempt to deal with this, but not bothering.
            abandon = .TRUE.
            EXIT
          ENDIF
        ENDDO

        IF ( ierr/=0 .AND. .NOT.abandon ) THEN
!         Assume that the error is that end of file was reached.
!         Now check that a tdef line can be found.
          found = .FALSE.
          DO i=1,nline
            SELECT CASE ( inLine(i)(1:4) )
              CASE ( 'tdef', 'TDEF' )
                found = .TRUE.
!               Split line into 5 space-delimited words.
!                CALL getFields( inLine(i),' ',.TRUE.,field,nfield )
                CALL getFields( LEN(field),.TRUE.,.TRUE.,.FALSE.,' ',inLine(i),errFlag  &
                               ,nfMaxUse=.TRUE.,nfield=nfield,fieldList=field )
                IF ( nfield /= 5 ) THEN
                  WRITE(*,*)'ERROR: did not find 5 fields in tdef line.'
                  abandon = .TRUE.
                ENDIF
                EXIT  ! leave do loop
            END SELECT
          ENDDO  !  i (lines)
          IF ( .NOT. found ) THEN
            WRITE(*,*)'ERROR: did not find tdef line.'
            abandon = .TRUE.
          ELSE
!           Edit the tdef line.
            WRITE(inLine(i),"('tdef ',i7,' linear ',a,tr1,a)") ntCtlNeed(iout),TRIM(field(4)),TRIM(field(5))
          ENDIF
          IF ( .NOT. abandon ) THEN
!           Rewrite ctl file.
            nline = nline - 1
!           Close file (2nd arg=TRUE), then reopen for write.
            CALL openFile(1,.TRUE.,unit,'write',formatAsc,outCtlFile(iout),'replace')
            DO i=1,nline
              WRITE(unit,"(a)") TRIM(inLine(i))
            ENDDO
          ENDIF  !  NOT abandon
        ENDIF

!       Close file to release unit.
        CALL closeFile( unit,formatAsc )

      ELSE  !  NOT found

        WRITE(*,*)'WARNING: rewrite_ctl: previous ctl file cannot be found.'
        WRITE(*,*)'File may have been moved or deleted.'
        abandon = .TRUE.

      ENDIF  !  found
    ENDIF  !  ntCtl

    IF ( abandon ) THEN
      WRITE(*,*)'Error raised - abandoning attempt to rewrite ctl file.'
      WRITE(*,*)'file: ',TRIM(outCtlFile(iout))
      WRITE(*,*)'Number of times in tdef was ',ntCtl(iout),', should be ',ntCtlNeed(iout)
      outWarnCtl = .TRUE.
    ENDIF

  ENDIF

  END SUBROUTINE rewrite_ctl
!################################################################################
!################################################################################
!################################################################################
! subroutine reset_outval
! Reset output variables, including the output storage variable, after output
! has been written.
! Snapshot and time-averaged variables are set to zero. Time-accumulations are
! only reset at the end of a "section" of the run.

  SUBROUTINE reset_outval( iout )

  USE inout, ONLY :   &
!  imported arrays with intent(in)
    nvarOut,pointsOut,pointsOutLand,taccumVar,varNlev,varPos,varStartPos,varType  &
!  imported arrays with intent(inout)
   ,outval  &
!  imported arrays with intent(out)
   ,outStep,outStepSamPer

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
    endSec
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout     !  counter of output profile

  INTEGER ::  &!  local SCALARS
    iposOne,ivar,jvar,np    !  work

!-------------------------------------------------------------------------------
! Reset counters.
  outStep(iout) = 0
  outStepSamPer(iout) = 0

! Reset parts of the storage space (outval).

  DO ivar=1,nvarOut(iout)
    jvar = varPos(iout,ivar)  !  position in output
    iposOne = varStartPos(jvar)  !  position of first value (of first level) in output
!   Get number of x,y points for this variable.
!
    SELECT CASE ( varType(jvar) )
      CASE ( 'LA','PF','SC','SN','SO','TI','TY' )
        np = pointsOutLand(iout)
      CASE ( 'RP' )
        np = 1
      CASE ( 'RG', 'SI' )
        np = pointsOut(iout)
      CASE default
        WRITE(*,*)'ERROR: reset_outval: varType=',varType(jvar)
        WRITE(*,*)'No code for this varType.'
        WRITE(*,*)'Stopping in reset_outval'
        STOP
    END SELECT

    IF ( taccumVar(jvar) ) THEN
!     Only reset at end of a section of run.
      IF ( endSec ) outval(iposOne:iposOne+varNlev(jvar)*np-1) = 0.0
    ELSE
!     Snapshot or time-average variable.
      outval(iposOne:iposOne+varNlev(jvar)*np-1) = 0.0
    ENDIF
  ENDDO

  END SUBROUTINE reset_outval
!################################################################################
!################################################################################
!################################################################################
!
! function gbmTileDiag
! Internal procedure in module output_mod.
! Get gridbox mean of a tile variable.


  FUNCTION gbmTileDiag( inval,tile_frac,tileMask ) RESULT(outval)

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts,ntiles  &
!  imported arrays with intent(in)
   ,tile_index,tile_pts

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL ::   &!  function result
    outval(land_pts)  !  function result

  REAL, INTENT(in) ::  &!  in arrays
    inval(land_pts,ntiles)     &!  input tile field.
   ,tile_frac(land_pts,ntiles)  !  tile fractions

  LOGICAL, INTENT(in), OPTIONAL :: tileMask(ntiles)  !  mask indicating which
!                         tiles are to be included in average. If not given,
!                         all tiles are used.
!                         Introduced for canopy snow diagnostics.

  INTEGER ::  &!  local SCALARS
    i,p,t  !  work/loop counters

  LOGICAL :: tmask(ntiles) !  mask indicating which tiles are to be included in average.

!-------------------------------------------------------------------------------
! Set mask.
  IF ( PRESENT( tileMask ) ) THEN
    tMask(:) = tileMask(:)
  ELSE
    tMask(:) = .TRUE.
  ENDIF

! Initialise the average.
  outval(:) = 0.0

  DO t=1,ntiles
    IF ( tMask(t) ) THEN
      DO i=1,tile_pts(t)
        p = tile_index(i,t)
        outval(p) = outval(p) + tile_frac(p,t)*inval(p,t)
      ENDDO
    ENDIF
  ENDDO

  END FUNCTION gbmTileDiag
!################################################################################
!################################################################################
!################################################################################
!
! function get_fsnow
! Internal procedure in module output_mod.
! Calculate a snow-covered fraction consisent with albedo formulation.
! This is useful only as long as this code agrees with the form of that used
! to calculate tile albedo (e.g. in subroutine tile_albedo)!


  FUNCTION get_fsnow( tile_frac ) RESULT(outval)

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts,ntiles  &
!  imported arrays with intent(in)
    ,tile_index,tile_pts

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     ntype

  USE p_s_parms, ONLY :  &
!  imported arrays with intent(in)
     z0_tile

  USE prognostics, ONLY :  &
!  imported arrays with intent(in)
     snowdepth

  USE snow_param, ONLY :  &
!  imported arrays with intent(in)
     canSnowTile

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     can_model,l_aggregate,l_point_data

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL ::   &!  function result
    outval(land_pts)  !  function result

! Array arguments with intent(in)
  REAL, INTENT(in) :: tile_frac(land_pts,ntiles)  !  tile fractions

! Local scalar variables.
  INTEGER :: j,l,n

! Local array variables.
  REAL :: fsnow(land_pts,ntype)  !  weight given to snow albedo
  REAL :: snowD(land_pts)        !  snow depth (m)
  REAL :: z0(land_pts)           !  roughness length (m)

!-------------------------------------------------------------------------------

! Calculate the snow-covered fraction for each type.
! This is the weight given to snow albedo when calculating tile albedo - must
! be consistent with form used in subroutine TILE_ALBEDO.

  fsnow(:,:) = 0.0
  DO N=1,NTYPE

    IF (.NOT. L_AGGREGATE) THEN
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        SNOWD(L) = SNOWDEPTH(L,N)
        Z0(L) = Z0_TILE(L,N)
      ENDDO
    ENDIF

    IF ( l_point_data ) THEN
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        IF ( SNOWD(L) .GT. 0.) FSNOW(L,N) = 1.0 - EXP( -50.0*SNOWD(L) )
      ENDDO
    ELSE
      DO J=1,TILE_PTS(N)
        L = TILE_INDEX(J,N)
        IF ( SNOWD(L) .GT. 0.) FSNOW(L,N) = SNOWD(L) / ( SNOWD(L) + 10.*Z0(L) )
      ENDDO
    ENDIF

  ENDDO   !  ntype

! Calculate gridbox mean.
  outval(:) = gbmTileDiag( fsnow,tile_frac )

  END FUNCTION get_fsnow
!###############################################################################
!###############################################################################
!###############################################################################
! function soilMoistDiag
! Internal procedure in module output_mod.
! Calculate a soil moisture diagnostic.

  FUNCTION soilMoistDiag( depthArg,diagType,swet ) RESULT(outval)

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts,sm_levels,soil_pts  &
!  imported arrays with intent(in)
   ,soil_index

  USE c_densty, ONLY :  &
!  imported scalar parameters
    rho_water

  USE soil_param, ONLY : dzsoil

  USE p_s_parms, ONLY :  &
!  imported arrays with intent(in)
    vsmcSat => smvcst,vsmcWilt => smvcwt

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
   i,iz,j        !  work/loop counters

  REAL ::   &!  function result
    outval(land_pts)  !  function result
!             This is a soil moisture (kg m-2).

  REAL, INTENT(in) ::  &!  in arrays
    swet(land_pts,sm_levels)   !  input soil layer field.
!                   This is a soil wetness (i.e. soil moisture as fraction of saturation)

  REAL, INTENT(in) ::   &! in scalars
    depthArg     !  depth over which the diagnostic to be calculated (m)
!                      If < 0, this means use the full soil depth.

  REAL ::  &!  local scalars
   dz        &!  a layer thickness (m)
  ,zbot      &!  depth to the bottom of a soil layer (m)
  ,zDiag     &!  depth for diagnostic (m)
  ,ztop       !  depth to the top of a soil layer (m)

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    diagType       !  the the type of diagnostic required
!                         'total': the total soil moisture
!                         'avail': the soil moisture available for evaporation (i.e. above wilting point)

  REAL ::  &!  local arrays
   vsmcOff(land_pts,sm_levels)   !  volumetric soil moisture content used as offset (dimensionless)

!-----------------------------------------------------------------------

! Initialise. This value is kept at non-soil points (land ice points).
  outval(:) = 0.0

! Get depth over which the diagnostic is calculated.
! If the input depth is <0, use the full depth of soil.
  zDiag = depthArg
  IF ( zDiag < 0.0 ) zDiag = SUM( dzSoil(:) )

! Calculate offset soil moisture.
  SELECT CASE ( diagType )
    CASE ( 'avail' )
      vsmcOff(:,:) = vsmcWilt(:,:)
    CASE ( 'total' )
      vsmcOff(:,:) = 0.0
    CASE default
      WRITE(*,*)'ERROR: soilMoistDiag: do not recognise diagType=',TRIM(diagType)
      STOP
  END SELECT

! Loop over soil points.
  DO j=1,soil_pts
    I = SOIL_INDEX(J)

    zbot = 0.0
    DO iz=1,sm_levels

      ztop = zbot
      zbot = zbot + dzSoil(iz)
      IF ( ztop >= zDiag ) EXIT   !  diagnostic is finished for this column

!     Get depth of this layer to consider.
      dz = dzSoil(iz)
      IF ( zbot > zDiag ) dz = zdiag - ztop

!     Add to diagnostic. Restrict values to be >=0.
      outval(i) = outval(i) + RHO_WATER * dz *  &
                      MAX( 0.0, (swet(i,iz)*vsmcSat(i,iz) - vsmcOff(i,iz) ) )

    ENDDO  !  iz
  ENDDO   !  j (soil_pts)

  END FUNCTION soilMoistDiag
!###############################################################################
!###############################################################################
!###############################################################################
! function frozenDepth
! Internal procedure in module output_mod.
! Calculates depth of frozen or thawed ground, according to input.
! Note that it is assumed that the subsurface temperature field (tsoil) is just that,
! whereas at present the top "subsurface layers" of JULES can contain snow.

  FUNCTION frozenDepth( tsoil,diagType ) RESULT( outval )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts,sm_levels

  USE soil_param, ONLY :  &
!  imported arrays with intent(in)
    dzSoil

   USE c_0_dg_c, ONLY :  &
!   imported scalar parameters
     tm
!-------------------------------------------------------------------------------
 IMPLICIT NONE

 REAL ::   &!  function result
    outval(land_pts)  !  function result

  REAL, INTENT(in) ::  &!  in arrays
    tsoil(land_pts,sm_levels)   !  temperature of each subsurface layer (K)

  INTEGER ::  &!  local scalars
    ip,iz,kk   !  work/loop counters

  CHARACTER(len=*), INTENT(in) ::  &! in scalars
    diagType   !  type of diagnostic to calculate
!                'frozen': calculate depth of frozen soil at surface
!                'unfrozen': calculate depth of unfrozen soil at surface

!-------------------------------------------------------------------------------

! Initialise.
  outval(:) = 0.0

  SELECT CASE ( diagType )

!-------------------------------------------------------------------------------
    CASE ( 'frozen' )

!     Get frozen depth from surface.
      DO ip=1,land_pts
        IF ( tsoil(ip,1) < tm ) THEN
!         Find first unfrozen layer.
          kk = sm_levels + 1
          DO iz=2,sm_levels
            IF ( tsoil(ip,iz) >= tm ) THEN
              kk = iz
              EXIT
            ENDIF
          ENDDO
          IF ( kk == sm_levels+1 ) THEN
!           No unfrozen layer was found - whole column is frozen.
            outval(ip) = SUM( dzSoil(:) )
          ELSE
!           Interpolate to estimate depth of zero degC isotherm.
            outval(ip) = SUM( dzSoil(1:kk-1) ) + dzSoil(kk) *  &
                      ( tm-tsoil(ip,kk-1) ) / ( tsoil(ip,kk)-tsoil(ip,kk-1) )
          ENDIF
        ENDIF  !  top layer frozen
      ENDDO  !  ip (land points)

!-------------------------------------------------------------------------------
    CASE ( 'unfrozen' )

!     Get unfrozen depth from surface.
      DO ip=1,land_pts
        IF ( tsoil(ip,1) >= tm ) THEN
!         Find first frozen layer.
          kk = sm_levels + 1
          DO iz=2,sm_levels
            IF ( tsoil(ip,iz) < tm ) THEN
              kk = iz
              EXIT
            ENDIF
          ENDDO
          IF ( kk == sm_levels+1 ) THEN
!           No frozen layer was found - whole column is unfrozen.
            outval(ip) = SUM( dzSoil(:) )
          ELSE
!           Interpolate to estimate depth of zero degC isotherm.
            outval(ip) = SUM( dzSoil(1:kk-1) ) + dzSoil(kk) *  &
                     ( tm-tsoil(ip,kk-1) ) / ( tsoil(ip,kk)-tsoil(ip,kk-1) )
          ENDIF
        ENDIF  !  top layer frozen
      ENDDO  !  ip (land points)
!-------------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'ERROR: frozenDepth: do not recognise diagType=',TRIM(diagType)
      STOP

  END SELECT

  END FUNCTION frozenDepth
!###############################################################################
!###############################################################################
!###############################################################################

! subroutine init_out_count
! Internal procedure in module output_mod.
! Driver for init_out_count_sub1 (which reinitialises output counters).

  SUBROUTINE init_out_count( iout,outActiveTime,outActiveDate )

  USE inout, ONLY : &
!  imported scalar parameters
     periodOneFile  &
!  imported arrays with intent(in)
    ,outDateFlag,outFilePer,outPer,outSamPer   &
!  imported arrays with intent(inout)
    ,outFileStep,outFirstWrite,outStep,outStepSamPer,outWriteCount

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     nspin,spinUp
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    iout            &!  number of output profile
   ,outActiveDate   &!  date (yyyymmdd) when profile is first active
   ,outActiveTime    !  time of day (s) when profile is first active

!----------------------------------------------------------------------------
! Since effectively identical code is needed for setting outStep and outFileStep,
! this is held in init_out_count_sub1 and called separately for each variable.

! outStep (counter within output period) and outStepSamPer.
  CALL init_out_count_sub1( outActiveTime,outActiveDate,outStep(iout),outPer(iout),'outStep'  &
                           ,outStepSamPer(iout),outSamPer(iout) )

! outFileStep (counter of times in output file)
  CALL init_out_count_sub1(  outActiveTime,outActiveDate,outFileStep(iout),outFilePer(iout),'outFileStep' )

!-------------------------------------------------------------------------------

! Set flag showing that no output has been written in this section.
  outFirstWrite(iout) = .TRUE.
! Reset this flag at the start of the main section after spin up, if there was output during spin up and
! all output is to go to the same file.
  IF ( .NOT.spinUp .AND. outDateFlag(iout)==0 .AND. nspin>0 .AND. outFilePer(iout)==periodOneFile ) outFirstWrite(iout) = .FALSE.
! Also reset if we are still in spin up, but all output is to go to one file (or all output
! in spin up to one file).
  IF ( spinUp .AND. outDateFlag(iout)==0 .AND. nspin>0 .AND.  &
         ( outFilePer(iout)==periodOneFile .OR. outFilePer(iout)==-8 ) ) &
            outFirstWrite(iout) = .FALSE.

! Reset counter.
  IF ( outFirstWrite(iout) ) outWriteCount(iout)=0

  END SUBROUTINE init_out_count

!##########################################################################################
!##########################################################################################
!##########################################################################################
! subroutine init_out_count_sub1
! Internal procedure in module output_mod.
! (Re)initialise output counters at the start of output from a new 'section' of the run,
!  if required.

  SUBROUTINE init_out_count_sub1( outActiveTime,outActiveDate,counter,period,callType,counter2,samPer )

  USE inout, ONLY :  &
!  imported scalar parameters
    periodOneFile

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    dateToBits,timeDate,timeDate_diff

  USE timeConst, ONLY : &
    iSecInDay,iSecInHour,iSecInMin

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in scalars
    outActiveDate   &!  date (yyyymmdd) when profile is first active
   ,outActiveTime   &!  time of day (s) when profile is first active
   ,period           !  outPer or outFilePer, according to callType

  INTEGER, INTENT(in), OPTIONAL ::  &!  optional in scalars
    samPer            !  sampling period (number of timesteps; present if callType=outStep)

  INTEGER, INTENT(out) ::  &!  out scalars
    counter          !  outStep or outFileStep, according to callType

  INTEGER, INTENT(out), OPTIONAL ::  &!  optional out scalars
    counter2          !  outStepSamPer (present if callType=outStep)

  INTEGER ::  &!  local SCALARS
    dateMonth,dateOne,dateReady,dateYear,ndy,nhr,nmin  &!  work
   ,nperToReady     &!  number of periods between first active time and time when output is ready
   ,nsec,nstep      &!  work
   ,outActiveDay    &!  day of month from outActiveDate
   ,outDay          &!  days from 00H 1st of month until first output complete
   ,periodUse       &!  period used in calculations. This may be the input period,
!                         or may be the sampling period.
   ,timeReady        !  time of day when first output value is complete

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType   !   flag indicating what variable is being processed
!                    Values are 'outStep' or 'outfileStep'.

!----------------------------------------------------------------------------

! Check that callType is recognised.
  SELECT CASE ( callType )
    CASE ( 'outFileStep', 'outStep' )
    CASE default
      WRITE(*,*)'ERORR: init_out_count_sub1: do not recognise callType=',TRIM(callType)
      STOP
  END SELECT

! Check that optional argument is only present if callType=outStep.
  SELECT CASE ( callType )
    CASE ( 'outStep' )
      IF ( .NOT. PRESENT( counter2 ) ) THEN
        WRITE(*,*)'ERROR: init_out_count_sub1: optional argument counter2 must be'
        WRITE(*,*)'provided for callType=',TRIM(callType)
        STOP
      ENDIF
    CASE ( 'outFileStep' )
      IF ( PRESENT( counter2 ) ) THEN
        WRITE(*,*)'ERROR: init_out_count_sub1: optional argument counter2 should not be'
        WRITE(*,*)'provided for callType=',TRIM(callType)
        STOP
      ENDIF
  END SELECT

! Check we have all or none of optional arguments.
  IF ( PRESENT(counter2) .NEQV. PRESENT(samPer) ) THEN
    WRITE(*,*)'ERORR: init_out_count_sub1: all or none of the optional arguments should be provided.'
    WRITE(*,*)'present(counter2)=',PRESENT(counter2),' present(samPer)=',PRESENT(samPer)
    STOP
  ENDIF

! For "special" file periods (period<0, e.g. monthly or annual), the file counter is not used
! since new files are triggered by events such as the end of a month (rather than an accumulating counter).
! However, for these "special" periods, we do need to initialise the counter for the number of timesteps
! in the output (outStep), since this is also used to decide when to add to time accumulations.
! All counters are also reset after output is written or a file is closed. The counters
! are used for weights when adding to time accumulation, but any errors there are corrected by normalisation.
! So we don't need to do much for the special file periods, but we set counter=0 here anyway.
  counter = 0
  IF ( PRESENT(counter2) ) counter2=0   !  counter2 is always initialised as zero

! If we did start to use outFileStep for the special periods, we would proably need to decide here
! whether counter is to be reset.  e.g. all data to a single file (outFilePer=periodOneFile) would
! only require that outFileStep is set at the start of the first output from the run (not at each section).
! Since outStep (counter within output period) is reset at the start of every "section" of the run,

! Decide what period is used in calculations. Generally this is the input period, but for "special" periods
! we may want to use the sampling period - see below.
  periodUse = period
  IF ( period<0 .AND. callType=='outStep' ) periodUse=samPer
!-------------------------------------------------------------------------------

! Set counter so that outStep=period at time of 1st output.
! e.g. 3 hourly output with start time 04H, first output will be at 06H
! e.g. 2 day output with start day 1, first output will be at OH on day 3 (i.e. end of day 2)
! This is done so that output is in phase with days.
! eg 1 hour averages, and the output starts at say 15H, we ensure the output is
! in phase with days by initialising the counter so that the averages cover 0-2H, 2-4H, 4-6H etc.
! This may mean that the first averaging period is not complete - dealt with elsewhere.
! For special periods (e.g. monthly), the counter is set so that outStep=outSamPer at the time when
! any time accumulation is first added to, i.e. so that mod(outStep,outSamPer)=0 indicates that
! the time accumulation is to be incremented.

  IF ( callType=='outStep' .OR. (callType=='outFileStep' .AND. period>0) ) THEN

!XX Does this really need to be so complicated?! Can't it be simplified into a few "mod"-type operations?

!   First, work out when output will first be generated. That is, when the first output value
!   will be complete (the time of the first snapshot, or the end of the first time average).

    IF ( periodUse*NINT(timeStep) <= iSecInDay ) THEN

!      Period is less than or equal to one day.
!      One day is known to be a multiple of the period.
!      The first output value is complete at the first time in the
!      current day (and after current time) when the time of day is a multiple of the output period.
       nperToReady = CEILING( REAL(outActiveTime) / REAL(periodUse*NINT(timeStep)) )
       timeReady = nperToReady * periodUse * NINT(timeStep)
!      Express the difference in time from current time to timeReady as a number of timesteps.
       nstep = ( timeReady - outActiveTime ) / NINT(timeStep)
       counter = periodUse - nstep

    ELSE

!     Period is more than one day (and is previously tested to ensure it is <=30days).
!     Period is known to be a multiple of one day.
!     The first output value is complete at the first date in the current month (after the current date)
!     when the day of month is a multiple of the output period. If current month is Feb, and period>28 days,
!     we may return a date in March.
!XX     NB This is unnecessary! Matching up to day of month generally only means anything in the first month
!XX     - therafter the output will be output of phase with day of month (not that it matters). Only
!XX        worthwhile keeping in pahse if l_360 and period is a factor of 30. Otherwise, just start output
!XX        1 period from current/start date (or could do at 0H on next day).
      CALL dateToBits( outActiveDate,outActiveDay,dateMonth,dateYear,.FALSE.,'No message')
      nperToReady = CEILING( REAL(outActiveDay) / REAL(periodUse*NINT(timeStep)/iSecInDay) )
      outDay = nperToReady * periodUse * NINT(timeStep) / iSecInDay   !  days from 00H 1st of month until first output complete
!     Get date and time at which first output is complete.
      dateOne = dateYear*10000 + dateMonth*100 + 1   !  date at start of current month
      CALL timeDate( 0,dateOne,outDay*iSecInDay,'sec',.FALSE.,timeReady,dateReady,'Nowt' )
!     Get time difference between now and time when output is ready.
      CALL timeDate_diff( outActiveTime,outActiveDate,timeReady,dateReady,.FALSE.,'nowt' &
                           ,nsec,nmin,nhr,ndy )
      nstep = (ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec) / NINT(timeStep)
      counter = periodUse - nstep

    ENDIF !  per

!   Subtract one, so that the first increment will recover the value set above.
    counter = counter - 1

  ENDIF   !  callType and period

  END SUBROUTINE init_out_count_sub1

!###############################################################################
!###############################################################################

! function loadSnowLayers
! Reverse the order of the last 2 dimensions of a snow layer diagnostic (so that
! each tile can be output as a variable of shape(land_pts,nsmax)), and set
! values where there is an empty layer to a given value.

  FUNCTION loadSnowLayers( fillValue,zrev,nsnow,inval )  &
                           RESULT( outval )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts,nsmax,ntiles

  USE grid_utils, ONLY :  &
!  imported procedures
     reverseCols

  IMPLICIT NONE

! Function result.
  REAL :: outval(land_pts,ntiles*nsmax)

! Scalar arguments with intent(in)
  REAL, INTENT(in) :: fillValue     !  value used where there is no snow
  LOGICAL, INTENT(in) :: zrev       !  flag to reverse order of levels

! Array arguments with intent(in)
  INTEGER, INTENT(in) :: nsnow(land_pts,ntiles)    !  number of layers with snow
  REAL, INTENT(in) :: inval(land_pts,ntiles,nsmax) !  a snow layer variable

! Local scalars.
  INTEGER :: i,j1,j2,t   !  counters and work

! Local arrays.
  REAL :: val(nsmax)     !  work
!-------------------------------------------------------------------------------
  DO i=1,land_pts
    j1 = 1
    DO t=1,ntiles

!     Initialise with the given value, which is kept if there is no snow layer.
      val(:) = fillValue

      IF ( nsnow(i,t) > 0 ) THEN
!       Get nsmax values.
        val(:) = inval(i,t,:)
!       Replace value in empty layers.
        IF ( nsnow(i,t) < nsmax ) val(nsnow(i,t)+1:) = fillValue
      ENDIF

!     Load into the output array.
      j2 = j1 + nsmax - 1
      outval(i,j1:j2) = val(:)
!     Reverse order of levels, if requested.
      IF ( zrev ) outval(i:i,j1:j2)= reverseCols( outval(i:i,j1:j2) )

!     Increment ready for next location.
      j1 = j2 + 1

    ENDDO
  ENDDO

  END FUNCTION loadSnowLayers

!################################################################################
!################################################################################

  SUBROUTINE create_annotation( iout,callType,nvarCtl  &
                ,optionsLine,titleLine,zdefLine,annot,tdefAnnot,varCtlLine,varNcAtt )

! Create metadata and annotation that is useful for a GrADS ctl file,
! netCDF attributes, etc. As this information can be used in more than one
! location, it's easier to create it once and share.

  USE ancil_info, ONLY : &
!  imported scalars with intent(in)
     dim_cs1,nsmax,nsoil=>sm_levels,ntiles

  USE inout, ONLY :  &
!  imported scalar parameters
     formatBin,formatLen,formatNc,periodOneFile  &
!  imported scalars with intent(in)
    ,gradsNc,outEndian,outFormat,outGrads,runID,usePseudo,yrevOut  &
    ,zrevOutSoil,zrevOutSnow  &
!  imported arrays with intent(in)
    ,havePFT,haveSnow,haveSoil,haveTile  &
    ,haveType,nlevMaxCtl,nvarOut,outCompress,outDateFlag,outFilePer  &
    ,outGridNxy,outLLorder,outName,outPer  &
    ,outTemplate,pointsOut,pointsFlag,rgProfile,rpProfile,snapProfile  &
    ,taccumVar,tmeanProfile,tmeanVar,varDesc,varName  &
    ,varNum,varPos,varType,varUnitsList

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     nnvg,npft,ntype

  USE nvegparm, ONLY : &
!  imported arrays with nitent(in)
     nvgName

  USE pftparm, ONLY :  &
!  imported arrays with intent(in)
     pftName

  USE soil_param, ONLY :  &
!  imported arrays with intent(in)
     dzsoil

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in
     nspin,spinUp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_aggregate

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timeStep  &
!  imported arrays with intent(in)
    ,dateMainRun,dateSpin

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------

  INTEGER, INTENT(in) :: iout  !  current output profile

  CHARACTER(len=*), INTENT(in) :: callType  !  type (location) of call
!                  This is used to indicate what type of data are to be output.
!         compress: means this call is for the
!                  supplementary data file that describes the mapping used
!                  when compressing output (e.g. the supplementary file referred
!                  to in the GrADS pdef entry) .Note this is FALSE when
!                  just writing a ctl file that contains a pdef entry.
!                  Here it is assumed that this data file is binary (GrADS).
!         diag: means this call is related to the writing of diagnostic output
!         dump: means this call is related to the writing of a dump (restart) file

!-------------------------------------------------------------------------------
! Scalar arguments with intent(out)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(out) :: nvarCtl !  number of variables, as represented in ctl file
  CHARACTER(len=*), INTENT(out) :: optionsLine  !  GrADS options line
  CHARACTER(len=*), INTENT(out) :: titleLine    !  GrADS title line
  CHARACTER(len=*), INTENT(out) :: zdefLine     !  GrADS zdef line

!-------------------------------------------------------------------------------
! Array arguments with intent(out)
!-------------------------------------------------------------------------------
  CHARACTER(len=*), INTENT(out) :: annot(:)    !  annotation lines
  CHARACTER(len=*), INTENT(out) :: tdefAnnot(:)!  annotation lines related to
!                                                   GrADS tdef
  CHARACTER(len=*), INTENT(out) :: varCtlLine(:)  !  details of variables
  CHARACTER(len=*), INTENT(out) :: varNcAtt(:,:)  !  variable attributes for
!             netCDF files. This has been rather hacked at present - in particular
!             note that snow layer variables can have ntiles entries in varNcAtt.
!             On exit, both varCtlLine and varNcAtt are expected to have nvarCtl
!             entries (lines).

!-------------------------------------------------------------------------------
! Local scalar variables
!-------------------------------------------------------------------------------
  INTEGER :: dz,i,ia,iline,ivar,iz,iz1,iz2,jvar,kvar,nlevCtl   !  work/loop counters
  INTEGER :: nvarRep  !  used if a variable is represented by more than one variable!
  INTEGER :: nzMax    !  number of z levels needed

  LOGICAL :: haveScVar    !  T means there is a soil carbon pool variable to be output
!                            Only T if there is more than 1 pool.
  LOGICAL :: useSoilLevs  !  T means zdef will use soil levels, F means uses an index
! For netCDF files, which use the name=>Name syntax in the ctl file, a robust length
! for name is 2*LEN(varname)+4 - but this is currently jolly long thanks to sufLen,
! and much longer than is needed in most cases (which don't use sufLen).
! On writing we use adjustl(name) to get entries that start on left and align, but
! this looks lousy with long gaps if name is long.
! So....this hardwired len for name is an unplesant hack that works most of the time...
! but if you use a long suffix it might fail (particularly for netCDF output) - in
! that case increase the length!
! An alternative would be to calculate and store all names, identify the longest,
! and then output allowing for this length, but that would require more work
! than I fancy right now (but nothing too hard)!
!  CHARACTER(len=2*LEN(varname)+4) :: name  !  work
  CHARACTER(len=LEN(varname)+2) :: name  !  work
  CHARACTER(len=formatLen) :: fileFormat !  format of output file
  CHARACTER(len=7) :: axesCtl    !  axes for netCDF ctl file (long enough
!                                   for "t,z,y,x")
  CHARACTER(len=20) :: axes      !  axes as written to netCDF file attribute
  CHARACTER(len=200) :: cwork,cwork1    !  work
  CHARACTER(len=200) :: desc     !  a description of a variable

!-------------------------------------------------------------------------------
! Local array variables
!-------------------------------------------------------------------------------
  REAL :: zsoil(nsoil)  !  representative depths for each soil layer (m)

!-------------------------------------------------------------------------------
! Check the cal type is recognised.
!-------------------------------------------------------------------------------
  SELECT CASE ( callType )
    CASE ( 'compress', 'diag', 'dump' )
!     OK
    CASE default
      WRITE(*,*)'ERROR: create_annotation: do nto recognise callType='  &
                  ,TRIM(callType)
      STOP
  END SELECT

!-------------------------------------------------------------------------------
! Initialise annotation as empty.
!-------------------------------------------------------------------------------
  ia = 0
  tdefAnnot(:) = ''
  annot(:) = ''
  optionsLine = ''
  varCtlLine(:) = ''
  varNcAtt(:,:) = ''

!-------------------------------------------------------------------------------
! Sort out how to represent vertical levels.
!-------------------------------------------------------------------------------
  SELECT CASE ( callType )
    CASE ( 'compress' )
      nzMax = 1
    CASE ( 'diag' )
      nzMax = nlevMaxCtl(iout)
    CASE ( 'dump' )
      nzMax = 1  !  cxyz
  END SELECT

!-------------------------------------------------------------------------------
! Establish what types of variable we have in this profile - this affects how
! z levels are represented and whether we write names of surface types.
! Also work out how many variables will be represented in ctl file.
!-------------------------------------------------------------------------------
  haveScVar = .FALSE.
  useSoilLevs = .FALSE.
  SELECT CASE ( callType )
    CASE ( 'diag' )
      fileFormat = outFormat
      nvarCtl = nvarOut(iout)
      DO ivar=1,nvarOut(iout)
        SELECT CASE ( varType(varPos(iout,ivar)) )
          CASE ( 'SC' )
            IF ( dim_cs1 > 1 ) haveScVar = .TRUE.
          CASE ( 'SN' )
!           Snow variables are represented separately for each tile, so
!           increment nvarCtl. This is not the case for netCDF files that are
!           not to be read by GrADS (and for which the ctl file cannot work).
!           Make sure attributes etc are consistent!
            IF ( outGrads ) nvarCtl = nvarCtl + ntiles - 1
        END SELECT
      ENDDO

!     If we have a soil variable, assume that this will be used to describe
!     levels, but then check if there is any other type of multi-level
!     variable (which would mean we would not use soil levels).
!     I'm assuming that a "multi-level" variable type has more than one
!     level (e.g. not ntype=1).
!      IF ( haveSoilVar ) useSoilLevs = .TRUE.
!      IF ( useSoilLevs ) THEN
!        DO ivar=1,nvarOut(iout)
!          SELECT CASE ( varType(varPos(iout,ivar)) )
!            CASE ( 'PF','SN','TI','TY' )
!              useSoilLevs = .FALSE.
!          END SELECT
!        ENDDO
!        IF ( haveScVar ) useSoilLevs = .FALSE.
!      ENDIF
!     Override useSoilLevs - it can be confusing if levels change between output files
!     only because different variables are present. The code is left...for now.
      useSoilLevs = .FALSE.

    CASE ( 'compress' )

      IF ( outGrads .OR. outFormat==formatNc ) THEN
        fileFormat = formatBin
      ELSE
        fileFormat = outFormat
      ENDIF
      nvarCtl = 3

  END SELECT  !  callType

!-------------------------------------------------------------------------------
! Create OPTIONS and TITLE lines for ctl.
!-------------------------------------------------------------------------------
  IF ( fileFormat == formatBin ) optionsLine = TRIM(outEndian)
  SELECT CASE ( callType )
    CASE ( 'compress' )
      titleLine = 'JULES run ' // TRIM(runID)    &
           // ', output compression (pdef) data for output profile '  &
           // TRIM(outname(iout))
    CASE ( 'diag' )
      IF ( outTemplate(iout) ) optionsLine = TRIM(optionsLine) // ' template'
      IF ( yrevOut .AND. .NOT.outCompress(iout) ) optionsLine =  &
               TRIM(optionsLine) // ' yrev'
      titleLine = 'JULES run ' // TRIM(runID) // ', output profile '  &
                    // TRIM(outname(iout))
    CASE ( 'dump' )
!      cxyz dump
  END SELECT

!-------------------------------------------------------------------------------
! Create a ZDEF line for ctl file. cxyz don;t bother if not writing a ctl
!-------------------------------------------------------------------------------
  IF ( .NOT. useSoilLevs ) THEN
!   Describe vertical levels in terms of an index 1,2,...
    IF ( nzMax == 1 ) THEN
      WRITE(zdefLine,"('zdef   1 linear 1 1')")
    ELSE
      WRITE(zdefLine,"('zdef ',i3,' linear 1 1')") nlevMaxCtl(iout)
    ENDIF
  ELSE
!   Describe vertical levels in terms of soil layers.
    WRITE(zdefLine,"(a)")'# zdef based on soil layer midpoints'
    DO iz=1,nsoil
      zsoil(iz) = -1.0 * ( SUM(dzsoil(1:iz-1)) + 0.5*dzsoil(iz) )
    ENDDO
    iz1=1; iz2=nsoil; dz=1
    IF ( zrevOutSoil ) THEN; iz1=nsoil; iz2=1; dz=-1; ENDIF
    WRITE(zdefLine,"('zdef ',i3,' levels ',(15f7.2))") nsoil,(zsoil(iz), iz=iz1,iz2,dz)
  ENDIF  !  useSoilLevs

!-------------------------------------------------------------------------------
! Create tdef and associated information.
!-------------------------------------------------------------------------------
  IF ( callType == 'diag' ) THEN
!   Add a comment if times will not be correctly represented.
    i=0
    IF ( outFilePer(iout)==periodOneFile .AND. outDateFlag(iout)==0 .AND. nspin>0  &
!                      one file and output through spin up
         .AND. .NOT. (nspin==1 .AND. dateSpin(2)==dateMainRun(1)) )   &
!                      but not one cycle of spin up and main run to follow
      THEN
      i=i+1
      tdefAnnot(i) = '# Dates cannot be correctly represented for this file, since'
      i=i+1
      tdefAnnot(i) = '# the same times are repeated during spin up.'
    ENDIF
    IF ( outFilePer(iout)==-8 .AND. spinUp .AND. nspin/=1 ) THEN
      i=i+1
      tdefAnnot(i) = '# Dates cannot be correctly represented for this file, since'
      i=i+1
      tdefAnnot(i) = '# the same times are repeated during spin up.'
    ENDIF
    IF ( outDateFlag(iout) == -2 ) THEN
      i=i+1
      tdefAnnot(i) =  &
     '# This file contains data for the start time of the run (initial state) only.'
    ENDIF
  ENDIF  !  callType

!-------------------------------------------------------------------------------
! Some annotation to give extra details.
!-------------------------------------------------------------------------------
  SELECT CASE ( callType )
    CASE ( 'diag' )

      IF ( .NOT. rpProfile(iout) ) THEN

        IF ( rgProfile(iout) ) THEN
          ia=ia+1; annot(ia) = '# This is output on the routing grid.'
        ENDIF

        SELECT CASE ( pointsFlag(iout,1) )
          CASE ( 0 )
            ia=ia+1; annot(ia) = '# All points in model grid are included in output.'
          CASE ( 1 )
            ia=ia+1; annot(ia) = '# This output is from a subarea of the model grid.'
          CASE ( 2 )
            ia=ia+1; annot(ia) = '# The points to be output were given in a list.'
        END SELECT

        SELECT CASE ( pointsFlag(iout,2) )
          CASE ( 0 )
            ia=ia+1; annot(ia) = '# The output grid is the same as the model grid.'
          CASE ( 1 )
            ia=ia+1; annot(ia) = '# The output grid is the same as the input grid.'
          CASE ( 2 )
            ia=ia+1; annot(ia) = '# The output grid was described via the run control file.'
          CASE ( 3 )
            ia=ia+1; annot(ia) = '# The output grid is the smallest rectangle that contained all'
            ia=ia+1; annot(ia) = '# points. The rectangle and order of points was determined by'
            IF ( outLLorder(iout) ) THEN
              ia=ia+1; annot(ia) = '# the latitude and longitude of each point.'
            ELSE
              ia=ia+1; annot(ia) = '# the row and column position in the input grid of each point.'
          ENDIF
        END SELECT

        IF ( .NOT. outCompress(iout) ) THEN
          IF ( pointsOut(iout) /= outGridNxy(iout,1)*outGridNxy(iout,2) ) THEN
            ia=ia+1
            WRITE(annot(ia),"(a,i5)")'# Number of undefined (padding) points='  &
                           ,outGridNxy(iout,1)*outGridNxy(iout,2)-pointsOut(iout)
          ENDIF
        ELSE
          ia=ia+1; annot(ia) = '# Output is a vector of model points that can be "scattered"'
          ia=ia+1; annot(ia) = '# across a (potentially) larger grid.'
          ia=ia+1; WRITE(annot(ia),"(a,i6)")'# Number of points in vector=',pointsOut(iout)
        ENDIF

      ELSE

!       rpOnlyProfile=TRUE.
        ia=ia+1; annot(ia) = '# All output variables are for points on routing grid.'

      ENDIF  ! rpProfile

!-------------------------------------------------------------------------------
!   If there are soil layer variables, describe order of levels.
!-------------------------------------------------------------------------------
    IF ( haveSoil(iout) ) THEN
      IF ( zrevOutSoil ) THEN
        ia=ia+1; annot(ia) = '# Soil layer variables are written bottom to top.'
      ELSE
        ia=ia+1; annot(ia) = '# Soil layer variables are top to bottom.'
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!   If there is a variable for more than one soil C pool, describe this.
!-------------------------------------------------------------------------------
    IF ( haveScVar ) THEN
      ia=ia+1; annot(ia) = '# Soil carbon pools: DPM, RPM, biomass, humus'
    ENDIF

!-------------------------------------------------------------------------------
!   If there are snow layer variables, describe order of levels.
!-------------------------------------------------------------------------------
    IF ( haveSnow(iout) ) THEN
      IF ( zrevOutSnow ) THEN
        ia=ia+1; annot(ia) = '# Snow layer variables are written top to bottom.'
      ELSE
        ia=ia+1; annot(ia) = '# Snow layer variables are bottom to top.'
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!   Write names of surface types. For now, assume we can get all these into single lines!
!-------------------------------------------------------------------------------
    IF ( havePFT(iout) .OR. haveTile(iout) .OR. haveType(iout) ) THEN
!     Get names and numbers into a string. For now, assume we can get all this into one line!
      WRITE(cwork,"(a)") '# PFTs:'
      DO iz=1,npft
        WRITE(cwork1,"(tr1,i2,tr1,a)") iz,TRIM(pftName(iz))
        cwork = TRIM(cwork) // TRIM(cwork1)
      ENDDO
      ia=ia+1; annot(ia) = TRIM( cwork )
    ENDIF
    IF ( haveTile(iout) .OR. haveType(iout) ) THEN
!     Write names of non-veg surface types, having already written PFT names.
!     Get names and numbers into a string. For now, assume we can get all this into one line!
      WRITE(cwork,"(a)") '# Non-veg types:'
      DO iz=1,nnvg
        WRITE(cwork1,"(tr1,i2,tr1,a)") npft+iz,TRIM(nvgName(iz))
        cwork = TRIM(cwork) // TRIM(cwork1)
      ENDDO
      ia=ia+1; annot(ia) = TRIM( cwork )
    ENDIF

    IF ( callType == 'diag' ) THEN  !  cxyz WRONG as we're in that IF already
!-------------------------------------------------------------------------------
!     Describe time-averaging annotation for variables.
!-------------------------------------------------------------------------------
      IF ( snapProfile(iout) ) THEN
        ia=ia+1
        annot(ia) = '# S denotes snapshot value at given time.'
      ENDIF

      IF ( tmeanProfile(iout) ) THEN
        IF ( outPer(iout)>0 ) THEN
          ia=ia+1
          WRITE(annot(ia),"(a,i6,a)") '# M denotes backward time average over '  &
            ,outPer(iout)*NINT(timeStep),'s, ending at given time.'
          ia=ia+1
          WRITE(annot(ia),"(a,i6,a)") '# A denotes time-accumulation*timestep (s), over '  &
             ,outPer(iout)*NINT(timeStep),'s, ending at given time.'
        ELSE
          ia=ia+1
          annot(ia) = '# M denotes backward time average over month/year, ending at given time.'
          ia=ia+1
          annot(ia) = '# A denotes time-accumulation*timestep (s), ending at given time.'
        ENDIF
      ENDIF
    ENDIF  !  callType

  END SELECT

!-------------------------------------------------------------------------------
! Write a warning if snow layer variables cannot be handled correctly by GrADS.
!-------------------------------------------------------------------------------
  IF ( haveSnow(iout) .AND. outFormat==formatNc .AND. .NOT.gradsNc ) THEN
    ia=ia+1; annot(ia) = '########################################&
                         &################################'
    ia=ia+1; annot(ia) = '#### NB gradsNc=F. This ctl file cannot &
                         &represent snow layer variables.'
    ia=ia+1; annot(ia) = '########################################&
                         &################################'
  ENDIF


!-------------------------------------------------------------------------------
! Get details of each variable.
!-------------------------------------------------------------------------------
  iline = 0
  SELECT CASE ( callType )

    CASE ( 'diag' )

      DO ivar=1,nvarOut(iout)
        jvar = varPos(iout,ivar)
        kvar = varNum(jvar)    !  position in master list
        cwork = 'S'  ! snapshot variable
        IF ( taccumVar(jvar) ) cwork = 'A'
        IF ( tmeanVar(jvar) ) cwork = 'M'

!       Get number of levels to use and axes for netCDF.
        nvarRep = 1
        SELECT CASE ( varType(jvar) )
          CASE ( 'LA', 'RG', 'RP', 'SI' )
            nlevCtl = 0
            axes = 'TYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,x'
            ELSE
              axesCtl = 't,y,x'
            ENDIF
          CASE ( 'PF' )
            nlevCtl = npft
            axes = 'TZYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,z,x'
            ELSE
              axesCtl = 't,z,y,x'
            ENDIF
          CASE ( 'SC' )
            nlevCtl = dim_cs1
            axes = 'TZYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,z,x'
            ELSE
              axesCtl = 't,z,y,x'
            ENDIF
          CASE ( 'SN' )
            IF ( .NOT. usePseudo ) THEN
!             For snow layer variables, treat each tile as a separate variable.
!             This is for ease of use with GrADS, which generally can't cope with all
!             these dimensions!
              nvarRep = ntiles
              nlevCtl = nsmax
              axes = 'TZYX'
              IF ( outCompress(iout) ) THEN
                axesCtl = 't,z,x'
              ELSE
                axesCtl = 't,z,y,x'
              ENDIF
            ELSE
!             netCDF output that is not required to be readable by GrADS (gradsNc=F)
!             Write a single variable.
!             nlevCtl (number of levels as in ctl is irrelevant as ctl will not work!)
              nlevCtl = nsmax
              axes = 'TPZYX'
!             Although ctl will not work for SN variables, it can still be used
!             for other types of variable. So just claim that axes are something
!             reasonable for SN variables.
              IF ( outCompress(iout) ) THEN
                axesCtl = 't,z,x'    !   actually t,p,z,x
              ELSE
                axesCtl = 't,z,y,x'  !   actualy t,p,z,y,x
              ENDIF
            ENDIF   !  usePseudo
          CASE ( 'SO' )
            nlevCtl = nsoil
            axes = 'TZYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,z,x'
            ELSE
              axesCtl = 't,z,y,x'
            ENDIF
          CASE ( 'TI' )
            nlevCtl = ntiles
            axes = 'TZYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,z,x'
            ELSE
              axesCtl = 't,z,y,x'
            ENDIF
          CASE ( 'TY' )
            nlevCtl = ntype
            axes = 'TZYX'
            IF ( outCompress(iout) ) THEN
              axesCtl = 't,z,x'
            ELSE
              axesCtl = 't,z,y,x'
            ENDIF
          CASE default
            WRITE(*,*)'ERROR: create_annotation #2: varType=',varType(jvar)
            WRITE(*,*)'No code for this varType.'
            WRITE(*,*)'Stopping in create_annotation'
            STOP
        END SELECT  !  varType

!       If not netCDF output, replaces axesCtl with the default '99' code for units.
        IF ( fileFormat /= formatNc ) axesCtl = '99'

        DO iz=1,nvarRep
          iline = iline + 1

          IF ( varType(jvar) /= 'SN' ) THEN
            name = varName(jvar)
            desc = varDesc(jvar)
          ELSEIF ( usePseudo ) THEN
            name = varName(jvar)
!           Append warning to desc.
            desc = varDesc(jvar)
            desc = TRIM(desc) // ' ## SEE gradsNC=F WARNING ABOVE!'
          ELSE
!           For snow layer variables, add tile number to name, add type to description.
!           Becoming cumbersome as bits have been bolted on!
            IF ( nvarRep < 10 ) THEN
              WRITE(cwork1,"(i1.1)") iz
            ELSE
              WRITE(cwork1,"(i2.2)") iz
            ENDIF
            name = TRIM(varName(ivar)) // TRIM(cwork1)
            IF ( l_aggregate ) THEN
              cwork1 = 'aggregate tile'
            ELSEIF ( iz <= npft ) THEN
              cwork1 = TRIM( pftName(iz) ) // ' tile'
            ELSE
              cwork1 = TRIM( nvgName(iz-npft) ) // ' tile'
            ENDIF
            desc = TRIM(varDesc(jvar)) // ' ' // TRIM(cwork1)
          ENDIF

!         For netCDF files, GrADS at v1.9b4 seems to need the name=>new_name
!         syntax, even if name=new_name!
          IF ( fileFormat == formatNc ) THEN
            name = TRIM(name) // '=>' // TRIM(name)
            varNcAtt(iline,1) = axes
            varNcAtt(iline,2) = varUnitsList(kvar)
            varNcAtt(iline,3) = name
            varNcAtt(iline,4) = desc
          ENDIF

          WRITE(varCtlLine(iline),"(a,tr1,i2,tr1,a,tr2,a1,tr2,a,tr1,a1,a,a1)")  &
                        ADJUSTL(name),nlevCtl,ADJUSTL(axesCtl),cwork(1:1)  &
                       ,TRIM(desc),'(',TRIM(varUnitsList(kvar)),')'

        ENDDO  !  iz

      ENDDO  !  ivar

    CASE ( 'compress' )

!     This always writes the same variables, and file format is not netCDF.
      i=0
      i=i+1; WRITE(varCtlLine(i),"('ind  0  -1,40,4  index (position in vector)')")
      i=i+1; WRITE(varCtlLine(i),"('wt   0  99       weights (1 at points in vector, else 0)')")
      i=i+1; WRITE(varCtlLine(i),"('rot  0  99       wind rotation values (all undef)')")

  END SELECT  !  callType


  END SUBROUTINE create_annotation

!###############################################################################
!###############################################################################
!###############################################################################
  END MODULE output_mod
