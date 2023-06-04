!###############################################################################
!###############################################################################
!
! subroutine init_drive
! Read details of driving data and initialise.
!
!###############################################################################
!###############################################################################
SUBROUTINE Init_Drive

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY :  &!
!  imported arrays with intent(out)
     z1_uv,z1_tq

  USE drive_io_vars, ONLY :  &
!   imported scalar parameters
      ndriveExtraMax  &
!   imported scalars with intent(out)
     ,byteSwapDrive,diffFracConst  &
     ,driveDataPer,driveDataStep,driveDataStepInit,driveDate,driveDateInit  &
     ,driveEndTime,driveFile,driveFilePer,driveFileStep,driveFileTemplate  &
     ,driveFormat,driveResetStep,driveResetStepPrev  &
     ,driveTemplateT,driveTemplateDate,driveTemplateUnits,driveTemplateTime  &
     ,driveTemplateV,driveTime,driveTimeInit  &
     ,ioPrecipType,io_rad_type  &
     ,ndriveDataTime,ndriveExtra,ndriveFileTime,ndriveUnit  &
     ,ndriveHeaderField,ndriveHeaderFile,ndriveHeaderTime  &
     ,ndriveVar,ndriveVarIn,ndriveVarMax  &
     ,nfieldDriveFile,noNewLineDrive,notNextDrive  &
     ,tForCRain,tForCRain,tForSnow,useWGen    &
!   imported arrays with intent(out)
     ,driveFileDate,driveFileName,driveFileTime  &
     ,driveVarInterp,driveVarInSort  &
     ,driveTimeIndex,driveUnit,driveFileOnUnit  &
     ,driveVarNameSDF,driveVarUse

  USE file_utils, ONLY : &
!   imported procedures
     closeFile,fileUnit,findTag,openFile

  USE soil_param, ONLY :  &
!  imported scalars with intent(out)
    conFrac

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatBin,formatNc,formatLen,formatPP,periodAnn,periodMon,periodOneFile  &
    ,jinUnit,tagAscBin,tagNc  &
!  imported scalars with intent(in)
    ,echo,nxIn,nyIn

  USE jules_netcdf, ONLY :  &
!  imported procedures
     checkNcType  &
!  imported scalars with intent(out)
    ,ncTypeDrive

  USE misc_utils, ONLY :  &
!  imported procedures
     check_template,read_list

  USE spin_mod, ONLY :  &
!   imported scalars with intent(in)
     spinUp

  USE switches, ONLY :  &
!   imported scalars with intent(in)
     l_360,l_point_data,routeOnly,l_imogen

  USE time_loc, ONLY :  &
!   imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :  &
!   imported procedures
     chhmmss_to_s,timeDate,timeDate_cmp

  USE timeConst, ONLY : &
     iSecInDay

  USE update_mod, ONLY :  &
!   imported procedures
     data_init,calc_reset_step

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    i,ierr,ipos,iunit,ivar   &!  work
   ,ndriveVarNeed            &!  number of driving variables for chosen configuration
!                                  Only used for diagnostic message to screen.
   ,t                        &!  work
   ,tmpDataStepMax            !  work

  INTEGER ::  &!  local ARRAYS
    tmpFlag(ndriveVarMax)   !  values of driveVarFlag as read in

  REAL ::  &!  local scalars
    z1_uv_val,z1_tq_val    ! values of z1_uv and z1_tq that are used at all points on grid

  LOGICAL ::  &!  local SCALARS
    haveAve     &!  T means that one or more variable in input data is a time-average
!                     that places extra restriction on timestep
   ,readList    &!  T means read a list of file names, F means read a single file name
   ,templateT,templateV  &!  work
   ,ioWindSpeed          &!  T means that the windspeed is input
!                            F means 2 components of wind are input
   ,sdfFile              &!  T means data will come from an SDF (self-describing file)
   ,useDiffRad            !  T means diffuse radiation is input
!                            F means diffuse radiation is set to a constant
!                            fraction (diffFracConst)



  LOGICAL ::  &! local ARRAYS
    needTime(-1:2)       !  flag indicating which times of data are required
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  CHARACTER(len=8) ::  &!  local SCALARS
    cdriveTime   !  time string

  CHARACTER(len=LEN(driveVarNameSDF)) ::  &!  local ARRAYS
    driveVarName(ndriveVarMax)   &!  names of all possible forcing variables
   ,driveVarFileName(ndriveVarMax)  &!  the part of the file name that changes between
!                  files that hold different variables, if driveTemplateV=TRUE.
   ,tmpName(ndriveVarMax)        &!  work: names of forcing variables, as read in.
   ,tmpNameSDF(ndriveVarMax)     &!  work: names of forcing variables for SDF files
   ,tmpNameFile(ndriveVarMax)     !  work: names used in file names, as read in

  CHARACTER(len=LEN(driveVarInterp)) ::  &!  local ARRAYS
    tmpInterp(ndriveVarMax)         !  work: values of driveVarInterp, as read in

  CHARACTER(len=150) ::  &!  local SCALARS
    listFile  !  work

  CHARACTER(len=70) ::  &!  local arrays
    varDesc(ndriveVarMax)  !  a description of each possible driving variable
!------------------------------------------------------------------------------------------


  if (echo) WRITE(*,"(50('-'),/,a)") 'init_drive'

!-------------------------------------------------------------------------------
! Initialise.
!-------------------------------------------------------------------------------
  byteSwapDrive = .FALSE.        !   data will not be byteswapped
  driveTemplateT = .FALSE.       !   no time-templated file names
  driveTemplateV = .FALSE.       !   no var-templated file names

!-------------------------------------------------------------------------------
! Hardwire so that weather generator is not available in this version.
! Just changing this line will not successfully activate the weather generator
!  - other variables have been commented out too and the model WILL fail!
!-------------------------------------------------------------------------------
  useWgen = .FALSE.

!-------------------------------------------------------------------------------
! Read information from run control file.
!-------------------------------------------------------------------------------

! Locate the start of this section in file.
  CALL findTag( jinUnit,'init_drive','>INIT_DRIVE' )

  READ(jinUnit,*) driveDataPer  !DSM Nao utilizado
  driveDataPer=timestep  !DSM o CCATT estah fornecendo os dados a cada timestep

  READ(jinUnit,*) ndriveFileTime,driveFilePer
! Establish if a list of files is to be read.
  READ(jinUnit,*) readList

!-------------------------------------------------------------------------------
! Allocate space for file names and times.
!-------------------------------------------------------------------------------
  CALL allocate_arrays( 'init_drive 1' )

!-------------------------------------------------------------------------------
! Read the name of file listing file names, or the name of the only file.
!-------------------------------------------------------------------------------
! If more than one file is indicated, insist that a list is read.
  IF ( ndriveFileTime>1 .AND. .NOT.readList ) THEN
    WRITE(*,*)'ERROR: init_drive: names of more than one file must be read from a list file.'
    WRITE(*,*)'Change readList to TRUE.'
    STOP
  ENDIF

  IF ( readList ) THEN
    READ(jinUnit,*) listFile
  ELSE
    READ(jinUnit,*) driveFileName(1)
  ENDIF

!-------------------------------------------------------------------------------
  READ(jinUnit,*) driveFileDate(1),cdriveTime
  READ(jinUnit,*) driveEndTime
  READ(jinUnit,*) driveFormat

  READ(jinUnit,*) ioPrecipType,l_point_data
  READ(jinUnit,*) tForSnow
  READ(jinUnit,*) tForCRain,conFrac
  READ(jinUnit,*) io_rad_type,ioWindSpeed
  READ(jinUnit,*) useDiffRad,diffFracConst

! Read a single pair of height values.
  READ(jinUnit,*) z1_uv_val,z1_tq_val

! Read number of "extra" variables.
  READ(jinUnit,*) ndriveExtra

!-------------------------------------------------------------------------------
! Check timesteps for weather generator.
!-------------------------------------------------------------------------------
  IF ( useWgen ) THEN
!   Weather generator is not needed if timestep >= data period.
    IF ( timeStep >= driveDataPer ) THEN
      WRITE(*,*)'ERROR: timeStep is >= interval (period) of driving data '
      WRITE(*,*)'No need for weather generator. Set useWgen to FALSE.'
      WRITE(*,*)'Stopping in init_drive'
      STOP
    ENDIF
!   Weather generator is coded for daily data.
    IF ( driveDataPer /= iSecInDay ) THEN
      WRITE(*,*)'ERROR: weather generator requires daily input data.'
      WRITE(*,*)'Stopping in init_drive'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Check values of flags.
!-------------------------------------------------------------------------------
  IF ( .NOT. routeOnly ) THEN
    IF ( ioPrecipType<1 .OR. ioPrecipType>4 ) THEN
      WRITE(*,*)'ERROR: init_drive_var: ioPrecipType must be in range 1 to 4.'
      STOP
    ENDIF
    IF ( l_point_data .AND. ( ioPrecipType==3 .OR. ioPrecipType==4 ) ) THEN
      WRITE(*,*)'ERROR: init_drive_var: l_point_data cannot be used with ioPrecipType=3 or 4'
      STOP
    ENDIF
    IF ( ECHO .AND. ( ioPrecipType/=1 .AND. ioPrecipType/=5 ) ) &
       WRITE(*,*)'NOTE: ioPrecipType=1 or 5 means tForSnow is not used.'
    IF ( ECHO .AND. ( ioPrecipType==3 .OR. ioPrecipType==4 ) )  &
       WRITE(*,*)'NOTE: ioPrecipType=3 or 4: tForCRain is not used.'
    IF ( ECHO .AND. l_point_data ) &
       WRITE(*,*)'NOTE: l_point_data: tForCRain is ignored and there is no convective precip.'
    IF ( io_rad_type<1 .OR. io_rad_type>3 ) THEN
      WRITE(*,*)'ERROR: init_drive_var: io_rad_type must be in range 1 to 3.'
      STOP
    ENDIF
    IF ( useDiffRad .AND. io_rad_type > 2 ) THEN
      WRITE(*,*)'ERROR: init_drive_var: useDiffRad requires io_rad_type=1 or 2.'
      WRITE(*,*)'This made coding easier, but could be circumvented.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Override certain flags if weather generator is used.
!-------------------------------------------------------------------------------
  IF ( useWgen ) THEN
    IF ( ioPrecipType /= 5 ) THEN
      WRITE(*,*)'Resetting ioPrecipType to 5 for use with weather generator.'
      ioPrecipType = 5
    ENDIF
    IF ( io_rad_type /= 1 ) THEN
      WRITE(*,*)'Resetting io_rad_type to 1 for use with weather generator.'
      io_rad_type = 1
    ENDIF
    IF ( .NOT. ioWindSpeed ) THEN
      WRITE(*,*)'Resetting ioWindSpeed to TRUE for use with weather generator.'
      ioWindSpeed = .TRUE.
    ENDIF
    IF ( l_point_data ) THEN
      WRITE(*,*)'ERROR: l_point_data not allowed with weather generator.'
      WRITE(*,*)'Mainly because it seems unlikely. And I''m in a hurry.'
      WRITE(*,*)'Stopping in INIT_DRIVE.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Check temperatures used with precip are consistent. Not used with point data.
!-------------------------------------------------------------------------------
  IF ( .NOT.l_point_data .AND.  &
       ( ioPrecipType==1 .OR. ioPrecipType==2 .OR. ioPrecipType==5) ) THEN
    IF ( tForCRain <= tForSnow ) THEN
      WRITE(*,*)'ERROR: init_drive: tForCRain <= tForSnow.'
      WRITE(*,*)'ForSnow=',tForSnow,' tForCRain=',tForCRain
      WRITE(*,*)'Precipitation is snow when T<=tForSnow.'
      WRITE(*,*)'Rainfall is convective when T>=tForCRain.'
      WRITE(*,*)'We require tForCRain > tForSnow, so that precip is liquid at tForCRain.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Check that number of "extra" variables is within range.
!-------------------------------------------------------------------------------
  IF ( ndriveExtra > ndriveExtraMax ) THEN
    WRITE(*,*)'ERROR:  ndriveExtra > ndriveExtraMax.'
    WRITE(*,*)'Stopping in INIT_DRIVE.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Note that the model and i/o specifications have already determined the number
! of variables to be input, but we read all variables provided. Any mistake is
! (hopefully) detected later.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.
!-------------------------------------------------------------------------------
  SELECT CASE ( driveFormat )

    CASE ( formatAsc,formatBin,formatPP )
      sdfFile = .FALSE.
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_drive',tagAscBin,preInit=.TRUE. )
      READ(jinUnit,*) byteSwapDrive
      READ(jinUnit,*) nfieldDriveFile
      READ(jinUnit,*) ndriveHeaderFile,ndriveHeaderTime,ndriveHeaderField
      READ(jinUnit,*) noNewLineDrive

!     Read details of variables. Read 4 values from a blank-delimited list.
      CALL read_list( jinUnit,4,ndriveVarMax,'>VARS','>ENDVARS',' ','init_drive'  &
                     ,ndriveVarIn,cvar1=tmpName,cvar1Pos=1,ivar=tmpFlag,ivarPos=2  &
                     ,cvar2=tmpInterp,cvar2Pos=3,cvar3=tmpNameFile,cvar3Pos=4 )

    CASE ( formatNc )
      sdfFile = .TRUE.
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_drive',tagNc,preInit=.TRUE. )
      READ(jinUnit,*) ncTypeDrive
!     Check that type of files is recognised.
      CALL checkNcType( ncTypeDrive,'init_drive' )

!     Read details of variables. Read 4 values from a blank-delimited list.
      CALL read_list( jinUnit,4,ndriveVarMax,'>VARS','>ENDVARS',' ','init_drive'  &
                     ,ndriveVarIn,cvar1=tmpName,cvar1Pos=1,cvar2=tmpNameSDF,cvar2Pos=2  &
                     ,cvar3=tmpNameFIle,cvar3Pos=3,cvar4=tmpInterp,cvar4Pos=4 )

    CASE default
      WRITE(*,*)'ERROR: init_drive: no code for driveFormat=',TRIM(driveFormat)
      STOP

  END SELECT
!-------------------------------------------------------------------------------
! If necessary, read list of file names and times.
!-------------------------------------------------------------------------------
  IF ( readList ) THEN
!   One or more non-time-templated files. Read details from another file.
    IF ( echo ) WRITE(*,*)'Reading names of drive files from file.'
    iunit = fileUnit( formatAsc )   !  get unit
    CALL openFile( 1,.FALSE.,iunit,'read',formatAsc,listFile,'old' )
    READ(iunit,*)   !  header
    DO i=1,ndriveFileTime
      READ(iunit,*,iostat=ierr) driveFileName(i),driveFileDate(i),cdriveTime
      IF ( ierr /= 0 ) THEN
        IF ( ierr < 0 ) WRITE(*,*)'End of file before found ',ndriveFileTime,' names.'
        WRITE(*,*)'ERROR: init_drive: error reading list of file names.'
        STOP
      ENDIF
      driveFileTime(i) = chhmmss_to_s( cdriveTime,'init_drive' )
    ENDDO
    !DSM CALL closeFile( iunit,formatAsc )
    if (iunit/=jinUnit) CALL closeFile( iunit,formatAsc ) !DSM
  ELSE
!   Have already read time, now convert to integer.
    driveFileTime(1) = chhmmss_to_s( cdriveTime,'init_drive' )
  ENDIF  !  readList

!-------------------------------------------------------------------------------
! Finished reading input for this section.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Check values.
!-------------------------------------------------------------------------------
  IF ( driveFormat /= formatAsc ) noNewLineDrive=.FALSE.

  IF ( noNewLineDrive ) THEN
    IF ( ndriveHeaderField /= 0 ) THEN
      WRITE(*,*)'ERROR: init_drive: ndriveHeaderField must be zero for noNewLineDrive.'
      STOP
    ENDIF
    IF ( nxIn*nyIn>1 ) THEN
      WRITE(*,*)'ERROR: init_drive: driving data on one line is only allowed if'
      WRITE(*,*)'there is only one point in input file.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Make a series of checks that the user-specified driving data `timestep' is
! valid.
!-------------------------------------------------------------------------------
! For now, demand that drive data period is 1 day or less (and not a "special" period).
  IF ( driveDataPer<0 .OR. driveDataPer>iSecInDay ) THEN
    WRITE(*,*)'ERROR: init_drive'
    WRITE(*,*)'Period of driving data must be <= one day, and not a'
    WRITE(*,*)'"special" period (e.g. monthly).'
    STOP
  ENDIF

! Check that data period is a multiple of timestep.
  IF ( MOD( driveDataPer,NINT(timeStep) ) /= 0 ) THEN
     WRITE(*,*)'ERROR: init_drive: period of driving data must be a multiple of model timestep'
     WRITE(*,*) 'driveDataPer, NINT(timeStep) = ',driveDataPer, NINT(timeStep)
     STOP
  ENDIF

! Check that drive data are in phase with days (data refer to same times each day).
  IF ( MOD(iSecInDay,driveDataPer) /= 0 ) THEN
    WRITE(*,*)'ERROR: init_drive: drive data must be in phase with days.'
    STOP
  ENDIF

! Convert data period from seconds to timesteps.
  driveDataPer = driveDataPer / NINT(timeStep)

!-------------------------------------------------------------------------------
! Establish if file names contain templating, and are consistent in this.
!-------------------------------------------------------------------------------

  DO i=1,ndriveFileTime

    CALL check_template( driveDataPer,driveFilePer,NINT(timeStep),driveFileTime(i)  &
                        ,driveFileDate(i),.FALSE.,.FALSE.,driveFileName(i),'init_drive'  &
                        ,templateT,templateV,driveTemplateUnits )

    IF ( templateT ) driveTemplateT = .TRUE.
    IF ( templateV ) driveTemplateV = .TRUE.

!   Stop if the template does not match our expectations.

    IF ( driveTemplateT .NEQV. templateT ) THEN
      WRITE(*,*)'ERROR: init_drive: driveTemplateT .NEQV. templateT'
      WRITE(*,*)'File names are inconsistent.'
      IF ( templateT ) THEN
        WRITE(*,*)'File name includes time template strings, but previous&
                    & files did not'
      ELSE
        WRITE(*,*)'File name does not have time template strings, but&
                    & previous files did.'
      ENDIF
      WRITE(*,*)'i=',i,' driveFileName=',TRIM(driveFilename(i))
      STOP
    ENDIF

    IF ( driveTemplateV .NEQV. templateV ) THEN
      WRITE(*,*)'ERROR: init_drive: driveTemplateV .NEQV. templateV'
      WRITE(*,*)'File names are inconsistent.'
      IF ( templateV ) THEN
        WRITE(*,*)'File name includes variable name template strings, but previous&
                    & files did not'
      ELSE
        WRITE(*,*)'File name does not have variable name template strings, but&
                    & previous files did.'
      ENDIF
      WRITE(*,*)'i=',i,' driveFileName=',TRIM(driveFilename(i))
      STOP
    ENDIF

  ENDDO  !  files

!-------------------------------------------------------------------------------
! Insist that only one file name was given if time-templating is to be used.
!-------------------------------------------------------------------------------
  IF ( driveTemplateT .AND. ndriveFileTime>1 ) THEN
    WRITE(*,*)'ERROR: init_drive: driveTemplateT .AND. ndriveFileTime>1'
    WRITE(*,*)'File name includes time template strings, but ndriveFileTime>1.'
    WRITE(*,*)'To use time templating, set ndriveFileTime=1.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! If noNewLineDrive, all data for a given time are read via a single read  -
! cannot access more than one file.
!-------------------------------------------------------------------------------
  IF ( driveTemplateV .AND. noNewLineDrive ) THEN
    WRITE(*,*)'ERROR: init_drive: driveTemplateV .AND. noNewLineDrive'
    WRITE(*,*)'Driving data on one line cannot be used with'
    WRITE(*,*)'file variable name templating (driveTemplateV must be F).'
    WRITE(*,*)'i.e. all data for a given time must be in one file.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Extra checks for time templating.
!-------------------------------------------------------------------------------

  IF ( driveTemplateT ) THEN
!-------------------------------------------------------------------------------
!   Check that the file period is acceptable.
!   File period is only used for time-templated files.
!-------------------------------------------------------------------------------
    IF ( driveFilePer < 0 ) THEN
!     Check that the period is recognised.
!     Note that we don't check if the file period is a multiple of data period.
      SELECT CASE ( driveFilePer )
        CASE ( periodAnn, periodMon )
!         OK, nothing to do.
        CASE ( periodOneFile )
!         All data in one file - this makes no sense for time templating.
          WRITE(*,*)'ERROR: init_drive: driveFilePer=',driveFilePer
          WRITE(*,*)'All data in one file does  not make sense for time templating.'
          STOP
        CASE default
          WRITE(*,*)'ERROR: init_drive: do not recognise driveFilePer=',driveFilePer
          STOP
      END SELECT
    ELSE
!     driveFilePer >= 0.
!     Check file period is a multiple of data period (which is already known to be
!     a multiple of timestep length, and to be a factor of one day).
      IF ( MOD( driveFilePer,driveDataPer*NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: init_drive: file period must be a multiple of data period.'
        STOP
      ENDIF
!     Further restriction for time template - insist that a single file cannot
!     hold data for more than 24 hours - makes life easier.
      IF ( driveFilePer > iSecInDay ) THEN
        WRITE(*,*)'ERROR: init_drive: files too infrequent.'
        WRITE(*,*)'For these template files, a single file cannot &
                   &hold data for more than 24 hours.'
        WRITE(*,*)'This restriction does not apply to files with a &
               &special" period (e.g. monthly).'
        STOP
      ENDIF
!     Convert file period to number of timesteps.
      driveFilePer = driveFilePer / NINT( timeStep )
    ENDIF  !  driveFilePer

!   Save some values.
    driveFileTemplate = driveFileName(1)
    driveTemplateDate = driveFileDate(1)
    driveTemplateTime = driveFileTime(1)

  ENDIF  !  driveTemplateT

!-------------------------------------------------------------------------------
! Check that file name or template, is (vaguely) acceptable.
! Note that this means that a file name for use without templating cannot
! contain the character "%".
!-------------------------------------------------------------------------------
  DO i=1,ndriveFileTime
    CALL check_template( driveDataPer,driveFilePer,NINT(timeStep),driveFileTime(i)  &
                        ,driveFileDate(i),.FALSE.,.TRUE.,driveFileName(i),'init_drive'  &
                        ,templateT,templateV,driveTemplateUnits )
  ENDDO

!-------------------------------------------------------------------------------
  IF ( .NOT. driveTemplateT ) THEN
!-------------------------------------------------------------------------------
!   Set the time/date for the last file (one more than we actually have) to be 100
!   years after start of penultimate file - easiest option for now (just don't
!   have hundreds of years of data....).
!-------------------------------------------------------------------------------
    i = ndriveFileTime+1
    driveFileName(i) = 'dummy_file'
    CALL timeDate( driveFileTime(i-1),driveFileDate(i-1),1200,'mon',l_360  &
                  ,driveFileTime(i),driveFileDate(i),'init_drive')
!-------------------------------------------------------------------------------
!   Check that the files are listed in chronological order.
!-------------------------------------------------------------------------------
    DO i=2,ndriveFileTime
      t = i - 1
      IF ( .NOT. timeDate_cmp( driveFileTime(i),driveFileDate(i),'>'  &
                      ,driveFileTime(t),driveFileDate(t),'init_drive' )  ) THEN
        WRITE(*,*)'ERROR: init_drive: files are not listed in chronological order.'
        WRITE(*,*)'Problem files and times:'
        WRITE(*,*) t,TRIM(driveFileName(t)),driveFileDate(t),driveFileTime(t)
        WRITE(*,*) i,TRIM(driveFileName(i)),driveFileDate(i),driveFileTime(i)
        STOP
      ENDIF
    ENDDO

  ENDIF  !  driveTemplateT
!-------------------------------------------------------------------------------

! Process the variables.
  CALL init_drive_var( ioWindSpeed,sdfFile,useDiffRad  &
                      ,tmpFlag,tmpInterp,tmpName,tmpNameSDF  &
                      ,tmpNameFile,driveVarName,driveVarFileName,varDesc )

! Count number of variables that are required.
  ndriveVarNeed = COUNT( driveVarUse(:) )

!-------------------------------------------------------------------------------
! If time interpolation is indicated, but data period equals timestep, change
! interpolation flag to indicate that no interpolation is required.
!-------------------------------------------------------------------------------
  IF ( driveDataPer == 1 ) THEN
    DO ivar=1,ndriveVarIn
      ipos = driveVarInSort(ivar)
      SELECT CASE ( driveVarInterp(ipos) )
        CASE ( 'i' )
          WRITE(*,*)'driveDataPer=timestep, so changing flag i to nf'
          driveVarInterp(ipos) = 'nf'
        CASE ( 'b' )
          WRITE(*,*)'driveDataPer=timestep, so changing flag b to nb'
          driveVarInterp(ipos) = 'nb'
        CASE ( 'c' )
          WRITE(*,*)'driveDataPer=timestep, so changing flag c to nc'
          driveVarInterp(ipos) = 'nc'
        CASE ( 'f' )
          WRITE(*,*)'driveDataPer=timestep, so changing flag f to nf'
          driveVarInterp(ipos) = 'nf'
      END SELECT
    ENDDO  !  ivar
  ENDIF

!-------------------------------------------------------------------------------
! Check that interpolation flag is recognised. Also look for time-averaged
! input, and work out how many and what time levels of input data are required.
! e.g. needTime(0:1)=T indicates that, when model is running, at a drive data
! timestep, we need the current and next drive data to be able to calculate /
! interpolate values for all model timesteps before the next drive read.
!-------------------------------------------------------------------------------
  haveAve = .FALSE.
  needTime(:) = .FALSE.
  DO ivar=1,ndriveVarIn
    ipos = driveVarInSort(ivar)
    SELECT CASE ( driveVarInterp(ipos) )
      CASE ( 'b' )
        needTime(0:2) = .TRUE.
        haveAve = .TRUE.
      CASE ( 'c' )
        needTime(-1:2) = .TRUE.
        haveAve = .TRUE.
      CASE ( 'f' )
        needTime(-1:1) = .TRUE.
        haveAve = .TRUE.
      CASE ( 'i' )
        needTime(0:1) = .TRUE.
      CASE ( 'nb' )
        needTime(1) = .TRUE.
      CASE ( 'nc' )
        needTime(0:1) = .TRUE.
      CASE ( 'nf' )
        needTime(0) = .TRUE.
      CASE default
        WRITE(*,*)'ERROR: init_drive: do not recognise driveVarInterp: ',TRIM(driveVarInterp(ipos))
        STOP
    END SELECT

!-------------------------------------------------------------------------------
!   If weather generator is used, check that flag indicates forward or backward averages.
!   The flag is not used to indicate interpolation, just timing.
!-------------------------------------------------------------------------------
    IF ( useWGen ) THEN
      SELECT CASE ( driveVarInterp(ipos) )
        CASE ( 'nb', 'nf' )
!         These are acceptable.
        CASE default
          WRITE(*,*)'ERROR: flag not acceptable for use with weather generator.'
          WRITE(*,*)'Flag=',TRIM(driveVarInterp(ipos))
          WRITE(*,*)'Only acceptable flags for use with weather generator are nb and nf.'
        STOP
      END SELECT
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Check that timestep allows interpolation of averages to be conservative.
!-------------------------------------------------------------------------------
  IF ( haveAve .AND. MOD(driveDataPer,2)/=0 ) THEN
    WRITE(*,*)'ERROR: init_drive: period of drive data is not a multiple of 2*timestep'
    WRITE(*,*)'The algorithm for interplation of time-averaged input data requires'
    WRITE(*,*)'that the data period comprises an even number of model timesteps.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Work out how many time levels of data to hold.  If a given time is required
! for one variable, all variables are held at that time.
!-------------------------------------------------------------------------------
  driveTimeIndex(:) = -9
! Get earliest time required.
  DO t=-1,2
    IF ( needTime(t) ) THEN
      driveTimeIndex(1) = t
      EXIT
    ENDIF
  ENDDO
! Get latest time required.
  DO t=2,-1,-1
    IF ( needTime(t) ) THEN
      driveTimeIndex(2) = t
      EXIT
    ENDIF
  ENDDO
  ndriveDataTime = driveTimeIndex(2) - driveTimeIndex(1) + 1
  if (echo) print*,'driveTimeIndex=',driveTimeIndex(:)

!-------------------------------------------------------------------------------
! Allocate space for input data.
!-------------------------------------------------------------------------------
  IF ( echo ) WRITE(*,*)'Allocating for drive data. ndriveVarIn=',ndriveVarIn  &
       ,' ndriveVar=',ndriveVar,' driveTimeIndex=',driveTimeIndex(:)
  CALL allocate_arrays( 'init_drive 3' )

!-------------------------------------------------------------------------------
! Get file unit(s). Open a scratch file on each to reserve it. Not needed for SDFs.
! In fact we could probably proceed without "reserving" units in this
! way...but this seems to be working for now.
!-------------------------------------------------------------------------------
! Initialise driveUnit. This value should NOT be a possible netCDF ID, so use -1.
  driveUnit(:) = -1

  IF ( driveFormat /= formatNc ) THEN
    DO i=1,ndriveUnit
      driveUnit(i) = fileUnit( driveFormat )
      driveFileOnUnit(i) = 'scratch_file'
      CALL openFile( 1,.FALSE.,driveUnit(i),'write','scratch',driveFileOnUnit(i),'new' )
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Set the reference level heights for wind and air temperature/humidity to
! user-specified constant values.
!-------------------------------------------------------------------------------
  IF ( .NOT. routeOnly ) THEN
    z1_uv(:,:) = z1_uv_val
    z1_tq(:,:) = z1_tq_val
  ENDIF

!--------------------------------------------------------------------------------
! If requested, report to stdout the results of reading INIT_DRIVE options.
!--------------------------------------------------------------------------------
  IF ( echo ) THEN

!   General options.
    IF ( .NOT. routeOnly ) THEN

      IF ( useWGen ) THEN
        WRITE(*,*)'A weather generator will be used to disaggregate daily input values.'
        WRITE(*,*)'NB Precipitation should be the daily amount (kg m-2)'
      ENDIF

      SELECT CASE ( ioPrecipType )
        CASE ( 1,5 )
          WRITE(*,*) 'Total precipitation is input.'
          WRITE(*,*) '  Convective snow is assumed to be zero.'
        CASE ( 2 )
          WRITE(*,*) 'Rainfall and snowfall are input separately.'
          WRITE(*,*) '  Convective snow is assumed to be zero.'
        CASE ( 3 )
          WRITE(*,*) 'Three components of precipitation are input:'
          WRITE(*,*) '  Large-scale rainfall'
          WRITE(*,*) '  Convective rainfall'
          WRITE(*,*) '  Large-scale (assumed to equal total) snowfall'
          WRITE(*,*) '  Convective snow is assumed to be zero.'
        CASE ( 4 )
          WRITE(*,*) 'Four components of precipitation are input:'
          WRITE(*,*) ' Large-scale and convective rainfall.'
          WRITE(*,*) ' Large-scale and convective snowfall.'
      END SELECT

      SELECT CASE ( ioPrecipType )
        CASE ( 1,5 )
          WRITE(*,*) 'Snowfall when temperature <=',tForSnow,' (K)'
          IF ( .NOT. l_point_data ) WRITE(*,*)  &
             'Convective rainfall when temperature >=',tForCRain,' (K)'
      END SELECT

      IF ( .NOT. l_point_data ) WRITE(*,*)  &
            'Fractional coverage for convective precipitation=',conFrac

!      IF ( useWgen ) THEN
!        WRITE(*,*)'Duration of convective rainfall events (hr)=',durConvRain
!        WRITE(*,*)'Duration of large-scale rainfall events (hr)=',durLSRain
!        WRITE(*,*)'Duration of large-scale snowfall events (hr)=',durLSSnow
!        WRITE(*,*)'Amplitude of diurnal temperature variation (K)=',tAmplConst
!      ENDIF

      SELECT CASE ( io_rad_type )
        CASE ( 1 )
          WRITE(*,*) 'Downward fluxes of shortwave and longwave radiation are input.'
        CASE ( 2 )
          WRITE(*,*) 'The downward shortwave and net (all-wavelength) downward radiation fluxes are input.'
        CASE ( 3 )
          WRITE(*,*) 'Net downward fluxes of shortwave and longwave radiation are input.'
      END SELECT

      IF ( useDiffRad ) THEN
        WRITE(*,*)'Diffuse radiation is input.'
      ELSE
        WRITE(*,*)'A constant fraction of diffuse radiation is assumed: ',diffFracConst
      ENDIF

      IF ( ioWindSpeed ) THEN
        WRITE(*,*)'The wind speed is input.'
      ELSE
        WRITE(*,*)'The two components of the horizontal wind (u,v) are input.'
      ENDIF

    ELSE  !  routeOnly

      WRITE(*,*)'Surface and subsurface components of runoff are input.'

    ENDIF  !  routeOnly

    WRITE(*,*) ndriveVarNeed,' driving variables are required.'
    WRITE(*,*) ndriveVarIn,' variables are read in.'

    SELECT CASE ( driveFormat )
      CASE ( formatAsc )
        WRITE(*,*)'Driving data are stored in ASCII files.'
      CASE ( formatBin )
        WRITE(*,*)'Driving data are stored in binary files.'
        IF ( byteSwapDrive ) THEN
          WRITE(*,"(50('#'))")
          WRITE(*,*) 'WARNING: all driving (meteorological) data will be byteswapped.'
          WRITE(*,*) 'This means that the byte order ("endianness") of the data'
          WRITE(*,*) 'will be be reversed. This is necessary when transferring files'
          WRITE(*,*) 'between different types of computers.'
          WRITE(*,*) 'However, if this flag is incorrectly set to TRUE, the result'
          WRITE(*,*) 'will be garbled data!'
          WRITE(*,"(50('#'))")
        ENDIF
      CASE ( formatNc )
        WRITE(*,*)'Driving data are stored in netCDF files.'
      CASE ( formatPP )
        WRITE(*,*)'Driving data are stored in PP format files.'
    END SELECT

    WRITE(*,*) ndriveUnit,' driving data files are open at any one time'

    IF ( driveFormat /= formatNc ) WRITE(*,*)'There are ',nfieldDriveFile,' variables in each file.'

    IF ( driveTemplateT ) THEN
      WRITE(*,*)'Files are named following time-templating conventions.'
      IF ( driveEndTime ) THEN
        WRITE(*,*)'Convention uses end time of file.'
      ELSE
        WRITE(*,*)'Convention uses start time of file.'
      ENDIF
    ENDIF

!--------------------------------------------------------------------------------
!   Report to stdout more detailed information for each variable.
!--------------------------------------------------------------------------------
    CALL init_drive_var_detail( driveVarname,varDesc )

  ENDIF   !  echo

!-------------------------------------------------------------------------------
! Initialise driving data. Read in any data values that should have already been
! read by the time of the first timestep (i.e. if it wasn't the start of the run,
! they would have been read by this time).
! The 2nd dummy argument is dataStepMax, which at present is only really needed
! for vegetation data (i.e. not here), but making it an optional argument has
! its disadvantages - so tmpDataStepMax is passed here (but not expected to change!).
! Admittedly this is rather untidy.....
!
! This is unnecessary if IMOGEN is being used
!-------------------------------------------------------------------------------
   tmpDataStepMax = driveDataPer
   IF( .NOT. l_imogen ) THEN
     CALL data_init( driveDataPer,driveFilePer,driveTemplateDate,driveTemplateTime  &
      ,tmpDataStepMax,driveDataStep,driveDataStepInit  &
      ,driveDate,driveDateInit  &
      ,driveFile,driveFileStep,driveResetStep,driveResetStepPrev  &
      ,driveTime,driveTimeInit  &
      ,driveFileDate,driveFileTime,driveFileTemplate,driveTemplateUnits,driveFileName  &
      ,driveTimeIndex,.FALSE.,driveEndTime,templateT,notNextDrive,'drive' )
   ENDIF

!   print*,'+data_init driveDateInit,driveTimeInit,driveDataStepInit=',driveDateInit,driveTimeInit,driveDataStepInit

!-------------------------------------------------------------------------------
! Deal with the possibility that driving data may need to be "recycled" to cope
! with spin up.
! Arguments: zero=timestep.
!-------------------------------------------------------------------------------
  IF ( spinUp ) CALL calc_reset_step( 0,driveResetStep,driveResetStepPrev,'drive' )

END SUBROUTINE init_drive






!################################################################################
!################################################################################
! subroutine init_drive_var
! Subroutine to establish what driving data are provided, and what are required.
!################################################################################
!################################################################################
SUBROUTINE init_drive_var( ioWindSpeed,sdfFile,useDiffRad  &
                          ,tmpFlag,tmpInterp,tmpName,tmpNameSDF  &
                          ,tmpNameFile,driveVarName,driveVarFileName,varDesc )

!-------------------------------------------------------------------------------
  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE drive_io_vars, ONLY :  &
!   imported scalar parameters
      ndriveExtraMax  &
!  imported scalars with intent(in)
     ,driveFormat,driveTemplateV  &
     ,ndriveExtra,ndriveVarIn,ndriveVarMax  &
     ,nfieldDriveFile,ioPrecipType,io_rad_type,useWgen  &
!  imported scalars with intent(out)
     ,iposDiffRad,iposExtra,iposLWD,iposLWN,iposPrecip,iposPrecipCR  &
     ,iposPrecipCS,iposPrecipLR,iposPrecipLS  &
     ,iposPrecipTR,iposPrecipTS,iposPstar,iposQ,iposRN,iposSubSurfRoff   &
     ,iposSurfRoff,iposSWD,iposSWN,iposT,iposU,iposV,iposWind  &
     ,iposOzone  &
     ,ndriveUnit,ndriveVar   &
!  imported arrays with intent(out)
     ,driveVarInterp,driveUnitUse,driveVarNameSDF,driveVarNameUnit,driveVarInSort  &
     ,driveVarFlag,driveVarPos,driveVarPosIn,driveVarStash,driveVarUse

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatBin,formatNc,formatPP

  USE misc_utils, ONLY :  &
!  imported procedures
     checkVarPos,repeatVal,varList

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     routeOnly,l_o3_damage
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    ioWindSpeed           &!  T means windspeed is input, F means 2 components of wind are input
   ,sdfFile               &!  T means data will come from an SDF (self-describing file)
   ,useDiffRad             !  T means diffuse radiation is input
!                             F means diffuse radiation is set to a constant fraction

  CHARACTER(len=*), INTENT(in) ::  &!  in ARRAYS
    tmpInterp(ndriveVarMax)    &!  values of driveVarInterp, as read in
   ,tmpName(ndriveVarMax)      &!  names of forcing variables, as read in
   ,tmpNameSDF(ndriveVarMax)   &!  names of forcing variables for SDF files
   ,tmpNameFile(ndriveVarMax)   !  names used in file names, as read in

  INTEGER, INTENT(inout) ::  &!  inout arrays
    tmpFlag(ndriveVarMax)  !  IN: values of driveVarFlag as read in
!                             OUT - not needed, but space reused in this routine.

  CHARACTER(len=*), INTENT(out) ::  &!  out arrays
    driveVarName(ndriveVarMax)      &!  names of all possible forcing variables
   ,driveVarFileName(ndriveVarMax)  &!  names of driving variables as used in file names
   ,varDesc(ndriveVarMax)            !  description of each variable

  INTEGER ::  &!  local SCALARS
    dinpos               &!  indicates next available location in driveDataIn array
   ,dpos                 &!  indicates next available location in driveData array
   ,i,ipos,ivar,jpos,jvar !  work

  INTEGER :: &!  local arrays
    ival(1)                 !  work

  LOGICAL ::  &!  local SCALARS
    checkNames,checkPos,errFound   !  work

  LOGICAL ::  &! local ARRAYS
    done(ndriveVarMax)         &!  work
   ,foundVar(ndriveVarMax)      !  work - currently redundant, but an argument
!--------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Setup list of all possible variables. There are ndriveVarMax of these.  These
! are essentially parameters, but the list is easier to maintain if setup as
! executable statements here.
! There are other possibilities for STASH code for some variables -
!   e.g. STASH 10 is q after timestep.
!-------------------------------------------------------------------------------
  i = 0
  i=i+1; iposDiffRad=i; driveVarName(i)='diff_rad'
  driveVarStash(i) = 0  !  not known
  varDesc(i) = 'Diffuse radiation'

  i=i+1; iposLWD=i; driveVarName(i)='lw_down'
  driveVarStash(i) = 2207  !  PP code 205
  varDesc(i) = 'Downward longwave radiation'

  i=i+1; iposLWN=i; driveVarName(i)='lw_net'
  driveVarStash(i) = 2201   !  PP code 187
  varDesc(i) = 'Net downward longwave radiation'

  i=i+1; iposPrecip=i; driveVarName(i)='precip'
  driveVarStash(i) = 5216   !  PP code 90
  varDesc(i) = 'Total precipitation rate'

  i=i+1; iposPrecipCR=i; driveVarName(i)='precipCR'
  driveVarStash(i) = 5205   !  PP code 98
  varDesc(i) = 'Convective rainfall rate'

  i=i+1; iposPrecipCS=i; driveVarName(i)='precipCS'
  driveVarStash(i) = 5206   !  PP code 117
  varDesc(i) = 'Convective snowfall rate'

  i=i+1; iposPrecipLR=i; driveVarName(i)='precipLR'
  driveVarStash(i) = 4203   !  PP code 99
  varDesc(i) = 'Large-scale rainfall rate'

  i=i+1; iposPrecipLS=i; driveVarName(i)='precipLS'
  driveVarStash(i) = 4204   !  PP code 118
  varDesc(i) = 'Large-scale snowfall rate'

  i=i+1; iposPrecipTR=i; driveVarName(i)='precipTR'
  driveVarStash(i) = 5214   !  PP code 97
  varDesc(i) = 'Total rainfall rate'

  i=i+1; iposPrecipTS=i; driveVarName(i)='precipTS'
  driveVarStash(i) = 5215
  varDesc(i) = 'Total snowfall rate'

  i=i+1; iposPstar=i; driveVarName(i)='pstar'
  driveVarStash(i) = 1
  varDesc(i) = 'Surface pressure'

  i=i+1; iposQ=i; driveVarName(i)='q'
  driveVarStash(i) = 3237   !  STASH 3237 is q at 1.5m
  varDesc(i) = 'Specific humidity'

  i=i+1; iposRN=i; driveVarName(i)='rad_net'
  driveVarStash(i) = 0   !  STASH code not known
  varDesc(i) = 'Net downward all-wavelength radiation.'

  i=i+1; iposSubSurfRoff=i; driveVarName(i)='subSurfRoff'
  driveVarStash(i) = 8235
  varDesc(i) = 'Subsurface runoff rate'

  i=i+1; iposSurfRoff=i; driveVarName(i)='surfRoff'
  driveVarStash(i) = 8234
  varDesc(i) = 'Surface runoff rate'

  i=i+1; iposSWD=i; driveVarName(i)='sw_down'
  driveVarStash(i) = 1235  !  PP code 200
  varDesc(i) = 'Downward shortwave radiation'

  i=i+1; iposSWN=i; driveVarName(i)='sw_net'
  driveVarStash(i) = 1201  !  PP code 186
  varDesc(i) = 'Net downward shortwave radiation'

  i=i+1; iposT=i; driveVarName(i)='t'
  driveVarStash(i) = 3236  !  STASH 3236 is T at 1.5m
  varDesc(i) = 'Air temperature'

  i=i+1; iposU=i; driveVarName(i)='u'
  driveVarStash(i) = 3225   !  STASH 3225 is u at 10m
  varDesc(i) = 'Zonal wind speed'

  i=i+1; iposV=i; driveVarName(i)='v'
  driveVarStash(i) = 3226  !  STASH 3226 is v at 10m
  varDesc(i) = 'Meridional wind speed'

  i=i+1; iposWind=i; driveVarName(i)='wind'
  driveVarStash(i) = 3249  !  STASH 3249 is wind at 10m
  varDesc(i) = 'Wind speed'

  i=i+1; iposOzone=i; driveVarName(i)='ozone'
  driveVarStash(i) = 0  !  STASH not known
  varDesc(i) = 'Surface ozone concentration'

! The possible "extra" variables.
  DO ivar=1,ndriveExtraMax
    i = i + 1
    IF ( ivar == 1 ) iposExtra = i
    WRITE(driveVarName(i),"(a5,i2.2)") 'extra',ivar
    driveVarStash(i) = 0  !  not known
    WRITE(varDesc(i),"(a16,i2)") 'extra variable #',ivar
  ENDDO

!-------------------------------------------------------------------------------
! Set flags indicating if a variable is needed for the selected configuration,
! and if it must be read in or can be derived from other variables.

! We differentiate between:
! a) the model configuration: i.e. what "physics" or "processes" are active
!      This determines what driving variables are required (ultimately -
!      they may be derived).
! b) the input/output (i/o) configuration: i.e. what input variables are input

! The model configuration is completely described by other variables (i.e. various
! flags) - it is not described by the choice of input variables. Rather, the
! necessary variables follow from the model configuration. For a given model
! configuration, there may be more than one possible i/o configurations, i.e.
! there may be more than one set of input variables that can provide the required
! input to the process descriptions.
!
! e.g. the model configuration determines if we need liquid precipitation.
! If we do need liquid precip, we can get this in more than one way (more than
! one i/o configuration) - e.g. liquid precip could be read in directly, or it
! could be calculated from total precip.
!
! At the end of this subroutine, variables will fall into one of 6 states:
! 1) driveVarUse=T, driveVarPosIn>0, driveVarPos>0
!      Variable is required by model processes, will be read in, and is stored
!      in driveData (temporally interpolated/disaggregated data).
!      Example: air temperature
! 2) driveVarUse=T, driveVarPosIn=0, driveVarPos>0.
!      Variable is required by model proceses, but is not read in - it is
!      calculated from others and then stored in driveData (temporally
!      interpolated/disaggregated data)
!      Example: Air temperature if using a "weather generator". The diurnal
!      temperature range can be read in, then the diurnal course of air
!      temperature calculated by the generator.
! 3) driveVarUse=T, driveVarPosIn=0, driveVarPos=0
!      Variable is required, but is not read in, nor is it stored in driveData
!      (temporally interpolated/disaggregated data)- rather, it is
!      calculated from others and loaded directly (not via driveDataIn).
!      Examples: 1) solid precipitation can be calculated from total precipitation
!      and a temperature threshold.  2) downward radiation can be calculated from
!      net radiation.
! 4) driveVarUse=F, driveVarPosIn>0, driveVarPos>=0
!      Variable is not directly required by model proceses, but will be read in.
!      It is not stored in driveData (temporally interpolated/disaggregated data).
!      This variable is used to derive others, but is not itself needed by the model.
!      Examples: (1) When using a weather generator, the diurnal temperature range
!      is read in, but is is not stored in driveData - the disaggregated temperature
!      is stored instead. (2) total liquid precip is read in, and is stored in
!      driveDataIn and driveData.
! 5) driveVarUse=F, driveVarPosIn=0, driveVarPos>0
!      Variable is not directly required by model proceses, and is not read in,
!      but is stored in driveData (temporally interpolated/disaggregated data).
!      This variable is used to derive others, but is not itself needed by the model.
!      Example: When using a weather generator, the total solid precipitation is not read
!      in, nor needed itself, but it is calculated from total precip and held disaggregated
!      in time (before being loaded into large-scale solid precip).
! 6) driveVarUse=F, driveVarInPos=0, driveVarPos=0
!      Variable is not used in any way.
!
! Much of this is currently overkill, and may look unneccesarily complicated (maybe
! it is overcomplicated!), but it is an attempt to provide moderate "future proofing"
! for situations that may arise in the "medium term"!?! In the short term: my
! apologies.
! My thinking: by holding all possible variables in arrays dim(ndriveVarMax),
! a "simple" set of comparisons can identify if and how a variable is being used.
! These arrays have to be set up correctly, but once set up, can be used elsewhere
! - maintenance is easier than if we have several IFs, at various points in code,
! to establish how a variable is to be used.

! At present we have 2 possible main model configurations:
! 1) NOT routeOnly, needing 11 variables:
!      sw_down,lw_down,T,q,u,v,pstar,ls_rain,ls_snow,con_rain,con_snow
!    These can be provided by various i/o configurations, and possibly using
!    a weather generator.
!    Diffuse radiation (diff_rad) is a bit different - we actually only need
!    the diffuse fraction (diff_frac), which can be got via 2 differnet routes.
! 2) routeOnly, needing 2 variables:
!      surf_roff, sub_surf_roff
!-------------------------------------------------------------------------------
! Initialise.
  driveVarUse(:) = .FALSE.     !  variable not needed
  driveVarFlag(:) = 0          !  position in file = 0
!   Note that driveVarFlag/=0 is used to test if a variable has been found in
!   the input list.
  driveVarPos(:) = 0           !  position in driveData
  driveVarPosIn(:) = 0         !  position in driveDataIn
!   Note that driveVarPosIn>0 is used to test if a variable is to be input
  driveVarFileName(:) = '-'
  driveVarNameSDF(:) = '-'

! Initialise counters. These are used to indicate the next available
! locations in the driveDataIn driveData arrays.
  dinpos = 0
  dpos = 0

  IF ( .NOT. routeOnly ) THEN

!   Variables that are always required, and must be read in (can't be derived).
!   These variables are held in both the driveDataIn and driveData arrays.

!-------------------------------------------------------------------------------
!   Surface pressure.
!   This is always required and is read in unless a weather generator is used.
!-------------------------------------------------------------------------------
    ipos = iposPstar
    driveVarUse(ipos) = .TRUE.
    dpos=dpos+1; driveVarPos(ipos) = dpos           !  next location in driveData
    IF ( .NOT. useWgen ) THEN
      dinpos=dinpos+1
      driveVarPosIn(ipos) = dinpos   !  next location in driveDataIn
    ENDIF

!   Specific humidity of air.
    ipos = iposQ
    driveVarUse(ipos) = .TRUE.
    dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
    dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData

!   Air temperature.
    ipos = iposT
    driveVarUse(ipos) = .TRUE.
    dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos   !  next location in driveDataIn
    dpos=dpos+1; driveVarPos(ipos) = dpos           !  next location in driveData

!-------------------------------------------------------------------------------
!   Radiation.
!   Currently we always require downward fluxes, but these can be
!   derived from net fluxes. This case is rather unusual in that the
!   net fluxes are stored in the variables used for downward fluxes, until
!   the downward fluxes are calculated in subroutine CONTROL.
!   So although we require the downward fluxes, we may not save them
!   at this level of code - waits until CONTROL.
!-------------------------------------------------------------------------------

!   Set driveVarUse for variables that are required.
    driveVarUse(iposLWD) = .TRUE.
    driveVarUse(iposSWD) = .TRUE.

!   Use i/o flag to work out what variables should be input.
    SELECT CASE ( io_rad_type )
      CASE ( 1 )
!       Downward fluxes are input (and so saved in driveDataIn and driveData).
!       Longwave.
        ipos = iposLWD
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Shortwave.
        ipos = iposSWD
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos        !  next location in driveData
      CASE ( 2 )
!       Net (all-wavelength) downward flux and downward shortwave fluxes are input.
!       This is saved in driveData before being loaded into downward longwave flux.
        ipos = iposRN
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Shortwave.
        ipos = iposSWD
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos        !  next location in driveData
      CASE ( 3 )
!       Net downward fluxes are input (and so saved in driveDataIn).
!       These are saved in driveData before being loaded into downward flux.
!       Longwave.
        ipos = iposLWN
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos        !  next location in driveData
!       Shortwave.
        ipos = iposSWN
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE default
        WRITE(*,*)'ERROR: init_drive: no code for io_rad_type=',io_rad_type
        STOP
    END SELECT

!-------------------------------------------------------------------------------
!   Diffuse radiation.
!-------------------------------------------------------------------------------
    IF ( useDiffRad ) THEN
      ipos = iposDiffRad
      driveVarUse(ipos) = .TRUE.
      dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
      dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
    ENDIF


!-------------------------------------------------------------------------------
!   Precipitation.
!   Currently we always require four components, but these can be
!   derived from totals. None of the totals is needed (except to derive components).
!-------------------------------------------------------------------------------
!   Set driveVarUse for variables that are required.
    driveVarUse(iposPrecipCR) = .TRUE.
    driveVarUse(iposPrecipCS) = .TRUE.
    driveVarUse(iposPrecipLR) = .TRUE.
    driveVarUse(iposPrecipLS) = .TRUE.

!   Use i/o flag to work out what variables should be input.
    SELECT CASE ( ioPrecipType )
      CASE ( 1 )
!       Total precip is input and held disaggregated.
        ipos = iposPrecip
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE ( 2 )
!       Total rainfall and total snowfall are input separately, and held disaggregated.
!       Total rainfall.
        ipos = iposPrecipTR
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Total snowfall
        ipos = iposPrecipTS
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE ( 3 )
!       Three components of precip are input and held disaggregated.
!       Convective rainfall.
        ipos = iposPrecipCR
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Large-scale rainfall.
        ipos = iposPrecipLR
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Total snowfall.
        ipos = iposPrecipTS
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE ( 4 )
!       Four components of precip are input and held disaggregated.
!       Convective rainfall.
        ipos = iposPrecipCR
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Large-scale rainfall.
        ipos = iposPrecipLR
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Convective snowfall.
        ipos = iposPrecipCS
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Large-scale snowfall.
        ipos = iposPrecipLS
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE ( 5 )
!       The total precipitation is input. Three components are held disaggregated.
!       Total precip.
        ipos = iposPrecip
        dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
!       Convective rainfall.
        ipos = iposPrecipCR
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Large-scale rainfall.
        ipos = iposPrecipLR
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!       Total snowfall.
        ipos = iposPrecipTS
        dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
      CASE default
        WRITE(*,*)'ERROR: init_drive: no code for ioPrecipType=',ioPrecipType
        STOP
    END SELECT

!-------------------------------------------------------------------------------
!   Wind.
!   Currently we always require two components, but these can be
!   set from the total. The total itself is not needed.
!-------------------------------------------------------------------------------
!   Set driveVarUse for variables that are required.
    driveVarUse(iposU) = .TRUE.
    driveVarUse(iposV) = .TRUE.

!   Use i/o flag to work out what variables should be input.
    IF ( ioWindSpeed ) THEN
!     (Total) wind speed is input.
      ipos = iposWind
      dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
      dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
    ELSE  !  NOT ioWindspeed
!     Two components are input.
!     Zonal component.
      ipos = iposU
      dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
      dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
!     Meridional component.
      ipos = iposV
      dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
      dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
    ENDIF !  ioWindSpeed

!-------------------------------------------------------------------------------
!   Ozone concentration
!-------------------------------------------------------------------------------
    IF ( l_o3_damage ) THEN
      ipos = iposOzone
      driveVarUse(ipos) = .TRUE.
      dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
      dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
    ENDIF

!-------------------------------------------------------------------------------

  ELSE  !  routeOnly

!   Currently we require two variables: surface and subsurface runoff.

!   Surface runoff.
    ipos = iposSurfRoff
    driveVarUse(ipos) = .TRUE.
    dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos   !  next location in driveDataIn
    dpos=dpos+1; driveVarPos(ipos) = dpos           !  next location in driveData

!   Subsurface runoff.
    ipos = iposSubSurfRoff
    driveVarUse(ipos) = .TRUE.
    dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos   !  next location in driveDataIn
    dpos=dpos+1; driveVarPos(ipos) = dpos           !  next location in driveData

  ENDIF   !  routeOnly

!-------------------------------------------------------------------------------
! Set up the required number of "extra" variables.
!-------------------------------------------------------------------------------
  DO ivar=1,ndriveExtra
    ipos = iposExtra + ivar - 1
    driveVarUse(ipos) = .TRUE.
    dinpos=dinpos+1; driveVarPosIn(ipos) = dinpos  !  next location in driveDataIn
    dpos=dpos+1; driveVarPos(ipos) = dpos          !  next location in driveData
  ENDDO

!-------------------------------------------------------------------------------
! Check that there are no repeated names in the list of names read in.
!-------------------------------------------------------------------------------
  IF ( repeatVal( tmpName(1:ndriveVarIn) ) ) THEN
    WRITE(*,*)'ERROR: init_drive_var: repeated variable name.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Establish what variables have been provided for input, and get related
! information into main arrays.
!-------------------------------------------------------------------------------
  CALL varList( ndriveVarIn,sdfFile,tmpName,driveVarName  &
               ,errFound,foundVar  &
               ,varFlagIn=tmpFlag,varFileNameIn=tmpNameFile  &
               ,varInterpIn=tmpInterp,varNameSDFin=tmpNameSDF  &
               ,varFlag=driveVarFlag,varFileName=driveVarFileName  &
               ,varInterp=driveVarInterp,varNameSDF=driveVarNameSDF )

  IF ( errFound ) THEN
    WRITE(*,*)'ERROR: init_drive_var: error from varList'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Check that the required variables are supplied.  The 6 possible states are
! listed above.   This is a check for bugs in the code above, and that variables
! are correctly specified in the run control file.
!-------------------------------------------------------------------------------
  errFound = .FALSE.

  DO ivar=1,ndriveVarMax

!-------------------------------------------------------------------------------
!   First we check for coding errors.
!-------------------------------------------------------------------------------
    IF ( driveVarUse(ivar) ) THEN
!     Either driveVarPosIn>0, driveVarPos>0
!     OR     driveVarPosIn=0, driveVarPos>0
!     OR     driveVarPosIn=0, driveVarPos=0
!     Check if we have a valid state.
      IF ( .NOT. (driveVarPosIn(ivar)>0 .AND. driveVarPos(ivar)>0 ) .AND.  &
           .NOT. (driveVarPosIn(ivar)==0 .AND. driveVarPos(ivar)>0 ) .AND.  &
           .NOT. (driveVarPosIn(ivar)==0 .AND. driveVarPos(ivar)==0 ) ) THEN
        errFound = .TRUE.
        WRITE(*,*)'ERROR: init_drive_var: variable=',driveVarname(ivar)
        WRITE(*,*)'We do not have a valid combination for this variable.'
        WRITE(*,*)'Note that we have driveVarUse=.TRUE.'
        WRITE(*,*)'This is a bug in the code!'
      ENDIF

    ELSE   !  NOT driveVarUse

!     Either driveVarPosIn>0, driveVarPos>=0
!     OR     driveVarPosIn=0, driveVarPos=0
!     Check if we have a valid state.
      IF ( .NOT. (driveVarPosIn(ivar)>0 .AND. driveVarPos(ivar)>=0 ) .AND.  &
           .NOT. (driveVarPosIn(ivar)==0 .AND. driveVarPos(ivar)==0 ) .AND.  &
           .NOT. (driveVarPosIn(ivar)==0 .AND. driveVarPos(ivar)>0 ) ) THEN
        errFound = .TRUE.
        WRITE(*,*)'ERROR: init_drive_var: variable=',driveVarname(ivar)
        WRITE(*,*)'We do not have a valid combination for this variable.'
        WRITE(*,*)'Note that we have driveVarUse=.FALSE.'
        WRITE(*,*)'This is a bug in the code!'
      ENDIF

    ENDIF   !  driveVaruse

!-------------------------------------------------------------------------------
!   Now check for errors in input.
!-------------------------------------------------------------------------------
!   Check that all required input variables have been found.
    IF ( driveVarPosIn(ivar)>0 .AND. driveVarFlag(ivar)==0 ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: init_drive_var: variable=',driveVarname(ivar)
      WRITE(*,*)'  This is a required input, but has not been found.'
    ENDIF

!   Check that we don't have variables that are not required - this suggests that
!   the user does not understand what's required.
    IF ( driveVarFlag(ivar)/=0 .AND. driveVarPosIn(ivar)==0 ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: init_drive_var: variable=',driveVarname(ivar)
      WRITE(*,*)'  This variable is not required to be input, but has been found.'
    ENDIF

!   If a variable is to be read in (other than from a SDF file), we must have a valid location.
    IF ( driveFormat/=formatNc .AND. driveVarPosIn(ivar)>0 .AND. driveVarFlag(ivar)<1 ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: init_drive_var: variable=',driveVarname(ivar)
      WRITE(*,*)'  driveVarFlag must be > 0 (this is field number on file)'
      WRITE(*,*)'  This is specified in the run control file.'
      WRITE(*,*)'  Current status suggests variable is not listed in run control file.'
    ENDIF

  ENDDO

!-------------------------------------------------------------------------------
! If an error in the driving variables setup was detected, stop.
!-------------------------------------------------------------------------------
  IF ( errFound ) STOP

!-------------------------------------------------------------------------------
! Get the number of files that are to be open at any time.  This is one, unless
! variable-name-templating is used, in which case need to count.
!-------------------------------------------------------------------------------
  ndriveUnit = 1
  IF ( driveTemplateV ) THEN
!   Need to establish how many different files are used.
    ndriveUnit = 0
    done(:) = .FALSE.
!   Don't consider variables that are not to be read in.
    WHERE( driveVarFlag(:) < 1 ) done(:)=.TRUE.
    DO ivar=1,ndriveVarMax
      IF ( .NOT. done(ivar) ) THEN
!       This is a new file.
        ndriveUnit = ndriveUnit + 1
!       Set flag at this and all other variables that use same file.
        WHERE( driveVarFileName(:)==driveVarFileName(ivar) ) done(:)=.TRUE.
      ENDIF
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Allocate some more space.
!-------------------------------------------------------------------------------
  CALL allocate_arrays( 'init_drive 2' )

!-------------------------------------------------------------------------------
! For each driving variable, associate a unit.
!-------------------------------------------------------------------------------
  i = 0
  driveUnitUse(:) = 1    !  all variables use driveUnit(1)
  IF ( ndriveUnit > 1 ) THEN
    done(:) = .FALSE.
!   Don't consider variables that are not to be read in.
    WHERE( driveVarFlag(:) < 1 ) done(:)=.TRUE.
    DO ivar=1,ndriveVarMax
      IF ( .NOT. done(ivar) ) THEN
!       This is a new file.
        i = i + 1
!       Set unit at this and all other variables that use same file.
        WHERE( driveVarFileName(:)==driveVarFileName(ivar) )
          done(:)=.TRUE.
          driveUnitUse(:) = i
        endwhere
!       Save template subsitution string for this unit.
        driveVarNameUnit(i) = driveVarFileName(ivar)
      ENDIF
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Process the variable positions.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Get a sorted list, giving input variables according to position in input
! file(s).  This is needed for file formats for which the order matters (or at
! least can be made easier).
!-------------------------------------------------------------------------------
  done(:) = .TRUE.

! Only consider variables that are to be read in (flag>0). SDF variables have flag set >0.
  WHERE ( driveVarFlag(:) > 0 ) done(:)=.FALSE.
  DO ivar=1,ndriveVarIn
    ival(1:1) = MINLOC( driveVarFlag(:), .NOT.done(:) )
    i = ival(1)
    done(i) = .TRUE.
    driveVarInSort(ivar) = i
  ENDDO

!-------------------------------------------------------------------------------
! Set flags indicating what to check, depending upon file format.
!-------------------------------------------------------------------------------
  checkNames = .FALSE.
  checkPos = .FALSE.
  SELECT CASE ( driveFormat )
    CASE ( formatAsc )
!     ASCII. Need positions of variables.
      checkPos = .TRUE.
    CASE ( formatBin )
!     Binary. Need positions of variables.
      checkPos = .TRUE.
    CASE ( formatNc )
!     netCDF. Need names of variables.
      checkNames = .TRUE.
    CASE ( formatPP )
!     PP format. Need positions of variables.
      checkPos = .TRUE.
    CASE default
      WRITE(*,*)'ERROR: init_drive_var: no code for driveFormat=',TRIM(driveFormat)
      STOP
  END SELECT

!-------------------------------------------------------------------------------
! Check positions of variables within file (not done for SDFs).
!-------------------------------------------------------------------------------
  IF ( checkPos ) THEN
!   Check each file separately.
    DO i=1,ndriveUnit
      tmpFlag(:) = 0
      WHERE ( driveUnitUse(:) == i ) tmpFlag(:) = driveVarFlag(:)
!     Set flag to zero for variables that are not read in (even if indicated as
!     being in this file).
      WHERE( driveVarFlag(:) < 1 ) tmpFlag(:)=0
!     First argument to checkVarPos indicates how important the order of the
!     variables is. Zero means order is not important.
      ivar = checkVarPos( 0,tmpFlag(:),' init_drive_var: tmpFlag'  &
                         ,nfield=nfieldDriveFile,geZero=.TRUE. )
      IF ( ivar < 0 ) THEN
        WRITE(*,*)'ERROR: init_drive: error from checkvarPos.'
        IF ( ndriveUnit > 1 ) THEN
          WRITE(*,*)'Error COULD be caused by mistakenly listing the same'
          WRITE(*,*)'file name and field number for different variables.'
          WRITE(*,*)'Error indicated for substitution=',driveVarNameUnit(i)
        ENDIF
        WRITE(*,*)'If error was repeated use of same varPos, but you do in fact want to reuse the'
        WRITE(*,*)'same data for more than one variable, comment out this stop!'
        STOP
      ENDIF
    ENDDO  !  units
  ENDIF

!------------------------------------------------------------------------------
! If necessary, check that names of variables are unique.
! Note that there is not yet any code to allow variables in different files to
! have the same name - a fairly unlikely situation anyway.
!-------------------------------------------------------------------------------
  IF ( checkNames ) THEN
    DO ivar=1,ndriveVarIn
      ipos = driveVarInSort(ivar)
      DO jvar=ivar+1,ndriveVarIn
        jpos = driveVarInSort(jvar)
        IF ( driveVarNameSDF(ipos) == driveVarNameSDF(jpos) ) THEN
          WRITE(*,*)'ERROR: init_drive_var: repeated driveVarNameSDF: ',TRIM(driveVarNameSDF(ipos))
          WRITE(*,*)'This means that a variable in a SDF file is used for more'
          WRITE(*,*)'than one driving variable. Possible, but unlikely.'
          WRITE(*,*)'Also unlikely and not coded for - same name is used in >1 file.'
          STOP
        ENDIF
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Count number of variables stored in driveData array.
!-------------------------------------------------------------------------------
  ndriveVar = 0
  DO ivar=1,ndriveVarMax
    IF ( driveVarPos(ivar) > 0 ) ndriveVar = ndriveVar + 1
  ENDDO

!-------------------------------------------------------------------------------
! Check for error: positions 1:ndriveVar should be used.
!-------------------------------------------------------------------------------
  IF ( MAXVAL(driveVarPos(:)) /= ndriveVar ) THEN
    WRITE(*,*)'ERROR: init_drive_var: maxval(driveVarPos) /= ndriveVar.'
    WRITE(*,*)'This suggests a coding error. Expect driveVarPos to be 1 to ndriveVar'
    WRITE(*,*)'driveVarPos=',driveVarPos(:)
    WRITE(*,*)'ndriveVar=',ndriveVar
    STOP
  ENDIF

END SUBROUTINE init_drive_var




!################################################################################
!################################################################################
! subroutine init_drive_var_detail
! Writes message to screen describing certain details of driving data.
! Possibly could be combined with varInfo.
!################################################################################
!################################################################################
SUBROUTINE init_drive_var_detail( driveVarName,varDesc )

  USE drive_io_vars, ONLY :  &
!   imported scalars with intent(in)
      driveFormat,drivetemplateV,ndriveVarMax  &
!   imported arrays with intent(in)
     ,driveUnitUse,driveVarFlag,driveVarInterp,driveVarNameSDF,driveVarNameUnit  &
     ,driveVarPos,driveVarPosIn,driveVarStash,driveVarUse

  USE inout, ONLY :  &
!   imported scalar parameters
      formatNc,formatPP

  IMPLICIT NONE

!-------------------------------------------------------------------------------
  CHARACTER(len=*), INTENT(in) ::  &!  in arrays
    driveVarName(ndriveVarMax)   &!  names of all possible forcing variables
   ,varDesc(ndriveVarMax)         !  description of each variable

  INTEGER ::  &!  scalars
    ipos   &!  work
   ,ivar    !  loop counter
!-------------------------------------------------------------------------------


  DO ivar=1,ndriveVarMax
!-------------------------------------------------------------------------------
!   For variables either required by model or used to derive another variable,
!   write details of variable setup to stdout.
!-------------------------------------------------------------------------------
    IF ( driveVarUse(ivar) .OR. (.NOT.driveVarUse(ivar) .AND.driveVarFlag(ivar)>0) ) THEN

      WRITE(*,*) TRIM(driveVarName(ivar)),'  (',TRIM(varDesc(ivar)),')'

      IF ( driveVarUse(ivar) .AND. driveVarPosIn(ivar)==0 ) THEN
         WRITE(*,*)'....is not read in, but is calculated from other variables.'
      ELSEIF ( .NOT.driveVarUse(ivar) .AND. driveVarFlag(ivar)>0 ) THEN
          WRITE(*,*)'....is only used to derive another variable (not used itself).'
      ENDIF

      IF ( driveVarFlag(ivar) > 0 ) THEN
!       Variable is read in.
        IF ( driveFormat == formatPP ) WRITE(*,*)'....uses STASH code=',driveVarStash(ivar)
        IF ( driveFormat /= formatNc ) THEN
          WRITE(*,*) '....is stored in file field #',driveVarFlag(ivar)
        ELSE
          WRITE(*,*) '....variable is ',TRIM(driveVarNameSDF(ivar))
        ENDIF
!       Variable-name templating.
        IF ( driveTemplateV ) THEN
          ipos = driveUnitUse(ivar)
          WRITE(*,*)'....comes from file ',TRIM(driveVarNameUnit(ipos))
        ENDIF

!       Time interpolation.
        SELECT CASE ( driveVarInterp(ivar) )
          CASE ( 'b' )
            WRITE(*,*)'....data are backward averages, 3 time levels needed for interpolation'
          CASE ( 'c' )
            WRITE(*,*)'....data are centred averages, 4 time levels needed for interpolation'
          CASE ( 'f' )
            WRITE(*,*)'....data are forward averages, 3 time levels needed for interpolation'
          CASE ( 'i' )
            WRITE(*,*)'....data are instantaneous values, 2 time levels needed for interpolation'
          CASE ( 'nb' )
            WRITE(*,*)'....data are backward averages, no interpolation applied'
          CASE ( 'nc' )
            WRITE(*,*)'....data are centred averages, no interpolation applied, 2 time levels stored'
          CASE ( 'nf' )
            WRITE(*,*)'....data are forward averages, no interpolation applied'
          CASE default
            WRITE(*,*)'WARNING: init_drive_var_detail: no code for interpFlag=',TRIM(driveVarInterp(ivar))
        END SELECT

        WRITE(*,*)'....is stored in driveDataIn location #',driveVarPosIn(ivar)

      ENDIF   !   driveVarFlag

      IF ( driveVarPos(ivar) > 0 )  &
         WRITE(*,*)'....is stored in driveData location #',driveVarPos(ivar)

    ENDIF  !  required or used to derive
  ENDDO

END SUBROUTINE init_drive_var_detail



!################################################################################
!################################################################################
!################################################################################
