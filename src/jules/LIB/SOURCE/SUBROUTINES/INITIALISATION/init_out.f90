!###############################################################################
!###############################################################################
! This file contains:
!   SUBROUTINE init_out.  Driver for initialisation of output.
!###############################################################################
!###############################################################################

  SUBROUTINE init_out(MYNUM,hfilout)

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY :   &
!   imported scalars with intent(in)
     sm_levels

  USE file_utils, ONLY :  &
!   imported procedures
     fileUnit,findTag

  USE init_out_map_mod, ONLY :  &
!   imported procedures
     init_out_map,init_out_map_land,init_out_map_compress

  USE inout, ONLY : compressgridFile,coord,coordLL,dumpFormat,dumpFreq,dumpStatus,dumpWarn  &
          ,echo,formatAsc,formatBin,formatNc,gradsNc,mapOut,mapOutCompress,mapOutLand  &
          ,nlevMax,nlevMaxCtl,nout  &
          ,ntCtl,ntCtlNeed,ntOutFilePer,ntOutPer  &
          ,numMonth,nvar,nvarOut,nvarOutTot,nxyMax,outActivePrev  &
          ,outCtlFile,outDataFile,outDate,outDateFlag,outDir  &
          ,outEndian   &
          ,outFilePer,outFirstActive,outFirstSection,outFirstWrite   &
          ,outGrads,outGridNxy,outLen,outName  &
          ,outPer,outRangeX,outRangeY,periodAnn,periodMon  &
          ,outAreaLL,outGridXY,outLLorder,outSamPer,outStatus  &
          ,outFormat,outTime,outUnit,outVal,outWarnCtl,outWarnEarly,outWarnUnder &
          ,outWriteCount,periodOneFile,pointsFlag,pointsOut,pointsOutLand,pointsOutMax  &
          ,redoTdef,rgProfile,rpProfile,runID  &
          ,snapProfile,jinUnit,sufLen,taccumVar,outTemplate,tmeanProfile,tmeanVar,undefOut  &
          ,varname,outCompress,useCompressGrid,varDesc,varDesclist,varNameList,varNlev,varNum  &
          ,varPos,varStartPos,varType,varTypeList,yrevOut  &
          ,zrevOutSnow,zrevOutSoil

  USE misc_utils, ONLY :   &
!   imported procedures
     getWord

  USE offline_diag, ONLY :  &
!  imported scalars with intent(out)
     useCiDiag,useGstomDiag,useRdcDiag,useRFlowDiag,useRoffInfDiag  &
    ,useRRunDiag,useSnowGMeltDiag,useWfluxDiag,useWfluxSfcDiag

  USE soil_param, ONLY :    &
!  imported scalars with intent(in)
     zsmc,zst  &
!  imported arrays with intent(in)
    ,dzsoil

  USE route_mod, ONLY :   &
!   imported scalars with intent(in)
     routeTimeStep

  USE spin_mod, ONLY:   &
!   imported scalars with intent(in)
      nspin

  USE switches, ONLY :    &
!   imported scalars with intent(in)
     l_360,route,routeOnly

  USE time_loc, ONLY :    &
!   imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :    &
!   imported scalars with intent(in)
     s_to_chhmmss

  USE timeConst, ONLY : &
     iSecInHour

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  local SCALARS
    nvarMax = 500   !  the maximum possible number of output variables (summed over all profiles)

  INTEGER ::   &!  local SCALARS
    i1,i2      &!  loop counter/work
   ,iout       &!  loop counter
   ,ivar,jvar,kvar !  work/loop counter

  INTEGER ::  &!  local ARRAYS
    tmpNum(nvarMax)  !  work: location of variable in list

  LOGICAL ::  &!  local SCALARS
    dirExist  &!  T if output directory exists
   ,found     &!  T if a variable is recognised
   ,useTemplate  !  T means GrADS output will use template option when possible

  LOGICAL ::  &!  local ARRAYS
    tmpTaccum(nvarMax)   &!  work: flag for time accumulation
   ,tmpTmean(nvarMax)     !  work: flag for time average

  CHARACTER(len=1) ::  &!  local SCALARS
    char1    !   work

  CHARACTER(len=100) ::  &!  local SCALARS
    inLine    !  a line read from file

  CHARACTER(len=LEN(varTypeList)) ::  &!  local arrays
    tmpvarType(nvarMax)   !  work: type of variable

  CHARACTER(len=LEN(varNameList)) ::  &!  local SCALARS
    tmpName       !  work: name of a selected variable

  CHARACTER(len=LEN(varName)) ::  &!  local arrays
    tmpvarname(nvarMax)   !  work: name used for a selected variable

  CHARACTER(len=sufLen) ::  &!  local ARRAYS
    tmpSuffix(nvarMax)    !  work: suffix for variable descriptions and names

  CHARACTER(len=100) ::  &!  local ARRAYS
    tmpLoc(nvarMax,2)    !  work: grid location fields

  CHARACTER(len=10) ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
!                                 ! that is required because some compilers can not cope
!                                 ! with recursive write statements

  INTEGER :: MYNUM
  CHARACTER(len=256) ::  hfilout

!-----------------------------------------------------------------------
  if (echo) WRITE(*,"(50('-'),/,a)") 'init_out'

!------------------------------------------------------------------------
! Initialisation.
!------------------------------------------------------------------------
  dumpWarn     = .FALSE.
  outWarnCtl   = .FALSE.
  outWarnEarly = .FALSE.
  outWarnUnder = .FALSE.
  pointsOutMax = 0
  tmpTaccum(:) = .FALSE.
  tmpTmean(:)  = .FALSE.
!------------------------------------------------------------------------
! Set switches to indicate that none of the "offline" diagnostics are
! required.
!------------------------------------------------------------------------
  useCiDiag       = .FALSE.
  useGstomDiag    = .FALSE.
  useRdcDiag      = .FALSE.
  useRFlowDiag    = .FALSE.
  useRoffInfDiag  = .FALSE.
  useRRunDiag     = .FALSE.
  useSnowGMeltDiag = .FALSE.
  useWfluxDiag    = .FALSE.
  useWfluxSfcDiag = .FALSE.

!------------------------------------------------------------------------
! Locate the start of this section in input file.
!------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_out','>INIT_OUT' )

!------------------------------------------------------------------------
! Read general details.
!------------------------------------------------------------------------
  READ(jinUnit,*) runID
  READ(jinUnit,*) outDir

!!DSM{ ----- Utilizando a configurcao HISTORY do RAMSIN
  IF (len(trim(hfilout(index(hfilout,'/',BACK=.TRUE.)+1:100))) > 9) THEN
     PRINT*, '********'
     PRINT*, 'ERRO... O nome do arquivo history (em HFILOUT) deve ter no maximo 5 caracteres'
     PRINT*, '********'
     STOP
  ENDIF
  !}

  !{--- Configurando o nome do history (dump) do JULES ---
  !outDir=trim(hfilout(1:index(hfilout,'/',BACK=.TRUE.)-1))
  !write(runID,'(a,i4.4)') trim(hfilout(index(hfilout,'/',BACK=.TRUE.)+1:100))//'-',MYNUM

  write(outDir,'(a,i4.4)') trim(hfilout(1:index(hfilout,'/',BACK=.TRUE.)-1))//'/p',MYNUM
  runID=trim(hfilout(index(hfilout,'/',BACK=.TRUE.)+1:100))
!DSM <BRAMS05 nao possui history>  CALL SYSTEM('mkdir -p '//outDir//' 2>/dev/null')
  !}
!DSM}

  READ(jinUnit,*) dumpFreq
dumpFreq=0  !BRAMS_05 nao possui history

  READ(jinUnit,*) dumpFormat
  READ(jinUnit,*) dumpStatus

  READ(jinUnit,*) nout
  READ(jinUnit,*) outFormat
  READ(jinUnit,*) gradsNc
  READ(jinUnit,*) outStatus
  READ(jinUnit,*) yrevOut
  READ(jinUnit,*) zrevOutSoil,zrevOutSnow
  READ(jinUnit,*) numMonth
  READ(jinUnit,*) useTemplate
  READ(jinUnit,*) undefOut
  READ(jinUnit,*) zsmc,zst
  READ(jinUnit,*) outEndian




!------------------------------------------------------------------------
! Deal with dump variables first, then other output, since there may not
! be any other output.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! Check we recognise frequency for dumps.
!------------------------------------------------------------------------
  SELECT CASE ( dumpFreq )
    CASE ( 0:4 )   !  valid
    CASE default
      WRITE(*,*)'ERROR: init_out: Invalid value for dumpFreq. Value=',dumpFreq
      STOP
  END SELECT

!------------------------------------------------------------------------
! Look to simplify dumpFreq if there is no spin up.
!------------------------------------------------------------------------
  IF ( nspin==0 .AND. dumpFreq==3 ) dumpFreq=2

!------------------------------------------------------------------------
! Check format of dump files.
!------------------------------------------------------------------------
  SELECT CASE ( dumpFormat )
    CASE ( formatAsc,formatNc )
!     Acceptable.
    CASE default
      WRITE(*,*)'ERROR: init_out: dumpFormat not allowed.'
      WRITE(*,*)'Dumps must be ASCII or netCDF.'
      STOP
  END SELECT

! In this version, don't allow netCDF dumps if routing is selected.
! I'm in a hurry and this simplifies grid information.
  IF ( dumpFormat==formatNc .AND. route ) THEN
    WRITE(*,*)'ERROR: init_out: netCDF dumps not allowed with routing.'
    WRITE(*,*)'Sorry - more code is required.'
    STOP
  ENDIF

!------------------------------------------------------------------------
! Check the file status for dump files.
!------------------------------------------------------------------------
  SELECT CASE ( dumpStatus )
    CASE ( 'new', 'replace' )
    CASE default
      WRITE(*,*)'ERROR: init_out: dumpStatus not allowed.'
      WRITE(*,*)'Acceptable values are: new, replace'
      STOP
  END SELECT

!------------------------------------------------------------------------
! Check that the output directory exists, if it is needed (for dumps or
! other output).
!------------------------------------------------------------------------
  IF ( dumpFreq/=0 .OR. nout>1 ) THEN
    IF ( LEN_TRIM(outdir) == 0 ) THEN
      WRITE(*,*)'ERROR: init_out: Invalid (empty) specification of output directory.'
      WRITE(*,*) 'Enter "." to use current directory.'
      STOP
    ENDIF
    INQUIRE (file=outDir, exist=dirExist )
!------------------------------------------------------------------------
!   If using an Intel compiler (tested for v9.0), you probably need to
!   remove the previous line and replace with the next line.
!    INQUIRE (directory=outDir, exist=dirExist )
!------------------------------------------------------------------------
    IF ( .NOT. dirExist ) THEN
      WRITE(*,*)'ERROR: init_out: Output directory does not exist.'
      WRITE(*,*)'outDir:',TRIM(outDir)
      STOP
    ENDIF
  ENDIF

!------------------------------------------------------------------------
! Issue a warning if no output is selected, then leave this subroutine.
!------------------------------------------------------------------------
  IF ( nout < 1 ) THEN
    if (echo) WRITE(*,"(50('#'),/,a,/,50('#'))")'WARNING: init_out: No output has been requested!'
    RETURN
  ENDIF

!------------------------------------------------------------------------
! Check the file status.
!------------------------------------------------------------------------
  SELECT CASE ( outStatus )
    CASE ( 'new', 'replace' )
    CASE default
      WRITE(*,*)'ERROR: outStatus not allowed.'
      WRITE(*,*)'Acceptable values are: new, replace'
      WRITE(*,*)'Stopping in init_out'
      STOP
  END SELECT

!------------------------------------------------------------------------
! Check that the output style is recognised.
!------------------------------------------------------------------------
  SELECT CASE (outFormat)
    CASE ( formatAsc,formatBin,formatNC )
!     OK, no problem.
    CASE default
      WRITE(*,*)'Format for output=',TRIM(outFormat),' .Not recognised.'
      WRITE(*,*)'Stopping in init_out'
      STOP
  END SELECT

!------------------------------------------------------------------------
! Set flag indicating that output is potentially useable with GrADS.
!------------------------------------------------------------------------
  outGrADS = .FALSE.
  SELECT CASE (outFormat)
    CASE ( formatBin )
      outGrads = .TRUE.
    CASE ( formatNc )
      IF ( gradsNc ) outGrads = .TRUE.
  END SELECT

! For clarity (and in fact used as a shorthand later on!), don't indicate
! GrADS-readable netCDF files if not netCDF output.
  IF ( outFormat /= formatNc ) gradsNc = .FALSE.

!------------------------------------------------------------------------
! Ignore (reset) certain options if not using GrADS output.
!------------------------------------------------------------------------
  IF ( outGrads ) THEN
    redoTdef = .TRUE.
  ELSE
    redoTdef = .FALSE.
    useTemplate = .FALSE.
  ENDIF

  IF ( .NOT. routeOnly ) THEN
!------------------------------------------------------------------------
!   Check that zsmc and zt are within soil depth.
!------------------------------------------------------------------------
    IF ( SUM( dzsoil(1:sm_levels) ) < zsmc ) THEN
      WRITE(*,*)'ERROR: init_out: zsmc is below bottom of soil column'
      WRITE(*,*)'Stopping in init_out'
      STOP
    ENDIF
    IF ( SUM( dzsoil(1:sm_levels) ) < zst ) THEN
      WRITE(*,*)'ERROR: init_out: zt is below bottom of soil column'
      WRITE(*,*)'Stopping in init_out'
      STOP
    ENDIF
  ENDIF  !  routeOnly

!------------------------------------------------------------------------
! Get a list of available variables.
!------------------------------------------------------------------------
  CALL init_out_varList

!------------------------------------------------------------------------
! Allocate memory for arrays dimensioned with nout.
!------------------------------------------------------------------------
  CALL allocate_arrays( 'init_out 1' )

!------------------------------------------------------------------------
! Initialise.
!------------------------------------------------------------------------
  nvarOutTot = 0
  outLen     = 0

  nlevMax(:)       = -1
  nlevMaxCtl(:)    = -1
  ntCtl(:)         = -1
  ntCtlNeed(:)     = -1
  ntOutFilePer(:)  = -1
  ntOutPer(:)      = -1
  nvarOut(:)       =  0
  nxyMax(:)        = -1
  outWriteCount(:) =  0
  pointsOutLand(:) =  0

  outActivePrev(:)   = .FALSE.
  outFirstActive(:)  = .TRUE.
  outFirstSection(:) = .TRUE.
  outFirstWrite(:)   = .TRUE.
  rgProfile(:)       = .FALSE.
  rpProfile(:)       = .FALSE.
  snapProfile(:)     = .FALSE.
  outTemplate(:)     =  useTemplate
  tmeanProfile(:)    = .FALSE.
  outCompress(:)     = .FALSE.

  outCtlFile(:)  = ''
  outDataFile(:) = ''

!------------------------------------------------------------------------
! Read the profiles.  Lines are read into temporary space until details
! (such as number of variables and size of output grid) can be calculated.
!------------------------------------------------------------------------
  DO iout=1,nout

!   Locate start of this profile.
    CALL findTag( jinUnit,'init_out >NEWPROF','>NEWPROF',.TRUE. )

!   Read name for profile.
    READ(jinUnit,*) outName(iout)

    IF ( LEN_TRIM(outName(iout)) == 0 ) THEN
      WRITE(*,*)'ERROR: Missing name for output profile.'
      WRITE(*,*)'Error for output profile #',iout
      WRITE(*,*)'Stopping in init_out'
      STOP
    ENDIF

!------------------------------------------------------------------------
!   Process dates and times.
!------------------------------------------------------------------------
    CALL init_out_time( iout )

!------------------------------------------------------------------------
!   Read details of output grid and the output mapping.  These are read,
!   but not much more done until we know what variables are selected.
!------------------------------------------------------------------------
    CALL init_out_map( iout,1,tmpLoc,tmpSuffix )

!------------------------------------------------------------------------
!   Read the list of selected variables.
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
!   Locate the start of the output variables section.
!------------------------------------------------------------------------
    CALL findTag( jinUnit,'init_out >VARS','>VARS',.TRUE. )

!------------------------------------------------------------------------
!   Read a list of selected variables and other details.  Everything is
!   read as a character string, since different types of variable have
!   different numbers of fields in this line.   A line containing a
!   variable is formatted as,
!
!   flag name [name2 xloc yloc desc]
!   where flag is a character indicating instantaneous or time-average value
!       name is the (short) name of the variable
!       name2 is the name to be used for the variable
!   and the following are only present for point variables on the routing grid
!       xloc,yloc give a grid location
!       desc is a description (e.g. flow at outlet).
!   Fields are separated by one or more blanks.
!
!   List ends at >ENDVARS.
!------------------------------------------------------------------------
    IF ( iout==1 ) ivar = 0
    DO
      READ(jinUnit,"(a)") inLine
      IF ( inLine == '>ENDVARS' ) EXIT  !  end of list
      IF ( LEN_TRIM(inLine) == 0 ) CYCLE   !  skip blank lines
      ivar = ivar + 1
      IF ( ivar>nvarMax ) THEN
        WRITE(*,*)'ERROR: init_out: too many variables! Increase nvarMax.'
        WRITE(*,*)'Stopping in init_out'
      ENDIF
!     Increment counters of total # of variables and variables in this profile.
      nvarOutTot = nvarOutTot + 1
      nvarOut(iout) = nvarOut(iout) + 1

!------------------------------------------------------------------------
!     Use the first field to establish whether variable is to be
!     time-averaged/accumulated.
!------------------------------------------------------------------------
      char1 = getWord( inline,1 )    !   only takes the first non-blank character
      SELECT CASE ( char1 )
        CASE ( 'A' )
          tmpTaccum(ivar) = .TRUE.
        CASE ( 'M' )
          tmpTmean(ivar) = .TRUE.
        CASE ( 'S' )
!         Nothing to do.
        CASE default
          WRITE(*,*)'ERROR: do not recognise (or cannot find) time flag :',char1
          WRITE(*,*)'in line:',TRIM(inline)
          WRITE(*,*)'NB One possible reason is that the list of variables'
          WRITE(*,*)'does not end with >ENDVARS.'
          WRITE(*,*)'Stopping in init_out'
          STOP
      END SELECT

!------------------------------------------------------------------------
!     Find out what variable is requested.  The name is the second field
!     (separated by blanks).
!------------------------------------------------------------------------
      tmpName = getWord( inline,2 )
      IF ( tmpName == '' ) THEN
        WRITE(*,*) 'ERROR: init_out: line:',TRIM(inLine)
        WRITE(*,*)'Can''t find variable name.'
        WRITE(*,*)'Stopping in init_out'
        STOP
      ENDIF

!------------------------------------------------------------------------
!     Check that this name is recognised and establish type of variable.
!------------------------------------------------------------------------
      found = .FALSE.
      DO jvar=1,nvar
        IF ( tmpName == varNameList(jvar) ) THEN
          found=.TRUE.
          tmpNum(ivar) = jvar
          tmpVarType(ivar) = varTypeList(jvar)
          EXIT
        ENDIF
      ENDDO
      IF ( .NOT. found ) THEN
        WRITE(*,*)'ERROR: init_out: tmpName=',tmpName
        WRITE(*,*)'Do not recognise the variable called ',TRIM(tmpName)
        STOP
      ENDIF
!------------------------------------------------------------------------
!     Get the name to be used for output, if none provided, reuse the
!     variable name.
!------------------------------------------------------------------------
      tmpvarname(ivar) = getWord( inline,3 )
      IF ( tmpvarname(ivar)=='-' .OR. LEN_TRIM(tmpvarname(ivar)) == 0 ) tmpvarname(ivar)=tmpName

!-----------------------------------------------------------------------
!     If routing variables are selected for a point, change type of
!     variable to RP.  RP variables require further fields: 2 giving
!     location, 1 giving name used in annotation.
!     It is possible to select a single point via pointsOut=1, but we do
!     not include that case here - since it's limited to one point per
!     profile.
!------------------------------------------------------------------------
      IF ( tmpVarType(ivar) == 'RG' ) THEN
        tmpLoc(ivar,1) = getWord( inline,4 )
        tmpLoc(ivar,2) = getWord( inline,5 )
        tmpSuffix(ivar) = getWord( inline,6 )
!       We need to get either 3 values or none.
        IF ( tmpLoc(ivar,1) /= '' ) THEN
          IF ( tmpSuffix(ivar)=='' ) THEN
            WRITE(*,*) 'ERROR: init_out: line:',TRIM(inLine)
            WRITE(*,*) 'We have more than 3 fields, so assuming this is a routing variable at a point.'
            WRITE(*,*) 'But can''t find all of fields #4-6.'
            WRITE(*,*) '4 & 5 give location, 6 is a name used in annotation.'
            STOP
          ENDIF
          tmpVarType(ivar) = 'RP'    !  this is a routing variable at a point

!------------------------------------------------------------------------
!         Check that the last characters of locations are N and E
!         respectively - they should be lat and lon.
!------------------------------------------------------------------------
          i1 = LEN_TRIM( tmpLoc(ivar,1) )
          i2 = LEN_TRIM( tmpLoc(ivar,2) )
          IF ( tmpLoc(ivar,1)(i1:i1)/='N' .OR. tmpLoc(ivar,2)(i2:i2)/='E' ) THEN
            WRITE(*,*) 'ERROR: init_out: line:',TRIM(inLine)
            WRITE(*,*) 'For a routing variable at a point, fields 4 and 5 should be lat and lon.'
            WRITE(*,*) 'They must end with N and E respectively.'
            STOP
          ENDIF
        ENDIF
      ENDIF  !  routing var

!     If this is a routing variable, set type of profile.
      IF ( tmpVarType(ivar) == 'RG' ) rgProfile(iout) = .TRUE.
      IF ( tmpVarType(ivar) == 'RP' ) rpProfile(iout) = .TRUE.

!     Check that a profile doesn't contain variables of incompatible types.
      IF ( rgProfile(iout) .AND. rpProfile(iout) ) THEN
        WRITE(*,*)'ERROR: init_out: profile cannot contain routing diagnostics'
        WRITE(*,*)'both on the grid and at listed points.'
        WRITE(*,*)'Error for profile #',iout,' called ',TRIM(outName(iout))
        WRITE(*,*)'Error raised at variable=',TRIM(tmpName)
        STOP
      ENDIF
      IF ( ( rgProfile(iout) .OR. rpProfile(iout) ) .AND.  &
           ( tmpVarType(ivar)/='RG' .AND.  tmpVarType(ivar)/='RP' ) ) THEN
        WRITE(*,*)'ERROR: init_out: a profile containing channel routing diagnostics'
        WRITE(*,*)'cannot also contain non-routing variables.'
        WRITE(*,*)'Error for profile #',iout,' called ',TRIM(outName(iout))
        WRITE(*,*)'Error raised at variable=',TRIM(tmpName)
        STOP
      ENDIF

    ENDDO  !  read a line

  ENDDO  !  iout (profiles)

!------------------------------------------------------------------------
! At this point we have finished reading necessary information from the
! input file.
!------------------------------------------------------------------------

!------------------------------------------------------------------------
! Check that profile names are unique.
!------------------------------------------------------------------------
  DO iout=1,nout-1
    IF ( ANY( outName(iout+1:nout)==outName(iout) ) ) THEN
      WRITE(*,*)'ERROR: Repeated output profile name: ',TRIM(outName(iout))
      WRITE(*,*)'Stopping in init_out'
      STOP
    ENDIF
  ENDDO

!------------------------------------------------------------------------
! Allocate further space now that number of variables is known.
!------------------------------------------------------------------------
  CALL allocate_arrays( 'init_out 2' )

!-------------------------------------------------------------------------------
! Load variables names etc into final storage.
!-------------------------------------------------------------------------------
  ivar = 0
  DO iout=1,nout
    DO jvar=1,nvarOut(iout)
      ivar = ivar + 1
      varPos(iout,jvar) = ivar
      varNum(ivar) = tmpNum(ivar)
      taccumVar(ivar) = tmpTaccum(ivar)
      tmeanVar(ivar) = tmpTmean(ivar)
      varname(ivar) = tmpvarname(ivar)
      varType(ivar) = tmpVarType(ivar)
      varDesc(ivar) = varDescList(varNum(ivar))
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Check that the requested variables are consistent with model
! configuration, e.g. don't ask for output from sub-model X if X has not been
! selected.
!-------------------------------------------------------------------------------
  CALL init_out_check

!------------------------------------------------------------------------
! Finish mappings.
!------------------------------------------------------------------------
  DO iout=1,nout
    CALL init_out_map( iout,2,tmpLoc,tmpSuffix )
  ENDDO

!------------------------------------------------------------------------
! Get mappings for land points.
!------------------------------------------------------------------------
! First, allocate space.
  CALL allocate_arrays( 'init_out 3' )
! Initialise.
  mapOutLand(:,:,:) = 0

! Get mappings.
  CALL init_out_map_land( nvarMax,tmpVarType )

!------------------------------------------------------------------------
! Further processing of variables.
!------------------------------------------------------------------------
  CALL init_out_var( nvarMax,tmpLoc )

!------------------------------------------------------------------------
! Allocate space used to store the output and initialise.
!------------------------------------------------------------------------
  CALL allocate_arrays( 'init_out 4' )
  outval(:) = 0.0

!------------------------------------------------------------------------
! Check sampling period and convert to timesteps.  Also work out how many
! sampling periods are in an output period.
!
! NB It would make sense to do this in the same procedure as the output
! period (outPer) is processed - currently that is in init_out_time.
!------------------------------------------------------------------------
  DO iout=1,nout

    IF ( tmeanProfile(iout) ) THEN

      IF ( outSamPer(iout) == 0 ) outSamPer(iout) = NINT(timeStep)
      IF ( outSamPer(iout) < 0 ) THEN
        WRITE(*,*)'ERROR: sampling period for output must be > 0.'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out'
        STOP
      ENDIF
      IF ( MOD(outSamPer(iout),NINT(timeStep))/=0 ) THEN
        WRITE(*,*)'ERROR: sampling period for output must be a multiple of timestep'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out'
        STOP
      ENDIF
      IF ( outPer(iout)>0 .AND. MOD(outPer(iout)*NINT(timeStep),outSamPer(iout))/=0 ) THEN
        WRITE(*,*)'ERROR: sampling period for output must be a factor of output period'
        WRITE(*,*)'Error for output profile #',iout
        WRITE(*,*)'Stopping in init_out'
        STOP
      ENDIF

!------------------------------------------------------------------------
!     Convert from seconds to # of timesteps.
!------------------------------------------------------------------------
      outSamPer(iout) = outSamPer(iout) / NINT(timeStep)

!------------------------------------------------------------------------
!     Work out how many sampling periods are in an output period. Only
!     needed for time average or accumulation.  This is estimated in
!     cases of "special" periods - e.g. months - and any error in estimate
!     is accounted for when time-average is calculated.
!
!     NOTE The estimation/assumption of constant month/year length is only
!     acceptable as long as the threshold for number of data that need to
!     be present before an average is calculated is set "not too high".
!     e.g., assume 30-day month, but have a 28-day Feb. If threshold is
!     set at 99.5%=29.85 days, the Feb average will always be set to
!     missing, since it can never get enough data. A 90% threshold will
!     be calculated when there are 27 or more days of data. If a high
!     threshold (e.g. 100% is required), ntOutPer will have to be
!     calculated separately at the start of each average, when we know
!     what month it is and therefore how long it is.  This is probably
!     not too tricky, except for any complicated cases if current date
!     refers to another month, e.g. at spin up?).
!------------------------------------------------------------------------
      SELECT CASE ( outPer(iout) )
        CASE ( periodAnn )
!         Estimate number of sampling periods in a year - assume 365 days.
          ntOutPer(iout) = ( 8760*iSecInHour ) / ( outSamPer(iout) * NINT(timeStep) )

        CASE ( periodMon )
!         Estimate number of sampling periods in a month - assume 30 days.
          ntOutPer(iout) = ( 720*iSecInHour ) / ( outSamPer(iout) * NINT(timeStep) )

        CASE ( 1: )
!         Constant output period.
          ntOutPer(iout) = outPer(iout) / outSamPer(iout)

        CASE default
          WRITE(*,*)'ERROR: init_out: no code for outPer(iout)=',outPer(iout)
          WRITE(*,*)'Stopping in init_out'
          STOP
      END SELECT

    ELSE

!     NOT tmeanProfile. Set to 1 timestep. This is not used but needs to pass certain tests.
      outSamPer(iout) = 1

    ENDIF  !  tmeanProfile

  ENDDO  !  iout

!------------------------------------------------------------------------
! Get file unit for each profile.  A scratch file is opened for each
! unit to reserve it (the unit is then in use at subsequent calls to
! fileUnit).
! In fact we could probably proceed without "reserving" units in this
! way...but this seems to be working for now.
! For netCDF files we don't need to reserve a unit, since a
! netCDF ID will later be given, but we need to set something here so that
! an attempt to close this output unit doesn't inadvertently close
! a netCDF file that is in use elsewhere! For netCDF, indicate an
! "impossible" unit (i.e. <0) which will never be used for i/o.
!------------------------------------------------------------------------
  IF ( outFormat == formatNc ) THEN
    outUnit(:) = -1
  ELSE
    DO iout=1,nout
!     Indicate ASCII format - just avoid netCDF!
      outUnit(iout) = fileUnit( formatAsc )
      PRINT*,'init_out: reserving scratch on ',outUnit(iout)
      OPEN(outUnit(iout),status='scratch')
    ENDDO
  ENDIF

!------------------------------------------------------------------------
! Check that "outEndian" type is OK.
!------------------------------------------------------------------------
  IF ( outFormat==formatBin .AND. outEndian/='big_endian' .AND. outEndian/='little_endian' ) THEN
    WRITE(*,*)'ERORR: init_out: outEndian must be big_endian or little_endian'
    WRITE(*,*)'Stopping in init_out'
    STOP
  ENDIF

!------------------------------------------------------------------------
! Write compression mapping (e.g. GrADS pdef) data files, if needed.
!------------------------------------------------------------------------
  CALL init_out_map_compress

!###############################################################################

!-----------------------------------------------------------------------
! Output to screen describing output.
!------------------------------------------------------------------------
  IF ( echo ) THEN

    WRITE(*,"(/70('-')/,a)") 'Dump (restart) file information:'
    SELECT CASE ( dumpFreq )
      CASE ( 0 )
        WRITE(*,*) 'dumpfreq=0: No model dumps will be written.'
      CASE ( 1 )
        WRITE(*,*) 'dumpfreq=1: Only the final state will be written to a dump.'
      CASE ( 2 )
        WRITE(*,*) 'dumpfreq=2: Initial and final states will be written to dumps.'
      CASE ( 3 )
        WRITE(*,*) 'dumpfreq=3: Dumps will be written for initial and final states, also after spin up.'
      CASE ( 4 )
        WRITE(*,*) 'dumpfreq=4: Dumps will be written for initial and final states, also at end of calendar years.'
    END SELECT

    WRITE(*,"(2a)") 'Dump file format will be ',dumpFormat

    WRITE(*,"(/70('-')/,a)") 'Selected output'
    WRITE(*,"(2a)") 'File format will be ',outFormat
    IF ( ANY(outTemplate(:)) ) WRITE(*,*)'Template ctl files will be used where possible.'
    IF ( l_360 ) WRITE(*,"(10('*'),a,10('*'))")  &
          ' WARNING: Using 360-day calendar: some visualisation software will struggle with this.'
    IF ( yrevOut .AND. ANY(.NOT.outCompress(:)) )  &
           WRITE(*,*)'Output will be written in N to S order.'

    DO iout=1,nout

      WRITE(*,"(/'Profile #',i2,': ',a)") iout,outName(iout)
      IF ( outFormat /= formatNc ) WRITE(*,*)'Output unit=',outUnit(iout)
      SELECT CASE ( outPer(iout) )
        CASE ( periodAnn );  WRITE(*,*)'Annual data.'
        CASE ( periodMon );  WRITE(*,*)'Monthly data.'
        CASE default;  WRITE(*,*)'Output period (s)=',outPer(iout)*NINT(timeStep)
      END SELECT
      SELECT CASE ( outFilePer(iout) )
        CASE ( periodOneFile );  WRITE(*,*)'All times to one file.'
        CASE ( -8 );  WRITE(*,*)'All times to one file, but spin-up times to a separate file'
        CASE ( -7 );  WRITE(*,*)'All times to one file, but each spin-up cycle in separate file'
        CASE ( periodAnn );  WRITE(*,*)'One file per year.'
        CASE ( periodMon );  WRITE(*,*)'One file per month.'
        CASE default; WRITE(*,*)'One file per ',outFilePer(iout)*NINT(timeStep),' s.'
      END SELECT

      IF ( outTemplate (iout) ) THEN
        WRITE(*,*) 'Template ctl files will be written.'
      ENDIF

      SELECT CASE ( outDateFlag(iout) )
        CASE ( 1: )
          time_hms = s_to_chhmmss( outTime(iout,1) )
          WRITE(*,*)'Start date and time for output:',outDate(iout,1),time_hms
          time_hms = s_to_chhmmss( outTime(iout,2) )
          WRITE(*,*)'End date and time for output:',outDate(iout,2),time_hms
        CASE ( 0 )
!         Dates may be misleading for this case (if nspin>0), so don't give them.
          WRITE(*,*)'Output is generated throughout the run (incl. any spin up).'
        CASE ( -1 )
          WRITE(*,*)'Output is generated at all times after spin up.'
          WRITE(*,*)'There is no output during spin up.'
        CASE ( -2 )
          WRITE(*,*)'Only the initial state will be output.'
      END SELECT

      IF ( .NOT. rpProfile(iout) ) THEN

        SELECT CASE ( pointsFlag(iout,1) )
          CASE ( 0 )
            WRITE(*,*)'All points in model grid will be included in output.'
          CASE ( 1 )
            WRITE(*,*)'A subarea of the model grid will be output, given by...'
            IF ( outAreaLL(iout) ) THEN
              WRITE(*,*)'latitude=',outRangeY(iout,1),' to ',outRangeY(iout,2)
             WRITE(*,*)'longitude=',outRangeX(iout,1),' to ',outRangeX(iout,2)
             ELSE
              WRITE(*,*)'column=',outRangeX(iout,1),' to ',outRangeX(iout,2)
              WRITE(*,*)'row=',outRangeY(iout,1),' to ',outRangeY(iout,2)
            ENDIF
          CASE ( 2 )
            WRITE(*,*)'Points to be output were given in a list.'
            IF ( coord(iout) ) THEN
              IF ( coordLL(iout) ) THEN
                WRITE(*,*)'Coordinates in list were latitude and longitude.'
              ELSE
                WRITE(*,*)'Coordinates in list were (x,y) on input grid.'
              ENDIF
            ENDIF
        END SELECT

        SELECT CASE ( pointsFlag(iout,2) )
          CASE ( 0 )
            WRITE(*,*)'The output grid is the same as the model grid.'
          CASE ( 1 )
            WRITE(*,*)'The output grid is the chosen subarea.'
          CASE ( 2 )
            WRITE(*,*)'The output grid described via the run control file.'
          CASE ( 3 )
            WRITE(*,*)'The output grid is the smallest rectangle that contained all output points'
            WRITE(*,*)'...the rectangle and order of points being determined by'
            IF ( outLLorder(iout) ) THEN
              WRITE(*,*)'...the latitude and longitude of each point.'
            ELSE
              WRITE(*,*)'...the row and column position in the input grid of each point.'
            ENDIF
          CASE ( 4 )
            WRITE(*,*)'The output grid is the same as the input grid.'
          CASE ( 5 )
            WRITE(*,*)'The selected points are simply written as a vector.'
        END SELECT

      ENDIF   !  rpProfile

      IF ( .NOT. outCompress(iout) ) THEN
        WRITE(*,*)'Size of output grid (x,y)=',outGridNxy(iout,1),outGridNxy(iout,2)
        WRITE(*,*)'Lat/lon of point (1,1) (SW corner)=',outGridXY(iout,2),outGridXY(iout,1)
        IF ( pointsOut(iout) /= outGridNxy(iout,1)*outGridNxy(iout,2) ) THEN
          WRITE(*,*)'....in which number of "defined" points=',pointsOut(iout)
          WRITE(*,*)'....and the number of padding points='  &
                         ,outGridNxy(iout,1)*outGridNxy(iout,2)-pointsOut(iout)
        ENDIF
      ELSE
        WRITE(*,*)'Output is a vector of model points that can be "scattered"'
        WRITE(*,*)'across a (potentially) larger grid.'
        WRITE(*,*)'....number of points in vector=',pointsOut(iout)
        WRITE(*,*)'....size of output grid (x,y)=',outGridNxy(iout,1),outGridNxy(iout,2)
        WRITE(*,*)'....lat/lon of point (1,1) (SW corner)=',outGridXY(iout,2),outGridXY(iout,1)
        WRITE(*,*)'The mapping between vector and grid is described by the data in file '  &
                ,TRIM(compressGridFile(iout))
        IF ( useCompressGrid(iout) /= iout ) WRITE(*,"(tr10,a,i2)")  &
                '....as written for profile #',useCompressGrid(iout)
        IF ( yrevOut ) WRITE(*,*)  &
           'NB This file is written in S to N order - i.e. yrevOut is ignored.'
      ENDIF

      IF ( rpProfile(iout) ) THEN
        WRITE(*,*) 'This grid is for selected points on the routing grid.'
      ELSEIF ( rgProfile(iout) ) THEN
        WRITE(*,*) 'This profile is for variables on the routing grid.'
      ENDIF

      WRITE(*,"(a,(/15i6))")'mapOut(points used from model grid)=',mapOut(iout,1:pointsOut(iout),1)
      IF ( .NOT. outCompress(iout) ) THEN
        WRITE(*,"(a,(/15i6))")'mapOut(points in output grid)=',mapOut(iout,1:pointsOut(iout),2)
      ELSE
        WRITE(*,"(a,(/15i6))")'mapOut(points in output grid)=',mapOutCompress(iout,1:pointsOut(iout))
      ENDIF

      IF ( pointsOutLand(iout) > 0 ) THEN
        WRITE(*,"(a)")'Mapping for land points:'
        WRITE(*,"(a,(/15i6))")'mapOutLand(1)=',mapOutLand(iout,1:pointsOutLand(iout),1)
        IF ( .NOT. outCompress(iout) ) WRITE(*,"(a,(/15i6))")'mapOutLand(2)='  &
                                ,mapOutLand(iout,1:pointsOutLand(iout),2)
      ENDIF

      WRITE(*,*)'Number of variables selected=',nvarOut(iout)
      IF ( .NOT. rpProfile(iout) ) THEN
        WRITE(*,*)'Variable #,varPos,varStartPos,nlev,tmean and taccum flags,name:'
      ELSE
        WRITE(*,*)'Variable #,varPos,varStartPos,tmean and taccum flags,grid index,name:'
      ENDIF
      DO ivar=1,nvarOut(iout)
        kvar = varPos(iout,ivar)   !  position in list of selected variables
        IF ( .NOT. rpProfile(iout) ) THEN
          WRITE(*,"(i3,tr1,i4,tr1,i7,tr1,i4,2(tr2,l1),tr1,a)")    &
                 ivar,varPos(iout,ivar),varStartPos(kvar),varNlev(kvar),tmeanVar(kvar)  &
                ,taccumVar(kvar),TRIM(varname(kvar))
        ELSE
          WRITE(*,"(i3,tr1,i4,tr1,i7,2(tr1,l1),tr1,i6,tr1,a)")  &
       &        ivar,varPos(iout,ivar),varStartPos(kvar),tmeanVar(kvar)   &
       &       ,taccumVar(kvar),mapOut(iout,ivar,1),TRIM(varname(kvar))
        ENDIF
      ENDDO

      WRITE(*,"(50('-'))")
    ENDDO  !  iout

  ENDIF  !  echo

!-------------------------------------------------------------------------------
! Warning messages.
!-------------------------------------------------------------------------------
! Warn if netCDF files might not be GrADS-readable.
  IF ( outFormat==formatNc .AND. .NOT.gradsNc ) THEN
    WRITE(*,"(/,a)")'WARNING: netCDF output but gradsNc=F'
    WRITE(*,*)'This means that the netCDF output might not be readable by GrADS.'
    WRITE(*,*)'e.g. snow layer variables will have 2 "z" dimensions.'
    WRITE(*,*)'This is NOT an issue if you don''t want to USE GrADS!'
  ENDIF

! A warning message to cover for the lack of more complete processing for
! diagnostics on the routing grid!
  IF ( ( ANY( rgProfile(:) ) .OR. ANY( rpProfile(:) ) ) .AND. &
       routeTimeStep /= 1 ) THEN
    WRITE(*,"(/,a)") 'WARNING: Take care when selecting timesteps for the'
    WRITE(*,*)'output of routing diagnostics. The code checks that all'
    WRITE(*,*)'timesteps for output are "sensible", e.g. are multiples'
    WRITE(*,*)'of the model timestep - but this is done assuming that the'
    WRITE(*,*)'"main" model timestep equals that of routing model, which'
    WRITE(*,*)'isn''t true for this run.'
    WRITE(*,*)'IN PARTICULAR, ensure that output and sampling period are'
    WRITE(*,*)'multiples of the timestep of the ROUTING model - otherwise'
    WRITE(*,*)'output can be misleading - e.g. if a time average is sampled'
    WRITE(*,*)'on timesteps when routing is not called, routing variables'
    WRITE(*,*)'will not have been updated since last routing timestep.'
  ENDIF
!-------------------------------------------------------------------------------

  END SUBROUTINE init_out

!###############################################################################
!###############################################################################
