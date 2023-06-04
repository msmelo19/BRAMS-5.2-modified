! Read parameters for PDM.
! This is currently very simple, as PDM only has scalar parameters, but much
! more code has been left, commented out, in anticipation of spatially-varying
! parameters! Of course it will be out of date by then...

  SUBROUTINE init_pdm()

!XX  USE ancil_info, ONLY :  &
!XX!  imported scalars with intent(in)
!XX     land_pts

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile,fileUnit,findTag,openFile  &
!  imported arrays with intent(inout)
    ,irecPrev

  USE inout, ONLY :  &
!  imported scalar parameters
     formatAsc,formatBin,formatLen,formatNc,formatPP,jinUnit,tagAscBin,tagNc  &
!  imported scalars with intent(in)
    ,echo,nxIn,nyIn  &
!  imported arrays with intent(in)
    ,mapInLand

!XX  USE jules_netcdf, ONLY :  &
!XX!  imported scalars with intent(in)
!XX     ncType
!XX
!XX  USE misc_utils, ONLY :  &
!XX!  imported procedures
!XX     checkVarPos,read_list,repeatVal,varInfo,varValue

  USE c_pdm, ONLY :  &
!  imported scalars with intent(out)
     b_pdm,dz_pdm

!XX  USE readwrite_mod, ONLY :  &
!XX!  imported procedures
!XX     readvar2dcomp
!XX
  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_pdm

!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar parameters.
!XX!-------------------------------------------------------------------------------
!XX  INTEGER, PARAMETER :: nvarMax = 4 !  the number of possible variables to be set
!XX!                            This may be larger than the number needed for any
!XX!                            particular configuration.
!XX!                            Possible variables:
!XX!                            fexp,ti_mean,ti_sig,ti_skew
!XX
!XX!-------------------------------------------------------------------------------
!XX! Local scalar variables.
!XX!-------------------------------------------------------------------------------
!XX  INTEGER ::        &
!XX       fieldNum     &!  field number in file
!XX      ,i,ierr       &!  loop counters/work
!XX      ,inUnit       &!  unit used to connect to input file
!XX      ,iorder       &!  work
!XX      ,ipos         &!  work
!XX!-------------------------------------------------------------------------------
!XX!  Every possible variable (nvarMax) must have an "iposX" variables declared here
!XX      ,iposFexp     &!  position in list of all possible variables of fexp
!XX      ,iposTi_mean  &!  position in list of all possible variables of ti_mean
!XX      ,iposTi_sig   &!  position in list of all possible variables of ti_sig
!XX      ,iposTi_skew  &!  position in list of all possible variables of ti_skew
!XX!-------------------------------------------------------------------------------
!XX      ,ivar,jvar    &!  work
!XX      ,l            &!  work
!XX      ,nfieldFile   &!  number of fields per time in a file
!XX      ,nheaderField &!  number of header records before each field in file
!XX      ,nheaderFile  &!  number of header records at start of file
!XX      ,nheaderT     &!  number of headers at start of each time
!XX      ,nlineField   &!  work
!XX      ,nvar      &!  number of variables required for chosen configuration
!XX      ,nvarIn       &!  number of variables that are read in
!XX!                       Often this =nvar, but may be < nvar if we can derive a field
!XX!                       from other fields. or set by some other assumption.
!XX      ,readT        &!  time level to be read from file
!XX      ,useIndex     &!  index in irecPrev
!XX      ,dumZ          !  z level
!XX
!XX LOGICAL ::                    &
!XX        checkNames,checkPos    &!  work
!XX       ,errFound               &!  flag indicating an error
!XX       ,readFile               &!  flag indicating if another file is to be read
!XX       ,summary                 !  work
!XX
!XX  CHARACTER(len=formatLen) ::  &
!XX       fileFormat               !  format of file
!XX
!XX  CHARACTER(len=150) ::     &
!XX       fileName             !  the name of a file
!XX
!XX!-------------------------------------------------------------------------------
!XX! Local array variables
!XX!-------------------------------------------------------------------------------
!XX  INTEGER ::                   &
!XX       ival(nvarMax)           &!  work
!XX      ,varFlag(nvarMax)        &!  flag indicating source of data.
!XX!                       At present this must either be zero (variable not used)
!XX!                       or >0 (location of first level of each variable in input file, in terms of
!XX!                       field number, or number of xy planes. Set to 1 for SDF files).
!XX      ,varInSort(nvarMax)      &!  a list of which variables in master list are to be read in,
!XX!                                  sorted into ascending order of location in file.
!XX!                                  Values 1:varIn give a list of variables that are to be
!XX!                                  read from file (in terms of location in masterlist of all possible
!XX!                                  variables), sorted into ascending order of location
!XX!                                  of data in file.
!XX      ,varFlagTmp(nvarMax)     &!  work: flags as read in
!XX      ,varStashCode(nvarMax)    !  STASH code for each variable
!XX
!XX  REAL :: tmpval(land_pts)      !  workspace, large enough for all levels of a variable
!XX  REAL :: varConst(nvarMax)     !  work: value of each variable used if a spatially constant
!XX!                                    value is to be used
!XX  REAL :: varConstTmp(nvarMax)  !  work
!XX
!XX  LOGICAL ::               &!  arrays
!XX       done(nvarMax)       &!  work
!XX       ,useVar(nvarMax)    &!  flag indicating if a variable is required for particular model setup
!XX       ,foundVar(nvarMax)   !  flag indicating if a variable is listed in run control file
!XX
!XX  CHARACTER(len=100) ::       &
!XX       varDesc(nvarMax)       &!  description of each variable
!XX       ,varName(nvarMax)      &!  names of all possible variables - used to identify
!XX!                                    variables in this routine
!XX       ,varNameTmp(nvarMax)   &!  work: variable names as read in
!XX       ,varNameSDF(nvarMax)   &!  names of variables in input file (SDF=self-describing file)
!XX       ,varNameSDFTmp(nvarMax) !  work: names as read in

!###############################################################################

! If data are not needed, nothing to do.
  IF ( .NOT. l_pdm ) RETURN

!-------------------------------------------------------------------------------

  WRITE(*,"(50('-'),/,a)") 'init_pdm'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_pdm','>INIT_PDM' )

!-------------------------------------------------------------------------------

! Read scalar parameters.
  READ(jinUnit,*) dz_pdm
  READ(jinUnit,*) b_pdm

!XX! Establish where data are to be read from.
!XX  READ(jinUnit,*) readFile
!XX  READ(jinUnit,*) fileFormat
!XX  READ(jinUnit,*) fileName
!XX
!XX! Set file format if reading from run control file.
!XX  IF ( .NOT. readFile ) fileFormat = formatAsc
!XX
!XX!-------------------------------------------------------------------------------
!XX! Set up list of all possible variables.
!XX  i = 0
!XX
!XX  i=i+1; iposFexp=i; varName(i)='fexp'
!XX  varStashcode(i) = 0   !  STASH code- not known!
!XX  varDesc(i)='Decay factor for saturated conductivity'
!XX
!XX  i=i+1; iposTi_mean=i; varName(i)='ti_mean'
!XX  varStashcode(i) = 0   !  STASH code- not known!
!XX  varDesc(i)='Mean value of topographic index'
!XX
!XX  i=i+1; iposTi_sig=i; varName(i)='ti_sig'
!XX  varStashcode(i) = 0   !  STASH code- not known!
!XX  varDesc(i)='Standard deviation of topographic index'
!XX
!XX  i=i+1; iposTi_skew=i; varName(i)='ti_skew'
!XX  varStashcode(i) = 0   !  STASH code- not known!
!XX  varDesc(i)='Skew of topographic indexti_skew'
!XX
!XX! At present, all variables are required.
!XX! In future we may select what variables are required.
!XX  useVar(:) = .TRUE.
!XX! ti_skew not currently allowed.
!XX  useVar(iposTi_skew) = .FALSE.
!XX
!XX! Count number of variables required for chosen configuration.
!XX  nvar = COUNT( useVar(:) )
!XX
!XX! Initialise flag as FALSE (variable not found in input list).
!XX  foundVar(:) = .FALSE.
!XX!-------------------------------------------------------------------------------
!XX! Only read parameters for the file format indicated.
!XX
!XX  SELECT CASE ( fileFormat )
!XX
!XX    CASE ( formatAsc,formatBin,formatPP )
!XX!     Locate the information in run control file.
!XX      CALL findTag( jinUnit,'init_pdm',tagAscBin,preInit=.TRUE. )
!XX      READ(jinUnit,*) nheaderFile,nheaderField
!XX
!XX!     Read variable names, flags and constVals from a blank-delimited list.
!XX      CALL read_list( jinUnit,3,nvarMax,'>VARS','>ENDVARS',' ','init_pdm'  &
!XX                       ,nvarIn,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2  &
!XX                       ,rvar=varConstTmp,rvarPos=3 )
!XX
!XX    CASE ( formatNc )
!XX!     Locate the information in run control file.
!XX      CALL findTag( jinUnit,'init_pdm',tagNc,preInit=.TRUE. )
!XX
!XX!     Read variable names, flags, constVals and names in file, from blank-delimited list.
!XX      CALL read_list( jinUnit,4,nvarMax,'>VARS','>ENDVARS',' ','init_pdm'  &
!XX                     ,nvarIn,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2  &
!XX                     ,rvar=varConstTmp,rvarPos=3,cvar2=varNameSDFTmp,cvar2Pos=4 )
!XX
!XX    CASE default
!XX      WRITE(*,*)'ERROR: init_pdm: no code for fileFormat=',TRIM(fileFormat)
!XX      STOP
!XX
!XX  END SELECT
!XX
!XX
!XX!-------------------------------------------------------------------------------
!XX
!XX! Check that there are no repeated names in input list.
!XX  DO ivar=1,nvarIn
!XX    DO jvar=ivar+1,nvarIn
!XX      IF ( varNameTmp(ivar) == varNameTmp(jvar) ) THEN
!XX        WRITE(*,*)'ERROR: init_pdm: repeated varNameTmp: ',TRIM(varNameTmp(ivar))
!XX        STOP'init_pdm'
!XX      ENDIF
!XX    ENDDO
!XX  ENDDO
!XX
!XX! Work out what variables are indicated.
!XX  DO ivar=1,nvarIn
!XX    SELECT CASE ( varNameTmp(ivar) )
!XX      CASE ( 'fexp' )
!XX        ipos = iposFexp
!XX      CASE ( 'ti_mean' )
!XX        ipos = iposTi_mean
!XX      CASE ( 'ti_sig' )
!XX        ipos = iposTi_sig
!XX      CASE ( 'ti_skew' )
!XX        ipos = iposTi_skew
!XX      CASE default
!XX        WRITE(*,*)'ERROR: init_pdm: do not recognise variable=',TRIM(varNameTmp(ivar))
!XX        STOP
!XX    END SELECT
!XX
!XX!   Save information for this variable in "master" list.
!XX    foundVar(ipos) = .TRUE.
!XX    varFlag(ipos) = varFlagTmp(ivar)
!XX    varConst(ipos) = varConstTmp(ivar)
!XX    varNameSDF(ipos) = varnameSDFTmp(ivar)
!XX
!XX  ENDDO
!XX
!XX! Check that all necessary variables are provided.
!XX! If a variable is provided but is not needed, reset flag.
!XX  errFound = .FALSE.
!XX  DO ivar=1,nvarMax
!XX    IF ( useVar(ivar) ) THEN
!XX
!XX      IF ( .NOT. foundVar(ivar) ) THEN
!XX        errFound = .TRUE.
!XX        WRITE(*,*)'ERROR: init_pdm: variable not supplied: ',TRIM(varName(ivar))
!XX        WRITE(*,*) TRIM(varDesc(ivar))
!XX      ELSE
!XX
!XX!       Variable is in the list.
!XX!       Determine source of data.
!XX
!XX        SELECT CASE ( varFlag(ivar) )
!XX
!XX          CASE ( 1: )
!XX!           A location is given. Nothing more to do.
!XX          CASE ( 0 )
!XX!           Variable is in the list, but we don't know how to get it.
!XX            errFound = .TRUE.
!XX            WRITE(*,*)'ERROR: init_pdm: varFlag=0'
!XX            WRITE(*,*)'Variable=',TRIM(varName(ivar))
!XX            WRITE(*,*)'Don''t know where to get data.'
!XX          CASE ( -1 )
!XX!           A constant value will be used. Nothing more to do.
!XX          CASE default
!XX!           Not allowed at present, but left in case of future use.
!XX            WRITE(*,*)'ERROR: init_pdm: varFlag<-1. Must be >= -1.'
!XX            WRITE(*,*)'Variable=',TRIM(varName(ivar))
!XX            errFound = .TRUE.
!XX
!XX        END SELECT
!XX
!XX      ENDIF  !  allVarPos
!XX
!XX    ELSE     !     NOT useVar
!XX
!XX      IF ( foundVar(ivar) ) THEN
!XX        WRITE(*,*) TRIM(varName(ivar)),' is provided but is not needed. Ignoring.'
!XX        foundVar(ivar) = .FALSE.
!XX      ENDIF
!XX
!XX    ENDIF   !  useVar
!XX  ENDDO
!XX
!XX! Stop if an error was detected.
!XX  IF ( errFound ) STOP
!XX
!XX! Calculate number of variables to be read in (not including spatially constant).
!XX  nvarIn = COUNT( varFlag(:) > 0 & varUse(:) )
!XX
!XX!-------------------------------------------------------------------------------
!XX
!XX! If no fields are to be read from external file, reset flags.
!XX  IF ( nvarIn==0 .AND. readFile ) THEN
!XX    readFile = .FALSE.
!XX    fileFormat = formatAsc
!XX    WRITE(*,*)'WARNING: init_pdm: No fields to be read from external file, but'
!XX    WRITE(*,*)'readFile=TRUE.'
!XX    WRITE(*,*)'All values will be set via flag<0 (e.g. spatially constant).'
!XX    WRITE(*,*)'Setting readFile to F.'
!XX  ENDIF
!XX
!XX! Get a sorted list, giving input variables according to position in input file(s).
!XX! This is needed for reading from run control file (and makes ASCII reading easier)
!XX! and is always done.
!XX! Note that we now insist that variables from run control file are listed in the
!XX! correct order, but this sorting code has been left - easier for now!
!XX  done(:) = .TRUE.
!XX  varInSort(:) = 0
!XX! Only consider variables that are to be read in (flag>0). SDF variables have flag set >0.
!XX  WHERE ( varFlag(:) > 0 ) done(:)=.FALSE.
!XX  DO ivar=1,nvarIn
!XX    ival(1:1) = MINLOC( varFlag(:), .NOT.done(:) )
!XX    i = ival(1)
!XX    done(i) = .TRUE.
!XX    varInSort(ivar) = i
!XX  ENDDO
!XX
!XX!-------------------------------------------------------------------------------
!XX! Some output to screen.
!XX  IF ( echo ) THEN
!XX    DO ivar=1,nvarMax
!XX      IF ( useVar(ivar) ) CALL varInfo( varName(ivar)  &
!XX                      ,varFlag(ivar),varStashcode(ivar),varNameSDF(ivar)  &
!XX                      ,fileFormat,constVal=varConst(ivar),varDesc=varDesc(ivar) )
!XX    ENDDO
!XX  ENDIF
!XX!-------------------------------------------------------------------------------
!XX
!XX! Check we have an acceptable file format, and set flags indicating what other
!XX! variables to check.
!XX  checkNames = .FALSE.
!XX  checkPos = .FALSE.
!XX  SELECT CASE ( fileFormat )
!XX    CASE ( formatAsc )
!XX!     ASCII. Need positions of variables.
!XX      checkPos = .TRUE.
!XX    CASE ( formatBin )
!XX!     GrADS. For now, need positions of variables.
!XX      checkPos = .TRUE.
!XX    CASE ( formatNc )
!XX!     netCDF. Need names of variables.
!XX      checkNames = .TRUE.
!XX    CASE ( formatPP )
!XX!     PP format. For now, need positions of variables.
!XX      checkPos = .TRUE.
!XX    CASE default
!XX      WRITE(*,*)'ERROR: init_pdm: no code for fileFormat=',TRIM(fileFormat)
!XX      STOP'init_pdm'
!XX  END SELECT
!XX
!XX  IF ( checkPos ) THEN
!XX!   Check positions of variables within file.
!XX!   Check for repeated locations.
!XX!
!XX!   First argument to checkVarPos indicates how important the order of the variables is.
!XX!   0 means order is not important. If reading from run control file, insist that variables
!XX!   are 1,2,3,.... - indicated by iorder=2. Note that for run control file, we interpret
!XX!   varFlag as variable number, not field number, since this is required for current code
!XX!   to pass checkVarPos with iorder=2 (want vars to appear adjacent, which wouldn't be if
!XX!   had multi-level data and gave non-consecutive field numbers. Admittedly, this is
!XX!   becoming a bit of a faff, and needs to be rewritten, one day.).
!XX!   In all cases, values <1 are ignored.
!XX    iorder = 0
!XX    IF ( .NOT. readFile ) iorder = 2
!XX
!XX!   We need to load file locations in sorted order.
!XX    ival(:) = 0
!XX    DO ivar=1,nvarIn
!XX      ival(ivar) = varFlag( varInSort(ivar) )
!XX    ENDDO
!XX
!XX    ierr = checkVarPos( iorder,ival(:),' init_pdm: ival' )
!XX    IF ( ierr < 0 ) THEN
!XX      WRITE(*,*)'ERROR: init_pdm: error from checkvarPos.'
!XX      WRITE(*,*)'If error was repeated use of same varPos, but you do in fact want to reuse the'
!XX      WRITE(*,*)'same data for more than one variable, comment out this stop!'
!XX      STOP
!XX    ENDIF
!XX  ENDIF
!XX
!XX  IF ( checknames ) THEN
!XX!   Check names only for variables that are to be read in (varFlag>0).
!XX    DO ivar=1,nvarMax
!XX      DO jvar=ivar+1,nvarMax
!XX        IF ( ( varFlag(ivar)>0 .AND. varFlag(jvar)>0 ) .AND.  &
!XX             ( varNameSDF(ivar)==varNameSDF(jvar) ) ) THEN
!XX          WRITE(*,*)'ERROR: init_pdm: repeated varNameSDF: ',TRIM(varNameSDF(ivar))
!XX          WRITE(*,*)'If you really do want to use the same variable from file to set values'
!XX          WRITE(*,*)'of more than one FORTRAN variable, comment out this check!'
!XX          STOP'init_pdm'
!XX        ENDIF
!XX      ENDDO
!XX    ENDDO
!XX  ENDIF
!XX!------------------------------------------------------------------------------
!XX
!XX! Open file.
!XX  IF ( readFile ) THEN
!XX    inUnit = fileUnit( fileFormat )  !  get unit
!XX    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_pdm',ncType )
!XX  ELSE
!XX    WRITE(*,*)'Reading PDM characteristics from the run control file.'
!XX    WRITE(*,*)'Data must be in fields 1,2,3,.... (i.e. sequential from 1).'
!XX    inUnit = jinUnit
!XX    CALL findTag( inUnit,'init_pdm char','>DATA',preInit=.TRUE. )
!XX  ENDIF
!XX
!XX!-------------------------------------------------------------------------------
!XX
!XX! Read data.
!XX! All levels of a variable are read into temporary space, then loaded into final variables.
!XX! Only data for the required points are returned.
!XX! Only need to do anything for variables that are to be read in, since currently all (required)
!XX! variables must be read in (no option to set as constant).
!XX! For non-SDF files, variables are read in the order in which they appear in file - this is
!XX! done so that values can be read from the run control file (which is open on jinUnit and
!XX! therefore can't be backspaced).
!XX
!XX  DO ivar=1,nvarIn
!XX
!XX!   Get location in master list of next variable to be read.
!XX    jvar = varInSort(ivar)
!XX
!XX    IF ( useVar(jvar) .AND. varFlag(jvar)>0 ) THEN
!XX
!XX      IF ( inUnit==jinUnit .AND. nxIn*nyIn==1 .AND. fileFormat==formatAsc) THEN
!XX!       If the input grid is a single point and reading from run control file, expect no
!XX!       new line between fields (eg all on one line).
!XX!       Calling readVar means we could cope with headers in the run control file.
!XX!       But there's no need since annotation is already simple in this case.
!XX        READ(jinUnit,*) tmpval(:)
!XX
!XX      ELSE
!XX
!XX!-------------------------------------------------------------------------------
!XX! Simplifying assumptions regarding input file. Otherwise have to read these in.
!XX!-------------------------------------------------------------------------------
!XX         readT      = 1         !  Time level to read from file.
!XX         dumZ       = 1         !  'z' level to read from file.
!XX!         nfieldFile = fieldNum  !  # of fields in file. Setting to field needed is OK while readT=1.
!XX         nfieldFile = nvar
!XX         nheaderT   = 0         !  No headers at top of each time.
!XX         nlineField = 0         !  Will not attempt to read ASCII line-by-line.

!XX!        Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!XX!        but need to keep index within bounds.
!XX         useIndex = inUnit
!XX         IF ( fileFormat == formatNC ) useIndex = 1

!XX
!XX!        Read data.
!XX         CALL readVar2dComp( readT,dumZ,varFlag(jvar),varStashcode(jvar),irecPrev(inUnit)  &
!XX              ,nfieldFile,nheaderFile,nheaderT  &
!XX              ,nheaderField,nlineField,nxIn,nyIn  &
!XX              ,inUnit,varNameSDF(jvar)  &
!XX              ,mapInLand(:),(/(i,i=1,land_pts)/)  &
!XX              ,fileFormat,tmpval(:),'init_pdm','init_pdm',ncType )
!XX
!XX      ENDIF
!XX
!XX!     Copy data into final variable.
!XX      IF ( jvar == iposFexp ) THEN
!XX        fexp(:) = tmpval(:)
!XX      ELSEIF ( jvar == iposTi_mean ) THEN
!XX        ti_mean(:) = tmpval(:)
!XX      ELSEIF ( jvar == iposTi_sig ) THEN
!XX        ti_sig(:) = tmpval(:)
!XX! Ti_skew not in code yet.
!XX!      ELSEIF ( jvar == iposTi_skew ) THEN
!XX!        ti_skew(:) = tmpval(:)
!XX      ELSE
!XX        WRITE(*,*)'ERROR: init_pdm: no code to load jvar=',jvar
!XX        STOP
!XX      ENDIF
!XX
!XX    ENDIF   !  read variable
!XX
!XX
!XX  ENDDO  !  variables
!XX
!XX! Close file.
!XX  CALL closeFile( inUnit,fileFormat )
!XX
!XX!-------------------------------------------------------------------------------
!XX! Now deal with variables that are to be set as spatially constant.
!XX! This involves a list of variables that is similar to that used above - it would
!XX! be much better to avoid this replication by dealing with all variables inside
!XX! one loop, regard less of whether a constant value is to be used.
!XX! But for now it's easier to do this....
!XX  DO ivar=1,nvarMax
!XX    IF ( varFlag(ivar) == -1 ) THEN
!XX      IF ( ivar == iposFexp ) THEN
!XX        fexp(:) = varConst(ivar)
!XX      ELSEIF ( ivar == iposTi_mean ) THEN
!XX        ti_mean(:) = varConst(ivar)
!XX      ELSEIF ( ivar == iposTi_sig ) THEN
!XX        ti_sig(:) = varConst(ivar)
!XX! Ti_skew not in code yet.
!XX      ELSEIF ( ivar == iposTi_skew ) THEN
!XX!        ti_skew(:) = varConst(ivar)
!XX      ELSE
!XX        WRITE(*,*)'ERROR: init_pdm: no code to load ivar=',ivar
!XX        STOP
!XX      ENDIF
!XX    ENDIF
!XX  ENDDO
!XX!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Optional writing of fields to screen.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN

!   Scalar parameters.
    write(*,*)'dz_pdm=',dz_pdm
    write(*,*)'b_pdm=',b_pdm

!XX    DO i=1,2
!XX!     i = 1 writes full field
!XX!     i = 2 writes summary statistics for whole field
!XX!     Don't bother summarizing a single point, single level field (i.e. 1 value).
!XX      if ( nxIn*nyIn==1 .AND. i==2 ) exit
!XX      IF ( i == 1 ) THEN
!XX        summary = .FALSE.
!XX      ELSE
!XX        summary = .TRUE.
!XX        WRITE(*,*)'### NB The ranges below include any ice points. ###'
!XX      ENDIF
!XX
!XX      DO ivar=1,nvarMax
!XX
!XX        IF ( useVar(ivar) ) THEN
!XX
!XX          IF ( ivar == iposFexp ) THEN
!XX            CALL varValue( summary,fexp,varFormat='f6.1',varName=varName(ivar) )
!XX
!XX          ELSEIF ( ivar == iposTi_mean ) THEN
!XX            CALL varValue( summary,ti_mean,varFormat='f6.2',varName=varName(ivar) )
!XX
!XX          ELSEIF ( ivar == iposTi_sig ) THEN
!XX            CALL varValue( summary,ti_sig,varFormat='f6.2',varName=varName(ivar) )
!XX
!XX!          ELSEIF ( ivar == iposTi_skew ) THEN
!XX!            CALL varValue( summary,ti_skew,varFormat='f6.2',varName=varName(ivar) )
!XX
!XX          ELSE
!XX            WRITE(*,*)'varName=',TRIM( varName(ivar) )
!XX            WRITE(*,*)'INIT_PDM: no code to echo this variable to screen.'
!XX            WRITE(*,*)'Not important - this is just diagnostic output.'
!XX
!XX          ENDIF
!XX
!XX        ENDIF  !  useVar
!XX
!XX      ENDDO  !  ivar
!XX    ENDDO  !  i
!XX
  ENDIF  !  echo

  END SUBROUTINE init_pdm
!###############################################################################
!###############################################################################
!###############################################################################
