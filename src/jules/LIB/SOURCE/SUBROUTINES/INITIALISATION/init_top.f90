! Read parameters for TOPMODEL.

SUBROUTINE init_top()

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts

  USE c_topog, ONLY :   &
!  imported scalars with intent(out)
     ti_max,ti_wetl,zw_max

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

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType

  USE misc_utils, ONLY :  &
!  imported procedures
     checkVarPos,checkVars,read_list,repeatVal,varInfo,varList,varValue

  USE top_pdm, ONLY :   &
!  imported arrays with intent(out)
     fexp,ti_mean,ti_sig

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar2dcomp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_top

!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: nvarMax = 4 !  the number of possible variables to be set
!                            This may be larger than the number needed for any
!                            particular configuration.
!                            Possible variables:
!                            fexp,ti_mean,ti_sig,ti_skew

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::        &
       i,ierr       &!  loop counters/work
      ,inUnit       &!  unit used to connect to input file
      ,iorder       &!  work
!-------------------------------------------------------------------------------
!  Every possible variable (nvarMax) must have an "iposX" variables declared here
      ,iposFexp     &!  position in list of all possible variables of fexp
      ,iposTi_mean  &!  position in list of all possible variables of ti_mean
      ,iposTi_sig   &!  position in list of all possible variables of ti_sig
      ,iposTi_skew  &!  position in list of all possible variables of ti_skew
!-------------------------------------------------------------------------------
      ,ivar,jvar    &!  work
      ,nfieldFile   &!  number of fields per time in a file
      ,nheaderField &!  number of header records before each field in file
      ,nheaderFile  &!  number of header records at start of file
      ,nheaderT     &!  number of headers at start of each time
      ,nlineField   &!  work
      ,nvar      &!  number of variables required for chosen configuration
      ,nvarFound    &!  number of variables found in list
      ,nvarIn       &!  number of variables that are read in
!                       Often this =nvar, but may be < nvar if we can derive a field
!                       from other fields. or set by some other assumption.
      ,readT        &!  time level to be read from file
      ,useIndex     &!  index in irecPrev
      ,dumZ           !  z level

 LOGICAL ::                    &
       checkNames,checkPos    &!  work
      ,errFound               &!  flag indicating an error
      ,readFile               &!  flag indicating if another file is to be read
      ,sdfFile                &!  TRUE if dealing with an SDF (self-describing file)
      ,summary                 !  work

  CHARACTER(len=formatLen) ::  &
       fileFormat               !  format of file

  CHARACTER(len=150) ::     &
       fileName             !  the name of a file

!-------------------------------------------------------------------------------
! Local array variables
!-------------------------------------------------------------------------------
  INTEGER ::                   &
       ival(nvarMax)           &!  work
      ,varFlag(nvarMax)        &!  flag indicating source of data.
!                       At present this must either be zero (variable not used)
!                       or >0 (location of first level of each variable in input file, in terms of
!                       field number, or number of xy planes. Set to 1 for SDF files).
      ,varInSort(nvarMax)      &!  a list of which variables in master list are to be read in,
!                                  sorted into ascending order of location in file.
!                                  Values 1:varIn give a list of variables that are to be
!                                  read from file (in terms of location in masterlist of all possible
!                                  variables), sorted into ascending order of location
!                                  of data in file.
      ,varFlagTmp(nvarMax)     &!  work: flags as read in
      ,varStashCode(nvarMax)   &!  STASH code for each variable
      ,varUse(nvarMax)          !  flag indicating if a variable
!                        is required for current setup, and how it is used
!                        0: not required
!                        1: needed. Can be set indirectly (from other variables).
!                        2: only needed indirectly to set other variables.


  REAL :: tmpval(land_pts)      !  workspace, large enough for all levels of a variable
  REAL :: varConst(nvarMax)     !  work: value of each variable used if a spatially constant
!                                    value is to be used
  REAL :: varConstTmp(nvarMax)  !  work

  LOGICAL ::               &!  arrays
       done(nvarMax)       &!  work
      ,foundVar(nvarMax)    !  flag indicating if a variable is listed in run control file

  CHARACTER(len=100) ::       &
        varDesc(nvarMax)       &!  description of each variable
       ,varName(nvarMax)      &!  names of all possible variables - used to identify
!                                    variables in this routine
       ,varNameTmp(nvarMax)   &!  work: variable names as read in
       ,varNameSDF(nvarMax)   &!  names of variables in input file (SDF=self-describing file)
       ,varNameSDFTmp(nvarMax) !  work: names as read in

!###############################################################################

! If data are not needed, nothing to do.
  IF ( .NOT. l_top ) RETURN

!-------------------------------------------------------------------------------

  WRITE(*,"(50('-'),/,a)") 'init_top'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_top','>INIT_TOP' )

!-------------------------------------------------------------------------------

! Read scalar parameters.
  READ(jinUnit,*) zw_max
  READ(jinUnit,*) ti_max
  READ(jinUnit,*) ti_wetl

! Establish where data are to be read from.
  READ(jinUnit,*) readFile
  READ(jinUnit,*) fileFormat
  READ(jinUnit,*) fileName

! Set file format if reading from run control file.
  IF ( .NOT. readFile ) fileFormat = formatAsc

!-------------------------------------------------------------------------------
! Set up list of all possible variables.
  i = 0

  i=i+1; iposFexp=i; varName(i)='fexp'
  varStashcode(i) = 0   !  STASH code- not known!
  varDesc(i)='Decay factor for saturated conductivity'

  i=i+1; iposTi_mean=i; varName(i)='ti_mean'
  varStashcode(i) = 0   !  STASH code- not known!
  varDesc(i)='Mean value of topographic index'

  i=i+1; iposTi_sig=i; varName(i)='ti_sig'
  varStashcode(i) = 0   !  STASH code- not known!
  varDesc(i)='Standard deviation of topographic index'

  i=i+1; iposTi_skew=i; varName(i)='ti_skew'
  varStashcode(i) = 0   !  STASH code- not known!
  varDesc(i)='Skew of topographic indexti_skew'

! At present, all variables are required.
! In future we may select what variables are required.
  varUse(:) = 1
! ti_skew not currently allowed.
  varUse(iposTi_skew) = 0
! Count number of variables required for chosen configuration.
  nvar = COUNT( varUse(:) > 0 )

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.

  SELECT CASE ( fileFormat )

    CASE ( formatAsc,formatBin,formatPP )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_top',tagAscBin,preInit=.TRUE. )
      READ(jinUnit,*) nheaderFile,nheaderField

!     Read variable names, flags and constVals from a blank-delimited list.
      CALL read_list( jinUnit,3,nvarMax,'>VARS','>ENDVARS',' ','init_top'  &
                       ,nvarFound,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2  &
                       ,rvar=varConstTmp,rvarPos=3 )
      sdfFile = .FALSE.

    CASE ( formatNc )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_top',tagNc,preInit=.TRUE. )

!     Read variable names, flags, constVals and names in file, from blank-delimited list.
      CALL read_list( jinUnit,4,nvarMax,'>VARS','>ENDVARS',' ','init_top'  &
                     ,nvarFound,cvar1=varNameTmp,cvar1Pos=1,ivar=varFlagTmp,ivarPos=2  &
                     ,rvar=varConstTmp,rvarPos=3,cvar2=varNameSDFTmp,cvar2Pos=4 )
      sdfFile = .TRUE.

    CASE default
      WRITE(*,*)'ERROR: init_top: no code for fileFormat=',TRIM(fileFormat)
      STOP

  END SELECT


!-------------------------------------------------------------------------------

! Check that there are no repeated names in input list.
  IF ( repeatVal( varNameTmp(1:nvarFound) ) ) THEN
    WRITE(*,*)'ERROR: init_top: repeated variable names.'
    STOP
  ENDIF

! Work out what variables are indicated.
  CALL varList( nvarFound,sdfFile,varNameTmp,varName  &
               ,errFound,foundVar  &
               ,varFlagIn=varFlagTmp,varConstIn=varConstTmp  &
               ,varNameSDFin=varNameSDFtmp  &
               ,varFlag=varFlag,varConst=varConst  &
               ,varNameSDF=varNameSDF )
  IF ( errFound ) THEN
    WRITE(*,*)'ERROR: init_top: error returned by varList'
    STOP
  ENDIF

! Check that all necessary variables are provided.
  CALL checkVars( nvarMax,.TRUE.,.FALSE.,varUse,foundVar,varDesc,varName  &
                 ,nvarIn,varFlag,errFound )
! Stop if an error was detected.
  IF ( errFound ) THEN
    WRITE(*,*)'ERROR: init_top: an error was found by checkVars.'
    STOP
  ENDIF

! Calculate number of variables to be read in (not including spatially constant).
  nvarIn = COUNT( varFlag(:) > 0 .AND. varUse(:) > 0 )

!-------------------------------------------------------------------------------

! If no fields are to be read from external file, reset flags.
  IF ( nvarIn==0 .AND. readFile ) THEN
    readFile = .FALSE.
    fileFormat = formatAsc
    WRITE(*,*)'WARNING: init_top: No fields to be read from external file, but'
    WRITE(*,*)'readFile=TRUE.'
    WRITE(*,*)'All values will be set via flag<0 (e.g. spatially constant).'
    WRITE(*,*)'Setting readFile to F.'
  ENDIF

! Get a sorted list, giving input variables according to position in input file(s).
! This is needed for reading from run control file (and makes ASCII reading easier)
! and is always done.
! Note that we now insist that variables from run control file are listed in the
! correct order, but this sorting code has been left - easier for now!
  done(:) = .TRUE.
  varInSort(:) = 0
! Only consider variables that are to be read in (flag>0). SDF variables have flag set >0.
  WHERE ( varFlag(:) > 0 ) done(:)=.FALSE.
  DO ivar=1,nvarIn
    ival(1:1) = MINLOC( varFlag(:), .NOT.done(:) )
    i = ival(1)
    done(i) = .TRUE.
    varInSort(ivar) = i
  ENDDO

!-------------------------------------------------------------------------------
! Some output to screen.
  IF ( echo ) THEN
    DO ivar=1,nvarMax
      IF ( varUse(ivar) > 0 ) CALL varInfo( varName(ivar)  &
                      ,varFlag(ivar),varStashcode(ivar),varNameSDF(ivar)  &
                      ,fileFormat,constVal=varConst(ivar),varDesc=varDesc(ivar) )
    ENDDO
  ENDIF
!-------------------------------------------------------------------------------

! Check we have an acceptable file format, and set flags indicating what other
! variables to check.
  checkNames = .FALSE.
  checkPos = .FALSE.
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )
!     ASCII. Need positions of variables.
      checkPos = .TRUE.
    CASE ( formatBin )
!     GrADS. For now, need positions of variables.
      checkPos = .TRUE.
    CASE ( formatNc )
!     netCDF. Need names of variables.
      checkNames = .TRUE.
    CASE ( formatPP )
!     PP format. For now, need positions of variables.
      checkPos = .TRUE.
    CASE default
      WRITE(*,*)'ERROR: init_top: no code for fileFormat=',TRIM(fileFormat)
      STOP'init_top'
  END SELECT

  IF ( checkPos ) THEN
!   Check positions of variables within file.
!   Check for repeated locations.
!
!   First argument to checkVarPos indicates how important the order of the variables is.
!   0 means order is not important. If reading from run control file, insist that variables
!   are 1,2,3,.... - indicated by iorder=2. Note that for run control file, we interpret
!   varFlag as variable number, not field number, since this is required for current code
!   to pass checkVarPos with iorder=2 (want vars to appear adjacent, which wouldn't be if
!   had multi-level data and gave non-consecutive field numbers. Admittedly, this is
!   becoming a bit of a faff, and needs to be rewritten, one day.).
!   In all cases, values <1 are ignored.
    iorder = 0
    IF ( .NOT. readFile ) iorder = 2

!   We need to load file locations in sorted order.
    ival(:) = 0
    DO ivar=1,nvarIn
      ival(ivar) = varFlag( varInSort(ivar) )
    ENDDO

    ierr = checkVarPos( iorder,ival(:),' init_top: ival' )
    IF ( ierr < 0 ) THEN
      WRITE(*,*)'ERROR: init_top: error from checkvarPos.'
      WRITE(*,*)'If error was repeated use of same varPos, but you do in fact want to reuse the'
      WRITE(*,*)'same data for more than one variable, comment out this stop!'
      STOP
    ENDIF
  ENDIF

  IF ( checknames ) THEN
!   Check names only for variables that are to be read in (varFlag>0).
    DO ivar=1,nvarMax
      DO jvar=ivar+1,nvarMax
        IF ( ( varFlag(ivar)>0 .AND. varFlag(jvar)>0 ) .AND.  &
             ( varNameSDF(ivar)==varNameSDF(jvar) ) ) THEN
          WRITE(*,*)'ERROR: init_top: repeated varNameSDF: ',TRIM(varNameSDF(ivar))
          WRITE(*,*)'If you really do want to use the same variable from file to set values'
          WRITE(*,*)'of more than one FORTRAN variable, comment out this check!'
          STOP'init_top'
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!------------------------------------------------------------------------------

! Open file.
  IF ( readFile ) THEN
    inUnit = fileUnit( fileFormat )  !  get unit
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_top',ncType )
  ELSE
    WRITE(*,*)'Reading topmodel characteristics from the run control file.'
    WRITE(*,*)'Data must be in fields 1,2,3,.... (i.e. sequential from 1).'
    inUnit = jinUnit
    CALL findTag( inUnit,'init_top char','>DATA',preInit=.TRUE. )
  ENDIF

!-------------------------------------------------------------------------------

! Read data.
! All levels of a variable are read into temporary space, then loaded into final variables.
! Only data for the required points are returned.
! Only need to do anything for variables that are to be read in, since currently all (required)
! variables must be read in (no option to set as constant).
! For non-SDF files, variables are read in the order in which they appear in file - this is
! done so that values can be read from the run control file (which is possibly open on stdIn and
! therefore can't be backspaced).

  DO ivar=1,nvarIn

!   Get location in master list of next variable to be read.
    jvar = varInSort(ivar)

    IF ( varUse(jvar)>0 .AND. varFlag(jvar)>0 ) THEN

      IF ( inUnit==jinUnit .AND. nxIn*nyIn==1 .AND. fileFormat==formatAsc) THEN
!       If the input grid is a single point and reading from run control file, expect no
!       new line between fields (eg all on one line).
!       Calling readVar means we could cope with headers in the run control file.
!       But there's no need since annotation is already simple in this case.
        READ(jinUnit,*) tmpval(:)

      ELSE

!-------------------------------------------------------------------------------
! Simplifying assumptions regarding input file. Otherwise have to read these in.
!-------------------------------------------------------------------------------
         readT      = 1         !  Time level to read from file.
         dumZ       = 1         !  'z' level to read from file.
         nfieldFile = nvar      !  # of fields in file
         nheaderT   = 0         !  No headers at top of each time.
         nlineField = 0         !  Will not attempt to read ASCII line-by-line.

!        Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!        but need to keep index within bounds.
         useIndex = inUnit
         IF ( fileFormat == formatNC ) useIndex = 1

!        Read data.
         CALL readVar2dComp( readT,dumZ,varFlag(jvar),varStashcode(jvar),irecPrev(useIndex)  &
              ,nfieldFile,nheaderFile,nheaderT  &
              ,nheaderField,nlineField,nxIn,nyIn  &
              ,inUnit,varNameSDF(jvar)  &
              ,mapInLand(:),(/(i,i=1,land_pts)/)  &
              ,fileFormat,tmpval(:),'init_top','init_top',ncType )

      ENDIF

!     Copy data into final variable.
      IF ( jvar == iposFexp ) THEN
        fexp(:) = tmpval(:)
      ELSEIF ( jvar == iposTi_mean ) THEN
        ti_mean(:) = tmpval(:)
      ELSEIF ( jvar == iposTi_sig ) THEN
        ti_sig(:) = tmpval(:)
! Ti_skew not in code yet.
!      ELSEIF ( jvar == iposTi_skew ) THEN
!        ti_skew(:) = tmpval(:)
      ELSE
        WRITE(*,*)'ERROR: init_top: no code to load jvar=',jvar
        STOP
      ENDIF

    ENDIF   !  read variable


  ENDDO  !  variables

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!-------------------------------------------------------------------------------
! Now deal with variables that are to be set as spatially constant.
! This involves a list of variables that is similar to that used above - it would
! be much better to avoid this replication by dealing with all variables inside
! one loop, regard less of whether a constant value is to be used.
! But for now it's easier to do this....
  DO ivar=1,nvarMax
    IF ( varFlag(ivar) == -1 ) THEN
      IF ( ivar == iposFexp ) THEN
        fexp(:) = varConst(ivar)
      ELSEIF ( ivar == iposTi_mean ) THEN
        ti_mean(:) = varConst(ivar)
      ELSEIF ( ivar == iposTi_sig ) THEN
        ti_sig(:) = varConst(ivar)
! Ti_skew not in code yet.
      ELSEIF ( ivar == iposTi_skew ) THEN
!        ti_skew(:) = varConst(ivar)
      ELSE
        WRITE(*,*)'ERROR: init_top: no code to load ivar=',ivar
        STOP
      ENDIF
    ENDIF
  ENDDO
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Optional writing of fields to screen.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN

!   Scalar parameters.
    write(*,*)'zw_max=',zw_max
    write(*,*)'ti_max=',ti_max
    write(*,*)'ti_wetl=',ti_wetl

    DO i=1,2
!     i = 1 writes full field
!     i = 2 writes summary statistics for whole field
!     Don't bother summarizing a single point, single level field (i.e. 1 value).
      if ( nxIn*nyIn==1 .AND. i==2 ) exit
      IF ( i == 1 ) THEN
        summary = .FALSE.
      ELSE
        summary = .TRUE.
        WRITE(*,*)'### NB The ranges below include any ice points. ###'
      ENDIF

      DO ivar=1,nvarMax

        IF ( varUse(ivar) > 0 ) THEN

          IF ( ivar == iposFexp ) THEN
            CALL varValue( summary,fexp,varFormat='f6.1',varName=varName(ivar) )

          ELSEIF ( ivar == iposTi_mean ) THEN
            CALL varValue( summary,ti_mean,varFormat='f6.2',varName=varName(ivar) )

          ELSEIF ( ivar == iposTi_sig ) THEN
            CALL varValue( summary,ti_sig,varFormat='f6.2',varName=varName(ivar) )

!          ELSEIF ( ivar == iposTi_skew ) THEN
!            CALL varValue( summary,ti_skew,varFormat='f6.2',varName=varName(ivar) )

          ELSE
            WRITE(*,*)'varName=',TRIM( varName(ivar) )
            WRITE(*,*)'INIT_TOP: no code to echo this variable to screen.'
            WRITE(*,*)'Not important - this is just diagnostic output.'

          ENDIF

        ENDIF  !  varUse

      ENDDO  !  ivar
    ENDDO  !  i

  ENDIF  !  echo


!-------------------------------------------------------------------------------
! Calculate fitting parameters.
  CALL CALC_FIT_FSAT

  END SUBROUTINE init_top
!###############################################################################
!###############################################################################
!###############################################################################
