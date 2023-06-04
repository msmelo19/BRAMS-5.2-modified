!###############################################################################
!###############################################################################
!###############################################################################
! subroutine init_veg
! Driver to read parameters for vegetation surface types (PFTs).

  SUBROUTINE init_veg()

  USE switches, ONLY :&
!  imported scalars with intent(in)
     routeOnly

  IMPLICIT NONE

!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

! Read values for variables that are a function of PFT only.
  CALL init_veg_pft()

! Read values for prescribed variables that may be functions of time and/or space.
  CALL init_veg_vary()

  END SUBROUTINE init_veg
!###############################################################################
!###############################################################################
!###############################################################################

! subroutine init_veg_pft
! Read vegetation parameter values that are fn(PFT) only.
! These may later be replaced by values that vary with time or location.

  SUBROUTINE init_veg_pft

  USE c_z0h_z0m, ONLY :  &
!  imported scalar parameters
     z0h_z0m

  USE file_utils, ONLY:  &
!  imported procedures
     closeFile,fileUnit,findTag,openFile

  USE inout, ONLY :   &
!  imported scalar parameters
     formatAsc,formatLen,jinUnit  &
!  imported scalars with intent(in)
    ,echo  &
!  imported scalars with intent(out)
    ,npftInFile  &
!  imported arrays with intent(out)
    ,pftUse

  USE misc_utils, ONLY :  &
!  imported procedures
     repeatVal

  USE nstypes, ONLY :  &
!  imported scalar parameters
     npft

  USE pftparm, ONLY :   &
!  imported arrays with intent(in)  &
     pftname  &
!  imported arrays with intent(out)
    ,albsnc_max,albsnc_min,albsnf_max,alpha,alnir,alpar,a_wl,a_ws,b_wl,c3,catch0  &
    ,dcatch_dlai,dz0v_dh,dgl_dm,dgl_dt,dqcrit,eta_sl  &
    ,fd,fsmc_of,f0,glmin,g_leaf_0,infil_f,kext,kpar,neff &
    ,nl0,nr_nl,ns_nl,omega,omnir,orient   &
    ,r_grow,rootd_ft,sigl,tleaf_of,tlow,tupp,emis_pft,fl_o3_ct,dfp_dcuo

  USE prognostics, ONLY :   &
!  imported arrays with intent(out)
     canht_ft,lai

  USE snow_param, ONLY :  &
!  imported arrays with intent(out)
     canSnowTile

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     can_model,l_aggregate

  USE ancil_info,  ONLY : land_pts  !DSM

!-------------------------------------------------------------------------------

  INTEGER, PARAMETER ::  &!  local SCALARS
    nvarMax=45  !  maximum possible number of parameter fields to read.
!                    In fact all fields must now be read.
!                    This does NOT include the name of PFT.

  INTEGER ::  &!  local SCALARS
    i            &!  loop counter
   ,ierr         &!  error flag
   ,ierrSum      &!  error flag
   ,inUnit       &!  unit used to connect to input file
   ,j            &!  loop counter
   ,n            &!  counter
   ,nvar          !  number of parameter fields to be read in

  INTEGER ::  &!  local ARRAYS
    snowCanPFT(npft)     ! flag indicating if canopy snow model is to be
!                            used for each PFT
!                            1=use, else=don't use

  REAL, ALLOCATABLE ::  &!  local ARRAYS
    tmpval(:,:)          !  work

  LOGICAL ::  &!  local SCALARS
    errFlag   &!  flag for error
   ,readFile  &!  flag indicating if another file is to be read
   ,rowWise    !  T if data for each type are arranged across rows
!              !  F if data for each type are arranged down a column

  LOGICAL ::  &!  local arrays
    foundPFT(npft)   !  work

  CHARACTER(len=formatLen) ::  &! local SCALARS
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &! local SCALARS
    fileName   !  the name of a file

  CHARACTER(len=LEN(pftName)), ALLOCATABLE ::  &!  local arrays
    tmpName(:)   !  work - used to hold name of PFT
!-------------------------------------------------------------------------------
  if (echo) WRITE(*,"(50('-'),/,a)") 'init_veg_pft'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_veg_pft','>INIT_VEG_PFT' )

! Establish where data are to be read from.
  READ(jinUnit,*) readFile
  READ(jinUnit,*) fileName

! Establish how many surface types are in the file.
  READ(jinUnit,*) npftInFile

!-------------------------------------------------------------------------------

! Ensure that these values are reasonable.
  IF ( npftInFile < npft ) THEN
    WRITE(*,*)'Need npft=',npft,' have only ',npftInFile
    WRITE(*,*)'Stopping in init_veg_pft'
    STOP
  ENDIF

!-----------------------------------------------------------------------

! Open file.
  fileFormat = formatAsc
  IF ( readFile ) THEN
!   Get unit
    inUnit = fileUnit( fileFormat )
!   First arg is not needed for formatted file, so using 1.
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old' )
  ELSE
    IF ( echo ) WRITE(*,*)'Reading parameters for veg types from the run control file.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_veg','>DATA',preInit=.TRUE. )
  ENDIF

! Decide how many parameter fields are to be read in. Now there's no choice.
  nvar = nvarMax   !   excluding name of PFT

! I couldn't decide whether data were best arranged in rows, or columns, so I've coded both but
! am only using one for now.
! Allocate space for all types, including those not ultimately kept.
  ierrSum = 0
  ALLOCATE( tmpval(nvar,npftInFile), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( tmpName(npftInFile), stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) THEN
    WRITE(*,*)'ERROR: init_veg_pft: could not allocate memory.'
    WRITE(*,*)'Stopping in init_veg_pft.'
    STOP
  ENDIF

! For now, assume rowwise or columnwise data.
  rowWise = .FALSE.

! Read the data.
  IF ( rowWise ) THEN
!   Data are arranged as all variables for type 1, all variables for type 2, etc...
    DO i=1,npftInFile
      READ(inUnit,*) tmpName(i),tmpval(:,i)
    ENDDO
  ELSE
!   Data are arranged as variable 1 for all types, variable 2 for all types, etc...
    READ(inUnit,*) tmpName(:)
    DO i=1,nvar
      READ(inUnit,*) tmpval(i,:)
    ENDDO
  ENDIF

! Check that names of PFTs are unique.
  IF ( repeatVal( tmpName(:) ) ) THEN
    WRITE(*,*)'ERROR: init_veg_pft: repeated name of PFT.'
    WRITE(*,*)'All PFTs must have unique names.'
    IF ( readFile ) THEN
      WRITE(*,*)'Error in file ',TRIM(fileName)
    ELSE
      WRITE(*,*)'Error reading from standard input.'
    ENDIF
    STOP
  ENDIF

! Extract the data for the required types.
  foundPFT(:) = .FALSE.
  DO j=1,npftInFile
    DO i=1,npft
      IF ( .NOT.foundPFT(i) .AND. pftName(i)==tmpName(j) ) THEN
        foundPFT(i) = .TRUE.
        pftUse(i) = j
!       Use the counter n to say which bit of tmpval to use. makes it easier to alter order!
        n = 0
        n=n+1; c3(i) = NINT( tmpval(n,j) )
        n=n+1; canht_ft(:,i) = tmpval(n,j)
        n=n+1; lai(:,i) = tmpval(n,j)
        n=n+1; catch0(i) = tmpval(n,j)
        n=n+1; dcatch_dlai(i) = tmpval(n,j)
        !DSM n=n+1; dz0v_dh(i) = tmpval(n,j)
        n=n+1; dz0v_dh(:,i) = tmpval(n,j)
        n=n+1; z0h_z0m(i) = tmpval(n,j)
        n=n+1; infil_f(i) = tmpval(n,j)
        n=n+1; rootd_ft(i) = tmpval(n,j)
        n=n+1; snowCanPFT(i) = NINT( tmpval(n,j) )
        n=n+1; albsnc_max(i) = tmpval(n,j)
        n=n+1; albsnc_min(i) = tmpval(n,j)
        n=n+1; albsnf_max(i) = tmpval(n,j)
        n=n+1; kext(i) = tmpval(n,j)
        n=n+1; kpar(i) = tmpval(n,j)
        n=n+1; orient(i) = NINT( tmpval(n,j) )
        n=n+1; alpha(i) = tmpval(n,j)
        n=n+1; alnir(i) = tmpval(n,j)
        n=n+1; alpar(i) = tmpval(n,j)
        n=n+1; omega(i) = tmpval(n,j)
        n=n+1; omnir(i) = tmpval(n,j)
        n=n+1; a_wl(i) = tmpval(n,j)
        n=n+1; a_ws(i) = tmpval(n,j)
        n=n+1; b_wl(i) = tmpval(n,j)
        n=n+1; eta_sl(i) = tmpval(n,j)
        n=n+1; g_leaf_0(i) = tmpval(n,j)
        n=n+1; dgl_dm(i) = tmpval(n,j)
        n=n+1; dgl_dt(i) = tmpval(n,j)
        n=n+1; glmin(i) = tmpval(n,j)
        n=n+1; dqcrit(i) = tmpval(n,j)
        n=n+1; fd(i) = tmpval(n,j)
        n=n+1; f0(i) = tmpval(n,j)
        n=n+1; fsmc_of(i) = tmpval(n,j)
        n=n+1; neff(i) = tmpval(n,j)
        n=n+1; nl0(i) = tmpval(n,j)
        n=n+1; nr_nl(i) = tmpval(n,j)
        n=n+1; ns_nl(i) = tmpval(n,j)
        n=n+1; r_grow(i) = tmpval(n,j)
        n=n+1; sigl(i) = tmpval(n,j)
        n=n+1; tleaf_of(i) = tmpval(n,j)
        n=n+1; tlow(i) = tmpval(n,j)
        n=n+1; tupp(i) = tmpval(n,j)
        n=n+1; emis_pft(i) = tmpval(n,j)
        n=n+1; fl_o3_ct(i) = tmpval(n,j)
        n=n+1; dfp_dcuo(i) = tmpval(n,j)
      ENDIF
    ENDDO  !  i  (pft)
  ENDDO   !  j (pft in file )

  IF ( ANY( .NOT.foundPFT(:) ) ) THEN
    WRITE(*,*)'ERROR: init_veg_pft: could not find one or more PFTs.'
    WRITE(*,*)'Could not find ',COUNT( .NOT.foundPFT(:) ),' of ',npft,' PFTs.'
    WRITE(*,*)'Could not find PFTs called:'
    DO i=1,npft
      IF ( .NOT. foundPFT(i) ) WRITE(*,*) TRIM(pftName(i))
    ENDDO
    WRITE(*,*)'NB Names are sensitive to case!'
    IF ( .NOT. readFile ) THEN
      WRITE(*,*)'Reading parameters for veg types from the run control file.'
    ELSE
      WRITE(*,*)'Reading parameters from ',TRIM(fileName)
    ENDIF
    STOP
  ENDIF

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

  DEALLOCATE( tmpName,tmpval, stat=ierr )
  IF ( ierr /= 0 ) WRITE(*,*) 'WARNING: init_veg_pft: could not deallocate memory.'
!------------------------------------------------------------------------------

! For now, don't use the location of PFTs in this file to indicate their
! location in all other files. Instead, assume that all other files will
! present the PFTs in the correct order for this run, and will have only npft
! PFTs represented.
  npftInFile = npft
  pftUse(:) = (/ (i, i=1,npft) /)
!-------------------------------------------------------------------------------

! Set switch for canopy snow model. This can be TRUE only at PFT tiles.
  canSnowTile(:) = .FALSE.
  IF ( .NOT. l_aggregate .AND. can_model==4 ) THEN
    DO i=1,npft
      IF ( snowCanPft(i) == 1 ) canSnowTile(i) = .TRUE.
    ENDDO
  ENDIF

! Check that glmin is >0.
! This ensures that wt_ext in subroutine soil_evap cannot become a NaN (which
! it would if gs=glmin and gsoil=0), or blow up, and might well be required
! elsewhere too.
  errFlag = .FALSE.
  DO i=1,npft
    IF ( glmin(i) < 1.0e-10 ) THEN
      errFlag = .TRUE.
      glmin(i) = 1.0e-10
    ENDIF
  ENDDO
  IF (errFlag ) THEN
     WRITE(*,*)'WARNING: one or more values of glmin were increased.'
     WRITE(*,*)'Very small values can cause model to blow up or NaNs.'
  ENDIF

!------------------------------------------------------------------------------
! Write to screen.
  IF ( echo ) THEN
    WRITE(*,*)'Depending upon options, some of these parameters may not be used.'
    WRITE(*,*)'canht_ft and lai may be overwritten in next section.'
!    write(*,"(a,10(tr1,a))") 'pftName=',(/ (trim(pftName(i)), i=1,npft ) /)  !  not working!
    WRITE(*,"(a,10(tr1,a))") 'pftName=',pftName(:)
    WRITE(*,*) 'c3=',c3(:)
    !DSM <vide Joe> WRITE(*,*) 'canht_ft=',canht_ft(1,:)
    !DSM <vide Joe> WRITE(*,*) 'lai=',lai(1,:)
    WRITE(*,*) 'catch0=',catch0(:)
    WRITE(*,*) 'dcatch_dlai=',dcatch_dlai(:)
    !DSM WRITE(*,*) 'dz0v_dh=',dz0v_dh(:)
    WRITE(*,*) 'dz0v_dh=',dz0v_dh(1,:)
    WRITE(*,*) 'z0h_z0m=',z0h_z0m(1:npft)
    WRITE(*,*) 'infil_f=',infil_f(:)
    WRITE(*,*) 'rootd_ft=',rootd_ft(:)
    WRITE(*,*) 'snowCanPFT=',snowCanPFT(:)
    WRITE(*,*) 'albsnc_max=',albsnc_max(:)
    WRITE(*,*) 'albsnc_min=',albsnc_min(:)
    WRITE(*,*) 'albsnf_max=',albsnf_max(:)
    WRITE(*,*) 'kext=',kext(:)
    WRITE(*,*) 'kpar=',kpar(:)
    WRITE(*,*) 'orient=',orient(:)
    WRITE(*,*) 'alpha=',alpha(:)
    WRITE(*,*) 'alnir=',alnir(:)
    WRITE(*,*) 'alpar=',alpar(:)
    WRITE(*,*) 'omega=',omega(:)
    WRITE(*,*) 'omnir=',omnir(:)
    WRITE(*,*) 'a_wl=',a_wl(:)
    WRITE(*,*) 'a_ws=',a_ws(:)
    WRITE(*,*) 'b_wl=',b_wl(:)
    WRITE(*,*) 'eta_sl=',eta_sl(:)
    WRITE(*,*) 'g_leaf_0=',g_leaf_0(:)
    WRITE(*,*) 'dgl_dm=',dgl_dm(:)
    WRITE(*,*) 'dgl_dt=',dgl_dt(:)
    WRITE(*,*) 'glmin=',glmin(:)
    WRITE(*,*) 'dqcrit=',dqcrit(:)
    WRITE(*,*) 'fd=',fd(:)
    WRITE(*,*) 'f0=',f0(:)
    WRITE(*,*) 'fsmc_of=',fsmc_of(:)
    WRITE(*,*) 'neff=',neff(:)
    WRITE(*,*) 'nl0=',nl0(:)
    WRITE(*,*) 'nr_nl=',nr_nl(:)
    WRITE(*,*) 'ns_nl=',ns_nl(:)
    WRITE(*,*) 'r_grow=',r_grow(:)
    WRITE(*,*) 'sigl=',sigl(:)
    WRITE(*,*) 'tleaf_of=',tleaf_of(:)
    WRITE(*,*) 'tlow=',tlow(:)
    WRITE(*,*) 'tupp=',tupp(:)
    WRITE(*,*) 'emis_pft=',emis_pft(:)
    WRITE(*,*) 'fl_o3_ct=',fl_o3_ct(:)
    WRITE(*,*) 'dfp_dcuo=',dfp_dcuo(:)
  ENDIF

  END SUBROUTINE init_veg_pft
!###############################################################################
!###############################################################################
!###############################################################################
! subroutine init_veg_vary
! Read details of prescribed vegetation variables that are a function of space and/or time.

  SUBROUTINE init_veg_vary

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE file_utils, ONLY:  &
!  imported procedures
    closeFile,fileUnit,findTag,openFile

  USE inout, ONLY :   &
!  imported scalar parameters
    formatAsc,formatBin,formatLen,formatNc,formatPP,npftInFile  &
   ,periodAnn,periodMon,periodOneFile,jinUnit  &
   ,tagAscBin,tagNc  &
!  imported scalars with intent(in)
   ,echo,nxIn,nyIn

  USE misc_utils, ONLY:  &
!  imported procedures
     check_template

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     npft

  USE pftparm, ONLY :   &
!  imported arrays with intent(in)  &
     pftname

  USE prognostics, ONLY :  &
!   imported arrays with intent(out)
     canht_ft,lai

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     spinUp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
    l_360,l_phenol,l_veg_compete

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
    timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    chhmmss_to_s,dateToBits,timeDate,timeDate_cmp,timeDate_diff

  USE timeConst, ONLY : &
    iSecInDay

  USE update_mod, ONLY :  &
!  imported procedures
    data_init,calc_reset_step,veg_update

  USE veg_io_vars, ONLY :  &
!  imported scalar parameters
    nvegVarMax,vegVarFlagPT,vegVarFlagPTX,vegVarFlagPX  &
!  imported scalars with intent(out)
   ,nfieldVegFile,notNextVeg,nvegDataTime,nvegFileTime,nvegFileVar  &
   ,nvegHeaderField,nvegHeaderFile,nvegHeaderTime,nvegVar  &
   ,noNewLineVeg,varNumCanht,varNumLAI,varNumRootd,vegClim  &
   ,vegDataPer,vegDataStep,vegDataStepInit,vegDataStepMax,vegDate,vegDateInit  &
   ,vegEndTime,vegFile,vegFilePer,vegFileStep,vegFileTemplate,vegFormat  &
   ,vegResetStep,vegResetStepPrev,vegTemplateDate,vegTemplateT,vegTemplateTime  &
   ,vegtemplateUnits,vegTemplateV,vegTime  &
   ,vegTimeInit,vegUpdatePer,vegUpdateStep,vegUpdateStepMax,vegUpdateStepInit,vegVaryT  &
!  imported arrays with intent(out)
   ,vegDataIn,vegFileDate,vegFileName  &
   ,vegFileTime,vegVarFlag,vegVarInterp  &
   ,vegVarName,vegVarNameFile,vegTimeIndex,vegUnit,vegUnitFile  &
   ,vegVarPos,vegVarStash
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &! local SCALARS
    i,ierr,ierrSum,iunit,ivar,j,jvar  &!  work
   ,nsec,nmin,nhr,ndy,nmon,nyr    &!  work
   ,t,tmpMonth,tmpYear     &!  work
   ,vegDataDay         !  day of month to which veg data refer to (monthly data only)

  INTEGER ::  &!  local ARRAYS
    ival(1)                   &!  work
   ,tmpPos(nvegVarMax)         !  work: temporary version of vegVarPos

  LOGICAL ::  &!  local SCALARS
    haveAve     &!  T means that one or more variable in input data is a time-average
!                     that places an extra restriction on timestep
   ,readList    &!  T means read a list of file names, F means read a single file name
   ,templateT,templateV   !  work

  LOGICAL ::  &! local ARRAYS
    done(nvegVarMax)  &!  work
   ,needTime(-1:2)     !  flag indicating which times of data are required
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  CHARACTER(len=8)  ::  &!  local SCALARS
    cvegTime   !  time string

  CHARACTER(len=LEN(vegVarName)) ::  &!  local ARRAYS
    vegName(nvegVarMax)     &!  work: names of veg variables
!                                     Only used to identify variables in this routine
   ,tmpName(nvegVarMax)     &!  work: temporary version of vegName
   ,tmpNameFile(nvegVarMax) &!  work: temporary version of vegVarNameFile
   ,tmpVegName(nvegVarMax)   !  work: temporary version of vegName

  CHARACTER(len=LEN(vegVarFlag)) ::  &!  local arrays
    tmpFlag(nvegVarMax)     !  work: temporary version of vegVarFlag

  CHARACTER(len=LEN(vegVarInterp)) ::  &!  local ARRAYS
    tmpInterp(nvegVarMax)         !  work: temporary version of vegVarInterp

  CHARACTER(len=150) ::  &!  local SCALARS
    line     &!  work
   ,listFile  !  work
!-------------------------------------------------------------------------------

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_veg_vary'

!-------------------------------------------------------------------------------
! Initialise.
!-------------------------------------------------------------------------------
  vegVaryT = .FALSE.           !   no time-varying veg
  vegTemplateT = .FALSE.       !   no time-templated file names
  vegTemplateV = .FALSE.       !   no var-templated file names
! Initialise varNumXX as showing that variables are not provided.
! Doing here so values can always be used, even if not later reset (when a variable is used).
  varNumCanht = -1
  varNumLAI = -1
  varNumRootd = -1

! Other initialisation.
  nfieldVegFile = 1  !  untested, aimed at netCDF - for f(pft,t) only I think we need nfieldVegFile>1 for tmpMap in veg_read!
! The following values are used for vegVaryT=FALSE (i.e. veg does not vary with time),
! and serve as initial values in all cases. (In fact some/all are later reset for vegVaryT=F.)
  vegFile = 1           !  use the first (only) file.
  vegFileStep = 1       !  use the first (only) time in file.
  vegDataStep = 0       !  vegDataIn will be allocated for times 0:0
  vegDataStepMax = -1   !  <=vegDataStep so first increment triggers read of data for vegVaryT=F
  vegUpdateStep = 0
  vegUpdateStepMax = -1   !  <=vegUpdateStep so we trigger an update for vegVaryT=F when vegUpdateStep is incremented

! Veg fields are prognostic with competing vegetation.
  IF ( l_veg_compete ) RETURN

!-------------------------------------------------------------------------------
! Read information from run control file.
!-------------------------------------------------------------------------------
! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_veg_vary','>INIT_VEG_VARY' )

! Read the number of time/space-varying veg variables.
  READ(jinUnit,*) nvegVar

  IF ( nvegVar > 0 ) THEN
    IF ( echo ) WRITE(*,*)'There are ',nvegVar,' prescribed veg fields that vary with space/time.'
  ELSE
    IF ( echo ) WRITE(*,*)'There are no prescribed veg fields that vary with space/time.'
    RETURN
  ENDIF

! Read details of data.
  READ(jinUnit,*) vegDataPer,vegUpdatePer
  READ(jinUnit,*) nvegFileTime,vegFilePer
  READ(jinUnit,*) vegClim
! Establish if a list of files is to be read.
  READ(jinUnit,*) readList

!-------------------------------------------------------------------------------

! Note that if vegVaryT=FALSE (set later), we may still have nvegFileTime>1 at this point.

! Allocate space for some of the veg variables.
  CALL allocate_arrays( 'init_veg_vary 1' )

! If more than one file is indicated, insist that a list is read.
  IF ( nvegFileTime>1 .AND. .NOT.readList ) THEN
    WRITE(*,*)'ERROR: init_veg: names of more than one file must be read from a list file.'
    WRITE(*,*)'Change readList to TRUE.'
    STOP
  ENDIF

! Read the name of file listing file names, or the name of the only file.
  IF ( readList ) THEN
    READ(jinUnit,*) listFile
  ELSE
    READ(jinUnit,*) vegFileName(1)
  ENDIF

  READ(jinUnit,*) vegFileDate(1),cvegTime

! Read naming convention and format of file.
  READ(jinUnit,*) vegEndTime
  READ(jinUnit,*) vegFormat

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.
!-------------------------------------------------------------------------------

  SELECT CASE ( vegFormat )

    CASE ( formatAsc,formatBin,formatPP )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_veg_vary',tagAscBin,preInit=.TRUE. )
      READ(jinUnit,*) nfieldVegFile
      READ(jinUnit,*) nvegHeaderFile,nvegHeaderTime,nvegHeaderField
      READ(jinUnit,*) noNewLineVeg
      DO ivar=1,nvegVar
        READ(jinUnit,*) tmpVegName(ivar),tmpFlag(ivar),tmpPos(ivar),tmpInterp(ivar),tmpNameFile(ivar)
      ENDDO

    CASE ( formatNc )
!     Locate the information in run control file.
      CALL findTag( jinUnit,'init_veg_vary',tagNc,preInit=.TRUE. )

      DO ivar=1,nvegVar
        READ(jinUnit,*) tmpVegName(ivar),tmpFlag(ivar),tmpInterp(ivar),tmpName(ivar),tmpNameFile(ivar)
      ENDDO

    CASE default
      WRITE(*,*)'ERROR: init_veg_vary: no code for vegFormat=',TRIM(vegFormat)
      STOP

  END SELECT

!-------------------------------------------------------------------------------
! If necessary, read list of file names and times.
!-------------------------------------------------------------------------------
  IF ( readList ) THEN
!   Read details from another file.
    IF ( echo ) WRITE(*,*)'Reading names of veg files from file.'
    iunit = fileUnit( formatAsc )   !  get unit
    CALL openFile( 1,.FALSE.,iunit,'read',formatAsc,listFile,'old' )
    READ(iunit,*)   !  header
    DO i=1,nvegFileTime
      READ(iunit,*,iostat=ierr) vegFileName(i),vegFileDate(i),cvegTime
      IF ( ierr /= 0 ) THEN
        IF ( ierr < 0 ) WRITE(*,*)'End of file before found ',nvegFileTime,' names.'
        WRITE(*,*)'ERROR: error reading list of file names.'
        WRITE(*,*)'Stopping in init_veg_vary'
        STOP
      ENDIF
      vegFileTime(i) = chhmmss_to_s( cvegTime,'init_veg_vary' )
    ENDDO
    !DSM CALL closeFile( iunit,formatAsc )
    if (iunit /= jinUnit) CALL closeFile( iunit,formatAsc )     !DSM
  ELSE
!   Have already read name and time, now convert time to integer.
    vegFileTime(1) = chhmmss_to_s( cvegTime,'init_veg_vary' )
  ENDIF  !  readList

!-------------------------------------------------------------------------------
! This routine does not read anything more from the run control file.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Check that flag values are recognised, and establish if there are time-varying fields.
!-------------------------------------------------------------------------------
  DO ivar=1,nvegVar
    SELECT CASE ( tmpFlag(ivar) )
      CASE ( vegVarFlagPT, vegVarFlagPTX )
        vegVaryT = .TRUE.
      CASE ( vegVarFlagPX )
!       Nothing to do.
      CASE default
        WRITE(*,*)'ERROR: init_veg_vary: do not recognise value of tmpFlag.'
        WRITE(*,*)'ivar=',ivar,' tmpFlag=',tmpFlag(ivar)
        STOP
    END SELECT
  ENDDO

!-------------------------------------------------------------------------------
! For now, insist that all fields have the same dimensions, e.g. all f(PFT,t) or
! all f(PFT,t,x).
!-------------------------------------------------------------------------------
  DO ivar=1,nvegVar
    DO jvar=ivar+1,nvegVar
      IF ( tmpFlag(ivar) /= tmpFlag(jvar) ) THEN
        WRITE(*,*)'ERROR: init_veg_vary: veg fields do not have consistent number of varying dimensions.'
        WRITE(*,*)'e.g. not all f(PFT,t) or f(PFT,t,x)'
        WRITE(*,*)'Have ',tmpFlag(ivar),' and ',tmpFlag(jvar)
        STOP
      ENDIF
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Don't allow fields to be functions of space if there's only one point in the input grid.
! (This also allows later test for f(x) to guarantee that # points>1).
!-------------------------------------------------------------------------------
  IF ( nxIn*nyIn == 1 ) THEN
    IF ( ANY( tmpFlag(1:nvegVar) == vegVarFlagPX ) .OR.  &
         ANY( tmpFlag(1:nvegVar) == vegVarFlagPTX ) ) THEN
      WRITE(*,*)'ERROR: init_veg_vary: veg fields cannot be regarded as varying'
      WRITE(*,*)'with position if there is only one point in the input grid.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
! Process the variable positions.
! Sort the list of input variables according to position in input file(s).
! This is needed for file formats for which the order matters (or at least can be made easier).
!-------------------------------------------------------------------------------
  done(:) = .TRUE.
  done(1:nvegVar) = .FALSE.
  vegVarPos(:) = 0
  IF ( vegFormat == formatNc ) tmpPos(:) = 1
  DO ivar=1,nvegvar
    ival(1:1) = MINLOC( tmpPos(1:nvegvar), .NOT.done(1:nvegvar) )
    i = ival(1)
    done(i) = .TRUE.
    vegName(ivar) = tmpVegName(i)
    vegVarFlag(ivar) = tmpFlag(i)
    vegVarInterp(ivar) = tmpInterp(i)
    SELECT CASE ( vegFormat )
      CASE ( formatAsc,formatBin )
        vegVarPos(ivar) = tmpPos(i)
        vegVarNameFile(ivar) = tmpNameFile(i)
      CASE ( formatNc )
        vegVarName(ivar) = tmpName(i)
        vegVarNameFile(ivar) = tmpNameFile(i)
    END SELECT
  ENDDO

!-------------------------------------------------------------------------------
! Establish if file names contain templating, and are consistent in this.
!-------------------------------------------------------------------------------
  DO i=1,nvegFileTime

    CALL check_template( vegDataPer,vegFilePer,NINT(timeStep),vegFileTime(i)  &
                        ,vegFileDate(i),vegClim,.FALSE.,vegFileName(i),'init_veg'  &
                        ,templateT,templateV,vegTemplateUnits )

    IF ( templateT ) vegTemplateT = .TRUE.
    IF ( templateV ) vegTemplateV = .TRUE.

!   Stop if the template does not match our expectations.

    IF ( vegTemplateT .NEQV. templateT ) THEN
      WRITE(*,*)'ERROR: init_veg: vegTemplateT .NEQV. templateT'
      WRITE(*,*)'File names are inconsistent.'
      IF ( templateT ) THEN
        WRITE(*,*)'File name includes time template strings, but previous&
                    & files did not'
      ELSE
        WRITE(*,*)'File name does not have time template strings, but&
                    & previous files did.'
      ENDIF
      WRITE(*,*)'i=',i,' vegFileName=',TRIM(vegFilename(i))
      STOP
    ENDIF

    IF ( vegTemplateV .NEQV. templateV ) THEN
      WRITE(*,*)'ERROR: init_veg: vegTemplateV .NEQV. templateV'
      WRITE(*,*)'File names are inconsistent.'
      IF ( templateV ) THEN
        WRITE(*,*)'File name includes variable name template strings, but previous&
                    & files did not'
      ELSE
        WRITE(*,*)'File name does not have variable name template strings, but&
                    & previous files did.'
      ENDIF
      WRITE(*,*)'i=',i,' vegFileName=',TRIM(vegFilename(i))
      STOP
    ENDIF

  ENDDO  !  files

!-------------------------------------------------------------------------------
! Insist that only one file name was given if time-templating is to be used.
!-------------------------------------------------------------------------------
  IF ( vegTemplateT .AND. nvegFileTime>1 ) THEN
    WRITE(*,*)'ERROR: init_veg_vary: vegTemplateT .AND. nvegFileTime>1'
    WRITE(*,*)'File name includes time template strings, but nvegFileTime>1.'
    WRITE(*,*)'To use time templating, set nvegFileTime=1.'
    STOP
  ENDIF


!-------------------------------------------------------------------------------
! Check time-related variables.
!-------------------------------------------------------------------------------
  IF ( vegVaryT ) THEN

!   Check that data period is >0 or monthly.
    IF ( vegDataPer<=0 .AND. vegDataPer/=periodMon ) THEN
       WRITE(*,*)'ERROR: period of veg data must be >0 OR monthly'
       WRITE(*,*)'Stopping in init_veg_vary'
       STOP
    ENDIF

    IF ( vegDataPer > 0 ) THEN
!     Check that data period is a multiple of timestep.
      IF ( MOD( vegDataPer,NINT(timeStep) ) /= 0 ) THEN
         WRITE(*,*)'ERROR: period of veg data must be a multiple of model timestep'
         WRITE(*,*)'Stopping in init_veg_vary'
         STOP
      ENDIF

!     Check that veg data are in phase with days (data refer to same times each day).
      IF ( vegDataPer<=iSecInDay ) THEN
        IF ( MOD(iSecInDay,vegDataPer) /= 0 ) THEN
          WRITE(*,*)'ERROR: veg data must be in phase with days.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
      ELSE
        IF ( MOD(vegDataPer,iSecInDay) /= 0 ) THEN
          WRITE(*,*)'ERROR: veg data must be in phase with days.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
      ENDIF

!     For climatological files, the data period must fit into one year.
!     Periods <=1 day are already known to fit into one day.
!     All "special" cases currently fit too.
!     This only leaves periods>1day (which are known to be multiples of 1 day), which
!     generally do not fit into all years - unless using 360-day calendar.
      IF ( vegClim .AND. ( vegDataPer > iSecInDay ) ) THEN
        IF ( .NOT.l_360 .OR. MOD(360, vegDataPer/iSecInDay) /= 0 ) THEN
          WRITE(*,*)'ERROR: init_veg_vary: climatological data must have a period that is a'
          WRITE(*,*)'factor of 1 year. This restriction is not important for runs of less than one year.'
          STOP
        ENDIF
      ENDIF

!     Convert data period from seconds to timesteps.
      vegDataPer = vegDataPer / NINT(timeStep)
    ENDIF  !  vegDataPer

!-------------------------

!   Change per=0 to equal timestep length.
    IF ( vegUpdatePer == 0 ) vegUpdatePer=NINT(timeStep)

!   Check that update period is >0 or monthly.
    IF ( vegUpdatePer<=0 .AND. vegUpdatePer/=periodMon ) THEN
      WRITE(*,*)'ERROR: period for update of veg data must be >0 OR monthly'
      WRITE(*,*)'Stopping in init_veg_vary'
      STOP
    ENDIF

!   Check that update period fits other constraints.
    IF ( vegUpdatePer > 0 ) THEN
!     Check that update period is a multiple of timestep.
      IF ( MOD( vegUpdatePer,NINT(timeStep) ) /= 0 ) THEN
        WRITE(*,*)'ERROR: period for update of veg data must be a multiple of model timestep'
        WRITE(*,*)'Stopping in init_veg_vary'
        STOP
      ENDIF

!     Check that update veg data are in phase with days (data refer to same
!     times each day), and with the data frequency.
      IF ( vegUpdatePer <= iSecInDay ) THEN
        IF ( MOD(iSecInDay,vegUpdatePer) /= 0 ) THEN
          WRITE(*,*)'ERROR: veg update must be in phase with days.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
      ELSE
        IF ( MOD(vegUpdatePer,iSecInDay) /= 0 ) THEN
          WRITE(*,*)'ERROR: veg update must be in phase with days.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
      ENDIF

!     At present vegDataPer is always >0, but start to future-proof by testing anyway.
      IF ( vegDataPer>0 ) THEN
        IF ( MOD(vegDataPer*NINT(timeStep),vegUpdatePer) /= 0 ) THEN
           WRITE(*,*)'ERROR: vegUpdatePer must be a factor of vegDataPer.'
           WRITE(*,*)'Stopping in init_veg_vary'
           STOP
        ENDIF
      ENDIF

!     Convert update period from seconds to timesteps.
      vegUpdatePer = vegUpdatePer / NINT(timeStep)

    ELSE
!     vegUpdatePer <= 0
!     Check that data frequency is not greater than that of update (only monthly
!     option at present - so we cannot read more than once a month).
      IF ( vegDataPer >= 0 ) THEN
        WRITE(*,*)'ERROR: monthly update of veg fields but data are more frequent (vegDataPer).'
        WRITE(*,*)'Stopping in init_veg_vary'
        STOP
      ENDIF
      IF ( vegDataPer /= periodMon ) THEN
        WRITE(*,*)'ERROR: init_veg_vary: monthly updates can only be used with monthly data.'
        STOP
      ENDIF

    ENDIF  !  vegUpdatePer

!-------------------------------------------------------------------------------
!   Extra checks for time templating.
!-------------------------------------------------------------------------------

    IF ( vegTemplateT ) THEN

!-------------------------------------------------------------------------------
!     Check that the file period is acceptable.
!     File period is only used for time-templated files.
!-------------------------------------------------------------------------------
      IF ( vegFilePer < 0 ) THEN
!       Check that the period is recognised.
!       Note that we don't check if the file period is a multiple of data period.
        SELECT CASE ( vegFilePer )
          CASE ( periodAnn, periodMon )
!           OK, nothing to do.
          CASE ( periodOneFile )
!           All data in one file - this makes no sense for time templating.
            WRITE(*,*)'ERROR: init_veg_vary: vegFilePer=',vegFilePer
            WRITE(*,*)'All data in one file does  not make sense for time templating.'
            STOP
          CASE default
            WRITE(*,*)'ERROR: init_veg_vary: do not recognise vegFilePer=',vegFilePer
            STOP
        END SELECT
      ELSE
!       vegFilePer >= 0.
!       Check file period is a multiple of data period (which is already known to be
!       a multiple of timestep length, and to be a factor of one day).
        IF ( MOD( vegFilePer,vegDataPer*NINT(timeStep) ) /= 0 ) THEN
          WRITE(*,*)'ERROR: file period must be a multiple of data period.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
!       Further restriction for time template - insist that a single file cannot
!       hold data for more than 24 hours - makes life easier.
        IF ( vegFilePer>iSecInDay ) THEN
          WRITE(*,*)'ERROR: init_veg: files too infrequent.'
          WRITE(*,*)'For these template files, a single file cannot &
                     &hold data for more than 24 hours.'
          WRITE(*,*)'This restriction does not apply to files with a &
                 &special" period (e.g. monthly).'
          STOP
        ENDIF
!       Convert file period to number of timesteps.
        vegFilePer = vegFilePer / NINT(timeStep)
      ENDIF  !  vegFilePer
    ENDIF  !  vegTemplateT

!-------------------------------------------------------------------------------
!   If monthly data, establish what day of month is used.
!-------------------------------------------------------------------------------
    IF ( vegDataPer == periodMon ) THEN
      CALL dateToBits( vegFileDate(1),vegDataDay,tmpMonth,tmpYear,l_360  &
                      ,'init_veg_vary: month date' )
      IF ( .NOT. l_360 ) THEN
!       If day of month > 28, this can be difficult - so avoid.
        IF ( vegDataDay > 28 ) THEN
          WRITE(*,*)'ERROR: init_veg_vary: vegDataDay > 28.'
          WRITE(*,*)'The date of the first veg file has day of month > 28.'
          WRITE(*,*)'This is difficult, because not all months have this number of days!'
          WRITE(*,*)'To use this date, see notes in code at this point.'
          WRITE(*,*)'If data are timestamped for end of month, use 00H on'
          WRITE(*,*)'1st of next month.'
          STOP
        ENDIF
      ELSE   !  l_360
        IF ( vegDataDay > 30 ) THEN
          WRITE(*,*)'ERROR: init_veg_vary: vegDataDay > 30.'
          WRITE(*,*)'The date of the first veg file has day of month > 30.'
          WRITE(*,*)'This cannot be used with 360-day calendar.'
          STOP
        ENDIF
      ENDIF  !  l_360
!     Notes on how to use day>28: More code is required! We could assume that
!     day>28indicates "use last day of month" or, less usefully, we
!     could assume that it meant "use this number of days in from end of month"(!?!).
!     Then set vegDataDay to <1, and whenever vegDataDay is used, would have to test
!     for values <1 and put in extra code. All seems bit of a faff for a fairly
!     unlikely occurence.
      WRITE(*,*)'Day of month on which monthly veg data are valid on/timestamped=' &
                 ,vegDataDay
    ENDIF

  ELSE

!-------------------------------------------------------------------------------
!   Not time-varying veg.
!   Ignore/reset options that are not needed in this case.
!-------------------------------------------------------------------------------
!   Insist that there are no time headers.
    IF ( nvegHeaderTime /= 0 ) WRITE(*,*)'WARNING: setting nvegHeaderTime to zero, since veg fields are not a function of time'
    nvegHeaderTime = 0
!   Only need 1 "time" of files.
    IF ( nvegFileTime /= 1 ) THEN
!     Code will work, but insist that user sets nvegFileTime=1.
      WRITE(*,*)'ERROR: init_veg_vary: NOT vegVaryT AND nvegFileTime/=1'
      WRITE(*,*)'There are no time-varying fields, so set nvegFileTime=1.'
      STOP
    ENDIF

  ENDIF  !  vegVaryT
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------

  IF ( vegtemplateT ) THEN
!   Save some values.
    vegFileTemplate = vegFileName(1)
    vegTemplateDate = vegFileDate(1)
    vegTemplateTime = vegFileTime(1)
  ENDIF

!-------------------------------------------------------------------------------
! Check that file name or template, is (vaguely) acceptable.
! Note that this means that a file name for use without templating cannot
! contain the character "%".
!-------------------------------------------------------------------------------

  DO i=1,nvegFileTime
    CALL check_template( vegDataPer,vegFilePer,NINT(timeStep),vegFileTime(i)  &
                        ,vegFileDate(i),vegClim,.TRUE.,vegFileName(i),'init_veg'  &
                        ,templateT,templateV,vegTemplateUnits )
  ENDDO

!-------------------------------------------------------------------------------
  IF ( .NOT. vegTemplateT ) THEN
!-------------------------------------------------------------------------------
!   Set the time/date for the last file (one more than we actually have) to be
!   100 years after start of penultimate file - easiest option for now (just
!   don't have hundreds of years of data....).
!   If we have climatological files, set date to one year after first file. If
!   this is before some of the other files listed, an error will be raised (correctly).
!-------------------------------------------------------------------------------
    i = nvegFileTime+1
    j = i - 1   !  reference file
    vegFileName(i) = 'dummy_file'
    nmon = 1200   !  10 years
    IF ( vegClim ) THEN
      j = 1
      nmon = 12
    ENDIF
    CALL timeDate( vegFileTime(j),vegFileDate(j),nmon,'mon',l_360  &
                    ,vegFileTime(i),vegFileDate(i),'init_veg_vary_time')

!-------------------------------------------------------------------------------
!   Check that the files are listed in chronological order.
!-------------------------------------------------------------------------------
    DO i=2,nvegFileTime
      j = i - 1
      IF ( .NOT. timeDate_cmp( vegFileTime(i),vegFileDate(i),'>'  &
                      ,vegFileTime(j),vegFileDate(j),'init_veg_vary' )  ) THEN
        WRITE(*,*)'ERROR: init_veg_vary: files are not listed in chronological order.'
        WRITE(*,*)'Problem files and times:'
        WRITE(*,*) j,TRIM(vegFileName(j)),vegFileDate(j),vegFileTime(j)
        WRITE(*,*) i,TRIM(vegFileName(i)),vegFileDate(i),vegFileTime(i)
        STOP
      ENDIF
    ENDDO

!-------------------------------------------------------------------------------
!   Check that the dates given for climatological files do not cover more than
!   one year. They can cover two calendar years.
!   At present this is a necessary but not sufficient condition that
!   things will work!
!-------------------------------------------------------------------------------
    IF ( vegClim ) THEN
      CALL timeDate_diff( vegFileTime(1),vegFileDate(1)  &
           ,vegFileTime(nvegFileTime),vegFileDate(nvegFileTime)  &
           ,l_360,'init_veg_vary'  &
           ,nsec,nmin,nhr,ndy,nmon,nyr )
      i = nsec+nmin+nhr+ndy+nmon
      IF ( nyr>1 .OR. ( nyr==1 .AND. i>0 ) ) THEN
        WRITE(*,*)'ERROR: init_veg_vary: bad dates for climatological files.'
        WRITE(*,*)'Dates of files must not cover more than one year.'
        WRITE(*,*)'They can cover two calendar years.'
        STOP
      ENDIF
    ENDIF


  ENDIF  !  vegTemplateT

!  do i=1,nvegFileTime+1
!    print*, i,trim(vegFileName(i)),' date,time=',vegFileDate(i),vegFileTime(i)
!  enddo

!-------------------------------------------------------------------------------

! Do further processing/checking.
  CALL init_veg_vary2( vegName )

!-------------------------------------------------------------------------------
! Don't allow prescribed time-varying LAI if phenology has been switched on.
! This possibly isn't the best place to be testing such things, but is done here for now!
!-------------------------------------------------------------------------------
  IF ( l_phenol ) THEN
    SELECT CASE ( vegVarFlag(varNumLAI) )
      CASE ( vegVarFlagPT, vegVarFlagPTX )
        WRITE(*,*)'ERROR: init_veg_vary: time-varying LAI should not be prescribed'
        WRITE(*,*)'when the phenology model (l_phenol) is also active.'
        STOP
    END SELECT
  ENDIF

!-------------------------------------------------------------------------------
! Sort out noNewLineVeg.
!-------------------------------------------------------------------------------
! noNewLineVeg only makes sense for ASCII files.
  IF ( vegFormat /= formatAsc ) noNewLineVeg=.FALSE.

  IF ( noNewLineVeg ) THEN
!   Only allow noNewLineVeg if either (1) there's only one point in input grid
!   or (2) variables are only functions of PFT and time (not space).
!   Previously we have tested that variables that are f(x) only come with nxIn*nyIn>1.
    IF ( ANY(vegVarFlag(:)==vegVarFlagPX) .OR. ANY(vegVarFlag(:)==vegVarFlagPTX) ) THEN
      WRITE(*,*)'ERROR: init_veg_vary: veg data on one line  (really every field'
      WRITE(*,*)'does not need to be on a separate line) is only allowed if'
      WRITE(*,*)'veg fields are not functions of position.'
      STOP
    ENDIF
    IF ( vegTemplateV ) THEN
      WRITE(*,*)'ERROR: init_veg_vary: veg data on one line  (really every field'
      WRITE(*,*)'does not need to be on a separate line) cannot be used with'
      WRITE(*,*)'variable name templating - little/no point!'
      STOP
    ENDIF
    IF ( nvegHeaderField /= 0 ) THEN
      WRITE(*,*)'ERROR: nvegHeaderField must be zero for noNewLineVeg.'
      WRITE(*,*)'Stopping in init_veg_vary'
      STOP
    ENDIF
    IF ( nfieldVegFile < npftInFile*nvegVar ) THEN
      WRITE(*,*)'ERROR: init_veg_vary: nfieldVegFile < npftInFile*nvegVar'
      WRITE(*,*)'It appears that there are not enough values in file.'
      WRITE(*,*)'nfieldVegfile=number of variables*number of PFTs in file.'
      STOP
    ENDIF
  ENDIF

!-------------------------------------------------------------------------------
  IF ( vegVaryT ) THEN

!   If time interpolation is indicated, but data period equals timestep, change interpolation flag
!   to indicate that no interpolation is required.
    IF ( vegDataPer == 1 ) THEN
      DO ivar=1,nvegvar
        SELECT CASE ( vegVarInterp(ivar) )
          CASE ( 'i' )
            WRITE(*,*)'vegDataPer=timestep, so changing flag i to nf'
            vegVarInterp(ivar) = 'nf'
          CASE ( 'b' )
            WRITE(*,*)'vegDataPer=timestep, so changing flag b to nb'
            vegVarInterp(ivar) = 'nb'
          CASE ( 'c' )
            WRITE(*,*)'vegDataPer=timestep, so changing flag c to nc'
            vegVarInterp(ivar) = 'nc'
          CASE ( 'f' )
            WRITE(*,*)'vegDataPer=timestep, so changing flag f to nf'
            vegVarInterp(ivar) = 'nf'
        END SELECT
      ENDDO  !  ivar
    ENDIF

!   Check that interpolation flag is recognised. Also look for time-averaged input, and
!   work out how many and what time levels of input data are required.
!   eg needTime(0:1)=T indicates that, when model is running, at a veg data timestep, we
!   need the current and next veg data
!   to be able to calculate/interpolate values for all model timesteps before the next veg read.
    haveAve = .FALSE.
    needTime(:) = .FALSE.
    DO ivar=1,nvegvar
      SELECT CASE ( vegVarInterp(ivar) )
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
          WRITE(*,*)'ERROR: do not recognise vegVarInterp: ',TRIM(vegVarInterp(ivar))
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
      END SELECT
    ENDDO

!   Monthly updates (which currently implies monthly data) must use interpolation
!   flags 'nb' or 'nf' - i.e. no interpolation.
    IF ( vegUpdatePer == periodMon ) THEN
      DO ivar=1,nvegvar
        SELECT CASE ( vegVarInterp(ivar) )
          CASE ( 'nb', 'nf' )
          CASE default
            WRITE(*,*)'ERROR: init_veg_vary: vegVarInterp invalid.'
            WRITE(*,*)'Monthly updates must use interpolation flags nb or nf,'
            WRITE(*,*)'i.e. no interpolation allowed - none needed.'
            STOP
        END SELECT
      ENDDO
    ENDIF

!   At present, "special" data periods cannot be used with either of the 'centred'
!   average flags - there's no code to decide when e.g. we're halfway through month.
    IF ( vegDataPer < 0 ) THEN
      DO ivar=1,nvegvar
        SELECT CASE ( vegVarInterp(ivar) )
          CASE ( 'nc', 'c' )
            WRITE(*,*)'ERROR: init_veg_vary: vegVarInterp invalid.'
            WRITE(*,*)'Data with "special" period (monthly, annual) cannot use'
            WRITE(*,*)'interpolation flags nc or c - i.e. cannot be centred averages'
            WRITE(*,*)'- unless you write some more code. Sorry.'
            STOP
        END SELECT
      ENDDO

    ENDIF

!   Check that timestep allows interpolation of averages to be conservative.
    IF ( haveAve ) THEN
      IF ( vegDataPer > 0 ) THEN
        IF ( MOD(vegDataPer,2)/=0 ) THEN
          WRITE(*,*)'ERROR: init_veg_vary: period of veg data is not an even number of timesteps.'
          WRITE(*,*)'The algorithm for interpolation of time-averaged input data requires'
          WRITE(*,*)'that the data period comprises an even number of model timesteps.'
          WRITE(*,*)'Stopping in init_veg_vary'
          STOP
        ENDIF
      ELSE
!       These special cases may give "irregular data" (eg monthly), meaning that the time
!       interval over which a datum is valid may vary between periods (e.g. Jan v. Feb).
!       Keeping track of this to ensure that we always have the correct data available and
!       that interpolation will conserve, requires more code than has yet been written.
!       So we refuse to proceed unless we are confident we have equal intervals - which is
!       the case if l_360=TRUE, since (at present) all special cases indicate periods that
!       are multiples of 1 day (eg month, year), and the model timestep is known to be a
!       factor of 1 day. A faff.
        IF ( .NOT. l_360 ) THEN
          WRITE(*,*)'ERROR: init_veg_vary: option not supported.'
          WRITE(*,*)'If the conserving interpolation is chosen, the current code cannot'
          WRITE(*,*)'deal with irregularly spaced data (e.g. monthly).'
          WRITE(*,*)'Either change interpolation flag, or write some new code which'
          WRITE(*,*)'needs to ensure conservation, and track when "future" data need'
          WRITE(*,*)'to be read.'
          STOP
        ELSE
          IF ( vegDataPer == periodMon ) THEN
!           Check that one month (30 days) consists of an even number of timesteps.
            IF ( MOD(30.0*REAL(iSecInDay),2*timeStep)/=0 ) THEN
              WRITE(*,*)'ERROR: init_veg_vary: period of veg data is not an even number of timesteps.'
              WRITE(*,*)'The algorithm for interpolation of time-averaged input data requires'
              WRITE(*,*)'that the data period comprises an even number of model timesteps.'
              WRITE(*,*)'Stopping in init_veg_vary'
              STOP
            ENDIF
          ELSE
            WRITE(*,*)'ERROR: init_veg_vary: no code for vegDatPer=',vegDataPer
            WRITE(*,*)'More code required!'
            STOP
          ENDIF
        ENDIF  !  l_360
      ENDIF    !  vegDataPer
    ENDIF      !  haveAve

!   Work out how many time levels of data to hold.
!   If a given time is required for one variable, all variables are held at that time.
    vegTimeIndex(:) = -9
!   Get earliest time required.
    DO t=-1,2
      IF ( needTime(t) ) THEN
        vegTimeIndex(1) = t
        EXIT
      ENDIF
    ENDDO
!   Get latest time required.
    DO t=2,-1,-1
      IF ( needTime(t) ) THEN
        vegTimeIndex(2) = t
        EXIT
      ENDIF
    ENDDO
    nvegDataTime = vegTimeIndex(2) - vegTimeIndex(1) + 1
    PRINT*, 'nvegDataTime=',nvegDataTime

  ELSE

!   NOT vegVaryT
    vegTimeIndex(:) = 0   !  a single time level of data will be held

  ENDIF  !  vegVaryT

!-------------------------------------------------------------------------------
! Allocate space for input data.
!-------------------------------------------------------------------------------
  IF ( echo ) WRITE(*,*)'Allocating for veg data. nvegvar=',nvegvar,' vegTimeIndex=',vegTimeIndex(:)
  CALL allocate_arrays( 'init_veg_vary 2' )

!-------------------------------------------------------------------------------
! Get file unit(s). Open a scratch file on each to reserve it. Not needed for netCDF files.
! In fact we could probably proceed without "reserving" units in this
! way...but this seems to be working for now.
!-------------------------------------------------------------------------------

! Initialise vegUnit. This value should NOT be a possible netCDF ID, so use -1.
  vegUnit(:) = -1
  IF ( vegFormat /= formatNc ) THEN
    DO i=1,nvegFileVar
      vegUnit(i) = fileUnit( vegFormat )
      vegUnitFile(i) = 'scratch'
      CALL openFile( 1,.FALSE.,vegUnit(i),'scratch','scratch',vegUnitFile(i),'scratch' )
    ENDDO
  ENDIF

!--------------------------------------------------------------------------------
! Some output to screen.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN
    WRITE(*,*) nvegVar,' veg variables vary with time or time and position.'
    IF ( vegVaryT ) THEN
      IF ( vegDataPer > 0 ) THEN
        WRITE(*,*)'Period of vegetation data = ',vegDataPer*timeStep,' seconds.'
      ELSEIF ( vegDataPer == periodMon ) THEN
        WRITE(*,*)'Period of vegetation data = 1 month.'
      ENDIF
      IF ( vegUpdatePer > 0 ) THEN
        WRITE(*,*)'Vegetation data are updated every ',vegUpdatePer*timeStep,' seconds.'
      ELSEIF ( vegUpdatePer == periodMon ) THEN
        WRITE(*,*)'Vegetation data are updated once a month.'
      ENDIF
    ENDIF
    IF ( vegTemplateT ) THEN
      WRITE(*,*)'Files are named following time-templating conventions.'
      IF ( vegEndTime ) THEN
        WRITE(*,*)'Convention uses end time of file.'
      ELSE
        WRITE(*,*)'Convention uses start time of file.'
      ENDIF
    ENDIF
    IF ( vegFormat /= formatNc ) WRITE(*,*)'There are ',nfieldVegFile,' fields in each file.'
    DO ivar=1,nvegVar
      IF ( ivar == varNumCanht ) THEN
        line = 'Canopy height'
      ELSEIF ( ivar == varNumLAI ) THEN
        line = 'LAI'
      ELSEIF ( ivar == varNumRootd ) THEN
        line = 'Rootdepth'
      ENDIF
      IF ( vegFormat == formatPP ) THEN
        line =  TRIM(line) // ' uses STASH code='
        i = LEN_TRIM(line)
        WRITE(line(i:i+3),"(i4)") vegVarStash(ivar)
      ENDIF
      IF ( vegFormat /= formatNc ) THEN
        line = TRIM(line) // ' starts at field #'
        i = LEN_TRIM(line)
        WRITE(line(i:i+2),"(i3)") vegVarPos(ivar)
      ELSE
        line = TRIM(line) // ' is called ' // TRIM(vegVarName(ivar)) // ' in input'
      ENDIF
      SELECT CASE ( vegVarFlag(ivar) )
        CASE ( vegVarFlagPT )
          line = TRIM(line) // ' and is a function of PFT and time'
        CASE ( vegVarFlagPTX )
          line = TRIM(line) // ' and is a function of PFT, time and location'
        CASE ( vegVarFlagPX )
          line = TRIM(line) // ' and is a function of PFT and location'
      END SELECT
      WRITE(*,*) TRIM(line)
      SELECT CASE ( vegVarFlag(ivar) )
        CASE ( vegVarFlagPT,vegVarFlagPTX )
          line  = '......with interpolation flag='
          i = LEN_TRIM(line)
          WRITE(line(i+1:),"(a)") vegVarInterp(ivar)
          WRITE(*,*) TRIM(line)
      END SELECT

    ENDDO  !  variables
    WRITE(*,*)'NB All other vegetation fields are functions of PFT only.'
  ENDIF  !  echo

!-------------------------------------------------------------------------------
! Initialise veg data.
!-------------------------------------------------------------------------------
  IF ( vegVaryT ) THEN

!   Read in past values that are needed on first timestep.
    CALL data_init( vegDataPer,vegFilePer,vegTemplateDate,vegTemplateTime  &
        ,vegDataStepMax,vegDataStep,vegDataStepInit  &
        ,vegDate,vegDateInit  &
        ,vegFile,vegFileStep,vegResetStep,vegResetStepPrev  &
        ,vegTime,vegTimeInit  &
        ,vegFileDate,vegFileTime,vegFileTemplate,vegTemplateUnits,vegFileName  &
        ,vegTimeIndex,vegClim,vegEndTime,templateT,notNextVeg,'veg'  &
        ,vegUpdatePer,vegUpdateStepMax,vegUpdateStep,vegUpdateStepInit )

!   Deal with the possibility that veg data may need to be "recycled" to cope with spin up.
!   Arguments: zero=timestep.
    IF ( spinUp ) CALL calc_reset_step( 0,vegResetStep,vegResetStepPrev,'veg' )

  ELSE

!   Read in a spatial field.
!   1st argument (a_step)=0, 2nd (next)=.TRUE. to read next data, 3rd should read
!   data into the only time level of vegDataIn.
    CALL veg_update( 0,.TRUE.,vegTimeIndex(1) )

  ENDIF

!-------------------------------------------------------------------------------
! If veg is not time-varying, don't need units or variables again.
! Also do optional print out.
!-------------------------------------------------------------------------------
  IF ( .NOT. vegVaryT ) THEN

!   Close files and deallocate.
    DO i=1,nvegFileVar
      !DSM CALL closeFile( vegUnit(i),vegFormat )
      if (vegUnit(i) /= jinUnit) CALL closeFile( vegUnit(i),vegFormat )  !DSM
    ENDDO
    ierrSum = 0
    DEALLOCATE( vegVarPos,vegVarFlag,vegVarInterp,vegVarName,vegVarNameFile, stat=ierr ); ierrSum=ierrSum+ierr
    DEALLOCATE( vegFileName,vegFileDate,vegFileTime, stat=ierr ); ierrSum=ierrSum+ierr
    DEALLOCATE( vegDataIn,vegUnit,vegUnitFile, stat=ierr ); ierrSum=ierrSum+ierr
    IF ( ierrSum/=0 ) WRITE(*,"(50('#'),/,a,50('#'))")'WARNING: init_veg_vary: could not deallocate'

!   Optional print to screen for fields that are not functions of time.
    IF ( echo ) THEN

!     First, print entire fields.
      DO ivar=1,nvegVar
        IF ( ivar == varNumCanht ) THEN
          DO i=1,npft
            WRITE(*,*)'Canopy height for ',TRIM(pftName(i))
            WRITE(*,*) canht_ft(:,i)
          ENDDO
        ELSEIF ( ivar == varNumLAI ) THEN
          DO i=1,npft
            WRITE(*,*)'LAI for ',TRIM(pftName(i))
            WRITE(*,*) lai(:,i)
          ENDDO
        ENDIF
      ENDDO  !  ivar

!     Now, print summary info. Note that this will include locations where frac=0.
      DO ivar=1,nvegVar
        IF ( ivar == varNumCanht ) THEN
          DO i=1,npft
            WRITE(*,"(a,2(a,f6.2))") 'Range of canopy height for '  &
               ,TRIM(pftName(i)),MINVAL(canht_ft(:,i)),' to ',MAXVAL(canht_ft(:,i))
          ENDDO
        ELSEIF ( ivar == varNumLAI ) THEN
          DO i=1,npft
            WRITE(*,"(a,2(a,f6.2))") 'Range of LAI for '  &
              ,TRIM(pftName(i)),MINVAL(lai(:,i)),' to ',MAXVAL(lai(:,i))
          ENDDO
        ENDIF
      ENDDO  !  ivar

    ENDIF   !  echo

  ENDIF   !  vegVaryT

  END SUBROUTINE init_veg_vary
!################################################################################
!################################################################################
!################################################################################
! subroutine init_veg_vary2
! Further processing of prescribed vegetation fields.

  SUBROUTINE init_veg_vary2( vegName )

  USE inout, ONLY :  &
!  imported scalar parameters
    formatAsc,formatBin,formatNc

  USE veg_io_vars, ONLY :  &
!  imported scalar parameters
     nvegVarMax,vegVarFlagPTX,vegVarFlagPX  &
!  imported scalars with intent(in)
    ,nfieldVegFile,nvegFileVar,nvegVar,nvegvar,nvegVarMax,vegFormat  &
      ,vegVarNameFile,vegTemplateV  &
      ,vegVarPos,nvegFileVar,nvegVar,nvegvar,nvegVarMax,nfieldVegFile  &
!  imported scalars with intent(out)
    ,varNumCanht,varNumLAI,varNumRootd  &
!  imported arrays with intent(in)
    ,vegVarName,vegVarNameFile  &
!  imported arrays with intent(inout)
    ,vegVarPos  &
!  imported arrays with intent(out)
    ,nvegFileVar,nvegVar,nvegvar,nvegVarMax,nfieldVegFile  &
    ,vegVarNameFile,vegVarStash

  USE misc_utils, ONLY :  &
!  imported procedures
    checkVarPos,repeatVal
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    ivar,jvar       !  work

  LOGICAL ::  &!  local SCALARS
    checkNames,checkPos  !  work

  CHARACTER(len=*), INTENT(in) ::  &!  in ARRAYS
    vegName(nvegVarMax)      !  names of veg variables

!--------------------------------------------------------------------------------

! Check that there are no repeated names.
  IF ( repeatVal( vegName(1:nvegvar) ) ) THEN
    WRITE(*,*)'ERROR: init_veg_vary2: repeated vegName'
    WRITE(*,*)'Stopping in init_veg_vary2'
    STOP
  ENDIF

! Establish what variables have been provided for input.
! Also set required PP code.
  DO ivar=1,nvegvar
    SELECT CASE ( vegName(ivar) )
      CASE ( 'canht' )
        varNumCanht = ivar
        vegVarStash(ivar) = 218  !  PP code 1393
      CASE ( 'lai' )
        varNumLAI = ivar
        vegVarStash(ivar) = 217  !  PP code 1392
      CASE ( 'rootd' )
        varNumRootd = ivar
        vegVarStash(ivar) = 8219 ! not sure this is what we want!
      CASE default
        WRITE(*,*)'ERROR: init_veg_vary2: do not recognise the variable:',TRIM(vegName(ivar))
!       Error message is only useful while it is kept up to date!
        WRITE(*,*)'Allowable names are : canht, lai, rootd.'
        STOP
    END SELECT
  ENDDO

! At this stage, varNumXX is set for all veg variables that will be read from file. For other
! variables (not read from file, and possibly not even required for this run), varNumXX retains its
! initial value (currently -1).

!--------------------------------------------------------------------------------

! Check we have an acceptable file format, and set flags indicating what other variables to check.
  checkNames = .FALSE.
  checkPos = .FALSE.
  SELECT CASE ( vegFormat )
    CASE ( formatAsc )
!     ASCII. Need positions of variables.
      checkPos = .TRUE.
    CASE ( formatBin )
!     GrADS. For now, need positions of variables.
      checkPos = .TRUE.
    CASE ( formatNc )
!     netCDF. Need names of variables.
      checkNames = .TRUE.
    CASE default
      WRITE(*,*)'ERROR: init_veg_vary2: no code for vegFormat=',TRIM(vegFormat)
      WRITE(*,*)'Stopping init_veg_vary2'
      STOP
  END SELECT

! Check positions of variables within file, unless each from separate file.
! Now only checking for doublers, and for variables beyond range in file.
! Should also check if vegTemplateV=T, as could have more than one var per
! file - but not at present (see notes later).
  IF ( checkPos .AND. .NOT.vegTemplateV ) THEN  !  XX should also check if vegTemplateV, since could have more than one var per file.
    IF ( ANY(vegVarPos(1:nvegvar) < 1) ) THEN
      WRITE(*,*)'ERROR: position for variable < 1'
      WRITE(*,*)'Stopping init_veg_vary2'
      STOP
    ENDIF
    IF ( ANY(vegVarPos(1:nvegvar) > nfieldVegFile) ) THEN
      WRITE(*,*)'ERROR: init_veg_vary2: field for variable exceeds number in file'
      WRITE(*,*)'Number in file=',nfieldVegFile,' positions=',vegVarPos(1:nvegvar)
      WRITE(*,*)'NB Each PFT is considered a separate field.'
    ENDIF
!   Check for doublers (meaning repeated values).
!   First argument to checkVarPos indicates how important the order of the variables is.
!   0 means order is not importamt
    ivar = checkVarPos( 0,vegVarPos(1:nvegvar),' init_veg_vary2: vegVarPos' )
    IF ( ivar < 0 ) THEN
      WRITE(*,*)'ERROR: init_veg_vary2: error from checkvarPos.'
      WRITE(*,*)'If error was repeated use of same varPos, but you do in fact want to reuse the'
      WRITE(*,*)'same data for more than one variable, comment out this stop!'
      STOP
    ENDIF
  ENDIF

  IF ( vegFormat/=formatNc .AND. nfieldVegFile<nvegvar ) THEN
!   Assuming that vegTemplateV mkeans all files have same number of variables!
    WRITE(*,*)'ERROR: init_veg_vary2: nfieldVegFile < nvegvar'
    WRITE(*,*)'nfieldVegFile=',nfieldVegFile,' nvegvar=',nvegvar
    WRITE(*,*)'Number of variables in a veg file is less than number of required variables.'
    STOP
  ENDIF

! Get the number of files that are to be open at any time.
! This is one, unless variable-name-templating is used, in which case a separate file is used
! for each variable (also see note below regarding relaxing this assumption).
  nvegFileVar = 1
  IF ( vegTemplateV ) nvegFileVar=nvegvar

!------------------------------------------------------------------------------
! If necessary, check that names of variables are unique.
! Also check names used in file names.
  DO ivar=1,nvegvar
    DO jvar=ivar+1,nvegvar
      IF ( checkNames .AND. ( vegVarName(ivar) == vegVarName(jvar) ) ) THEN
        WRITE(*,*)'ERROR: init_veg_vary2: repeated vegVarName: ',TRIM(vegVarName(ivar))
        WRITE(*,*)'Stopping in init_veg_vary2'
        STOP
      ENDIF
      IF ( vegTemplateV .AND. ( vegVarNameFile(ivar) == vegVarNameFile(jvar) ) ) THEN
!       This check is consistent with the current assumption that all variables are
!       either in one file (then should use vegTemplateV) or one variable per file.
!       If want to allow >1 file with >1 var per file, have to write new code
!       along the lines used for driving data (see driveVarNameUnit and the likes).
        WRITE(*,*)'ERROR: repeated vegVarNameFile: ',TRIM(vegVarNameFile(ivar))
        WRITE(*,*)'Stopping in init_veg_vary2'
        STOP
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE init_veg_vary2


!################################################################################
!################################################################################
!################################################################################

