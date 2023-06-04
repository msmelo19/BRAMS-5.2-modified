!###############################################################################
!###############################################################################
  SUBROUTINE init_hgt()

! Initialise the surface height field.

!-------------------------------------------------------------------------------

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts,ntiles

  USE c_elevate, ONLY :  &
!  imported arrays with intent(out)
     surf_hgt

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
     varInfo,varValue

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar3dcomp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_aggregate,routeOnly

!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Local scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: stashCode = 0  !  STASH code for tile height field (not known)

!-------------------------------------------------------------------------------
! Local scalars.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    fieldNum     &!  field number in file of first field in frac
   ,i            &!  loop counter
   ,inUnit       &!  unit used to connect to input file
   ,nfieldFile   &!  number of fields per time in a file
   ,nheaderField &!  number of header records before each field in file
   ,nheaderFile  &!  number of header records at start of file
   ,nheaderT     &!  number of headers at start of each time
   ,nlineField   &!  work
   ,readT        &!  time level to be read from file
   ,useIndex      !  index in irecPrev

  LOGICAL ::  &
    levels       &!  work
   ,readFile     &!  flag indicating if another file is to be read
   ,summary      &!  work
   ,zeroHeight    !  flag indicating if all heights are to be set to zero

  CHARACTER(len=formatLen) ::  &
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &
    fileName    &!  the name of a file
   ,varNameSDF   !  name of variable in input netCDF file

!###############################################################################

!-------------------------------------------------------------------------------
! If data are not needed, nothing to do.
!-------------------------------------------------------------------------------
  IF ( routeOnly ) RETURN

  if (echo) WRITE(*,"(50('#'),/,a)") 'init_hgt'

!-------------------------------------------------------------------------------
! If aggregate tiles are being used, set tile heights to zero.
!-------------------------------------------------------------------------------
  IF ( l_aggregate ) THEN
    surf_hgt(:,:) = 0.0
!   Nothing more to do, so leave this routine.
    RETURN
  ENDIF

!-------------------------------------------------------------------------------

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_hgt','>INIT_HGT' )

! Read switch to set all heights to zero.
  READ(jinUnit,*) zeroHeight

! If all heights are to be zero, do this now.
  IF ( zeroHeight ) THEN
    if (echo) write(*,*)'INIT_HGT: zeroHeight=TRUE, all tile heights will be set to zero.'
    surf_hgt(:,:) = 0.0
!   Nothing more to do, so leave this routine.
    RETURN
  ENDIF

!-------------------------------------------------------------------------------
! Establish where data are to be read from, and other details of format.  Only
! read parameters for the file format indicated.
!-------------------------------------------------------------------------------
  READ(jinUnit,*) readFile

  IF ( readFile ) THEN
!-------------------------------------------------------------------------------
!   An external file will be read.
!-------------------------------------------------------------------------------
    READ(jinUnit,*) fileFormat
    READ(jinUnit,*) fileName

    SELECT CASE ( fileFormat )

      CASE ( formatAsc,formatBin,formatPP )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_hgt',tagAscBin )
        READ(jinUnit,*) nheaderFile,nheaderField
        READ(jinUnit,*) fieldNum

      CASE ( formatNc )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_grid',tagNc,preInit=.TRUE. )
        READ(jinUnit,*) varNameSDF

      CASE default
        WRITE(*,*)'ERROR: init_hgt: no code for fileFormat=',TRIM(fileFormat)
        STOP

    END SELECT

  ELSE   !  NOT readFile
!-------------------------------------------------------------------------------
!   Data will be read from run control file.
!   The first field will be read, no headers expected. Field number is
!   redundant for stdIn, but is used to set nfieldFile.
!-------------------------------------------------------------------------------
    fileFormat = formatAsc
    nheaderFile = 0
    nheaderField = 0
    fieldNum = 1

  ENDIF

!-------------------------------------------------------------------------------
! Open file.
!-------------------------------------------------------------------------------
  IF ( readFile ) THEN
    inUnit = fileUnit( fileFormat )  !  get unit
!   Use first arg to openFile to set recl (unformatted file) to be enough for a single value.
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_hgt',ncType )
!   Summarise source of data.
    IF ( echo ) CALL varInfo( 'surf_hgt',fieldNum,stashCode,varNameSDF  &
                             ,fileFormat,varDesc='Tile surface height' )
  ELSE
    WRITE(*,*)'Reading tile height data from the run control file.'
    WRITE(*,*)'Data must be the first field encountered.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_hgt','>DATA',preInit=.TRUE. )
  ENDIF

!-------------------------------------------------------------------------------
! Read data.
!-------------------------------------------------------------------------------
  IF ( inUnit==jinUnit .AND. nxIn*nyIn==1 ) THEN
!   If only need to read one space point from run control file, expect no new
!   line between fields (eg all on one line).  Calling readVar means we could
!   cope with headers in the run control file!  But there's no need since
!   annotation is already simple in this case.
    READ(jinUnit,*) surf_hgt(:,:)
  ELSE
!   Simplifying assumptions regarding input file. Otherwise have to read these in.
    readT      = 1                !  time level to read from file
    nfieldFile = fieldNum+ntiles-1 !  # of fields in file. Set to last level of required field - OK while readT=1
    nheaderT   = 0                !  no headers at top of each time
    nlineField = 0                !  0 means will not attempt to read ASCII line-by-line

!   Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!   but need to keep index within bounds.
    useIndex = inUnit
    IF ( fileFormat == formatNC ) useIndex = 1

    CALL readVar3dComp( readT,fieldNum,stashCode,.FALSE.,irecPrev(useIndex),nfieldFile,nheaderFile  &
                       ,nheaderT,nheaderField,nlineField,nxIn,nyIn,ntiles,.FALSE. &
                       ,inUnit,varNameSDF  &
                       ,mapInLand(:),(/(i,i=1,land_pts)/)  &
                       ,fileFormat,surf_hgt(:,:),'init_hgt','init_hgt',ncType )
  ENDIF

!-------------------------------------------------------------------------------
! Close file if it is not the run control file
!-------------------------------------------------------------------------------
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!-------------------------------------------------------------------------------
! Optional writing of fields to screen.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN

    DO i=1,2
!     i = 1 writes full field, level by level
!     i = 2 writes summary statistics for whole field
!     Don't bother summarizing a single point, single level field (i.e. 1 value).
      IF ( nxIn*nyIn==1 .AND. i==2 ) EXIT
      IF ( i == 1 ) THEN
        summary = .FALSE.
        levels = .TRUE.
      ELSE
        summary = .TRUE.
        levels = .FALSE.
        WRITE(*,*)'### NB The ranges below include any ice points. ###'
      ENDIF

      CALL varValue( levels,summary,surf_hgt,varFormat='f6.1',varName='surf_hgt' )

    ENDDO  !  i

  ENDIF  !  echo

  END SUBROUTINE init_hgt
!###############################################################################
!###############################################################################
!###############################################################################


