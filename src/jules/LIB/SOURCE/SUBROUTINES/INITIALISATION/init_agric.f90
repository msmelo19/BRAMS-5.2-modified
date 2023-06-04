!###############################################################################
!###############################################################################
! Read fractional cover of agriculture.
!###############################################################################
!###############################################################################
SUBROUTINE Init_Agric()

  USE ancil_info, ONLY : &
!  imported scalars with intent(in)
     land_pts

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile,fileUnit,findTag,openFile &
!  imported scalars with intent(in)
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

  USE readwrite_mod, ONLY :  &
!  imported procedures
    readvar2dComp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_triffid,routeOnly

  USE trifctl, ONLY :  &
!  imported arrays with intent(out)
     frac_agr
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, parameter ::  &!  scalar parameters
    stashCodeAgric = 0   !  STASH code for agricultural fraction - not known!

  INTEGER ::  &!  local SCALARS
    fieldNum     &!  field number in file
   ,i            &!  loop counter
   ,inUnit       &!  unit used to connect to input file
   ,nfieldFile   &!  number of fields per time in a file
   ,nheaderField &!  number of header records before each field in file
   ,nheaderFile  &!  number of header records at start of file
   ,nheaderT     &!  number of headers at start of each time
   ,nlineField   &!  work
   ,readT        &!  time level to be read from file
   ,readZ        &!  'z' level to be read from file
   ,useIndex      !  index in irecPrev

  LOGICAL ::  &!  local SCALARS
    readFile   !  flag indicating if another file is to be read

  CHARACTER(len=formatLen) ::  &! local SCALARS
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &! local SCALARS
    fileName  &!  the name of a file
   ,varName    !  name of variable in input file

!------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! If not using TRIFFID, nothing to do.
!-------------------------------------------------------------------------------
  IF ( .NOT.l_triffid .OR. routeOnly ) RETURN

  WRITE(*,"(50('-'),/,a)") 'init_agric'

!-------------------------------------------------------------------------------
! Locate the start of this section in input file and establish where data are
! to be read from.
!-------------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_agric','>INIT_AGRIC' )

  READ(jinUnit,*) readFile

!-------------------------------------------------------------------------------
! Only read parameters for the file format indicated.
!-------------------------------------------------------------------------------
  IF ( readFile ) THEN

!   An external file will be read.
    READ(jinUnit,*) fileFormat
    READ(jinUnit,*) fileName

    SELECT CASE ( fileFormat )

      CASE ( formatAsc,formatBin,formatPP )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_agric',tagAscBin,preInit=.TRUE. )
        READ(jinUnit,*) nheaderFile,nheaderField
        READ(jinUnit,*) fieldNum

      CASE ( formatNc )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_agric',tagNc,preInit=.TRUE. )
        READ(jinUnit,*) varName

      CASE default
        WRITE(*,*)'ERROR: init_agric: no code for fileFormat=',TRIM(fileFormat)
        STOP

    END SELECT

!-------------------------------------------------------------------------------
!   Data will be read from run control file.  The first field will be
!   read, no headers expected. Field number is redundant for stdIn, but is used
!   to set nfieldFile.
!-------------------------------------------------------------------------------
  ELSE  !  NOT readFile
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
!   Use first arg to openFile to set recl (unformatted file) to be enough for a
!   single value.
    !DSM CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_agric',ncType )
  ELSE
    WRITE(*,*)'Reading fractional cover of agriculture from the run control file.'
    WRITE(*,*)'Data must be the first field encountered.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_agric','>DATA',preInit=.TRUE. )
  ENDIF

!-------------------------------------------------------------------------------
! Simplifying assumptions regarding input file. Otherwise have to read these in.
!-------------------------------------------------------------------------------
  readT      = 1         !  Time level to read from file.
  readZ      = 1         !  'z' level to read from file.
  nfieldFile = fieldNum  !  # of fields in file. Setting to field needed is OK while readT=1.
  nheaderT   = 0         !  No headers at top of each time.
  nlineField = 0         !  Will not attempt to read ASCII line-by-line.

! Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
! but need to keep index within bounds.
  useIndex = inUnit
  if ( fileFormat == formatNC ) useIndex = 1

! Read data.
  !DSM CALL readVar2dComp( readT,readZ,fieldNum,stashCodeAgric,irecPrev(useIndex)  &
  !DSM                    ,nfieldFile,nheaderFile,nheaderT  &
  !DSM                    ,nheaderField,nlineField,nxIn,nyIn,inUnit  &
  !DSM                    ,varName  &
  !DSM                    ,mapInLand(:),(/(i,i=1,land_pts)/),fileFormat,frac_agr(:)  &
  !DSM                    ,'init_agric','init_agric',ncType )

  frac_agr=0.0   !DSM Colocando provisoriamente que nao existe area de agricultura
!-------------------------------------------------------------------------------
! Close file if it not the JULES control file
!-------------------------------------------------------------------------------
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!------------------------------------------------------------------------------
! Write to stdout agricultural fraction for each gridbox.
!-------------------------------------------------------------------------------
  IF ( echo ) WRITE(*,"('frac_agr=',8f5.2,(10f5.2))") frac_agr(:)

END SUBROUTINE init_agric

