!#######################################################################
!#######################################################################
!#######################################################################
! module file_utils
!
! Contains various procedures to open a file, allocate
! a unit to a file etc. Also saves irecPrev.
!
!-------------------------------------------------------------------------------
  MODULE file_utils

  USE jules_netcdf, ONLY :  &
!   imported procedures
     check_nc_dims

  USE inout, ONLY :  &
!   imported scalar parameters
     formatAsc,formatBin,formatNc,formatPP,stdIn  &
!   imported scalars with intent(in)
    ,echo

  USE netcdf, ONLY :  &
!  imported procedures
     nf90_close,nf90_create,nf90_inquire,nf90_open  &
!  imported scalar parameters
    ,nf90_clobber,nf90_noClobber,nf90_noErr,nf90_noWrite

  USE rwErr_mod, ONLY : &
!  imported procedures
     rwErr

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER ::   &!  scalar parameters
     & nUnitMax=50    !  maximum unit used to connect to a file

  INTEGER, SAVE ::   &!  arrays
     & irecPrev(nUnitMax)  !  the last record number used for input from
!                               a given unit. Only used for sequential
!                               access files (ASCII and PP).

! Note: irecPrev(i) gives a record number for the file open on unit i.
! netCDF files don't use this "unit", rather they are given a netCDF
! ID number, and they also don't use irecPrev. However, calls to
! readVar include an argument of the likes of irecPrev(unit). We don't
! need this for netCDF files, and we don't want this to go out of
! bounds if the netCDF ID is large - this is currently achieved by
! passing irecPrev(1) for all netCDF files....which we get away with
! as long as irecPrev is not then updated when reading from netCDF!
!
! Potentially there is a similar problem with irecOut(:) - but that's
! not a problem at present as there is no output to netCDF files!

  CONTAINS

!#######################################################################
!#######################################################################
!#######################################################################
! function fileUnit
! Internal procedure in module fileUtils.
! Assign a unit for connecting next file.
! Loops over a prescibed range of units and selects first unit tha
! is not already in use.
! Must be called immediately prior to a file being opened, otherwise
! intervening code may open another file on the same unit.
! Also requires that files are closed once finished with, otherwise
! unit remains unavailable.
! If fileFormat=netCDF, does nothing (since netCDF interface does not
! need to be told a file unit).

  FUNCTION fileUnit( fileFormat ) RESULT(unit)

  USE inout, ONLY :  &
!  imported scalar parameters
    formatNc,stdIn,stdOut

  IMPLICIT NONE
!-------------------------------------------------------------------------------

  INTEGER :: unit   !  function result

  INTEGER :: i          !  local scalars

  LOGICAL :: fileOpen   !  local scalars

  CHARACTER(len=*), INTENT(in) ::  &!  IN scalars
    fileFormat           !  format of file to be opened
!-----------------------------------------------------------------------

  unit = 0

  IF ( fileFormat /= formatNc ) THEN

    DO i=1,nUnitMax

!     Avoid units for standard i/o, hard-wired!
      IF ( i==stdIn .OR. i==stdOut ) CYCLE

!     Find out if anything is connected to this unit.
      INQUIRE (unit=i, opened=fileOpen )
      IF ( .NOT. fileopen ) THEN
!       This unit is free, so exit.
        unit = i
        EXIT
      ENDIF

    ENDDO

    IF ( unit < 1 ) THEN
      WRITE(*,*)'All allowed units are in use. '
      WRITE(*,*)'Searched up to nUnitMax=',nUnitMax
      WRITE(*,*)'No unit free to open another file.'
      WRITE(*,*)'Increase nUnitMax.'
      WRITE(*,*)'Stopping in fileUnit'
      STOP
    ENDIF

  ENDIF  !  fileFormat

  END FUNCTION fileUnit

!###############################################################################
!###############################################################################
!###############################################################################
! subroutine closeFile
! Internal procedure in module fileUtils.
! Close a file connected to a given unit.

  SUBROUTINE closeFile( unit,fileFormat )

  IMPLICIT NONE

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) :: unit   !  the unit to be closed
!                                  For netCDF files, unit is the netCDF ID.
  CHARACTER(len=*), INTENT(in) :: fileFormat       !  format of file

! Local scalar variables.
  INTEGER :: ierr   !  error flag


!-----------------------------------------------------------------------
  SELECT CASE ( fileFormat )

    CASE ( formatNc )
!     Check if there is a file with this netCDF ID. If the inquiry
!     returns an error, assume that there is no such file open and
!     so nothing to close.
      ierr = nf90_inquire( unit )
      IF ( ierr==nf90_noErr ) THEN
        IF ( echo ) WRITE(*,*)'Closing file with netCDF ID=',unit
        ierr = nf90_close( unit )
!       If error, call error routine, but don't stop.
        IF ( ierr /= nf90_noErr ) CALL rwErr( formatNc,.FALSE.  &
              ,ierr=ierr,ncID=unit,errMess1='closeFile nf90_close' )
      ENDIF

    CASE default
!     Used for all file formats other than netCDF.
!     Note that no file is closed if the unit=standard in.
      IF ( unit /= stdIn ) THEN
        IF ( echo ) WRITE(*,*)'Closing file open on unit=',unit
        CLOSE( unit, iostat=ierr )

!       Ignore any error. This allows this procedure to be called
!       without worrying whether a file is actually connected to the unit.
!        if ( ierr/=0 ) then
!          write(*,*)'WARNING: closeFile: unit=',unit,' iostat=',ierr
!          inquire( unit, exist=existArg, opened=openedArg )
!          if ( .NOT. existArg ) write(*,*)'This unit does not exist.'
!          if ( .NOT. openedArg ) write(*,*)'No file is open on this unit.'
!        endif

      ENDIF

  END SELECT

  END SUBROUTINE closeFile

!#######################################################################
!#######################################################################
!#######################################################################

! subroutine commentLine
! Internal procedure in module fileUtils.

  SUBROUTINE commentLine(unit,ierr)

! Reads the first character of a line of text (from formatted input).
! If the character indicates a comment line (! or #), reads more lines
! until finds one that is not a comment. Leaves file in position to read
! first non-comment line. Has to read from an external file so that
! can use backspace.

  IMPLICIT NONE
!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: unit    !   unit via which to read file

!-------------------------------------------------------------------------------
! Scalar arguments with intent(out)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(out) :: ierr   !   return status (error code)

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  CHARACTER(len=1) :: char1  !  the first character of each non-blank line

! Loop over lines in input.
  DO

!   Read fist character on line.
    READ(unit,*,iostat=ierr) char1

    IF ( ierr /= 0 ) EXIT

    IF ( char1=='!' .OR. char1=='#' ) THEN
!     This is a comment line, so cycle to top.
      CYCLE
    ELSE
!     This is not a comment line. Reposition file at start of line.
      BACKSPACE unit
      EXIT

    ENDIF

  ENDDO

  END SUBROUTINE commentLine

!#######################################################################
!#######################################################################
!#######################################################################

! subroutine findTag
! Internal procedure in module fileUtils.
! Reads the first characters of a line of text (from formatted input).
! If the string matches the input string, exits. If not, reads more
! lines until a match is found. Stops at any error (eg end of file).
! Includes option to restrict search to current "section" of file (this
! is hardwired to expect each section to begin ">INIT").

  SUBROUTINE findTag( unit,message,tag,preInit )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::   &!  in scalars
    unit    !   unit via which to read file

  CHARACTER(len=*), INTENT(in) ::  &!  IN scalars
    message  &!  message written at error
   ,tag       !  the string to search for

  LOGICAL, INTENT(in), OPTIONAL ::  &!  optional in scalars
    preInit   !  TRUE indicates that the tag must be found before
!                  any tag starting '>INIT'. This is used to
!                  ensure that information is read from the same
!                  section of the file (i.e. it's not that we've skipped
!                  to a much later section of the file).

  INTEGER ::  &!  scalars
    ierr      &!  error code
   ,clen       !  work

  LOGICAL ::  &!  local scalars
    initTag     &!  TRUE if a line begins with ">INIT".
   ,preInitVal   !  local version of preInit

  CHARACTER(len=MAX(LEN(tag),5)) ::  &!  local scalars
    instring   !  length must be enough for tag and for ">INIT"

  CHARACTER(len=5) :: fmt  !  format used for read
!-------------------------------------------------------------------------------

! Deal with optional argument.
  preInitVal = .FALSE.
  IF ( PRESENT(preInit) ) preInitVal = preInit

  clen = LEN( instring )
  IF ( clen < 10 ) THEN
    WRITE(fmt,"('(a',i1,')')") clen
  ELSEIF ( clen < 100 ) THEN
    WRITE(fmt,"('(a',i2,')')") clen
  ELSE
    WRITE(*,*)'findTag: search string=',TRIM(instring),' too long.'
    WRITE(*,*)'search string length must be <100 characters.'
    WRITE(*,*)'Stopping in findTag'
    STOP
  ENDIF

  DO     !   do until succesful search or read error

    READ(unit,fmt=fmt,iostat=ierr) instring
    IF ( ierr /= 0 ) THEN
      WRITE(*,*) 'ERROR: findTag: ',TRIM(message)
      WRITE(*,*)'findTag: Error looking for ',TRIM(tag)
      WRITE(*,*)'i.e. was looking for this text in file'
      IF ( ierr < 0 ) THEN
        WRITE(*,*)'findTag: Error: end of file (iostat=',ierr,')'
      ELSE
        WRITE(*,*)'findTag: Error reading from file. iostat=',ierr
      ENDIF
      STOP
    ENDIF

!   Find out if this line starts ">INIT".
    initTag = .FALSE.
    IF ( instring(1:5) == '>INIT' ) initTag = .TRUE.

    IF ( instring /= tag ) THEN
      IF ( preInitVal .AND. initTag ) THEN
        WRITE(*,*) 'findTag: ERROR: ',TRIM(message)
        WRITE(*,*) 'Error looking for tag=',TRIM(tag)
        WRITE(*,*) 'Did not find tag in this section.'
        STOP
      ELSE
        CYCLE !  cycle to top to read next line
      ENDIF
    ELSE
      EXIT  !  successful search
    ENDIF

  ENDDO

  END SUBROUTINE findTag

!#######################################################################
!#######################################################################
!#######################################################################
! subroutine openFile
! Internal procedure in module fileUtils.
! Open a file with specified attributes.

  SUBROUTINE openFile( nval,closeFirst,unit,action,fileFormat,filename,status  &
                      ,ncCallType,ncFileType )

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: nval !  number of data values per record. Only used
!                                for unformatted files.

  LOGICAL, INTENT(in) :: closeFirst  !  TRUE means first close any file already
!                            connected on unit before opening new file.

  CHARACTER(len=*), INTENT(in) ::  &!  SCALARS
    action          &!
   ,fileFormat      &!  format of file
   ,filename,status

!-------------------------------------------------------------------------------
! Scalar arguments with intent(inout)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(inout) :: unit   !  The unit to be used.
!               For netCDF files, unit is the netCDF ID, and has intent(inout).
!               For other values of fileFormat, this is the unit that
!                 is used to connect to the file, and has intent(in).

!-------------------------------------------------------------------------------
! Optional scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  CHARACTER(len=*), INTENT(in), OPTIONAL ::  &!  optional in scalars
    ncCallType    &!  flag indicating what type of data are to be read - a hack
!                     for netCDF
   ,ncFileType     !  indicates what "type" of netCDF files are to be read

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    fStatus   &!  creation mode argument for netCDf files
   ,ierr      &!  error code
   ,reclOne    !  record length for a single word.

  REAL ::  &
    tmpReal   !  work

!-----------------------------------------------------------------------
! If requested, close any existing file.
! closeFile does not raise an error if no file is open on this unit.
!-----------------------------------------------------------------------
  IF ( closeFirst ) CALL closeFile( unit,fileFormat )

!-----------------------------------------------------------------------
! Optional print to screen.
!-----------------------------------------------------------------------
  IF ( echo ) THEN
    IF ( fileFormat /= 'scratch' ) THEN
      WRITE(*,*) 'Opening ',TRIM(filename)
      WRITE(*,*) 'fileFormat=',fileFormat
      IF ( fileFormat /= formatNc ) WRITE(*,*) ' action=',action,' status=',status
    ELSE
      WRITE(*,*) 'Opening scratch file.'
    ENDIF
    IF ( fileFormat /= formatNc ) WRITE(*,*)'unit=',unit
  ENDIF

!-----------------------------------------------------------------------
! Open file.
!-----------------------------------------------------------------------

  SELECT CASE ( fileFormat )

!--------------------------------------------------
    CASE ( formatAsc )
      OPEN(unit,file=filename,action=action,status=status)
!     Reset irecPrev.
      irecPrev(unit) = 0

!--------------------------------------------------
    CASE ( formatBin )
!     Establish length of record - this is processor dependent.
!     FOR NOW, assuming we want to write a 'default' real value!
      INQUIRE( iolength=reclOne ) tmpReal
      IF ( echo ) WRITE(*,*)'recl=',reclOne*nval
      OPEN( unit,file=filename,form='unformatted',access='direct'  &
                ,recl=reclOne*nval,action=action,status=status )
!     Reset irecPrev.
      irecPrev(unit) = 0

!--------------------------------------------------
    CASE ( formatNc )

      IF ( status == 'old' ) THEN

!      Open file for read only - currently OK for all 'old' files.
       fstatus = nf90_noWrite
       ierr = nf90_open( filename,fStatus,unit )

!      If error, call error routine, indicating a stop.
       IF ( ierr /= 0 ) CALL rwErr( formatNc,.TRUE.  &
           ,ierr=ierr,fileName=fileName,ncID=unit,errMess1='openFile nf90_open' )

!       If the "type" of netCDF call has been provided, check that
!       that the dimensions in the file are as expected.
        IF ( PRESENT(ncCallType) )   &
             CALL check_nc_dims( unit,ncCallType,ncFileType )

      ELSE

!       Create a new netCDF file.
!       Get value for file status.
        IF ( status == 'new' ) THEN
          fStatus = NF90_NOCLOBBER
        ELSE
          fStatus = NF90_CLOBBER
        ENDIF
 
!      Open (create) the file.
       ierr = nf90_create( filename,fStatus,unit )

!      If error, call error routine, indicating a stop.
       IF ( ierr /= 0 ) CALL rwErr( formatNc,.TRUE.  &
              ,ierr=ierr,fileName=fileName,ncID=unit,errMess1='createNcFile nf90_create' )

      ENDIF  !  status

      IF ( echo ) WRITE(*,*)'netCDF ID=',unit

!--------------------------------------------------
    CASE ( formatPP )
      OPEN(unit,file=filename,form='unformatted',access='sequential',action=action,status=status)
!     Reset irecPrev.
      irecPrev(unit) = 0

!--------------------------------------------------
    CASE ( 'scratch' )
      OPEN( unit,status='scratch' )

!--------------------------------------------------
    CASE default
      WRITE(*,*)'ERORR: openFile: do not recognise fileFormat=',TRIM(fileFormat)
      STOP

  END SELECT

  END SUBROUTINE openFile

!#######################################################################
!#######################################################################
!#######################################################################

  END MODULE file_utils
