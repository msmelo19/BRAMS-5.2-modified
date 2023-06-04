! Read TRIFFID parameters.

  SUBROUTINE init_trif()

  USE file_utils, ONLY :  &
!   imported procedures
     closeFile,fileUnit,findTag,openFile

  USE inout, ONLY :  &
!  imported scalar parameters
    formatAsc,formatLen,jinUnit  &
!  imported scalars with intent(in)
   ,echo

  USE misc_utils, ONLY :  &
!  imported procedures
     repeatVal

  USE nstypes, ONLY :   &
!  imported scalars with intent(in)
    npft

  USE pftparm, ONLY :  &
!  imported arrays with intent(in)
     pftName

  USE switches, ONLY :  &
!  imported scalars with intent(in)
    l_phenol,l_triffid,routeOnly

  USE trif, ONLY :  &
!   imported arrays with intent(out)
     crop,g_area,g_grow,g_root,g_wood,lai_max,lai_min

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  SCALARS
    nvar=7  !  # of parameter fields to read (excluding trifName)

  INTEGER ::  &!  local SCALARS
    i            &!  loop counter
   ,ierr         &!  error flag
   ,ierrSum      &!  error flag
   ,inUnit       &!  unit used to connect to input file
   ,j            &!  loop counter
   ,ntypeInFile   !  number of surface types (all PFTs) in the input file

  REAL, ALLOCATABLE ::  &!  local ARRAYS
    tmpval(:,:)          !  work

  LOGICAL ::  &!  local SCALARS
    readFile  &!  flag indicating if another file is to be read
   ,rowWise    !  T if data for each type are arranged across rows
!              !  F if data for each type are arranged down a column

  LOGICAL ::  &!  local arrays
    foundPFT(npft)   !  work

  CHARACTER(len=formatLen) ::  &! local SCALARS
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &! local SCALARS
    fileName   !  the name of a file

  CHARACTER(len=LEN(pftName)), ALLOCATABLE ::  &!  local arrays
    tmpName(:)   !  work - used to hold PFT names
!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_trif'

! Nothing to do if neither triffid not phenology are selected.
  IF ( .NOT.l_triffid .AND. .NOT.l_phenol ) RETURN

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_trif','>INIT_TRIF' )

! Establish where data are to be read from, and how many surface types there are in the file.
  READ(jinUnit,*) readFile
  READ(jinUnit,*) fileName
  READ(jinUnit,*) ntypeInFile

!-----------------------------------------------------------------------

! Open file.
  fileFormat = formatAsc
  IF ( readFile ) THEN
!   Get unit
    inUnit = fileUnit( fileFormat )
!   First arg is not needed for formatted file, so using 1.
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old' )
  ELSE
    WRITE(*,*)'Reading parameters for veg types from the run control file.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_trif','>DATA',preInit=.TRUE. )
  ENDIF

! Allocate space for all types, including those not ultimately kept.
  ierrSum = 0
  ALLOCATE( tmpval(nvar,ntypeInFile), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( tmpName(ntypeInFile), stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) THEN
    WRITE(*,*)'ERROR: init_veg: could not allocate memory.'
    STOP
  ENDIF


! I couldn't decide whether data were best arranged in rows, or columns, so I've coded both but
! am only using one for now.
! For now, assume rowwise or columnwise data.
  rowWise = .FALSE.

! Read the data.
  IF ( rowWise ) THEN
!   Data are arranged as all variables for type 1, all variables for type 2, etc...
    DO i=1,ntypeInFile
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
    WRITE(*,*)'ERROR: init_trif: repeated name of PFT.'
    WRITE(*,*)'All PFTs must have unique names.'
    IF ( readFile ) THEN
      WRITE(*,*)'Error in file ',TRIM(fileName)
    ELSE
      WRITE(*,*)'Error reading from standard input.'
    ENDIF
    STOP
  ENDIF

! Extract the data for the required PFTs.
  foundPFT(:) = .FALSE.
  DO j=1,nTypeInFile
    DO i=1,npft
      IF ( .NOT.foundPFT(i) .AND. pftName(i)==tmpName(j) ) THEN
        foundPFT(i) = .TRUE.
        crop(i) = NINT( tmpval(1,j) )
        g_area(i) = tmpval(2,j)
        g_grow(i) = tmpval(3,j)
        g_root(i) = tmpval(4,j)
        g_wood(i) = tmpval(5,j)
        lai_max(i) = tmpval(6,j)
        lai_min(i) = tmpval(7,j)
      ENDIF
    ENDDO  !  i  (pft)
  ENDDO   !  j (type in file )

  IF ( ANY( .NOT.foundPFT(:) ) ) THEN
    WRITE(*,*)'ERROR: init_trif: could not find one or more PFTs.'
    WRITE(*,*)'Could not find ',COUNT( .NOT.foundPFT(:) ),' of ',npft,' PFTs.'
    WRITE(*,*)'Could not find PFTs called:'
    DO i=1,npft
      IF ( .NOT. foundPFT(i) ) WRITE(*,*) TRIM(pftName(i))
    ENDDO
    WRITE(*,*)'NB Names are sensitive to case!'
    IF ( .NOT. readFile ) THEN
      WRITE(*,*)'Reading parameters from the run control file.'
    ELSE
      WRITE(*,*)'Reading parameters from ',TRIM(fileName)
    ENDIF
    STOP
  ENDIF

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

  DEALLOCATE( tmpval, tmpName, stat=ierr )
  IF ( ierr /= 0 ) WRITE(*,*) 'WARNING: init_veg: could not deallocate memory.'
!------------------------------------------------------------------------------

! Write to screen.
  IF ( echo ) THEN
!    write(*,"(a,10(tr1,a))") 'pftName=',(/ (trim(pftName(i)), i=1,npft ) /)  !  not working!
    WRITE(*,"(a,10(tr1,a))") 'pftName=',pftName(:)
    WRITE(*,*) 'crop=',crop(:)
    WRITE(*,*) 'g_area=',g_area(:)
    WRITE(*,*) 'g_grow=',g_grow(:)
    WRITE(*,*) 'g_root=',g_root(:)
    WRITE(*,*) 'g_wood=',g_wood(:)
    WRITE(*,*) 'lai_max=',lai_max(:)
    WRITE(*,*) 'lai_min=',lai_min(:)
  ENDIF

  END SUBROUTINE init_trif
