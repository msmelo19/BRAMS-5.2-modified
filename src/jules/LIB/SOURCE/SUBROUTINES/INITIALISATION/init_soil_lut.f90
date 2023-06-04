!###############################################################################
  SUBROUTINE init_soil_lut( land_pts,nsoilVar,sm_levels,albSoilConst,soilLUTfile  &
                          ,soilType )

!-------------------------------------------------------------------------------
! Read a look-up table (LUT) of soil characteristics for a number of soil types,
! and then set spatial soil variables.
!-------------------------------------------------------------------------------

  USE file_utils, ONLY : &
!  imported procedures
     closeFile,commentLine,fileUnit,openFile

  USE inout, ONLY : &
!  imported scalar parameters
     formatAsc,formatLen

  USE p_s_parms, ONLY : &
!  imported arrays with intent(out)
     albsoil,b,hcap,hcon,satcon,sathh,smvccl,smvcst,smvcwt

  USE soil_param, ONLY :  &
!  imported arrays with intent(out)
     dzsoil

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: land_pts   !  number of land points
  INTEGER, INTENT(in) :: nsoilVar   !  number of soil variables required for chosen configuration
!         At present this is 2 larger than number read here, because soilType and
!         soilAlbConst are included.
  INTEGER, INTENT(in) :: sm_levels  !  number of soil layers
  REAL, INTENT(in) :: albSoilConst  !  prescribed soil albedo
  CHARACTER(len=*), INTENT(in) :: soilLUTfile   !  name of the
!               file containing look-up table with characteristics of  each soil type

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: soilType(land_pts) ! soil type at each gridpoint

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER :: i        !  loop counter  
  INTEGER :: ierr     !  error value
  INTEGER :: ierrSum  !  error counter
  INTEGER :: inUnit   !  unit used to connect to input file
  INTEGER :: isoil    !  work
  INTEGER :: iz       !  loop counter  
  INTEGER :: j        !  loop counter  
  INTEGER :: nSoil    !  number of soil types
  INTEGER :: nzSoilIn !  number of soil layers used when preparing input data
  LOGICAL :: found    !  work
  CHARACTER(len=formatLen) :: fileFormat   !  format of file

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  REAL :: dzSoilIn(sm_levels)          ! soil layer thicknesses (m)
  INTEGER, ALLOCATABLE :: soilNum(:)  !  soil type numbers
  REAL, ALLOCATABLE :: soilChar(:,:,:) !  soil characteristics for each soil and
!                                              each layer

!-------------------------------------------------------------------------------
! This code is not generic, and it's much easier if we prescribe the order in
! which the soil characteritics appear in the file.
! Stop if number of variables isn't as expected (a minimum test).
!-------------------------------------------------------------------------------
  IF ( nsoilVar /= 10 ) THEN
    WRITE(*,*)'ERROR: init_soil_lut: nsoilVar not as expected.'
    WRITE(*,*)'nsoilVar=',nsoilVar
    WRITE(*,*)'Stopping in init_soil_lut'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Open ASCII LUT file on external unit.
!-------------------------------------------------------------------------------
  fileFormat = formatAsc
  inUnit = fileUnit( fileFormat )  !  get unit
  CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,soilLUTfile,'old' )

!-------------------------------------------------------------------------------
! Skip header lines (starting with # or !).
!-------------------------------------------------------------------------------
  CALL commentLine( inUnit,ierr )
  IF ( ierr /= 0 ) THEN
    WRITE(*,*)'ERROR: setSoilChar: error while skipipng header lines.'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Read details of soil layers that were assumed when creating input file.
!-------------------------------------------------------------------------------
  READ(inUnit,*) nzSoilIn
  IF ( nzSoilIn /= sm_levels ) THEN
    WRITE(*,*)'ERROR: readSoilChar: nzSoilIn /= sm_levels.'
    WRITE(*,*)'The number of model soil layers does not match that used'
    WRITE(*,*)'when creating input data.'
    WRITE(*,*)'soilLUTfile=',TRIM(soilLUTfile)
    WRITE(*,*)'Stopping in init_soil_lut'
    STOP
  ENDIF
  READ(inUnit,*) dzSoilIn(:)
  DO iz=1,nzSoilIn
    IF ( ABS( dzSoil(iz)-dzSoilIn(iz) ) > EPSILON(dzSoil) ) THEN
      WRITE(*,*)'ERROR: readSoilChar: dzSoil /= dzSoilIn.'
      WRITE(*,*)'Soil layer thickness does not match that used'
      WRITE(*,*)'when creating input data.'
      WRITE(*,*)'Model layers=',dzSoil(:)
      WRITE(*,*)'Input file layers=',dzSoilIn(:)
      WRITE(*,*)'soilLUTfile=',TRIM(soilLUTfile)
      WRITE(*,*)'Stopping in readSoilChar'
    STOP
    ENDIF
  ENDDO

!-------------------------------------------------------------------------------
! Read number of soil types and allocate memory.
!-------------------------------------------------------------------------------
  READ(inUnit,*) nsoil
  WRITE(*,*)'Number of defined soil types=',nsoil
  ierrSum = 0
  ALLOCATE(soilNum(nsoil), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE(soilChar(nsoil,nsoilVar-2,sm_levels), stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) THEN
    WRITE(*,*)'ERROR: readSoilChar: error allocating memory.'
    WRITE(*,*)'nsoil=',nsoil,' ierrSum=',ierrSum
    WRITE(*,*)'Stopping in readSoilChar'
  ENDIF

!-------------------------------------------------------------------------------
! Read characteristics for each soil type.
!-------------------------------------------------------------------------------
  DO i=1,nsoil
    READ(inUnit,*) soilNum(i)
    DO iz=1,sm_levels
      READ(inUnit,*) soilChar(i,:,iz)
    ENDDO
  ENDDO

!-------------------------------------------------------------------------------
! Close file.
!-------------------------------------------------------------------------------
  CALL closeFile( inUnit,fileFormat )

!-------------------------------------------------------------------------------
! Now set the characteristics at each point.
! Note that each gridbox has a single soil type, but each type has
! characteristics for each layer.
!-------------------------------------------------------------------------------

  DO i=1,land_pts

!   Find soil type in list.
    found = .FALSE.
    DO j=1,nsoil
      IF ( soilType(i) == soilNum(j) ) THEN
        isoil = j
        found = .TRUE.
      ENDIF
    ENDDO
    IF ( .NOT. found ) THEN
      WRITE(*,*)'ERORR: init_soil_lut: soil type not found.'
      WRITE(*,*)'Land point #',i,' soil number=',soilType(i)
      STOP
    ENDIF

!-------------------------------------------------------------------------------
!   Set soil characteristics at this point.
!   Note that this assumes that the data were provided in this order.
!-------------------------------------------------------------------------------
    sathh(i,:) = soilChar(isoil,1,:)
    b(i,:) = soilChar(isoil,2,:)
    hcap(i,:) = soilChar(isoil,3,:)
    hcon(i,1:sm_levels) = soilChar(isoil,4,:)
    satcon(i,1:sm_levels) = soilChar(isoil,5,:)
    smvccl(i,:) = soilChar(isoil,6,:)
    smvcst(i,:) = soilChar(isoil,7,:)
    smvcwt(i,:) = soilChar(isoil,8,:)

  ENDDO

!-------------------------------------------------------------------------------
! Set soil albedo.
! In this implementation, this is a constant.
!-------------------------------------------------------------------------------
  albsoil(:) = albSoilConst

!-------------------------------------------------------------------------------
! Deallocate memory.
!-------------------------------------------------------------------------------
  ierrSum = 0
  DEALLOCATE(soilNum, stat=ierr ); ierrSum=ierrSum+ierr
  DEALLOCATE(soilChar, stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) THEN
    WRITE(*,*)'WARNING: readSoilChar: error deallocating memory.'
    WRITE(*,*)'ierrSum=',ierrSum
  ENDIF

  END SUBROUTINE init_soil_lut
!###############################################################################
!###############################################################################

