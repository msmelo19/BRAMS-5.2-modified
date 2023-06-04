! Read parameters for non-veg surface types.

  SUBROUTINE init_nonveg()

  USE c_z0h_z0m, ONLY : z0h_z0m
  USE file_utils, ONLY: closeFile,fileUnit,findTag,openFile
  USE inOut, ONLY : echo,formatAsc,formatLen,jinUnit

  USE misc_utils, ONLY :  &
!  imported procedures
     repeatVal

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     nnvg,npft,        &
!  imported scalars with intent(out)
     ice,lake,soil,urban,urban_canyon,urban_roof

  USE nvegparm

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     routeOnly

  USE switches_urban, ONLY :  &
!  imported scalars with intent(out)
     l_urban2t,l_moruses

! The following is required for the URBAN-2T & MORUSES code and will be removed
! in a future release
  USE urban_param, ONLY :  &
     !  imported scalars with intent(out)
     albsnc_c, albsnc_rf, albsnf_c, albsnf_rf, catch_c, catch_rf,        &
     gs_c, gs_rf, infil_c, infil_rf, z0_c, z0_rf, z0h_z0m_c, z0h_z0m_rf, &
     ch_c, ch_rf, vf_c, vf_rf, emis_c, emis_rf

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  SCALARS
    nvar=10       !  # of parameter fields to read (not including name of type)

  INTEGER ::  &!  local SCALARS
    i            &!  loop counter
   ,ierr         &!  error flag
   ,ierrSum      &!  error flag
   ,inUnit       &!  unit used to connect to input file
   ,j            &!  work/loop counter
   ,nnvgInFile    !  number of non-veg surface types in the input file
!                      In fact this is the total number of surface types in the file.

  REAL, ALLOCATABLE ::  &!  local ARRAYS
    tmpval(:,:)          !  work

  LOGICAL ::  &!  local SCALARS
    readFile  &!  flag indicating if another file is to be read
   ,rowWise    !  T if data for each type are arranged across rows
!              !  F if data for each type are arranged down a column

  LOGICAL ::  &!  local arrays
    foundNvg(nnvg)   !  work

  CHARACTER(len=LEN(nvgName)), PARAMETER ::  &!  scalar parameters
! If using Intel compiler (9.0 tested), you probably have to remove the previous line
! and replace with the next line. Make sure length is sufficient - see length of
! nvgName in module nvegparm.
!  character(len=20), parameter ::  &!  scalar parameters
    iceName = 'ice'     &!  value of nvgName that indicates the ice type
   ,lakeName = 'lake'   &!  value of nvgName that indicates the lake type
   ,soilName = 'soil'   &!  value of nvgName that indicates the soil type
   ,urbanName = 'urban' &!  value of nvgName that indicates the urban type
   ,canyonName = 'urban_canyon' &!  value of nvgName that indicates canyon type
   ,roofName   = 'urban_roof'    !  value of nvgName that indicates roof type

  CHARACTER(len=formatLen) ::  &! local SCALARS
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &! local SCALARS
    fileName   !  the name of a file

  CHARACTER(len=LEN(nvgName)), ALLOCATABLE ::  &!  local arrays
    tmpName(:)   !  work - used to hold names of types
!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_nonveg'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_nonveg','>INIT_NONVEG' )

! Establish where data are to be read from.
  READ(jinUnit,*) readFile
  READ(jinUnit,*) fileName

! Establish how many surface types are in the file.
  READ(jinUnit,*) nnvgInFile

! Ensure that these values are reasonable.
  IF ( nnvgInFile < nnvg ) THEN
    WRITE(*,*)'Need nnvg=',nnvg,' have only ',nnvgInFile
    WRITE(*,*)'Stopping in init_nonveg'
    STOP
  END IF

!-----------------------------------------------------------------------

! Open file.
  fileFormat = formatAsc
  IF ( readFile ) THEN
!   Get unit
    inUnit = fileUnit( fileFormat )
!   First arg is not needed for formatted file, so using 1.
    CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old' )
  ELSE
    IF ( echo ) WRITE(*,*)'Reading parameters for non-veg types from the run control file.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_nonveg','>DATA', preInit=.TRUE. )
  END IF

! I couldn't decide whether data were best arranged in rows, or columns, so I've coded both but
! am only using one for now.
! Allocate space for all types, including those not ultimately kept.
  ierrSum = 0
  ALLOCATE( tmpval(nvar,nnvgInFile), stat=ierr ); ierrSum=ierrSum+ierr
  ALLOCATE( tmpName(nnvgInFile), stat=ierr ); ierrSum=ierrSum+ierr
  IF ( ierrSum /= 0 ) THEN
    WRITE(*,*)'ERROR: init_nonveg: could not allocate memory.'
    STOP
  END IF

! For now, assume rowwise or columnwise data.
  rowWise = .FALSE.

! Read the data.
  IF ( rowWise ) THEN
!   Data are arranged as all variables for type 1, all variables for type 2, etc...
    DO i=1,nnvgInFile
      READ(inUnit,*) tmpName(i),tmpval(:,i)
    END DO
  ELSE
!   Data are arranged as variable 1 for all types, variable 2 for all types, etc...
    READ(inUnit,*) tmpName(:)
    DO i=1,nvar
      READ(inUnit,*) tmpval(i,:)
    END DO
  END IF

! Check that names of surface types are unique.
  IF ( repeatVal( tmpName(:) ) ) THEN
    WRITE(*,*)'ERROR: init_nonveg: repeated name of surface type.'
    WRITE(*,*)'All surface types must have unique names.'
    IF ( readFile ) THEN
      WRITE(*,*)'Error in file ',TRIM(fileName)
    ELSE
      WRITE(*,*)'Error reading from standard input.'
    END IF
    STOP
  END IF

! Extract the data for the required types.
  foundNvg(:) = .FALSE.
  DO j=1,nnvgInFile
    DO i=1,nnvg
      IF ( .NOT.foundNvg(i) .AND. nvgName(i)==tmpName(j) ) THEN
        foundNvg(i) = .TRUE.
        albsnc_nvg(i) = tmpval(1,j)
        albsnf_nvg(i) = tmpval(2,j)
        catch_nvg(i) = tmpval(3,j)
        gs_nvg(i) = tmpval(4,j)
        infil_nvg(i) = tmpval(5,j)
        z0_nvg(i) = tmpval(6,j)
        z0h_z0m(npft+i) = tmpval(7,j)   !  assuming nvg follow PFTs
        ch_nvg(i) = tmpval(8,j)
        vf_nvg(i) = tmpval(9,j)
        emis_nvg(i) = tmpval(10,j)
      END IF
    END DO  !  i  (nvg)
  END DO  !  j (nvg in file)

  IF ( ANY( .NOT.foundNvg(:) ) ) THEN
    WRITE(*,*)'ERROR: init_veg_nvg: could not find one or more non-veg types.'
    WRITE(*,*)'Could not find ',COUNT( .NOT.foundNvg(:) ),' of ',nnvg,' types.'
    WRITE(*,*)'Could not find types called:'
    DO i=1,nnvg
      IF ( .NOT. foundNvg(i) ) WRITE(*,*) TRIM(nvgName(i))
    END DO
    WRITE(*,*)'NB Names are sensitive to case!'
    IF ( .NOT. readFile ) THEN
      WRITE(*,*)'Reading parameters for types from the run control file.'
    ELSE
      WRITE(*,*)'Reading parameters from ',TRIM(fileName)
    END IF
    STOP
  END IF

! Close file if it is not the run control file
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

  DEALLOCATE( tmpval,tmpName, stat=ierr )
  IF ( ierr /= 0 ) WRITE(*,*) 'WARNING: init_nonveg: could not deallocate memory.'

!------------------------------------------------------------------------------


! Check that we have soil type, and also check if lake and ice types available.
! We need the soil type since this is used for the bare soil under PFTs.

! Note that the non-veg surface types must come after the PFTs since
! (a) PHYSIOL refers to type soil-npft
! (b) INIT_NONVEG refers to z0h_z0m(npft+i)

  ice          = -1
  lake         = -1
  soil         = -1
  urban        = -1
  urban_canyon = -1
  urban_roof   = -1

  DO i=1,nnvg
    SELECT CASE ( nvgName(i) )
      CASE ( iceName )
        ice = npft + i
      CASE ( lakeName )
        lake = npft + i
      CASE ( soilName )
        soil = npft + i
      CASE ( urbanName )
        urban = npft + i
      CASE ( canyonName )
        urban_canyon = npft + i
      CASE ( roofName )
        urban_roof = npft + i
      END SELECT
  END DO

  IF ( soil < 1 ) THEN
    WRITE(*,*)'ERROR: init_nonveg: no surface type "soil"'
    WRITE(*,*)'There must be a soil surface type, called ',TRIM(soilName)
    STOP
  END IF

!------------------------------------------------------------------------------
! Do checks on urban tiles
  IF ( urban > 0 .AND. urban_canyon > 0 ) THEN
    WRITE(*,*) 'ERROR: init_nonveg: cannot have both urban surface types'
    WRITE(*,*) TRIM(urbanName),' and ',TRIM(canyonName)
    STOP
  END IF

! The urban_roof has to be specified to use the two-tile urban schemes. Check
! that the roof has the correct counterparts
  IF ( urban_roof > 0 ) THEN
    IF ( urban < 0 .AND. urban_canyon < 0 ) THEN
      WRITE(*,*) 'ERROR: init_nonveg: cannot have a roof without a canyon'
      WRITE(*,*) '  URBAN-2T: As urban surface type ',TRIM(roofName),     &
         ' is present there must'
      WRITE(*,*) '  also be an urban surface type called ', TRIM(canyonName)
      WRITE(*,*) '  Alternatively, to use URBAN-1T scheme remove ',       &
         'urban surface type ',TRIM(roofName)
      STOP
    END IF
  ELSE ! Check that if not present then neither is the canyon or MORUSES
    IF ( urban_canyon > 0 ) THEN
      WRITE(*,*) 'ERROR: init_nonveg: cannot have a canyon without a roof'
      WRITE(*,*) '  URBAN-2T: As urban surface type ',TRIM(canyonName),     &
         ' is present there must'
      WRITE(*,*) '  also be an urban surface type called ', TRIM(roofName)
      WRITE(*,*) '  Alternatively, to use URBAN-1T scheme rename ',         &
         TRIM(canyonName),' --> ',TRIM(urbanName)
      STOP
    END IF
    IF ( l_moruses ) THEN
      WRITE(*,*) 'ERROR: init_nonveg: MORUSES has no urban roof surface type'
      WRITE(*,*) 'There must be an urban surface type, called ',TRIM(roofName)
      STOP
    END IF
  END IF

! Either "urban_canyon" or "urban" can be specified in run control file as
! urban fraction as a whole is specified and then is split between canyon &
! roof depending on canyon fraction in init_urban.
  IF ( urban_roof > 0 ) THEN
    l_urban2t = .TRUE.
    IF ( urban_canyon < 0 ) urban_canyon = urban
    IF ( l_moruses ) THEN
      WRITE(*,*) 'WARNING: init_nonveg: Using urban surface scheme - MORUSES'
    ELSE
      WRITE(*,*) 'WARNING: init_nonveg: Using urban surface scheme - URBAN-2T'
    END IF
  END IF

!------------------------------------------------------------------------------
! Write to screen.
  IF ( echo ) THEN
    WRITE(*,*)'Soil is surface type #',soil
    IF ( lake > 0 ) THEN
      WRITE(*,*) 'Lake (inland water) is surface type #',lake
    ELSE
      WRITE(*,*) 'No lake (inland water) type has been indicated.'
    END IF
    IF ( ice > 0 ) THEN
      WRITE(*,*) 'Land ice is surface type #',ice
    ELSE
      WRITE(*,*) 'No land ice type has been indicated.'
    END IF
    IF ( urban > 0 ) THEN
      WRITE(*,*) 'Urban is surface type #',urban
    ELSE
      WRITE(*,*) 'No urban type has been indicated (URBAN-1T).'
    END IF
    IF ( urban_canyon > 0 ) THEN
      WRITE(*,*) 'Urban canyon is surface type #',urban_canyon
      IF ( urban > 0 ) WRITE(*,*)'  (Will be the same as urban surface type #)'
    ELSE
      WRITE(*,*) 'No canyon type has been indicated (URBAN-2T or MORUSES).'
    END IF
    IF ( urban_roof > 0 ) THEN
      WRITE(*,*) 'Urban roof is surface type #',urban_roof
    ELSE
      WRITE(*,*) 'No roof type has been indicated (URBAN-2T or MORUSES).'
    END IF
    WRITE(*,*)'Depending upon options, some of these parameters may not be used.'
    WRITE(*,"(a,10(tr1,a))") 'nvgName=',nvgName(:)
    WRITE(*,*) 'albsnc_nvg=',albsnc_nvg(:)
    WRITE(*,*) 'albsnf_nvg=',albsnf_nvg(:)
    WRITE(*,*) 'catch_nvg=',catch_nvg(:)
    WRITE(*,*) 'gs_nvg=', gs_nvg(:)
    WRITE(*,*) 'infil_nvg=', infil_nvg(:)
    WRITE(*,*) 'z0_nvg=', z0_nvg(:)
    WRITE(*,*) 'z0h_z0m=', z0h_z0m(npft+1:npft+nnvg)
    WRITE(*,*) 'ch_nvg=', ch_nvg(:)
    WRITE(*,*) 'vf_nvg=', vf_nvg(:)
    WRITE(*,*) 'emis_nvg=', emis_nvg(:)
  END IF

! The following is required for the URBAN-2T & MORUSES code and will be removed
! in a future release

  IF ( l_urban2t ) THEN
    albsnc_c   = albsnc_nvg(urban_canyon-npft)
    albsnc_rf  = albsnc_nvg(urban_roof-npft)

    albsnf_c   = albsnf_nvg(urban_canyon-npft)
    albsnf_rf  = albsnf_nvg(urban_roof-npft)

    catch_c    = catch_nvg(urban_canyon-npft)
    catch_rf   = catch_nvg(urban_roof-npft)

    gs_c       = gs_nvg(urban_canyon-npft)
    gs_rf      = gs_nvg(urban_roof-npft)

    infil_c    = infil_nvg(urban_canyon-npft)
    infil_rf   = infil_nvg(urban_roof-npft)

    z0_c       = z0_nvg(urban_canyon-npft)
    z0_rf      = z0_nvg(urban_roof-npft)

    z0h_z0m_c  = z0h_z0m(urban_canyon)
    z0h_z0m_rf = z0h_z0m(urban_roof)

    ch_c       = ch_nvg(urban_canyon-npft)
    ch_rf      = ch_nvg(urban_roof-npft)

    vf_c       = vf_nvg(urban_canyon-npft)
    vf_rf      = vf_nvg(urban_roof-npft)

    emis_c     = emis_nvg(urban_canyon-npft)
    emis_rf    = emis_nvg(urban_roof-npft)
  END IF

END SUBROUTINE init_nonveg

