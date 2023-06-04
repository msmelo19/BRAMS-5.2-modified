SUBROUTINE init_imogen
  
  USE alloc_mod, ONLY : allocate_arrays
  
  USE switches, ONLY : l_imogen,l_360
  
  USE file_utils, ONLY : closeFile,fileUnit,findTag,openFile
    
  USE inout, ONLY : formatAsc,jinUnit
  
  USE time_loc, ONLY : dateRun,timeRun,timeStep

  USE time_mod, ONLY : dateToBits,timeToBits
  
  USE spin_mod, ONLY : nspin
  
  USE timeConst, ONLY : iSecInDay
  
! Import namelists and associated variables
  USE imogen_anlg_vals
  USE imogen_run
  
  USE imogen_time
  
  USE imogen_progs, ONLY :                                        &
    CO2_PPMV,CO2_CHANGE_PPMV,DTEMP_O,FA_OCEAN,SEED_WG,SEED_RAIN
  
  USE imogen_constants, ONLY : n_imogen_land,drive_month
  
  USE imogen_map, ONLY : get_imogen_map,sgindinv
  
  USE imogen_clim
  
  USE imogen_io_vars, ONLY : YR_EMISS,C_EMISS
  
  IMPLICIT NONE
  
  INTEGER :: funit  ! File unit number for I/O
  INTEGER :: i,j,l,n ! Loop counters
  
  CHARACTER(len=180) ::                                           &
    namelistFile,                                                 &
                ! File containing the IMOGEN namelists
    imogenOrderFile,                                              &
                ! File containing the IMOGEN land mask/point order
    FILE_CLIM   ! WORK File containing initialising climatology.
    
  INTEGER ::                                                      &
    year,month,day,hours,mins,secs
       
    
! Nothing to do if IMOGEN is not enabled
  IF(.NOT. l_imogen) RETURN

!-----------------------------------------------------------------
! Read the values defined in the control file
!-----------------------------------------------------------------
  WRITE(*,"(50('-'),/,a)") 'init_imogen'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_imogen','>INIT_IMOGEN' )

  READ(jinUnit,*) imogenOrderFile
  READ(jinUnit,*) namelistFile

!-----------------------------------------------------------------
! Check that the setup of JULES time variables is compatible with
! IMOGEN and set up IMOGEN time variables
!-----------------------------------------------------------------
! -> Force 360 day year
  IF( .NOT. l_360 ) THEN
    WRITE(*,*)'ERROR: INIT_IMOGEN: 360 day year must be used (l_360 must be TRUE)'
    STOP
  ENDIF

! -> Check no spinup has been requested (a full IMOGEN experiment
!    consists of 3 JULES runs with the first 2 acting as a really
!    long spinup)
  IF( nspin /= 0 ) THEN
    WRITE(*,*)'ERROR: INIT_IMOGEN: Spinup has been requested'
    STOP
  ENDIF
  
! -> Run must start at 00:00 on 1st Jan for some year
  CALL dateToBits(dateRun(1),day,month,year,l_360,'DateToBits failed on start date')
  CALL timeToBits(timeRun(1),secs,mins,hours,'TimeToBits failed on start time')
! Since secs,mins,hours are all integers and >= 0, we can check for midnight
! by checking that the sum of them is 0
  IF((secs + mins + hours) /= 0 .OR. day /= 1 .OR. month /= 1) THEN
    WRITE(*,*)'ERROR: INIT_IMOGEN: Run must start at 00:00 on 1st Jan for some year'
    STOP
  ENDIF
! Store the start year in the IMOGEN variables
  YEAR1   = year
  IYSTART = year
  IYEAR   = year

! -> Set the IMOGEN end year variable - although it doesn't make
!    much sense to have an IMOGEN run ending in the middle of
!    a year, it won't mess anything up, so we allow it
  CALL dateToBits(dateRun(2),day,month,year,l_360,'DateToBits failed on start date')
  IYEND = year
  
! -> Set steps per day and check it is not more than the maximum
!    number of timesteps per day allowed by IMOGEN
!    Note that we already know that iSecInDay MOD timestep is 0
!    since this is enforced in init_time
  STEP_DAY = iSecInDay / NINT(timestep)
  IF(STEP_DAY > NSDMAX) THEN
    WRITE(*,*)'ERROR: INIT_IMOGEN: Too many timesteps per day'
    STOP
  ENDIF
  
!-----------------------------------------------------------------
! Check that the grid matches the IMOGEN grid (HadCM3 grid)
! This restriction may be lifted at some point in the future
!-----------------------------------------------------------------
  

!-----------------------------------------------------------------    
! Find corresponding land sites from the imogen grid 'sgind'.
!-----------------------------------------------------------------
  CALL GET_IMOGEN_MAP(imogenOrderFile)
    
!-----------------------------------------------------------------
! Allocate imogen arrays
!-----------------------------------------------------------------
  CALL allocate_arrays( 'init_imogen' )
  
!-----------------------------------------------------------------
! Read in IMOGEN namelists
!-----------------------------------------------------------------
  funit = fileUnit( formatAsc )
  CALL openFile(1,.FALSE.,funit,'read',formatAsc,namelistFile,'old')
  READ(funit,ANLG_VALS)
  READ(funit,RUN)
  CALL closeFile(funit,formatAsc)
  
! Weather generator is not available at present
  IF(WGEN) THEN
    PRINT *,'Weather generator not available at present'
    STOP
!        CALL WGEN_MAIN(SEED_WG,DIR_CLIM,PRECIP_WG,TMIN_WG,
!     &            TMAX_WG,SW_WG,RH15M_WG)
  ENDIF
  
!-----------------------------------------------------------------
! Check that the configurations proposed are valid.
! This contains the runs that are allowed on the "decision" tree
! given under directory "plots" and called "imogen.jpg"
!-----------------------------------------------------------------
  CALL IMOGEN_CHECK(                                              &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,            &
    OCEAN_FEED,ANLG,ANOM                                          &
  )

!-----------------------------------------------------------------
! This subroutine only allows runs that have been checked.
! Other runs can be made (editing out the "STOP" command),
! but at one's peril!
!-----------------------------------------------------------------
  CALL IMOGEN_CONFIRMED_RUN(                                      &
    C_EMISSIONS,INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,            &
    OCEAN_FEED,.TRUE.,WGEN,ANLG,ANOM                              &
  )
  
!-----------------------------------------------------------------
! Set up the various initial conditions. 
!-----------------------------------------------------------------
  CO2_PPMV = CO2_INIT_PPMV
    
  IF( INCLUDE_CO2 ) THEN
    CO2_CHANGE_PPMV = 0.0
  ENDIF
        
  IF( ANLG ) THEN
    DTEMP_O(:)=0.0

    IF(INCLUDE_CO2 .AND. OCEAN_FEED .AND. C_EMISSIONS) THEN
      FA_OCEAN(:)=0.0
    ENDIF 
  ENDIF 
 
! Initiate seeding values for subdaily rainfall
  SEED_RAIN(1) = 9465
  SEED_RAIN(2) = 1484
  SEED_RAIN(3) = 3358
  SEED_RAIN(4) = 8350
! Initiate seeding values for the weather generator 
  IF(WGEN) THEN
    SEED_WG(1) = 5810
    SEED_WG(2) = 5575
    SEED_WG(3) = 5817
    SEED_WG(4) = 9119
  ENDIF

!------------------------------------------------------------------------
! Now get the emissions/CO2 concentrations.
! At present only coded for reading in a file of emission/CO2 concentrat
! imogen projects of the future climate, with carbon cycle feedbacks) or
! in a file of CO2 concentrations for "Hydrology 20th Century" simulatio
!------------------------------------------------------------------------
  IF(C_EMISSIONS .AND. INCLUDE_CO2 .AND. ANOM .AND. ANLG) THEN
    funit = fileUnit( formatAsc )
    CALL openFile(1,.FALSE.,funit,'read',formatAsc,FILE_SCEN_EMITS,'old')
    DO N=1,NYR_EMISS
      READ(funit,*) YR_EMISS(N),C_EMISS(N)
    ENDDO
    CALL closeFile(funit,formatAsc)
  ENDIF

!-----------------------------------------------------------------------
! Read in monthly control climate data.
!-----------------------------------------------------------------------
  DO J=1,MM
    FILE_CLIM = TRIM(DIR_CLIM) // DRIVE_MONTH(J)
        
    funit = fileUnit( formatAsc )
    CALL openFile(1,.FALSE.,funit,'read',formatAsc,FILE_CLIM,'old')
    READ(funit,*) LONGMIN_CLIM,LATMIN_CLIM,LONGMAX_CLIM,          &
                  LATMAX_CLIM
                  
    DO I=1,n_imogen_land
      IF (SGINDINV(I) > 0) THEN
        L = SGINDINV(I)
        READ(funit,*) LONG(L),LAT(L),T_CLIM(L,J),                 &
                      RH15M_CLIM(L,J),UWIND_CLIM(L,J),            &
                      VWIND_CLIM(L,J),LW_CLIM(L,J),SW_CLIM(L,J),  &
                      DTEMP_CLIM(L,J),RAINFALL_CLIM(L,J),         &
                      SNOWFALL_CLIM(L,J),PSTAR_HA_CLIM(L,J),      &
                      F_WET_CLIM(L,J)
      ELSE
        READ(funit,*)
      ENDIF
    ENDDO
    CALL closeFile(funit,formatAsc)
  ENDDO


  RETURN

END SUBROUTINE init_imogen
