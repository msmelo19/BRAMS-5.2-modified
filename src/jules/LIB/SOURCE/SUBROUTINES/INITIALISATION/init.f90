!#######################################################################
!#######################################################################
! subroutine init
! Subroutine to read input data and initialise the model.

  SUBROUTINE init(mxp,myp,npatch,nzg,nstyp,dtlong,iyear1,imonth1,idate1,itime1,timmax &
                  ,patch_area,glon,glat,veg_fracarea,leaf_class,stom_resist,soil_water,soil_text,slmsts &
                  ,veg_lai,temp2,soil_energy1,soil_energy2  &
                  ,veg_ht,slcpd,nvtyp_teb  &
                  ,slbs,slcons,slden,Wwilt,PHIsat,alfa_vG  &
                  ,MYNUM,hfilout,hfilin,runtype,zts2J,slz)

!<in sfclyr_jules>    USE init_grid_mod, ONLY :  &
!  imported procedures
!<in sfclyr_jules>       init_grid

  USE initial_mod, ONLY :  &
!  imported procedures
     dump_io,init_ic

  USE inout, ONLY :  &
!  imported scalars with intent(in)
     dumpFreq,stdIn,formatAsc  &
!  imported scalars with intent(out)
    ,jinUnit, echo

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     spinUp

  USE veg_io_vars, ONLY :  &
!  imported scalars with intent(in)
     vegVaryT

  USE file_utils, ONLY : fileUnit,openFile,closeFile

  !DSM{
  USE ancil_info, ONLY :  &!
!  imported arrays with intent(out)
     z1_uv,z1_tq,land_pts,land_index,row_length,frac

   USE nstypes, ONLY : npft

  USE p_s_parms, ONLY : z0_tile
  
  USE pftparm
  USE prognostics, ONLY : canht_ft, lai
  
  !DSM}

! Import procedures needed to read arguments
! This line is only required for the NAG compiler
!  USE F90_UNIX_ENV, ONLY : IARGC,GETARG

!-----------------------------------------------------------------------

  IMPLICIT NONE

! Declare local variables
  CHARACTER(len=250) :: jinFilename

!-----------------------------------------------------------------------
      integer, intent(in)  :: mxp,myp,npatch,nzg,nstyp,iyear1,imonth1,idate1,itime1
      real,    intent(in)  :: dtlong,timmax
      
      real,    intent(in)  :: patch_area(mxp,myp,npatch),glon(mxp,myp)  &
                             ,glat(mxp,myp),veg_fracarea(mxp,myp,npatch) &
                             ,leaf_class(mxp,myp,npatch),stom_resist(mxp,myp,npatch)    &
                             , soil_water(nzg,mxp,myp) &
                             ,soil_text(nzg,mxp,myp),slmsts(nstyp)              &
                             ,veg_lai(mxp,myp,npatch)    &
                             ,temp2(mxp,myp),soil_energy1(nzg,mxp,myp),soil_energy2(nzg,mxp,myp) &
                             ,slz(nzg)
      
      REAL, INTENT(IN)     :: zts2J(mxp,myp) 
      
      LOGICAL              :: tem_in
      
      REAL :: sst(mxp,myp)!DSM
      INTEGER                             :: nvtyp_teb
      REAL, DIMENSION (nstyp+nvtyp_teb)   :: veg_ht 
      REAL, dimension(nstyp)              :: slbs,slcpd,slcons,slden,Wwilt,PHIsat,alfa_vG

      INTEGER :: MYNUM
      CHARACTER(len=256) :: hfilout,hfilin
      CHARACTER(len=16)  :: runtype
      
      INTEGER   :: i,j,n,l  !DSM


! To enable passing of the control file as an argument rather than
! piping standard input, un-comment the code below. This makes
! attaching certain debuggers easier, particularly those integrated with
! IDES (e.g. Photran for Eclipse).
! This code causes problems compiling with the Oracle Studio (Sun)
! and Portland Group compilers, so is disabled by default.
!
! This code is written with the expectation that either:
!  * One argument is supplied that contains a path to a formatted file
!  * No argument is given and we read from standard input that should be
!    redirected to a formatted file
!-----------------------------------------------------------------------
!  IF ( IARGC() > 0 ) THEN
! The first argument will be the filename
!    CALL GETARG(1, jinFilename)
!    WRITE(*,"(50('-'),/,a)")'Reading model control file from file '  &
!                         // TRIM(jinFilename)                        &
!                         // '...'

! Get a unit for the control file specified and open it
!    jinUnit = fileUnit(formatAsc)
!    CALL openFile(1,.FALSE.,jinUnit,'READ',formatAsc,jinFilename,'old')
!  ELSE
    !DSM WRITE(*,*)'If this program appears to be waiting for input,'
    !DSM WRITE(*,*)'that''s because it is intended to be used with input'
    !DSM WRITE(*,*)'from a file, and expects a particular format.'
    !DSM WRITE(*,*)'Use with standard input redirected,'
    !DSM WRITE(*,*)'e.g. xjules < jules.in'
    !DSM WRITE(*,"(50('-'),/,a)")'Reading model control file from stdin...'

    !<in sfclyr_jules>  jinUnit = stdIn
!  ENDIF

!-----------------------------------------------------------------------
! Read model options and misc other values.
!-----------------------------------------------------------------------
  !<in sfclyr_jules> CALL init_opts(mxp,myp)

!-----------------------------------------------------------------------
! Read date, time and location information
!-----------------------------------------------------------------------
  !<in sfclyr_jules> CALL init_time(....)

!-----------------------------------------------------------------------
! Read details of model grid and allocate arrays.
!-----------------------------------------------------------------------
  !<in sfclyr_jules> CALL init_grid(....)

!-----------------------------------------------------------------------
! Read fractional coverages.
!-----------------------------------------------------------------------
  CALL init_frac(mxp,myp,npatch,patch_area,veg_fracarea,leaf_class)

!-----------------------------------------------------------------------
! Read soil layer thicknesses, and hydraulic and thermal characteristics.
!-----------------------------------------------------------------------
  CALL init_soil(mxp,myp,npatch,nzg,soil_text &
                 ,slbs,slcpd,slcons,slden,Wwilt,PHIsat,alfa_vG,slmsts,nstyp,slz)

!-----------------------------------------------------------------------
! Read TOPMODEL parameters.
!-----------------------------------------------------------------------
  CALL init_top

!-----------------------------------------------------------------------
! Read PDM parameters.
!-----------------------------------------------------------------------
  CALL init_pdm

!-----------------------------------------------------------------------
! Read surface heights
!-----------------------------------------------------------------------.
  CALL init_hgt

!-----------------------------------------------------------------------
! Read PFT parameters.
!-----------------------------------------------------------------------
  CALL init_veg

!-----------------------------------------------------------------------
! Read parameters for non-veg types.
!-----------------------------------------------------------------------
  CALL init_nonveg

!-----------------------------------------------------------------------
! Read urban parameters and calculate disp and ztm from MacDonald (1998)
! formulation
!-----------------------------------------------------------------------
  CALL init_urban

!-----------------------------------------------------------------------
! Read snow parameters.
!-----------------------------------------------------------------------
  CALL init_snow

!-----------------------------------------------------------------------
! Read TRIFFID parameters.
!-----------------------------------------------------------------------
  CALL init_trif

!-----------------------------------------------------------------------
! Read agricultural fraction.
!-----------------------------------------------------------------------
  CALL init_agric

!-----------------------------------------------------------------------
! Read miscellaneous surface and carbon/veg parameters.
!-----------------------------------------------------------------------
  CALL init_misc

!-----------------------------------------------------------------------
! Read runoff routing parameters.
!-----------------------------------------------------------------------
  CALL init_route
  
!-----------------------------------------------------------------------
! Read IMOGEN parameters
!-----------------------------------------------------------------------
  CALL init_imogen

!-----------------------------------------------------------------------
! Read details of meteorological driving data.
!-----------------------------------------------------------------------
  CALL init_drive
  
!-----------------------------------------------------------------------
! Read the initial state.
!-----------------------------------------------------------------------
  CALL init_ic(mxp,myp,npatch,nzg,nstyp,patch_area   &
               ,veg_fracarea,leaf_class,stom_resist,soil_water,soil_text &
               ,slmsts,veg_lai,temp2,soil_energy2 & 
               ,veg_ht,slcpd,nvtyp_teb,MYNUM,hfilout,hfilin,runtype)



  
  !DSM{
  !--- - Colocando z1_uv e z1_tq em funcao da altura do segundo nivel do CCATT-BRAMS {  
  DO l=1,land_pts
    j = ( land_index(l)-1 ) / row_length + 1
    i = land_index(l) - (j-1)*row_length
    

    z1_uv(i,j)=0.0
    DO n=1,npft
      !--- Calculando o dz0v_dh em funcao da formulacao do LEAF ---
      IF (dz0v_dh(l,n)<0.) THEN  ! se configurar no jules.in com dz0v_dh negativo
         dz0v_dh(l,n)=1.-0.91*exp(-0.0075*lai(l,n))
      ENDIF

      if (frac(l,n)<0.0 .or. frac(l,n)>1.0) STOP "frac no initialized..."
      !DSM z1_uv(i,j) = z1_uv(i,j) + frac(l,n) * (canht_ft(l,n) - dz0v_dh(n)/0.264)
      z1_uv(i,j) = z1_uv(i,j) + frac(l,n) * (canht_ft(l,n) - dz0v_dh(l,n)/0.264)
      !write(*,'(a,3i,10f8.2)') "SDFG=",i,j,n, z1_uv(i,j), canht_ft(l,n), frac(l,n), dz0v_dh(n)
    ENDDO
    
    z1_uv(i,j)=zts2J(i,j) - z1_uv(i,j)
  ENDDO
  
  !print*,"SSFSFF=",z1_uv
  z1_tq(:,:)=z1_uv(:,:)
  !DSM}

!-----------------------------------------------------------------------
! Setup the output.
!-----------------------------------------------------------------------
  CALL init_out(MYNUM,hfilout)

  if (echo) WRITE(*,"(60('-'),/,a)")'Finished reading control file.'

!-----------------------------------------------------------------------
! Temporary (i.e. for this version) initialisation of variables.
!-----------------------------------------------------------------------
  !{DSM
  sst(:,:)= (soil_energy1(nzg,:,:)-334000)/4186 + 273.15
  !DSM}
  CALL INIT_VARS_TMP(sst,mxp,myp)

!-----------------------------------------------------------------------
! Set index arrays and initialise other variables.
!-----------------------------------------------------------------------
  CALL INIT_PARMS

!-----------------------------------------------------------------------
! Save initial state if spinning up. Arrays are allocated here.
!-----------------------------------------------------------------------
  IF ( spinUp ) CALL spin_check

!-----------------------------------------------------------------------
! Write an initial dump (if requested).
! If veg fields vary with time, wait until these are definitely updated.
!-----------------------------------------------------------------------
  IF ( .NOT.vegVaryT .AND. dumpFreq>1 ) CALL dump_io( .TRUE., dumpTypeArg='init' )

!-----------------------------------------------------------------------
! Deallocate space that was only needed during intialisation.
!-----------------------------------------------------------------------
  CALL deallocate_arrays( 'init_end' )

!-----------------------------------------------------------------------
! Close the jin file if we were reading from a file other than stdin
!-----------------------------------------------------------------------
  IF ( jinUnit /= stdIn ) THEN
    CALL closeFile(jinUnit, formatAsc)
! Set jinUnit to a value that cannot possibly be a file unit to avoid future collisions
    jinUnit = -1
  END IF

  if (echo) WRITE(*,"(70('#'),/,a,/,70('#'))") 'Initialisation complete. Start of run.'

  END SUBROUTINE init

!#######################################################################
!#######################################################################
