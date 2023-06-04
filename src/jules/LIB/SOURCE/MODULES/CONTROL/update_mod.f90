! module update_mod
! Module containing procedures relevant to updating of forcing data (meteorology and vegetation).
! The meteorology vars are implemented in a relatively flexible way and using
! code similar to that used elsewhere for other variables. The same cannot be
! said for the veg vars.

!-------------------------------------------------------------------------------

  MODULE update_mod

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType,ncTypeDrive

  CONTAINS

!################################################################################
!################################################################################

! subroutine drive_update
! Internal procedure in update_mod.
! Update driving data. This may involve reading data from file.


  !DSM SUBROUTINE drive_update( tstep,next,istepArg )
  SUBROUTINE drive_update( mxp,myp,npatch,nzg,rshort_diffuse,rshort,rlong,temps,ups,vps,pps,rvs,pcpgl   &
                          ,tstep,next,istepArg )  !DSM

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     nx=>row_length,ny=>rows,land_pts,land_index

  USE drive_io_vars, ONLY :  &
!  imported scalars with intent(in)
      diffFracConst,driveDataPer,driveEndTime,driveFilePer  &
     ,driveFileTemplate,driveResetStep,driveResetStepPrev  &
     ,driveTemplateDate,driveTemplateTime,driveTemplateUnits,ioPrecipType  &
     ,io_rad_type,tForCRain,tForSnow  &
!     ,iposExtra  & !  uncomment this line to access the optional "extra" variables
     ,iposDiffRad,iposLWD,iposLWN,iposPrecip,iposPrecipCR,iposPrecipCS  &
     ,iposPrecipLR,iposPrecipLS,iposPrecipTR,iposPrecipTS,iposPstar,iposQ,iposRN &
     ,iposSubSurfRoff,iposSurfRoff,iposSWD,iposSWN,iposT,iposU,iposV,iposWind  &
     ,iposOzone  &
!     ,pstarConst,tAmplConst  &
     ,useWgen  &
!  imported scalars with intent(inout)
     ,driveDataStep,driveDate,driveDateInit,driveFile,driveFileStep  &
     ,driveTemplateT,driveTime,driveTimeInit,notNextDrive  &
!  imported arrays with intent(in)
     ,driveVarPos,driveVarUse  &
!  imported arrays with intent(inout)
     ,driveData,driveFileDate,driveFileName,driveFileTime

  USE fluxes, ONLY :  &
!  imported arrays with intent(inout)
     sub_surf_roff,surf_roff

  USE forcing, ONLY : con_rain,con_snow,pstar,qw_1  &
                     ,ls_rain,ls_snow,lw_down,sw_down,tl_1  &
                     ,u_0,u_1,v_0,v_1

  USE ozone_vars, ONLY : o3

  USE grid_utils, ONLY :  &
!  imported procedures
     getLandValues

  USE misc_utils, ONLY :  &
!  imported procedures
     find_file

  USE surf_param, ONLY : &
!  imported arrays with intent(out)
     diff_frac

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_360,l_point_data,routeOnly,l_imogen

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
     timeStep,date,time

  USE time_mod, ONLY :  &
!  imported procedures
     timeDate,dateToBits

  USE u_v_grid, ONLY :  &
!  imported arrays with intent(inout)
     u_0_p,u_1_p,v_0_p,v_1_p
     
  USE timeConst, ONLY : &
    iSecInDay
    
  USE imogen_drive_vars, ONLY : &
    PSTAR_OUT,WIND_OUT,CONV_RAIN_OUT,CONV_SNOW_OUT,LS_RAIN_OUT,  &
    LS_SNOW_OUT,SW_OUT,LW_OUT,QHUM_OUT,T_OUT

!-------------------------------------------------------------------------------

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) ::  tstep  !   timestep number.
!                                    0 indicates a call during initialisation.

  LOGICAL, INTENT(in) ::  next  !  T means assume that the requested data are
!             held in the next time level, either
!             of the current file or the next file if appropriate
!             F means search file details to locate the data

!-------------------------------------------------------------------------------
! Optional scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in), OPTIONAL :: istepArg  !  current value of the time index
!        Generally istep=driveDataStep, in which case istepArg is not provided.
!        However, at initialisation, it may be necessary to use a different value of istep.

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    i,istep,j,k,l       &! loop counters/work
   ,tmpDate,tmpTime     &!  work
   ,dateDay,dateMonth,dateYear,insd  ! Date parts for indexing IMOGEN arrays

  LOGICAL ::  &
    nextData  &!  a local version of next. Values equals that of next, unless
!                      we need to "reset" the time of data.
   ,resetStep  !  TRUE when data time is to be reset, to deal with spin up

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  INTEGER :: precipMask(nx,ny)   !  used to partition precipitation
!       0=(large-scale) snowfall,  1=large-scale rainfall,  2=convective rainfall

  REAL :: diff_rad(nx,ny)   ! diffuse radiation (W m-2)
  REAL :: wrkDrvData(nx,ny) ! Work array for holding data from driveData().


!DSM{
  integer, intent(in)  :: mxp,myp,npatch,nzg
  real,    intent(in)  :: rshort_diffuse(mxp,myp),rshort(mxp,myp),rlong(mxp,myp),temps(mxp,myp)  &
                          ,ups(mxp,myp),vps(mxp,myp)        &    
                          ,pps(mxp,myp),rvs(mxp,myp),pcpgl(mxp,myp)

!DSM}

!--------------------------------------------------------------------------------

! At present this code cannot deal with "special" periods (e.g. monthly data).
! If these are required, veg_update might provide soem useful help as to how to do.
  IF ( driveDataPer < 0 ) THEN
    WRITE(*,*)'ERROR: drive_update: driveDataPer < 0'
    WRITE(*,*)'There''s no code for "special" periods (e.g. monthly data).'
    WRITE(*,*)'Driving data must have a "regular" period (driveDataPer>0).'
    STOP
  ENDIF
  
  
!-------------------------------------------------------------------------------
! If IMOGEN is enabled, just copy the correct timestep of the
! current climatology into the driving variables and we're done
!-------------------------------------------------------------------------------
  IF( l_imogen ) THEN
    CALL dateToBits( date,dateDay,dateMonth,dateYear,l_360,'DateToBits failed' )
! Get the timestep in the day that we are on
    insd = (time / timestep) + 1

    DO l=1,land_pts
      j = (land_index(l) - 1) / nx + 1
      i = land_index(l) - (j-1) * nx
      
      pstar(i,j)    = PSTAR_OUT(l,dateMonth,dateDay,insd)
      u_1(i,j)      = WIND_OUT(l,dateMonth,dateDay,insd)
      v_1(i,j)      = 0.0
      u_0(i,j)      = 0.0
      v_0(i,j)      = 0.0
      con_rain(i,j) =                                             &
        CONV_RAIN_OUT(l,dateMonth,dateDay,insd) / float(iSecInDay)
      con_snow(i,j) =                                             &
        CONV_SNOW_OUT(l,dateMonth,dateDay,insd) / float(iSecInDay)
      ls_rain(i,j)  =                                             &
          LS_RAIN_OUT(l,dateMonth,dateDay,insd) / float(iSecInDay)
      ls_snow(i,j)  =                                             &
          LS_SNOW_OUT(l,dateMonth,dateDay,insd) / float(iSecInDay)
      sw_down(i,j)  = SW_OUT(l,dateMonth,dateDay,insd)
      lw_down(i,j)  = LW_OUT(l,dateMonth,dateDay,insd)
      qw_1(i,j)     = QHUM_OUT(l,dateMonth,dateDay,insd)
      tl_1(i,j)     = T_OUT(l,dateMonth,dateDay,insd)
      u_0_p(i,j)    = u_0(i,j)
      v_0_p(i,j)    = v_0(i,j)
      u_1_p(i,j)    = u_1(i,j)
      v_1_p(i,j)    = v_1(i,j)
    ENDDO
    
    RETURN
  ENDIF
  

! Initialise.
  resetStep = .FALSE.

! Increment counter.
  driveDataStep = driveDataStep + 1

! Deal with resetting time and date for spin up..
! driveResetStep is used to reset driving data for spin up.
! Note that we test on both driveResetStep and driveResetStepPrev, because a test on only
! the former fails if we need to reset data on the first timestep of a new cycle of spin up.
! In that case, newTime has already recalculated driveResetStep before this routine has
! had a chance to reset the data - but we can catch that case using the previous value (driveResetStepPrev).
!XX  Is there a problem here if we have tstep=0 (to indicate "special" call, but this also happens to be driveResetPer?)
  IF ( tstep==driveResetStep .OR. tstep==driveResetStepPrev ) resetStep = .TRUE.

! Decide if it is time to read new driving data.
  IF ( driveDataStep>driveDataPer .OR. resetStep ) THEN

!   Read new data.

!   Initialise.
    nextData = next

!   If switch set, force code to search for appropriate time.
    IF ( notNextDrive ) THEN
      nextData = .FALSE.
      notNextDrive = .FALSE.  !  reset the switch
    ENDIF

!   Reset counter.
    driveDataStep = 1
!    print*,'Read new data, reset driveDataStep to',driveDataStep

!   Get date and time of the next driving data.
    CALL timeDate( driveTime,driveDate,driveDataPer*NINT(timeStep),'sec',l_360,tmpTime,tmpDate,'drive_update' )
    driveTime = tmpTime
    driveDate = tmpDate
!    print*,'updated driveTime,date=',driveTime,driveDate

    IF ( resetStep ) THEN
!     Resetting time and date of data, for spin up.
!     Reuse initial values.
      driveTime = driveTimeInit
      driveDate = driveDateInit
!      print*,'drive_update: Have reset to use initial values, driveTime,date=',driveTime,driveDate
!     Set switch to trigger a search for these data.
      nextData = .FALSE.
    ENDIF

!   Work out time level and file that hold the required data.
    CALL find_file( driveTime,driveDate,.FALSE.,driveEndTime,nextData,driveDataPer,NINT(timeStep)  &
                   ,driveFile,driveFileStep  &
                   ,driveTemplateT,driveFilePer,driveFileTime,driveFileDate  &
                   ,driveFileName,driveFileTemplate,driveTemplateDate  &
                   ,driveTemplateTime,driveTemplateUnits,'drive_update' )
!-------------------------------------------------------------------------------

!    PRINT*,'To drive_read: Next drive time,date=',driveDate,driveTime
!    PRINT*,'driveFileStep=',driveFileStep,' driveFile=',driveFile

!   Read driving data.
    !DSM CALL drive_read
    CALL drive_read(mxp,myp,npatch,nzg,rshort_diffuse,rshort,rlong,temps,ups,vps,pps,rvs,pcpgl)
!    print*,'+ drive_read'

    IF ( .NOT. useWGen ) THEN
!     Interpolate data in time.
!     This calculates values for all timesteps until the next data are read.
      CALL drive_tinterp
!      print*,'+ drive_tinterp'
    ELSE
!     Call the weather generator.
!     Set values that are not input in this version.
!!      tAmpl(:,:) = tAmplConst
!!      pstarVal(:,:) = pstarConst
!      print*,'maxval(driveDataIn)=',maxval(driveDataIn)

!      print*,'driveDataIn(driveVarPosIn(iposT),1,1,0)=',driveDataIn(driveVarPosIn(iposT),1,1,0)
!      print*,'driveDataIn(driveVarPosIn(iposPrecip),1,1,0)=',driveDataIn(driveVarPosIn(iposPrecip),1,1,0)

!     CALL day_calc( nx,ny,driveDataPer,timestep,latitude,longitude  &
!               ,tAmpl  &
!               ,driveDataIn(driveVarPosIn(iposLWD),:,:,0)     &
!               ,driveDataIn(driveVarPosIn(iposPrecip),:,:,0)  &
!               ,pstarVal  &
!               ,driveDataIn(driveVarPosIn(iposQ),:,:,0)  &
!               ,driveDataIn(driveVarPosIn(iposSWD),:,:,0)     &
!               ,driveDataIn(driveVarPosIn(iposT),:,:,0)       &
!               ,driveDataIn(driveVarPosIn(iposWind),:,:,0)    &
!               ,driveData(driveVarPos(iposPrecipCR),:,:,:)    &
!               ,driveData(driveVarPos(iposPrecipLR),:,:,:)    &
!               ,driveData(driveVarPos(iposPrecipTS),:,:,:)    &
!               ,driveData(driveVarPos(iposLWD),:,:,:)         &
!               ,driveData(driveVarPos(iposPstar),:,:,:)       &
!               ,driveData(driveVarPos(iposQ),:,:,:)           &
!               ,driveData(driveVarPos(iposSWD),:,:,:)         &
!               ,driveData(driveVarPos(iposT),:,:,:)           &
!               ,driveData(driveVarPos(iposWind),:,:,:) )

!      print*,'+DAY driveData(driveVarPos(iposT),1,1,:)=',driveData(driveVarPos(iposT),1,1,:)
!      print*,'+DAY CR,LR,TS(1,1)=',driveData(driveVarPos(iposPrecipCR),1,1,1),  &
!         driveData(driveVarPos(iposPrecipLR),1,1,1),driveData(driveVarPos(iposPrecipTS),1,1,1)
    ENDIF

  ENDIF  !  driveDataStep

!-------------------------------------------------------------------------------

! Load the driving data into the final variables.

! Find out what time level of driveData to use.
  istep = driveDataStep
  IF ( PRESENT(istepArg) ) istep=istepArg
  IF ( istep<LBOUND(driveData,4) .OR. istep>UBOUND(driveData,4)  ) THEN
    WRITE(*,*)'ERROR: drive_update: istep is out of range.'
    WRITE(*,*)'Value=',istep,' allowed range=',LBOUND(driveData,4),' to ',UBOUND(driveData,4)
    IF ( PRESENT(istepArg) ) THEN
      WRITE(*,*)'The most likely way for this to happen is if the actual argument corresponding'
      WRITE(*,*)'to the optional dummy argument istepArg is wrongly specified.'
    ENDIF
    STOP
  ENDIF

!  print*,'LOADING drive data with driveDataStep=',driveDataStep,' istep=',istep
!  print*,'driveData(driveVarPos(iposT),:,:,:)=',driveData(driveVarPos(iposT),:,:,:)

!-------------------------------------------------------------------------------

! If a variable is required (driveVarUse=TRUE), load it.
! Some of the testing below is currently unnecessary, since earlier tests have
! ensured that a variable is present if necessary, e.g. if ioPrecipType=2 we must
! have driveVarUse(iposPrecipTL)=TRUE.
! But we test anyway, for moderate robustness!
! Also prevent values going out of range because of interpolation - this is only a
! possible problem with certain interpolation options (b,c,f), but is always enforced.

!   Variables that must be held in driveData - can load directly.
    IF ( driveVarUse(iposPstar) ) pstar(:,:) = driveData(driveVarPos(iposPstar),:,:,istep)

    IF ( driveVarUse(iposQ) ) THEN
      qw_1(:,:) = driveData(driveVarPos(iposQ),:,:,istep)
      qw_1(:,:) = MAX( qw_1, 0.0 )
    ENDIF

    IF ( driveVarUse(iposT) ) tl_1(:,:) = driveData(driveVarPos(iposT),:,:,istep)

    IF ( driveVarUse(iposU) ) THEN
      IF ( driveVarPos(iposU) > 0 ) THEN
        u_1(:,:) = driveData(driveVarPos(iposU),:,:,istep)
      ELSE
!       Assume that ioWindSpeed=T, so we can use total wind speed.
        u_1(:,:) = driveData(driveVarPos(iposWind),:,:,istep)
      ENDIF
    ENDIF

!-------------------------------------------------------------------------------
!   Optional variables (can be read in, or derived from another variable).

!   Precipitation variables

!   To ensure correct partitioning of precip with respect to temperature, create a mask field.
!   Assuming that temperature is available.
!   Only needed for some input choices (ioPrecipType=1 or 2).
    IF ( ioPrecipType==1 .OR. ioPrecipType==2 ) THEN
      precipMask(:,:) = 1    !   large-scale rainfall
      IF ( ioPrecipType == 1 ) WHERE ( tl_1(:,:) <= tForSnow ) precipMask(:,:) = 0   !  large-scale snowfall
      IF ( .NOT. l_point_data ) WHERE ( tl_1(:,:) >= tForCRain )  &
                      precipMask(:,:) = 2   !  convective rainfall
    ENDIF

!   Initialise.
    con_rain(:,:) = 0.0
    con_snow(:,:) = 0.0
    ls_rain(:,:) = 0.0
    ls_snow(:,:) = 0.0

!-------------------------------------------------------------------------------
!#### NOTE ####
!  Do not replace lines such as:
!     wrkDrvData(:,:) = driveData(driveVarPos(iposPrecip),:,:,istep)
!     where ( precipMask(:,:) == 1 ) ls_rain(:,:) = wrkDrvData(:,:)
!  with the apparently simpler equivalent
!    where ( precipMask(:,:) == 1 ) ls_rain(:,:) =  driveData(driveVarPos(iposPrecip),:,:,istep)
!  as the latter gives errors with some compilers.
!-------------------------------------------------------------------------------

!   Convective rainfall.
    IF ( driveVarUse(iposPrecipCR) ) THEN
!     If using point data, leave convective as zero.
      IF ( .NOT. l_point_data ) THEN
        IF ( ioPrecipType == 1 ) THEN
!         Total precip is available.
          wrkDrvData(:,:) = driveData(driveVarPos(iposPrecip),:,:,istep)
          WHERE ( precipMask(:,:) == 2 ) con_rain(:,:) = wrkDrvData(:,:)
        ELSEIF ( ioPrecipType == 2 ) THEN
!         Total rainfall is available.
          wrkDrvData(:,:) = driveData(driveVarPos(iposPrecipTR),:,:,istep)
          WHERE ( precipMask(:,:) == 2 ) con_rain(:,:) = wrkDrvData(:,:)
        ELSE
!         Conv rainfall is available.
          con_rain(:,:) = driveData(driveVarPos(iposPrecipCR),:,:,istep)
        ENDIF
        con_rain(:,:) = MAX( con_rain, 0.0 )
      ENDIF
    ENDIF

!   Convective snowfall.
    IF ( driveVarUse(iposPrecipCS) ) THEN
!     Only non-zero with ioPrecipType=4.
      IF ( ioPrecipType == 4 ) THEN
!       Convective snowfall was input.
        con_snow(:,:) = driveData(driveVarPos(iposPrecipCS),:,:,istep)
      ENDIF
      con_snow(:,:) = MAX( con_snow, 0.0 )
    ENDIF

!   Large-scale rainfall.
    IF ( driveVarUse(iposPrecipLR) ) THEN
      IF ( ioPrecipType == 1 ) THEN
!       Total precipitation was input.
        wrkDrvData(:,:) = driveData(driveVarPos(iposPrecip),:,:,istep)
        WHERE ( precipMask(:,:) == 1 ) ls_rain(:,:) = wrkDrvData(:,:)
      ELSEIF( ioPrecipType == 2 ) THEN
!       Total rainfall was input.
        wrkDrvData(:,:) = driveData(driveVarPos(iposPrecipTR),:,:,istep)
        WHERE ( precipMask(:,:) == 1 ) ls_rain(:,:) = wrkDrvData(:,:)
      ELSE
!       Large-scale rainfall is available.
        ls_rain(:,:) = driveData(driveVarPos(iposPrecipLR),:,:,istep)
      ENDIF
      ls_rain(:,:) = MAX( ls_rain, 0.0 )
    ENDIF

!   Large-scale snowfall.
    IF ( driveVarUse(iposPrecipLS) ) THEN
      SELECT CASE ( ioPrecipType )
        CASE ( 1 )
!         Total precipitation was input.
          wrkDrvData(:,:) = driveData(driveVarPos(iposPrecip),:,:,istep)
          WHERE ( precipMask(:,:) == 0 ) ls_snow(:,:) = wrkDrvData(:,:)
        CASE ( 2, 5 )
!         Total snowfall is available. Assume all of this is large-scale.
          ls_snow(:,:) = driveData(driveVarPos(iposPrecipTS),:,:,istep)
        CASE ( 3 )
!         Total snowfall was input. Assume all of this is large-scale.
          ls_snow(:,:) = driveData(driveVarPos(iposPrecipTS),:,:,istep)
        CASE default
!         Large-scale snowfall is available.
          ls_snow(:,:) = driveData(driveVarPos(iposPrecipLS),:,:,istep)
      END SELECT
      ls_snow(:,:) = MAX( ls_snow, 0.0 )
    ENDIF

!-------------------------------------------------------------------------------
!   Radiation variables.
!   Here we load any input (net flux) variables into the downward flux variables.
!   If the inputs are net fluxes, the calculation of downward fluxes is currently done in CONTROL.

!   Longwave.
    IF ( driveVarUse(iposLWD) ) THEN
      IF ( io_rad_type == 1 ) THEN
!       Downward fluxes are provided.
        lw_down(:,:) = driveData(driveVarPos(iposLWD),:,:,istep)
      ELSEIF ( io_rad_type == 2 ) THEN
!       Net downward (all-wavelength) and downward shortwave fluxes are provided.
!       Hold all-wavelength flux in lw_down.
        lw_down(:,:) = driveData(driveVarPos(iposRN),:,:,istep)
      ELSEIF ( io_rad_type == 3 ) THEN
!       Net downward fluxes are provided.
!       Hold net fluxes in downward flux variables.
        lw_down(:,:) = driveData(driveVarPos(iposLWN),:,:,istep)
      ENDIF
    ENDIF

!   Shortwave.
    IF ( driveVarUse(iposSWD) ) THEN
      IF ( io_rad_type == 1 ) THEN
!       Downward fluxes are provided.
        sw_down(:,:) = driveData(driveVarPos(iposSWD),:,:,istep)
      ELSEIF ( io_rad_type == 2 ) THEN
!       Net downward (all-wavelength) and downward shortwave fluxes are provided.
        sw_down(:,:) = driveData(driveVarPos(iposSWD),:,:,istep)
      ELSEIF ( io_rad_type == 3 ) THEN
!       Net downward fluxes are provided.
!       Hold net fluxes in downward flux variables.
        sw_down(:,:) = driveData(driveVarPos(iposSWN),:,:,istep)
      ENDIF
      sw_down(:,:) = MAX( sw_down, 0.0 )
    ENDIF

!   Downward radiation should be >=0....but this is not checked because
!   it's bit of a faff to establish whether we are holding downward or net fluxes at present,
!   and anyway only solar is likely to ever be <0!

!-------------------------------------------------------------------------------
!   Fraction of diffuse radiation.
!-------------------------------------------------------------------------------
    IF ( driveVarUse(iposDiffRad) ) THEN
!     Note that driveVarUse(iposDiffRad)=T is equivalent to useDiffRad=T.
!     Diffuse radiation has been input. Now calculate fraction.
!     For simplicity this option is only allowed if io_rad_type=1 or 2, meaning
!     we already know sw_down.

!     First load diffuse radiation.
      diff_rad(:,:) = driveData(driveVarPos(iposDiffRad),:,:,istep)

      k = 0
      DO i=1,nx
        DO j=1,ny
          k = k + 1
          IF ( sw_down(i,j) > 1.0 ) THEN
            diff_frac(k) = diff_rad(i,j) / sw_down(i,j)
            diff_frac(k) = MIN( 1.0, diff_frac(k) )
          ELSE
            diff_frac(k) = 0.0
          ENDIF
        ENDDO
      ENDDO
    ELSE
!     Use a constant fraction.
      diff_frac(:) = diffFracConst
    ENDIF

!--------------------------------------------------------------------------------
!   Wind variables. u-component was dealt with above.
    IF ( driveVarUse(iposV) ) THEN
      v_1(:,:) = 0.0
      IF ( driveVarPos(iposV) > 0 ) v_1(:,:) = driveData(driveVarPos(iposV),:,:,istep)
    ENDIF

!-------------------------------------------------------------------------------

! Get runoff values.
! These need to be extracted from the full model grid onto land points.
! Prevent values going out of range because of interpolation - this is only a
! possible problem with certain interpolation options.

  IF ( driveVarUse(iposSurfRoff) ) THEN
    surf_roff(:) = getLandValues( driveData(driveVarPos(iposSurfRoff),:,:,istep) )
    surf_roff(:) = MAX( surf_roff(:), 0.0 )
  ENDIF

  IF ( driveVarUse(iposSubSurfRoff) ) THEN
    sub_surf_roff(:) = getLandValues( driveData(driveVarPos(iposSubSurfRoff),:,:,istep) )
    sub_surf_roff(:) = MAX( sub_surf_roff(:), 0.0 )
  ENDIF


! Read the ozone field from driving data only if required
! We initialise the field in any case
  o3(:) = 0.0
  IF ( driveVarUse(iposOzone) ) THEN
! We need to map the ozone field onto land points only
    o3(:) = getLandValues( driveData(driveVarPos(iposOzone),:,:,istep) )
  END IF


!-------------------------------------------------------------------------------
! Get data from optional "extra" variables here - you have to add your own code!
! There are ndriveExtra variables available. They have been interpolated in time
! (if necessary) and the values for this timestep are stored in
! driveData(i,:,:,istep), where i=driveVarPos(iposExtra),
! driveVarPos(iposExtra+1),...,driveVarPos(iposExtra+ndriveExtra-1).
! For example, if ndriveExtra=2, and the values are to be loaded into new
! variables called co2 and ozone, this could be achieved with the following code:
! i = driveVarPos(iposExtra)
! co2(:,:) = driveData(i,:,:,istep)
! i = i + 1
! ozone(:,:) = driveData(i,:,:,istep)
! iposExtra is available in module drive_io_vars, and can be USEd by
! uncommenting the line at the top of this routine.
! The user has to make the destination variables (e.g. co2 and ozone) available
! to this routine - either by arguments or via modules.
!--------------------------------------------------------------------------------

  IF ( .NOT. routeOnly ) THEN

!-------------------------------------------------------------------------------
!   Set up information on U, V and T grids (assume that att grids are the same)
!-------------------------------------------------------------------------------
    DO I=1,NX
      DO J=1,NY
        U_0_P(I,J)=U_0(I,J)
        V_0_P(I,J)=V_0(I,J)
        U_1_P(I,J)=U_1(I,J)
        V_1_P(I,J)=V_1(I,J)
      ENDDO
    ENDDO
   ENDIF

!  print*,'drive_update final tl_1(1,1)=',tl_1(1,1)
!  print*,'drive_update final qw_1(1,1)=',qw_1(1,1)
!  print*,'drive_update final ls_rain(1,1)==',ls_rain(1,1)

  END SUBROUTINE drive_update
!###############################################################################
!###############################################################################
! subroutine drive_read
! Internal procedure in update_mod.
! Read driving data from file, making sure the correct files are open.

 !DSM SUBROUTINE drive_read( )
 SUBROUTINE drive_read(mxp,myp,npatch,nzg,rshort_diffuse,rshort,rlong,temps,ups,vps,pps,rvs,pcpgl)

  USE drive_io_vars, ONLY :  &
!   imported scalars with intent(in)
     byteSwapDrive,driveFile,driveFileStep,driveFormat,driveTemplateV  &
    ,ndriveUnit,ndriveHeaderField,ndriveHeaderFile,ndriveHeaderTime  &
    ,ndriveVarIn,nfieldDriveFile,noNewLineDrive  &
!   imported arrays with intent(in)
    ,driveFileName,driveTimeIndex,driveUnit,driveFileOnUnit,driveUnitUse  &
    ,driveVarFlag,driveVarNameSDF,driveVarNameUnit,driveVarPosIn,driveVarInSort,driveVarStash &
!   imported arrays with intent(inout)
    ,driveDataIn

  USE file_utils, ONLY : closeFile,irecPrev,openFile
  USE inout, ONLY : formatNc,mapIn,npoints,nxIn,nyIn
  USE misc_utils, ONLY : replaceTemplate   !  procedures
  USE readwrite_mod, ONLY : readvar2d      !  procedures

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    dInPos       &!  position in driveDataIn to use
   ,i            &!  loop counter
   ,ipos         &!  work
   ,ivar         &!  loop counter
   ,nfieldFile   &!  number of fields per time in a file
   ,nlineField   &!  work
   ,readField    &!  number of field to read
   ,readNyIn     &!  y size of input grid
   ,readZ        &!  'z' level to be read from file
   ,t,unit       &!  work
   ,useIndex     !  index in irecPrev

  INTEGER ::  &!  local ARRAYS
    tmpMap(nfieldDriveFile,2)

  REAL ::  &!  local SCALARS
    tmpval(nfieldDriveFile,1)  !  work - used for input if all data are on one line

  CHARACTER(len=LEN(driveFileName)) ::  &!  local SCALARS
    fileName  !  the name of the required file

!DSM{
  integer, intent(in)  :: mxp,myp,npatch,nzg
  real,    intent(in)  ::   rshort_diffuse(mxp,myp),rshort(mxp,myp),rlong(mxp,myp),temps(mxp,myp)   &
                           ,ups(mxp,myp),vps(mxp,myp)          &
                           ,pps(mxp,myp),rvs(mxp,myp),pcpgl(mxp,myp)

!DSM}
!--------------------------------------------------------------------------------

!  PRINT*,'******** in drive_read with driveFileStep=',driveFileStep

!--------------------------------------------------------------------------------
! Move the existing driving data back one time (e.g. so new data become old data).
!--------------------------------------------------------------------------------
  driveDataIn(:,:,:,:) = EOSHIFT( driveDataIn(:,:,:,:),1,dim=4 )

!--------------------------------------------------------------------------------
! Ensure that the correct file is open.
!--------------------------------------------------------------------------------
  fileName = driveFileName(driveFile)
  DO i=1,ndriveUnit
!--------------------------------------------------------------------------------
!   Get file name if using template. 2nd argument (doTime)=.FALSE. to indicate
!   that we are not changing time template strings.
!--------------------------------------------------------------------------------
    IF ( driveTemplateV ) fileName=replaceTemplate( driveFileName(driveFile)  &
                           ,.FALSE.,'drive_read',replString=driveVarNameUnit(i) )
    IF ( driveFileOnUnit(i) /= fileName ) THEN
!     Close any file that is open on this unit, and open next file.
      !DSM CALL openFile( 1,.TRUE.,driveUnit(i),'read',driveFormat,fileName,'old'  &
      !DSM               ,'drive_read',ncTypeDrive )
      driveFileOnUnit(i) = fileName
    ENDIF
  ENDDO

!--------------------------------------------------------------------------------
! Read driving data into the latest time level.
!--------------------------------------------------------------------------------
  t = driveTimeIndex(2)

! NB at present, ascii file fields will be 1,2,3,...

  IF ( .NOT. noNewLineDrive ) THEN

    DO ivar=1,ndriveVarIn
      ipos = driveVarInSort(ivar)     !  location in master list
      dInpos = driveVarPosIn(ipos)    !  location to use in driveDataIn
!     Get unit to use for this variable.
      unit = driveUnit(driveUnitUse(ipos))
!     Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!     but need to keep index within bounds.
      useIndex = unit
      IF ( driveFormat == formatNc ) useIndex = 1
      readZ = 1   !   'z' level to read from file

      !DSM CALL readVar2d( driveFileStep,readZ,driveVarFlag(ipos)  &
      !DSM                ,driveVarStash(ipos),irecPrev(useIndex)  &
      !DSM                ,nfieldDriveFile,ndriveHeaderFile,ndriveHeaderTime  &
      !DSM                ,ndriveHeaderField,0,nxIn,nyIn,unit  &
      !DSM                ,driveVarNameSDF(ipos)  &
      !DSM                ,mapIn(:),(/(i,i=1,npoints)/),driveFormat,driveDataIn(dInPos,:,:,t)  &
      !DSM                ,byteSwapDrive,'drive_read','drive_read',ncTypeDrive )

      !{DSM
      
      !print*,'ipos=',ipos,dInpos ,unit,driveVarFlag(ipos),driveVarStash(ipos),irecPrev(useIndex)
 
      !--- Vide arquivo LIB/SOURCE/SUBROUTINES/INITIALISATION/init_drive.f90 para saber o numero 
      !--- correspondente a cada variavel.
      IF     (driveVarStash(ipos)==0) THEN !--- diff_rad = Diffuse radiation
         driveDataIn(dInPos,:,:,t)=rshort_diffuse(:,:)
      ELSEIF     (driveVarStash(ipos)==1235) THEN !--- sw_down = Downward shortwave radiation
         driveDataIn(dInPos,:,:,t)=rshort(:,:)
      ELSEIF (driveVarStash(ipos)==2207) THEN !--- lw_down = Downward longwave radiation
         driveDataIn(dInPos,:,:,t)=rlong(:,:)
      ELSEIF (driveVarStash(ipos)==5216) THEN !--- precip = Total precipitation rate
         driveDataIn(dInPos,:,:,t)=pcpgl(:,:)
      ELSEIF (driveVarStash(ipos)==3236) THEN !--- t = Air temperature
         driveDataIn(dInPos,:,:,t)=temps(:,:)
      ELSEIF (driveVarStash(ipos)==3225) THEN !--- u = Zonal wind speed
         driveDataIn(dInPos,:,:,t)=ups(:,:)
      ELSEIF (driveVarStash(ipos)==3226) THEN !--- v = Meridional wind speed
         driveDataIn(dInPos,:,:,t)=vps(:,:)
      ELSEIF (driveVarStash(ipos)==1) THEN !--- pstar = Surface pressure
         driveDataIn(dInPos,:,:,t)=pps(:,:)
      ELSEIF (driveVarStash(ipos)==3237) THEN !--- q = Specific humidity
         driveDataIn(dInPos,:,:,t)=rvs(:,:)

!write(*,'(a,8f12.4)') 'wwwwe1====>(5,12) = ',t,rshort(5,12),rlong(5,12),pcpgl(5,12),temps(5,12),ups(5,12),vps(5,12),pps(5,12),rvs(5,12)
!write(*,'(a,8f12.4)') 'wwwwe2====>(11,22) = ',rshort(11,22),rlong(11,22),pcpgl(11,22),temps(11,22),ups(11,22),vps(11,22),pps(5,12),rvs(11,22)
      ENDIF 

      !DSM}

!      print*,'drive_read driveDataIn(ipos,1,1,t)=',driveDataIn(ipos,1,1,t)
    ENDDO  !  ivar

  ELSE

!--------------------------------------------------------------------------------
!   noNewLineDrive.
!   No new line between fields - all fields are arranged across line or
!   line(s) in the ASCII input file. nxIn=nyIn=1.
!   Read nfieldDriveFile values but as one field.
!     Note: In attempt to cope with files that have character timestamps at end of each line, I
!     thought about only reading as many fields as required to get data. That works OK if data are
!     on a single line. However, say data are on 2 lines, but all required fields are on line 1.
!     If we only read line 1, the file is in the wrong position for the next read.
!--------------------------------------------------------------------------------

!   Create fields to replace mapIn.
    tmpMap(:,1) = (/ (i, i=1,nfieldDriveFile) /)
    tmpMap(:,2) = tmpMap(:,1)
!   Set other values needed for reading.
    unit = driveUnit(1)
    readZ = 1         !   "z" level to read from file
    readField = 1     !   read field number 1
    nfieldFile = 1    !   file appears to only have one field per time
    nlineField = 0    !   will not attempt to read ASCII line-by-line
    readNyIn = 1      !   input has ny=1

    CALL readVar2d( driveFileStep,readZ,readField,driveVarStash(1),irecPrev(unit)  &
                   ,nfieldFile,ndriveHeaderFile,ndriveHeaderTime  &
                   ,ndriveHeaderField,nlineField,nfieldDriveFile,readNyIn,unit  &
                   ,driveVarNameSDF(1)  &
                   ,tmpMap(:,1),tmpMap(:,2),driveFormat,tmpval(1:nfieldDriveFile,1:1)  &
                   ,byteSwapDrive,'drive_read','drive_read',ncTypeDrive)

!--------------------------------------------------------------------------------
!   Extract the required fields.
!--------------------------------------------------------------------------------
    DO ivar=1,ndriveVarIn
      ipos = driveVarInSort(ivar)   !   location in master list
      dInpos = driveVarPosIn(ipos)  !  location to use in driveDataIn
      driveDataIn(dInPos,1,1,driveTimeIndex(2)) = tmpVal(driveVarFlag(ipos),1)
    ENDDO

  ENDIF  !  noNewLineDrive
!  print*, 'With new data driveDataIn(ivar=5)=',driveDataIn(5,1,1,:)
!  print*,'end drive_read'

  END SUBROUTINE drive_read
!################################################################################
!################################################################################
!################################################################################
! subroutine drive_tinterp
! Internal procedure in update_mod.
! Given driving data, calls tinterp to interpolate in time onto model timesteps.

  SUBROUTINE drive_tinterp()

  USE ancil_info, ONLY :  &
!   imported scalars with intent(in)
     nx=>row_length,ny=>rows

  USE drive_io_vars, ONLY :  &
!   imported scalars with intent(in)
     driveDataPer,ndriveVarIn  &
!   imported arrays with intent(in)
    ,driveDataIn,driveTimeIndex,driveVarInterp,driveVarPos,driveVarPosIn,driveVarInSort &
!   imported arrays with intent(out)
    ,driveData
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    dpos         &!  position in driveData to use
   ,dInPos       &!  position in driveDataIn to use
   ,ipos,ivar,ix,iy,t1,t2      !  work

  CHARACTER(len=1) ::  &!  local scalars
    disAggFlag      !  work: disaggregation flag passed to interpolation routine

!-------------------------------------------------------------------------------
  disAggFlag = '0'   !  disaggregation flag, always '0' for now
  t1 = driveTimeIndex(1)
  t2 = driveTimeIndex(2)

! Call tinterp with 1st and 2nd args (ntOut and nt)=driveDataPer so that values
! for all timesteps in the period of the data (driveDataPer) are calculated.
! 3rd arg (t) =1 to fill all values starting at first element.

  DO iy=1,ny
    DO ix=1,nx
      DO ivar=1,ndriveVarIn
        ipos = driveVarInSort(ivar)   !   location in master list
        dInpos = driveVarPosIn(ipos)  !  location to use in driveDataIn
        dpos = driveVarPos(ipos)  !  location to use in driveData
!        print*, 'drive_tinterp ivar=',ivar,' to interp with inval=',driveDataIn(dInPos,ix,iy,t1:t2)
        CALL tinterp( driveDataPer,driveDataPer,1,driveTimeIndex,driveDataIn(dInPos,ix,iy,t1:t2)  &
                         ,disAggFlag,driveVarInterp(ipos),driveData(dpos,ix,iy,:) )
!        print*, 'drive_tinterp: After tinterp, outval=',driveData(dpos,ix,iy,:)
      ENDDO
    ENDDO
  ENDDO

  END SUBROUTINE drive_tinterp

!###############################################################################
!###############################################################################
!###############################################################################
! subroutine tinterp
! Internal procedure in update_mod.
! Do temporal interpolation or disaggregation of input data.
! Output is a vector, but can be of length 1.
! Based on drv_finterp from GSWP2, by Dirmeyer, Guo and Lohmann.
! The "conserving interpolation" options ( 'b','c', and 'f') are, apparently,
! simplified versions of Sheng and Zwiers, 1998, Climate Dynamics, 14: 609-613.

! Notes:
! The available data has a period (frequency) of nt timesteps.
! We can request values at ntOut timesteps within this period, where
! ntOut is either 1 or nt.
! ntOut=1 means one value is requested. This is the value at timestep t (in the range 1:nt).
! ntOut=nt means that values are requested for all nt timesteps in the data period. In this case,
!     t should be 1 (since we are starting wit hthe value for t=1).

  SUBROUTINE tinterp( ntOut,nt,t,tindex,inval,disAggFlag,interpFlag,outval )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    nt      &!  the number of timesteps in period (interval) of input data
   ,ntOut   &!  the number of output values, i.e. the number of timesteps for which to calculate values
!                 This is currently restricted to be either 1 or nt, i.e. either calculate a single
!                 value (for a given time) or calculate values for all timesteps within data validity.
   ,t        !  the number of the first timestep (in the interval of length nt timesteps) to calculate
!                 a value for.
!                 e.g. 2 means calculate values starting at the 2nd timestep
!                 e.g. nt=10, ntOut=1,t=2
!                        nt=10   : period of input data is 10 timesteps
!                        ntOut=1 : values are required for a single timestep
!                        t=2     : values are required starting 2nd timestep in interval of 10 (and
!                                   since ntOut=1, this is the only timestep that values are required for)
!   ,tUse     !  the first position in outval to store data in
!                 Often tUse=t, but if we want to store the output starting in the first element
!                 of outval, we use tUse=1 (e.g. want to interpolate with time but return answer to first element).

   INTEGER ::  &!  local SCALARS
     jt        &!  loop counter
    ,jt1,jt2    !  bounds for loop (timesteps within the 1:nt timesteps of the data period)

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    tindex(2)   !  index of "relative" times of data

  REAL ::   &!  local SCALARS
    denom   &!  denominator of scaling factor for conserving interpolation
   ,numer    !  numerator of scaling factor for conserving interpolation

  REAL, INTENT(in) ::  &!  in ARRAYS
    inval( tindex(1):tindex(2) )   !  input data

  REAL, INTENT(out) ::  &!  out ARRAYS
    outval( ntOut )    !  data after interpolation/disaggregation

  REAL ::  &!  local SCALARS
    rnt     !  real version of nt

  REAL ::  &!  local ARRAYS
    fac(-1:1)  !  weights for previous, current and next values

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    disAggFlag    &!   type of temporal disaggregation to be performed
   ,interpFlag     !   type of temporal interpolation to be performed

!--------------------------------------------------------------------------------
!  print*,'tinterp: ntOut,nt,t,tindex=',ntOut,nt,t,tindex
!  print*,'tinterp: inval,disAggFlag,interpFlag=',inval,disAggFlag,interpFlag

! Checks.
! This is bit of a faff - alternative is not to check for errors and accept the consequences of any duff input.

  IF ( ntOut > nt ) THEN
    WRITE(*,*)'ERROR: tinterp: ntOut>nt: ntOut=',ntOut,' nt=',nt
    WRITE(*,*)'Cannot calculate values for more times than are in one data interval.'
    STOP
  ENDIF
  IF ( t > nt ) THEN
    WRITE(*,*)'ERROR: tinterp: t>nt: timestep number=',t,' nt=',nt
    WRITE(*,*)'Cannot calculate a value for a timestep after end of data validity.'
    STOP
  ENDIF

  IF ( ntOut>1 .AND. ntOut==nt ) THEN
!   Want to output data for more than one time, and this number is all times for which input is valid.
    IF ( t /= 1 ) THEN
      WRITE(*,*)'ERROR: tinterp: ntOut>1 .AND. ntOut==nt must use t=tUse=1'
      WRITE(*,*)'i.e. if requesting values for more than one time, must request all times'
      STOP
    ENDIF
!   Do all times.
    jt1 = 1
    jt2 = nt
  ELSEIF ( ntOut==1 .AND. nt==1 ) THEN
!   Input data is only valid for one timestep.
    IF (  t /= 1 ) THEN
      WRITE(*,*)'ERROR: tinterp: ntOut=1 .AND. nt=1'
      WRITE(*,*)'i.e. if there only a single value is available, it must be requested'
      STOP
    ENDIF
!   Do all times (equals one).
    jt1 = 1
    jt2 = jt1
  ELSEIF ( ntOut==1 .AND. nt>1 ) THEN
!   A single value is requested (from a longer 'possible' vector of interpolated  data).
    jt1 = t
    jt2 = jt1
  ENDIF

!-------------------------------------------------------------------------------
  rnt = REAL( nt )

!  print*,'TOP tinterp tindex=',tindex(:),' inval=',inval(:)
!  print*,'nt,ntOut,t=',nt,ntOut,t
!  print*,'jt1,jt2=',jt1,jt2

  SELECT CASE ( disAggFlag )
    CASE ( '0' )
!     No temporal disaggregation.
    CASE default
      WRITE(*,*)'ERROR: tinterp: do not recognise disaggregation flag=',TRIM(disAggFlag)
      STOP
  END SELECT

!--------------------------------------------------------------------------------


  SELECT CASE ( interpFlag )

    CASE ( 'i' )
!     Inputs are instantaneous values at given times.
      DO jt=jt1,jt2
        fac(0) = REAL(nt+1-jt)/REAL(nt)
!        fac(0) = real(nt-jt)/real(nt)
        fac(1) = 1.0 - fac(0)
        outval(jt-jt1+1) = inval(0)*fac(0) + inval(1)*fac(1)
!        print*, 'i jt=',jt,' fac=',fac(0:1),' inval=',inval(0:1)
!        print*, jt-jt1+1,' outval=',outval(jt-jt1+1)
      ENDDO

    CASE ( 'b' )
!     Inputs are backward time averages, i.e. time average ending at given time (in GSWP2 this is L)
!     We have previously tested that nt is even, so algorithm conserves.
      DO jt=jt1,jt2
        fac(0) = ( 2.0*rnt - ABS(REAL(2*jt-nt-1)) ) / (2.0*rnt)
        fac(-1) = MAX( 1.0-REAL(jt*2+nt-1)/(2.0*rnt), 0.0 )
        fac(1) = MAX( 1.0-REAL((nt+1-jt)*2+nt-1)/(2.0*rnt) , 0.0 )
        denom = 0.5*(inval(0)+inval(2)) + 3.0*inval(1)
        numer = 4.0*inval(1)
        IF ( denom > EPSILON(denom) ) THEN
          outval(jt-jt1+1) = (inval(0)*fac(-1)+inval(1)*fac(0)+inval(2)*fac(1)) * numer / denom
        ELSE
           outval(jt-jt1+1) = 0.0
        ENDIF
      ENDDO

    CASE ( 'c' )
!     Inputs are centred time averages, i.e. time average centred on given time (GSWP2 C)
!     We have previously tested that nt is even, so algorithm conserves.
      DO jt=jt1,jt2
        fac(0) = ( 2.0*rnt - (REAL(2*jt-1)) ) / (2.0*rnt)
        fac(-1) = 0.0
        fac(1) = 1.0 - fac(0)
        IF ( jt > nt/2 ) THEN
          denom = 0.5*(inval(0)+inval(2)) + 3.0*inval(1)
          numer = 4.0*inval(1)
        ELSE
          denom = 0.5*(inval(-1)+inval(1)) + 3.0*inval(0)
          numer = 4.0*inval(0)
        ENDIF
        IF ( denom > EPSILON(denom) ) THEN
           outval(jt-jt1+1) = ( inval(-1)*fac(-1) + inval(0)*fac(0) + inval(1)*fac(1) ) * numer / denom
        ELSE
           outval(jt-jt1+1) = 0.0
        ENDIF
      ENDDO

    CASE ( 'f' )
!     Inputs are forward time averages, i.e. time average starting at given time (GSWP2 N)
!     We have previously tested that nt is even, so algorithm conserves.
      DO jt=jt1,jt2
        fac(0) = ( 2.0*rnt - ABS(REAL(2*jt-nt-1)) ) / (2.0*rnt)
        fac(-1) = MAX( 1.0-REAL(jt*2+nt-1)/(2.0*rnt), 0.0 )
        fac(1) = MAX( 1.0-REAL((nt+1-jt)*2+nt-1)/(2.0*rnt) , 0.0)
        denom = 0.5*(inval(-1)+inval(1)) + 3.0*inval(0)
        numer = 4.0*inval(0)
        IF (denom > EPSILON(denom)) THEN
           outval(jt-jt1+1) = ( inval(-1)*fac(-1) + inval(0)*fac(0) + inval(1)*fac(1) ) * numer / denom
        ELSE
           outval(jt-jt1+1) = 0.0
        ENDIF
      ENDDO

    CASE ( 'nb' )
!     No interpolation, (current) value is valid over time interval ending at given time (not in GSWP2).
!     Use the 'next' time.
      outval(:) = inval(1)

    CASE ( 'nc' )
!     No interpolation, (current) value is valid over time interval centred on given time  (GSWP2 0)
      outval(:) = inval(0)     !   initialise as using "current" time
!     Decide what times will use the next data value.
!      jt = max( nt/2+1, ceiling(real(nt)/2.0) )   ! odd t=int(t/2)+1, even t=t/2+1
!     If period is an even number of timesteps, start using next value halfway through the interval.
!     If period is odd, start using next value only after halfway through the interval
!      (eg 3 timestep data, 2nd timestep starts only 1 hour after timestamp for current data,
!       so continue to use current data).
      jt = CEILING( REAL(nt)/2.0 ) + 1      !  timestep in 1:nt at which we switch to later value
      IF ( jt <= nt ) THEN
        IF ( ntOut > 1 ) THEN
          outval(jt:) = inval(1)
        ELSE
          outval(1) = inval(1)
        ENDIF
      ENDIF

    CASE ( 'nf' )
!     No interpolation, (current) value is valid over time interval starting at given time (GSWP2 default)
      outval(:) = inval(0)

    CASE default
      WRITE(*,*)'ERROR: tinterp: do not recognise interpolation flag=',TRIM(interpFlag)
      STOP

  END SELECT  !   interpFlag

!  print*,'bottom tinterp outval(:)=',outval(:)

  END SUBROUTINE tinterp
!###############################################################################
!###############################################################################
!################################################################################

! subroutine data_init
! Does various calculations involved in setting up the reading of time-varying
! data (e.g. driving data), such as working out where the first data to be read
! are held in the input file. Also calls routines to read any times of
! data (e.g. earlier times) that are required at the first timestep.
! Note that this routine is not called for data that do not vary with time
! (e.g. "static" vegetation fields) - although it might be better/more
! consistent if it was.

  SUBROUTINE data_init( dataPer,dataFilePer,templateDate,templateTime  &
    ,dataStepMax,dataStep,dataStepInit  &
    ,dataDate,dataDateInit  &
    ,dataFile,dataFileStep,resetStep,resetStepPrev  &
    ,dataTime,dataTimeInit  &
    ,dataFileDate,dataFileTime,fileTemplate,fileTemplateUnits,dataFilename  &
    ,dataTimeIndex,climatol,useEndTime,templateT,notNext,callType  &
    ,dataUpdatePer,dataUpdateStepMax,dataUpdateStep,dataUpdateStepInit )

  USE inout, ONLY :  &!
!   imported scalar parameters
     periodMon

  USE misc_utils, ONLY :  &!
!   imported procedures
     find_file

  USE spin_mod, ONLY : &
!  imported scalars with intent(in)
    nspin,spinUp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
    l_360

  USE time_loc, ONLY :   &
!  imported scalars with intent(in)
    timeStep  &
!  imported arrays with intent(in)
   ,dateRun,dateMainRun,dateSpin,timeRun

  USE time_mod, ONLY : &
!  imported procedures
    timeDate,timeDate_diff

  USE timeConst, ONLY : &
    iSecInDay,iSecInHour,iSecInMin

!  use veg_io_vars, only : &
!  imported scalars with intent(in)
!    vegDataDay

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL, PARAMETER :: initCall = .TRUE.   !   argument to data_init_vals,
!        TRUE when call is during initialisation (as now)

  INTEGER, INTENT(in) ::  &!  in SCALARS
    dataPer       &!  period of input data
!                       >0 is number of timesteps
!                      <=0 is a special case (eg monthly)
   ,dataFilePer   &!  period of data files
!                       >0 is number of timesteps
!                      <=0 is a special case (eg monthly)
   ,templateDate  &!  date associated with template files
   ,templateTime   !  time associated with template files

  INTEGER, INTENT(inout) ::  &!  inout scalars
    dataStepMax          !  Maximum value of dataStep (often equals dataPer).
!     In previous code, this was an optional argument (as driving data currently
!     does not use - only really used for monthly veg data), but that got a bit
!     untidy in here, so now it's not optional...and the downside is that although
!     this is an inout argument, driving code elsewhere assumes that the out value=in value!

  INTEGER, OPTIONAL, INTENT(in) ::  &!  optional in scalars
    dataUpdatePer

  INTEGER, OPTIONAL, INTENT(inout) ::  &!  optional inout scalars
    dataUpdateStep      &!
   ,dataUpdateStepInit  &!
   ,dataUpdateStepMax

  INTEGER, INTENT(out) ::  &!  SCALARS
!!  NOTE: The following XInit variables refer to the latest data that are read at the
!!        start of the run. They could all be recalculated rather than stored....but right now it's
!!        much easier if I just store them.
!   In fact...I think these variables refer to the first time read, not the last time!
    dataStepInit    &!  the value of dataStep used at the start of the run
   ,dataDateInit    &!  the value of dataDate used at the start of the run
   ,dataTimeInit     !  the value of dataTime used at the start of the run

  INTEGER, INTENT(out) ::  &!  out SCALARS
    dataStep       &!  counter of timesteps within a period of input data
!                         i.e. index of currently-used data time, in a vector of dataPer times
   ,dataDate       &!  date (yyyymmdd) of the last data that were read (furthest forward in time)
   ,dataTime       &!  time of day (s) of the last data that were read (furthest forward in time)
   ,dataFile       &!  the number (index) of the data file that is currently in use
   ,dataFileStep   &!  index of the time level last read from a data gfile (i.e. 1,2,3,..)
!                         This index refers to dataTime, dataDate.
   ,resetStep      &!  the timestep number (a_step) when time and date of data will be reset,
                         !    to account for spin up. This is not used if nspin<0 (model-determined spin up).
   ,resetStepPrev   !  previous value of resetStep

  INTEGER, INTENT(in) ::  &!  in arrays
    dataTimeIndex(2)       !  index showing range of times of data held
!          index(1) and index(2) are the earliest and latest time respectively
!          Values are relative to the current time when the data are read, and are
!          given in terms of the period of the input data:
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  INTEGER, INTENT(inout) ::  &!  inout arrays
    dataFileDate(:)       &!  the date of first data in each data file
   ,dataFileTime(:)        !  the time of first data in each data file

  INTEGER ::  &!  local SCALARS
    dt                           &!  work
   ,ndy,nhr,nmin,nsec            &!  work
   ,readDate,readTime            &!  date and time at which data will be read
!                                     There might not be data with exactly this timestamp.
   ,tmpDate,tmpDate2,tmpTime      !  work

  LOGICAL, INTENT(in) ::  &!  in scalars
    climatol   &!  switch indicating if files are climatological
   ,useEndTime &!  switch indicating file-naming convention if templateT=TRUE
   ,templateT   !  switch indicating if time templating is used for file names

  LOGICAL, INTENT(out) ::  &!  out SCALARS
    notNext   !    switch indicating if next data to be used are next in file

  LOGICAL ::  &!  local scalars
    stepMax   !  work

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType     &!  flag indicating what data this call refers to
!                  'drive' = driving (meteorological) data
!                  'veg' = time-varying vegetation data
   ,fileTemplate &!  template file name
   ,fileTemplateUnits  !  time units associated with template file name

  CHARACTER(len=*), INTENT(inout) ::  &!  inout arrays
    dataFilename(:)   !  names of data files

!--------------------------------------------------------------------------------

! At present, optional arguments should only be present for callType='veg', and then
! all arguments should be present - so no further testing is done.

  SELECT CASE ( callType )
    CASE ( 'drive','veg' )   !  acceptable values
    CASE default
      WRITE(*,*)'ERROR: data_init: do not recognise callType=',TRIM(callType)
      STOP
  END SELECT

! Initialise timestep numbers for "resetting" as something that cannot be reached (<0).
  resetStep = -1
  resetStepPrev = -1
! Initialise flag.
  stepMax = .FALSE.

! Check dataPer. Of the "special" values, only monthly has been coded for.
! Checking now saves having to check several times below.
  IF ( dataPer<0 .AND. dataPer/=periodMon ) THEN
    WRITE(*,*)'ERROR: data_init: dataPer<0 .AND. dataPer/=periodMon (=',periodMon,')'
    WRITE(*,*)'There is no code for this yet - sorry!'
    STOP
  ENDIF

!-------------------------------------------------------------------------------
! Get the time when data are first required.
!-------------------------------------------------------------------------------
! This time is not necessarily the time of the first data that will be required,
! since earlier data may be needed for interpolation (or later data if they are
! backward averages).
! Start by providing time at the start of the first timestep.
  readTime = timeRun(1)
  readDate = dateRun(1)
  IF ( spinUp ) readDate = dateSpin(1)

!-------------------------------------------------------------------------------
! Locate a time with data, at or before the time when data are first required.
! If before, this is the time closest to the time when data required.
!-------------------------------------------------------------------------------
  CALL data_time(  dataPer,dataFileDate(1),dataFileTime(1),dataTimeIndex(1)  &
                       ,'data_init',climatol  &
                       ,readDate,readTime  &
                       ,dataDate,dataTime )

!-------------------------------------------------------------------------------
! Find the file that holds the required data.
!-------------------------------------------------------------------------------
! First find the file.
  CALL find_file( dataTime,dataDate,climatol,useEndTime,.FALSE.,dataPer,NINT(timeStep)  &
                 ,dataFile,dataFileStep  &
                 ,templateT,dataFilePer,dataFileTime,dataFileDate  &
                 ,dataFileName,fileTemplate,templateDate,templateTime  &
                 ,fileTemplateUnits,'data_init' )
!  PRINT*,'+find_file dataFile=',dataFile,' dataFileStep=',dataFileStep
!  PRINT*,'Time/date of first data to be read=',dataTime,dataDate,' dataFile=',dataFile

dataFileStep=1   !DSM - Forcando este valor para nao identificar problema de data
  IF ( dataFileStep < 1 ) THEN
    WRITE(*,*)'ERROR: data_init: First required data are before start of available data.'
    WRITE(*,*)'Required time/date=',dataTime,dataDate,' first available=',dataFileTime(1),dataFileDate(1)
    STOP
  ENDIF

! Note: if the interpolation code is 'c' (no interpolation, value valid at times centred on given time),
! we may have demanded an earlier datum when it was not required. If time of first data read is more than
! halfway through the data period, the 'current' value (index=0) is not actually needed.
! This is the least of my worries, but I thought I'd point it out.

!  PRINT*,'later.... dataFile=',dataFile,' dataFileStep=',dataFileStep
!-------------------------------------------------------------------------------
! Save the date and time of the earliest data needed.
  dataTimeInit = dataTime
  dataDateInit = dataDate

!--------------------------------------------------

! Initialise dataStepMax.
  IF ( dataPer > 0 ) THEN
!   At present we expect dataStepMax to be passed in as dataPer, and not to change,
!   but I've added code here in case we later do decide that dataStepMax may vary here!
    dataStepMax = dataPer
  ELSEIF ( dataPer == periodMon ) THEN
    IF ( callType/='veg' ) THEN
      WRITE(*,*)'ERROR: data_init: monthly data must be veg!'
      WRITE(*,*)'Stopping in data_init'
      STOP
    ENDIF
!   This value is over-written by the first call to veg_update, so it only needs
!   to satisfy/trigger certain bits of code during initialisation! Helpful, eh?
!   The following is from veg_update, so if used here should implement as function.

!   At this point, dataDate gives the date of the data that will be read next,
!   but we may need to use an earlier month - depends on interpolation.
!   dataTimeIndex can be used to calculate how many times of data are used.
!   We then take dataDate and go back by the number of times (months) that are used.
!   Some of this code is close to or actually circular...but it's the best I can do for now.

!   Earlier we have restricted monthly veg to have interpFlags i,nb,nc or nf,
!   and this code is written assuming that holds - check via index.
    IF ( ANY(dataTimeIndex(:)<0) .OR. ANY(dataTimeIndex(:)>1) ) THEN
      WRITE(*,*)'ERROR: data_init: monthly data with unexpected datTimeIndex,'
      WRITE(*,*)'suggesting unexpected interpolation flag.'
      WRITE(*,*)'New code required.'
      STOP
    ENDIF
!   Now we know that we only need times 0 and 1, which turns out to mean that we only need
!   to know the length of the current month (actually from 1 month from data time). Start
!   from time of the last data read, and work back.
!   Work out required offset from dataTimeIndex(2).
    dt = 0
    IF ( dataTimeIndex(2) == 1 ) dt = 1
!   Now we need to work out the time/date of the latest veg, knowing earliest veg.
!   We do this (rather than starting from earliest veg) for consistency with veg_update.
!   Again, coded assuming we have restricted dataTimeIndex to 0 or 1.
    tmpDate = dataDateInit
    IF ( dataTimeIndex(2)==1 .AND. dataTimeIndex(1)==0 ) CALL timeDate(  &
                                   dataTimeInit,dataDateInit  &
                                  ,1,'mon',l_360,tmpTime,tmpDate,'data_init' )
!   tmpDate is now time of latest veg.
!   Take latest data time, and wind back by dt months.
    CALL timeDate( 0,tmpDate,-dt,'mon',l_360,tmpTime,tmpDate2,'data_init')
!   Now get date one month later than this (i.e. end of data validity).
!   I'm assuming (as elsewhere) that the same day of month can be used in both months.
    CALL timeDate( 0,tmpDate2,1,'mon',l_360,tmpTime,tmpDate,'data_init')
!   Get difference between times (the length of first month).
    CALL timeDate_diff( 0,tmpDate2,0,tmpDate,l_360,'data_init',nsec,nmin,nhr,ndy )
    dataStepMax = ( ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec ) / NINT(timeStep)
  ENDIF  !  dataPer
!  print*,'dataStepMax=',dataStepMax

!----------------------------------------------------------

! Calculate how far through the data interval the start time lies (i.e. calculate dataStep).

! We can use the time of the first data that are used (not necessarily the data closest
! in time to the start of the run), since the relationship of either to the start time
! is the same once modulo is taken (below).

! Get time difference between time of first data and start of run.
  CALL timeDate_diff( dataTimeInit,dataDateInit,timeRun(1),dateRun(1),l_360,'data_init',nsec,nmin,nhr,ndy)
  tmpTime = ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec

! Convert to number of timesteps.
  dataStepInit = tmpTime/ NINT(timeStep)
  dataStepInit = MODULO( dataStepInit, dataStepMax )
  dataStepInit = dataStepInit + 1   !  since tmpTime=0 corresponds to dataStepInit=1
! dataStepInit = 1 means that data will be read at start of first timestep. Indicate this as dataStepMax+1.
  IF ( dataStepInit == 1 ) THEN
    dataStepInit = dataStepMax + 1
!   Set flag to indicate that we want to have dataStepInit = max+1. This is
!   necessary because dataStepMax is currently an arbitrary value, and may well be changed
!   when veg_update reads data. If flag is set, we can later update dataStepInit.
!   Possibly now redundant...but still doing it.
    stepMax = .TRUE.
  ENDIF

! Initialise counter using this value.
  dataStep = dataStepInit
!-------------------------------------------------------------------------------
! Initialise counters for update of data (between data reads).
! If necessary, calculate how far between updates the start time lies.
! This is only necessary for fields that need not be updated every timestep.

  IF ( PRESENT( dataUpdatePer) ) THEN

!   Initialise dataUpdateStepMax.
    IF ( dataUpdatePer > 0 ) THEN
      dataUpdateStepMax = dataUpdatePer
    ELSE
!     Coded assuming monthly data and monthly updates.
      dataUpdateStepMax = dataStepMax
    ENDIF

!   Initialise dataUpdateStep.
    IF ( dataUpdatePer == 1 ) THEN
!     dataUpdateStepInit = 1 means that data will be updated at start of first timestep.
!     Indicate this as dataUpdateStepMax+1.
      dataUpdateStepInit = dataUpdateStepMax + 1
    ELSEIF ( dataUpdatePer > 1 ) THEN
!     Get time difference between earliest data and start of run.
!     We can use the time of the first data that are used (not necessarily the data closest
!     in time to the start of the run), since the relationship of either to the start time
!     is the same once modulo is taken (below).
      CALL timeDate_diff( dataTime,dataDate,timeRun(1),dateRun(1),l_360,'data_init',nsec,nmin,nhr,ndy)
      tmpTime = ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec

!     Convert to number of timesteps.
      dataUpdateStepInit = tmpTime/ NINT(timeStep)

      dataUpdateStepInit = MODULO( dataUpdateStepInit, dataUpdatePer )
      dataUpdateStepInit = dataUpdateStepInit + 1   !  since tmpTime=0 corresponds to dataUpdateStepInit=1
!     dataUpdateStepInit = 1 means that data will be updated at start of first
!     timestep. Indicate this as dataUpdateStepMax+1.
      IF ( dataUpdateStepInit == 1 ) dataUpdateStepInit = dataUpdateStepMax+1
    ELSE
!     Assuming monthly data and monthly updates.
      dataUpdateStepInit = dataStepInit
    ENDIF  !  dataUpdatePer

  ENDIF   !  present(dataUpdatePer)

!-------------------------------------------------------------------------------

  IF ( spinUp ) THEN

!   Calculate when we will need to "reset" the data for the start of another cycle of spin up.
!   This is not necessary for the case when the dates of the main run follow on from the
!   spin up - in that case the resetting is done at the start of a spin-up cycle (see newTime).
!   In fact, here we calculate a time BEFORE the start of the run - when we would have had to
!   reset had this not been the start of the run - and this is later incremented to calculate
!   when to reset in next cycle.
    IF ( nspin>0 .AND. dateSpin(1)==dateMainRun(1) ) THEN
!     Work out what timestep the initial data correspond to.
      CALL timeDate_diff( dataTimeInit,dataDateInit,timeRun(1),dateRun(1),l_360,'data_init'  &
                     ,nsec,nmin,nhr,ndy )
      nsec = ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec
!      print*,'nsec=',nsec,' from ',dataTimeInit,dataDateInit,' to ',timeRun(1),dateRun(1)
!XX   Debug test!
      IF ( MOD(nsec,NINT(timeStep)) /= 0 ) THEN
        WRITE(*,*)'ERROR: data_init: earliest data needed have not given an integer # of timesteps.'
        WRITE(*,*)'dataTimeInit=',dataTimeInit,' dataDateInit=',dataDateInit
        WRITE(*,*)'nsec=',nsec,' nsec-mod(nsec,nint(timeStep))*nint(timeStep)=',nsec-MOD(nsec,NINT(timeStep))*NINT(timeStep)
        WRITE(*,*)'Stopping in data_init'
        STOP
      ENDIF
!     The 1 on the next line can be understood in terms of the fact that nsec=0 corresponds to timestep #1.
      resetStep = 1 - nsec/NINT(timeStep)
!     Account for the fact that data are often needed, before their nominal time, for interpolation.
      IF ( dataPer > 0 ) THEN
        resetStep = resetStep - dataTimeIndex(2)*dataPer
      ELSEIF ( dataPer == periodMon ) THEN
!       We have to go back a further dataTimeIndex months.
        CALL timeDate( dataTimeInit,dataDateInit,-dataTimeIndex(2),'mon',l_360,tmpTime,tmpDate,'data_init: resetStep')
!       Calculate length of time back to this new date.
        CALL timeDate_diff( tmpTime,tmpDate,dataTimeInit,dataDateInit,l_360,'data_init: resetStep'  &
                     ,nsec,nmin,nhr,ndy )
        nsec = ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec
        resetStep = resetStep - nsec / timeStep
      ENDIF

    ENDIF   !  spin-up dates

  ENDIF  !  spinUp

!-------------------------------------------------------------------------------

! Read the data and calculate initial values.
! For now, use an IF clause to deal with optional arguments.
  IF ( PRESENT(dataUpdateStep) ) THEN
    CALL data_init_vals( dataDateInit,dataPer,dataStepMax,dataStepInit  &
                        ,dataTimeInit,dataStep,dataDate,dataTime  &
                        ,dataTimeIndex,initCall,notNext,callType  &
                        ,dataUpdateStep,dataUpdateStepInit )
  ELSE
    CALL data_init_vals( dataDateInit,dataPer,dataStepMax,dataStepInit  &
                        ,dataTimeInit,dataStep,dataDate,dataTime  &
                        ,dataTimeIndex,initCall,notNext,callType )
  ENDIF

! If flag set, update dataStep.
  IF ( stepMax ) THEN
    dataStepInit = dataStepMax + 1
    dataStep = dataStepInit
  ENDIF

  END SUBROUTINE data_init
!################################################################################
!################################################################################
! subroutine data_time
! Given a date and time, calculates the time to use for data.
! This is the closest time to that requested, with the condition that the
! returned time is not later than that requested.
! Expected use is as part of initialisation (via data_init).

  SUBROUTINE data_time( dataPer,dateFile,timeFile,timeOffset  &
                       ,callType,climatol  &
                       ,dateReq,timeReq  &
                       ,dateData,timeData )

  USE inout, ONLY :  &
!   imported scalar parameters
     periodMon

  USE switches, ONLY :   &
!   imported scalars with intent(in)
     l_360

  USE time_loc, ONLY :   &
!   imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :   &
!   imported procedures
     dateToBits,timeDate,timeDate_cmp

  USE timeConst, ONLY : &
!   imported scalar parameters
     iSecInDay

!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)

  INTEGER, INTENT(in) ::  &
    dataPer    &!    period of input data
!                     >0 is number of timesteps
!                    <=0 is a special case (eg monthly)
   ,dateFile   &!  date of start of first data file
   ,timeFile   &!  time of day (s) of start of first data file
   ,timeOffset  !  offset indicating how many earlier data are required for
!                    interpolation

  LOGICAL, INTENT(in) ::  &
    climatol  !  flag indicating is data are to be treated as "climatological"
!                  T=climatological

  CHARACTER(len=*), INTENT(in) ::  &
    callType      !  flag indicating what data this call refers to
!                  'drive' = driving (meteorological) data
!                  'veg' = time-varying vegetation data
!-------------------------------------------------------------------------------
! Scalar arguments with intent(inout)

  INTEGER, INTENT(inout) ::  &
    dateReq    &!  date for which data requested
   ,timeReq     !  time of day (s) for which data requested

!-------------------------------------------------------------------------------
! Scalar arguments with intent(out)

  INTEGER, INTENT(out) ::  &
    dateData   &!  date of data
   ,timeData    !  time of day (s) of data
!-------------------------------------------------------------------------------
! Local scalar variables.

  INTEGER ::  &
    i          &!  counter
   ,datePrev   &!  work: a date
   ,dayFile    &!  work: a day of month
   ,dayReq     &!  work: a day of month
   ,monthReq   &!  work: a month of year
   ,yearReq    &!  work: a year
   ,timePrev   &!  work: a time of day
   ,tmpDate    &!  work: a date
   ,tmpMonth   &!  work: a month of year
   ,tmpTime    &!  work: a time of day
   ,tmpYear     !  work: a year

!-------------------------------------------------------------------------------
! Adjust time to account for need for any earlier data for interpolation.
!-------------------------------------------------------------------------------
  IF ( dataPer > 0 ) THEN
    CALL timeDate( timeReq,dateReq,timeOffset*dataPer*NINT(timeStep),'sec'  &
                  ,l_360,tmpTime,tmpDate,'data_time')
  ELSEIF ( dataPer == periodMon ) THEN
    CALL timeDate( timeReq,dateReq,timeOffset,'mon'  &
                  ,l_360,tmpTime,tmpDate,'data_time')
  ELSE
    WRITE(*,*)'ERROR: date_time: no code for dataPer=',dataPer
    WRITE(*,*)'callType=',TRIM(callType)
    STOP
  ENDIF

  timeReq = tmpTime
  dateReq = tmpDate

!-------------------------------------------------------------------------------
! Find data time.
!-------------------------------------------------------------------------------

  IF ( dataPer > 0 ) THEN
!-------------------------------------------------------------------------------
!   Regular data.
!-------------------------------------------------------------------------------

!   Life isn't simple. Actually it is for periods <=1day, as these have been
!   "synchronised" with days. Longer periods are more difficult as finding a
!   reference point is harder - e.g. 10 day data will generally not have a
!   fixed relationship to start of month or year.
    IF ( dataPer*NINT(timeStep) <= isecInDay ) THEN
!     Period <= 1 day.
!     00H at the start of this day was a data time.
      timeData = INT( timeReq/(dataPer*NINT(timeStep)) ) * dataPer * NINT(timeStep)
      dateData = dateReq
    ELSE
!     Period > 1 day.
!     Tricky, without further assumptions. Use brute force and the only
!     fixed reference we have - namely the date and time given for the first
!     data. Advance this time until we get to the first time that is >
!     requested time, then use the previous time
!     Start point for search depends upon whether data are climatological or not.
      IF ( climatol )  THEN
!       We have previously ensured that we have l_360=TRUE and dataPer is a
!       factor of 1 year (so we can use same dates in each year).
!       Start at date of first file, with year changed to requested year-1 (minus -1
!       to avoid having to check if requested time is earlier or later in year).
        CALL dateToBits( dateReq,dayReq,monthReq,yearReq,l_360,'data_time' )
        CALL dateToBits( dateFile,dayFile,tmpmonth,tmpYear,l_360,'data_time' )
        dateData = (yearReq-1)*10000 + tmpMonth*100 + dayFile
        timeData = timeFile
      ELSE
!       Note that in this case we also check if time is before start of available
!       data - a check that is not done by this code for other data periods.
        timeData = timeFile
        dateData = dateFile
!       First check that required time is not before start of available data.
        IF ( timeDate_cmp( timeReq,dateReq,'<',timeFile,dateFile,'data_time' ) ) THEN
          WRITE(*,*)'ERROR: data_time: requested time is before start of data.'
          WRITE(*,*)'callType=',TRIM(callType)
          WRITE(*,*)'Requested time and date=',timeReq,dateReq
          STOP
        ENDIF
       ENDIF  !  climatol

       timePrev = -1
       datePrev = -1   !  so don't try to find dates <0!
       i = 0
       DO
        i = i + 1
!       An arbitrary guard against an infinite loop.
        IF ( i > 100000 ) EXIT
!       Check if time is > requested time.
        IF ( timeDate_cmp ( timeData,dateData,'>',timeReq,dateReq,'data_time' ) ) THEN
!         Use the previous time, as this was <= requested time.
          timeData = timePrev
          dateData = datePrev
          EXIT
        ENDIF
!       Advance time by data period.
        timePrev = timeData
        datePrev = dateData
        CALL timeDate( timePrev,datePrev,dataPer*NINT(timeStep),'sec',l_360  &
                      ,timeData,dateData,'data_time' )
      ENDDO
!     Check we found data.
      IF ( timePrev < 0 ) THEN
        WRITE(*,*)'ERROR: data_time: could not find suitable data.'
        WRITE(*,*)'callType=',TRIM(callType)
        WRITE(*,*)'Requested time and date=',timeReq,dateReq
        WRITE(*,*)'DataPer (s) = ',dataPer*NINT(timeStep)
        WRITE(*,*)'Time and date of first file=',timeFile,dateFile
        WRITE(*,*)'ERROR: data_time: max loop size exceeded'
        WRITE(*,*)'Check for input errors, or increase loop size!'
        STOP
      ENDIF
    ENDIF   !  datePer relative to 1 day

  ELSEIF ( dataPer == periodMon ) THEN

!-------------------------------------------------------------------------------
!   Monthly data.
!-------------------------------------------------------------------------------
!   Data are timestamped with a particular day of month.
!   Start with requested time and move back through months (to same day of month)
!   until we find a time <= requested time.
!   To get time (of day) of data, use that given by first file.
    timeData = timeFile
!   Work out starting date - must use given day of month.
!   Note that this day is known to exist in all months (i.e. <=28).
    CALL dateToBits( dateReq,dayReq,monthReq,yearReq,l_360,'data_time' )
    CALL dateToBits( dateFile,dayFile,tmpMonth,tmpYear,l_360,'data_time' )
!   First assume start in month given by dateReq.
    dateData = yearReq*10000 + monthReq*100 + dayFile
!   If required day/time are earlier in month, start one month earlier.
    IF ( timeDate_cmp( timeReq,dayReq,'<',timeFile,dayFile,'data_time' ) ) THEN
      CALL timeDate( timeData,dateData,-1,'mon',l_360,tmpTime,tmpDate,'data_time' )
      dateData = tmpDate
    ENDIF
!   Initialise for loop.
    datePrev = -1
    i = 0
    DO
      i = i + 1
!     An arbitrary guard against an infinite loop.
      IF ( i > 100000 ) EXIT
!     Check if time is <= requested time.
      IF ( timeDate_cmp ( timeData,dateData,'<=',timeReq,dateReq,'data_time' ) ) THEN
!       Set datePrev to >0 to show success.
        datePrev = 1
        EXIT
      ENDIF
!     Wind time back by 1 month.
      timePrev = timeData
      datePrev = dateData
      CALL timeDate( timePrev,datePrev,-1,'mon',l_360,timeData,dateData,'data_time' )
    ENDDO
!   Check we found data.
    IF ( datePrev < 0 ) THEN
      WRITE(*,*)'ERROR: data_time: could not find suitable data.'
      WRITE(*,*)'callType=',TRIM(callType)
      WRITE(*,*)'Requested time and date=',timeReq,dateReq
      WRITE(*,*)'dataPer = 1 month'
      WRITE(*,*)'Time and date of first file=',timeFile,dateFile
      WRITE(*,*)'ERROR: data_time: max loop size exceeded'
      WRITE(*,*)'Check for input errors, or increase loop size!'
      STOP
    ENDIF

  ELSE
!-------------------------------------------------------------------------------
    WRITE(*,*)'ERROR: date_time: no code for dataPer=',dataPer
    WRITE(*,*)'callType=',TRIM(callType)
    STOP

  ENDIF

  END SUBROUTINE data_time

!################################################################################
!################################################################################
! subroutine data_init_vals
! Internal procedure in update_mod.
! Set data to values as used at initial time.
! This is used to set the initial data (for time-varying fields), and also to reset at
! the start of a spin-up cycle (if dates of main run follow those for spin up).

  SUBROUTINE data_init_vals( dataDateInit,dataPer,dataStepMax,dataStepInit  &
                            ,dataTimeInit,dataStep,dataDate,dataTime &
                            ,dataTimeIndex,initCall,notNext,callType  &
                            ,dataUpdateStep,dataUpdateStepInit )

  USE inout, ONLY :  &!
!   imported scalar parameters
     periodMon

  USE switches, ONLY :  &
!   imported scalars with intent(in)
      l_360

  USE time_loc, ONLY :  &
!   imported scalars with intent(in)
     timeStep

  USE time_mod, ONLY :  &
!   imported procedures
     timeDate

  USE veg_io_vars, ONLY :  &
!   imported scalars with intent(out
     vegDataStepMax
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  LOGICAL, PARAMETER :: readAllSwitch = .FALSE. !  switch affecting data during spin up
!   TRUE means that a call here after initialisation (initCall=F) will read all
!       required times of data. This has the effect of instantaneously resetting
!       variables at the start of a spin-up cycle.
!   FALSE means that a call here after initialisation (initCall=F) will only read
!       the last required time of data. This can have the effect of more gradually
!       resetting variables at the start of a spin-up cycle, because we don't
!       immediately change all values used in interpolation (but if there is no
!       temporal interpolation, this switch can't have any effect!).
!  Before this functionality was added, the code used to behave as if readAllSwitch=TRUE.

  INTEGER, INTENT(in) ::  &!  in SCALARS
    dataDateInit   &!  dataDate at start of run
   ,dataPer        &!  period of data
!                        >0 is period in terms of number of model timesteps
!                       <=0 indicate "special" cases (e.g. monthly)
   ,dataStepInit   &!  dataStep at start of run
   ,dataTimeInit    !  dataTime at start of run

  INTEGER, INTENT(inout) ::  &!  inout SCALARS
    dataStep          &!  counter of timesteps when a datum is used
   ,dataStepMax        !  maximum value of dataStep (often equals dataPer)

  INTEGER, OPTIONAL, INTENT(in) ::  &!  optional in scalars
    dataUpdateStepInit

  INTEGER, OPTIONAL, INTENT(inout) ::  &!  optional inout scalars
    dataUpdateStep

  INTEGER, INTENT(out) ::  &!  out SCALARS
    dataDate   &!  date of latest data read
   ,dataTime    !  time of latest data read

  INTEGER ::  &!  local SCALARS
    i,t,t2    &!  work
   ,dataDateNew,dataTimeNew   !  work

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    dataTimeIndex(:)       !  indices of required data times

  LOGICAL, INTENT(in) :: initCall    !  TRUE when this subroutine is called during initialisation

  LOGICAL, INTENT(out) ::  &!  out scalars
    notNext    !  used to force a search for a file with data for the appropriate time,
!                          rather than assuming that "next" file can be used.

  LOGICAL :: readAll  !  TRUE all times of data will be read
                      !  FALSE only latest time of data will be read

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType    !  flag indicating what data this call refers to
!                  'drive' = driving (meteorological) data
!                  'veg' = time-varying vegetation data

!-------------------------------------------------------------------------------

! At present, optional arguments should only be present for callType='veg', and
! then all arguments should be present - so no further testing is done.

! If more than one time level of data is required (for interpolation), read any
! earlier times. Also, if data are not read on the first timestep (because that
! is not time for reading), read all data required.
!-------------------------------------------------------------------------------

! Reset time of data to initial value. Decrement by one period, so that next
! increment returns required values.
  IF ( dataPer > 0 ) THEN
    CALL timeDate( dataTimeInit,dataDateInit,-1*dataPer*NINT(timeStep),'sec',l_360  &
                  ,dataTime,dataDate,'data_init_vals' )
  ELSEIF ( dataPer == periodMon ) THEN
    CALL timeDate( dataTimeInit,dataDateInit,-1,'mon',l_360  &
                  ,dataTime,dataDate,'data_init_vals' )
  ENDIF

! Work out the latest data to read.
  t2 = dataTimeIndex(2)
! If data are to be read on first timestep, no need to read the last time level here.
  IF ( dataStepInit > dataStepMax ) t2 = t2 - 1

! Set switch affecting how many data are read.
  readAll = .TRUE.
  IF ( .NOT.initCall .AND. .NOT.readAllSwitch ) readAll = .FALSE.

!-------------------------------------------------------------------------------
! If we only hold a single time of data (dataTimeIndex(1)=dataTimeIndex(2)), and
! this is to be read on the first timestep, we currently have t2=dataTimeIndex(1)-1,
! and we don't want to read any data here. So only read if t2>=dataTimeIndex(1).

  IF ( dataTimeIndex(1) <= t2 ) THEN

!   Loop over times of data.
    DO t=dataTimeIndex(1),t2

!     Decide whether to read this time or not.
!     If readAll=FALSE, only the latest time is read.
      IF ( readAll .OR. t==t2 ) THEN

!       Read data.

!       Reset data step so that the next increment (e.g. in drive_update) triggers
!       more data to be read.
        dataStep = dataStepMax

!       Call (e.g.) drive_update with first argument=0 (zeroth or 'special' timestep)
!       and optional argument specified to indicate which time level of data to
!       use at t=0. Argument 'next' should be false, to ensure files are searched
!       for correct time.
        i = dataStepInit
        IF ( i == dataStepMax+1 ) i=dataStepMax 
        SELECT CASE ( callType )
          CASE ( 'drive' )
            Print*, "Desabilitado esta chamada"
            STOP
            !DSM CALL drive_update( 0,.FALSE.,i )
          CASE ( 'veg' )
            CALL veg_update( 0,.FALSE.,i )
!           For "irregular" (eg monthly) data, we need to update dataStepMax
!           before veg_update is next called (in the next loop).
!           This assignment means that vegDataStepMax is updated here, and no longer needs
!           to be passed as an argument - but it's still passed until this code is rewritten
!           (or scrapped...).
            dataStepMax = vegDataStepMax
          CASE default
            WRITE(*,*)'ERROR: data_init_vals: do not recognise callType=',TRIM(callType)
            STOP
        END SELECT

     ELSE

!      Don't read data for this time.
!      Update the data time so that by end of loop we still have correct value.
!      This update wouldn't be necessary if data time was calculated differently
!      (i.e. accounting for readAll) above...but easier this way for now.
       
       IF ( dataPer > 0 ) THEN
         CALL timeDate( dataTime,dataDate,dataPer*NINT(timeStep),'sec',l_360  &
                       ,dataTimeNew,dataDateNew,'data_init_vals' )
       ELSEIF ( dataPer == periodMon ) THEN
         CALL timeDate( dataTime,dataDate,1,'mon',l_360  &
                       ,dataTimeNew,dataDateNew,'data_init_vals' )
       ENDIF
       dataTime = dataTimeNew
       dataDate = dataDateNew

     ENDIF  !  readAll

    ENDDO  !  t (times)
!-------------------------------------------------

  ELSE

!   No data were read here (because all data required will be read at start of
!   next timestep). Set switch so that when the next data are read, it is not
!   assumed that these are the "next" data in the file

    notNext = .TRUE.

  ENDIF
 
! Reset counters. Minus 1 so that first increment takes back to original value.
  dataStep = dataStepInit - 1
  IF ( PRESENT(dataUpdateStep) ) dataUpdateStep = dataUpdateStepInit - 1

  END SUBROUTINE data_init_vals
!###############################################################################
!###############################################################################
! subroutine calc_reset_step
! Internal procedure in update_mod.
! Calculates when data need to be changed because of spin up.

  SUBROUTINE calc_reset_step( a_step,resetStep,resetStepPrev,callType )

  USE switches, ONLY :  &
!   imported scalars with intent(in)
     l_360

  USE time_loc, ONLY :  &
!   imported scalars with intent(in)
     timeStep  &
!   imported arrays with intent(in)
    ,dateMainRun,dateSpin

  USE time_mod, ONLY :  &!
!   imported procedures
     timeDate_diff

  USE timeConst, ONLY : &
     iSecInDay,iSecInHour

  IMPLICIT NONE
!-------------------------------------------------------------------------------

  INTEGER, INTENT(in) :: &!  in SCALARS
    a_step   !  timestep number   XX not used now

  INTEGER, INTENT(inout) ::  &!  inout scalars
    resetStep      &!  the timestep number (a_step) when time and date of  data will be reset,
                         !    to account for spin up. This is not used if nspin<0 (model-determined spin up).
   ,resetStepPrev   !  previous value of resetStep

  INTEGER ::  &!  local SCALARS
    ndy,nhr,nmin,nsec  &!  work
   ,nstep               !  number of timesteps in a cycle of spin up

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType    !  flag indicating what data this call refers to
!                  'drive' = driving (meteorological) data
!                  'veg' = time-varying vegetation data
!-------------------------------------------------------------------------------

  SELECT CASE ( callType )
    CASE ( 'drive','veg' )   !  acceptable values
    CASE default
      WRITE(*,*)'ERROR: calc_reset_step: do not recognise callType=',TRIM(callType)
      STOP
  END SELECT

! Save the value of timestep when values were last reset.
  resetStepPrev = resetStep

! If needed, set flag to indicate that data will have to 'cycle' back
! to deal with the start of next spin up or main run.
! The flag is not used if nspin>0 and times for main run follow spin up - in that
! case the data are automatically reset at start of each cycle of spin up.

  resetStep = -1

! Set resetStep to 0 to indicate that we need to cycle back.
  IF ( dateSpin(1) == dateMainRun(1) ) THEN
!   Spin up is over same times as (or at least starts at same time as) main run.
!   We will need to "cycle" back to start of data.
    resetStep = 0
  ELSE
!   Main run dates follow on from those of spin up.
!   Assume that we will not be cycling  (i.e. assume this is end of spin up) -
!   this is corrected later if necessary.
!   Nothing to do!
  ENDIF

! Work out when we need to cycle - get timestep number.
  IF ( resetStep == 0 ) THEN
!   Calculate the number of timesteps in a cycle of spin up.
!   A long spin-up and truncation could result in a wrong answer here....
    CALL timeDate_diff( 0,dateSpin(1),0,dateSpin(2),l_360,'calc_reset_step'  &
                       ,nsec,nmin,nhr,ndy )
    nStep = ( ndy*iSecInDay + nhr*iSecInHour + nsec ) / NINT(timeStep)
!   Increment the timestep when data were last 'reset' by the length of spin up cycle.
    resetStep = resetStepPrev + nstep
    IF ( resetStep < 1 ) THEN
      WRITE(*,*)'ERROR: calc_reset_step: resetStep<1'
      WRITE(*,*)'This probably means that you are trying to spin up over a period that is'
      WRITE(*,*)'short relative to the period of some forcing data.'
      WRITE(*,*)'e.g. spin up over 1 day, vegetation data period=10 days'
      WRITE(*,*)'This problem is a shortcoming of the code!'
      WRITE(*,*)'Possible solutions include: longer spin up or improve this code!!'
      WRITE(*,*)'If this error is not correct for your application, comment out this stop!'
      WRITE(*,*)'callType=',callType,' a_step=',a_step
      WRITE(*,*)'resetStepPrev=',resetStepPrev,' resetStep=',resetStep
      STOP
    ENDIF

  ENDIF

  END SUBROUTINE calc_reset_step

!###############################################################################
!###############################################################################

! subroutine veg_update
! Internal procedure in update_mod.
! Update veg data, if required. This may involve reading data from file.
! Note that much of this follows drive_update, so bug fixes to one are likely also needed in the other.

  SUBROUTINE veg_update( tstep,next,istepArg )

  USE ancil_info, ONLY : &
!  imported scalars with intent(in)
    land_pts,ntiles,tile_index,tile_pts  &
!  imported arrays with intent(in)
   ,frac

  USE initial_mod, ONLY :  &
!  imported procedures
    dump_io

  USE inout, ONLY : &
!  imported scalar parameters
    periodMon  &
!  imported scalars with intent(in)
    ,dumpFreq

  USE misc_utils, ONLY :  &
!  imported procedures
    find_file

  USE p_s_parms, ONLY :  &
!  imported arrays with intent(in)
    satcon  &
!  imported arrays with intent(inout)
   ,catch_snow,catch,infil_tile,z0_tile

  USE prognostics, ONLY :  &
!   imported arrays with intent(inout)
     canht_ft,lai

  USE switches, ONLY :  &
!  imported scalars with intent(in)
    can_model,l_360,l_aggregate

  USE time_loc, ONLY :  &
!  imported scalars with intent(in)
    timeStep

  USE time_mod, ONLY :  &
!  imported procedures
    timeDate,timeDate_diff

  USE timeConst, ONLY : &
    iSecInDay,iSecInHour,iSecInMin

  USE veg_io_vars, ONLY :  &
!   imported scalars with intent(in)
     nvegVarMax,vegClim,vegDataPer,vegDataStepMax,vegEndTime,vegFilePer,vegFileTemplate  &
    ,vegResetStep,vegResetStepPrev  &
    ,vegTemplateDate,vegTemplateTime,vegTemplateUnits,vegUpdatePer  &
    ,vegVarFlagPT,vegVarFlagPTX,vegVarFlagPX  &
    ,varNumCanht,varNumLAI,varNumRootd,vegVaryT  &
!   imported scalars with intent(inout)
    ,notNextVeg,vegDataStep,vegDate,vegDateInit,vegFile,vegFileStep  &
    ,vegTemplateT,vegTime,vegTimeInit,vegUpdateStep,vegUpdateStepMax  &
!   imported arrays with intent(in)
    ,vegTimeIndex,vegVarFlag  &
!   imported arrays with intent(inout)
     ,vegFileDate,vegFileName,vegFileTime

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    tstep    !   timestep number. 0 indicates a call during initialisation.

  INTEGER, INTENT(in), OPTIONAL ::  &!  optional in SCALARS
    istepArg    !  the current value of the time index
!                 Generally istep=vegDataStep, in which case istepArg is not provided.
!                 However, at initialisation, it may be necessary to use a different value of istep.

  INTEGER ::  &!  local SCALARS
    dt                 &!  work: index value
   ,istep,ivar         &!  loop counters
   ,ndy,nhr,nmin,nsec  &!  work
   ,tmpDate,tmpDate2,tmpTime,varNum    !  work

  LOGICAL, INTENT(in) ::  &!  inout SCALARS
    next  !  T means assume that the requested data are held in the next time level, either
!              of the current file or the next file if appropriate
!            F means search file details to locate the data

  LOGICAL ::  &!  local SCALARS
    nextData  &!  a local version of next. Value equals that of next, unless
!                      we need to "reset" the time of data.
   ,resetStep  !  TRUE when data time is to be reset, to deal with spin up
!--------------------------------------------------------------------------------

! Check that the data period is one that is coded for.
  IF ( vegDataPer<0 .AND. vegDataPer/=periodMon ) THEN
    WRITE(*,*)'ERROR: veg_update: vegDataPer<0 .AND. vegDataPer/=periodMon (=',periodMon,')'
    WRITE(*,*)'There''s no code for this case!'
    STOP
  ENDIF

! Initialise.
  resetStep = .FALSE.

! Increment counter of steps within data interval (i.e. between reads).
  vegDataStep = vegDataStep + 1

! Deal with resetting time and date for spin up.
! vegResetStep is used to reset veg data for spin up.
! Note that we test on both vegResetStep and vegResetStepPrev, because a test on only
! the former fails if we need to reset data on the first timestep of a new cycle of spin up.
! In that case, newTime has already recalculated vegResetStep before this routine has
! had a chance to reset the data - but we can catch that case using the previous value (vegResetStepPrev).
!XX    Is there a problem here if we have tstep=0 (to indicate "special" call, but this also happens to be vegResetPer?)
  IF ( tstep==vegResetStep .OR. tstep==vegResetStepPrev ) resetStep = .TRUE.

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Decide if it is time to read new veg data.
  IF ( vegDataStep>vegDataStepMax .OR. resetStep ) THEN

!   Read new data.

!   Reset counter.
    vegDataStep = 1

!   Trigger an update of the data.
    vegUpdateStep = vegUpdateStepMax + 1
!XX If resetStep, need to update...but maybe with differnt counter. eg istepArg (but not given)

!   Initialise.
    nextData = next

!-------------------------------------------------------------------------------
    IF ( vegVaryT ) THEN

!     If switch set, force code to search for appropriate time.
      IF ( notNextVeg ) THEN
        nextData = .FALSE.
        notNextVeg = .FALSE.  !  reset the switch
      ENDIF

!     Get date and time of the next veg data.
      IF ( vegDataPer > 0 ) THEN
        CALL timeDate( vegTime,vegDate,vegDataPer*NINT(timeStep),'sec',l_360,tmpTime,tmpDate,'veg_update' )
      ELSEIF ( vegDataPer == periodMon ) THEN
        CALL timeDate( vegTime,vegDate,1,'mon',l_360,tmpTime,tmpDate,'veg_update' )
      ENDIF
      vegTime = tmpTime
      vegDate = tmpDate

      IF ( resetStep ) THEN
!       Resetting time and date of data,for spin up.
!       Reuse initial values.
        vegTime = vegTimeInit
        vegDate = vegDateInit
!       Set switch to trigger a search for these data.
        nextData = .FALSE.
      ENDIF

!     Work out time level and file that hold the required data.
      CALL find_file( vegTime,vegDate,vegClim,vegEndTime,nextData,vegDataPer,NINT(timeStep)  &
                     ,vegFile,vegFileStep  &
                     ,vegTemplateT,vegFilePer,vegFileTime,vegFileDate  &
                     ,vegFileName,vegFileTemplate,vegTemplateDate  &
                     ,vegTemplateTime,vegTemplateUnits,'veg_update' )

!     Get length of period of validity of data (i.e. set vegDataStepMax).
      vegDataStepMax = vegDataPer
      IF ( vegDataPer == periodMon ) THEN

!       The following code is very similar to some in data_init, so bugs in one may
!       well be in the other, and code might better be put into a function.

!       At this point, vegDate gives the date of the data that will be read next,
!       but we may need to use an earlier month - depends on interpolation.
!       vegTimeIndex can be used to calculate how many times of data are used.
!       We then take vegDate and go back by the number of times (months) that are used.
!       Some of this code is close to or actually circular...but it's the best I can do for now.

!       Earlier we have restricted monthly veg to have interpFlags i,nb,nc or nf,
!       and this code is written assuming that holds - check via index.
        IF ( ANY(vegTimeIndex(:)<0) .OR. ANY(vegTimeIndex(:)>1) ) THEN
          WRITE(*,*)'ERROR: data_init: monthly data with unexpected datTimeIndex,'
          WRITE(*,*)'suggesting unexpected interpolation flag.'
          WRITE(*,*)'New code required.'
          STOP
        ENDIF
!       Now we know that we only need times 0 and 1, which turns out to mean that we
!       only need to know the length of the current month (actually from 1 month from
!       data time). Start from time of the last data read, and work back.
!       Work out required offset from vegTimeIndex(2).
        dt = 0
        IF ( vegTimeIndex(2) == 1 ) dt = 1
!       Take latest data time, and wind back by dt months.
        CALL timeDate( 0,vegDate,-dt,'mon',l_360,tmpTime,tmpDate2,'veg_update')
!       Now get date one month later than this (i.e. end of data validity).
!       I'm assuming (as elsewhere) that the same day of month can be used in both months.
        CALL timeDate( 0,tmpDate2,1,'mon',l_360,tmpTime,tmpDate,'veg_update')
!       Get difference between times (the length of first month).
        CALL timeDate_diff( 0,tmpDate2,0,tmpDate,l_360,'veg_update',nsec,nmin,nhr,ndy )
        vegDataStepMax = ( ndy*iSecInDay + nhr*iSecInHour + nmin*iSecInMin + nsec ) / NINT(timeStep)

      ENDIF  !  periodMon

!     Get number of timesteps until veg is next updated.
      vegUpdateStepMax = vegUpdatePer
      IF ( vegUpdatePer == periodMon ) THEN
!       Monthly updates - assume monthly data.
        IF ( vegDataPer /= periodMon ) THEN
          WRITE(*,*)'ERROR: veg_update: monthly updates only allowed for monthly data.'
          STOP
        ENDIF
        vegUpdateStepMax = vegDataStepMax
      ENDIF
!XX   Need to reduce if necessary to fit with spin up?

!--------------------------------------------------
    ELSE

!     NOT vegVaryT. Veg properties are f(PFT,x).
!     We come here during initialisation only. Use the first (only) time in the only file.
!     Some or all of this was done in init_veg_vary, but current code has incremented
!     soem counters...so I'm resetting here.
      vegFile = 1
      vegFileStep = 1
      vegDataStep = 0        !  vegDataIn was allocated for times 0:0
      vegDataStepMax = 0
      vegUpdateStep = 0
      vegUpdateStepMax = 0   !  so we trigger an update when vegUpdateStep is incremented

    ENDIF  !  vegVaryT
!-------------------------------------------------------------------------------

!   Read veg data.
    CALL veg_read

  ENDIF  !  vegDataStep (time to read data)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! Increment counter of steps between updates.
  vegUpdateStep = vegUpdateStep + 1

! Decide if it is time to update veg data.
  IF ( vegUpdateStep >= vegUpdateStepMax+1 ) THEN

    vegUpdateStep = 1
!XX Should probably reset vegUpdateStepMax here - to cover next period?

!   Do any time interpolation and load the veg data into the final variables.

!   Check that the time level of vegDataIn requested is acceptable.
    IF ( vegVaryT ) THEN
      istep = vegDataStep
      IF ( PRESENT(istepArg) ) istep=istepArg
      IF ( istep<1 .OR. istep>vegDataStepMax ) THEN
        WRITE(*,*)'ERROR: veg_update: istep is out of range.'
        WRITE(*,*)'Value=',istep,' allowed range=1 to vegDataStepMax=',vegDataStepMax
        IF ( PRESENT(istepArg) ) THEN
          WRITE(*,*)'The most likely way for this to happen is if the actual argument corresponding'
          WRITE(*,*)'to the optional dummy argument istepArg is wrongly specified.'
        ENDIF
        STOP
      ENDIF
      vegDataStep = istep
    ENDIF  !  vegVaryT

!   All variables are optional - may or may not be read in.
    DO ivar=1,nvegVarMax

!     Establish what field this is.
      SELECT CASE ( ivar )
        CASE ( 1 )
          varNum = varNumCanht
        CASE ( 2 )
          varNum = varNumLAI
        CASE ( 3 )
          varNum = varNumRootd
        CASE default
          WRITE(*,*)'ERROR: veg_update: no code for ivar=',ivar
          STOP
      END SELECT

      IF ( varNum < 0 ) CYCLE  !  nothing to do for this variable

      SELECT CASE ( vegVarFlag(varNum) )
        CASE ( vegVarFlagPT )
!         Interpolate in time and copy input field (PFT) to (PFT,x)
!         or, in the case of rootd_ft, copy input f(PFT) to f(PFT).
          IF ( varNum == varNumRootd ) THEN
            CALL veg_move_p_to_p( varNum )
          ELSE
            CALL veg_move_p_to_px( varNum )
          ENDIF
        CASE ( vegVarFlagPTX )
!         Interpolate in time and copy input field (PFT,x) to (PFT,x).
          CALL veg_move_px_to_px( varNum,.TRUE. )
        CASE ( vegVarFlagPX )
!         Copy input field (PFT,x) to (PFT,x).
          CALL veg_move_px_to_px( varNum,.FALSE. )
        CASE default
          WRITE(*,*)'ERROR: veg_update: no code for vegVarFlag=',vegVarFlag(varNum)
          STOP
      END SELECT

    ENDDO  !  ivar

!   Calculate parameters (characteristics) for surface types.
!   Don't do this during initialisation (tstep<1) when tile_pts has not been set.
    IF ( tstep > 0 ) CALL SPARM ( LAND_PTS,NTILES,CAN_MODEL,L_AGGREGATE  &
                                 ,TILE_PTS,TILE_INDEX  &
                                 ,FRAC,CANHT_FT,LAI,SATCON,CATCH_SNOW,CATCH  &
                                 ,INFIL_TILE,Z0_TILE )

  ENDIF  !  vegUpdateStep ( time for update)

!-------------------------------------------------------------------------------
!  If an initial dump is required and vegetation characteristics vary with time,
!  the dump was not written during initialisation because (depending upon
!  time interpolation) the veg may not have been updated until the first timestep.
!  At start of 1st timestep, any updating is now complete, so write a dump now, if necessary.
!DSM <BRAMS05 nao possui history>      IF ( tstep==1 .AND. vegVaryT .AND. dumpFreq>1 ) CALL dump_io( .TRUE., dumpTypeArg='init' )

  END SUBROUTINE veg_update
!###############################################################################
!###############################################################################
! subroutine veg_read
! Internal procedure in update_mod.
! Read veg data from file, making sure the correct files are open.

  SUBROUTINE veg_read( )

  USE ancil_info, ONLY : &
!  imported scalars with intent(in)
    land_pts

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile,openFile  &
!  imported arrays with intent(inout)
    ,irecPrev

  USE inout, ONLY :  &
!  imported scalar parameters
    formatNc  &
!  imported scalars with intent(in)
   ,nxIn,nyIn  &
!  imported arrays with intent(in)
   ,mapInLand,pftUse

  USE misc_utils, ONLY :  &
!  imported procedures
    replaceTemplate

  USE nstypes, ONLY : &
!  imported scalars with intent(in)
    npft

  USE readwrite_mod, ONLY :  &
!  imported procedures
    readvar2dComp

  USE veg_io_vars, ONLY :  &
!  imported scalar parameters
     vegVarFlagPT  &
!  imported scalars with intent(in)
    ,nfieldVegFile,noNewLineVeg,nvegFileVar,nvegHeaderField,nvegHeaderFile  &
    ,nvegHeaderTime,nvegVar,vegFile,vegFileStep,vegFormat,vegTemplateV  &
!  imported arrays with intent(in)
    ,vegFileName,vegTimeIndex,vegUnit,vegVarFlag,vegVarName,vegVarNameFile  &
    ,vegvarPos,vegVarStash  &
!  imported arrays with intent(inout)
    ,vegDataIn,vegUnitFile

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    i,ipft,ivar   &!  work
   ,nfieldFile   &!  number of fields per time in a file
   ,nlineField   &!  work
   ,nxRead       &!  x size of input grid
   ,nyRead       &!  y size of input grid
   ,readField    &!  number of field to read
   ,readZ        &!  'z' level to be read from file
   ,t,unit       &!  work
   ,useIndex      !  index in irecPrev

  INTEGER ::  &!  local ARRAYS
    tmpMap(nfieldVegFile,2)

  REAL ::  &!  local SCALARS
    tmpval(nfieldVegFile)  !  work - used for input if all data are on one line

  CHARACTER(len=LEN(vegFileName)) ::  &!  local SCALARS
    fileName  !  work
!--------------------------------------------------------------------------------

!  PRINT*,'top veg_read'

! Move the existing veg data back one time (e.g. so new data become old data).
!  print*, 'Before shift vegDataIn(ivar=1)=',vegDataIn(1,:,:,:)
  vegDataIn(:,:,:,:) = EOSHIFT( vegDataIn(:,:,:,:),1,dim=4 )
!  print*, 'After shift vegDataIn(ivar=1)=',vegDataIn(1,:,:,:)
!  if ( nvegvar>1 ) print*, 'After shift vegDataIn(ivar=2)=',vegDataIn(2,:,:,:)

! Ensure that the correct file is open.
  fileName = vegFileName(vegFile)
  DO i=1,nvegFileVar
!   Get file name if using template. 2nd argument (doTime)=.FALSE. to indicate
!   that we are not changing time template strings.
    IF ( vegTemplateV ) fileName=replaceTemplate( vegFileName(vegFile)  &
                                  ,.FALSE.,'veg_read',replString=vegVarNameFile(i) )
    IF ( vegUnitFile(i) /= fileName ) THEN
!     Close any file that is open on this unit, and open next file.
      CALL openFile( 1,.TRUE.,vegUnit(i),'read',vegFormat,fileName,'old'  &
                    ,'veg_read',ncType )
      vegUnitFile(i) = fileName
    ENDIF
  ENDDO

! Read veg data into the latest time level.
  t = vegTimeIndex(2)

! NB at present, ascii file fields will be 1,2,3,...

  unit = vegUnit(1)

!-------------------------------------------------------------------------------
  IF ( .NOT. noNewLineVeg ) THEN

    DO ivar=1,nvegvar
      IF ( vegTemplateV ) unit=vegUnit(ivar)
!     Set index to use for irecPrev with netCDF files - this irecPrev isn't
!     changed, but need to keep index within bounds.
      useIndex = unit
      IF ( vegFormat == formatNc ) useIndex = 1

!     If input data are not f(x) (e.g. are f(pft,t) only), we are not reading
!     from the usual input grid, so need to use different values.
!     For now, just do in separate IF clauses.
!-------------------------------------------------------------------------------
      IF ( vegVarFlag(ivar) == vegVarFlagPT ) THEN
        nxRead = 1
        nyRead = 1
!       Create fields to replace mapIn.
        tmpMap(1,:) = 1

        DO ipft=1,npft
          nlineField = 0     !   will not attempt to read ASCII line-by-line
          CALL readVar2dComp( vegFileStep,pftUse(ipft),vegVarPos(ivar)+pftUse(ipft)-1  &
                       ,vegVarStash(ivar),irecPrev(useIndex),nfieldVegFile  &
                       ,nvegHeaderFile,nvegHeaderTime  &
                       ,nvegHeaderField,nlineField,nxRead,nyRead,unit  &
                       ,vegVarName(ivar)  &
                       ,tmpMap(1,1),tmpMap(1,2),vegFormat,vegDataIn(ivar,1:1,ipft,t)  &
                       ,'veg_read','veg_read',ncType)
        ENDDO

!-------------------------------------------------------------------------------
      ELSE

!       vegVarFlag /= vegVarFlagPT, so we read from the usual (spatial) input grid

        DO ipft=1,npft
          nlineField = 0     !   will not attempt to read ASCII line-by-line

          CALL readVar2dComp( vegFileStep,pftUse(ipft),vegVarPos(ivar)+pftUse(ipft)-1  &
                       ,vegVarStash(ivar),irecPrev(useIndex),nfieldVegFile  &
                       ,nvegHeaderFile,nvegHeaderTime  &
                       ,nvegHeaderField,nlineField,nxIn,nyIn,unit,vegVarName(ivar)  &
                       ,mapInLand(:),(/(i,i=1,land_pts)/),vegFormat  &
                       ,vegDataIn(ivar,:,ipft,t),'veg_read','veg_read',ncType)

        ENDDO  !  pft

      ENDIF   !  vegVarFlag

    ENDDO  !  ivar

!-------------------------------------------------------------------------------
  ELSE

!   noNewLineVeg.
!   No new line between fields - all fields are arranged across one or more lines(s) in the ASCII input file.
!   nxIn=nyIn=1. Read all values for all variables as one field.
!   nfieldVegFile=number of data per time=nvarInFile*npftInFile.
!   Create fields to replace mapIn.
    tmpMap(:,1) = (/ (i, i=1,nfieldVegFile) /)
    tmpMap(:,2) = tmpMap(:,1)
!    irecPrev(unit) = 2
!   Pretend all data in file are a single field, with nx=nfieldVegFile*npftInFile.
    readZ = 1         !   "z" level to read from file
    readField = 1     !   field number to read
    nfieldFile = readField    !   # of fields in file. Set to field number - OK while readT=1
    nlineField = 0    !   will not attempt to read ASCII line-by-line
    nyRead = 1        !   input has ny=1

    CALL readVar2dComp( vegFileStep,readZ,readField,vegVarStash(1),irecPrev(unit)  &
                   ,nfieldFile,nvegHeaderFile,nvegHeaderTime  &
                   ,nvegHeaderField,nlineField,nfieldVegFile,nyRead,unit  &
                   ,vegVarName(1)  &
                   ,tmpMap(:,1),tmpMap(:,2),vegFormat,tmpval(1:nfieldVegFile)  &
                   ,'veg_read','veg_read',ncType)

!    print*,'One line veg read ',tmpval(1:nfieldVegFile)
!   Extract the required fields.
    DO ivar=1,nvegvar
      DO ipft=1,npft
        i = vegVarPos(ivar) + pftUse(ipft) - 1
!        print*,'ivar=',ivar,' ipft=',ipft,' pftUse=',pftUse(ipft),' reading i=',i
        vegDataIn(ivar,1,ipft,vegTimeIndex(2)) = tmpVal(i)
      ENDDO
!      print*,'#####veg_read: noNewLine: ivar=',ivar,' vegDataIn=',vegDataIn(ivar,1,1:npft,vegTimeIndex(2))
    ENDDO

  ENDIF  !  noNewLineVeg

  END SUBROUTINE veg_read
!################################################################################
!###############################################################################
!###############################################################################
!###############################################################################
!###############################################################################
! subroutine veg_move_p_to_p
! Internal procedure in update_mod.
! Interpolate in time and put input field(PFT) onto output (PFT).

  SUBROUTINE veg_move_p_to_p( varNum )

  USE nstypes, ONLY : &
!  imported scalars with intent(in)
     npft

  USE pftparm, ONLY : &
!  imported arrays with intent(inout)
     rootd_ft

  USE veg_io_vars, ONLY :  &
!  imported scalars with intent(in)
     varNumRootd,vegDataStep,vegDataStepMax  &
!   imported arrays with intent(in)
    ,vegDataIn,vegTimeIndex,vegVarInterp

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    varNum   !   index of current variable

  INTEGER ::  &!  local SCALARS
    ipft,t1,t2  !  work

  REAL ::  &!  local ARRAYS
    vegData(npft)  !  work: data after time-interpolation
!-------------------------------------------------------------------------------

! Interpolate in time.
  t1 = vegTimeIndex(1)
  t2 = vegTimeIndex(2)

! Call tinterp with 1st arg (ntOut)=1, meaning only return values for 1 time,
! 2nd arg (nt)=vegDataStepMax, defining the period of validity of the input data.
! 3rd arg (t)=vegDataStep, meaning the one value is for this timestep within the period of validity of input data
  DO ipft=1,npft
    CALL tinterp( 1,vegDataStepMax,vegDataStep,vegTimeIndex,vegDataIn(varNum,1,ipft,t1:t2)  &
                           ,'0',vegVarInterp(varNum),vegData(ipft) )
  ENDDO

! Copy data to final variable.
  IF ( varNum == varNumRootd ) THEN

    rootd_ft(:) = vegData(:)

  ELSE

    WRITE(*,*)'ERROR: veg_move_p_to_p: no code for varNum=',varNum
    STOP

  ENDIF

!  print*,'veg_move_p_to_p: rootd_ft=',rootd_ft(:)

  END SUBROUTINE veg_move_p_to_p
!###############################################################################
!###############################################################################
!###############################################################################
! subroutine veg_move_p_to_px
! Internal procedure in update_mod.
! Interpolate in time and put input field(PFT) onto output (PFT,land_pts).

  SUBROUTINE veg_move_p_to_px( varNum )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts  &
!  imported arrays with intent(in)
   ,frac

  USE nstypes, ONLY : &
!  imported scalars with intent(in)
    npft

  USE prognostics, ONLY :  &
!  imported arrays with intent(inout)
    canht_ft,lai

  USE veg_io_vars, ONLY :  &
!   imported scalars with intent(in)
     varNumCanht,varNumLAI,vegDataStep,vegDataStepMax  &
!   imported arrays with intent(in)
    ,vegDataIn,vegTimeIndex,vegVarInterp

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    varNum   !   index of current variable

  INTEGER ::  &!  local SCALARS
    i,ipft,t1,t2  !  work

  REAL ::  &!  local ARRAYS
    vegData(npft)  !  work: data after time-interpolation
!-------------------------------------------------------------------------------

! Interpolate in time.
  t1 = vegTimeIndex(1)
  t2 = vegTimeIndex(2)

! Call tinterp with 1st arg (ntOut)=1, meaning only return values for 1 time,
! 2nd arg (nt)=vegDataStepMax, defining the period of validity of the input data.
! 3rd arg (t)=vegDataStep, meaning the one value is for this timestep within the period of validity of input data
  DO ipft=1,npft
    CALL tinterp( 1,vegDataStepMax,vegDataStep,vegTimeIndex,vegDataIn(varNum,1,ipft,t1:t2)  &
                           ,'0',vegVarInterp(varNum),vegData(ipft) )
  ENDDO

! Create a spatial field according to frac.
  IF ( varNum == varNumCanht ) THEN
    canht_ft(:,:) = -1.0
    DO ipft=1,npft
      DO i=1,land_pts
        IF ( frac(i,ipft) > EPSILON(frac(:,:)) ) canht_ft(i,ipft)=vegData(ipft)
      ENDDO
    ENDDO

  ELSEIF ( varNum == varNumLAI ) THEN
    lai(:,:) = -1.0
    DO ipft=1,npft
      DO i=1,land_pts
        IF ( frac(i,ipft) > EPSILON(frac(:,:)) ) lai(i,ipft)=vegData(ipft)
      ENDDO
!      print*,'final lai for pft=',ipft,lai(:,ipft)
    ENDDO

  ELSE
    WRITE(*,*)'ERROR: veg_move_p_to_px: no code for varNum=',varNum
    STOP

  ENDIF

  END SUBROUTINE veg_move_p_to_px
!###############################################################################
!###############################################################################
! subroutine veg_move_px_to_px
! Internal procedure in update_mod..
! Interpolate in time if necessary and put input field(PFT,land_pts ) onto output (PFT,land_pts).

  SUBROUTINE veg_move_px_to_px( varNum,interp )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts

  USE nstypes, ONLY : &
!  imported scalars with intent(in)
    npft

  USE prognostics, ONLY :  &
!  imported arrays with intent(inout)
    canht_ft,lai

  USE veg_io_vars, ONLY :  &
!   imported scalars with intent(in)
     varNumCanht,varNumLAI,vegDataStep,vegDataStepMax  &
!   imported arrays with intent(in)
    ,vegDataIn,vegTimeIndex,vegvarInterp

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    varNum   !   index of current variable

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    interp  !  T means do time interpolation

  INTEGER ::  &!  local SCALARS
    i,ipft,t1,t2  !  work

  REAL ::  &!  local ARRAYS
    vegData(land_pts,npft)  !  work: data after time-interpolation
!-------------------------------------------------------------------------------

  IF ( interp ) THEN
!   Interpolate in time.
    t1 = vegTimeIndex(1)
    t2 = vegTimeIndex(2)
    DO i=1,land_pts
      DO ipft=1,npft
        CALL tinterp( 1,vegDataStepMax,vegDataStep,vegTimeIndex,vegDataIn(varNum,i,ipft,t1:t2)  &
                             ,'0',vegVarInterp(varNum),vegData(i,ipft) )
      ENDDO
    ENDDO
  ELSE
!   Simply copy data.
    vegData(:,:) = vegDataIn(varNum,:,:,vegDataStep)  !  vegDataStep will be zero
  ENDIF

! Copy to final variable.
  IF ( varNum == varNumCanht ) THEN
    canht_ft(:,:) = vegData(:,:)

  ELSEIF ( varNum == varNumLAI ) THEN
    lai(:,:) = vegData(:,:)

  ELSE
    WRITE(*,*)'ERROR: veg_move_px_to_px: no code for varNum=',varNum
    STOP

  ENDIF

  END SUBROUTINE veg_move_px_to_px
!###############################################################################
!###############################################################################

  END MODULE update_mod
!###############################################################################
!###############################################################################
