! module drive_io_vars
! Module containing variables relevant to (meteorological) driving data and their input.
!
! Note: the meteorological driving data are updated every timestep (although, depending
!   upon the selected flags, the values may be constant over several timesteps).
!
! Note: There are references to "weather generator" variables in here, but the
!   weather generator is not available in this version.

!-------------------------------------------------------------------------------

  MODULE drive_io_vars

  USE inout, ONLY :  &
!  imported scalar parameters
     formatLen

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: ndriveExtraMax = 10   ! the number of possible "extra"
!     driving variables.
!     The "extra" variables are provided to make it relatively easy to access
!     novel driving data while developing the model. The standard model requires
!     various driving data, but a new development might require some new variables.
!     Rather than having to understand all the code required to formally add new
!     driving variables, the user can access them via a standard route.
!     If you make this much larger (so that number of digits increases, the
!     formats used for character variables (e.g. drivevarName in INIT_DRIVE) will
!     also have to be able to cope with more digits.
  INTEGER, PARAMETER :: ndriveVarMax = 22 + ndriveExtraMax   ! the number of
!     possible driving variables
!     These are the "normal" variables + ndriveExtraMax "extra" variables
!     This may be larger than the number needed for any particular configuration.
!     This number includes the ndriveExtraMax "extra" variables.
!     Each possible variable should have an iposXX variable set below - except
!     that the "extra" variables share a single ipos variable (iposExtra).

!  REAL, PARAMETER :: pstarConst = 1.0e5  !  constant value of surface pressure for
!                                           possible use with weather generator (Pa)

!-------------------------------------------------------------------------------
! Scalar variables.
!-------------------------------------------------------------------------------

!--------------------------------------------------
!!  NOTE: The following driveXInit variables refer to the latest driving data that were read at the
!!        start of the run. They could all be recalculated rather than stored....but right now it's
!!        much easier if I just store them.
!   In fact...I think these variables refer to the first time read, not the last time!

  INTEGER ::  &
    driveDataStepInit    &!  the value of driveDataStep used at the start of the run
   ,driveDateInit       &!  the values of driveDate used at the start of the run
   ,driveTimeInit        !  the values of driveTime used at the start of the run
!--------------------------------------------------

  INTEGER ::  &
    driveResetStep      &!  the timestep number (a_step) when time and date of driving data will be reset,
                         !    to account for spin up. This is not used if nspin<0 (model-determined spin up).
!XX Is this required, now that only (what was previously) nspin<0 case is allowed???
   ,driveResetStepPrev   !  previous value of driveResetStep

  INTEGER ::   &
    driveDataPer       &!  period of input driving data (# of timesteps)
   ,driveDataStep      &!  counter of timesteps within a period of input driving data
!                           i.e. index of currently-used timestep, in a vector of driveDataPer times
!                           i.e. counter of how long until next data are read
   ,driveDate          &!  date (yyyymmdd) of the last driving data that were read (furthest
!                              forward in time)
   ,driveTemplateDate  &!  date associated with time-templated driving files
   ,driveTemplateTime  &!  time of day (s) associated with time-templated driving files
   ,driveTime          &!  time of day (s) of the last driving data that were read (furthest
!                            forward in time)
   ,driveFile          &!  the number (index) of the driving data file that is currently in use
   ,driveFilePer       &!  period of driving data files (i.e. interval between files) (# of timesteps)
!                            >0  period
!                            -1 (periodMon) monthly files
!                            -2 (periodAnn) annual files
!                            -9 (periodOneFile) all data times in a single file
   ,driveFileStep      &!  index of the time level last read from a driving file (i.e. 1,2,3,..)
!                            This index refers to driveTime, driveDate.
   ,ndriveExtra        &!  the number of "extra" driving variables used in a run
   ,ndriveFileTime     &!  the number of driving data files, each for a separate time period
!                             If each variable is in a separate file, but the single file for each variable
!                             holds all times of data, ndriveFileTime=1.
   ,ndriveUnit         &!  the number of driving data files required at any one time
!                            =1 if all driving variables come from the same file
!                            =ndriveVarIn if each driving variable comes from a separate file
!                            Can be >1 AND <ndriveVarIn if some files contain more than one variable.
   ,ndriveHeaderField  &!  number of header records before each field (xy slice) in a driving data file
!                            For an ASCII file, this is number of lines.
!                            For a GrADS file, this is number of data.
   ,ndriveHeaderFile   &!  number of header records at start of a driving data file
!                            For an ASCII file, this is number of lines.
!                            For a GrADS file, this is number of data.
   ,ndriveHeaderTime   &!  number of header records before each time in a driving data file
!                            For an ASCII file, this is number of lines.
!                            For a GrADS file, this is number of data.
   ,ndriveVar          &!  the number of driving variables that are held in driveData array
!                            Often this equals ndriveVarIn, but can differ if a variable is to
!                            be dervied from other input variables.
   ,ndriveVarIn        &!  the number of driving variables that are read from file
!                            This can be /= ndriveVar if some of the required variables are derived
!                            from input variables, e.g. snowfall might be required, but could be
!                            calculated from total precipitation.
   ,nfieldDriveFile    &!  the number of fields (xy planes) per time in a driving
!                            data file (not including headers or timestamp)
   ,ndriveDataTime     &!  number of time levels of driving input data that are stored
   ,ioPrecipType       &!  flag indicating how precipitation is input
!           1 = total precipitation is read in
!           2 = values for rainfall and snowfall are read in
!           3 = values for large-scale rainfall, convective rainfall and large-scale
!                 (assumed to be total) snowfall are read in
!           4 = values for convective rainfall, large-scale rainfall, convective snowfall and
!                 large-scale snowfall are read in.
!           5 = as 1, but distinct in that it is used with the weather generator. The
!                 distinction is useful as the components are held separately in this case.
!               Not allowed in this version!
   ,io_rad_type         ! Flag indicating how radiation is input.
!                         1 downward fluxes provided
!                         2 net (all wavelength) longwave flux and downward
!                             shortwave flux are provided
!                         3 net downward fluxes are provided
!-------------------------------------------------------------------------------
! The following iposXX  are indices that locate the driving variables in the list
! of all possible driving variables. There should be ndriveVarMax indices, one
! for each possible driving variable.

  INTEGER ::  &
     iposDiffRad      &!  position of diffuse radiation
    ,iposExtra        &!  position of first "extra" variable
    ,iposLWD          &!  position of downward longwave radiation
    ,iposLWN          &!  position of net downward longwave radiation
    ,iposPrecip       &!  position of total precipitation
!      Precip components use suffixes: 1st letter C convective, L large-scale, T total (C+L)
!                                    : 2nd letter R rain (liquid), S snow (solid)
!      Total precip (iposPrecip, above) is the exception to these "rules"!
    ,iposPrecipCR     &!  position of convective rainfall
    ,iposPrecipCS     &!  position of convective snowfall
    ,iposPrecipLR     &!  position of large-scale rainfall
    ,iposPrecipLS     &!  position of large-scale snowfall
    ,iposPrecipTR     &!  position of total rainfall
    ,iposPrecipTS     &!  position of total snowfall
    ,iposPstar        &!  position of surface pressure
    ,iposQ            &!  position of specific humidity
    ,iposRN           &!  position of net downward all-wavelength radiation
    ,iposSubSurfRoff  &!  position of subsurface runoff rate
    ,iposSurfRoff     &!  position of surface runoff rate
    ,iposSWD          &!  position of downward shortwave radiation
    ,iposSWN          &!  position of net downward shortwave radiation
    ,iposT            &!  position of air temperature
    ,iposU            &!  position of zonal wind
    ,iposV            &!  position of meridional wind
    ,iposWind         &!  position of (total, horizontal) windspeed
    ,iposOzone         !  position of ozone concentration

  REAL ::  &
    diffFracConst  &!  a constant value for fraction of radiation that is diffuse
!   ,durConvRain &!  duration of convective rainfall events (hours)
!   ,durLSRain   &!  duration of large-scale rainfall events (hours)
!   ,durLSSnow   &!  duration of large-scale snowfall events (hours)
!   ,tAmplConst  &!  constant value of amplitude of diurnal
!                   temperature variation for possible use with weather generator (K)
   ,tForCRain   &!   air temperature (K) at or above which rainfall is taken to
!                     be convective (rather than large-scale) in nature.
!                     Used with precipType=1, 2 or 5.
   ,tForSnow     !   air temperature (K) at or below which precipitation is assumed to be snow
!                     Used with precipType=1, 4 or 5.

  LOGICAL ::  &
    byteSwapDrive     &!  Switch indicatng that byte order of driving data is to be reversed
!                           Only affects binary (non-SDF) files.
!                           T means swap order after reading data from file
!                           F means leave order as in input file
   ,driveEndTime      &!  Switch indicating file naming convention, if driveTemplateT=TRUE.
!                           T means files are named according to end (last) time in file
!                            e.g. hourly data in daily files, data for 1 Jan are in file timed 2 Jan.
!                           F means files are named according to first time in file (as in GrADS)
!                            e.g. hourly data in daily files, data for 1 Jan are in file timed 1 Jan.
   ,driveTemplateT    &!  T means that the driving data is split into files according to time,
!                           with file naming convention given by a template
   ,driveTemplateV    &!  T means that each driving variable (for a given time interval)
!                           is in a separate file, with the name of each file following a template
!                           At present this is the only way to get each var in a separate file
!                              - i.e. must use template.
   ,noNewLineDrive    &!  T means that all driving variables for each time are arranged across
!                           one or more lines in
!                           an ASCII input file, WITHOUT each variable starting on a new line.
!                           This is only allowed if there is only 1 x,y point in the input file.
!                         F means that each driving variable starts a new line in an ASCII input file.
    ,notNextDrive     &!  used to force the drive_update to search for a file with data
!                           for the appropriate time,
!                           rather than assuming that "next" file can be used.
!                           Needed when no data are read by init_drive, so that first call
!                           to drive_update
!                           is from the main program, where the assumption is that next=TRUE.
   ,useWGen            !  T means the weather generator is used to disaggregate in time
!                         F means no weather generator is used

  CHARACTER(len=formatLen) :: driveFormat  !  format (type) of driving data file

  CHARACTER(len=2) :: driveTemplateUnits   !  indicates what time units are used
!                                             to represent time templating

  CHARACTER(len=150) :: driveFileTemplate  !  template file name for driving data

!-------------------------------------------------------------------------------
! Array variables.
!-------------------------------------------------------------------------------
  INTEGER ::   &
    driveTimeIndex(2)        &!  index showing range of times of driving data held
!          index(1) and index(2) are the earliest and latest time respectively
!          Values are relative to the current time when the data are read, and are given
!          in terms of the period of the input data
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)
   ,driveVarFlag(ndriveVarMax)   &!  flag indicating where to find each driving variable in input
!                                 >0 gives the location (field number) in file
!                                 <-1 means derive variable from others
   ,driveVarPos(ndriveVarMax)    &!  location in driveData array where each stored driving variable
!                                      is stored.
!                                      <1 means variable is not stored - either it's not used or
!                                      it's read in but not stored  - e.g. because it is used to
!                                      derive another variable.
   ,driveVarPosIn(ndriveVarMax)  &!  location in driveDataIn array where each input driving variable
!                                      is stored. <1 means variable is not input.
   ,driveVarStash(ndriveVarMax)  &!  STASH code for each driving variable that is input
   ,driveVarInSort(ndriveVarMax) &!  a list of which variables in master list are to be read in,
!                                      sorted into ascending order of location in file.
!                                      Values 1:ndriveVarIn give a list of variables that are to be
!                                      read from file (in terms of location in master list of all possible
!                                      driving variables), sorted into ascending order of location
!                                      of data in file (which makes ASCII access quicker).
!                                      Hence a loop over driveVarInSort(1:ndriveVarIn) identifies
!                                      variables in the master list that are to be read in.
   ,driveUnitUse(ndriveVarMax)    !  the index in driveUnit to use for each input variable
!                                      This allows each input variable to be associated with a unit.

  INTEGER, ALLOCATABLE ::   &
    driveFileDate(:)  &!  the date of first data in each driving file (not used if driveTemplateT=T)
   ,driveFileTime(:)  &!  the time of first data in each driving file (not used if driveTemplateT=T)
   ,driveUnit(:)       !  the units used to connect to driving data files
!                           For a netCDF file, this is the netCDF ID.

  REAL, ALLOCATABLE ::  &
    driveData(:,:,:,:)  &!  driving data for current driving data period, after temporal
!                             interpolation/disaggregation
!                        !    Only variables that are read in are held (not any derived variables).
!                        !    If the period of data equals model timestep (so no time-interpolation) and
!                        !    no temporal disaggregation, driveData becomes a rather unnecessary
!                        !    copy of the input driveDataIn.
!                        !  Dimensions are: 1 variables, 2 x, 3 y, 4 times
   ,driveDataIn(:,:,:,:) !    driving data as read in (i.e. only those variables that are read in)
!                        !  Dimensions are: 1 variables, 2 x, 3 y, 4 times
!                           NOTE that driveData and driveDataIn are dimensioned with nx,ny. i.e. for full
!                           model grid. For runoff variables, which are only defined on land points, this
!                           may be much more space than is actually needed.

  LOGICAL :: driveVarUse(ndriveVarMax)     !  flag indicating if a variable is required for the chosen
!                                      model configuration. A variable can be needed (driveVarUse=T)
!                                      but not be input if it can be derived from
!                                      variables that are input.

  CHARACTER(len=2) ::  &
    driveVarInterp(ndriveVarMax)     !  code describing time interval over which input data applies and
!                                    !  how data are used/interpolated
!          'b' backward time average, i.e. time average ending at given time (in GSWP2 this is L)
!          'c' centred time average, i.e. time average centred on given time (GSWP2 C)
!          'f' forward time average, i.e. time average starting at given time (GSWP2 N)
!          'i' instantaneous value at given time (interpolation will be linear in time)(GSWP2 I)
!          'nb' no interpolation, value is valid over time interval ending at given time (not in GSWP2)
!          'nc' no interpolation, value is valid over time interval centred on given time  (GSWP2 0)
!          'nf' no interpolation, value is valid over time interval starting at given time (GSWP2 default)
!          Note that interpolation for b,c and f will generate values that, generally, line outside the range of the input averages.

  CHARACTER(len=25) :: allDriveVarName(ndriveVarMax)  !  names of all possible driving variables
  CHARACTER(len=25) :: driveVarNameSDF(ndriveVarMax)  !  names of driving variables (as used in
!                                           self-describing files,eg netCDF files)

  CHARACTER(len=30), ALLOCATABLE :: driveVarNameUnit(:)   !  names of driving
!                       variables as used in file names - i.e. that part
!                       of file name that is replaced if using template naming. There is
!                       one value for each unit (each file) open at any one time.


  CHARACTER(len=150), ALLOCATABLE :: driveFileName(:)    !  names of driving
!                                                        !  data files
!        If time-templating is used (driveTemplateT), there is only one name,
!        the template. If templating is used for variable name (driveTemplateV),
!        that part is left in driveFileName.
  CHARACTER(len=150), ALLOCATABLE :: driveFileOnUnit(:)  !  the name of the
!        driving data file currently open on each unit

  END MODULE drive_io_vars
!###############################################################################
!###############################################################################
