! module veg_io_vars
! Module containing variables relevant to updating prescribed vegetation fields that vary with time and/or space.
!
! Note: There are two timescales/periods that are relevant to time-varying vegetation data:
!  1) period/frequency of input data (vegDataPer) - e.g. read a new value once a month
!  2) period/frequency at which the value used in the model is updated - e.g. update daily.
!  Period #1 must be >= #2

!-------------------------------------------------------------------------------

  MODULE veg_io_vars

  USE inout, ONLY :  &
!  Imported scalar parameters
     formatLen

  IMPLICIT NONE

!-------------------------------------------------------------------------------

! Global (i.e. public) declarations.

! Global scalar parameters
  INTEGER, PARAMETER ::  &
     nvegVarMax = 3    !  Maximum possible number of prescribed time or position varying veg vars.
!                           Exists to ease dimensioning of arrays.
!                           The 3 are canht, LAI and rootdepth.
!                           Each variable should have a varNumXX entry below.

  CHARACTER(len=2), PARAMETER ::  &
     vegVarFlagPT = 't'      &!  value of vegVarFlag that indicates that variable is f(PFT,t) only
    ,vegVarFlagPTX = 'tx'    &!  value of vegVarFlag that indicates that variable is f(PFT,t,x)
    ,vegVarFlagPX = 'x'       !  value of vegVarFlag that indicates that variable is f(PFT,x) only

!-------------------------------------------------------------------------------
!  Global scalars

  INTEGER ::  &!  SCALARS
!!  NOTE: The following vegXInit variables refer to the latest veg data that were read at the
!!        start of the run. They could all be recalculated rather than stored....but right now it's
!!        much easier if I just store them.
!   In fact...I think these variables refer to the first time read, not the last time!
    vegDataStepInit    &!  the value of vegDataStep used at the start of the run
   ,vegDateInit        &!  the values of vegDate used at the start of the run
   ,vegTimeInit        &!  the values of vegTime used at the start of the run
   ,vegUpdateStepInit   !  the value of vegUpdateStep used at the start of the run


   INTEGER ::       &!
    nfieldVegFile   &!  the number of fields (xy planes) per time in the veg data file
                     !     (not including headers or timestamp). Each PFT is considered a separate field.
   ,nvegFileTime    &!  the number of veg data files, each for a separate time period
!                          If each variable is in a separate file, but the single file for each variable 
!                          holds all times of data, nvegFileTime=1.
   ,nvegFileVar     &!  the number of veg data files required at any one time
!                    !    =1 if all veg variables come from the same file
!                    !    =nvegVar if each veg variable comes from a separate file
   ,nvegVar         &!  the number of prescribed time or space-varying veg variables
   ,nvegDataTime    &!  number of time levels of veg input data that are stored
   ,nvegHeaderField  &!  number of header records before each field (xy slice) in a veg data file
                      !    For an ASCII file, this is number of lines.
                      !    For a GrADS file, this is number of data.
   ,nvegHeaderFile   &!  number of header records at start of a veg data file
                      !    For an ASCII file, this is number of lines.
                      !    For a GrADS file, this is number of data.
   ,nvegHeaderTime   &!  number of header records before each time in a veg data file
                      !    For an ASCII file, this is number of lines.
                      !    For a GrADS file, this is number of data.
!-------------------------------------------------------------------------------
! The following are indices that locate the input veg variables in variables such as vegVarIn.
! They do not say where the data are in the file - that is given by the values of vegVarPos.
! e.g. varNumLAI=2 with vegVarPos(2)=3 indicates that LAI is held in the 2nd position
!      of variables, and the data are the 3rd field in input files.
    ,varNumCanht        &!  position of canopy height
    ,varNumLAI          &!  position of LAI
    ,varNumRootd        &!  position of root depth
!-------------------------------------------------------------------------------
    ,vegDataPer   &!  period of input veg data (# of timesteps)
    ,vegDataStep  &!  counter of timesteps within a period of input veg data
!                       i.e. counter of how long until next data are read
    ,vegDataStepMax  &!  number of timesteps over which current input veg data remain valid
!        For most cases vegDataStepMax=vegDataPer, but special cases are:
!        1) monthly data - vegDatStepMax depends upon length of month
!        2) if time "jumps" due to spin up - interval over which veg data are used is shortened
    ,vegDate      &!  date (yyyymmdd) of the last veg data that were read (furthest forward in time)
    ,vegTemplateDate  &!  date associated with time-templated driving files
    ,vegTemplateTime  &!  time of day (s) associated with time-templated driving files
    ,vegUpdateStep      &!  counter of timesteps within veg update period (i.e. counter for time to next update)
    ,vegUpdateStepMax   &!  maximum possible value of vegUpdateStep
    ,vegTime      &!  time of day (s) of the last veg data that were read (furthest forward in time)
    ,vegTimeLev   &!  time level (index) corresponding to vegTime
    ,vegFile      &!  the number (index) of the veg data file that is currently in use
    ,vegFilePer   &!  period of veg data files (i.e. interval between files) (# of timesteps)
!                      >0  period
!                      -1 (periodMon) monthly files
!                      -2 (periodAnn) annual files
!                      -9 (periodOneFile) all data times in a single file
   ,vegFileStep       &!  index of the time level last read from a vegfile (i.e. 1,2,3,..)
!                         This index refers to vegTime, vegDate.
   ,vegResetStep      &!  the timestep number (a_step) when time and date of veg data will be reset,
                       !    to account for spin up. This is not used if nspin<0 (model-determined spin up).
!XX Is this required, now that only (what was previously) nspin<0 case is allowed???
   ,vegResetStepPrev  &!  previous value of vegResetStep
   ,vegUpdatePer       !  period for update of time-varying veg fields (# of timesteps)

  LOGICAL ::       &!
    noNewLineVeg   &!  T means that all veg variables for each time are arranged across one or more lines in
!                        an ASCII input file, WITHOUT each variable starting on a new line.
!                        This is only allowed if there is only 1 x,y point in the input file.
!                      F means that each veg variable starts a new line in an ASCII input file.
   ,notNextVeg     &!  Used to force veg_update to search for a file with data for the
!                        appropriate time, rather than assuming that "next" file can be used.
!                        Needed when no data are read by init_veg, so that first call to veg_update
!                        is from the main program, where the assumption is that next=TRUE.
   ,vegClim        &!  T means that veg data are "climatological", meaning
!                        the same data are reused for every year
   ,vegEndTime      &!  Switch indicating file naming convention, if vegTemplateT=TRUE.
!                           T means files are named according to end (last) time in file
!                            e.g. hourly data in daily files, data for 1 Jan are in file timed 2 Jan.
!                           F means files are named according to first time in file (as in GrADS)
!                            e.g. hourly data in daily files, data for 1 Jan are in file timed 1 Jan.
   ,vegTemplateT   &!  T means that the veg data is split into files according to time, with file naming
!                          convention given by a template
   ,vegTemplateV   &!  T means that each vegvariable (for a given time interval) is in a separate file,
!                          with the name of each file following a template
!                         At present this is the only way to get each var in a separate file - i.e. must use template.
   ,vegVarSeparate &!  T means that each veg variable (for a given time interval) is in a separate file
!                         F means there all variables (for a given time interval) are in one file
   ,vegVaryT        !  T means that one or more prescribed veg field varies with time

  CHARACTER(len=formatLen) ::  &!
    vegFormat       !  format (type) of veg data file

  CHARACTER(len=2) :: &!  SCALARS
    vegTemplateUnits   !  indicates what time units are used to represent time templating

  CHARACTER(len=150) ::  &!  scalars
    vegFileTemplate   !  template file name for veg data
!--------------------------------------------------------------------------------

! Global arrays

  INTEGER ::   &!
    vegTimeIndex(2)          !  index showing range of times of data held
!          index(1) and index(2) are the earliest and latest time respectively
!          Values are relative to the current time, and are given in terms of the period of the input data
!         -1 = previous datum
!          0 = current datum
!          1 = next datum
!          2 = next datum after next! (i.e. two ahead)

  INTEGER, ALLOCATABLE ::   &!
    vegFileDate(:)  &!  the date of first data in each veg file (not used if vegTemplateT=T)
   ,vegFileTime(:)  &!  the time of first data in each veg file (not used if vegTemplateT=T)
   ,vegUnit(:)      &!  the units used to connect to veg data files
   ,vegVarPos(:)    &!  the location (field or xy plane number) of each veg variable in the input file
   ,vegVarStash(:)   !  the STASH code (from the UM) to use for each veg variable

  REAL, ALLOCATABLE ::  &!
    vegDataIn(:,:,:,:)   !  veg data as read in (i.e. only those variables that are read in)
                         !    Dimensions are: 1 variables, 2 point, 3 PFT, 4 time.

  CHARACTER(len=2), ALLOCATABLE :: vegVarFlag(:)   !  flag indicating how variable is modelled
!        Length of this must match that of vegVarFlagPT above.
!        gfortran complains if try to set length using len(vegVarFlagPT) here.

  CHARACTER(len=2), ALLOCATABLE ::  &!
    vegVarInterp(:)     !  code describing time interval over which input data applies and
!                                    !  how data are used/interpolated
!          'b' backward time average, i.e. time average ending at given time (in GSWP2 this is L)
!          'c' centred time average, i.e. time average centred on given time (GSWP2 C)
!          'f' forward time average, i.e. time average starting at given time (GSWP2 N)
!          'i' instantaneous value at given time (interpolation will be linear in time)(GSWP2 I)
!          'nb' no interpolation, value is valid over time interval ending at given time (not in GSWP2)
!          'nc' no interpolation, value is valid over time interval centred on given time  (GSWP2 0)
!          'nf' no interpolation, value is valid over time interval starting at given time (GSWP2 default)
!          Note that interpolation for b,c and f will generate values that, generally,
!          lie outside the range of the input averages (but that conserve the average).

  CHARACTER(len=25), ALLOCATABLE ::  &!
    vegVarName(:)      &!  names of veg variables (as used in eg netCDF files)
   ,vegVarNameFile(:)   !  names of veg variables as used in file names (if each var from separate file)

  CHARACTER(len=150), ALLOCATABLE ::  &!
    vegFileName(:) &!  names of veg data files
!          If time-templating is used (vegTemplateT), there is only one name, the template.
!          If templating is used for variable name (vegTemplateV), that part is left in vegFileName.
   ,vegUnitFile(:)  !  the name of the veg data file currently open on each unit

  END MODULE veg_io_vars
!###############################################################################
!###############################################################################
