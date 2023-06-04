! module inout
! Contains details of input and output (including dumps).

  MODULE inout

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER ::  &
    formatLen = 3      &!  length of character variable used for file formats
   ,longLen = 90       &!  max possible length for variable descriptions
   ,nameLen = 15       &!  max possible length for variable names
   ,ntemplString = 16  &!  number of template strings
   ,optLen = 10        &!  length used for many model options
!-------------------------------------------------------------------------------
!  The following 'period' values should be distinct and less than zero.
!-------------------------------------------------------------------------------
   ,periodAnn = -2     &!  the value given to periods to indicate an annual period
   ,periodMon = -1     &!  the value given to periods to indicate a monthly period
   ,periodOneFile = -9 &!  the value given to file periods to indicate that a
!                            single file is used
!-------------------------------------------------------------------------------
   ,sufLen = 20         !  max possible length for name appended to variable's description

!-------------------------------------------------------------------------------
! Units used for standard input and output (i.e. the units that
! would be used if an explicit unit number was omitted from a
! format, e.g. write(*,*) ). Apparently these can vary between systems.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER ::  &
    !DSM stdIn = 5     &!  standard input
    stdIn = 77     &!  standard input      !DSM para ler do arquivo: jules.in
   ,stdOut = 6     !  standard output

!-------------------------------------------------------------------------------
! Unit used to read JULES input file
! This will be obtained using fileUnit if a file is specified as an argument
! Otherwise will default to stdIn to allow input piping
!-------------------------------------------------------------------------------
  INTEGER ::  &
    jinUnit = stdIn

!-------------------------------------------------------------------------------
! Strings used to indicate file formats.
! These are used both within the code and as filename extensions (except that
! formatBin uses extension ".gra").
!-------------------------------------------------------------------------------
  CHARACTER(len=formatLen), PARAMETER ::  &
    formatAsc = 'asc'  &!  indicates ASCII file format
   ,formatBin = 'bin'  &!  indicates gridded 'flat' binary file format
!                             GrADS ("standard", not netCDF) files are a subset of this.
!                             This is often used as a synonym for 'GrADS format'.
   ,formatNc = 'nc'    &!  indicates netCDF file format
   ,formatPP = 'pp'     !  indicates Met Office PP file format

  CHARACTER(len=7), PARAMETER ::  &
!   These "tag" variables are used to locate the start of sections of the run
!   control file that give information specific to a file format.
    tagAscBin = '>ASCBIN'  &!  tag for information for ASCII or binary
   ,tagNc = '>NC'           !  tag for netCDF information

!-------------------------------------------------------------------------------
! Array parameters.
!-------------------------------------------------------------------------------
  character(len=1), parameter ::  &
    metOChar(0:35) =  &!  characters in order used in Met Office file naming convention
                   (/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'  &
                     ,'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'  &
                     ,'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't'  &
                     ,'u', 'v', 'w', 'x', 'y', 'z' /)

  CHARACTER(len=3), PARAMETER ::  &
    templString(ntemplString) =  (/  &!  strings that are used to indicate templating
!           in file names. These strings are replaced to get the name of a file.
!           Note that only some of these are recognised by GrADS,
     '%tc'     &!  1-character decade, as used by Met Office (not in GrADS)
    ,'%y2'     &!  2-digit year (in GrADS)
    ,'%y4'     &!  4-digit year (in GrADS)
    ,'%yc'     &!  1-character year in decade, as used by Met Office (not in GrADS)
    ,'%m1'     &!  1- or 2-digit month (in GrADS)
    ,'%m2'     &!  2-digit month (in GrADS)
    ,'%mc'     &!  3-character month abbreviation (in GrADS)
    ,'%mm'     &!  1-character month, as used by Met Office (not in GrADS)
    ,'%d1'     &!  1- or 2-digit day of month (in GrADS)
    ,'%d2'     &!  2-digit day of month (in GrADS)
    ,'%dc'     &!  1-character day of month, as used by Met Office (not in GrADS)
    ,'%h1'     &!  1- or 2-digit hour of day (in GrADS)
    ,'%h2'     &!  2-digit hour of day (in GrADS)
    ,'%hc'     &!  1-character hour of day, as used by Met Office (not in GrADS)
    ,'%n2'     &!  2-digit minute of hour (in GrADS)
    ,'%vv'     &!  name of a variable (not in GrADS)
       /)

!-------------------------------------------------------------------------------

  INTEGER ::  &!  SCALARS
    dumpFreq    &!  Flag indicating how often the model state is to be written to a dump
!                     0: don't write dumps
!                     1: dump final state only
!                     2: dump initial and final states
!                     3: As 2 but also after end of spin up phase
!                     4: As 3 but also at end of every calendar year
   ,nout          &!  number of output profiles selected
   ,npoints       &!  number of points in the grid (model domain)
   ,nvar          &!  the number of variables available
   ,nvarOutTot    &!  total # of variables to be output, summed over all output profiles
   ,nxIn          &!  number of points in the x-direction (row length) in input data
   ,nyIn          &!  number of points in the y-direction (# of rows) in input data
   ,outLen        &!  extent of the output vector (outval)
   ,outLenWrite   &!  extent of the output vector used to write to file (outvalWrite)
   ,pointsOutMax  &!  maxVal(pointsOut) = max number of points in any output profile
   ,print_Step    &!  Number of timesteps between time messages to screen
   ,x1In          &!  x index in the input grid of first point in the model grid
   ,y1In           !  y index in the input grid of first point in the model grid

!-------------------------------------------------------------------------------

  INTEGER, ALLOCATABLE ::  &
    irecPrevOut(:)  &!  record number for unformatted output files
   ,nlevMax(:)      &!  maximum number of levels in any variable to be output
!                       It is the maximum value of varNlev for any variable in
!                       a given output profile.
   ,nlevMaxCtl(:)   &!  maximum number of levels in any variable to be output -
!                           as represented in GrADS ctl file
   ,ntCtl(:)        &!  the number of times as initially written to GrADS ctl tdef
   ,ntCtlNeed(:)    &!  the number of times that should be in GrADS ctl tdef
!                        (i.e. count of times written to files)
   ,ntOutFilePer(:) &!  number of times expected in a file (number of outper in outFilePer)
!                         In some cases this is estimated, rather than working out exactly,
!                         such as if spin up etc makes life more complicated.
   ,ntOutPer(:)     &!  expected number of values in a time accumulation used for time average.
!                         This is the number of sampling periods in the output period.
!                         In some cases this is estimated, rather than working out exactly.
   ,nvarOut(:)      &!  number of variables in each output profile
   ,nxyMax(:)       &!  maximum nx*ny of fields for each output profile (these are nx,ny
!                         for input fields on model domain, not output nxout etc)
!                         Redundant while I insist that a profile only has one output grid.
   ,outAccum(:)     &!  time accumulated since last output (s)
   ,outDate(:,:)    &!  start and finish dates for each output profile (see also outTime)
   ,outDateFlag(:)  &!  flag indicating when output is to be generated. See also outDate.
!                         Used to deal with special cases when outDate may be confusing -
!                         e.g. output through spin up.
!                         >0: output is generated only after spin up, between dates
!                           and times given by outDate and outTime.
!                         0: output throught run (spin up, if there is any, and main run)
!                           (including start of first timestep)
!                         -1: output at all times after spin up (and there is spin up)
!                         -2: output the initial state only (at the start of first timestep)
   ,outEndPos(:)    &!  the last location in output vector (outlen) used by each profile
   ,outFileAccum(:) &!  time accumulated since last output file (s)
   ,outGridNxy(:,:) &!  x any y extents of output grid for each profile.
!                         2nd dimension: 1=nx 2=ny
!                         If the output is to be "compressed" (see outCompress), these
!                         give the shape of the full, uncompressed output grid.
!                         Otherwise, these give the shape of the actual output
!                         including padding, so that their product for a profile=outNpWrite.
   ,outPer(:)       &!  the period (number of timesteps) for each output profile
!                        >0  period
!                         0  use model timestep for period
!                        -1 (periodMon) monthly output
!                        -2 (periodAnn) annual output
   ,outPerNunits(:) &!  outPer expressed in terms of the units chosen.
!                       e.g. if outPer=2days, outPerNunits=2.
!                       At present this is used to increment a time index for netCDF files,
!                       and the units themselves are not available globally - just
!                       in subroutine newOutFilename!
   ,outFilePer(:)   &! the period (# of timesteps) for output files
!                        >0  period
!                         0 means use model timstep for period
!                        -1 (periodMon) means new file each month
!                        -2 (periodAnn)  means new file each year
!                        -7 means as -9 but output from each spin-up cycle to a separate file
!                        -8 means as -9 but all spin-up output times to a single separate file
!                        -9 (periodOneFile) means all output times to one file
!                            (NB all output times, not necessarily all times)
   ,outFileStep(:) &!  counter of timesteps for output files
!                      Note that under certain circumstances this can be zero even after
!                      the first time of data have been written - for reasons that are lost
!                      in the mists of time (but might have something to do with dealing
!                      with the case when output is requested at all times, incl t=0).
!                      See outWriteCount for a definitive count of times in a file!
   ,outNpWrite(:)  &!  number of values in a single field of each output profile
!                        This is length of vector that is written, including any padding.
!                        If compressing output (see outCompress), outNpWrite=pointsOut and
!                        is generally smaller than the full, uncompressed grid (which is sized
!                        as the product of outGridNxy).
   ,outSamPer(:)   &!  sampling period (# of timesteps) for time-average output
   ,outStep(:)     &!  counter of timesteps in current output interval
!                        This is initialised in such a way as to ensure that output is
!                        generated at times that are "synchronous" with some important
!                        period - e.g. 3-hour averages are generated for 0-3H, 3-6H etc
!                        (and NOT 1-4H, 4-7H, which would not be "in synchrony" with the days.
   ,outStepSamPer(:)  &!  counter of number of times that have been added to time accumulation
!                           for time-average of time-accumulated output. If the sampling period
!                           is longer than a timestep, this counter is incremented less
!                           often than outStep.
!                           We avoid simply calculating outStepSamPer from outStep, so as
!                           to correctly count times in accumulation under all circumstances.
   ,outTime(:,:)      &!  start and finish times of day for each output profile
!                           Output starts at or after outDate(:,1),outTime(:,1) and ends at or after
!                           outDate(:,2),outTime(:,2).
   ,outTimeID(:,:)    &!  netCDF IDs of time variables
!                         2nd dimension: 1=time, 2=timestep
   ,outUnit(:   )     &!  unit used to connect file for each output profile
   ,outVarID(:)       &!  netCDF ID for each variable
   ,outWriteCount(:)  &!  counter of number of times that have been written to an output file
   ,pointsFlag(:,:)   &!  flag indicating what points are to be output and to what grid
!             pointsFlag(1) indicates what points in the model grid are to be output
!               0 = all points in the grid are to be output
!               1 = points in a given area are to be output
!               2 = points to be output will be listed
!             pointsFlag(2) indicates what the output grid will look like
!                Note that whatever the output grid, it can be compressed (see outCompress).
!               0 = the output grid is the model grid
!               1 = the output grid will be exactly as given via pointsFlag(1)=1.
!                     i.e. flag(1)=1 indicates area to be output
!                          flag(2)=1 will output that area exactly
!               2 = the shape of the output grid and mapping from model grid will be
!                     specified directly (but mapping can be altered by compression).
!                     Most likely will be used with pointsFlag(1)=2,
!                     and arguably is less useful than most other options because everything
!                     has to be prescribed. However, this does mean that it should be
!                     able to act as a "do-anything" option that can be used to deal with
!                     currently unforeseen or unlikely cases....if you set up the input correctly.
!               3 = the output grid is the smallest possible regular, rectangular
!                   grid that encompasses all selected points
!               4 = the output grid is the input grid
!               5 = the output grid is a vector. This can only be used if pointsFlag(1)=2,
!                     i.e. a list of locations was read in, in which case it simply writes
!                     each location to the next point in the vector,
!                     e.g. chosen point #1 is written to point #1 in output
!                     This saves having to define a trivial mapOut(2).
   ,pointsOut(:)      &!  number of (chosen or defined) points in each output grid
   ,pointsOutLand(:)  &!  number of (chosen or defined) land points in each output grid
   ,useCompressGrid(:)    &!  If compressed grid output (e.g. GrADS pdef) output
!                               is selected, this gives the profile number
!                               that generates the compression mapping data file to be used.
!                           Often useCompressGrid(iout)=iout, i.e. a profile can have its own
!                           compression mapping data file, but if there are several identical grids only
!                           a single mapping file is generated, e.g. useCompressGrid(1:3)=1.
   ,varNlev(:)     &!  # of levels for each output variable
!                        It is the product of all "non-horizontal" sizes, e.g. ntiles, or ntiles*nsmax.
!                        Currently this is likely redundant, since we can get nlevs knowing varType.
!                        But it could be useful if later allow subsetting (eg only take levs #2-4).
!                        Main use is in setting up and acessing storage space for output.
   ,varNum(:)      &!  position (in list of available variables) of each  variable selected for output
!                        Used less and less as I move to allocating more space so that each selected
!                        variable carrie smore information about itself only.
   ,varPos(:,:)    &!  position in the list of chosen output variables (e.g. in varNum) of each
!                        variable (For convenience only. This is just 1,2,3,..,nvarOutTot.)
   ,varStartPos(:)  ! position of first element of each variable in output vector (outval)

!-------------------------------------------------------------------------------
! Mappings.
! Note that all mappings use the default/GrADS order: left to right, bottom to top.
! If rows are to be reorder (yrevIn,yrevOut), the point numbers still refer to this GrADS ordering.
!-------------------------------------------------------------------------------
! Note: mapOut, mapOutLand,mapOutCompress are all allocated to have enough space for the
! largest grid required. This can be rather wasteful if one grid is much larger than the others,
! especially if a mapping is not required for a particular profile! In future, probably move
! to a vector with "pointers" to say which part of vector to use for a particular mapping - this
! would also allow easy reuse of the same mapping.

  INTEGER, ALLOCATABLE ::  &!
    mapIn(:)     &!  mapping between input grid and model grid
!                      For each point in the model grid, this is the point number
!                      to use from the input grid (in the order in which the input
!                      grid is presented)
   ,mapInLand(:) &!  mapping between input grid and model land points vector
!                      For each land point in the model grid, this is the point number
!                      to use from the input grid (in the order in which the input grid is presented)
!                      Using a separate map in this way saves land fields from first
!                      having to be mapped from the input grid to the model grid, and
!                      then onto the land points vector.
   ,mapOut(:,:,:)  &!  mapping between the model grid and the output variable for each output grid
!                        Note that the relevant model grid might be e.g. the routing grid.
!                        1st dimension - output profile
!                        2nd dimension - points
!                        3rd dimension - output points
!                        mapOut(:,:,1) is a list of the points to be output,
!                              given as point number in model grid (1:points).
!                        mapOut(:,:,2) is a list of the locations to which the
!                          points (given in mapOut(:,:,1)) are mapped in the
!                          output variable. If the output is not "compressed"
!                          (outCompress=FALSE), these locations are also the locations
!                          in the output grid.
!                        e.g. points=100, mapOut(1,1:2,1)=4,6, outCompress=FALSE
!                            mapOut(1,1:2,2)=10,34 will output the points #4,6
!                            from the model grid of 100 points to points #10,34
!                            in an output grid (extents set up elsewhere).
   ,mapOutCompress(:,:)  &!  mapping between the output variable and locations in the
!                              full, uncompressed output grid. Only used if output is
!                              "compressed" (outCompress=TRUE).
!                        1st dimension - output profile
!                        2nd dimension - points in full output grid
!                       e.g. mapOutCompress=1,10,27 means that the first three points held
!                          in the output vector represent values for points numbered 1, 10 and 27
!                          in the full, uncompressed grid. They are however the first three
!                          values that are actually written.
   ,mapOutLand(:,:,:)   !  mapping between the model land points vector and the output
!                            variable for each output grid
!                            1st dimension - output profile
!                            2nd dimension - land points
!                            3rd dimension - output points
!                            mapOutLand(:,:,1) is a list of the land points to be output,
!                            given as point number in model land points vector (1:land_pts).
!                            mapOutLand(:,:,2) is a list of the locations to which the
!                            points (given in mapOutLand(:,:,1))are mapped in the
!                            output variable. If the output is not "compressed"
!                           (outCompress=FALSE), these locations are also the locations
!                           in the output grid.
!           e.g. land_points=100, mapOutLand(1,1:2,1)=4,6, mapOutLand(1,1:2,2)=10,34
!                will output the points #4,6 from the model land points vector of 100 points
!                to points #10,34 in the output variable.

!###############################################################################
  REAL ::  &!  SCALARS
    undefOut  !  value used for missing data in output (GrADS undef)

  REAL, ALLOCATABLE ::  &!  ARRAYS
    outval(:)         &!  output data
   ,outGridDxy  (:,:) &!  gridbox size of output grid
!                           2nd dimension is 1:2, with values:
!                           1: gridbox size (degrees longitude)
!                           2: gridbox size (degrees latitude)
   ,outGridXY  (:,:)  &!  The location of the point (1,1) of the output grid
!                           2nd dimension is 1:2, with values:
!                           1: longitude of first column of gridpoints in output grid
!                           2: latitude of first row of gridpoints in output grid
   ,outRangeX(:,:)    &!  gives range of x (or lon) values that are to be output
!                           2nd dimension 1=start value  2=end value
   ,outRangeY(:,:)     !  gives range of y (or lat) values that are to be output
!                           2nd dimension 1=start value  2=end value

!###############################################################################
  LOGICAL ::   &!  SCALARS
    dumpWarn     &!  T means an error has been raised regarding a failure to generate the name for a dump
   ,echo         &!  T means that extra messages (and prints of fields) are written to
!                      to standard output throughout the run - can be lots of output.
   ,gradsNc      &!  T means the output will be netCDF that is readable by GrADS
!                    In particular, variables with 2 "levels" dimensions
!                    (e.g. snow layer avriables s(z,ntiles) ) will be
!                    represented as separarate variables
!                    (e.g. ntiles variables each with nz levels).
!                    F means ...all other cases, including netCDF files that
!                    need make no special effort to be readable by GrADS
!                    NB gradsNc must be set to FALSE if not using netCDF files - this is
!                    (perhaps naughtily) relied upon at various points in the code!
   ,outGrads     &!  T means output is potentially readable by GrADS
!                    This is either a GrADS "classic" binary file, or a netCDF file
!                    that has been made consistent with GrADS (gradsNc=TRUE).
!                    One implication is that a GrADS ctl file will be written.
!                    F means no need to ensure readable by GrADS (e.g. ASCII files)
   ,outWarnCtl   &!  T means that an error was raised while attempting to rewrite a GrADS ctl file
   ,outWarnEarly &!  T means that output was generated 'early' (e.g. monthly average
!                       calculated before end of month,
!                       but as near end of month as possible) or not at all (output on
!                       "regular" or short periods
!                       is abandoned) at the end of a section of the run
   ,outWarnUnder &!  T means that a time-average  was not calculated (but was set to
!                      undef value) because
!                      there were insufficient times of data (e.g. monthly average if
!                      run starts mid-month).
   ,numMonth     &!  T means that months in file names are represented as numbers
                  !  F means that months in file names are represented using characters
   ,readFracIC   &!  T means fractional cover (frac) is read as part of initial condition in init_ic
                  !  F means it is read separately by init_frac
   ,redoTdef     &!  T means that a GrADS ctl file will be rewritten if it is found that the number of
!                 !     times written to the data file(s) was not as "expected".
   ,usePseudo    &!  T means a pseudo layer dimension is represented
!                    F there is no pseudo layer dimension; z dimension used instead.
!                    If F, variables with both z and pseudo dimensions
!                       (eg snow layer variables) are represented as a
!                       separate variable for each pseudo level (e.g. tile)
!                       Generally FALSE, except for netCDF files with gradsNC=FALSE.
   ,yrevIn       &!  T means the order of the rows in the input data will be reversed
                  !     Points in the grid are numbered starting from the bottom left corner,
                  !     then going across rows. yrevIn indicates that the input data are not in this order.
                  !     So....F means input data are arranged N to S
                  !           T means input data are arranged S to N
   ,yrevOut      &!  As yrevIn but for output data
   ,zrevOutSnow  &!  T means the order of the snow layers will be reversed in the output data.
!                      F means the layers will be output starting with the surface layer
!                         (furthest from soil surface)
!                      T means the layers will be output starting with the deepest layer
!                         (closest to soil surface)
   ,zrevOutSoil   !   T means the order of the soil layers will be reversed in the output data.
!                      F means the layers will be output starting with the surface layer
!                      T means the layers will be output starting with the deepest layer


  LOGICAL, ALLOCATABLE ::  &!  ARRAYS
!   Each "type" of variables that potentially has >1 level should have a haveXX variable
!   declared below. These are used to determine dimensions required for netCDF files.
    havePFT(:)            &!  indicates if a profile includes PFT variables
   ,haveSCpool(:)         &!  indicates if a profile includes soil pool variables
   ,haveSnow(:)           &!  indicates if a profile includes snow layer variables
   ,haveSoil(:)           &!  indicates if a profile includes soil layer variables
   ,haveTile(:)           &!  indicates if a profile includes tile variables
   ,haveType(:)           &!  indicates if a profile includes surface type variables
   ,outActivePrev(:)      &!  indicates if an output profile was "active" on the
!                              previous timestep. Currently only use is to close
!                              netCDF files as soon as possible.
   ,outAreaLL(:)          &!  indicates if the area of the grid that is to be output is
!                            defined by a range of latitude and longitude.
!                            Only used with pointsFlag(1)=1.
!                            TRUE = range is defiend by latitude and longitude
!                            FALSE = range is defined by row and column numbers
   ,outCompress(:)     &!  T means output is written as a vector in which only model
!                            points are included. This facility was designed with GrADS' pdef
!                            in mind, but it may be possible to adapt for other packages.
!                          F means output is written as a rectangle or vector, with
!                            extra points added (as padding) to give the required grid shape.
   ,outFirstActive(:)  &!  T until first timestep in run when output profile is active
   ,outFirstSection(:) &!  T until first timestep in current "section" of the run at
!                            which profile is "active"
   ,outFirstWrite(:)   &!  T until first output in current "section" of the run has
!                            been written to file
!                            For special cases (e.g. all output to one file), this is
!                            only T until the first output in the FIRST section of the
!                            run for which output was requested.
   ,outLLorder(:)      &!  Flag indicating what coordinates are to be used when determining order
!                            of points in output grid. Only used if pointsFlag(2)=3 (rectangle)
!                            or 1 (select subarea)..
!                            T means latitude and longitude of points are used
!                            F means x,y position in input grid are used
!                            Note that if output is compressed (outCompress=TRUE),
!                            outLLorder does NOT affect the order in which the
!                            data are written, it only affects the ancillary information
!                            that say how the data were mapped (e.g. the GrADS pdef file).
   ,outTemplate(:)     &!  T means GrADS output will use template option for this profile
!                            A single ctl file will be used for all output (not that it's
!                            obvious why anyone would want more than one template ctl
!                            file...eg to cover different subsets of times?!).
   ,rgProfile(:)       &!  T if an output profile includes variables of type RG (routing grid), and these
!                            are not the single point type considered by rpProfile (q.v.)
!                            If TRUE, all variables in the profile are routing variables (types RG and RP).
!                            We could allow non-routing variables to mix with RG-type IF the routing
!                            grid was identical to the main model grid, but that hasn't been coded.
   ,rpProfile(:)       &!  T if an output profile only contains variables of type RP (routing vars at
!                            individual points). A routing variable becomes type RP when extra
!              arguments are given to indicate location on grid and pointsOut=1. Otherwise a
!              routing variable is considered to be type RG and the locations are given via
!              mapOut, as usual.
   ,snapProfile(:)     &!  T if an output profile includes snapshot (instantaneous) values
   ,taccumVar(:)       &!  T if a variable is to be accumulated over time
   ,tmeanProfile(:)    &!  T if an output profile includes time-averaged or time-accumulated values
   ,tmeanVar(:)         !  T if a variable is to be time averaged

!###############################################################################

  CHARACTER(len=7) :: dumpStatus   !  status used when writing dump files (eg new,replace)
  CHARACTER(len=7) :: outStatus   !  status for output files (eg old,new,replace) (not used for dumps)

  CHARACTER(len=formatLen) :: dumpFormat  !  format for dump files
  CHARACTER(len=formatLen) :: outFormat   !  format for output (diagnostic) files
  CHARACTER(len=10) :: runID    !  name of run
  CHARACTER(len=13) :: outEndian   !  for GrADS output, 'little_endian' or 'big_endian'
  CHARACTER(len=150) :: outDir   ! directory used for output

  CHARACTER(len=2), ALLOCATABLE :: &
    varType(:)        &!  the type of each output variable (for each selected variable)
   ,varTypeList(:)     !  the type of each output variable (in list of available vars)

  CHARACTER(len=10), ALLOCATABLE :: &
    outName(:)    !  name for each output profile

  CHARACTER(len=nameLen), ALLOCATABLE :: &
    varNameList(:)   !   name of each output variable (in list of available vars)

  CHARACTER(len=nameLen+sufLen), ALLOCATABLE :: &
    varName(:)       &!  final (unique in profile) name of each selected output variable
   ,varUnitsList(:)   !  list of units of each output variable

  CHARACTER(len=longLen), ALLOCATABLE :: &
    varDesc(:)      &!   description of each selected output variable
   ,varDescList(:)   !   description of each output variable (in list of available vars)

  CHARACTER(len=200), ALLOCATABLE :: &
    openedFileName(:)  &!  The name of an output data file that was opened earlier in
!                            the timestep (when endCall=F in output).
!                            This was introduced to deal with a particular situation - if output
!                            is to be generated
!                            at the start and end of a timestep (which happens when
!                            we want an initial condition),
!                            and this timestep is (eg) the first in a month and
!                            monthly files are to be used, the
!                            code would otherwise try to open the same new file
!                            twice over the timestep.
   ,outCtlFile(:)   &!   name of last opened GrADS ctl file for each output profile
   ,outDataFile(:)  &!   name of last opened data file for each output profile
   ,compressGridFile(:)    !  name of the compression mapping data file (e.g.
!                               GrADS pdef supplementary data file) used for
!                               each profile (if outCompress=TRUE)

!###############################################################################
! Variables describing input files with PFT data.

! Note that these are currently used in INIT_VEG_PFT and then ARE RESET there!!

  INTEGER ::   &!  SCALARS
    npftInFile  !  number of PFTs in each input file
!                      In fact this is the total number of surface types in the file.

  INTEGER, ALLOCATABLE ::  &!  ARRAYS
    pftUse(:)       !  the numbers/positions in the input file of the PFTs to use

!###############################################################################

! Variables that are used during initialisation of output.

! Array variables.

  real, allocatable :: coordList(:,:,:)  !  a list of coordinates of points selected for output

  logical, allocatable :: coord(:)  !  flag indicating if points are selected
!                                        via a list of coordinates
!                                     TRUE: coordinates,   FALSE: index values

  logical, allocatable :: coordLL(:)  !  flag indicating if coordinates of
!                                  selected points are (lat,lon) pairs (TRUE)
!                                  or (x,y) pairs (FALSE)

!-------------------------------------------------------------------------------

  END MODULE inout

