! Could separate into two generic types: readVar and readVarComp.
! eg readVar
!      readVar2d, readVar3d  (2 and 3d separated on # and/or rank of args, eg an nz arg)
!        readVar2dInt,readVar2dReal,readVar3dInt,readVar3dReal
!
! Or could even group together under readVar, then separate as:
! eg readVar
!      readVarFull, readVarComp (comp has an extra argument, probably gratuitous! - AND SO A BAD IDEA!)
!        readVarFull2d, readVarFull3d, readVarComp2d,readVarComp3d
!          readVarFull2dInt, ...etc...
!
! Note that ASCII-reading routines are very similar and largely consist of file positioning, which is common
! between int and real versions - might well be better to do all positioning via a single routine that only diverges at last possible moment to read int or real (certainly easier to maintain!).
! This would suggest that the organisation of routines might be something like readvar > readvarASCII > readvarASCIIReal. Maybe?
!
! Why do we need 2D and 3D versions? Can't we make 2D a subset of 3D? Probably difficulty is to be conformable arrays, but maybe could go part way?
! Possibly if 3D version could read a subset of the 3D field, e.g. iz=1,3?? Trouble is likely still array dims!

!----------------------------------------------------------------------
! module readWrite_Mod
!
! Contains routines to read and writes variable to/from a file.
! Specific routines are available for real or integer variables, and for different file formats.
!
! Note: irecPrev must not be changed (i.e. should be intent(in)) when reading
! from netCDF files. See note in FILE_UTILS.
!
!###############################################################################
!###############################################################################

  MODULE readWrite_mod

  USE inout, ONLY :  &
!  imported scalar parameters
    formatAsc,formatBin,formatNc,formatPP,stdIn,jinUnit

  USE grid_utils, ONLY :  &
!  imported procedures
     getXYpos,mapAtoB

  USE jules_netcdf, ONLY :  &
!  imported procedures
     readVar2dReal_ncVector_gswp2,readVar2dReal_nc_pilps2e  &
    ,readVar2dReal_nc_princet,readVar2dReal_nc_tseries  &
!  imported scalar parameters
    ,ncTypeGswp2,ncTypePilps2e,ncTypePrincet,ncTypeTseries,ncTypeWatch

  USE netcdf, ONLY :  &
!  imported scalar parameters
     nf90_noerr  &
!  imported procedures
    ,nf90_put_var

  USE rwErr_mod, ONLY :  &
!  imported procedures
     rwErr
!-------------------------------------------------------------------------------

  IMPLICIT NONE

! Interfaces to switch between integer and real variables.

  INTERFACE readVar2d
  MODULE PROCEDURE readVar2dInt,readVar2dReal
  END INTERFACE

  INTERFACE writeVar
  MODULE PROCEDURE writeVarInt,writeVarReal
  END INTERFACE

! Scalar parameters.
  INTEGER, PARAMETER :: lineLen = 10 !  number of values per line in ASCII output

  CONTAINS

!#######################################################################
!#######################################################################
! subroutine readVar2dInt
! Internal procedure in readwrite_mod..
! Module procedure with generic name readVar2d.
! Reads a 2D integer variable from file.

  SUBROUTINE readVar2dInt( t,z,field,fieldCode,irecPrev,nfieldFile,nheaderFile  &
                          ,nheaderT,nheaderField,nlineField,nxIn,nyIn  &
                          ,unit,varName,mapIn,mapOut  &
                          ,fileFormat,outval,errMess,ncCallType,ncTypeArg )

  IMPLICIT NONE

! NB outval has to be declared before any other declaration refers to it.
  INTEGER, INTENT(out) ::   &!  out ARRAYS
    outval(:,:)      !  the output array

  INTEGER, INTENT(in) ::  &!  in SCALARS
    field          &!  the number (index or location) of the field to be read
   ,fieldCode      &!  the field code (e.g. STASH) for the field (for PP files only)
   ,nfieldFile     &!  the number of fields per time in the file
   ,nheaderFile   &!  number of headers records at start of file
   ,nheaderT   &!  number of header records at start of each time
   ,nheaderField  &!  number of header records before each field (xy plane)
   ,nlineField     &!  the number of records per field (excl headers)
   ,nxIn       &!  x extent of input data grid
   ,nyIn       &!  y extent of input data grid
   ,t          &!  "time" level to be read from file
   ,unit       &!  unit to read from
   ,z           !  z level to be read from file. Only used for netCDF files.

  INTEGER, INTENT(inout) ::  &!  in SCALARS
    irecPrev    !  record number of last field read from this file

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    mapIn(SIZE(outval))    &!  points to select from input
   ,mapOut(SIZE(outval))    !  points to fill in output

  INTEGER ::  &!  local SCALARS
    irecOne   &!  record number of first record of field
   ,nx        &!  x extent of output array
   ,ny         !  y extent of output array

  INTEGER ::  &! local ARRAYS
    tmpval(nxIn*nyIn)      &!  all input data for one level
   ,tmpval2(SIZE(outval))   !  input data at the points selected for output

  CHARACTER(len=*), INTENT(in) ::   &!  in SCALARS
    errMess     &!  message printed on error
   ,fileFormat  &!  format of file
   ,ncCallType  &!  flag indicating what type of data to be read - a hack for netCDF
   ,ncTypeArg   &!  indicates "type" (format) of netCDF files
   ,varName      !  name of variable

!--------------------------------------------------------------------------------

! Get extents of output array.
  nx = SIZE( outval,1 )
  ny = SIZE( outval,2 )

! Check for (some) incorrect arguments.


! Note that we can have nx>nxIn or ny>nyIn, if a rectangular input grid has
! been converted to a vector model grid (e.g. by selecting land points only).

! Check that field number is within range.
  SELECT CASE ( fileFormat )
    CASE ( formatAsc,formatBin )
      IF ( field > nfieldFile ) THEN
        WRITE(*,*)'ERROR: readVar2dInt: field > nfieldFile.'
        WRITE(*,*)'Requested field number exceeds number said to be in the file.'
!       Call error routine, indicating a stop.
        CALL rwErr( fileFormat,.TRUE.  &
                   ,field=field,nfieldFile=nfieldFile,unit=unit,errMess1=errMess )
      ENDIF
  END SELECT

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )
!     ASCII input
!     Read the entire input field and then select points.
!     Note that if we know nlineField, we could now just read the selected points, as done for GrADS...
!     irecOne is the field number (accumulated from start of file, excl headers)
      irecOne = (t-1)*nfieldFile + field
      CALL readIntAscii( irecPrev,irecOne,nfieldFile  &
                        ,nheaderFile,nheaderT,nheaderField,nlineField  &
                        ,unit,tmpval,errMess )
!     Get the desired points.
      tmpval2(:) = mapAtoB( nxIn*nyIn,nx*ny,nx*ny,mapIn,mapOut,0,tmpval)
!     Convert vector to 2-D array.
      outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )

!-------------------------------------------------------------------------------
    CASE ( formatBin )
!     GrADS input
      WRITE(*,*)'SORRY! Code does not yet exist to read integer from GrADS!'
      WRITE(*,*)'Stopping in readVar2dInt'
      STOP
!-------------------------------------------------------------------------------

    CASE ( formatNc )
!     netCDF input
!     irecOne is the "time" level of required field.
!     In future may want to pass a "z" index too.
      irecOne = t
      WRITE(*,*)'SORRY! netCDF code does not yet exist for reading integer!'
      WRITE(*,*)'Stopping in readVar2dInt'
      STOP
!-------------------------------------------------------------------------------

    CASE ( formatPP )
!     PP file input
      WRITE(*,*)'SORRY! Code does not yet exist to read integer from PP file!'
      WRITE(*,*)'Stopping in readVar2dInt'
      STOP
!-------------------------------------------------------------------------------

    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in readVar2dInt'
      STOP

  END SELECT

  END SUBROUTINE readVar2dInt
!#######################################################################
!#######################################################################
! subroutine readVar2dReal
! Internal procedure in readwrite_mod.
! Module procedure with generic name readVar2d.
! Reads a 2D real variable from file.

  SUBROUTINE readVar2dReal( t,z,field,fieldCode,irecPrev,nfieldFile,nheaderFile  &
                           ,nheaderT,nheaderField,nlineField,nxIn,nyIn  &
                           ,unit,varName,mapIn,mapOut,fileFormat  &
                           ,outval,byteSwap,errMess,ncCallType,ncTypeArg )

  IMPLICIT NONE

! NB outval has to be declared before any other declaration refers to it.
  REAL, INTENT(out) ::   &!  out ARRAYS
    outval(:,:)      !  the output array

  INTEGER, INTENT(in) ::  &!  in SCALARS
    field          &!  the number (index or location) of the field to be read
   ,fieldCode      &!  the field code (e.g. STASH) for the field (for PP files only)
   ,nfieldFile     &!  the number of fields per time in the file
   ,nheaderField   &!  number of header records before each field (xy plane)
   ,nheaderFile    &!  number of headers records at start of file
   ,nheaderT       &!  number of header records at start of each time
   ,nlineField     &!  number of records per field (excl headers)
   ,nxIn           &!  x extent of input data grid
   ,nyIn           &!  y extent of input data grid
   ,t              &!  "time" level to be read from file
   ,unit           &!  unit to read from
   ,z               !  z level to be read from file. Only used for netCDF files.

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    mapIn(SIZE(outval))    &!  points to select from input
   ,mapOut(SIZE(outval))    !  points to fill in output

  INTEGER, INTENT(inout) ::  &!  in SCALARS
    irecPrev    !  record number of last field read from this file

  INTEGER ::  &!  local SCALARS
    i,irec    &!  work
   ,irecOne   &!  record number of first record of field
   ,ix,ixOut,iy,iyOut  &!  work
   ,np        &!  number of points for which values are required (=size(outval))
   ,npRange   &!  number of points in a section of data
   ,nx        &!  x extent of output array
   ,ny        &!  y extent of output array
   ,p1,p2      !  point numbers

  REAL ::  &! local ARRAYS
    tmpval(nxIn*nyIn)      &!  all input data for one level
   ,tmpval2(SIZE(outval))   !  input data at the points selected for output

  LOGICAL, INTENT(in) :: byteSwap  !  flag indicating if byte order ("endianness")
!                                       of data is to be reversed

  CHARACTER(len=*), INTENT(in) ::   &!  in SCALARS
    errMess    &!  message printed on error
   ,fileFormat &!  format of file
   ,ncCallType &! flag indicating what type of data to be read - a hack for netCDF
   ,ncTypeArg  &!  indicates "type" (format) of netCDF files
   ,varName     !  name of variable

!--------------------------------------------------------------------------------

! Get extents of output array.
  np = SIZE( outval )
  nx = SIZE( outval,1 )
  ny = SIZE( outval,2 )

! Check for (some) incorrect arguments.

! Note that we can have nx>nxIn if a rectangular input grid has
! been converted to a vector model grid (e.g. by selecting land points only).

! Check that field number is within range.
  SELECT CASE ( fileFormat )
    CASE ( formatAsc,formatBin,formatPP )
      IF ( field > nfieldFile ) THEN
        WRITE(*,*)'ERROR: readVar2dReal: field > nfieldFile.'
        WRITE(*,*)'Requested field number exceeds number said to be in the file.'
!       Call error routine, indicating a stop.
        CALL rwErr( fileFormat,.TRUE.  &
                   ,field=field,nfieldFile=nfieldFile,unit=unit,errMess1=errMess )
      ENDIF
  END SELECT

!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )

!     ASCII input
!     Read the entire input field and then select points.
!     Note that if we know nlineField, we could now just read the selected points,
!     as done for GrADS...
!     irecOne is the field number (accumulated from start of file, excl headers)
      irecOne = (t-1)*nfieldFile  + field
      CALL readRealAscii( irecPrev,irecOne,nfieldFile  &
                         ,nheaderFile,nheaderT,nheaderField,nlineField  &
                         ,unit,tmpval,errMess )
!     Extract the desired points.
      tmpval2(:) = mapAtoB( nxIn*nyIn,nx*ny,nx*ny,mapIn,mapOut,0.0,tmpval)
!     Convert vector to 2-D array.
      outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )
!     Update record counter
      irecPrev = irecOne

!-----------------------------------
    CASE ( formatBin )
!     GrADS input (direct access)
!     (Note this could be done similarly to how ASCII are read above - reading all points
!     into tmpval2 and then reshaping. But doing this way for now.)
!     Get record number of first point in field (after any headers).
      irecOne = nheaderFile +  &! file headers
            (t-1)*nfieldFile*(nxIn*nyIn+nheaderField) + (t-1)*nheaderT +   &! earlier times
            (field-1)*(nxIn*nyIn+nheaderField) +   &!  earlier fields
            nheaderField + 1
!     Read the selected points.
      IF ( byteSwap ) THEN
        DO i=1,nx*ny
          irec = irecOne + mapIn(i) - 1
!         Get x,y location in output grid.
          CALL getXYpos( mapOut(i),nx,ny,ixOut,iyOut )
          CALL readRealGrads( irec,unit,outval(ixOut:ixOut,iyout),errMess )
!         Reverse byte order.
          tmpval(1) = outval(ixOut,iyout)
          CALL native_4byte_real( tmpval(1),outval(ixOut,iyout) )
        ENDDO
      ELSE
        DO i=1,nx*ny
          irec = irecOne + mapIn(i) - 1
!         Get x,y location in output grid.
          CALL getXYpos( mapOut(i),nx,ny,ixOut,iyOut )
          CALL readRealGrads( irec,unit,outval(ixOut:ixOut,iyout),errMess )
        ENDDO
      ENDIF  !  byteSwap

!-----------------------------------
    CASE ( formatNc )
!     netCDF input
!     At present we deal with netCDF files from different sources (and therefore
!     with different structures) by using different but similar pieces of code!

       SELECT CASE ( ncTypeArg )

         CASE ( ncTypeGSWP2,ncTypeWatch )
!          GSWP2- and WATCH-specific code.

!           Establish the range of input points requested.
            p1 = MINVAL( mapIn(:) )
            p2 = MAXVAL( mapIn(:) )
            npRange = p2 - p1 + 1

!           If the number of points needed to cover the range of input points (npRange)
!           is comparable to the number of required points, read the full range of
!           points and then select the required points. Otherwise, read each point individually.
!           The thinking behind this approach is that we don't want to read any more
!           data than we have to, but there is some overhead (unknown because not
!           investigated) to the point-by-point approach, as currently coded (i.e. several
!           calls to netCDF library functions).

            IF ( npRange > 3*np ) THEN
!             Read only the selected points, point-by-point.
              DO i=1,nx*ny
                CALL readVar2dReal_ncVector_gswp2( ncCallType,ncTypeArg  &
                                    ,unit,varName,t,1,1,tmpval(1:1),mapIn(i) )
!               Get x,y location in output grid.
                CALL getXYpos( mapOut(i),nx,ny,ixOut,iyOut )
                outval(ixout,iyout) = tmpval(1)
              ENDDO
            ELSE
!             Read a range of points, then select the required points.
!             Note that this option is also used for a single point (nx*ny=1).
              CALL readVar2dReal_ncVector_gswp2( ncCallType,ncTypeArg,unit  &
                                         ,varName,t,1,npRange,tmpval(1:npRange),p1 )
!             Extract the desired points. Use an offset with mapIn to account for the
!             fact that mapIn is defined for the
!             full input grid, whereas we only have a subsection. Likely needs to be
!             revisited when we have more general
!             (non-vector) netCDF code here!
              tmpval2(:) = mapAtoB( npRange,nx*ny,nx*ny,mapIn-p1+1,mapOut,0.0,tmpval)
!             Convert vector to 2-D array.
              outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )
            ENDIF

!--------------------------------------------------
         CASE ( ncTypePilps2e )
!          Code for PILPS2e.

!          If the number of points requested is < a small number, read only the required points,
!          otherwise read the full grid and extract the required points. This probably doesn't
!          make much difference to time taken for PILPS because the grid is relatively small
!          (timings have not been investigated)...but I'm following the existing GSWP2 code
!          where the larger grid made i/o savings more substantial.

           IF ( np < 6 ) THEN
!            Read only the selected points, point-by-point.
             DO i=1,nx*ny
!              Get x,y location in input grid.
               CALL getXYpos( mapIn(i),nxIn,nyIn,ix,iy )
!              Read this datum.
               CALL readVar2dReal_nc_pilps2e( ncCallType  &
                                   ,unit,varName,1,ix,1,iy,t,tmpval(1:1) )
!              Get x,y location in output grid.
               CALL getXYpos( mapOut(i),nx,ny,ixOut,iyOut )
               outval(ixout,iyout) = tmpval(1)
             ENDDO
           ELSE
!            Read the full grid, then select the required points.
             CALL readVar2dReal_nc_pilps2e( ncCallType,unit,varName  &
                                               ,nxIn,1,nyIn,1,t,tmpval )
!            Extract the desired points.
             tmpval2(:) = mapAtoB( nxIn*nyIn,nx*ny,nx*ny,mapIn,mapOut,0.0,tmpval)
!            Convert vector to 2-D array.
             outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )
           ENDIF
!--------------------------------------------------
         CASE ( ncTypePrincet )
!          Princeton-specific code.

!          Following PILPS2e example.
!          If the number of points requested is < a small number, read only the required points,
!          otherwise read the full grid and extract the required points.

           IF ( np < 6 ) THEN
!            Read only the selected points, point-by-point.
             DO i=1,nx*ny
!              Get x,y location in input grid.
               CALL getXYpos( mapIn(i),nxIn,nyIn,ix,iy )
!              Read this datum.
               CALL readVar2dReal_nc_princet( ncCallType  &
                                   ,unit,varName,1,ix,1,iy,t,tmpval(1:1) )
!              Get x,y location in output grid.
               CALL getXYpos( mapOut(i),nx,ny,ixOut,iyOut )
               outval(ixout,iyout) = tmpval(1)
             ENDDO
           ELSE
!            Read the full grid, then select the required points.
             CALL readVar2dReal_nc_princet( ncCallType,unit,varName  &
                                               ,nxIn,1,nyIn,1,t,tmpval )
!            Extract the desired points.
             tmpval2(:) = mapAtoB( nxIn*nyIn,nx*ny,nx*ny,mapIn,mapOut,0.0,tmpval)
!            Convert vector to 2-D array.
             outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )
           ENDIF

!--------------------------------------------------
         CASE ( ncTypeTseries )
!          A "simple" time series at a point.

!          Read this datum.
           CALL readVar2dReal_nc_tseries( ncCallType,unit,varName,t,tmpval(1:1) )
!          Load into output variable.
           outval(1,1) = tmpval(1)

!--------------------------------------------------
         CASE default
           WRITE(*,*)'ERROR: readVar2dReal: No code for ncTypeArg=',TRIM(ncTypeArg)
           STOP

         END SELECT
!-------------------------------------------------------------------------------
    CASE ( formatPP )

!     PP format input
!     Read the entire input field and then select points.
!     irecOne is the field number (accumulated from start of file, excl headers)
      irecOne = (t-1)*nfieldFile  + field
      CALL readRealPP( fieldCode,nfieldFile,unit,errMess,irecOne,irecPrev,tmpval )

!     Extract the desired points.
      tmpval2(:) = mapAtoB( nxIn*nyIn,nx*ny,nx*ny,mapIn,mapOut,0.0,tmpval)
!     Convert vector to 2-D array.
      outval(:,:) = RESHAPE( tmpval2(:), (/ nx,ny /) )
!     Update record counter
      irecPrev = irecOne
!-------------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in readVar2dReal'
      STOP

  END SELECT

  END SUBROUTINE readVar2dReal
!#######################################################################
!#######################################################################
! subroutine readVar2dComp
! Internal procedure in readWrite_mod.
! Read a 2-D variable from file. Points in the horizontal are selected via a mapping.
! The output is a vector in which the x and y directions are compressed (e.g. full grid is
! compressed to land points), e.g. input A(nx,ny), output B(nn), where nn<=nx*ny and the
! points in B are selected by a mapping.

  SUBROUTINE readVar2dComp( t,z,field,fieldCode,irecPrev,nfieldFile,nheaderFile,nheaderT  &
                           ,nheaderField,nlineField,nxIn,nyIn  &
                           ,unit,varName,mapIn,mapOut  &
                           ,fileFormat,outval,errMess,ncCallType,ncTypeArg )

  IMPLICIT NONE

! NB outval has to be declared before any other declaration refers to it.
  REAL, INTENT(out) ::   &!  out ARRAYS
    outval(:)      !  the output vector

  INTEGER, INTENT(in) ::  &!  in SCALARS
    field         &!  the number (index or location) of the field to be read
   ,fieldCode     &!  the field code (e.g. STASH) for the field (for PP files only)
   ,nfieldFile    &!  the number of fields per time in the file
   ,nheaderFile   &!  number of headers records at start of file
   ,nheaderT      &!  number of header records at start of each time
   ,nheaderField  &!  number of header records before each field (xy plane)
   ,nlineField    &!  number of records per field (excl headers)
   ,nxIn          &!  x extent of input data grid
   ,nyIn          &!  y extent of input data grid
   ,t             &!  "time" level to be read from file
   ,unit          &!  unit to read from
   ,z              !  z level to be read from file. Only used for netCDF files.

  INTEGER, INTENT(inout) ::  &!  in SCALARS
    irecPrev    !  record number of last field read from this file

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    mapIn(SIZE(outval))    &!  points to select from input
   ,mapOut(SIZE(outval))    !  points to fill in output

  INTEGER ::  &!  local SCALARS
    i,irec    &!  work
   ,irecOne   &!  record number of first record of field
   ,j         &!  work
   ,npoints    !  number of points in output array

  REAL ::  &! local ARRAYS
    tmpval(nxIn*nyIn)   !  all input data for one level

  CHARACTER(len=*), INTENT(in) ::   &!  in SCALARS
    errMess     &!  message printed on error
   ,fileFormat  &!  format of file
   ,ncCallType  &!  flag indicating what type of data to be read - a hack for netCDF
   ,ncTypeArg   &!  indicates "type" (format) of netCDF files
   ,varName      !  name of variable

!--------------------------------------------------------------------------------

! Get number of points to be found.
  npoints = SIZE( outval,1 )

! Check for (some) incorrect arguments.
! Check that field number is within range.
  SELECT CASE ( fileFormat )
    CASE ( formatAsc,formatBin,formatPP )
      IF ( field > nfieldFile ) THEN
        WRITE(*,*)'ERROR: readVar2dComp: field > nfieldFile.'
        WRITE(*,*)'Requested field number exceeds number said to be in the file.'
!       Call error routine, indicating a stop.
        CALL rwErr( fileFormat,.TRUE.  &
                   ,field=field,nfieldFile=nfieldFile,unit=unit,errMess1=errMess )
      ENDIF
  END SELECT

!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )

    CASE ( formatAsc )
!     ASCII input
!     Read the entire input field and then select points.
!     Note that if we know nlineField, we could now just read the selected points,
!      as done for GrADS...
!     irecOne is the field number (accumulated from start of file, excl headers)
      irecOne = (t-1)*nfieldFile  + field
      CALL readRealAscii( irecPrev,irecOne,nfieldFile  &
                         ,nheaderFile,nheaderT,nheaderField,nlineField  &
                         ,unit,tmpval,errMess )

!     Get the desired points.
      outval(:) = mapAtoB( nxIn*nyIn,npoints,npoints,mapIn,mapOut,0.0,tmpval)

!-------------------------------------------------------------------------------
    CASE ( formatBin )
!     GrADS input (direct access)
!     Read the selected points.
!     Get record number of first point in field (after any headers).
      irecOne = nheaderFile +  &!  file headers
            (t-1)*nfieldFile*(nxIn*nyIn+nheaderField) + (t-1)*nheaderT +  &!  earlier times
            (field-1)*(nxIn*nyIn+nheaderField) +  &!  earlier fields
            nheaderField + 1
      DO i=1,npoints
        irec = irecOne + mapIn(i) - 1
        j = mapOut(i)
        CALL readRealGrads( irec,unit,outval(j:j),errMess )
      ENDDO

!-------------------------------------------------------------------------------
    CASE ( formatNc )
!     netCDF input

!     At present we deal with netCDF files from different sources (and therefore
!     with different structures) by using different but similar pieces of code!

       SELECT CASE ( ncTypeArg )

         CASE ( ncTypeGswp2,ncTypeWatch )
!         GSWP2- and WATCH-specific code.
!         Read the entire input field and then select points.
          CALL readVar2dReal_ncVector_gswp2( ncCallType,ncTypeArg,unit  &
                                            ,varName,t,z,nxIn,tmpval )
!         Get the desired points.
          outval(:) = mapAtoB( nxIn*nyIn,npoints,npoints,mapIn,mapOut,0.0,tmpval)

         CASE default
           WRITE(*,*)'ERROR: readVar2dComp: No code for ncTypeArg=',TRIM(ncTypeArg)
           STOP

         END SELECT
!-------------------------------------------------------------------------------

    CASE ( formatPP )
!     PP format input
!     Read the entire input field and then select points.
!     irecOne is the field number (accumulated from start of file, excl headers)
      irecOne = (t-1)*nfieldFile  + field
      CALL readRealPP( fieldCode,nfieldFile,unit,errMess,irecOne,irecPrev,tmpval )

!     Get the desired points.
      outval(:) = mapAtoB( nxIn*nyIn,npoints,npoints,mapIn,mapOut,0.0,tmpval)

!-------------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in readVar2dComp'
      STOP

  END SELECT

  END SUBROUTINE readVar2dComp
!###############################################################################
!###############################################################################
! subroutine readVar3dComp
! Internal procedure in readWrite_mod.
! Read a 3-D variable (eg x,y,z or x,y,type) from file. Points in the horizontal
! are selected via a mapping.
! In the output, x and y directions are compressed (e.g. full grid is compressed
! to land points).
! e.g. input A(nx,ny,nz), output B(nn,nz), where nn<=nx*ny and the x,y points
! in B are selected by a mapping.

  SUBROUTINE readVar3dComp( t,field,fieldCode,noChange,irecPrev,nfieldFile  &
                           ,nheaderFile,nheaderT  &
                           ,nheaderField,nlineField,nxIn,nyIn  &
                           ,nz,zrev,unit,varName,mapIn,mapOut,fileFormat,outval  &
                           ,errMess,ncCallType,ncTypeArg )

  IMPLICIT NONE

! NB outval has to be declared before any other declaration refers to it.
  REAL, INTENT(out) ::   &!  out ARRAYS
    outval(:,:)      !  the output array

  INTEGER, INTENT(in) ::  &!  in SCALARS
    field         &!  the number (index or location) of the first field to be read
!                       The first field is the first level of the variable.
   ,fieldCode     &!  the field code (e.g. STASH) for the field (for PP files only)
   ,nfieldFile    &!  the number of fields per time in the file
   ,nheaderFile   &!  number of headers records at start of file
   ,nheaderT      &!  number of header records at start of each time
   ,nheaderField  &!  number of header records before each field (xy plane)
   ,nlineField    &!  number of records per field (exl headers)
   ,nxIn          &!  x extent of input data grid
   ,nyIn          &!  y extent of input data grid
   ,nz            &!  z extent of input and output data grids
   ,t             &!  "time" level to be read from file
   ,unit           !  unit to read from

  INTEGER, INTENT(inout) ::  &!  in SCALARS
    irecPrev    !  record number of last field read from this file

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    mapIn(SIZE(outval,1))    &!  points to select from input
   ,mapOut(SIZE(outval,1))    !  points to fill in output

  INTEGER ::  &!  local SCALARS
    f,i,irec  &!  work
   ,irec1,irec2   &!  work: record numbers
   ,iz,izz,j      &!  work
   ,npoints        !  number of points in a single vertical level of output array

!  real ::  &! local ARRAYS
!    tmpval(nxIn*nyIn*nz)   &!  all input data
!   ,tmpval2(size(outval))   !  input data at the points selected for output

  LOGICAL, INTENT(in) ::  &!  in SCALARS
    noChange   &!  T means do not change the field number between reads
                !    This is used with ASCII files when field=0, which readReal
                !    interprets as 'read next field'. In this case, don't want
                !    the field number to change in loop over levels - always want zero.
   ,zrev    !  T means that the order of the levels is to be reversed before output

  CHARACTER(len=*), INTENT(in) ::   &!  in SCALARS
    errMess     &!  message printed on error
   ,fileFormat  &!  format of file
   ,ncCallType  &!  flag indicating what type of data to be read - a hack for netCDF
   ,ncTypeArg   &!  indicates "type" (format) of netCDF files
   ,varName      !  name of variable

!--------------------------------------------------------------------------------

! Get number of points to be found.
  npoints = SIZE( outval,1 )

!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )

    CASE ( formatAsc )
!     ASCII input

!     Get the field number of the field immediately before the first required
!     field (first level)(accumulated from start of file, excl headers).
!     An increment of 1 will then be the first required field.
      irec = (t-1)*nfieldFile  + field - 1

!     Loop over z levels. readVar2Dcomp selects the required points via mapping.
      DO iz=1,nz
!       Get field number.
        irec = irec + 1
!       Override if necessary.
        f = irec
        IF ( noChange ) f=field
!       Reverse order of levels, if necessary.
        izz = iz
        IF ( zrev ) izz = nz - iz + 1
        CALL readVar2dComp( t,iz,f,fieldCode,irecPrev,nfieldFile  &
                ,nheaderFile,nheaderT,nheaderField,nlineField,nxIn,nyIn  &
                ,unit,varName,mapIn,mapOut,fileFormat  &
                ,outval(:,izz),errMess,ncCallType,ncTypeArg)
      ENDDO

!-------------------------------------------------------------------------------
    CASE ( formatBin )

!     GrADS input (direct access)
!     (Note this could be done similarly to how ASCII are read above - reading all points
!     into tmpval2 and then reshaping. If we were reading entire field, we would use
!     repeated calls to readVar2dComp - as done for ASCII - but as it is, we can go direct
!     to readRealGrads.)
!

!     Get record number just before start of first field (first level of variable)
!     (before field headers).
      irec1 = nheaderFile +   &!  file headers
          (t-1)*nfieldFile*(nxIn*nyIn+nheaderField) + (t-1)*nheaderT +   &!  earlier times
          (field-1)*(nxIn*nyIn+nheaderField) !  earlier fields

      DO iz=1,nz
!       Get record number just before first point in this field (level) (at end
!       of field headers).
        irec2 = irec1 + (iz-1)*(nxIn*nyIn+nheaderField) + nheaderField
!       Reverse order of levels, if necessary.
        izz = iz
        IF ( zrev ) izz = nz - iz + 1
!       Read the selected points for this level.
        DO i=1,npoints
          irec = irec2 + mapIn(i)
!         Calculate what point to use in output array.
          j = mapOut(i)
          CALL readRealGrads( irec,unit,outval(j:j,izz),errMess )
        ENDDO
      ENDDO

!-------------------------------------------------------------------------------
    CASE ( formatNc )

!     netCDF input
!     Loop over z levels. readVar2Dcomp selects the required points via mapping.
      DO iz=1,nz
!       Reverse order of levels, if necessary.
        izz = iz
        IF ( zrev ) izz = nz - iz + 1
        CALL readVar2dComp( t,iz,iz,fieldCode,irecPrev,nfieldFile  &
                    ,nheaderFile,nheaderT,nheaderField,nlineField,nxIn,nyIn  &
                    ,unit,varName,mapIn,mapOut,fileFormat  &
                    ,outval(:,izz),errMess,ncCallType,ncTypeArg)
      ENDDO
!-------------------------------------------------------------------------------
    CASE( formatPP )
!     PP-format input.

!     Get the field number of the field immediately before the first required
!     field (first level)(accumulated from start of file).
!     An increment of 1 will then be the first required field.
      irec = (t-1)*nfieldFile  + field - 1

!     Loop over z levels. readVar2Dcomp selects the required points via mapping.
      DO iz=1,nz
!       Get field number.
        irec = irec + 1
!       Reverse order of levels, if necessary.
        izz = iz
        IF ( zrev ) izz = nz - iz + 1
        CALL readVar2dComp( t,iz,irec,fieldCode,irecPrev,nfieldFile,nheaderFile,nheaderT  &
                           ,nheaderField,nlineField,nxIn,nyIn,unit,varName  &
                           ,mapIn,mapOut,fileFormat  &
                           ,outval(:,izz),errMess,ncCallType,ncTypeArg)
      ENDDO

!-------------------------------------------------------------------------------
    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in readVar3dComp'
      STOP

  END SELECT

  END SUBROUTINE readVar3dComp
!#######################################################################
!#######################################################################
! subroutine readIntASCII
! Internal procedure in readWrite_mod.
! Read an integer variable from an ASCII file.
! Repositions file as necessary.
! It is assumed that any previous reads of this file read a whole field, so
! file is either positioned at start of file or at the end of a field.
! Will fail if previously read only part of a field that extends over
! more than one line.
! This routine could be much shorter if time headers were not allowed - then
! we could just loop over fields.

  SUBROUTINE readIntASCII( fPrevAcc,fNextAcc,nfieldFile,nheaderFile  &
                          ,nheaderT,nheaderField,nlineField  &
                          ,unit,val,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    nfieldFile    &!  the number of fields (per time) in the file
   ,nheaderFile   &!  the number of header lines at the start of the file
   ,nheaderT      &!  the number of header lines at the start of each time
   ,nheaderField  &!  the number of header lines at the start of each field
   ,nlineField    &!  the number of lines per field (excl headers) - this is needed
                   !    if need to backspace within a file
                   !  <1 means don't attempt to backspace, always move forward through file
   ,unit           !  unit to read from

  INTEGER, INTENT(inout) ::  &!  inout SCALARS
    fNextAcc   &!  the field number to be read. This is the number
!                    accumulated from start of file.
!                    0 means next field is read.
   ,fPrevAcc    !  the last field read from this file. This is the number
!                    accumulated from start of file.

  INTEGER, INTENT(out) ::   &!  out ARRAYS
    val(:)      !  the data that are read

  INTEGER ::  &!  local SCALARS
    f        &!  loop counter
   ,fNext    &!  the field (at t=tNext, i.e. not accumulated) to be read next
!                  (fields do not include headers)
   ,fPrev    &!  the field (at t=tPrev, i.e. not accumulated) that was last read
!                  (fields do not include headers)
   ,ierr     &!  error value
   ,inc      &!  the number of fields between fPrevAcc and fNextAcc
   ,jt       &!  work
   ,l        &!  loop counter
   ,nline    &!  work: number of lines to read or backspace
   ,nlineT    &!  the number of lines per time in file (incl headers)
   ,nlineFieldTot   &!  the number of lines per field per time (incl headers)
   ,t        &!  loop counter
   ,tNext    &!  the "time" level in the file of fNextAcc
   ,tPrev    &!  the "time" level in the file of fPrevAcc
   ,xfPrev    !  work: version of fPrev

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    errMess    !  message printed on error

!-------------------------------------------------------------------------------

! For now, I don't expect nlineField>0 (although the code "should" still work).
  IF ( nlineField > 0 ) THEN
    WRITE(*,*)'nlineField=',nlineField,' unit=',unit
    WRITE(*,*)'readIntASCII: nlineField>0: Unexpected at present!'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF

  IF ( nfieldFile < 1 ) THEN
!   Note that we need nfieldFile, even if we are attempting to read the
!   next record, since there is always a division with nfieldFile in denominator.
!   In future we might like to have a case (e.g. fNextAcc<0) which indicates
!   "REALLY read next record, not trying to avoid headers etc.".
    WRITE(*,*)'ERROR: readIntASCII: nfieldFile < 1'
    WRITE(*,*)'nfieldFile must be >= 0'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF

! Can't rewind standard input, so assume we will read next field.
! Since the file on jinUnit must behave like stdIn (in most instances it will be!),
! we also don't rewind if we are on jinUnit
! Ignore possibility of skipping any fields.
  IF ( unit==stdIn .OR. unit == jinUnit ) THEN
    fPrevAcc = 0
    fNextAcc = 0
  ENDIF

! fNextAcc=0 means read next field (still deals with headers).
  IF ( fNextAcc == 0 ) fNextAcc=fPrevAcc+1

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  IF ( nlineField >= 1 ) THEN

!   Knowing the number of lines of data per field, we can move forward or back
!   through the file.

    IF ( fNextAcc<fPrevAcc/2 .AND. fNextAcc>100 ) THEN
!   If next field is nearer start of file that current position, rewind and go
!   from start of file.
      if (unit/=77) REWIND( unit )
      fPrevAcc = 0
    ENDIF

!   Get the number of intervening fields.
    inc = fNextAcc - fPrevAcc - 1


!   Calculate number of lines per field, including headers.
    nlineFieldTot = nheaderField + nlineField
!   Calculate number of lines per time, including headers.
    nlineT = nheaderT + nfieldFile*nlineFieldTot

!   Work out where the last field read was in the file.
    tPrev = CEILING( REAL(fPrevAcc) / REAL(nfieldFile) )
    jt = MAX( 1, tPrev )   !   if fPrevAcc=0, want tPrev=0
    fPrev = fPrevAcc - (jt-1)*nfieldFile

!   Work out where the next field to read is in the file..
    tNext = CEILING( REAL(fNextAcc) / REAL(nfieldFile) )
    fNext = fNextAcc - (tNext-1)*nfieldFile

!   Calculate number of lines to be read or backspaced so as to reposition the
!   file as necessary.
    nline = 0
    IF ( inc < 0 ) THEN

!     Move back through the file.

!     Lines back to start of current time.
      IF ( tNext-tPrev < 0 ) nline = nline + fPrev*nlineFieldTot + nheaderT

      nline = nline  &
!        Lines back over intervening times, to end of required time.
          + (tPrev-tNext-1)*nlineT  &
!        Lines back over all later fields at required time, and the required field.
          + (nfieldFile-fNext+1)*nlinefieldTot &
!        Don't include any header for required field.
        - nheaderField

!     Reposition the file.
      DO l=1,nline
        BACKSPACE(unit,iostat=ierr)   !  backspace over a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: backspace',errMess2=errMess )
      ENDDO

    ELSE  !  inc>= 0
!     Move forward through file.
!     Note that inc=0 (meaning read next field) still may have to skip headers.

      IF ( inc > 0 ) THEN
!       Lines to end of current time.
        IF ( tNext-tPrev>0 .AND. tPrev>0 ) nline=nline+(nfieldFile-fPrev)*nlineFieldTot
!       Lines for all intervening times.
        nline = nline + (tNext-tPrev-1)*nlineT
      ENDIF

!     Lines for any file headers - needed if we are at start of file.
      IF ( tPrev==0 ) nline = nline + nheaderFile

      nline = nline  &
!       Lines for any time headers at required time.
        + nheaderT

      xfPrev = 1
      IF ( tPrev == tNext ) xfPrev=fPrev+1
      nline = nline  &
!       Lines for all earlier fields at required time.
        + (fNext-xfPrev)*nlineFieldTot  &
!       Lines for any field headers for required field.
        + nheaderField

!     Reposition the file.
      DO l=1,nline
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: reposition',errMess2=errMess )
      ENDDO

    ENDIF  !  inc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  ELSE  !  nlineField<1

!   Get the number of intervening fields.
    inc = fNextAcc - fPrevAcc - 1

    IF ( inc < 0 ) THEN
!     To move back through file, rewind and read from start.
      if (unit/=77) REWIND( unit )
      fPrevAcc = 0
    ENDIF

!   Work out where the last field read was in the file.
    tPrev = CEILING( REAL(fPrevAcc) / REAL(nfieldFile) )
    jt = MAX( 1, tPrev )   !   if fPrevAcc=0, want tPrev=0
    fPrev = fPrevAcc - (jt-1)*nfieldFile

!   Work out where the next field to read is in the file..
    tNext = CEILING( REAL(fNextAcc) / REAL(nfieldFile) )
    fNext = fNextAcc - (tNext-1)*nfieldFile

!   Note that inc=0 (meaning read next field) still may have to skip headers.

!   Read any file headers - needed if we are at start of file.
    IF ( tPrev==0 ) THEN
      DO l=1,nheaderFile
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: file header',errMess2=errMess )
      ENDDO
    ENDIF

!   Read to end of current time.
    IF ( tNext-tPrev>0 .AND. tPrev>0  ) THEN
      DO f=fPrev+1,nfieldFile
!       Read (skip) field headers
        DO l=1,nheaderField
          READ(unit,*,iostat=ierr)  !  read a line
!         If error, call error routine, indicating a stop.
          IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                    ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                    ,ierr=ierr,unit=unit,t=tNext &
                    ,errMess1='readIntAscii: current time field header',errMess2=errMess )
        ENDDO
!       Read a data field.
        READ(unit,*,iostat=ierr) val(:)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: current time data',errMess2=errMess )
      ENDDO
    ENDIF

!   Read all intervening times.
    DO t=tPrev+1,tNext-1
!     Read time headers.
      DO l=1,nheaderT
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii inter.time header',errMess2=errMess )

      ENDDO
      DO f=1,nfieldFile
!       Read (skip) field headers
        DO l=1,nheaderField
          READ(unit,*,iostat=ierr)  !  read a line
!         If error, call error routine, indicating a stop.
          IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                    ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                    ,ierr=ierr,unit=unit,t=tNext  &
                    ,errMess1='readIntAscii inter.field header',errMess2=errMess )
        ENDDO
!       Read a data field.
        READ(unit,*,iostat=ierr) val(:)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii inter.time data',errMess2=errMess )
      ENDDO
    ENDDO

!   Read any time headers at required time.
    IF ( tPrev < tNext ) THEN
      DO l=1,nheaderT
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: time header',errMess2=errMess )
      ENDDO
    ENDIF

!   Read all earlier fields at required time.
    xfPrev = 1
    IF ( tPrev == tNext ) xfPrev=fPrev+1
    DO f=xfPrev,fNext-1
!     Read field headers.
      DO l=1,nheaderField
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: prev.field header',errMess2=errMess )
      ENDDO
!     Read a data field.
      READ(unit,*,iostat=ierr) val(:)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                  ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
                  ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readIntAscii: prev.data',errMess2=errMess )
    ENDDO

!   Read any field headers for required field.
    DO l=1,nheaderField
      READ(unit,*,iostat=ierr)  !  read a line
!     If error, call error routine, indicating a stop.
      IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                ,ierr=ierr,unit=unit,t=tNext  &
                ,errMess1='readIntAscii: field header',errMess2=errMess )
    ENDDO

  ENDIF  !  nlineField
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! File is now correctly positioned to read requested field.
  READ(unit,*,iostat=ierr) val(:)

! If error, call error routine, indicating a stop.
  IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
            ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
            ,ierr=ierr,unit=unit,t=tNext  &
            ,errMess1='readIntAscii: data',errMess2=errMess )

! Update fPrevAcc to show latest field read.
  fPrevAcc = fNextAcc

  END SUBROUTINE readIntASCII
!#######################################################################
!#######################################################################

! subroutine readIntGrads
! Internal procedure in readWrite_mod.
! Read an integer variable from a GrADS file.
! NB  This reads a real variable that is then loaded into an integer, it is NOT reading
!     an integer from the file.

  SUBROUTINE readIntGrads( irec,unit,val,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    irec       &!  record number to be read
   ,unit        !  unit to read from

  INTEGER, INTENT(out) ::   &!  out ARRAYS
    val(:)      !  the data that are output

  INTEGER ::  &!  local SCALARS
    ierr

  REAL ::   &!  local ARRAYS
    rval( SIZE(val) )  !  the real data that are read in

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    errMess    !  message printed on error
!-------------------------------------------------------------------------------
  READ(unit,rec=irec,iostat=ierr) rval(:)

! If error, call error routine, indicating a stop.
  IF ( ierr /= 0 ) CALL rwErr( formatBin,.TRUE.,nval=SIZE(rval(:))  &
            ,irec=irec,ierr=ierr,unit=unit,errMess1='readIntGrads',errMess2=errMess )

  val(:) = NINT( rval(:) )

  END SUBROUTINE readIntGrads
!#######################################################################
!#######################################################################
! subroutine readRealASCII
! Internal procedure in readWrite_mod.
! Read a real variable from an ASCII file.
! Repositions file as necessary.
! It is assumed that any previous reads of this file read a whole field, so
! file is either positioned at start of file or at the end of a field.
! Will fail if previously read  only part of a field that extends over more than one line.
! This routine could be much shorter if time headers were not allowed - then
! we could just loop over fields.


  SUBROUTINE readRealASCII( fPrevAcc,fNextAcc,nfieldFile  &
                           ,nheaderFile,nheaderT,nheaderField,nlineField  &
                           ,unit,val,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    nfieldFile    &!  the number of fields (per time) in the file
   ,nheaderFile   &!  the number of header lines at the start of the file
   ,nheaderT      &!  the number of header lines at the start of each time
   ,nheaderField  &!  the number of header lines at the start of each field
   ,nlineField    &!  the number of lines per field (excl headers) - this is needed
                   !    if need to backspace within a file
!                       <1 means don't attempt to backspace, always move
!                       forward through file
   ,unit           !  unit to read from

  INTEGER, INTENT(inout) ::  &!  inout SCALARS
    fNextAcc   &!  the field number to be read. This is the number accumulated
!                    from start of file.
!                    0 means next field is read.
   ,fPrevAcc    !  the last field read from this file. This is the number
!                     accumulated from start of file.

  INTEGER ::  &!  local SCALARS
    f        &!  loop counter
   ,fNext    &!  the field (at t=tNext, i.e. not accumulated) to be read next
!                  (fields do not include headers)
   ,fPrev    &!  the field (at t=tPrev, i.e. not accumulated) that was last read
!                  (fields do not include headers)
   ,ierr     &!  error value
   ,inc      &!  the number of fields between fPrevAcc and fNextAcc
   ,jt       &!  work
   ,l        &!  loop counter
   ,nline    &!  work: number of lines to read or backspace
   ,nlineT    &!  the number of lines per time in file (incl headers)
   ,nlineFieldTot   &!  the number of lines per field per time (incl headers)
   ,t        &!  loop counter
   ,tNext    &!  the "time" level in the file of fNextAcc
   ,tPrev    &!  the "time" level in the file of fPrevAcc
   ,xfPrev    !  work: version of fPrev

  REAL, INTENT(out) ::   &!  out ARRAYS
    val(:)      !  the data that are read

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    errMess    !  message printed on error
!-------------------------------------------------------------------------------

! For now, I don't expect nlineField>0 (although the code "should" still work).
  IF ( nlineField > 0 ) THEN
    WRITE(*,*)'nlineField=',nlineField,' unit=',unit
    WRITE(*,*)'readRealASCII: nlineField>0: Unexpected at present!'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF

  IF ( nfieldFile < 1 ) THEN
!   Note that we need nfieldFile, even if we are attempting to read the
!   next record, since there is always a division with nfieldFile in denominator.
!   In future we might like to have a case (e.g. fNextAcc<0) which indicates
!   "REALLY read next record, not trying to avoid headers etc.".
    WRITE(*,*)'ERROR: readRealASCII: nfieldFile < 1'
    WRITE(*,*)'nfieldFile must be >= 0'
    WRITE(*,*) TRIM( errMess )
    STOP
  ENDIF

! Can't rewind standard input, so assume we will read next field.
! Ignore possibility of skipping any fields.
  IF ( unit == stdIn .OR. unit == jinUnit ) THEN
    fPrevAcc = 0
    fNextAcc = 0
  ENDIF

! fNextAcc=0 means read next field (still deals with headers).
  IF ( fNextAcc == 0 ) fNextAcc=fPrevAcc+1

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  IF ( nlineField >= 1 ) THEN

!   Knowing the number of lines of data per field, we can move forward or back
!   through the file.

    IF ( fNextAcc<fPrevAcc/2 .AND. fNextAcc>100 ) THEN
!   If next field is nearer start of file that current position, rewind and go
!   from start of file.
      if (unit/=77) REWIND( unit )
      fPrevAcc = 0
    ENDIF

!   Get the number of intervening fields.
    inc = fNextAcc - fPrevAcc - 1

!   Calculate number of lines per field, including headers.
    nlineFieldTot = nheaderField + nlineField
!   Calculate number of lines per time, including headers.
    nlineT = nheaderT + nfieldFile*nlineFieldTot

!   Work out where the last field read was in the file.
    tPrev = CEILING( REAL(fPrevAcc) / REAL(nfieldFile) )
    jt = MAX( 1, tPrev )   !   if fPrevAcc=0, want tPrev=0
    fPrev = fPrevAcc - (jt-1)*nfieldFile

!   Work out where the next field to read is in the file..
    tNext = CEILING( REAL(fNextAcc) / REAL(nfieldFile) )
    fNext = fNextAcc - (tNext-1)*nfieldFile

!   Calculate number of lines to be read or backspaced so as to reposition the
!   file as necessary.
    nline = 0
    IF ( inc < 0 ) THEN

!     Move back through the file.

!     Lines back to start of current time.
      IF ( tNext-tPrev < 0 ) nline = nline + fPrev*nlineFieldTot + nheaderT

      nline = nline  &
!        Lines back over intervening times, to end of required time.
          + (tPrev-tNext-1)*nlineT  &
!        Lines back over all later fields at required time, and the required field.
          + (nfieldFile-fNext+1)*nlinefieldTot &
!        Don't include any header for required field.
        - nheaderField

!     Reposition the file.
      DO l=1,nline
        BACKSPACE(unit,iostat=ierr)   !  backspace over a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: backspace',errMess2=errMess )
      ENDDO

    ELSE  !  inc>= 0
!     Move forward through file.
!     Note that inc=0 (meaning read next field) still may have to skip headers.

      IF ( inc > 0 ) THEN
!       Lines to end of current time.
        IF ( tNext-tPrev>0 .AND. tPrev>0 ) nline=nline+(nfieldFile-fPrev)*nlineFieldTot
!       Lines for all intervening times.
        nline = nline + (tNext-tPrev-1)*nlineT
      ENDIF

!     Lines for any file headers - needed if we are at start of file.
      IF ( tPrev==0 ) nline = nline + nheaderFile

      nline = nline  &
!       Lines for any time headers at required time.
        + nheaderT

      xfPrev = 1
      IF ( tPrev == tNext ) xfPrev=fPrev+1
      nline = nline  &
!       Lines for all earlier fields at required time.
        + (fNext-xfPrev)*nlineFieldTot  &
!       Lines for any field headers for required field.
        + nheaderField

!     Reposition the file.
      DO l=1,nline
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: reposition',errMess2=errMess )
      ENDDO

    ENDIF  !  inc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  ELSE  !  nlineField<1

!   Get the number of intervening fields.
    inc = fNextAcc - fPrevAcc - 1

    IF ( inc < 0 ) THEN
!     To move back through file, rewind and read from start.
      if (unit/=77) REWIND( unit )
      fPrevAcc = 0
    ENDIF

!   Work out where the last field read was in the file.
    tPrev = CEILING( REAL(fPrevAcc) / REAL(nfieldFile) )
    jt = MAX( 1, tPrev )   !   if fPrevAcc=0, want tPrev=0
    fPrev = fPrevAcc - (jt-1)*nfieldFile

!   Work out where the next field to read is in the file..
    tNext = CEILING( REAL(fNextAcc) / REAL(nfieldFile) )
    fNext = fNextAcc - (tNext-1)*nfieldFile

!   Note that inc=0 (meaning read next field) still may have to skip headers.

!   Read any file headers - needed if we are at start of file.
    IF ( tPrev==0 ) THEN
      DO l=1,nheaderFile
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: file header',errMess2=errMess )
      ENDDO
    ENDIF

!   Read to end of current time.
    IF ( tNext-tPrev>0 .AND. tPrev>0  ) THEN
      DO f=fPrev+1,nfieldFile
!       Read (skip) field headers
        DO l=1,nheaderField
          READ(unit,*,iostat=ierr)  !  read a line
!         If error, call error routine, indicating a stop.
          IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                   ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                   ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: current time field header',errMess2=errMess )
        ENDDO
!       Read a data field.
        READ(unit,*,iostat=ierr) val(:)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: current time data',errMess2=errMess )
      ENDDO
    ENDIF

!   Read all intervening times.
    DO t=tPrev+1,tNext-1
!     Read time headers.
      DO l=1,nheaderT
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: inter.time time header',errMess2=errMess )
      ENDDO
      DO f=1,nfieldFile
!       Read (skip) field headers
        DO l=1,nheaderField
          READ(unit,*,iostat=ierr)  !  read a line
!         If error, call error routine, indicating a stop.
          IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                   ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                   ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: inter time field header',errMess2=errMess )
        ENDDO
!       Read a data field.
        READ(unit,*,iostat=ierr) val(:)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.,nval=SIZE(val(:))  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: inter.time data',errMess2=errMess )
      ENDDO
    ENDDO

!   Read any time headers at required time.
    IF ( tPrev < tNext ) THEN
      DO l=1,nheaderT
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: time header',errMess2=errMess )
      ENDDO
    ENDIF

!   Read all earlier fields at required time.
    xfPrev = 1
    IF ( tPrev == tNext ) xfPrev=fPrev+1
    DO f=xfPrev,fNext-1
!     Read field headers.
      DO l=1,nheaderField
        READ(unit,*,iostat=ierr)  !  read a line
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
                 ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
                 ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: prev.field header',errMess2=errMess )
      ENDDO
!     Read a data field.
      READ(unit,*,iostat=ierr) val(:)
!     If error, call error routine, indicating a stop.
      IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
               ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
               ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: prev.data',errMess2=errMess )
    ENDDO

!   Read any field headers for required field.
    DO l=1,nheaderField
      READ(unit,*,iostat=ierr)  !  read a line
!     If error, call error routine, indicating a stop.
      IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
               ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
               ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: field header',errMess2=errMess )
    ENDDO

  ENDIF  !  nlineField
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

! File is now correctly positioned to read requested field.
  READ(unit,*,iostat=ierr) val(:)

! If error, call error routine, indicating a stop.
  IF ( ierr /= 0 ) CALL rwErr( formatAsc,.TRUE.  &
           ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile,nval=SIZE(val(:))  &
           ,ierr=ierr,unit=unit,t=tNext  &
                  ,errMess1='readRealAscii: data',errMess2=errMess )

! Update fPrevAcc to show latest field read.
  fPrevAcc = fNextAcc

  END SUBROUTINE readRealASCII
!#######################################################################
!#######################################################################
! subroutine readRealGrads
! Internal procedure in readWrite_mod.
! Read a real variable from a GrADS file.

  SUBROUTINE readRealGrads( irec,unit,val,errMess )

  IMPLICIT NONE

  INTEGER, INTENT(in) ::  &!  in SCALARS
    irec       &!  record number to be read
   ,unit        !  unit to read from

  INTEGER ::  &!  local SCALARS
    ierr

  REAL, INTENT(out) ::   &!  out ARRAYS
    val(:)      !  the data that are read

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    errMess    !  message printed on error
!-------------------------------------------------------------------------------

  READ(unit,rec=irec,iostat=ierr) val(:)

! If error, call error routine, indicating a stop.
  IF ( ierr /= 0 ) CALL rwErr( formatBin,.TRUE.,nval=SIZE(val(:))  &
            ,irec=irec,ierr=ierr,unit=unit,errMess1='readRealGrads',errMess2=errMess )

  END SUBROUTINE readRealGrads
!#######################################################################
!#######################################################################
! subroutine readRealPP
! Internal procedure in readWrite_mod.
! Read a real variable from a PP file (sequential).
! Repositions file as necessary.
! It is assumed that any previous reads of this file read a whole field, so
! file is either positioned at start of file or at the end of a field.


  SUBROUTINE readRealPP( fieldCode,nfieldFile,unit,errMess  &
                        ,fNextAcc,fPrevAcc  &
                        ,val )

  IMPLICIT NONE

!--------------------------------------------------
! Scalar parameters.
  INTEGER, PARAMETER ::  &
    niHeader = 45  &!  number of words in integer header
   ,nrHeader = 19   !  number of words in real header

!--------------------------------------------------
! Scalars with intent(in)
  INTEGER, INTENT(in) ::  &
    fieldCode     &!  the field code (e.g. STASH) for the required field
!                        (used as a check only). The PP and STASH codes are
!                        available from the headers, so either can be used
!                        At present the STASH code is used.
   ,nfieldFile    &!  the number of fields (per time) in the file
   ,unit           !  unit to read from

  CHARACTER(len=*), INTENT(in) ::  &
    errMess    !  message printed on error

!--------------------------------------------------
! Scalars with intent(inout)
  INTEGER, INTENT(inout) ::  &
    fNextAcc   &!  the field number to be read. This is the number accumulated
!                    from start of file.
   ,fPrevAcc    !  the last field read from this file. This is the number accumulated
!                    from start of file.

!--------------------------------------------------
! Arrays with intent(out)
  REAL, INTENT(out) ::   &
    val(:)      !  the data that are read

!--------------------------------------------------
! Local scalars.
  INTEGER ::  &
    f        &!  loop counter
   ,ierr     &!  error value
   ,nval     &!  size of data field in file (including "extra" data after "main" data)
   ,nxFile   &!  x extent of grid (from file header)
   ,nyFile   &!  y extent of grid (from file header)
   ,fieldCodeFile  !  field code from file header

  LOGICAL ::  &
    errFound  !  flag indicating an error
!--------------------------------------------------
! Local arrays.
  INTEGER :: &
    iHeader(niHeader)  !  integer header

  REAL :: &!  local arrays
    rHeader(nrHeader)  !  real header
!-------------------------------------------------------------------------------

  IF ( fNextAcc <= fPrevAcc ) THEN
!   To move back through file, rewind and read from start.
    if (unit/=77) REWIND( unit )
    fPrevAcc = 0
  ENDIF

! Loop over fields (including desired field).
  DO f=fPrevAcc+1,fNextAcc

!-------------------------------------------------------------------------------
!   Read headers.
!-------------------------------------------------------------------------------
    READ(unit,iostat=ierr) iHeader(:),rHeader(:)
!   If error, call error routine, indicating a stop.
    IF ( ierr /= 0 ) CALL rwErr( formatPP,.TRUE.  &
               ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
               ,fieldCode=fieldCode,unit=unit  &
               ,errMess1='readRealPP: header read',errMess2=errMess )

!   Get information from header.
    nyFile = iHeader(18)
    nxFile = iHeader(19)
    fieldCodeFile = iHeader(42)   !  42 is STASH code. 23 is PP code.
    nval = iHeader(15)

    errFound = .FALSE.

!   Check that field size is as expected - just check total size.
    IF ( nxFile*nyFile /= SIZE(val(:)) ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: readRealPP: field size not as expected.'
      WRITE(*,*)'Required field size=',SIZE(val(:)),' points.'
      WRITE(*,*)'PP header indicates nx,ny=',nxFile,nyFile,' giving nx*ny=',nxFile*nyFile
    ENDIF

!   If this is the required field, heck that field code is as expected.
    IF ( f==fNextAcc .AND. fieldCodeFile/=fieldCode ) THEN
      errFound = .TRUE.
      WRITE(*,*)'ERROR: readRealPP: field code not as expected.'
      WRITE(*,*)'Required field code=',fieldCode,' found ',fieldCodeFile
    ENDIF

!   Act on error. Call error handling procedure, indicating a stop.
    IF ( errFound ) CALL rwErr( formatPP,.TRUE.  &
               ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
               ,fieldCode=fieldCode,unit=unit  &
               ,errMess1='readRealPP: header error',errMess2=errMess )

!-------------------------------------------------------------------------------
!   Read a data field.
!-------------------------------------------------------------------------------
    CALL readRealPPData( nval,unit,ierr,val )
!   If error, call error routine, indicating a stop.
    IF ( ierr /= 0 ) CALL rwErr( formatPP,.TRUE.  &
               ,irec=fnextAcc,irecPrev=fprevAcc,nfieldFile=nfieldFile  &
               ,fieldCode=fieldCode,unit=unit  &
               ,errMess1='readRealPP: data read',errMess2=errMess )

!   Update fPrevAcc to show latest field read.
    fPrevAcc = fPrevAcc + 1

  ENDDO

  END SUBROUTINE readRealPP
!#######################################################################
!#######################################################################
! subroutine readRealPPData
! Internal procedure in readWrite_mod.
! Read a PP data record (not headers).

  SUBROUTINE readRealPPData( nval,unit,ierr,valOut )
! Read a record of nval values, a subset (or all of which) are returned.
! This is done in a subroutinbe to make it easier to read any "extra" data
! after the "main" data.

  IMPLICIT NONE

! Scalar arguments with intent(in)
  INTEGER, INTENT(in) ::  &
    nval    &!  the length of the record to be read
!                 This includes any "extra" data after the main data record.
   ,unit     !  the unit connected to the PP file

! Scalar arguments with intent(out)
  INTEGER, INTENT(out) ::  &
    ierr   !  error code from read

! Array arguments with intent(out)
  REAL, INTENT(out) ::  &
    valOut(:)    !  the output data (all or part of val)
!                     size(valout) is known to be <=size(val)

! Local arrays.
  REAL ::  &
    val(nval)   !  data record from file (including any "extra" data)
!-------------------------------------------------------------------------------

! Read the full data record.
  READ(unit,iostat=ierr) val(:)

! Extract required data.
  IF ( ierr == 0 ) valOut(1:SIZE(valOut)) = val(1:SIZE(valOut))

  END SUBROUTINE readRealPPData
!#######################################################################
!#######################################################################
! subroutine writeVarInt
! Internal procedure in readWrite_mod.
! Module procedure with generic name writeVar.

  SUBROUTINE writeVarInt( irec,unit,val,fileFormat,errMess )

  INTEGER, INTENT(in) ::  &!  in SCALARS
    irec      &!  record number to write, only used for unformatted
   ,unit       !  unit to write to

  INTEGER, INTENT(in) ::  &!  in ARRAYS
    val(:)  !  data

  INTEGER ::  &!  local SCALARS
    i,ierr,j1,j2   &!
   ,nline          &!  number of lines to write
   ,np              !  length of input data

  CHARACTER(len=*), INTENT(in) ::  &!  in SCALARS
    errMess     &!  message printed on error
   ,fileFormat   !  format (type) of output file
!   ,writeFormat

!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )
!     ASCII output

!      write(unit,fmt=writeFormat) val
!      WRITE(unit,*,iostat=ierr) val(:)   !  until get better format!

!     Get size of input data.
      np = SIZE( val(:) )
!     Work out how many lines to be written.
      nline = CEILING( REAL(np) / REAL(lineLen) )
      j2 = 0
      DO i=1,nline
        j1 = j2 + 1
        j2 = MIN( j1+lineLen-1, np )
        WRITE(unit,*,iostat=ierr) val(j1:j2)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( fileFormat,.TRUE.  &
                 ,ierr=ierr,unit=unit,nval=j2-j1+1  &
                 ,errMess1='writeVarReal: formatAsc',errMess2=errMess )
      ENDDO

    CASE ( formatBin )
!     GrADS output (direct access)
      WRITE(unit,rec=irec,iostat=ierr) val(:)
!     If error, call error routine, indicating a stop.
      IF ( ierr /= 0 ) CALL rwErr( fileFormat,.TRUE.  &
               ,ierr=ierr,irec=irec,nval=SIZE(val(:)),unit=unit  &
               ,errMess1='writeVarInt: formatBin',errMess2=errMess )

    CASE ( formatNc )
!     netCDF output
      WRITE(*,*)'SORRY! netCDF code does not yet exist!'
      WRITE(*,*)'Stopping in writeVarInt'
      STOP

    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in writeVarInt'
      STOP

  END SELECT

  END SUBROUTINE writeVarInt
!#######################################################################
!#######################################################################
! subroutine writeVarReal
! Internal procedure in readWrite_mod.
! Module procedure with generic name writeVar.

  SUBROUTINE writeVarReal( irec,varID,unit,ncCount,ncStart,val  &
                          ,fileFormat,errMess,mapArg )

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) ::  &
    irec      &!  record number to write, only used for unformatted
   ,varID     &!  netCDF ID of variable
   ,unit       !  unit to write to (or netCDF ID of file)

  CHARACTER(len=*), INTENT(in) ::  &
    errMess     &!  message printed on error
   ,fileFormat   !  format (type) of output file

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: ncCount(:)  !  count argument for netCDF variable
  INTEGER, INTENT(in) :: ncStart(:)  !  start argument for netCDF variable
  REAL, INTENT(in) :: val(:)         !  data

!-------------------------------------------------------------------------------
! Optional array arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in), OPTIONAL :: mapArg(:) ! map argument for netCDF variable

!-------------------------------------------------------------------------------
! Local scalar variables
!-------------------------------------------------------------------------------
  INTEGER ::  &
    i,idim,ierr,j1,j2   &!
   ,nline          &!  number of lines to write
   ,np              !  length of input data

!-------------------------------------------------------------------------------
! Local array variables
!-------------------------------------------------------------------------------
  INTEGER :: map(SIZE(ncCount))  !  map argument for nf90_put_var

!-------------------------------------------------------------------------------
  SELECT CASE ( fileFormat )
    CASE ( formatAsc )
!     ASCII output

!      write(unit,fmt=writeFormat) val
!      WRITE(unit,*,iostat=ierr) val(:)   !  until get better format!

!     Get size of input data.
      np = SIZE( val(:) )
!     Work out how many lines to be written.
      nline = CEILING( REAL(np) / REAL(lineLen) )
      j2 = 0
      DO i=1,nline
        j1 = j2 + 1
        j2 = MIN( j1+lineLen-1, np )
        WRITE(unit,*,iostat=ierr) val(j1:j2)
!       If error, call error routine, indicating a stop.
        IF ( ierr /= 0 ) CALL rwErr( fileFormat,.TRUE.  &
                 ,ierr=ierr,unit=unit,nval=j2-j1+1  &
                 ,errMess1='writeVarReal: formatAsc',errMess2=errMess )
      ENDDO

    CASE ( formatBin )
!     GrADS output (direct access)
      WRITE(unit,rec=irec,iostat=ierr) val(:)
!     If error, call error routine, indicating a stop.
      IF ( ierr /= 0 ) CALL rwErr( fileFormat,.TRUE.  &
               ,ierr=ierr,irec=irec,unit=unit,nval=SIZE(val(:))  &
               ,errMess1='writeVarReal: formatBin',errMess2=errMess )

    CASE ( formatNc )
!     netCDF output
!     Get map.
      IF ( PRESENT(mapArg) ) THEN
        map(:) = mapArg(:)
      ELSE
!       Create a map that will do nothing!
        map(1) = 1
        DO idim=2,SIZE(ncCount)
          map(idim) = map(idim-1) * ncCount(idim-1)
        ENDDO
      ENDIF

      ierr = nf90_put_var( unit,varID,val,start=ncStart,count=ncCount,map=map )
!     If error, call error routine, indicating a stop.
      IF ( ierr /= nf90_noerr ) CALL rwErr( fileFormat,.TRUE.  &
               ,ierr=ierr,ncID=unit,nval=SIZE(val(:)),varID=varID  &
               ,errMess1='writeVarReal: formatNc',errMess2=errMess  &
               ,ncCount=ncCount,ncStart=ncStart,ncMap=map )


    CASE default
      WRITE(*,*)'Do not recognise fileFormat=',TRIM(fileFormat)
      WRITE(*,*)'Stopping in writeVarReal'
      STOP

  END SELECT

  END SUBROUTINE writeVarReal

!###############################################################################
!###############################################################################

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SUBPROGRAM: native_4byte_real
!
!         AUTHOR: David Stepaniak, NCAR/CGD/CAS
! DATE INITIATED: 29 April 2003
!  LAST MODIFIED: 29 April 2003
!
!       SYNOPSIS: Converts a 32 bit, 4 byte, REAL from big Endian to
!                 little Endian, or conversely from little Endian to big
!                 Endian.
!
!    DESCRIPTION: This subprogram allows one to convert a 32 bit, 4 byte,
!                 REAL data element that was generated with, say, a big
!                 Endian processor (e.g. Sun/sparc, SGI/R10000, etc.) to its
!                 equivalent little Endian representation for use on little
!                 Endian processors (e.g. PC/Pentium running Linux). The
!                 converse, little Endian to big Endian, also holds.
!                 This conversion is accomplished by writing the 32 bits of
!                 the REAL data element into a generic 32 bit INTEGER space
!                 with the TRANSFER intrinsic, reordering the 4 bytes with
!                 the MVBITS intrinsic, and writing the reordered bytes into
!                 a new 32 bit REAL data element, again with the TRANSFER
!                 intrinsic. The following schematic illustrates the
!                 reordering process
!
!
!                  --------    --------    --------    --------
!                 |    D   |  |    C   |  |    B   |  |    A   |  4 Bytes
!                  --------    --------    --------    --------
!                                                             |
!                                                              -> 1 bit
!                                       ||
!                                     MVBITS
!                                       ||
!                                       \/
!
!                  --------    --------    --------    --------
!                 |    A   |  |    B   |  |    C   |  |    D   |  4 Bytes
!                  --------    --------    --------    --------
!                         |           |           |           |
!                         24          16          8           0   <- bit
!                                                                 position
!
!          INPUT: realIn,  a single 32 bit, 4 byte REAL data element.
!         OUTPUT: realOut, a single 32 bit, 4 byte REAL data element, with
!                 reverse byte order to that of realIn.
!    RESTRICTION: It is assumed that the default REAL data element is
!                 32 bits / 4 bytes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE native_4byte_real( realIn, realOut )
  IMPLICIT NONE
  REAL, INTENT(IN)                              :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
  REAL, INTENT(OUT)                             :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Local variables (generic 32 bit INTEGER spaces):

  INTEGER                                       :: i_element
  INTEGER                                       :: i_element_br
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer 32 bit of realIn to generic 32 bit INTEGER space:
  i_element = TRANSFER( realIn, 0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reverse order of 4 bytes in 32 bit INTEGER space:
  CALL MVBITS( i_element, 24, 8, i_element_br, 0  )
  CALL MVBITS( i_element, 16, 8, i_element_br, 8  )
  CALL MVBITS( i_element,  8, 8, i_element_br, 16 )
  CALL MVBITS( i_element,  0, 8, i_element_br, 24 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transfer reversed order bytes to 32 bit REAL space (realOut):
  realOut = TRANSFER( i_element_br, 0.0 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  END SUBROUTINE native_4byte_real

!###############################################################################
!###############################################################################

  END MODULE readWrite_mod
