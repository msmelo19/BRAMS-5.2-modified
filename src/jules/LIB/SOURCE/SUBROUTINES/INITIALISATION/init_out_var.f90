!###############################################################################
!###############################################################################
! subroutine init_out_var
! Further processing of output variables.

  SUBROUTINE init_out_var( nvarMax,tmpLoc )

  USE ancil_info, ONLY :  &
!   imported scalars with intent(in)
     dim_cs1,land_pts,nsmax,nsoil=>sm_levels,ntiles,nx=>row_length,ny=>rows

  USE inout, ONLY : formatNc,gradsNc,havePFT,haveSCpool,haveSnow,haveSoil  &
                   ,haveTile,haveType,mapOut,nlevMax,nlevMaxCtl,nout  &
                   ,nvarOut,nvarOutTot,nxyMax  &
                   ,outDateFlag,outEndPos,outLenWrite  &
                   ,outLen,outFormat,pointsOut,pointsOutLand,rpProfile  &
                   ,snapProfile,taccumVar,tmeanProfile,tmeanVar,usePseudo  &
                   ,varname,varNameList,varNlev,varNum,varPos,varStartPos,varType

  USE nstypes, ONLY :  &
!   imported scalars with intent(in)
     npft,ntype

  USE offline_diag, ONLY :  &
!  imported scalar parameters
     offDiag  &
!  imported scalars with intent(inout)
    ,useCiDiag,useGstomDiag,useRdcDiag,useRFlowDiag  &
    ,useRoffInfDiag,useRRunDiag  &
    ,useSnowGMeltDiag,useWfluxDiag,useWfluxSfcDiag

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     nxRoute,nyRoute

!-------------------------------------------------------------------------------
  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar arguments with intent(in)
!-------------------------------------------------------------------------------
  INTEGER, INTENT(in) :: nvarMax   !  size of array

!-------------------------------------------------------------------------------
! Array arguments with intent(in)
!-------------------------------------------------------------------------------
  CHARACTER(len=*), INTENT(in) :: tmpLoc(nvarMax,2)

!-------------------------------------------------------------------------------
! Local scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
    i,iout,ip,ivar,jvar,k,kvar   &!  work
   ,nlevCtl  &!  number of levels in a variable (for use in GrADS ctl file)
   ,np       &!  number of points in a variable (full model domain)
   ,npOut     !  number of points to be output for a variable

  LOGICAL :: errFlag,errFlagVar      !  work

  CHARACTER(len=LEN(varname)) :: varnameVal   !  work

!-------------------------------------------------------------------------------
! Local array variables.
!-------------------------------------------------------------------------------
  LOGICAL :: done(nvarOutTot)   !  T when a variable name has been changed to avoid clash


!-------------------------------------------------------------------------------

  WRITE(*,"(50('-'),/,a)") 'init_out_var'

!-------------------------------------------------------------------------------
! Set flags indicating type of time processing required for a profile.
! If only initial state is requested, change variables to snapshot.
!-------------------------------------------------------------------------------
  DO iout=1,nout

    errFlag = .FALSE.

    DO jvar=1,nvarOut(iout)
      ivar = varPos(iout,jvar)   !  position in list of chosen variables

      IF ( taccumVar(ivar) .OR. tmeanVar(ivar) ) THEN

!       Indicate that this profile includes time-averaged variables.
        tmeanProfile(iout) = .TRUE.

!       Reset flags if only initial state required.
        IF ( outDateFlag(iout) == -2 ) THEN
          taccumVar(ivar) = .FALSE.
          tmeanVar(ivar) = .FALSE.
          errFlag = .TRUE.
        ENDIF
  
      ELSE

!       Indicate that this profile includes snapshot variables.
        snapProfile(iout) = .TRUE.

      ENDIF

    ENDDO  !  jvar

!   Print a warning if flags were changed.
    IF ( errFlag ) THEN
      WRITE(*,*)'WARNING: Output profile #',iout,' - only the initial state is requested.'
      WRITE(*,*)'Time average or accumulation variables were changed to snapshots.'
      WRITE(*,*)'This may create repeated variables.'
    ENDIF

  ENDDO  !  iout

!-------------------------------------------------------------------------------
! Check that the same combination of variable and accumulation/mean/instantaneous value is not
! repeated in any one profile - i.e. that there are no doublers.
! Variables of type RP (selected points from routing grid) are identical only
! if the location is also identical.
!-------------------------------------------------------------------------------
  DO iout=1,nout
    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)
!     Compare with all later variables.
      DO i=ivar+1,nvarOut(iout)
        kvar = varPos(iout,i)  !  position in overall list
        IF ( varNum(jvar) == varNum(kvar) ) THEN   !  same var
          IF ( (taccumVar(jvar).EQV.taccumVar(kvar)) .AND. (tmeanVar(jvar).EQV.tmeanVar(kvar)) ) THEN
            IF ( .NOT. rpProfile(iout) ) THEN
              WRITE(*,*)'Repeated variable in output profile #',iout
              WRITE(*,*)'Variables #',ivar,i
              WRITE(*,*) 'Variable description (taccumVar,tmeanVar,name): ',taccumVar(jvar)  &
                        ,tmeanVar(jvar),varNameList(varNum(jvar))
              WRITE(*,*) 'init_out_var: repeated variable'
              STOP
            ELSEIF ( varType(jvar)=='RP' .AND. varType(kvar)=='RP' ) THEN
              IF ( mapOut(iout,ivar,1) == mapOut(iout,i,1) ) THEN
                WRITE(*,*)'ERROR: init_out_var: Repeated variable in output profile #',iout
                WRITE(*,*)'Variables #',ivar,i
                WRITE(*,*) 'Variable description (tmeanVar,name): ',tmeanVar(jvar),varNameList(varNum(jvar))
                WRITE(*,*) 'Variable description (location): ',TRIM(tmpLoc(jvar,1)),' ',TRIM(tmpLoc(jvar,2))
                WRITE(*,*) 'init_out_var'
                STOP
              ENDIF
!             At this point it is possible that we have point routing variables at different locations but
!             with identical descriptions. These will be differentiated in output by appending information
!             about the location to the descriptions.
            ENDIF
          ENDIF !  same time ave etc
        ENDIF   !  same var
      ENDDO     !  i (var)
    ENDDO  !  ivar
  ENDDO    !  iout

!-------------------------------------------------------------------------------
! Check if variable names are repeated. If names are repeated, but the variables in fact differ in
! time-processing, add suffix to each name to distiniguish them.  e.g. sthu snapshot and sthu
! time-average selected, by default these both get varname=sthu, but change these to sthuS and sthuM.
! The new name could now clash with the name of a different variable - this is checked later.
! This could all be done better, but is OK for now.
!-------------------------------------------------------------------------------

  done(:) = .FALSE.

  DO iout=1,nout

    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)  !  position in overall list
      varnameVal = varname(jvar)
!     Compare with all later variables.
      DO i=ivar+1,nvarOut(iout)
        kvar = varPos(iout,i)  !  position in overall list
        IF ( varnameVal == varname(kvar) ) THEN
!         Same name is given.
!         Add suffix to each name. Only change a name once!
          DO k=1,2
            ip = jvar
            IF ( k == 2 ) ip = kvar
            IF ( done(ip) ) CYCLE  !  name already changed
            np = LEN_TRIM(varname(ip))
            IF ( taccumVar(ip) ) THEN
              WRITE(*,*)'WARNING: Changing variable name to avoid clash.'
              WRITE(*,*) TRIM(varname(ip)),' => ',varname(ip)(1:np) // 'A'
              varname(ip) = varname(ip)(1:np) // 'A'
            ENDIF
            IF ( tmeanVar(ip) ) THEN
              WRITE(*,*)'WARNING: Changing variable name to avoid clash.'
              WRITE(*,*) TRIM(varname(ip)),' => ',varname(ip)(1:np) // 'M'
              varname(ip) = varname(ip)(1:np) // 'M'
            ENDIF
            IF ( .NOT.taccumVar(ip) .AND. .NOT.tmeanVar(ip) ) THEN
              WRITE(*,*)'WARNING: Changing variable name to avoid clash.'
              WRITE(*,*) TRIM(varname(ip)),' => ',varname(ip)(1:np) // 'S'
              varname(ip) = varname(ip)(1:np) // 'S'
            ENDIF
            done(ip) = .TRUE.
          ENDDO
        ENDIF   !  same name
      ENDDO     !  i (vars)
    ENDDO       !  ivar
  ENDDO         !  iout (profiles)

!-------------------------------------------------------------------------------
! In the interests of moderate rigour....check that we haven't changed the name of a variable
! to something that now clashes with the another name.
!-------------------------------------------------------------------------------
  DO iout=1,nout
    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)  !  position in overall list
!     Compare with all earlier variables.
      DO i=1,ivar-1
        kvar = varPos(iout,i)  !  position in overall list
        IF ( varname(jvar) == varname(kvar) ) THEN
          WRITE(*,*)'ERROR: init_out_var: repeated variable name: ',TRIM(varname(jvar))
          WRITE(*,*)'Output profile #',iout,' variables numbered ',i,ivar
          STOP
        ENDIF
      ENDDO
    ENDDO
  ENDDO 

!-------------------------------------------------------------------------------
! Work out if there is to be a dimension for "pseudo" layers.
! This dimension is only used for netCDF files that are not to be GrADS-readable.
! In all other cases, the "z" dimension is used for "pseudo".
! This is done in a separate loop so that the information is available for
! later loops.
!-------------------------------------------------------------------------------
  usePseudo = .FALSE.
  IF ( outFormat==formatNc .AND. .NOT. gradsNc )THEN
    DO iout=1,nout
      DO ivar=1,nvarOut(iout)
        SELECT CASE ( varType(varPos(iout,ivar)) )
          CASE ( 'SC','SN','PF','TI','TY' )
            usePseudo = .TRUE.
!           Exit as soon as detected a need for pseudo dimension.
            EXIT
        END SELECT
      ENDDO
    ENDDO
  ENDIF

!-------------------------------------------------------------------------------
! Establish how much space is needed to store each variable.
! Also find out what sizes (eg for work space) are needed for each profile.
! And what the maximum sizes for levels and pseudo dimensions in netCDF are.
!-------------------------------------------------------------------------------

  havePFT(:) = .FALSE.
  haveSCpool(:) = .FALSE.
  haveSnow(:) = .FALSE.
  haveSoil(:) = .FALSE.
  haveTile(:) = .FALSE.
  haveType(:) = .FALSE.

  DO iout=1,nout
    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)  !  position in list of chosen vars
      kvar = varNum( jvar )     !  position on list of all available vars
!     Save the location in outval of the first datum of this variable.
      varStartPos(jvar) = outLen + 1
!     Save number of levels for this variable, and update total space needed.
      SELECT CASE ( varType(jvar) )
        CASE ( 'LA' )
          varNlev(jvar) = 1
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + pointsOutLand(iout)
        CASE ( 'PF' )
          havePFT(iout) = .TRUE.
          varNlev(jvar) = npft
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + npft*pointsOutLand(iout)
        CASE ( 'RP' )
          varNlev(jvar) = 1
          np = nxRoute * nyRoute
          npOut = 1
          outLen = outLen + 1
        CASE ( 'RG' )
          varNlev(jvar) = 1
          nlevCtl = varNlev(jvar)
          np = nxRoute * nyRoute
          npOut = pointsOut(iout)
          outLen = outLen + pointsOut(iout)
        CASE( 'SC' )
          haveSCPool(iout) = .TRUE.
          varNlev(jvar) = dim_cs1
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + dim_cs1*pointsOutLand(iout)
        CASE ( 'SI' )
          varNlev(jvar) = 1
          nlevCtl = varNlev(jvar)
          np = nx*ny
          npOut = pointsOut(iout)
          outLen = outLen + pointsOut(iout)
        CASE ( 'SN' )
          haveSnow(iout) = .TRUE.
!         Combine tile and snow layers in varNlev.
          varNlev(jvar) = ntiles*nsmax
          nlevCtl = nsmax  !  since a ctl not used/useful if both pseudo and z dims exist
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + ntiles*nsMax*pointsOut(iout)
        CASE( 'SO' )
          haveSoil(iout) = .TRUE.
          varNlev(jvar) = nsoil
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + nsoil*pointsOutLand(iout)
        CASE ( 'TI' )
          haveTile(iout) = .TRUE.
          varNlev(jvar) = ntiles
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + ntiles*pointsOutLand(iout)
        CASE ( 'TY' )
          haveType(iout) = .TRUE.
          varNlev(jvar) = ntype
          nlevCtl = varNlev(jvar)
          np = land_pts
          npOut = pointsOutLand(iout)
          outLen = outLen + ntype*pointsOutLand(iout)
        CASE default
          WRITE(*,*)'ERROR: init_out_var: varType=',varType(jvar)
          WRITE(*,*)'No code for this varType.'
          STOP
      END SELECT

!     Check if this variable is largest so far.
      IF ( varNlev(jvar) > nlevMax(iout) ) nlevMax(iout) = varNlev(jvar)
      IF ( nlevCtl > nlevMaxCtl(iout) ) nlevMaxCtl(iout) = nlevCtl
      IF ( np > nxyMax(iout) ) nxyMax(iout) = np

!     Save values that relate to size of output for this profile.
      IF ( ivar == nvarOut(iout) ) THEN
!       Save the last storage location used by each profile.
        outEndPos(iout) = varStartPos(jvar) + npOut*varNlev(jvar) - 1
!       Check if this profile needs the largest amount of space in output vector.
        IF ( iout == 1 ) THEN
          outLenWrite = outEndPos(iout)
        ELSE
          IF ( outEndPos(iout)-outEndPos(iout-1) > outLenWrite )  &
                  outLenWrite = outEndPos(iout)-outEndPos(iout-1)
        ENDIF
      ENDIF

    ENDDO  !  ivar

  ENDDO    !  iout

!-------------------------------------------------------------------------------
! Set flags and count up space required for "offline" diagnostics.
! We only need to allocate this space once (not separately for every profile that
! requires the same data).
! However, it's not an efficient use of space if we only want a few points from a large field,
! as at present we save the entire field!
! NB varNameLists used here must, of course, agree with those in init_out_varlist.
!
! But....we're not allowed these "offline diagnostics" in "official" JULES, since
! they have to be loaded in the "physics" routines which have to be compatible
! with the UM. The switch "offDiag" indicates if these diagnostics are available.
! If the diags are not available but have been requested, raise an error.
!-------------------------------------------------------------------------------

  errFlag = .FALSE.

  DO iout=1,nout
    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)

      errFlagvar = .FALSE.

      SELECT CASE ( varNameList(varNum(jvar)) )

        CASE ( 'ciP' )
          useCiDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'gstomP' )
          useGstomDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'rdcP' )
          useRdcDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'rflow' )
          useRflowDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'rInflow' )
!         This is derived from the river (out)flow diagnostic, so indicate that that is needed.
          useRflowDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'surfRoffInf' )
          useRoffInfDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'runoffR' )
          useRrunDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'snowGroundMeltT' )
          useSnowGMeltDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'snowTotMeltT' )
!         Note that if snowTotMeltT is selected, we set useSnowGMeltDiag to get
!         snow_ground contribution.
          useSnowGMeltDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'wFlux' )
          useWfluxDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

        CASE ( 'wFluxSfc' )
          useWfluxSfcDiag = .TRUE.
          IF ( .NOT. offDiag ) errFlagVar = .TRUE.

      END SELECT

      IF ( errFlagVar ) THEN
         errFlag = .TRUE.
         WRITE(*,*)'SORRY ',TRIM(varNameList(varNum(jvar)))  &
              ,' is not available in this version of code. See subroutine init_out_var.'
      ENDIF

    ENDDO  !  ivar
  ENDDO  !  iout

!-------------------------------------------------------------------------------
! If an error has been raised from the offline diagnostics, stop.
!-------------------------------------------------------------------------------
  IF ( errFlag ) THEN
    WRITE(*,*) 'ERROR: init_out_var: offline diagnostics not available.'
    WRITE(*,*) 'These "offline" diagnostics are disabled in "official" versions of JULES.'
    STOP
!   If you would like to use these diagnostics, contact the JULES developers!
  ENDIF

  END SUBROUTINE init_out_var
!###############################################################################
!###############################################################################
