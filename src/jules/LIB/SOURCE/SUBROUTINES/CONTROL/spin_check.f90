!###############################################################################
!###############################################################################
! subroutine spin_check
! Assess whether model has spun up to within given tolerance.
!
! Note that at present if any point has not spun up, all points continue
! to be modelled - it is not just the "not spun up" points that are
! modelled. This is essential if there is any communication between
! gridpoints (e.g. river flow depends upon upstream conditions), (and is
! by far the easiest option if variables are on different grids - e.g. if
! routing is on a coarser grid, would have to work out which model
! gridboxes are covered by each routing box that has not spun up), but
! otherwise we could move to just modelling those points that are still
! spinning up (or do not use any "communicating" variables, such as river
! flow, to assess spin up).
!
!###############################################################################
!###############################################################################
SUBROUTINE Spin_Check

  USE alloc_mod, ONLY :  &
!  imported procedures
     allocate_arrays

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts,sm_levels

  USE grid_utils, ONLY :  &
!  imported procedures
     getXYPos

  USE prognostics, ONLY :  &
!  imported arrays with intent(in)
     routeStore,smcl,t_soil

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     npRoute,nxRoute,nyRoute  &
!  imported arrays with intent(in)
    ,routeIndex

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     iposRouteStore,iposSmcl,iposTsoil,ispin,nspin,nspinVarMax,spinFail  &
!  imported scalars with intent(inout)
     ,npSpinMax,nspinFinal,nspinVar,nzSpinMax,spinEnd,spinUp  &
!  imported arrays with intent(in)
     ,spinTol,spinTolPercent,spinVarName  &
!  imported arrays with intent(inout)
     ,spinValOld,spinvar,spinVarNp,spinVarNz

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    ip,ivar,ix,iy,jvar,np,nz        !  work/loop counters

  INTEGER :: spinPoint(npSpinMax)   !  list of points that are not yet fully spun-up.
!                                      At present only for output to screen.

  REAL ::  &!  local ARRAYS
    newVal(nzSpinMax)     &!  current values of spin up fields at a point
   ,toler(nzSpinMax)       !  tolerances

  LOGICAL ::  &!  local SCALARS
    firstCall=.TRUE.  !  F except on first call to this routine

  LOGICAL ::  &!  local ARRAYS
    spinDoneVar(nspinVar,npSpinMax)  !  T when a variable is spun up at this point

!###############################################################################

  IF ( firstCall ) THEN
!-------------------------------------------------------------------------------
!   On first call, process chosen variables, allocate space and save initial values.
!-------------------------------------------------------------------------------
    firstCall = .FALSE.

!-------------------------------------------------------------------------------
!   Set sizes of each chosen variable.
!-------------------------------------------------------------------------------
    DO ivar=1,nspinVarMax
      IF ( spinVar(ivar) ) THEN
        IF ( ivar == iposRouteStore ) THEN
          spinVarNp(ivar) = npRoute
          spinVarNz(ivar) = 1
        ELSEIF ( ivar==iposSmcl .OR. ivar==iposTsoil ) THEN
          spinVarNp(ivar) = land_pts
          spinVarNz(ivar) = sm_levels
        ELSE
          WRITE(*,*)'ERROR: spin_check: no code for name=',TRIM(spinVarName(ivar))
        ENDIF
      ENDIF  !  spinVar
    ENDDO  !  ivar

!-------------------------------------------------------------------------------
!   Save maximum size.
!-------------------------------------------------------------------------------
    npSpinMax = MAXVAL( spinVarNp(:),spinVar(:) )
    nzSpinMax = MAXVAL( spinVarNz(:),spinVar(:) )

!-------------------------------------------------------------------------------
!   Allocate space.
!-------------------------------------------------------------------------------
    CALL allocate_arrays( 'spin_check' )

!-------------------------------------------------------------------------------
!   Save initial values.
!-------------------------------------------------------------------------------
    spinEnd = .FALSE.
    spinValOld(:,:,:) = 0.0
    jvar = 0
    DO ivar=1,nspinVarMax
      IF ( spinVar(ivar) ) THEN
        np = spinVarNp(ivar)
        nz = spinVarNz(ivar)
        jvar = jvar + 1

        IF ( ivar == iposRouteStore ) THEN
          spinValOld(jvar,1:np,1) = RESHAPE( routeStore(:,:), (/ npRoute /) )
        ELSEIF ( ivar == iposSmcl ) THEN
          spinValOld(jvar,1:np,1:nz) = smcl(:,:)

        ELSEIF ( ivar == iposTsoil ) THEN
          spinValOld(jvar,1:np,1:nz) = t_soil(:,:)

        ELSE
          WRITE(*,*)'ERROR: spin_check: no code for name=',TRIM(spinVarName(ivar))
        ENDIF

      ENDIF !  spinVar
    ENDDO  !  ivar

!-------------------------------------------------------------------------------
! At all subsequent calls, check for spin up.
!-------------------------------------------------------------------------------
  ELSE  !  NOT firstCall

!   Reset flag.
    spinDoneVar(:,:) = .FALSE.

!-------------------------------------------------------------------------------
!   Test all variables (even if some have spun up - changes to other variables
!   could plausibly have knock-on effects).
!-------------------------------------------------------------------------------
    jvar = 0
    DO ivar=1,nspinVarMax
      IF ( spinVar(ivar) ) THEN
        np = spinVarNp(ivar)
        nz = spinVarNz(ivar)
        jvar = jvar + 1

!-------------------------------------------------------------------------------
!       Test all points (even if some may previously have spun up for this 
!       variable - changes in other variables may have impacted on this 
!       apparently spun-up variable).
!-------------------------------------------------------------------------------
        DO ip=1,np

!-------------------------------------------------------------------------------
!         Get the current values of this variable at this point.
!-------------------------------------------------------------------------------
          IF ( ivar == iposRouteStore ) THEN
!           Use routeIndex to extract "active" routing points.
            CALL getXYPos( routeIndex(ip),nxRoute,nyRoute,ix,iy )
            newVal(1) = routeStore(ix,iy)

          ELSEIF ( ivar == iposSmcl ) THEN
            newVal(1:nz) = smcl(ip,1:nz)  

          ELSEIF ( ivar == iposTsoil ) THEN
            newVal(1:nz) = t_soil(ip,1:nz)    

          ELSE
            WRITE(*,*)'ERROR: spin_check: no code for name=',TRIM(spinVarName(ivar))
          ENDIF
      
!-------------------------------------------------------------------------------
!         Get max allowed change.
!-------------------------------------------------------------------------------
          IF ( spinTolPercent(ivar) ) THEN
!           Take absolute value (since <0 will always fail the spin up test).
            toler(1:nz) = ABS( spinValOld(jvar,ip,1:nz) ) * spinTol(ivar)/100.0
          ELSE
            toler(1:nz) = spinTol(ivar)
          ENDIF

!-------------------------------------------------------------------------------
!         Test if all layers have passed spin-up.
!-------------------------------------------------------------------------------
          IF ( ALL( ABS( newval(1:nz)-spinValOld(jvar,ip,1:nz) ) <= toler(1:nz) ) )  &
            spinDoneVar(jvar,ip) = .TRUE.

!-------------------------------------------------------------------------------
!         In all cases (spun-up or not), save these values for next use as 'old'
!         values.  Even if this variable has spun up at this point, other 
!         variables may not have fully spun up yet, which in turn could impact 
!         on this apparently spun up variable.
!-------------------------------------------------------------------------------
          spinValOld(jvar,ip,1:nz) = newval(1:nz)

        ENDDO       !  ip (points)

!       Set mask to "done" at all points not needed for this variable.
        spinDoneVar(jvar,spinVarNp(ivar)+1:) = .TRUE.

      ENDIF !  spinVar
    ENDDO   !  ivar

!-------------------------------------------------------------------------------
!   Test if all variables have spun up at all points.
!-------------------------------------------------------------------------------
    IF ( ALL( spinDoneVar(:,:) ) ) THEN
!     Spin up has finished.
      spinEnd = .TRUE.
      spinUp = .FALSE.
      WRITE(*,"(a,i4,a,/,50('#'))")'All points spun up after ',ispin,' cycles.'
!     Save the number of cycles used.
      nspinFinal = ispin
    ELSE
!-------------------------------------------------------------------------------
!     One or more points have not fully spun up yet.
!-------------------------------------------------------------------------------
      IF ( ispin==nspin ) WRITE(*,*)'Model has not fully spun up after max number (=',nspin,') of spin up cycles'
      jvar = 0
      DO ivar=1,nspinVarMax
        IF ( spinVar(ivar) ) THEN
          jvar = jvar + 1
          np = COUNT( .NOT.spinDoneVar(jvar,:) )
          IF ( np > 0 ) THEN
            WRITE(*,*) TRIM(spinVarName(ivar)),' not fully spun up at ',np,' points.'
          ELSE
            WRITE(*,*) TRIM(spinVarName(ivar)),' has fully spun up at all points.'
          ENDIF
          IF ( ispin == nspin ) THEN
!           Get a list of points that are not spun up for this variable.
            np = 0
            DO ip=1,spinVarNp(ivar)
              IF ( .NOT. spinDoneVar(jvar,ip) ) THEN
                np = np + 1
                spinPoint(np) = ip
              ENDIF
            ENDDO
            IF ( np < spinVarNp(ivar) ) THEN
              WRITE(*,"(a,(20i6))")'Points that have not spun-up are:',spinPoint(1:np)
            ELSE
              WRITE(*,"(tr4,a)")'i.e. No points have completely spun up for this variable.'
            ENDIF
          ENDIF   !  ispin=nspin
        ENDIF  !  spinVar
      ENDDO  !  ivar

!-------------------------------------------------------------------------------
!     If model has not fully spun up, then either stop the simulation or report
!     that spin-up is imperfect and continue with simulation.
!-------------------------------------------------------------------------------
      IF ( ispin == nspin ) THEN
        IF ( spinFail ) THEN
!         Terminate run, as gracefully as possible.
          WRITE(*,"(/,a)") 'Model has failed to spin up sufficiently.'
          WRITE(*,*)'This run will be terminated.'
          CALL jules_final( 'spinFail' )
          WRITE(*,*)'Stopping in spinUp_check: spinFail: run terminated'
          STOP
        ELSE
          WRITE(*,"(a,/,50('#'))")'No further spin up will be attempted.'
          spinEnd = .TRUE.
          spinUp = .FALSE.
!         Save the number of cycles used.
          nspinFinal = ispin
        ENDIF
      ENDIF  !  ispin

    ENDIF  !  all( spinDoneVar)

  ENDIF  !  firstCall

END SUBROUTINE Spin_Check
  
