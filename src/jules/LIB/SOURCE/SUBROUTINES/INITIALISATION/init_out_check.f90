!##########################################################################################
!##########################################################################################
! subroutine init_out_check
! Check that the requested output variables are available - i.e. that they are consistent with
! the model setup. This is particularly important if a FORTRAN variable is only allocated
! if a model option is selected. A variable that is available (e.g. allocated) but that is
! not used by the model setup may be allowed - this will likely become increasingly rare if
! we move to only allocating space for variables that are needed for a given configuration.

  SUBROUTINE init_out_check

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     nsmax

  USE inout, ONLY :  &
!  imported scalars with intent(in)
     echo,nout  &
!  imported arrays with intent(in)
    ,nvarOut,varNameList,varNum,varPos,varType

  USE route_mod, ONLY :  &
!   imported scalar parameters
      routeTypeTRIP  &
!   imported scalars with intent(in)
     ,routeType   

  USE switches , ONLY :  &
!  imported scalars with intent(in)
     can_model,l_spec_albedo,route,routeOnly

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    iout,ivar,j,jvar   !  loop counters/work

  LOGICAL ::  &!  local scalars
    errFlag   !  T means an error has been detected
!-------------------------------------------------------------------------------
  errFlag = .FALSE.

  DO iout=1,nout
    DO ivar=1,nvarOut(iout)
      jvar = varPos(iout,ivar)  !  position in list of selected variables
      j = varNum(jvar)          !  position in list of available variables

!-------------------------------------------------------------------------------
!     First look for incorrect types of variable and raise an error if found.
!     These catch-all cases mean we will detect an error even if a test for a
!     specific variable has not been included below.
!-------------------------------------------------------------------------------

!     If only doing routing, no non-routing variables are allowed.
!     In fact, the input runoff rates should be allowed.
!      IF ( routeOnly .AND. varType(jvar)/='RP' .AND. varType(jvar)/='RG' ) THEN
!        errFlag = .TRUE.
!        WRITE(*,*)'ERROR: init_out_check: only routing variables can be output.'
!        WRITE(*,*)'Non-routing variables are not available because routeOnly=TRUE.'
!        WRITE(*,*)'Incorrect request for: ',TRIM( varNameList(j) )
!      ENDIF

!     If not doing routing, no routing variables are allowed.
      IF ( .NOT. route .AND. ( varType(jvar)=='RP' .OR. varType(jvar)=='RG' ) ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: no routing variables can be output.'
        WRITE(*,*)'Routing variables are not available because route=FALSE.'
        WRITE(*,*)'Incorrect request for: ',TRIM( varNameList(j) )
      ENDIF

!     If not using multi-layer snow model, can't get those diagnostics.
      IF ( .NOT. routeOnly .AND. nsmax<1 .AND. varType(jvar)=='SN' ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: no snow layer  variables can be output,'
        WRITE(*,*)'because the multi-layer snow model is not selected (nsmax<1)'
        WRITE(*,*)'Incorrect request for: ',TRIM( varNameList(j) )
      ENDIF

!-------------------------------------------------------------------------------
!     Now test for specific variables - to give more information.
!     For some, a warning is raised, but not a fatal error.
!     Obviously, variable names here must be the same as used in init_out_varlist!
!     It is possible to output these variables (because variables exist) but
!     the values are meaningless if the appropriate sub-model is not selected.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!     Warnings for snow variables.
!-------------------------------------------------------------------------------

!     can_model=4 variables.
      IF ( can_model /= 4 .AND. echo ) THEN
        IF ( varNameList(j)=='snomltSubHtf' ) THEN
          WRITE(*,*)'WARNING: snomltSubHtf selected for output but can_model/=4'
          WRITE(*,*)'This field is not modelled - will be zero.'
        ELSEIF ( varNameList(j)=='snowGrCan' ) THEN
          WRITE(*,*)'WARNING: snowGrCan selected for output but can_model/=4'
          WRITE(*,*)'This field is not modelled.'
        ELSEIF ( varNameList(j)=='snowGrCanMeltT' ) THEN
          WRITE(*,*)'WARNING: snowGrCanMeltT selected for output but can_model/=4'
          WRITE(*,*)'This field is not modelled.'
        ELSEIF ( varNameList(j)=='snowGrCanT' ) THEN
          WRITE(*,*)'WARNING: snowGrCanT selected for output but can_model/=4'
          WRITE(*,*)'This field is not modelled.'
        ENDIF
      ENDIF   !  can_model and echo

      IF ( ( varNameList(j)=='rgrainT' .OR. varNameList(j)=='rgrainL' )  &
           .AND. .NOT.l_spec_albedo .AND. echo ) THEN
        WRITE(*,*)'WARNING: ',TRIM(varNameList(j)),' selected for output but l_spec_albedo=F.'
        WRITE(*,*)'This field is not modelled.'
      ENDIF

!     Diagnostics that require snow layer variables (but are not layer variables themselves).
      IF ( nsmax < 1 ) THEN
        SELECT CASE ( varNameList(j) )
          CASE ( 'snowIceT', 'snowIceTot', 'snowLiqT', 'snowLiqTot' )
            errFlag = .TRUE.
            WRITE(*,*)'ERROR: init_out_check: variable requires snow layers, but'
            WRITE(*,*)'the multi-layer snow model is not selected (nsmax<1).'
            WRITE(*,*)'variable=',TRIM( varNameList(j) )
        END SELECT
      ENDIF

!-------------------------------------------------------------------------------
!     Warnings for routing variables.
!-------------------------------------------------------------------------------
!     Above, we have tested that routing variables are only available when routing is done.
!     Now we ensure that a variable is only available if the correct routing model is used.
      IF ( varNameList(j)=='rflow' .AND. route .AND. routeType/=routeTypeTRIP ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: rflow selected for output but TRIP routing not selected.'
      ENDIF

      IF ( varNameList(j)=='rInflow' .AND. route .AND. routeType/=routeTypeTRIP ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: rInflow selected for output but TRIP routing not selected.'
      ENDIF

      IF ( varNameList(j)=='rstore' .AND. route .AND. routeType/=routeTypeTRIP ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: rstore selected for output but TRIP routing not selected.'
      ENDIF

      IF ( varNameList(j)=='runoffR' .AND. route .AND. routeType/=routeTypeTRIP ) THEN
        errFlag = .TRUE.
        WRITE(*,*)'ERROR: init_out_check: runoffR selected for output but TRIP routing not selected.'
      ENDIF

!-------------------------------------------------------------------------------
    ENDDO !  ivar
  ENDDO   !  iout

  IF ( errFlag ) THEN
    WRITE(*,*)'ERROR: init_out_check: variable(s) not available - see above.'
    STOP
  ENDIF

  END SUBROUTINE init_out_check

!###############################################################################
!###############################################################################
