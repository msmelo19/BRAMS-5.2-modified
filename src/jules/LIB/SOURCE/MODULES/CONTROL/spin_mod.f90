! module spin_Mod
! Stores variables related to spin up.

  MODULE spin_Mod

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar parameters.
!-------------------------------------------------------------------------------
  INTEGER, PARAMETER :: nspinVarMax = 3  !  max possible number of variables
!                       used to test for spin up (routeStore,smcl,t_soil)

!-------------------------------------------------------------------------------
! Scalar variables.
!-------------------------------------------------------------------------------
  INTEGER ::  &
!  Each possible variables used in spin up must have an associated iposXX variables below
   iposRouteStore  &!  position of routeStore in list of spin up variables
  ,iposSmcl   &!  position of smcl in list of spin up variables
  ,iposTsoil  &!  position of t_soil in list of spin up variables
  ,ispin      &!  the number of the current spin up cycle
  ,npSpinMax  &!  an extent (number of grid points) allocated to spin up variables
  ,nspin      &!  maximum number of cycles of spin up to be performed
  ,nspinFinal &!  the final number of spin up cycles performed
  ,nspinVar   &!  the number of spinup variables selected
  ,nzSpinMax   !  an extent (number of levels) allocated to spin up variables

  LOGICAL ::  &
    spinEnd   &!  T at the start of the first timestep after spin up (i.e. end
!                 of spin up) XX still needed??
   ,spinFail  &!  If nspin<0 (model assesses when spin up complete),
!              !     T means run ends if model has not fully spun-up within
!                    abs(nspin) cycles of spin-up
!              !     F means run continues anyway.
   ,spinUp     !  TRUE while model is spinning up

!-------------------------------------------------------------------------------
! Array variables.
!-------------------------------------------------------------------------------
  INTEGER :: spinVarNp(nspinVarMax)  !  number of space points in each spin up variable
  INTEGER :: spinVarNz(nspinVarMax)  !  number of levels in each spin up variable

  REAL :: spinTol(nspinVarMax) ! maximum change over a cycle of spin-up
!                            that allow model to be considered spun up

  REAL, ALLOCATABLE :: spinValOld(:,:,:)  !  values of spin-up fields at end of previous cycle
!     Space is allocated only for variables that are to be used in determining spin up.


  LOGICAL :: spinTolPercent(nspinVarMax)  !  T if spinTol is a percentage change.
!                                            F if spinTol is to be used as given (not a %)
  LOGICAL :: spinVar(nspinVarMax)         !  T if a variable is used to determine whether
!                                             model has spun-up

!  logical, allocatable :: spinDone(:)  !  T when a point is fully spun up

  CHARACTER(len=11) :: spinVarName(nspinVarMax)  !  names used to identify spin up variables

  END MODULE spin_Mod
