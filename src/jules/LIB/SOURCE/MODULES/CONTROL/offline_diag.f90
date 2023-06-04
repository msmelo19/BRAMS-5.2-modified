! module offline_diag
! Module containing declarations of variables that are used to save diagnostics
! from "physics-level" code. This code is only for use in offline (i.e. not UM)
! JULES. It is kept separate for ease of removing from the full JULES code.

!###############################################################################

  MODULE offline_diag
!###############################################################################

  IMPLICIT NONE

!-------------------------------------------------------------------------------
! Scalar parameters.
!-------------------------------------------------------------------------------
  LOGICAL, PARAMETER ::  &
    offDiag = .FALSE.     !  switch indicating if the version of JULES code being
!                            used allows these "offline diagnostics".
!                            This is FALSE, unless you have added the relevant
!                            code.
! These "offline diagnostics" are disabled in the "official" JULES code, since
! they have to be loaded in the "physics" routines which have to be compatible
! with the UM. All the "offline diagnostic" code is left in the control level of
! this code, but has been removed from the physics/UM code - so the diagnostic
! cannot be set. If you would like to use these diagnostics, contact the JULES
! developers!

!-------------------------------------------------------------------------------
! It is expected that every diagnostic covered by this code will have
! 1) a logical variable: useXXDiag
! 2) space to save the data: XXDiag
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Scalar variables.
!-------------------------------------------------------------------------------
  LOGICAL ::  &
    useCiDiag        &!  switch for ci diagnostic
   ,useGstomDiag     &!  switch for gstom diagnostic
   ,useRdcDiag       &!  switch for rdc diagnostic
   ,useRflowDiag     &!  switch for rflow diagnostic
   ,useRoffInfDiag   &!  switch for roffInf diagnostic
   ,useRrunDiag      &!  switch for runoffR diagnostic
   ,useSnowGMeltDiag &!  switch for snowGMelt diagnostic
   ,useWfluxDiag     &!  switch for wflux diagnostic
   ,useWfluxSfcDiag   !  switch for wfluxSfc diagnostic

!-------------------------------------------------------------------------------
! Array variables.
!-------------------------------------------------------------------------------
  REAL, ALLOCATABLE ::  &
    ciDiag(:,:)          &!  ci diagnostic
   ,gstomDiag(:,:)       &!  gstom diagnostic
   ,rdcDiag(:,:)         &!  rdc diagnostic
   ,rflowDiag(:,:)       &!  rflow diagnostic (river outflow, on routing grid)
   ,roffInfDiag(:)       &!  roffInf diagnostic (infiltration excess runoff)
   ,rrunDiag(:,:)        &!  runoffR diagnostic (runoff on routing grid)
   ,snowGMeltDiag(:,:)   &!  snowGMelt diagnostic (melt of snow beneath canopy)
   ,wfluxDiag(:,:)       &!  wflux diagnostic (soil vertical water flux)
   ,wfluxSfcDiag(:)       !  wfluxSfc diagnostic (vertical water flux at soil surface)

!###############################################################################

  END MODULE offline_diag
!###############################################################################
