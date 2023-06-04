! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module setting of 
! Variables that specify numbers of surface types.
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE nstypes

  IMPLICIT NONE

!-----------------------------------------------------------------------
! In the UM, these values are set as parameters in nstypes.h
! We don't want parameters for JULES - they are just set to default
! values that are the same as in the UM:

! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Urban
!     7 - Water
!     8 - Soil
!     9 - Ice

!-----------------------------------------------------------------------

! scalars
  INTEGER ::                                                      &
   nnvg  = 4                                                      &
                  !  Number of non-vegetation surface types
  ,npft  = 5                                                      &
                  !  Number of plant functional types
  ,ntype = 9                                                      &
                  !  Number of surface types
  ,urban = 6                                                      &
                  ! Index of the surface type 'Urban'
  ,lake  = 7                                                      &
                  !  Index of the surface type 'Lake'
  ,soil  = 8                                                      &
                  !  Index of the surface type 'Soil'
  ,ice   = 9                                                      &
                  !  Index of the surface type 'Ice'

! For use when l_urban2T == .TRUE.
  ,urban_canyon = 6                                               &
                  ! This is used instead of urban
  ,urban_roof   = 10                                              &
                  ! NB: Not same as in UM (== ice!) as ice tile is hijacked
  ,nurb         = 0
                  ! Number of extra urban tiles types

!-----------------------------------------------------------------------
! Declare a namelist so that we can set the variables from IO
!-----------------------------------------------------------------------
  NAMELIST /jules_nstypes/ npft,nnvg,lake,soil,ice,urban

END MODULE nstypes
