! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module with setting of 
! ratio of roughness lengths
  
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.
 
MODULE c_z0h_z0m

  USE max_dimensions, ONLY : ntype_max
  
  IMPLICIT NONE

  REAL, ALLOCATABLE ::                                            &
    z0h_z0m(:)            ! Ratio of roughness length for heat
!                         to roughness length for momentum
!                         for each surface type.

!-----------------------------------------------------------------------
! Create a fixed length counterpart for namelist
!-----------------------------------------------------------------------
  REAL ::                                                         &
    z0h_z0m_io(ntype_max)

  NAMELIST /jules_z0h_z0m/ z0h_z0m_io

END MODULE c_z0h_z0m
