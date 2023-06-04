!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

program main

  use ModMemory, only: &
       CreateMemory,   & ! subroutine
       DumpMemory,     & ! subroutine
       DestroyMemory     ! subroutine

  use ModTimeStamp, only: &
       CreateTimeStamp,    & ! subroutine
       DestroyTimeStamp      ! subroutine

  use ParLib, only: &
       parf_init_mpi, & ! subroutine
       parf_barrier, &  ! subroutine
       parf_exit_mpi    ! subroutine

  use ModOneProc, only: &
       OneProc

  ! Main:
  !   starts MPI
  !   initializes memory use and execution time instrumentation
  !   dispatches processes
  !   invoking master/slave processes or full model process
  !   dumps and destroys instrumentation
  !   finishes MPI

  implicit none

  ! process rank and size (local variables)

  integer :: nmachs_in
  integer :: mchnum_in
  integer :: master_num_in
  character(len=*), parameter :: h="**(main)**"
  character(len=8) :: c0, c1
  
  ! execution time instrumentation

  include "tsNames.h"

  integer :: ierr

  ! enroll MPI; get number of processes, process Id, master Id

  call parf_init_mpi(mchnum_in, nmachs_in, master_num_in)

  ! initialize memory instrumentation

  if (mchnum_in==master_num_in) then
     call CreateMemory(mchnum_in, nmachs_in, master_num_in)
  end if

  ! initialize execution time instrumentation

  call CreateTimeStamp(nmachs_in, mchnum_in, names)

  ! dispatch processes

  call OneProc(nmachs_in, mchnum_in, master_num_in)

  ! finishes execution

  call parf_barrier(0)

  call DestroyTimeStamp()

  if (mchnum_in==master_num_in) then
     call DumpMemory("Fim")
     call DestroyMemory()
  end if

  if (mchnum_in==master_num_in) write(*,"(a)") " ****** BRAMS execution ends ******"

  call parf_exit_mpi()
end program main
