!###############################################################################
!###############################################################################
! program jules
! Main program unit for JULES.

PROGRAM jules

  USE inout, ONLY :  &
!  imported scalars with intent(out)
    print_step

  USE output_mod, ONLY :  &
!  imported procedures
    output

  USE spin_mod, ONLY :  &
!  imported scalars with intent(inout)
    ispin,spinUp

  USE time_loc, ONLY :   &
!  imported scalars with intent(out)
    date,time,newYear,endYear

  USE time_mod, ONLY :   &
!  imported procedures
    s_to_chhmmss

  USE trifctl, ONLY :  &
!  imported scalars with intent(out)
    asteps_since_triffid

  USE update_mod, ONLY :  &
!  imported procedures
   drive_update,veg_update

  USE veg_io_vars, ONLY :  &
!  imported scalars with intent(out)
    vegVaryT
    
  USE switches, ONLY : l_imogen

!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ::   &!  local SCALARS
    a_step      !  timestep counter

  LOGICAL ::  &!  local SCALARS
    endRun    !  TRUE at end of last timestep in run

  CHARACTER(len=10) ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
!                                 ! that is required because some compilers can not cope
!                                 ! with recursive write statements

print*, "O fonte jules.f90 nao eh utilizado no modelo acoplado"
STOP
!-------------------------------------------------------------------------------
! Call the initialisation routine.
!-------------------------------------------------------------------------------
  CALL init

!-------------------------------------------------------------------------------
! Loop over timesteps.
! Note that the number of timesteps is of unknown length at the start of run,
! if the model is to determine when it has spun up.
!-------------------------------------------------------------------------------

  a_step = 0
  DO    !  timestep

!-------------------------------------------------------------------------------
!   Generate output if this is required at the start of a time period (i.e.
!   rarely).  The logical argument (endCall) is always false at this call.
!-------------------------------------------------------------------------------
    CALL output( a_step,.FALSE. )

!-------------------------------------------------------------------------------
!   Increment timestep counters.
!-------------------------------------------------------------------------------
    a_step = a_step + 1
    ASTEPS_SINCE_TRIFFID=ASTEPS_SINCE_TRIFFID+1

    IF ( MOD(a_step,print_step)==0 .OR. a_step==1 ) THEN
      time_hms = s_to_chhmmss( time )
      IF ( spinUp ) THEN
        WRITE(*,"(a,i7,tr2,a,i8,tr1,a,a,i3)") 'timestep:',a_step  &
                   ,'Start date and time: '  &
                    ,date,time_hms,' Spin up cycle: ',ispin
      ELSE
        WRITE(*,"(a,i7,tr2,a,i8,tr1,a)") 'timestep:',a_step  &
                  ,'Start date and time: ',date,time_hms
      ENDIF
    ENDIF
   
!-------------------------------------------------------------------------------
!   Update the IMOGEN climate if required
!-------------------------------------------------------------------------------
    IF( l_imogen .AND. newYear ) CALL imogen_update_clim

!-------------------------------------------------------------------------------
!   Update meteorological data.  2nd argument (next) is TRUE to indicate that
!   the "next" data in file are to be used.
!-------------------------------------------------------------------------------
    !DSM CALL drive_update( a_step,.TRUE. )
    print*, "Desativada a chamada da subrotina drive_update em jules.f90"
    STOP

!-------------------------------------------------------------------------------
!   Update prescribed vegetation fields.  2nd argument (next) is TRUE to
!   indicate that the "next" data in file are to be used.
!-------------------------------------------------------------------------------
    IF ( vegVaryT ) CALL veg_update( a_step,.TRUE. )

!-------------------------------------------------------------------------------
!   Call the main model routine.
!-------------------------------------------------------------------------------
    CALL control (a_step)
    
!-------------------------------------------------------------------------------
!   Update IMOGEN carbon if required
!-------------------------------------------------------------------------------
    IF( l_imogen .AND. endYear ) CALL imogen_update_carb

!-------------------------------------------------------------------------------
!   Generate output. This call (at the end of the timestep) is expected to
!   generate most of the output for a run.  The logical argument (endCall) is
!   always true at this call.
!-------------------------------------------------------------------------------
    CALL output( a_step,.TRUE. )

!-------------------------------------------------------------------------------
!   Update time.
!-------------------------------------------------------------------------------
    CALL newTime( a_step,endRun, 0.,0.)

!-------------------------------------------------------------------------------
!   Check for end of run.
!-------------------------------------------------------------------------------
    IF ( endRun ) EXIT

  ENDDO  !  timestep loop

!-------------------------------------------------------------------------------
! End of model run.
!-------------------------------------------------------------------------------
  time_hms = s_to_chhmmss( time )
  WRITE(*,"(/,a,i7,a,i8,tr1,a,/)") 'End of model run after ',a_step  &
      ,' timesteps. Date and time: ',date,time_hms

!-----------------------------------------------------------------------
! Call routine to dump final state, deallocate memory etc.
!-----------------------------------------------------------------------
  CALL jules_final( 'final' )

  WRITE(*,"(/,a)")'End.'

END PROGRAM jules

!###############################################################################
!###############################################################################
