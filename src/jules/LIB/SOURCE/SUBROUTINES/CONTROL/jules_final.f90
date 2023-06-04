!#################################################################################
!#################################################################################
!
! subroutine jules_final
! Routine that is called at the end of a run.
! Dumps the final state, deallocates memory etc.
!
!#################################################################################
!#################################################################################
SUBROUTINE jules_final( callType )

  USE file_utils, ONLY :  &
!  imported procedures
     closeFile

  USE initial_mod, ONLY :  &
!  imported procedures
     dump_io

  USE inout, ONLY :  &
!  imported scalar parameters
     formatNc  &
!  imported scalars with intent(in)
    ,dumpFreq,nout,outFormat  &
!  imported arrays with intent(in)
    ,outActivePrev,outUnit  &
!  imported scalars with intent(inout)
    ,dumpWarn  &
!  imported scalars with intent(out)
    ,nout,outWarnCtl,outWarnEarly,outWarnUnder,redoTdef

  USE output_mod, ONLY : &
!  imported procedures
     rewrite_ctl
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local scalars
    iout   !  loop counter

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callType   !  flag indicating what type of call this is

!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Check that callType is recognised.
!-----------------------------------------------------------------------
  SELECT CASE ( callType )
    CASE ( 'final','spinFail' )   !  nothing to do
    CASE default
      WRITE(*,*)'WARNING: jules_final: do not recognise callType=',TRIM(callType)
!     Proceed anyway!
  END SELECT

!-----------------------------------------------------------------------
! Check that final GrADS ctl files have correct tdef.
!-----------------------------------------------------------------------
  IF ( redoTdef ) THEN   !  only T for GrADS output
    DO iout=1,nout
      CALL rewrite_ctl( .TRUE., iout )
    ENDDO
  ENDIF

!-----------------------------------------------------------------------
! Close any open netCDF output files.
! outActivePrev will be TRUE if an output stream was active on the last timestep.
!-----------------------------------------------------------------------
  IF ( outFormat == formatNc ) THEN
    DO iout=1,nout
      IF ( outActivePrev(iout) ) CALL closeFile( outUnit(iout),outFormat )
    ENDDO
  ENDIF

!-----------------------------------------------------------------------
! Write the final dump, unless no dumps are requested.
!-----------------------------------------------------------------------
!DSM <BRAMS05 nao possui history>    IF ( dumpFreq /= 0 ) THEN
!DSM <BRAMS05 nao possui history>      WRITE(*,*)'Writing dump...'
!DSM <BRAMS05 nao possui history>      SELECT CASE ( callType )
!DSM <BRAMS05 nao possui history>        CASE ( 'final' )
!DSM <BRAMS05 nao possui history>          CALL dump_io( .TRUE., dumpTypeArg='final' )
!DSM <BRAMS05 nao possui history>        CASE ( 'spinFail' )
!DSM <BRAMS05 nao possui history>          CALL dump_io( .TRUE., dumpTypeArg='spinFail' )
!DSM <BRAMS05 nao possui history>      END SELECT
!DSM <BRAMS05 nao possui history>    ENDIF

!-----------------------------------------------------------------------
! Warnings.
!-----------------------------------------------------------------------
  IF ( outWarnEarly ) THEN
    WRITE(*,"(/,70('#'),/,a)")'WARNING: output was generated early or not at all the end of one or more sections.'
    WRITE(*,"(a,/,70('#'))")'See above for details - look for outWarnEarly.'
  ENDIF
  IF ( outWarnUnder ) THEN
    WRITE(*,"(/,70('#'),/,a)")'WARNING: one or more time average or accumulation was set to undef (missing data)'
    WRITE(*,"(a)") 'value because insufficient times were present.'
    WRITE(*,"(a)")'See above for details - look for outWarnUnder.'
    IF ( outWarnEarly ) WRITE(*,"(a)") 'This may be related to outWarnEarly raised above.'
    WRITE(*,"(70('#'))")
  ENDIF
  IF ( outWarnCtl ) THEN
    WRITE(*,"(/,70('#'),/,a)")'WARNING: error raised by rewrite_ctl - possibly could not find file.'
    WRITE(*,"(a,/,70('#'))")'See above for details.'
  ENDIF
  IF ( dumpWarn ) THEN
    WRITE(*,"(/,70('#'),/,a)")'WARNING: a dump file could not be created with the desired name.'
    WRITE(*,*) 'A dump file may have been created with a "non-standard" name.'
    WRITE(*,"(a,/,70('#'))")'See above for details.'
  ENDIF

!-----------------------------------------------------------------------
! Deallocate array space
!-----------------------------------------------------------------------
  CALL deallocate_arrays( 'end' )

!-----------------------------------------------------------------------
! More warnings!
!-----------------------------------------------------------------------
  SELECT CASE ( callType )
    CASE ( 'spinFail' )
      WRITE(*,"(/,a)")'WARNING: jules_final: end of model run because of failure to spin up.'
      WRITE(*,*)'Output might not have been finished as would otherwise have been expected.'
  END SELECT

END SUBROUTINE jules_final
!###############################################################################
!###############################################################################
