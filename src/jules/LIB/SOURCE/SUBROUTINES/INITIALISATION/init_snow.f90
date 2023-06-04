! Read snow parameters.

  SUBROUTINE init_snow()

  USE c_0_dg_c, ONLY :  &!
!   imported scalar parameters
     tm

  USE file_utils, ONLY :  &!
!   imported procedures
     findTag

  USE inout, ONLY :  &!
!   imported scalars with intent(in)
     echo,jinUnit

  USE rad_param, ONLY :  &!
!   imported scalars with intent(out)
     amax,dtland,kland,maskd,r0,rmax,tcland  &!
!   imported scalars with intent(out)
    ,snow_ggr

  USE ancil_info, ONLY :   &!
!   imported scalar with intent(out)
     nsmax

  USE snow_param, ONLY :   &!
!   imported scalars with intent(out)
     rho_snow_const,snow_hcap,snow_hcon,snowInterceptFact, &!
     snowLoadLAI,snowUnloadFact,snowliqcap,      &
!   imported array with intent(out)
     dzsnow

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     routeOnly
!-------------------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::  &!  local SCALARS
    inUnit     !  unit used to connect to input file

!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_snow'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_snow','>INIT_SNOW' )

  inUnit = jinUnit

! Read the data.
  READ(inUnit,*) dzsnow(1:nsmax)
  READ(inUnit,*) rho_snow_const
  READ(inUnit,*) snow_hcap,snow_hcon
  READ(inUnit,*) snowliqcap
  READ(inUnit,*) r0,rmax
  READ(inUnit,*) snow_ggr(:)
  READ(inUnit,*) amax(1:2)
  READ(inUnit,*) dtland,kland  !  kland is NOT expected to have dtland in denominator
  READ(inUnit,*) maskd
! Values used for can_model=4.
  READ(inUnit,*) snowLoadLAI,snowInterceptFact,snowUnloadFact

!------------------------------------------------------------------------------------------

! MOSES code divides kland by dtland. Note this blows up if dtland=0.
  kland = kland / dtland
  tcland = tm - dtland

!------------------------------------------------------------------------------

! Write to screen.
  IF ( echo ) THEN
    WRITE(*,*) 'dzsnow=',dzsnow
    WRITE(*,*) 'rho_snow_const=',rho_snow_const
    WRITE(*,*) 'snow_hcap=',snow_hcap,' snow_hcon=',snow_hcon
    WRITE(*,*) 'snowliqcap=',snowliqcap
    WRITE(*,*) 'r0=',r0,' rmax=',rmax
    WRITE(*,*) 'snow_ggr=',snow_ggr(:)
    WRITE(*,*) 'amax=',amax(:)
    WRITE(*,*) 'dtland,kland=',dtland,kland
    WRITE(*,*) 'maskd=',maskd
    WRITE(*,*) 'snowLoadLAI=',snowLoadLAI,' snowInterceptFact=',snowInterceptFact
    WRITE(*,*) 'snowUnloadFact=',snowUnloadFact
  ENDIF

  END SUBROUTINE init_snow
