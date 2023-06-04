! Read miscellaneous surrface and carbon/veg parameters.

  SUBROUTINE init_misc()

  USE aero, ONLY : &
!   imported scalars with intent(out)
      co2_mmr

  USE csmin, ONLY :  &
!   imported scalars with intent(out)
      cs_min

  USE file_utils, ONLY:  &
!   imported procedures
      findTag

  USE inout, ONLY :  &
!   imported scalars with intent(in)
      echo,jinUnit

  USE seed, ONLY :  &
!   imported scalars with intent(out)
      frac_min,frac_seed

  USE sigm, ONLY :  &
!   imported scalars with intent(out)
      pow

  USE surf_param, ONLY :  &
!   imported scalars with intent(out)
      beta1,beta2,fwe_c3,fwe_c4,hleaf,hwood,kaps,kaps_roth,q10_leaf,q10_soil

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     routeOnly

  IMPLICIT NONE

!------------------------------------------------------------------------------------------

! If data are not needed, nothing to do.
  IF ( routeOnly ) RETURN

  if (echo) WRITE(*,"(50('-'),/,a)") 'init_misc'

! Locate the start of this section in input file.
  CALL findTag( jinUnit,'init_misc','>INIT_MISC' )

  READ(jinUnit,*) hleaf,hwood
  READ(jinUnit,*) beta1,beta2
  READ(jinUnit,*) fwe_c3,fwe_c4
  READ(jinUnit,*) q10_leaf
  READ(jinUnit,*) kaps
  READ(jinUnit,*) kaps_roth(:)
  READ(jinUnit,*) q10_soil
  READ(jinUnit,*) cs_min
  READ(jinUnit,*) co2_mmr
  READ(jinUnit,*) frac_min, frac_seed
  READ(jinUnit,*) pow

  IF ( echo ) THEN
    WRITE(*,*) 'hleaf=',hleaf,' hwood=',hwood
    WRITE(*,*) 'beta1=',beta1,' beta2=',beta2
    WRITE(*,*) 'fwe_c3=',fwe_c3,' fwe_c4=',fwe_c4
    WRITE(*,*) 'q10_leaf=',q10_leaf
    WRITE(*,*) 'kaps=',kaps
    WRITE(*,*) 'kaps_roth=',kaps_roth(:)
    WRITE(*,*) 'q10_soil=',q10_soil
    WRITE(*,*) 'cs_min=',cs_min
    WRITE(*,*) 'co2_mmr=',co2_mmr
    WRITE(*,*) 'frac_min=',frac_min,' frac_seed=',frac_seed
    WRITE(*,*) 'pow=',pow
  ENDIF

  END SUBROUTINE init_misc
