#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine FTSA -------------------------------------------------
!
!    It calculates (true) surface albedos for P234.
!      Release 2.8 of the UM allows for separate surface
!    albedos for direct and diffuse light over sea, calculating the
!    former from the formula of Briegleb and Ramanathan (1982, J. Appl.
!    Met., 21, 1160-1171) and passes the albedos out in different form
!    Land & ice albedos are still the same for direct and diffuse beams.
!                                              William Ingram 25/9/92
!    Suitable for single column model use.
!
!    Author: William Ingram
!
!
!    It conforms to programming standard A of UMDP 4, version 2.
!    It contains ! comments, but otherwise conforms to the FORTRAN 77
!    standard with no features deprecated by 8X.
!
!   Logical components covered : P233
!     (ancillary calculations for the shortwave scheme)
!
!   Project task : P23
!
!    Offline documentation is in UMDP 23, sections "True surface albedo
!    specification" and "Modifications to the radiation scheme to
!    accommodate the leads model"
!
!
SUBROUTINE ftsa (                                                 &
  land, flandg, aice, tstar, tstar_sice, cosz, s, s_sea,          &
  alpham,alphac,alphab,dtice,l_moses_ii,l_ssice_albedo,           &
  l_mod_barker_albedo,                                            &
  l_sice_meltponds,l_sice_scattering,l_sice_hadgem1a,             &
  dt_bare,dalb_bare_wet,pen_rad_frac,beta,version,                &
  l1, l2, sa_land, sa_sice, saos)

USE c_0_dg_c

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

CHARACTER                                                         &
                                    !, INTENT(IN) ::
     version*3                      ! Version of radiation scheme
     
INTEGER                                                           &
                                    !, INTENT(IN) ::
     l1,                                                          &
                                    ! Full field dimension
     l2                             ! Number of points to treat
     
LOGICAL                                                           &
                                    !, INTENT(IN) ::
     land(l1)                                                     &
                                    ! Land-sea mask (land .TRUE.)
    ,l_moses_ii                                                   &
                                    ! .TRUE. if MOSES II land
!                                     surface is selected.
!         Switch on the effect of snow on sea-ice albedo
    ,l_ssice_albedo                                               &
    ,l_mod_barker_albedo                                          &
                                    ! Use modified Barker
!                                     albedo (open sea).
    ,l_sice_meltponds                                             &
                       ! switch on seaice albedo meltponds
    ,l_sice_scattering                                            &
                       ! switch on seaice albedo internal scatter
    ,l_sice_hadgem1a   ! switch for HadGEM1 bug correction
!                        NOT USED IN JULES
!                        Remains to maintain a consistent
!                        subroutine interface for the UM

REAL                                                              &
                                    !, INTENT(IN) ::
     flandg(l1),                                                  &
                                    ! Land fraction
     aice(l1),                                                    &
                                    ! Sea-ice fraction
     tstar(l1),                                                   &
                                    ! Surface temperature
     tstar_sice(l1),                                              &
                                    ! Seaice surface temperature
     cosz(l1),                                                    &
                                    ! cos(solar zenith angle)
     s(l1)                                                        &
                                    ! Snow amount (mass/area)
    ,s_sea(l1)                                                    &
                       ! Snow amount on sea ice (mass/area of ice)
!                       
!     Constants used to determine the albedo of sea-ice:
! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
!
    ,alpham                                                       &
                       ! "M" for "melting"
! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
!
    ,alphac                                                       &
                       ! "C" for "cold"
! Albedo of snow-free sea-ice if l_ssice_albedo
    ,alphab                                                       &
                       ! "B" for "bare"
!
! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
    ,dtice                                                        &
! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
    ,dt_bare                                                      &
! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
    ,dalb_bare_wet                                                &
! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
    ,pen_rad_frac                                                 &
! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
   ,beta

REAL                                                              &
     !, INTENT(OUT)
     sa_land(l1),                                                 &
                                    ! Surface Albedos for Land.
!                                     (Not output for MOSESII).
     sa_sice(l1),                                                 &
                                    ! Surface Albedos for seaice
     saos(l1,2)                     ! and Ice, and for Open Sea,
!    respectively, with zeroes for safety where no value applies
!    FTSA has no dynamically allocated workspace, no EXTERNAL calls
!    and no significant structure - just one loop and some IF blocks


INTEGER j                           ! Loops over points
REAL dsa                            ! Deep-snow albedo (alphasubD)
REAL bx                                                           &
                                    ! 1 - COSZ
! Temperature at which (snow on) sea-ice reaches its "cold" value
    ,tcice                                                        &
! Slope and intercept of temperature dependence of the albedo of
! (snow on) sea-ice
    ,ice1,ice2

!     Local parameters
REAL dtland, kland, tcland, adifc, fcatm,                         &
     snow_albedo,                                                 &
                                    ! Snow albedo
     maskd                          ! Masking depth (S in 3.6.1)
!     Note that the same masking depth is always used, both for land,
!     regardless of vegetation cover, and for sea-ice. This assumption
!     of constancy may be doubtful.
PARAMETER ( maskd = 0.2 )

!     Basic quantities for land CSSA calculations:
PARAMETER ( dtland = 2.,                                          &
                                    ! delta(T) in 3.6.2
     fcatm = 0.3 )                  ! Fraction by which deep-snow
!     albedo changes (from "cold" value towards snow-free value) at TM
!     From these, 2 constants precalculated for efficiency in 3.6.2:
PARAMETER ( kland = 0.3/dtland,                                   &
     tcland = tm-dtland )

PARAMETER ( adifc = 0.06 )          ! Surface albedo of ice-free
!                                     sea for the diffuse beam
!     derive 3 constants from the basic quantities (supplied in the
!     namelist RUNCNST) for sea-ice CSSA calculations:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('FTSA',zhook_in,zhook_handle)
tcice = tm - dtice
ice1 = (alpham-alphac) / dtice
ice2 = alpham - tm*ice1

DO j=1, l2
  sa_land(j)=0.0
  sa_sice(j)=0.0
END DO

DO j=1, l2
  IF (flandg(j) <  1.0) THEN
    IF ((version /= '03C').AND.(version /= '03Z').AND.            &
                               (version /= '03A') ) THEN
       saos(j,1) = 0.05 / ( 1.1 * cosz(j)**1.4 + 0.15 )
    ELSE IF (l_mod_barker_albedo) THEN
       bx=1.0-cosz(j)
       saos(j, 1)= 0.0315 + 0.0421*bx**2                          &
             + 0.128*bx**3 - 0.04*bx**4                           &
             + (3.12/(5.68) + (0.074*bx/(1.0)))*bx**5
    ELSE
       saos(j, 1)=0.026/(cosz(j)**1.7+0.065)                      &
          +0.15*(cosz(j)-0.1)*(cosz(j)-0.5)*(cosz(j)-1.0)
    END IF
    saos(j,2) = adifc
!         Note that the following will add in ICE1*(TSTAR-TFS) to CSSA
!         if AICE#0 when it should be - even if only very small: for
!         large enough TSTAR this will give very large surface heating
!         and even negative atmospheric heating.  Check whether this
!         could occur.
    IF ( aice(j)  ==  0. ) THEN
       sa_sice(j) = 0.
    ELSE
       IF (l_ssice_albedo) THEN
         IF (s_sea(j) >  0.0) THEN   ! snow on sea ice
                                     ! Cox et al., Clim Dyn,1999
           IF (tstar_sice(j) >  tcice) THEN
             snow_albedo=ice2+ice1*tstar_sice(j)
           ELSE
             snow_albedo=alphac
           END IF
           sa_sice(j)=alphab                                      &
           +(snow_albedo-alphab)*(1.0-EXP(-maskd*s_sea(j)))
         ELSE           ! no snow so bare ice only
           IF(l_sice_meltponds) THEN
              ! bare ice, temperature dep. (Ebert &Curry,1993)
              IF (tstar_sice(j) > (tm - dt_bare))THEN
                sa_sice(j) = alphab+(tstar_sice(j)-tm+dt_bare)    &
                                   *dalb_bare_wet
              ELSE      ! surface is dry
                sa_sice(j)=alphab
              END IF     ! end melt ponds
           ELSE         ! just use bare ice albedo
             sa_sice(j)=alphab
           END IF       ! l_sice_meltponds
           IF(l_sice_scattering) THEN
               !Semtner modification dry albedo for internal
               ! scattering (Semnter, 1976)
               sa_sice(j)=sa_sice(j)+beta*(1.0-sa_sice(j))        &
                              *pen_rad_frac
           END IF       ! l_sice_scattering
         END IF         !  any snow on ice
       ELSE             ! default to operational NWP scheme
                        !  3.5.1:
       IF ( tstar_sice(j)  <   tcice ) THEN
          sa_sice(j) = alphac
        ELSE
          sa_sice(j) = ice1 * tstar_sice(j) + ice2
       END IF
       END IF
    END IF
 END IF

#if defined(A01_3C)||defined(A01_3Z)
  IF (sa_sice(j) <  0.0) THEN

     sa_sice(j)=alpham

  END IF
#endif
END DO


IF (lhook) CALL dr_hook('FTSA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ftsa
#endif
