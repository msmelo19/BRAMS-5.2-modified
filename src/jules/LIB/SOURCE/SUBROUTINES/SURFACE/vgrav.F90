#if defined(L08_1A)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine VGRAV -------------------------------------------------

! Purpose: To calculate the gravitational sedimentation velocity of
!          tracer particles according to Stoke's law, including the
!          Cunningham correction factor.

!Documentation: Ref. Pruppacher & Klett
!                    Microphysics of clouds & ppn    1978,1980 edns.

!-----------------------------------------------------------------------

SUBROUTINE vgrav(                                                 &
 row_length,rows,nlevs,diam,rhop,p,t,                             &
 vstokes,ccf,etaa)

USE dust_param, ONLY: accf, bccf, cccf
USE c_sulchm
USE c_0_dg_c
USE c_g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE


INTEGER row_length         !IN row length
INTEGER rows               !IN number of rows
INTEGER nlevs              !IN number of levels

REAL diam                  !IN particle diameter
REAL rhop                  !IN particles density
REAL p(row_length,rows,nlevs)!IN pressure
REAL t(row_length,rows,nlevs)!IN temperature

REAL vstokes(row_length,rows,nlevs) !OUT sedimentation velocity
REAL etaa(row_length,rows,nlevs)!OUT viscosity of air
REAL ccf(row_length,rows,nlevs) !OUT cunningham correction factor


! local variables

INTEGER ilev               !LOC loop counter for levels
INTEGER i                  !LOC loop counter
INTEGER j                  !LOC loop counter

REAL tc(row_length,rows)   !LOC temperature in deg C
REAL lamdaa(row_length,rows)!LOC mean free path of particle
REAL alphaccf(row_length,rows)!LOC

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Calculate viscosity of air (Pruppacher & Klett p.323)
IF (lhook) CALL dr_hook('VGRAV',zhook_in,zhook_handle)
DO ilev=1,nlevs
  DO j=1,rows
    DO i = 1,row_length
     tc(i,j)=t(i,j,ilev)-zerodegc
     IF (tc(i,j)  >=  0.) THEN
       etaa(i,j,ilev)=(1.718+0.0049*tc(i,j))*1.e-5
     ELSE
       etaa(i,j,ilev)=                                             &
            (1.718+0.0049*tc(i,j)-1.2e-5*tc(i,j)*tc(i,j))*1.e-5
     END IF

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NLEVS

DO ilev=1,nlevs
  DO j=1,rows
    DO i=1,row_length

! Calculate mean free path of particle (Pruppacher & Klett p.323)
      lamdaa(i,j)=mfp_ref*pref_mfp*t(i,j,ilev)/                    &
      (p(i,j,ilev)*tref_mfp)
! Calculate Cunningham correction factor(Pruppacher & Klett p.361)
      alphaccf(i,j)=accf+bccf*EXP(cccf*diam*.5/lamdaa(i,j))
      ccf(i,j,ilev)=(1.+alphaccf(i,j)*lamdaa(i,j)/(.5*diam))
! Calculate sedimentation velocity (Pruppacher & Klett p.362)
      vstokes(i,j,ilev)=ccf(i,j,ilev)*(diam*diam*g*rhop)/          &
             (18.*etaa(i,j,ilev))

    END DO !ROW_LENGTH
  END DO !ROWS
END DO !NLEV

IF (lhook) CALL dr_hook('VGRAV',zhook_out,zhook_handle)
RETURN
END SUBROUTINE vgrav
#endif
