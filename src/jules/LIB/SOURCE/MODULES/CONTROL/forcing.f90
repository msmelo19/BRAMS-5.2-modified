! Module containing all of the driving (atmospheric forcing) variables.

  module forcing

!-------------------------------------------------------------------------------

! The forcing variables.
  real, allocatable ::  &!
    QW_1(:,:)          &!  Total water content (Kg/Kg)
   ,TL_1(:,:)          &! Ice/liquid water temperature (k)
   ,U_0(:,:)           &! W'ly component of surface current (m/s)
   ,V_0(:,:)           &! S'ly component of surface current (m/s)
   ,U_1(:,:)           &! W'ly wind component (m/s)
   ,V_1(:,:)           &! S'ly wind component (m/s)
   ,PSTAR(:,:)         &! Surface pressure (Pascals)
   ,LS_RAIN(:,:)       &! Large-scale rain (kg/m2/s)
   ,CON_RAIN(:,:)      &! Convective rain (kg/m2/s)
   ,LS_SNOW(:,:)       &! Large-scale snowfall (kg/m2/s)
   ,CON_SNOW(:,:)      &! Convective snowfall (kg/m2/s)
   ,SW_DOWN(:,:)       &! Surface downward SW radiation (W/m2)
   ,LW_DOWN(:,:)        ! Surface downward LW radiation (W/m2)

!-------------------------------------------------------------------------------

  end module forcing
