  MODULE imogen_run
  
    IMPLICIT NONE
    
    REAL ::                                                       &
      CO2_INIT_PPMV=286.03
                  ! Initial CO2 concentration (ppmv)
        
    CHARACTER(len=180) ::                                         &
      FILE_SCEN_EMITS = 'emits_HADCM3.dat',                       &
                  ! If used, file containing CO2 emissions in G
      FILE_NON_CO2_VALS='',                                       &
                  ! If used, file containing non-CO2 values
      FILE_SCEN_CO2_PPMV = 'co2_vals.dat'
                  ! If used, file containing CO2 values
                  
    LOGICAL ::                                                    &
      ANLG = .TRUE.,                                              &
                  ! If true, then use the GCM analogue model
      ANOM = .TRUE.,                                              &
                  ! If true, then use the GCM analogue model
      C_EMISSIONS = .TRUE.,                                       &
                  ! If true, means CO2 concentration is calcula
      INCLUDE_CO2 = .TRUE.,                                       &
                  ! Are adjustments to CO2 values allowed?
      INCLUDE_NON_CO2 = .TRUE.,                                   &
                  ! Are adjustments to non-CO2 values allowed? 
      LAND_FEED = .FALSE.,                                        &
                  ! Are land feedbacks allowed on atmospheric C
      OCEAN_FEED = .FALSE.,                                       &
                  ! Are ocean feedbacks allowed on atmospheric 
      WGEN = .FALSE.
                  ! Is the weather generator switched on.
                  
    INTEGER ::                                                    &
      NYR_EMISS = 241
                  ! Number of years of emission data in file.
      
    NAMELIST /RUN/ CO2_INIT_PPMV,FILE_SCEN_EMITS,                 &
                   FILE_SCEN_CO2_PPMV,NYR_EMISS,C_EMISSIONS,      &
                   INCLUDE_CO2,INCLUDE_NON_CO2,LAND_FEED,         &
                   OCEAN_FEED,WGEN,ANOM,ANLG,FILE_NON_CO2_VALS
  
  
  
  END MODULE imogen_run
