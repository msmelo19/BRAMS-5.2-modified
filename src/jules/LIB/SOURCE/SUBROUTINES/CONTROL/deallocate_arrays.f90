!###############################################################################
!###############################################################################
! subroutine deallocate_arrays
! Subroutine to deallocate arrays in modules.
!###############################################################################
!###############################################################################
SUBROUTINE Deallocate_Arrays( callType )

  USE Aero

  USE ancil_info, ONLY :  &
!  imported arrays with intent(in)
    frac,ice_fract,ice_fract_ncat,land_index,land_mask,lice_index,soil_index  &
   ,tile_index,tile_pts,z1_uv,z1_tq,sice_pts_ncat,ssi_index,sea_index  &
   ,sice_index,sice_index_ncat,fssi,sea_frac,sice_frac,sice_frac_ncat

  USE c_elevate, ONLY :  &
!  imported arrays
     surf_hgt,z_land

  USE c_z0h_z0m
  USE Coastal

  USE drive_io_vars, ONLY :  &
!  imported arrays with intent(in)
     driveData,driveDataIn,driveFileDate,driveFileName,driveFileTime  &
    ,driveFileOnUnit,driveVarNameUnit,driveUnit

  USE Fluxes
  USE Forcing

  USE inout, ONLY :  &
!  imported arrays with intent(in)
     coord,coordList,coordLL,havePFT,haveSCpool,haveSnow,haveSoil,haveTile  &
    ,haveType,irecPrevOut,mapIn,mapInLand,mapOut,mapOutCompress,mapOutLand &
    ,nlevMax,nlevMaxCtl,ntCtl,ntCtlNeed  &
    ,ntOutFilePer,ntOutPer,nvarOut,nxyMax,openedFilename,outActivePrev,outAreaLL,outCtlFile  &
    ,outDataFile,outName,outNpWrite,compressGridFile,outGridDxy,outGridNxy  &
    ,outGridXy,pointsFlag,pointsOut,pointsOutLand   &
    ,outDate,outDateFlag,outEndPos,outFilePer,outFileStep,outFirstActive   &
    ,outFirstSection,outFirstWrite,outLLorder,outPer,outPerNunits,outRangeX,outRangeY  &
    ,outSamPer,outStep,outStepSamPer,outTime,outTimeID,outUnit   &
    ,outval,outVarID,outWriteCount,pftUse,rgProfile,rpProfile,snapProfile,taccumVar &
    ,tmeanVar,outTemplate,tmeanProfile,outCompress,useCompressGrid,varDesc  &
    ,varDescList,varName,varNameList    &
    ,varNlev,varNum,varPos,varStartPos,varType,varTypeList,varUnitsList

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error

  USE nvegparm

  USE offline_diag, ONLY :  &
!  imported arrays
     ciDiag,gstomDiag,rdcDiag,rflowDiag,roffInfDiag  &
    ,rrunDiag,snowGMeltDiag,wfluxDiag,wfluxSfcDiag

  USE Orog
  USE P_s_parms
  USE pftparm
  USE Prognostics

  USE route_mod, ONLY :  &
!  imported arrays
     mapInRoute,routeIndex,routeMask,routeNext,routeOrder,roffAccumLand

  USE Screen

  USE snow_param, ONLY :  &
!  imported arrays
     canSnowTile,ds,dzsnow

  USE surf_param, ONLY : diff_frac

  USE soil_param, ONLY : dzsoil

  USE spin_mod, ONLY :  &
!  imported arrays with intent(in)
     spinValOld

  USE time_loc, ONLY :  &
!  imported arrays with intent(in)
     latitude,longitude

  USE Top_pdm, ONLY :  &
!  imported arrays with intent(inout)
     a_fsat,a_fwet,c_fsat,c_fwet,drain,dun_roff,fch4_wetl,fexp,fsat,fwetl  &
    ,gamtot,qbase,qbase_zw,sthzw,ti_mean,ti_sig,zw,INLANDOUT_ATM

  USE trif
  USE Trifctl
  USE U_v_grid

  USE veg_io_vars, ONLY :  &
!  imported arrays with intent(in)
     vegDataIn,vegFileDate,vegFileName,vegFileTime,vegUnit,vegUnitFile  &
    ,vegVarFlag,vegVarInterp,vegVarName,vegvarNameFile,vegVarPos,vegVarStash

  USE urban_param, ONLY :                                                     &
     hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr
     
  USE ozone_vars, ONLY : o3, flux_o3_ft, fo3_ft
  
  USE imogen_clim, ONLY :                                         &
    T_CLIM,RAINFALL_CLIM,SNOWFALL_CLIM,RH15M_CLIM,UWIND_CLIM,     &
    VWIND_CLIM,DTEMP_CLIM,PSTAR_HA_CLIM,SW_CLIM,LW_CLIM,          &
    F_WET_CLIM,LAT,LONG,DCTOT
    
  USE imogen_drive_vars, ONLY :                                   &
    T_OUT,CONV_RAIN_OUT,CONV_SNOW_OUT,LS_RAIN_OUT,LS_SNOW_OUT,    &
    QHUM_OUT,WIND_OUT,PSTAR_OUT,SW_OUT,LW_OUT
    
  USE jules_mod, ONLY : snowdep_surf

!-------------------------------------------------------------------------------
  IMPLICIT NONE

! Scalar arguments with intent(in)
  CHARACTER(len=*), INTENT(in) ::  &
    callType    !  idenitifier for the call to this subroutine

! Local scalars
  INTEGER ::  &
    ierr,ierrSum    !  error codes

!-------------------------------------------------------------------------------

  ierrSum = 0

!-------------------------------------------------------------------------------
! First, specific calls that are intended to deallocate variables that are no
! longer needed, then the "main" call to deallocate everything at the end of the
! run.  At present, the main call does not need to deallocate or test variables
! covered by "specific" calls.
!-------------------------------------------------------------------------------

  SELECT CASE ( callType )

    CASE ( 'init_end' )
!-------------------------------------------------------------------------------
!     Deallocate variables that are no longer needed after initialisation stage.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(coord)) THEN; DEALLOCATE(coord,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(coordList)) THEN; DEALLOCATE(coordList,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(coordLL)) THEN; DEALLOCATE(coordLL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(mapInRoute)) THEN; DEALLOCATE(mapInRoute,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(routeMask)) THEN; DEALLOCATE(routeMask,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varDescList)) THEN; DEALLOCATE(varDescList,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varTypeList)) THEN; DEALLOCATE(varTypeList,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

    CASE ( 'end' )
!-------------------------------------------------------------------------------
!     This is the "main" call, to deallocate at the end of a run.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!     Index variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(land_index)) THEN; DEALLOCATE(LAND_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tile_index)) THEN; DEALLOCATE(TILE_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(soil_index)) THEN; DEALLOCATE(SOIL_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lice_index)) THEN; DEALLOCATE(LICE_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(land_mask)) THEN; DEALLOCATE(LAND_MASK,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(frac)) THEN; DEALLOCATE(FRAC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ice_fract)) THEN; DEALLOCATE(ICE_FRACT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ice_fract_ncat)) THEN; DEALLOCATE(ICE_FRACT_NCAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z1_uv)) THEN; DEALLOCATE(Z1_UV,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z1_tq)) THEN; DEALLOCATE(Z1_TQ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SSI_INDEX)) THEN; DEALLOCATE(SSI_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SEA_INDEX)) THEN; DEALLOCATE(SEA_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SICE_INDEX)) THEN; DEALLOCATE(SICE_INDEX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SICE_PTS_NCAT)) THEN; DEALLOCATE(SICE_PTS_NCAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SICE_INDEX_NCAT)) THEN; DEALLOCATE(SICE_INDEX_NCAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(FSSI)) THEN; DEALLOCATE(FSSI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SEA_FRAC)) THEN; DEALLOCATE(SEA_FRAC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SICE_FRAC)) THEN; DEALLOCATE(SICE_FRAC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SICE_FRAC_NCAT)) THEN; DEALLOCATE(SICE_FRAC_NCAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Runoff routing variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(routeIndex)) THEN;  DEALLOCATE( routeNext,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(routeNext)) THEN;  DEALLOCATE( routeNext,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(routeOrder)) THEN;  DEALLOCATE( routeOrder,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(roffAccumLand)) THEN;  DEALLOCATE( roffAccumLand,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Screen variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(q1p5m)) THEN; DEALLOCATE(Q1P5M,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(q1p5m_tile)) THEN; DEALLOCATE(Q1P5M_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(t1p5m)) THEN; DEALLOCATE(T1P5M,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(t1p5m_tile)) THEN; DEALLOCATE(T1P5M_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(u10m)) THEN; DEALLOCATE(U10M,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(v10m)) THEN; DEALLOCATE(V10M,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Plant and soil parameters
!-------------------------------------------------------------------------------
      IF (ALLOCATED(catch)) THEN; DEALLOCATE(CATCH,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(catch_snow)) THEN; DEALLOCATE(CATCH_SNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(cosz)) THEN; DEALLOCATE(COSZ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(diff_frac)) THEN; DEALLOCATE(DIFF_FRAC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

      IF (ALLOCATED(b)) THEN; DEALLOCATE(B,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dzsnow)) THEN; DEALLOCATE(DZSNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dzsoil)) THEN; DEALLOCATE(DZSOIL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ds)) THEN; DEALLOCATE(DS,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sathh)) THEN; DEALLOCATE(SATHH,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albsoil)) THEN; DEALLOCATE(ALBSOIL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(hcon)) THEN; DEALLOCATE(HCON,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(smvccl)) THEN; DEALLOCATE(SMVCCL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(smvcst)) THEN; DEALLOCATE(SMVCST,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(smvcwt)) THEN; DEALLOCATE(SMVCWT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(hcap)) THEN; DEALLOCATE(HCAP,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(satcon)) THEN; DEALLOCATE(SATCON,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(infil_tile)) THEN; DEALLOCATE(INFIL_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z0_tile)) THEN; DEALLOCATE(Z0_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sthu)) THEN; DEALLOCATE(STHU,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sthf)) THEN; DEALLOCATE(STHF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(soil_clay)) THEN; DEALLOCATE(SOIL_CLAY,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Surface type variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(z0h_z0m)) THEN; DEALLOCATE( z0h_z0m ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Veg surface type variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(albsnc_max)) THEN; DEALLOCATE(albsnc_max,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albsnc_min)) THEN; DEALLOCATE(albsnc_min,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albsnf_max)) THEN; DEALLOCATE(albsnf_max,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(alpar)) THEN; DEALLOCATE(alpar,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(alpha)) THEN; DEALLOCATE(alpha,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(alnir)) THEN; DEALLOCATE(alnir,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(a_wl)) THEN; DEALLOCATE(a_wl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(a_ws)) THEN; DEALLOCATE(a_ws,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(b_wl)) THEN; DEALLOCATE(b_wl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(c3)) THEN; DEALLOCATE(c3,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(catch0)) THEN; DEALLOCATE(catch0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dcatch_dlai)) THEN; DEALLOCATE(dcatch_dlai,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dgl_dm)) THEN; DEALLOCATE(dgl_dm,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dgl_dt)) THEN; DEALLOCATE(dgl_dt,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dqcrit)) THEN; DEALLOCATE(dqcrit,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dz0v_dh)) THEN; DEALLOCATE(dz0v_dh,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(emis_pft)) THEN; DEALLOCATE(emis_pft,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(eta_sl)) THEN; DEALLOCATE(eta_sl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(f0)) THEN; DEALLOCATE(f0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fd)) THEN; DEALLOCATE(fd,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fsmc_of)) THEN; DEALLOCATE(fsmc_of,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf_0)) THEN; DEALLOCATE(g_leaf_0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(glmin)) THEN; DEALLOCATE(glmin,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(infil_f)) THEN; DEALLOCATE(infil_f,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(kext)) THEN; DEALLOCATE(kext,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(kpar)) THEN; DEALLOCATE(kpar,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(neff)) THEN; DEALLOCATE(neff,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nl0)) THEN; DEALLOCATE(nl0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nr_nl)) THEN; DEALLOCATE(nr_nl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ns_nl)) THEN; DEALLOCATE(ns_nl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(omega)) THEN; DEALLOCATE(omega,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(omnir)) THEN; DEALLOCATE(omnir,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(orient)) THEN; DEALLOCATE(orient,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(r_grow)) THEN; DEALLOCATE(r_grow,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rootd_ft)) THEN; DEALLOCATE(rootd_ft,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sigl)) THEN; DEALLOCATE(sigl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tleaf_of)) THEN; DEALLOCATE(tleaf_of,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tlow)) THEN; DEALLOCATE(tlow,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tupp)) THEN; DEALLOCATE(tupp,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fl_o3_ct)) THEN; DEALLOCATE(fl_o3_ct,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dfp_dcuo)) THEN; DEALLOCATE(dfp_dcuo,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pftName)) THEN; DEALLOCATE(pftName,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Non-veg surface type variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(albsnc_nvg)) THEN; DEALLOCATE(albsnc_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albsnf_nvg)) THEN; DEALLOCATE(albsnf_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(catch_nvg)) THEN; DEALLOCATE(catch_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ch_nvg)) THEN; DEALLOCATE(ch_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(emis_nvg)) THEN; DEALLOCATE(emis_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gs_nvg)) THEN; DEALLOCATE(gs_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(infil_nvg)) THEN; DEALLOCATE(infil_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vf_nvg)) THEN; DEALLOCATE(vf_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z0_nvg)) THEN; DEALLOCATE(z0_nvg,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nvgName)) THEN; DEALLOCATE(nvgName,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     New urban scheme MORUSES variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(hgt)) THEN; DEALLOCATE(hgt,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(hwr)) THEN; DEALLOCATE(hwr,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(wrr)) THEN; DEALLOCATE(wrr,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(disp)) THEN; DEALLOCATE(disp,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ztm)) THEN; DEALLOCATE(ztm,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albwl)) THEN; DEALLOCATE(albwl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(albrd)) THEN; DEALLOCATE(albrd,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(emisw)) THEN; DEALLOCATE(emisw,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(emisr)) THEN; DEALLOCATE(emisr,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!     Height above gridbox mean.
      IF (ALLOCATED(surf_hgt)) THEN; DEALLOCATE(surf_hgt,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
!     Land height
      IF (ALLOCATED(z_land)) THEN; DEALLOCATE(z_land,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     TRIFFID variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(crop)) THEN; DEALLOCATE(crop ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_area)) THEN; DEALLOCATE(g_area ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_grow)) THEN; DEALLOCATE(g_grow ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_root)) THEN; DEALLOCATE(g_root ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_wood)) THEN; DEALLOCATE(g_wood ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lai_max)) THEN; DEALLOCATE(lai_max ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lai_min)) THEN; DEALLOCATE(lai_min ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Forcing variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(qw_1)) THEN; DEALLOCATE(QW_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tl_1)) THEN; DEALLOCATE(TL_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(u_0)) THEN; DEALLOCATE(U_0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(v_0)) THEN; DEALLOCATE(V_0,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(u_1)) THEN; DEALLOCATE(U_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(v_1)) THEN; DEALLOCATE(V_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pstar)) THEN; DEALLOCATE(PSTAR,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ls_rain)) THEN; DEALLOCATE(LS_RAIN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(con_rain)) THEN; DEALLOCATE(CON_RAIN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ls_snow)) THEN; DEALLOCATE(LS_SNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(con_snow)) THEN; DEALLOCATE(CON_SNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sw_down)) THEN; DEALLOCATE(SW_DOWN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lw_down)) THEN; DEALLOCATE(LW_DOWN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

! Ozone forcing variable
      IF (ALLOCATED(o3)) THEN; DEALLOCATE(o3,STAT=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(flux_o3_ft)) THEN; DEALLOCATE(flux_o3_ft,STAT=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fo3_ft)) THEN; DEALLOCATE(fo3_ft,STAT=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Variables connected to i/o of driving data (drive_io_vars).
!-------------------------------------------------------------------------------
      IF (ALLOCATED(driveFileDate)) THEN; DEALLOCATE(driveFileDate,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveFileName)) THEN; DEALLOCATE(driveFileName,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveFileTime)) THEN; DEALLOCATE(driveFileTime,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveData)) THEN; DEALLOCATE(driveData,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveDataIn)) THEN; DEALLOCATE(driveDataIn,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveUnit)) THEN; DEALLOCATE(driveUnit,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveFileOnUnit)) THEN; DEALLOCATE(driveFileOnUnit,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(driveVarNameUnit)) THEN; DEALLOCATE(driveVarNameUnit,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Prognostics variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(lai)) THEN; DEALLOCATE(LAI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(canht_ft)) THEN; DEALLOCATE(CANHT_FT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(smcl)) THEN; DEALLOCATE(SMCL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(t_soil)) THEN; DEALLOCATE(T_SOIL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rgrain)) THEN; DEALLOCATE(RGRAIN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rho_snow_grnd)) THEN; DEALLOCATE(RHO_SNOW_GRND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(routeStore)) THEN;  DEALLOCATE( routeStore,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snow_tile)) THEN; DEALLOCATE(SNOW_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(soot)) THEN; DEALLOCATE(SOOT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar_tile)) THEN; DEALLOCATE(TSTAR_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(canopy)) THEN; DEALLOCATE(CANOPY,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(canopy_gb)) THEN; DEALLOCATE(CANOPY_GB,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(cs)) THEN; DEALLOCATE(CS,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ti)) THEN; DEALLOCATE(TI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z0msea)) THEN; DEALLOCATE(Z0MSEA,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gs)) THEN; DEALLOCATE(GS,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gc)) THEN; DEALLOCATE(GC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(smc)) THEN; DEALLOCATE(SMC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(di)) THEN; DEALLOCATE(DI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(di_ncat)) THEN; DEALLOCATE(DI_NCAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snow_grnd)) THEN; DEALLOCATE(SNOW_GRND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snow_mass)) THEN; DEALLOCATE(SNOW_MASS,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snow_mass_sea)) THEN; DEALLOCATE(SNOW_MASS_SEA,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nsnow)) THEN; DEALLOCATE(NSNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rgrainl)) THEN; DEALLOCATE(RGRAINL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sice)) THEN; DEALLOCATE(SICE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sliq)) THEN; DEALLOCATE(SLIQ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snowdepth)) THEN; DEALLOCATE(SNOWDEPTH,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tsnow)) THEN; DEALLOCATE(TSNOW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snowdep_surf)) THEN; DEALLOCATE(snowdep_surf,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      
!-------------------------------------------------------------------------------
!     Forcing fluxes
!-------------------------------------------------------------------------------
      IF (ALLOCATED(alb_tile)) THEN; DEALLOCATE(ALB_TILE ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar)) THEN; DEALLOCATE(TSTAR,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(e_sea)) THEN; DEALLOCATE(E_SEA,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fqw_1)) THEN; DEALLOCATE(FQW_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fsmc)) THEN; DEALLOCATE(fsmc,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ftl_1)) THEN; DEALLOCATE(FTL_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ftl_tile)) THEN; DEALLOCATE(FTL_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(le_tile)) THEN; DEALLOCATE(LE_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(h_sea)) THEN; DEALLOCATE(H_SEA,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(taux_1)) THEN; DEALLOCATE(TAUX_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tauy_1)) THEN; DEALLOCATE(TAUY_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fqw_tile)) THEN; DEALLOCATE(FQW_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fqw_ice)) THEN; DEALLOCATE(FQW_ICE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ftl_ice)) THEN; DEALLOCATE(FTL_ICE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ecan)) THEN; DEALLOCATE(ECAN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(esoil_tile)) THEN; DEALLOCATE(ESOIL_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sea_ice_htf)) THEN; DEALLOCATE(SEA_ICE_HTF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_ht_flux)) THEN; DEALLOCATE(SURF_HT_FLUX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_htf_tile)) THEN; DEALLOCATE(surf_htf_tile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sice_mlt_htf)) THEN; DEALLOCATE(SICE_MLT_HTF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snomlt_surf_htf)) THEN; DEALLOCATE(SNOMLT_SURF_HTF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(land_albedo)) THEN; DEALLOCATE(LAND_ALBEDO ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(latent_heat)) THEN; DEALLOCATE(LATENT_HEAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ei)) THEN; DEALLOCATE(EI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ei_tile)) THEN; DEALLOCATE(EI_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ecan_tile)) THEN; DEALLOCATE(ECAN_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(esoil)) THEN; DEALLOCATE(ESOIL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ext)) THEN; DEALLOCATE(EXT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snowmelt)) THEN; DEALLOCATE(SNOWMELT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(melt_tile)) THEN; DEALLOCATE(MELT_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(hf_snow_melt)) THEN; DEALLOCATE(HF_SNOW_MELT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(radnet_tile)) THEN; DEALLOCATE(RADNET_TILE ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sw_tile)) THEN; DEALLOCATE(SW_TILE ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(emis_tile)) THEN; DEALLOCATE(EMIS_TILE ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snow_melt)) THEN; DEALLOCATE(SNOW_MELT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snomlt_sub_htf)) THEN; DEALLOCATE(SNOMLT_SUB_HTF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sub_surf_roff)) THEN; DEALLOCATE(SUB_SURF_ROFF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_roff)) THEN; DEALLOCATE(SURF_ROFF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tot_tfall)) THEN; DEALLOCATE(TOT_TFALL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_ht_store)) THEN; DEALLOCATE(surf_ht_store,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(anthrop_heat)) THEN; DEALLOCATE(anthrop_heat,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Aerosol variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(co2_3d)) THEN; DEALLOCATE(CO2_3D,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rho_cd_modv1)) THEN; DEALLOCATE(RHO_CD_MODV1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rho_aresist)) THEN; DEALLOCATE(RHO_ARESIST,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(aresist)) THEN; DEALLOCATE(ARESIST,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resist_b)) THEN; DEALLOCATE(RESIST_B,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rho_aresist_tile)) THEN; DEALLOCATE(RHO_ARESIST_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(aresist_tile)) THEN; DEALLOCATE(ARESIST_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resist_b_tile)) THEN; DEALLOCATE(RESIST_B_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(r_b_dust)) THEN; DEALLOCATE(R_B_DUST,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(cd_std_dust)) THEN; DEALLOCATE(CD_STD_DUST,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(u_s_std_tile)) THEN; DEALLOCATE(U_S_STD_TILE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Orographic roughness variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(sil_orog_land)) THEN; DEALLOCATE(SIL_OROG_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ho2r2_orog)) THEN; DEALLOCATE(HO2R2_OROG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(h_blend_orog)) THEN; DEALLOCATE(H_BLEND_OROG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(z0m_eff)) THEN; DEALLOCATE(Z0M_EFF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!     Grid-change variables
      IF (ALLOCATED(u_0_p)) THEN; DEALLOCATE(U_0_P,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(v_0_p)) THEN; DEALLOCATE(V_0_P,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(u_1_p)) THEN; DEALLOCATE(U_1_P,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(v_1_p)) THEN; DEALLOCATE(V_1_P,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dtrdz_charney_grid_1)) THEN; DEALLOCATE(DTRDZ_CHARNEY_GRID_1,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Triffid/veg variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(g_leaf_acc)) THEN; DEALLOCATE(G_LEAF_ACC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(npp_ft_acc)) THEN; DEALLOCATE(NPP_FT_ACC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_w_ft_acc)) THEN; DEALLOCATE(RESP_W_FT_ACC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_s_acc)) THEN; DEALLOCATE(RESP_S_ACC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf_phen_acc)) THEN; DEALLOCATE(G_LEAF_PHEN_ACC,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gpp)) THEN; DEALLOCATE(GPP,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(npp)) THEN; DEALLOCATE(NPP,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_p)) THEN; DEALLOCATE(RESP_P,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf)) THEN; DEALLOCATE(G_LEAF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf_phen)) THEN; DEALLOCATE(G_LEAF_PHEN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gpp_ft)) THEN; DEALLOCATE(GPP_FT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(npp_ft)) THEN; DEALLOCATE(NPP_FT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_p_ft)) THEN; DEALLOCATE(RESP_P_FT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_s)) THEN; DEALLOCATE(RESP_S,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_w_ft)) THEN; DEALLOCATE(RESP_W_FT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lai_phen)) THEN; DEALLOCATE(LAI_PHEN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(c_veg)) THEN; DEALLOCATE(C_VEG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(cv)) THEN; DEALLOCATE(CV,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf_day)) THEN; DEALLOCATE(G_LEAF_DAY,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(g_leaf_dr_out)) THEN; DEALLOCATE(G_LEAF_DR_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lit_c)) THEN; DEALLOCATE(LIT_C,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(lit_c_mn)) THEN; DEALLOCATE(LIT_C_MN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(npp_dr_out)) THEN; DEALLOCATE(NPP_DR_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_w_dr_out)) THEN; DEALLOCATE(RESP_W_DR_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(resp_s_dr_out)) THEN; DEALLOCATE(RESP_S_DR_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(frac_agr)) THEN; DEALLOCATE(FRAC_AGR,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     TOPMODEL and PDM variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(fexp)) THEN; DEALLOCATE(FEXP,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gamtot)) THEN; DEALLOCATE(GAMTOT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ti_mean)) THEN; DEALLOCATE(TI_MEAN,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ti_sig)) THEN; DEALLOCATE(TI_SIG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fsat)) THEN; DEALLOCATE(FSAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fwetl)) THEN; DEALLOCATE(FWETL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(zw)) THEN; DEALLOCATE(ZW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(dun_roff)) THEN; DEALLOCATE(DUN_ROFF,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(qbase)) THEN; DEALLOCATE(QBASE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(qbase_zw)) THEN; DEALLOCATE(QBASE_ZW,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(fch4_wetl)) THEN; DEALLOCATE(FCH4_WETL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(drain)) THEN; DEALLOCATE(drain,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(sthzw)) THEN; DEALLOCATE(sthzw,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(a_fsat)) THEN; DEALLOCATE(a_fsat,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(c_fsat)) THEN; DEALLOCATE(c_fsat,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(a_fwet)) THEN; DEALLOCATE(a_fwet,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(c_fwet)) THEN; DEALLOCATE(c_fwet,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(inlandout_atm)) THEN; DEALLOCATE(inlandout_atm,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Coastal tiling variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(fland)) THEN; DEALLOCATE(FLAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(flandg)) THEN; DEALLOCATE(FLANDG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar_land)) THEN; DEALLOCATE(TSTAR_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar_sea)) THEN; DEALLOCATE(TSTAR_SEA,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar_sice)) THEN; DEALLOCATE(TSTAR_SICE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tstar_ssi)) THEN; DEALLOCATE(TSTAR_SSI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(taux_land)) THEN; DEALLOCATE(TAUX_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(taux_ssi)) THEN; DEALLOCATE(TAUX_SSI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tauy_land)) THEN; DEALLOCATE(TAUY_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tauy_ssi)) THEN; DEALLOCATE(TAUY_SSI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vshr_land)) THEN; DEALLOCATE(VSHR_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vshr_ssi)) THEN; DEALLOCATE(VSHR_SSI,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_ht_flux_land)) THEN; DEALLOCATE(SURF_HT_FLUX_LAND,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(surf_ht_flux_sice)) THEN; DEALLOCATE(SURF_HT_FLUX_SICE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     ANCIL variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(tile_pts)) THEN; DEALLOCATE( tile_pts ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pftUse)) THEN; DEALLOCATE( pftUse ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Snow variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(canSnowTile)) THEN; DEALLOCATE( canSnowTile ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Location variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(latitude)) THEN; DEALLOCATE(LATITUDE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(longitude)) THEN; DEALLOCATE(LONGITUDE,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Spin up variables.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(spinValOld)) THEN; DEALLOCATE( spinValOld ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Input/output variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(havePFT)) THEN; DEALLOCATE(havePFT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(haveSCpool)) THEN; DEALLOCATE(haveSCpool,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(haveSnow)) THEN; DEALLOCATE(haveSnow,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(haveSoil)) THEN; DEALLOCATE(haveSoil,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(haveTile)) THEN; DEALLOCATE(haveTile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(haveType)) THEN; DEALLOCATE(haveType,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(irecPrevOut)) THEN; DEALLOCATE(irecPrevOut,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(mapIn)) THEN; DEALLOCATE(mapIn,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(mapInLand)) THEN; DEALLOCATE(mapInLand,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(mapOut)) THEN; DEALLOCATE(mapOut,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      if (ALLOCATED(mapOutCompress)) THEN; DEALLOCATE(mapoutCompress,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(mapOutLand)) THEN; DEALLOCATE(mapOutLand,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nlevMax)) THEN; DEALLOCATE(nlevMax,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nlevMaxCtl)) THEN; DEALLOCATE(nlevMaxCtl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ntCtl)) THEN; DEALLOCATE(ntCtl,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ntCtlNeed)) THEN; DEALLOCATE(ntCtlNeed,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ntOutFilePer)) THEN; DEALLOCATE(ntOutFilePer,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(ntOutPer)) THEN; DEALLOCATE(ntOutPer,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nvarOut)) THEN; DEALLOCATE(nvarOut,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(nxyMax)) THEN; DEALLOCATE(nxyMax,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(openedFileName)) THEN; DEALLOCATE(openedFileName,stat=ierr)
                   ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outCtlFile)) THEN; DEALLOCATE(outCtlFile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outActivePrev)) THEN; DEALLOCATE(outActivePrev,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outAreaLL)) THEN; DEALLOCATE(outAreaLL,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outDataFile)) THEN; DEALLOCATE(outDataFile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outDate)) THEN; DEALLOCATE(outDate,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outDateFlag)) THEN; DEALLOCATE(outDateFlag,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outEndPos)) THEN; DEALLOCATE(outEndPos,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outFirstActive)) THEN; DEALLOCATE(outFirstActive,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outFirstSection)) THEN; DEALLOCATE(outFirstSection,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outFirstWrite)) THEN; DEALLOCATE(outFirstWrite,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outFilePer)) THEN; DEALLOCATE(outFilePer,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outFileStep)) THEN; DEALLOCATE(outFileStep,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outLLorder)) THEN; DEALLOCATE(outLLorder,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outName)) THEN; DEALLOCATE(outName,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outNpWrite)) THEN; DEALLOCATE(outNpWrite,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outPer)) THEN; DEALLOCATE(outPer,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outPerNunits)) THEN; DEALLOCATE(outPerNunits,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outSamPer)) THEN; DEALLOCATE(outSamPer,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outStep)) THEN; DEALLOCATE(outStep,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outStepSamPer)) THEN; DEALLOCATE(outStepSamPer,stat=ierr)
                          ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outTime)) THEN; DEALLOCATE(outTime,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outTimeID)) THEN; DEALLOCATE(outTimeID,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outUnit)) THEN; DEALLOCATE(outUnit,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outval)) THEN; DEALLOCATE(outval,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outWriteCount)) THEN; DEALLOCATE(outWriteCount,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(compressGridFile)) THEN; DEALLOCATE(compressGridFile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outGridDxy)) THEN; DEALLOCATE(outGridDxy,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outGridNxy)) THEN; DEALLOCATE(outGridNxy,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outGridXY)) THEN; DEALLOCATE(outGridXY,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outRangeX)) THEN; DEALLOCATE(outRangeX,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outRangeY)) THEN; DEALLOCATE(outRangeY,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pointsFlag)) THEN; DEALLOCATE(pointsFlag,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pointsOut)) THEN; DEALLOCATE(pointsOut,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(pointsOutLand)) THEN; DEALLOCATE(pointsOutLand,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rgProfile)) THEN; DEALLOCATE(rgProfile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rpProfile)) THEN; DEALLOCATE(rpProfile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

      IF (ALLOCATED(snapProfile)) THEN; DEALLOCATE(snapProfile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(taccumVar)) THEN; DEALLOCATE(taccumVar,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outTemplate)) THEN; DEALLOCATE(outTemplate,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tmeanProfile)) THEN; DEALLOCATE(tmeanProfile,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(tmeanVar)) THEN; DEALLOCATE(tmeanVar,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outCompress)) THEN; DEALLOCATE(outCompress,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(useCompressGrid)) THEN; DEALLOCATE(useCompressGrid,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varDesc)) THEN; DEALLOCATE(varDesc,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varName)) THEN; DEALLOCATE(varName,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varNameList)) THEN; DEALLOCATE(varNameList,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(outVarID)) THEN; DEALLOCATE(outVarID,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varNlev)) THEN; DEALLOCATE(varNlev,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varNum)) THEN; DEALLOCATE(varNum,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varPos)) THEN; DEALLOCATE(varPos,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varStartPos)) THEN; DEALLOCATE(varStartPos,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varType)) THEN; DEALLOCATE(varType,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(varUnitsList)) THEN; DEALLOCATE(varUnitsList,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Output variables that are only allocated if diagnostic is requested.
!-------------------------------------------------------------------------------
      IF (ALLOCATED(ciDiag)) THEN; DEALLOCATE(ciDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(gstomDiag)) THEN; DEALLOCATE(gstomDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rdcDiag)) THEN; DEALLOCATE(rdcDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rflowDiag)) THEN; DEALLOCATE(rflowDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(roffInfDiag)) THEN; DEALLOCATE(roffInfDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(rrunDiag)) THEN; DEALLOCATE(rrunDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(snowGMeltDiag)) THEN; DEALLOCATE(snowGmeltDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(wfluxDiag)) THEN; DEALLOCATE(wfluxDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(wfluxSfcDiag)) THEN; DEALLOCATE(wfluxSfcDiag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF

!-------------------------------------------------------------------------------
!     Varying veg variables (veg_io_vars)
!-------------------------------------------------------------------------------
      IF (ALLOCATED(vegDataIn)) THEN; DEALLOCATE(vegDataIn ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarPos)) THEN; DEALLOCATE(vegVarPos ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarStash)) THEN; DEALLOCATE(vegVarStash ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarFlag)) THEN; DEALLOCATE(vegVarFlag ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarInterp)) THEN; DEALLOCATE(vegVarInterp ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarName)) THEN; DEALLOCATE(vegVarName ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegVarNameFile)) THEN; DEALLOCATE(vegVarNameFile ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegFileName)) THEN; DEALLOCATE(vegFileName ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegFileDate)) THEN; DEALLOCATE(vegFileDate ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegFileTime)) THEN; DEALLOCATE(vegFileTime ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegUnit)) THEN; DEALLOCATE(vegUnit ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(vegUnitFile)) THEN; DEALLOCATE(vegUnitFile ,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      
!-------------------------------------------------------------------------------
! IMOGEN variables
!-------------------------------------------------------------------------------
      IF (ALLOCATED(T_CLIM)) THEN; DEALLOCATE(T_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(RAINFALL_CLIM)) THEN; DEALLOCATE(RAINFALL_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SNOWFALL_CLIM)) THEN; DEALLOCATE(SNOWFALL_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(RH15M_CLIM)) THEN; DEALLOCATE(RH15M_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(UWIND_CLIM)) THEN; DEALLOCATE(UWIND_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(VWIND_CLIM)) THEN; DEALLOCATE(VWIND_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(DTEMP_CLIM)) THEN; DEALLOCATE(DTEMP_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(PSTAR_HA_CLIM)) THEN; DEALLOCATE(PSTAR_HA_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SW_CLIM)) THEN; DEALLOCATE(SW_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(LW_CLIM)) THEN; DEALLOCATE(LW_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(F_WET_CLIM)) THEN; DEALLOCATE(F_WET_CLIM,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      
      IF (ALLOCATED(T_OUT)) THEN; DEALLOCATE(T_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(CONV_RAIN_OUT)) THEN; DEALLOCATE(CONV_RAIN_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(CONV_SNOW_OUT)) THEN; DEALLOCATE(CONV_SNOW_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(LS_RAIN_OUT)) THEN; DEALLOCATE(LS_RAIN_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(LS_SNOW_OUT)) THEN; DEALLOCATE(LS_SNOW_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(QHUM_OUT)) THEN; DEALLOCATE(QHUM_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(WIND_OUT)) THEN; DEALLOCATE(WIND_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(PSTAR_OUT)) THEN; DEALLOCATE(PSTAR_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(SW_OUT)) THEN; DEALLOCATE(SW_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(LW_OUT)) THEN; DEALLOCATE(LW_OUT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      
      IF (ALLOCATED(LAT)) THEN; DEALLOCATE(LAT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(LONG)) THEN; DEALLOCATE(LONG,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      IF (ALLOCATED(DCTOT)) THEN; DEALLOCATE(DCTOT,stat=ierr); ierrSum=ierrSum+ierr; ENDIF
      
!-------------------------------------------------------------------------------
!   Anything other that the two cases above is an error.
!-------------------------------------------------------------------------------
    CASE default

      WRITE(*,*)'ERROR: deallocate_arrays: no code for callType=',TRIM(callType)
      STOP

  END SELECT

  IF ( ierrSum /= 0 ) CALL allocate_Error( 'dealloc',ierrSum,'deallocate_arrays:'//TRIM(callType) )

END SUBROUTINE deallocate_arrays
!###############################################################################
!###############################################################################

