! module alloc_mod
!
! Contains:
!   subroutine allocate_arrays

!###############################################################################

  MODULE alloc_mod

  IMPLICIT NONE

  CONTAINS

!###############################################################################

! subroutine allocate_arrays
! Internal procedure in module alloc_mod.
! Subroutine to allocate space for arrays, and (for some) initialise.
! Remember to also deallocate in deallocate_arrays!
! Note that this routine is in a module to give explicit interface for optional argument.

  SUBROUTINE ALLOCATE_ARRAYS( callName,np )

  USE Aero, ONLY :  &
!  imported arrays with intent(out)
     aresist,ARESIST_TILE,CD_STD_DUST,co2_3d,R_B_DUST,RESIST_B,RESIST_B_TILE  &
    ,RHO_ARESIST,RHO_ARESIST_TILE,RHO_CD_MODV1,U_S_STD_TILE

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     co2_dim_len,co2_dim_row,land_pts,land_pts_trif,n_rows,nice,npft_trif  &
    ,nsmax,ntiles,row_length,rows,sm_levels,ssi_pts  &
!  imported arrays with intent(inout)
   ,frac,ice_fract,ice_fract_ncat,land_index,land_mask,lice_index,soil_index  &
   ,tile_index,tile_pts,z1_uv,z1_tq,sice_pts_ncat,ssi_index,sea_index  &
   ,sice_index,sice_index_ncat,fssi,sea_frac,sice_frac,sice_frac_ncat  &
   ,dim_cs1

  USE dust_param, ONLY :  &
!  imported scalar parameters
     ndiv

  USE C_elevate, ONLY :  &
!  imported arrays
     surf_hgt,z_land

  USE c_z0h_z0m, ONLY :  &
!  imported arrays with intent(inout)
     z0h_z0m

  USE Coastal

  USE drive_io_vars, ONLY :  &
!  imported scalars with intent(in)
     driveDataPer,ndriveFileTime,ndriveUnit,ndriveVar,ndriveVarIn  &
!  imported arrays with intent(in)
    ,driveTimeIndex  &
!  imported arrays with intent(out)
    ,driveData,driveDataIn,driveFileDate,driveFileName,driveFileTime  &
    ,driveFileOnUnit,driveVarNameUnit,driveUnit

  USE Fluxes

  USE Forcing

  USE inout, ONLY :  &
!  imported scalars with intent(in)
     nout,npoints,nvar,nvarOutTot,outLen,outLenWrite,pointsOutMax  &
!  imported arrays with intent(out)
    ,coordList,havePFT,haveSCpool,haveSnow,haveSoil,haveTile,haveType  &
    ,irecPrevOut,mapIn,mapInLand  &
    ,mapOut,mapOutCompress,mapOutLand,nlevMax,nlevMaxCtl,ntCtl,ntCtlNeed  &
    ,ntOutFilePer,ntOutPer,nvarOut,nxyMax,openedFilename,outActivePrev,outAreaLL,outCtlFile  &
    ,outDataFile,outName,compressGridFile  &
    ,outGridDxy,outGridNxy,outGridXY,outLLorder,outNpWrite,outRangeX,outRangeY  &
    ,outVarID,pointsFlag,pointsOut,pointsOutLand   &
    ,outDate,outDateFlag,outEndPos,outFilePer,outFileStep,outFirstActive   &
    ,outFirstSection,outFirstWrite,outPer,outPerNunits,outSamPer,outStep  &
    ,outStepSamPer,outTime,outTimeID,outUnit   &
    ,outval,outWriteCount,pftUse,rgProfile,rpProfile,snapProfile,taccumVar &
    ,tmeanVar,outTemplate,tmeanProfile,outCompress,useCompressGrid,varDesc  &
    ,varDescList,varName,varNameList    &
    ,varNlev,varNum,varPos,varStartPos,varType,varTypeList,varUnitsList

  USE misc_utils, ONLY :  &
!  imported procedures
     allocate_error

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
    nnvg,npft,ntype

  USE nvegparm, ONLY :  &
!  imported arrays with intent(inout)
     albsnc_nvg,albsnf_nvg,catch_nvg,gs_nvg,infil_nvg,z0_nvg,ch_nvg,vf_nvg,emis_nvg,nvgName

  USE offline_diag, ONLY :  &
!  imported scalars with intent(in)
     useCiDiag,useGstomDiag,useRdcDiag,useRflowDiag,useRoffInfDiag  &
    ,useRrunDiag,useSnowGMeltDiag,useWfluxDiag,useWfluxSfcDiag  &
!  imported arrays with intent(out)
    ,ciDiag,gstomDiag,rdcDiag,rflowDiag,roffInfDiag  &
    ,rrunDiag,snowGmeltDiag,wfluxDiag,wfluxSfcDiag

  USE Orog
  USE P_s_parms
  USE pftparm
  USE Prognostics

  USE route_mod, ONLY :  &
!  imported scalars with intent(in)
     npRoute,nxRoute,nyRoute  &
!  imported arrays with intent(inout)
    ,mapInRoute,routeIndex,routeMask,routeNext,routeOrder,roffAccumLand

  USE Screen

  USE snow_param, ONLY :  &
!  imported arrays with intent(out)
     canSnowTile,ds,dzsnow

  USE surf_param, ONLY : diff_frac

  USE soil_param, ONLY : dzsoil

  USE spin_mod, ONLY :  &
!  imported scalars with intent(in)
     npSpinMax,nspinVar,nzSpinMax  &
!  imported arrays with intent(inout)
    ,spinValOld

  USE switches, ONLY :   &
! imported scalars with intent(in)
    l_phenol,l_triffid,routeOnly

  USE time_loc, ONLY :  &
!  imported arrays with intent(inout)
     latitude,longitude

  USE Top_pdm, ONLY :  &
!  imported arrays with intent(out)
     a_fsat,a_fwet,c_fsat,c_fwet,drain,dun_roff,fch4_wetl,fexp,fsat,fwetl  &
    ,gamtot,qbase,qbase_zw,sthzw,ti_mean,ti_sig,zw,INLANDOUT_ATM

  USE trif, ONLY :  &
!  imported arrays with intent(inout)
     crop,g_area,g_grow,g_root,g_wood,lai_max,lai_min

  USE Trifctl
  USE U_v_grid

  USE veg_io_vars, ONLY : &
!  imported scalars with intent(in)
     nvegFileTime,nvegFileVar,nvegVar  &
!  imported arrays with intent(in)
    ,vegTimeIndex  &
!  imported arrays with intent(out)
    ,vegDataIn,vegFileDate,vegFileName,vegFileTime,vegUnit,vegUnitFile  &
    ,vegVarFlag,vegVarInterp,vegVarName,vegvarNameFile,vegVarPos,vegVarStash

  USE switches_urban, ONLY : l_urban2t, l_moruses
  USE urban_param, ONLY :                                                     &
     hgt, hwr, wrr, disp, ztm, albwl, albrd, emisw, emisr
     
  USE ozone_vars, ONLY : o3, flux_o3_ft, fo3_ft
  
  USE imogen_time, ONLY : mm, md, nsdmax
  
  USE imogen_clim, ONLY :                                         &
    T_CLIM,RAINFALL_CLIM,SNOWFALL_CLIM,RH15M_CLIM,UWIND_CLIM,     &
    VWIND_CLIM,DTEMP_CLIM,PSTAR_HA_CLIM,SW_CLIM,LW_CLIM,          &
    F_WET_CLIM,LAT,LONG,DCTOT
    
  USE imogen_drive_vars, ONLY :                                   &
    T_OUT,CONV_RAIN_OUT,CONV_SNOW_OUT,LS_RAIN_OUT,LS_SNOW_OUT,    &
    QHUM_OUT,WIND_OUT,PSTAR_OUT,SW_OUT,LW_OUT
    
  USE jules_mod, ONLY : snowdep_surf

!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, INTENT(in), OPTIONAL :: &!  optional in scalars
    np          !  used to provide extra information for certain values of callName

  CHARACTER(len=*), INTENT(in) ::  &!  in scalars
    callName    !   identifies where this code has been called from

  INTEGER ::  &!  local scalars
    ierr,ierrSum  &!  used for error checking
   ,iout           !  work

!  INTEGER :: nsmax_tmp

!-----------------------------------------------------------------------

! Initialise.
  ierrSum = 0

  SELECT CASE ( callName )
!-----------------------------------------------------------------------
    CASE ( 'init_opts' )
      ALLOCATE( nvgName(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( pftName(npft), stat=ierr ); ierrSum=ierrSum+ierr

!-----------------------------------------------------------------------
    CASE ( 'init_grid 1' )
!     Note that these variables may be deallocated and then reallocated
!     by subroutine init_grid_realloc.
      ALLOCATE( mapIn(npoints), STAT=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( flandg(row_length,rows), STAT=ierr); ierrSum=ierrSum+ierr
      ALLOCATE( land_mask(row_length,rows), STAT=ierr ); ierrsum=ierrsum+ierr
      ALLOCATE( latitude(row_length,rows), STAT=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( longitude(row_length,rows), STAT=ierr ); ierrSum=ierrSum+ierr
!-----------------------------------------------------------------------
    CASE ( 'init_grid 2' )

!     "Main" call.
!     Note that we will already have stopped if land_pts=0.

!     Index variables
      ALLOCATE( LAND_INDEX(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!     Coastal tiling variables
      ALLOCATE( FLAND(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!     Input mapping.
      ALLOCATE( mapInLand(land_pts), stat=ierr ); ierrSum=ierrSum+ierr

!     Runoff components.
      ALLOCATE( SUB_SURF_ROFF(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
      ALLOCATE( SURF_ROFF(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!     Variables that are not required if only doing routing.
      IF ( .NOT. routeOnly ) THEN

!       Index variables
        ALLOCATE( TILE_INDEX(LAND_PTS,NTYPE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SOIL_INDEX(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LICE_INDEX(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FRAC(LAND_PTS,NTYPE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ICE_FRACT(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ICE_FRACT_NCAT(ROW_LENGTH,ROWS,NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Z1_UV(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Z1_TQ(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SSI_INDEX(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SEA_INDEX(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SICE_INDEX(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SICE_PTS_NCAT(NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SICE_INDEX_NCAT(ROW_LENGTH*ROWS,NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FSSI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Screen variables
        ALLOCATE( Q1P5M(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Q1P5M_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( T1P5M(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( T1P5M_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( U10M(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( V10M(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Plant and soil parameters
        ALLOCATE( ALBSOIL(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CATCH(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CATCH_SNOW(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( COSZ(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DIFF_FRAC(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

        ALLOCATE( B(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DZSOIL(SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SATHH(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( HCON(LAND_PTS,0:SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SMVCCL(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SMVCST(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SMVCWT(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( HCAP(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SATCON(LAND_PTS,0:SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( INFIL_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Z0_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( STHU(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( STHF(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SOIL_CLAY(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Surface type variables.
        ALLOCATE( z0h_z0m(ntype), stat=ierr ); ierrSum=ierrSum+ierr

!       Veg surface type variables.
        ALLOCATE( albsnc_max(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( albsnc_min(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( albsnf_max(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( alpha(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( alnir(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( alpar(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( a_wl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( a_ws(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( b_wl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( catch0(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( c3(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dcatch_dlai(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dgl_dm(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dgl_dt(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dqcrit(npft), stat=ierr ); ierrSum=ierrSum+ierr
        !DSM ALLOCATE( dz0v_dh(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dz0v_dh(land_pts,npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( emis_pft(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( eta_sl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( fd(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( fsmc_of(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( f0(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( glmin(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( g_leaf_0(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( infil_f(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( kext(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( kpar(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( neff(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( nl0(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( nr_nl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( ns_nl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( omega(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( omnir(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( orient(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( r_grow(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( rootd_ft(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( sigl(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( tleaf_of(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( tlow(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( tupp(npft), stat=ierr ); ierrSum=ierrSum+ierr
!       Ozone damage parameters
        ALLOCATE( fl_o3_ct(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( dfp_dcuo(npft), stat=ierr ); ierrSum=ierrSum+ierr

!       Non-veg surface type variables.
        ALLOCATE( albsnc_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( albsnf_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( catch_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( emis_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( gs_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( infil_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( z0_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( ch_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr
        ALLOCATE( vf_nvg(nnvg), stat=ierr ); ierrSum=ierrSum+ierr

!       Height above mean grid-box
        ALLOCATE( surf_hgt(land_pts,ntiles), stat=ierr ); ierrSum=ierrSum+ierr
!       Land height
        ALLOCATE( Z_LAND(ROW_LENGTH,ROWS), stat=ierr ); ierrSum=ierrSum+ierr

!       TRIFFID variables - only needed if TRIFFID and/or phenology is selected.
        IF ( l_triffid .OR. l_phenol ) THEN
          ALLOCATE( crop(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( g_area(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( g_grow(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( g_root(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( g_wood(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( lai_max(npft), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( lai_min(npft), stat=ierr ); ierrSum=ierrSum+ierr
        ENDIF

!       Forcing variables
        ALLOCATE( QW_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TL_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( U_0(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( V_0(ROW_LENGTH,N_ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( U_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( V_1(ROW_LENGTH,N_ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( PSTAR(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LS_RAIN(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CON_RAIN(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LS_SNOW(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CON_SNOW(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SW_DOWN(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LW_DOWN(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Ozone variables
        ALLOCATE( o3(land_pts), STAT=ierr ); ierrsum = ierrsum + ierr
        ALLOCATE( flux_o3_ft(land_pts,npft), STAT=ierr ); ierrsum = ierrsum + ierr
        ALLOCATE( fo3_ft(land_pts,npft), STAT=ierr ); ierrsum = ierrsum + ierr

!       Prognostics variables
        ALLOCATE( LAI(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CANHT_FT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SMCL(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( T_SOIL(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RGRAIN(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RHO_SNOW_GRND(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOW_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SOOT(ROW_LENGTH*ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TSTAR_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CANOPY(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CANOPY_GB(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CS(LAND_PTS,DIM_CS1), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Z0MSEA(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( GS(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( GC(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SMC(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DI_NCAT(ROW_LENGTH,ROWS,NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOW_GRND(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOW_MASS(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOW_MASS_SEA(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( NSNOW(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOWDEPTH(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOWDEP_SURF(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
!       Snow layer variables.
!       Some compilers seem to be unhappy passing these if they have not been
!       allocated, others don't mind. So now allocating even if nsmax=0.
        ALLOCATE( DZSNOW(nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DS(LAND_PTS,NTILES,nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RGRAINL(LAND_PTS,NTILES,nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SICE(LAND_PTS,NTILES,nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SLIQ(LAND_PTS,NTILES,nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TSNOW(LAND_PTS,NTILES,nsmax), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       (Forcing) fluxes
        ALLOCATE( ALB_TILE(LAND_PTS,NTILES,4) )
        ALLOCATE( TSTAR(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( E_SEA(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FQW_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( fsmc(land_pts,npft), STAT=ierr ); ierrsum=ierrsum+ierr
        ALLOCATE( FTL_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FTL_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LE_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( H_SEA(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUX_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUY_1(ROW_LENGTH,N_ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FQW_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FQW_ICE(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FTL_ICE(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ECAN(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ESOIL_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SEA_ICE_HTF(ROW_LENGTH,ROWS,NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SURF_HT_FLUX(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( surf_htf_tile(land_pts,ntiles), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SICE_MLT_HTF(ROW_LENGTH,ROWS,NICE), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOMLT_SURF_HTF(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LAND_ALBEDO(ROW_LENGTH,ROWS,4) )
        ALLOCATE( LATENT_HEAT(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( EI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( EI_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ECAN_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ESOIL(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( EXT(LAND_PTS,SM_LEVELS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOWMELT(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( MELT_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( HF_SNOW_MELT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RADNET_TILE(LAND_PTS,NTILES) )
        ALLOCATE( SW_TILE(LAND_PTS,NTILES) )
        ALLOCATE( EMIS_TILE(LAND_PTS,NTILES) )
        ALLOCATE( SNOW_MELT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SNOMLT_SUB_HTF(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TOT_TFALL(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( surf_ht_store(land_pts,ntiles), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( anthrop_heat(land_pts,ntiles), STAT=IERR ); IERRSUM=IERRSUM+ierr

        IF ( IERRSUM == 0 ) THEN
          alb_tile(:,:,:) = 0.0
          ftl_tile(:,:) = 0.0
          le_tile(:,:) = 0.0
          fqw_tile(:,:) = 0.0
          esoil_tile(:,:) = 273.15
          ecan_tile(:,:) = 0.0
          ei_tile(:,:) = 0.0
          melt_tile(:,:) = 0.0
          radnet_tile(:,:) = 0.0
          sw_tile(:,:) = 0.0
          emis_tile(:,:) = 0.0
          snow_grnd(:,:) = 0.0
          surf_htf_tile(:,:) = 0.0
          surf_ht_store(:,:) = 0.0
          anthrop_heat(:,:) = 0.0
        ENDIF

!       Aerosol variables
        ALLOCATE( CO2_3D(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RHO_CD_MODV1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RHO_ARESIST(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ARESIST(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESIST_B(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RHO_ARESIST_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ARESIST_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESIST_B_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( R_B_DUST(ROW_LENGTH,ROWS,NDIV), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CD_STD_DUST(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( U_S_STD_TILE(LAND_PTS,NTILES), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Orographic roughness variables
        ALLOCATE( SIL_OROG_LAND(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( HO2R2_OROG(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( H_BLEND_OROG(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( Z0M_EFF(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Grid-change variables
        ALLOCATE( U_0_P(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( V_0_P(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( U_1_P(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( V_1_P(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DTRDZ_CHARNEY_GRID_1(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Triffid variables
        ALLOCATE( G_LEAF_ACC(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( G_LEAF_PHEN_ACC(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( GPP(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( NPP(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_P(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( G_LEAF(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( G_LEAF_PHEN(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( GPP_FT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( NPP_FT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_P_FT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_S(LAND_PTS,DIM_CS1), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_W_FT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LAI_PHEN(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( C_VEG(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( CV(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( G_LEAF_DAY(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( G_LEAF_DR_OUT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LIT_C(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( LIT_C_MN(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( NPP_DR_OUT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_W_DR_OUT(LAND_PTS,NPFT), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( RESP_S_DR_OUT(LAND_PTS,5), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FRAC_AGR(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       TOPMODEL and PDM variables
        ALLOCATE( A_FSAT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( A_FWET(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( C_FSAT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( C_FWET(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DRAIN(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( DUN_ROFF(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FCH4_WETL(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FEXP(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FSAT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( FWETL(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( GAMTOT(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( QBASE(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( QBASE_ZW(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( STHZW(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TI_MEAN(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TI_SIG(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( ZW(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( INLANDOUT_ATM(LAND_PTS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       Coastal tiling variables
        ALLOCATE( TSTAR_LAND(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TSTAR_SEA(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TSTAR_SICE(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TSTAR_SSI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUX_LAND(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUX_SSI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUY_LAND(ROW_LENGTH,N_ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( TAUY_SSI(ROW_LENGTH,N_ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( VSHR_LAND(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( VSHR_SSI(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SURF_HT_FLUX_LAND(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr
        ALLOCATE( SURF_HT_FLUX_SICE(ROW_LENGTH,ROWS), STAT=IERR ); IERRSUM=IERRSUM+ierr

!       ANCIL variables.
        ALLOCATE( tile_pts(ntype), stat=ierr ); ierrSum=ierrSum+ierr

!       Variables related to i/o of PFT data.
        ALLOCATE( pftUse(npft), stat=ierr ); ierrSum=ierrSum+ierr

!       Snow variables.
        ALLOCATE( canSnowTile(ntiles), stat=ierr ); ierrSum=ierrSum+ierr

      ENDIF   !  routeOnly

!------------------------------------------------------------------------------
!   For two-tile urban schemes
!------------------------------------------------------------------------------
    CASE ( 'init_urban' )
      IF ( l_urban2t ) THEN
        ALLOCATE( wrr(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
        IF ( ierrSum == 0.0 ) wrr(:) = 0.0
        IF ( l_moruses ) THEN
          ALLOCATE( hgt(land_pts)  , stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( hwr(land_pts)  , stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( disp(land_pts) , stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( ztm(land_pts)  , stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( albwl(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( albrd(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( emisw(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
          ALLOCATE( emisr(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
          IF  ( ierrSum == 0.0 ) THEN
            hgt(:)   = 0.0
            hwr(:)   = 0.0
            albwl(:) = 0.0
            albrd(:) = 0.0
            emisw(:) = 0.0
            emisr(:) = 0.0
            ztm(:)   = 0.0
            disp(:)  = 0.0
          ENDIF
        ENDIF
      ENDIF

!-------------------------------------------------------------------------------
    CASE ( 'init_drive 1' )
      ALLOCATE( driveFileName(ndriveFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr
!     Allocate an extra space in fileDate and time, so that we can always compare with file i+1.
      ALLOCATE( driveFileDate(ndriveFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( driveFileTime(ndriveFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_drive 2' )
      ALLOCATE( driveUnit(ndriveUnit), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( driveFileOnUnit(ndriveUnit), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( driveVarNameUnit(ndriveUnit), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_drive 3' )
      ALLOCATE( driveDataIn(ndriveVarIn,row_length,rows,driveTimeIndex(1):driveTimeIndex(2)), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( driveData(ndriveVar,row_length,rows,driveDataPer), stat=ierr ); ierrSum=ierrSum+ierr

      IF ( ierrSum == 0 ) THEN
        driveData(:,:,:,:) = 0.0
        driveDataIn(:,:,:,:) = 0.0
      ENDIF

!-------------------------------------------------------------------------------
    CASE ( 'init_route 1' )
      ALLOCATE( mapInRoute(nxRoute*nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( routeMask(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_route 2' )
      ALLOCATE( routeIndex(npRoute), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( routeNext(npRoute), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( routeOrder(npRoute), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( routeStore(nxRoute,nyRoute), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( roffAccumLand(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
!-------------------------------------------------------------------------------
    CASE ( 'init_veg_vary 1' )
      ALLOCATE( vegVarPos(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegVarFlag(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegVarInterp(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegVarName(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegVarNameFile(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegVarStash(nvegVar), stat=ierr ); ierrSum=ierrSum+ierr
!     Allocate an extra space in certain fields, so that we can always compare with file i+1.
      ALLOCATE( vegFileName(nvegFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegFileDate(nvegFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegFileTime(nvegFileTime+1), stat=ierr ); ierrSum=ierrSum+ierr
!-------------------------------------------------------------------------------
    CASE ( 'init_veg_vary 2' )
      ALLOCATE( vegDataIn(nvegvar,land_pts,npft,vegTimeIndex(1):vegTimeIndex(2)), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegUnit(nvegFileVar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( vegUnitFile(nvegFileVar), stat=ierr ); ierrSum=ierrSum+ierr

      IF ( ierrSum == 0 ) vegDataIn(:,:,:,:) = 0.0
!-------------------------------------------------------------------------------
    CASE ( 'init_out 1' )

!     Integers.
      ALLOCATE( irecPrevOut(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( nlevMax(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( nlevMaxCtl(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( ntCtl(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( ntCtlNeed(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( ntOutFilePer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( ntOutPer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( nvarOut(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( nxyMax(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outGridNxy(nout,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( pointsFlag(nOut,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( pointsOut(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( pointsOutLand(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outDate(nout,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outDateFlag(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outEndPos(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outFilePer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outFileStep(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outNpWrite(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outPer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outPerNunits(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outRangeX(nOut,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outRangeY(nOut,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outSamPer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outStep(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outStepSamPer(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outTime(nout,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outTimeID(nout,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outUnit(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outWriteCount(nOut), stat=ierr ); ierrSum=ierrSum+ierr
!      allocate( snowCount(nOut,points), stat=ierr ); ierrSum=ierrSum+ierr
!      allocate( snowCountTile(nOut,points,ntiles), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( useCompressGrid(nOut), stat=ierr ); ierrSum=ierrSum+ierr

!     Reals.
      ALLOCATE( outGridDxy(nout,2), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outGridXY(nout,2), stat=ierr ); ierrSum=ierrSum+ierr

!     Logicals.
      ALLOCATE( havePFT(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( haveSCpool(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( haveSnow(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( haveSoil(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( haveTile(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( haveType(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outActivePrev(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outAreaLL(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outFirstActive(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outFirstSection(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outFirstWrite(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outLLorder(nout), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( rgProfile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( rpProfile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( snapProfile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outTemplate(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( tmeanProfile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outCompress(nOut), stat=ierr ); ierrSum=ierrSum+ierr

!     Characters.
      ALLOCATE( openedFileName(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outCtlFile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outDataFile(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( outName(nOut), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( compressGridFile(nOut), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_out 2' )

!     Integers.
      ALLOCATE( outVarID(nvarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varNlev(nvarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varNum(nvarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varStartPos(nvarOutTot), stat=ierr ); ierrSum=ierrSum+ierr

!     The following are allocated enough space to deal with the profile that has the most variables.
      iout = MAXVAL( nvarOut(1:nout) )  !  profile with most variables
      ALLOCATE( varPos(nout,iout), stat=ierr ); ierrSum=ierrSum+ierr

!     Characters and logicals.
      ALLOCATE( varDesc(nVarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varName(nVarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varType(nvarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( taccumVar(nVarOutTot), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( tmeanVar(nVarOutTot), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_out 3' )

      ALLOCATE( mapOutLand(nOut,pointsOutMax,2), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_out 4' )
      WRITE(*,*) 'Length of output vector=outLen=',outLen
      WRITE(*,*) 'Length of output vector used for file i/o=outLenWrite=',outLenWrite
      ALLOCATE( outVal(outLen), stat=ierr ); ierrSum=ierrSum+ierr

!     Allocate space used to store any extra "offline" diagnostics.
      IF ( useCiDiag ) ALLOCATE( ciDiag(land_pts,npft), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useGstomDiag ) ALLOCATE( gstomDiag(land_pts,npft), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useRdcDiag ) ALLOCATE( rdcDiag(land_pts,npft), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useRflowDiag ) ALLOCATE( rflowDiag(nxRoute,nyRoute), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useRoffInfDiag ) ALLOCATE( roffInfDiag(land_pts), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useRRunDiag ) ALLOCATE( rrunDiag(nxRoute,nyRoute), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useSnowGMeltDiag ) ALLOCATE( snowGMeltDiag(land_pts,ntiles), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useWfluxDiag ) ALLOCATE( wfluxDiag(land_pts,sm_levels), STAT=ierr ); ierrSum=ierrSum+ierr
      IF ( useWfluxSfcDiag ) ALLOCATE( wfluxSfcDiag(land_pts), STAT=ierr ); ierrSum=ierrSum+ierr

!     Initialise.
      IF ( ierrSum == 0 ) THEN
        IF ( useCiDiag ) ciDiag(:,:) = 0.0
        IF ( useGstomDiag ) gstomDiag(:,:) = 0.0
        IF ( useRdcDiag ) rdcDiag(:,:) = 0.0
        IF ( useRflowDiag ) rflowDiag(:,:) = 0.0
        IF ( useRoffInfDiag ) roffInfDiag(:) = 0.0
        IF ( useSnowGMeltDiag ) snowGMeltDiag(:,:) = 0.0
        IF ( useWfluxDiag ) wfluxDiag(:,:) = 0.0
        IF ( useWfluxSfcDiag ) wfluxSfcDiag(:) = 0.0
      ENDIF

!-------------------------------------------------------------------------------
    CASE ( 'init_out_map allocate_mapOut' )
      IF ( PRESENT(np) ) THEN
        ALLOCATE( mapOut(nout,np,2), stat=ierr ); ierrSum=ierrSum+ierr
      ELSE
        WRITE(*,*)'ERROR: allocate_arrays: argument np must be given for callName='  &
                 ,TRIM(callName)
        STOP
      ENDIF

!-------------------------------------------------------------------------------
    CASE ( 'init_out_map allocate_mapOutCompress' )
      IF ( PRESENT(np) ) THEN
        ALLOCATE( mapOutCompress(nout,np), stat=ierr ); ierrSum=ierrSum+ierr
      ELSE
        WRITE(*,*)'ERROR: allocate_arrays: argument np must be given for callName='  &
                 ,TRIM(callName)
        STOP
      ENDIF
!-------------------------------------------------------------------------------
    CASE ( 'init_out_map allocate_coordList' )
      IF ( PRESENT(np) ) THEN
        ALLOCATE( coordList(2,nout,np), stat=ierr ); ierrSum=ierrSum+ierr
      ELSE
        WRITE(*,*)'ERROR: allocate_arrays: argument np must be given for callName='  &
                 ,TRIM(callName)
        STOP
      ENDIF

!-------------------------------------------------------------------------------
    CASE ( 'init_out_varlist' )
      ALLOCATE( varNameList(nvar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varDescList(nvar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varTypeList(nvar), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( varUnitsList(nvar), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'spin_check' )
!     Spin up variables.
      ALLOCATE( spinValOld(nspinVar,npSpinMax,nzSpinMax), stat=ierr ); ierrSum=ierrSum+ierr

!-------------------------------------------------------------------------------
    CASE ( 'init_parms' )
!     Sea and sea-ice variables.
      ALLOCATE( SEA_FRAC(SSI_PTS), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( SICE_FRAC(SSI_PTS), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE( SICE_FRAC_NCAT(SSI_PTS,NICE), stat=ierr ); ierrSum=ierrSum+ierr
      
!-------------------------------------------------------------------------------
    CASE ( 'init_imogen' )
! Imogen variables
      ALLOCATE(T_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(RAINFALL_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(SNOWFALL_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(RH15M_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(UWIND_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(VWIND_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(DTEMP_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(PSTAR_HA_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(SW_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(LW_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(F_WET_CLIM(land_pts,mm), stat=ierr ); ierrSum=ierrSum+ierr
      
      ALLOCATE(T_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(CONV_RAIN_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(CONV_SNOW_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(LS_RAIN_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(LS_SNOW_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(QHUM_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(WIND_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(PSTAR_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(SW_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(LW_OUT(land_pts,mm,md,nsdmax), stat=ierr ); ierrSum=ierrSum+ierr
      
      ALLOCATE(LAT(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(LONG(land_pts), stat=ierr ); ierrSum=ierrSum+ierr
      ALLOCATE(DCTOT(land_pts), stat=ierr ); ierrSum=ierrSum+ierr

!     Initialise.
      IF ( ierrSum == 0 ) THEN
        T_CLIM(:,:)        = 0.0
        RAINFALL_CLIM(:,:) = 0.0
        SNOWFALL_CLIM(:,:) = 0.0
        RH15M_CLIM(:,:)    = 0.0
        UWIND_CLIM(:,:)    = 0.0
        VWIND_CLIM(:,:)    = 0.0
        DTEMP_CLIM(:,:)    = 0.0
        PSTAR_HA_CLIM(:,:) = 0.0
        SW_CLIM(:,:)       = 0.0
        LW_CLIM(:,:)       = 0.0
        F_WET_CLIM(:,:)    = 0.0

        T_OUT(:,:,:,:)         = 0.0
        CONV_RAIN_OUT(:,:,:,:) = 0.0
        CONV_SNOW_OUT(:,:,:,:) = 0.0
        LS_RAIN_OUT(:,:,:,:)   = 0.0
        LS_SNOW_OUT(:,:,:,:)   = 0.0
        QHUM_OUT(:,:,:,:)      = 0.0
        WIND_OUT(:,:,:,:)      = 0.0
        PSTAR_OUT(:,:,:,:)     = 0.0
        SW_OUT(:,:,:,:)        = 0.0
        LW_OUT(:,:,:,:)        = 0.0

        LAT(:)   = 0.0
        LONG(:)  = 0.0
        DCTOT(:) = 0.0
      ENDIF

!-------------------------------------------------------------------------------
    CASE default
       WRITE(*,*)'ERROR: allocate_arrays: do not recognise callName=',TRIM(callName)
       STOP
!-------------------------------------------------------------------------------
  END SELECT

! Check for error.
  IF ( ierr/=0 .OR. ierrSum /= 0 ) THEN
    IF ( ierrSum == 0 ) ierrSum = ierr
    CALL allocate_error( 'alloc',ierrSum,'allocate_arrays: '//TRIM(callName) )
    STOP
  ENDIF

  END SUBROUTINE allocate_arrays
!###############################################################################
!###############################################################################

  END MODULE alloc_mod
!###############################################################################

