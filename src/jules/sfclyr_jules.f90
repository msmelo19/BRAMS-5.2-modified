!############################# Change Log ##################################
! Interface entre BRAMS e JULES
! 

SUBROUTINE sfclyr_jules(mzp,mxp,myp,ia,iz,ja,jz,jdim)

  !--- Modulos do BRAMS ---
  use node_mod, only: &
       MYNUM,         & ! INTENT(OUT)
       i0,            & ! INTENT(IN)  
       j0	            ! INTENT(IN) 
  use io_params, only : frqanl, hfilout,frqhis,hfilin

  USE leaf_coms, ONLY : veg_ht,slcpd, nstyp, nvtyp_teb,slbs,slcons,slden,Wwilt,PHIsat &
                        ,alfa_vG,slmsts, dtll_factor,gzotheta 

  USE mem_leaf      

  use mem_carma,     ONLY: carma

  USE mem_jules      

  USE mem_basic,     ONLY: basic_g

  USE mem_grid ,     ONLY : dtlt,npatch, nzg, nzs, grid_g, time, dtlong, iyear1,imonth1, idate1 &
                            ,itime1, timmax,zt,dzt,ngrid,runtype   ! INTENT(IN)
                    
  USE rconstants,    ONLY : cpi,cp,alvl,vonk,g,p00,cpor

  USE mem_radiate,   ONLY : radiate_g,iswrtyp 

  USE mem_turb,      ONLY : turb_g
  
  USE nstypes,       ONLY : npft,ntype

  USE mem_cuparm,    ONLY: cuparm_g, nnqparm

  USE micphys,       ONLY: level

  USE mem_micro,     ONLY: micro_g
 
  USE chem1_list,    ONLY: CO2,chemical_mechanism
 
  USE mem_chem1,     ONLY: chem1_g, nsrc, chem1_src_vars, bioge, chem1_src_g,ntimes_src,chemistry

  USE mem_globaer,   ONLY: wtmol_air 

  !--- Modulos do JULES ---
  USE csigma, ONLY :  sbcon
  
  USE initial_mod, ONLY : dump_io

  USE init_grid_mod, ONLY :  init_grid

  USE inout,         ONLY : print_step,jinUnit,stdIn,outDir,runID,echo

  USE ancil_info,    ONLY : ntiles,sm_levels,land_pts,row_length,n_rows,land_index & ! INTENT(IN)
                            ,co2_dim_len,co2_dim_row,tile_pts,tile_index,frac

  USE p_s_parms, ONLY     : sthf,sthu   !   sthu may be used for sthuf

  USE prognostics,   ONLY : lai,gc,t_soil,tstar_tile,cs

  USE screen, ONLY : q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m

  USE output_mod,    ONLY : output
  
  USE spin_mod,      ONLY : ispin,spinUp

  USE time_loc,      ONLY : date,time_jules=>time

  USE coastal, ONLY :  tstar_sea

  USE time_mod,      ONLY : s_to_chhmmss

  USE trifctl,       ONLY : asteps_since_triffid, &
                            NPP,   &                   ! Net primary productivity (kg C/m2/s)
                            RESP_S, &                  ! Soil respiration (kg C/m2/s)
                            RESP_P, &                  ! Plant respiration (kg C/m2/s)
                            GPP                        ! Gross primary productivity (kg C/m2/s)

  USE update_mod,    ONLY : drive_update,veg_update

  USE veg_io_vars,   ONLY : vegVaryT

  USE fluxes,        ONLY : emis_tile, land_albedo,  &
                            LATENT_HEAT, &   !fluxo de calor latente do JULES
                            ftl_1, &       !fluxo de calor sensivel do JULES
                            TAUX_1, & !   W'ly component of surface wind stress (N/sq m)
                            TAUY_1    !   S'ly component of surface wind stress (N/sq m)

  USE switches, ONLY :  l_aggregate,routeOnly

  USE aero,          ONLY : co2_3d,u_s_std_tile

  IMPLICIT NONE
  
  INTEGER               :: nsoil, fat_dtlong
  
  LOGICAL  :: is_point
  
  INTEGER, INTENT(IN) :: mzp,mxp,myp,ia,iz,ja,jz,jdim
  
  !Local Variables
  INTEGER            :: ng,a_step,i ,j, k2,jd,l,n,ip,tB(ntype),qual_ip(ntype),it,k,mm

  LOGICAL            :: endRun    !  TRUE at end of last timestep in run

  LOGICAL              :: tem_in

  REAL               :: hcpi, zoverl,cx,wtol,psin,piv,ustar2,tstar2,rstar2,resp_s_tot,dz &
                        ,dz_inv,co2_flux,tempk1,tempk2,tempk3,fracliq,picpi,dens2,zts2,ths2 &
                        ,mass,AOD,diff_frac,rshort_diffuse(iz,jz)

  REAL, DIMENSION(:), ALLOCATABLE :: rlongupJ

  
  REAL, DIMENSION(:,:), ALLOCATABLE :: tile_frac   !  tile fractions

  
  CHARACTER(len=10)  ::  time_hms  ! Time in the form "hh:mm:ss H", used as a local variable
                                   ! that is required because some compilers can not cope
                                   ! with recursive write statements
  CHARACTER (LEN=80) :: veg(ntype)

  CHARACTER (LEN=2)  :: tB_str    

  REAL, DIMENSION(mxp,myp) :: rvs2,ups2,vps2,pp2,temp2,pcpgl

  INTEGER            :: nLAND,mxJ,myJ,nPEs,iam,ierr  ! =0 => possui ponto de continente na grade, =1=> nLAND=0

  LOGICAL            :: start=.TRUE.
  
  REAL, ALLOCATABLE :: zts2J(:,:)
  
  INTEGER, SAVE :: nts

  DATA wtol/1e-20/

  ng=ngrid

  mxJ=iz-ia+1
  myJ=jz-ja+1

  IF (MYNUM==2 .and. time==dtlong) WRITE(*,*) '-------- Executing with JULES... --------'
   
  ALLOCATE( zts2J(mxJ,myJ) )

!--- proximo do nascer e do por do sol o JULES serah chamado com maior frequencia, pois os processos de superficie variam mais rapidamente ---
!if (abs(radiate_g(ng)%cosz(4,4)) < 0.1) then
   fat_dtlong=1 
!elseif (abs(radiate_g(ng)%cosz(4,4)) < 0.2) then
!   fat_dtlong=3  ! no nascer e no por do sol o JULES serah chamado a cada timestep
!elseif (abs(radiate_g(ng)%cosz(4,4)) < 0.4) then
!   fat_dtlong=6  ! no nascer e no por do sol o JULES serah chamado a cada timestep
!else
!   fat_dtlong=10
!endif

   !--- Inicializando a grade do JULES ---
   IF (start) THEN     !DSM --- Semelhante a time==0, mas serve tambem para o history

      jinUnit = stdIn
      INQUIRE(FILE='./jules.in',EXIST=tem_in)
      IF (tem_in) THEN
         OPEN (77,FILE='./jules.in', STATUS='old')
      ELSE
         WRITE(*,*)'Not found ./jules.in control file.'
         STOP 
      ENDIF
 
      if ( echo ) WRITE(*,"(50('-'),/,a)")'Reading model control file (./jules.in)...'

      !-----------------------------------------------------------------------
      ! Read model options and misc other values.
      !-----------------------------------------------------------------------
      CALL init_opts(mxJ,myJ,nzg)

      !-----------------------------------------------------------------------
      ! Read date, time and location information
      !-----------------------------------------------------------------------
      CALL init_time(dtlong,iyear1,imonth1,idate1,itime1,timmax,time,frqhis,runtype)

      !-----------------------------------------------------------------------
      ! Read details of model grid and allocate arrays.
      !-----------------------------------------------------------------------

      CALL init_grid(mxJ,myJ,npatch,leaf_g(ng)%patch_area (   ia:iz, ja:jz, :),  &
                                        grid_g(ng)%glon         (   ia:iz, ja:jz   ),  &
                                        grid_g(ng)%glat         (   ia:iz, ja:jz   )      )
   ENDIF !(time==0)

   hcpi = .5 * cpi
   DO l=1,land_pts
      j = ( land_index(l)-1 ) / row_length + 1
      i = land_index(l) - ( j-1 ) * row_length

      IF (i<1 .or. i>mxJ .or. j<1 .or. j>myJ) THEN
         PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
         STOP
      ENDIF
      
      rvs2(i+ia-1,j+ja-1) = basic_g(ng)%rv(2,i+ia-1,j+ja-1)
      ups2(i+ia-1,j+ja-1) = basic_g(ng)%up(2,i+ia-1,j+ja-1)
      vps2(i+ia-1,j+ja-1) = basic_g(ng)%vp(2,i+ia-1,j+ja-1)
      
      picpi = (basic_g(ng)%pi0(2,i+ia-1,j+ja-1) + basic_g(ng)%pp(2,i+ia-1,j+ja-1)) * cpi
      pp2(i+ia-1,j+ja-1)=p00 * picpi ** cpor

      k2=int(grid_g(ng)%lpw(i+ia-1,j+ja-1))
      piv = hcpi * (basic_g(ng)%pi0(k2-1,i+ia-1,j+ja-1) + basic_g(ng)%pi0(k2,i+ia-1,j+ja-1)   &
                  + basic_g(ng)%pp(k2-1,i+ia-1,j+ja-1) + basic_g(ng)%pp(k2,i+ia-1,j+ja-1))
      
      temp2(i+ia-1,j+ja-1) = basic_g(ng)%theta(k2,i+ia-1,j+ja-1) * piv   != airtemp = tstar_tile
      
!DEBUG
!if (MYNUM==9600) then
!   write(15,*) 'time,i,j,temp2=',time,i+ia-1,j+ja-1,temp2(i+ia-1,j+ja-1)-273.15,basic_g(ng)%theta(2,i+ia-1,j+ja-1),basic_g(ng)%pi0(2,i+ia-1,j+ja-1)  !todos os ponto no oceano
!endif
!IF ( is_point(-44.86,-26.89,i+ia-1,j+ja-1,ng) ) THEN
!   write(15,*) 'time,i,j,lon,lat,temp2,theta,pi0,pp=>>>',time,i+ia-1,j+ja-1,grid_g(ng)%glon(i+ia-1,j+ja-1)  &
!               ,grid_g(ng)%glat(i+ia-1,j+ja-1), temp2(i+ia-1,j+ja-1)-273.15,basic_g(ng)%theta(2,i+ia-1,j+ja-1)  &
!               ,basic_g(ng)%pi0(2,i+ia-1,j+ja-1), basic_g(ng)%pp(2,i+ia-1,j+ja-1)
!ENDIF
      
      !--- Precipitacao total ---
      pcpgl(i+ia-1,j+ja-1)=0.
      IF (nnqparm(ng) > 0 .and. level >= 3) THEN
         pcpgl(i+ia-1,j+ja-1)=cuparm_g(ng)%conprr(i+ia-1,j+ja-1) + micro_g(ng)%pcpg(i+ia-1,j+ja-1)
      ELSEIF(nnqparm(ng) == 0 .and. level >= 3) THEN
         pcpgl(i+ia-1,j+ja-1)=micro_g(ng)%pcpg(i+ia-1,j+ja-1)
      ENDIF
      
      zts2J(i,j) = zt(2) * grid_g(ng)%rtgt(i+ia-1,j+ja-1)

   ENDDO


   IF (start) THEN     !DSM --- Semelhante a time==0, mas serve tambem para o history
      start=.false.

      CALL sfcdata

      a_step = 0
      
      !--- Inicializando o JULES ---
      CALL init( mxJ,myJ,npatch,nzg,nstyp,dtlong*fat_dtlong,iyear1,imonth1,idate1   &
                     ,itime1,timmax,leaf_g(ng)%patch_area   (   ia:iz, ja:jz, :)    &
                     ,grid_g(ng)%glon         (   ia:iz, ja:jz   )                  &
                     ,grid_g(ng)%glat         (   ia:iz, ja:jz   )                  &
                     ,leaf_g(ng)%veg_fracarea (   ia:iz, ja:jz, :)                  &
                     ,leaf_g(ng)%leaf_class   (   ia:iz, ja:jz, :)                  &
                     ,leaf_g(ng)%stom_resist  (   ia:iz, ja:jz, :)                  &
                     ,leaf_g(ng)%soil_water   (:, ia:iz, ja:jz, 2)                  &
                     ,leaf_g(ng)%soil_text    (:, ia:iz, ja:jz, 2)                  &
                     ,slmsts(:)                                                     &
                     ,leaf_g(ng)%veg_lai      (   ia:iz, ja:jz, :)                  &
                     ,temp2                   (   ia:iz, ja:jz   )                  &
                     ,leaf_g(ng)%soil_energy  (:, ia:iz, ja:jz, 1)                  &
                     ,leaf_g(ng)%soil_energy  (:, ia:iz, ja:jz, 2)                  &
                     ,veg_ht,slcpd,nvtyp_teb                                        &
                     ,slbs,slcons,slden,Wwilt,PHIsat,alfa_vG                        &
                     ,MYNUM,hfilout,hfilin,runtype, zts2J,slz                       &
               )
       
      !{--- Verificando opcoes do RAMSIN ---
      IF (NZG/=sm_levels) CALL erro_RAMSIN('NZG',len('NZG'),sm_levels,4)  ! 'nome da variavel descrita no RAMSIN', tamanho da variavel, valor definido em jules.in, valor recomendado

      !--- Abrindo os arquivos para imprimir as variaveis em ascii ---
      DO l=1,land_pts
         j = ( land_index(l)-1 ) / row_length + 1
         i = land_index(l) - ( j-1 ) * row_length
         IF (i<1 .or. i>mxJ .or. j<1 .or. j>myJ) THEN
            PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
            STOP
         ENDIF
         
         ustar2 = max(0.000001,sqrt( sqrt( (turb_g(ng)%sflux_u(i+ia-1,j+ja-1))**2 + (turb_g(ng)%sflux_v(i+ia-1,j+ja-1))**2 )))
         DO ip=1,npatch
            leaf_g(ng)%ustar(i+ia-1,j+ja-1,ip)=ustar2
         ENDDO
      ENDDO
      
      if (trim(runtype)=='HISTORY') CALL newTime( a_step,endRun,time,frqhis )  !DSM - Para ajustar a data da simulacao
   
   END IF !(if(start)
    
   !RMF! ELSEIF (time==dtlong .or. mod(time,dtlong*fat_dtlong)==0) THEN
      !{---Executando o JULES ---
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
         time_hms = s_to_chhmmss( time_jules )
      ENDIF


      !-------------------------------------------------------------------------------
      !   Update meteorological data.  2nd argument (next) is TRUE to indicate that 
      !   the "next" data in file are to be used.
      !-------------------------------------------------------------------------------

      !IF (co2_dim_len /= mxJ .or. co2_dim_row /= myJ) then
      !   PRINT*, 'ERRO... (co2_dim_len /= mxJ .or. co2_dim_row /= myJ):'     &
      !         , co2_dim_len,' /= ',mxJ, ' or ',co2_dim_row,' /= ',myJ
      !   STOP
      !ENDIF

      IF ( (chemical_mechanism .eq. 'RELACS_TUV' .or. chemical_mechanism .eq. 'CO2') &
            .and. CHEMISTRY == 0) THEN
         co2_3d(1:mxJ, 1:myJ) = chem1_g(CO2,ng)%sc_p (2, ia:iz, ja:jz) * 1.E-9
      ELSE
         co2_3d = 384.
      END IF
      
      !--------------------------- Calculo da radiacao difusa - Fornecido por Nilton Rosario -------------------------{
      
      !m=(cosz+0.50572*(6.07995+(90-acos(cosz)*180/pi))^(-1.6364))^-1 
      !z= ângulo zenital solar
      
      !AOD=profundidade óptica dos aerossóis no comprimento de onda 670 nm.
      !No rad_carma a AOD para 670 nm é o aotl(:,15)
      
      !Fração da radiação difusa (diff_frac)
      !m <= 1.1; diff_frac=0.0115*AOD^3 -0.1115*AOD^2+0.4693*AOD +0.1258
      !1.1<m <= 1.25; diff_frac=0.0129*AOD^3 -0.1235*AOD^2+0.4997*AOD +0.1304
      !1.25<m <= 1.4; diff_frac=0.0075*AOD^3 -0.1087*AOD^2+0.5035*AOD +0.1477
      !1.4<m <= 1.7; diff_frac=0.0052*AOD^3 -0.1031*AOD^2+0.5077*AOD +0.1795
      !1.7<m <= 2.0 diff_frac=0.0144*AOD^3 -0.1634*AOD^2+0.6207*AOD +0.1696
      !2.0<m <= 2.8 diff_frac=0.0166*AOD^3 -0.2237*AOD^2+0.7458*AOD +0.1851
      !m > 2.8 diff_frac=0.0736*AOD^3 -0.4631*AOD^2+1.0152*AOD +0.2005
      
      !rshort_diffuse=rshort*diff_frac
      !rshort_direct=rshort-rshort_diffuse
      
      DO l=1,land_pts
         j = ( land_index(l)-1 ) / row_length + 1
         i = land_index(l) - ( j-1 ) * row_length
         IF (i<1 .or. i>mxJ .or. j<1 .or. j>myJ .or. i+ia-1>iz .or. j+ja-1>jz) THEN
            PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
            STOP
         ENDIF
         
	 IF(iswrtyp == 4) THEN
	    
	    !- for CARMA only
	    mass=(( &
	   	   radiate_g(ng)%cosz(i+ia-1,j+ja-1)+0.50572* &
	   	   (6.07995+(90.-acos(radiate_g(ng)%cosz(i+ia-1,j+ja-1))* &
	   	   180/3.14159))**(-1.6364))**(-1))

            AOD = MAX(carma(ngrid)%aot(i+ia-1,j+ja-1,15),0.)
            
            IF ( mass <= 1.1 ) THEN
               diff_frac=0.0115*AOD**3 - 0.1115*AOD**2 + 0.4693*AOD + 0.1258
            ELSE IF (1.1 < mass .and. mass <= 1.25) THEN
               diff_frac=0.0129*AOD**3 - 0.1235*AOD**2 + 0.4997*AOD + 0.1304
            ELSE IF (1.25 < mass .and. mass <= 1.4) THEN
               diff_frac=0.0075*AOD**3 - 0.1087*AOD**2 + 0.5035*AOD + 0.1477
            ELSE IF (1.4 < mass .and. mass <= 1.7) THEN
               diff_frac=0.0052*AOD**3 - 0.1031*AOD**2 + 0.5077*AOD + 0.1795
            ELSE IF (1.7 < mass .and. mass <= 2.0) THEN
               diff_frac=0.0144*AOD**3 - 0.1634*AOD**2 + 0.6207*AOD + 0.1696
            ELSE IF (2.0 < mass .and. mass <= 2.8) THEN
               diff_frac=0.0166*AOD**3 - 0.2237*AOD**2 + 0.7458*AOD + 0.1851
            ELSE
               diff_frac=0.0736*AOD**3 - 0.4631*AOD**2 + 1.0152*AOD + 0.2005
            ENDIF	  
            
	ELSE
	
            diff_frac=0.  
	
	ENDIF
         !diff_frac=0.  !para rodar sem a parte difusa
	 
	 rshort_diffuse(i+ia-1, j+ja-1)=radiate_g(ng)%rshort(i+ia-1, j+ja-1)*diff_frac
        
      ENDDO

      CALL drive_update( mxJ,myJ,npatch,nzg                          &
                        ,rshort_diffuse       (ia:iz, ja:jz)         &
                        ,radiate_g(ng)%rshort (ia:iz, ja:jz)         &
                        ,radiate_g(ng)%rlong  (ia:iz, ja:jz)         &
                        ,temp2                (ia:iz, ja:jz)         &
                        ,ups2                 (ia:iz, ja:jz)         &
                        ,vps2                 (ia:iz, ja:jz)         &
                        ,pp2                  (ia:iz, ja:jz)         &
                        ,rvs2                 (ia:iz, ja:jz)         &
                        ,pcpgl                (ia:iz, ja:jz)         &
                        ,a_step,.TRUE.                               &
                       )


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
      !   Generate output. This call (at the end of the timestep) is expected to
      !   generate most of the output for a run.  The logical argument (endCall) is 
      !   always true at this call.
      !-------------------------------------------------------------------------------
      CALL output( a_step,.TRUE. )

      !-------------------------------------------------------------------------------
      !   Update time.
      !-------------------------------------------------------------------------------
      CALL newTime( a_step,endRun,time,frqhis )
      if (time==timmax-dtlong) CALL newTime( a_step,endRun,time+dtlong,frqhis )  !DSM - Para imprimir o ultimo history (dump)
   
      !--- Compatibilizando as variaveis do BRAMS com a do JULES ---
      CALL brams2jules(veg,ntype)
   
      dtll_factor = 1. / float(max(1,nint(dtlt/40.+.4)))

      !-------------------------------------------------------------------------------
      
      ALLOCATE(tile_frac(land_pts,ntiles))
      ALLOCATE(rlongupJ(land_pts))
      
      
      ! Calculate tile fractions.
      IF ( .NOT. routeOnly ) THEN
        tile_frac(:,:) = 0.0
        IF ( l_aggregate ) THEN
          tile_frac(:,1) = 1.0
        ELSE
          DO n=1,ntype
            DO j=1,tile_pts(n)
              i = tile_index(j,n)
              tile_frac(i,n) = frac(i,n)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      !-------------------------------------------------------------------------------
      
      !--- Acoplando o rlongup ---{
      call gbmTileDiag2( emis_tile * tstar_tile*tstar_tile*tstar_tile*tstar_tile,tile_frac,rlongupJ )
      rlongupJ = sbcon * rlongupJ
      !-----}
      
      !--- Escrevendo o vento de hora em hora ---
      IF (MOD(time,3600.*6+dtlong)==0 .or. time<=dtlong) THEN
         nts=1
      ELSE
         IF (MOD(time,3600.+dtlong)==0 .or. MOD(time,3600.*6+2*dtlong)==0 .or. time==2*dtlong)  nts=nts+1
      END IF

      DO l=1,land_pts
         j = ( land_index(l)-1 ) / row_length + 1
         i = land_index(l) - ( j-1 ) * row_length
         IF (i<1 .or. i>mxJ .or. j<1 .or. j>myJ) THEN
            PRINT*, "ERRO... conversao incorreta de l para i,j -> l,i,j=",l,i,j
            STOP
         ENDIF
         
         !--- LAI(land_pts,npft) para veg_lai(i,j,ip) ---
         DO n=1,npft
            READ(veg(n)(1:2),*) tB(n)  !caso este ponto (i,j) nao possua nenhum representante para este tile (k) serah utilizado o primeiro tipo definido em brams2jules

            qual_ip(n)=-999

            DO ip=1,npatch
               WRITE(tB_str,'(i2.2)') nint(leaf_g(ng)%leaf_class(i+ia-1,j+ja-1,ip))  ! Jogando o tipo do BRAMS em tB_str
               IF (index(veg(n),tB_str)/=0) THEN
                  READ(tB_str,*) tB(n)
                  qual_ip(n)=ip   
                  EXIT  ! jah encontrou para o menor patch (mais representativo)
               ENDIF
            ENDDO
            
            IF (qual_ip(n) /= -999 ) THEN               
               !--- Acoplando o lai ---
               leaf_g(ng)%veg_lai(i+ia-1,j+ja-1,qual_ip(n)) = lai(l,n)  !OK
               !!!!!veg_tai = veg_lai + sai(nveg) + dead_lai
            ENDIF
         ENDDO

         !--- gc(land_pts,ntiles) para 1/stom_resist(i,j,ip) ---
         DO n=1,ntiles
         
            READ(veg(n)(1:2),*) tB(n)  !caso este ponto (i,j) nao possua nenhum representante para este tile (k) serah utilizado o primeiro tipo definido em brams2jules

            qual_ip(n)=-999

            DO ip=1,npatch
               WRITE(tB_str,'(i2.2)') nint(leaf_g(ng)%leaf_class(i+ia-1,j+ja-1,ip))  ! Jogando o tipo do BRAMS em tB_str
               IF (index(veg(n),tB_str)/=0) THEN
                  READ(tB_str,*) tB(n)
                  qual_ip(n)=ip   
                  EXIT  ! jah encontrou para o menor patch (mais representativo)
               ENDIF
            ENDDO
         
            IF (qual_ip(n) /= -999 .and. qual_ip(n) /= 1) THEN    
               !--- Acoplando a stom_resist ---
               leaf_g(ng)%stom_resist(i+ia-1,j+ja-1,qual_ip(n)) = 1/max(0.000001,gc(l,n))  !OK

               !--- Acoplando a teveg ---
               leaf_g(ng)%veg_temp(i+ia-1,j+ja-1,qual_ip(n))=tstar_tile(l,n)
               !jules_g(ng)%u_s(i+ia-1,j+ja-1,qual_ip(n))=u_s_std_tile(l,n)  ! DSM 11/04/2014
            ENDIF

         ENDDO
         
         !--- Acoplando o rlongup ---{
         radiate_g(ng)%rlongup(i+ia-1,j+ja-1) = rlongupJ(l)
         !-----}
         
         !--- Acoplando o albedt ---{
         radiate_g(ng)%albedt(i+ia-1,j+ja-1)=(land_albedo(i,j,1)+land_albedo(i,j,2)+land_albedo(i,j,3)+land_albedo(i,j,4))/4
         !-----}
         
         !--- Escrevendo as variaveis do JULES de output ---
         jules_g(ng)%u10mj(i+ia-1,j+ja-1)=u10m(i,j)
         jules_g(ng)%v10mj(i+ia-1,j+ja-1)=v10m(i,j)
         jules_g(ng)%t2mj(i+ia-1,j+ja-1)=t1p5m(i,j)
         jules_g(ng)%rv2mj(i+ia-1,j+ja-1)=q1p5m(i,j)
         jules_g(ng)%csj(i+ia-1,j+ja-1)=cs(l,1)
         jules_g(ng)%u10mj1hr(nts,i+ia-1,j+ja-1)=u10m(i,j)
         jules_g(ng)%v10mj1hr(nts,i+ia-1,j+ja-1)=v10m(i,j)
	 jules_g(ng)%gpp(i+ia-1,j+ja-1)=gpp(l)
	 jules_g(ng)%resp_p(i+ia-1,j+ja-1)=resp_p(l)
	 jules_g(ng)%npp(i+ia-1,j+ja-1)=npp(l)
	 jules_g(ng)%resp_s(i+ia-1,j+ja-1)=sum(resp_s(l,:))
                  
         IF ( (chemical_mechanism .eq. 'RELACS_TUV' .or. chemical_mechanism .eq. 'CO2') &
             .and. CHEMISTRY >= 0) THEN
            resp_s_tot = SUM( resp_s(l,:) )
            dz        = grid_g(ng)%rtgt(i+ia-1,j+ja-1)/dzt(2) ! dzt=1/(z(k)-z(k-1))	
            dz_inv    = 1./dz
            co2_flux=(resp_s_tot - npp(l)) * 44./12. * dz_inv * 1E+9  ! kg[CO2]/(m^3 seg)
         
            DO it=1,ntimes_src(bioge) 
               chem1_src_g(it,bioge,CO2,ng)%sc_src(1,i+ia-1,j+ja-1) = co2_flux
            ENDDO
         END IF

         DO k=1,sm_levels
            
            !--- Acoplando a umidade do solo ---
            nsoil = nint(leaf_g(ng)%soil_text(k,i+ia-1,j+ja-1,2))
            leaf_g(ng)%soil_water(k,i+ia-1,j+ja-1,2:npatch)=sthu(l,sm_levels-k+1)*slmsts(nsoil)
         ENDDO
        
            !--- Acoplando fluxo de calor latente ---
            turb_g(ng)%sflux_r(i+ia-1,j+ja-1) = latent_heat(i,j)/alvl
            
            !--- Acoplando fluxo de calor sensivel ---
            turb_g(ng)%sflux_t(i+ia-1,j+ja-1) = ftl_1(i,j)/cp

            dens2 = (basic_g(ng)%dn0(1,i+ia-1,j+ja-1) + basic_g(ng)%dn0(2,i+ia-1,j+ja-1)) * .5

            !--- Acoplando fluxo sflux_u e sflux_v ---
            turb_g(ng)%sflux_u(i+ia-1,j+ja-1) = -1*taux_1(i,j)/dens2

            turb_g(ng)%sflux_v(i+ia-1,j+ja-1) = -1*tauy_1(i,j)/dens2
        
            !--- Acoplando ustar, tstar, rstar ---
            ustar2 = max(0.000001,sqrt( sqrt( (turb_g(ng)%sflux_u(i+ia-1,j+ja-1))**2 + (turb_g(ng)%sflux_v(i+ia-1,j+ja-1))**2 )))
            tstar2 = ftl_1(i,j)/(dens2 * cp * ustar2)
            rstar2 = latent_heat(i,j)/(dens2 * alvl * ustar2)

            DO ip=1,npatch
               leaf_g(ng)%ustar(i+ia-1,j+ja-1,ip)=ustar2
               leaf_g(ng)%tstar(i+ia-1,j+ja-1,ip)=tstar2
               leaf_g(ng)%rstar(i+ia-1,j+ja-1,ip)=rstar2
            ENDDO
            
            !--- fluxo sflux_w ---
            zts2 = zt(2) * grid_g(ng)%rtgt(i+ia-1,j+ja-1)
            ths2 = basic_g(ng)%theta(2,i+ia-1,j+ja-1)

            gzotheta = g * zts2 / ths2
            zoverl = gzotheta * vonk * tstar2 / (ustar2 * ustar2)

            IF (zoverl < 0.) THEN
               cx = zoverl * sqrt(sqrt(1. - 15. * zoverl))
            ELSE
               cx = zoverl / (1.0 + 4.7 * zoverl)
            ENDIF
            
            psin = sqrt((1.-2.86 * cx) / (1. + cx * (-5.39 + cx * 6.998 )))

            !--- Acoplando sflux_w ---
            ip=2
            turb_g(ng)%sflux_w(i+ia-1,j+ja-1) = (0.27 * max(6.25 * (1. - cx) * &
              psin,wtol)- 1.18 * cx * psin) * ustar2 * ustar2 * leaf_g(ng)%patch_area(i+ia-1,j+ja-1,ip)    

      ENDDO

  !RMF! ENDIF  !IF (time==0)

!DEBUG
!if (MYNUM==9600) then
!   write(15,*) 'time,temp2=',time,temp2-273.15  !todos os ponto no oceano
!   write(15,*) 'time,t2mj=',time,jules_g(ng)%t2mj-273.15  !todos os ponto no oceano
!endif
   


   call leaf_bcond(mxp,myp,nzg,nzs,npatch,jdim                &    
        ,leaf_g(ng)%soil_water   ,leaf_g(ng)%sfcwater_mass    &
        ,leaf_g(ng)%soil_energy  ,leaf_g(ng)%sfcwater_energy  &
        ,leaf_g(ng)%soil_text    ,leaf_g(ng)%sfcwater_depth   &
        ,leaf_g(ng)%ustar        ,leaf_g(ng)%tstar            &
        ,leaf_g(ng)%rstar        ,leaf_g(ng)%veg_albedo       &
        ,leaf_g(ng)%veg_fracarea ,leaf_g(ng)%veg_lai          &
        ,leaf_g(ng)%veg_tai                                   &
        ,leaf_g(ng)%veg_rough    ,leaf_g(ng)%veg_height       &
        ,leaf_g(ng)%patch_area   ,leaf_g(ng)%patch_rough      &
        ,leaf_g(ng)%patch_wetind ,leaf_g(ng)%leaf_class       &
        ,leaf_g(ng)%soil_rough   ,leaf_g(ng)%sfcwater_nlev    &
        ,leaf_g(ng)%stom_resist  ,leaf_g(ng)%ground_rsat      &
        ,leaf_g(ng)%ground_rvap  ,leaf_g(ng)%veg_water        &
        ,leaf_g(ng)%veg_temp     ,leaf_g(ng)%can_rvap         &
        ,leaf_g(ng)%can_temp     ,leaf_g(ng)%veg_ndvip        &
        ,leaf_g(ng)%veg_ndvic    ,leaf_g(ng)%veg_ndvif        )

   DEALLOCATE(tile_frac)
   DEALLOCATE(rlongupJ)
  
   RETURN
   
END SUBROUTINE sfclyr_jules


!*****************************************************************************
!--- Esta subrotina tem como finalidade imprimir mensagem caso o RAMSIN nao estiver de acordo
!--- Com o modelo JULES.
SUBROUTINE erro_RAMSIN(varName,tam,varJULES,varRecomend)
   IMPLICIT NONE
   CHARACTER (LEN=*),  INTENT(IN) :: varName
   INTEGER,             INTENT(IN) :: tam, varJULES,varRecomend
   
   PRINT*, '*********************************************************************************'
   PRINT*, '*********************************************************************************'
   PRINT*
   PRINT*, '    Valor setado em "jules.in" deve ser igual ao do RAMSIN, portanto setar: '
   PRINT*
   PRINT*, '               '//varName(1:tam),' = ',varJULES
   PRINT*
   
   IF (varJULES/=varRecomend) THEN
      PRINT*
      PRINT*,  '     Atencao!!! o valor recomendado eh: '//varName(1:tam),' = ',varRecomend
      PRINT*
      PRINT*
   ENDIF
   
   PRINT*, '*********************************************************************************'
   PRINT*, '*********************************************************************************'

   STOP
   RETURN
END SUBROUTINE erro_RAMSIN


!################################################################################
!################################################################################
!################################################################################
  SUBROUTINE gbmTileDiag2( inval,tile_frac,outval )

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
    land_pts,ntiles  &
!  imported arrays with intent(in)
   ,tile_index,tile_pts

!-------------------------------------------------------------------------------
  IMPLICIT NONE

  REAL ::   &!  function result
    outval(land_pts)  !  function result

  REAL, INTENT(in) ::  &!  in arrays
    inval(land_pts,ntiles)     &!  input tile field.
   ,tile_frac(land_pts,ntiles)  !  tile fractions

  INTEGER ::  &!  local SCALARS
    i,p,t  !  work/loop counters

  LOGICAL :: tmask(ntiles) !  mask indicating which tiles are to be included in average.

!-------------------------------------------------------------------------------
    tMask(:) = .TRUE.

! Initialise the average.
  outval(:) = 0.0

  DO t=1,ntiles
    IF ( tMask(t) ) THEN
      DO i=1,tile_pts(t)
        p = tile_index(i,t)
        outval(p) = outval(p) + tile_frac(p,t)*inval(p,t)
      ENDDO
    ENDIF
  ENDDO

  END SUBROUTINE  gbmTileDiag2


!*****************************************************************************
!--- Dado um lon e um lat, esta funcao tem como objetivo verificar se o ponto
!--- x,y corresponde a esta posicao. Utilizada para imprimir sempre o mesmo
!--- ponto independente do numero de processadores utilizado.
FUNCTION is_point(lon,lat,i,j,ng)
   USE mem_grid ,     ONLY : grid_g,deltaxn,deltayn

   IMPLICIT NONE
   
   LOGICAL  :: is_point
   REAL     :: lon,lat,delta
   INTEGER  :: i,j,ng
   is_point=.false.
   
   if (abs( grid_g(ng)%glon(i,j)-lon ) <= abs( (grid_g(ng)%glon(i,j)-grid_g(ng)%glon(i+1,j) )/2. ) .and. &
       abs( grid_g(ng)%glat(i,j)-lat ) <= abs( (grid_g(ng)%glat(i,j)-grid_g(ng)%glat(i,j+1) )/2. )            ) is_point=.true.
   RETURN
END FUNCTION is_point
