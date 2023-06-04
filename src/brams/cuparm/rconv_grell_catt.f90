!--------------------------------------------------------------------------------
! Cumulus Parameterization by G. Grell
!--------------------------------------------------------------------------------

subroutine CUPARM_GRELL_CATT(iens)

  use mem_basic         , only: basic_g
  use mem_tend          , only: tend
  use mem_cuparm        , only: confrq,cuparm_g,cuparm_g_sh
  use node_mod          , only: mynum,   &   ! INTENT(IN)
       mxp,     &   ! INTENT(IN)
       myp,     &   ! INTENT(IN)
       mzp,     &   ! INTENT(IN)
       ia,      &   ! INTENT(IN)
       iz,      &   ! INTENT(IN)
       ja,      &   ! INTENT(IN)
       jz,      &   ! INTENT(IN)
       i0,      &   ! INTENT(IN)  
       j0	     ! INTENT(IN) 
  use mem_grid          , only: time,    &   ! INTENT(IN)
       initial, &   ! INTENT(IN)
       dtlt,	 &   ! INTENT(IN)
       itime1,  &   ! INTENT(IN)
       ngrid,   &   ! INTENT(IN)
       grid_g,  &   ! INTENT(IN)  	  
       dtlongn, &   ! INTENT(IN)
       deltaxn, &   ! INTENT(IN)
       deltayn, &
       npatch


  use rconstants        , only: tkmin 
  !use extras            , only: extra3d,extra2d,na_EXTRA3D 
  use mem_turb          , only: turb_g
  use mem_micro         , only: micro_g
  use mem_scratch       , only: scratch
  use mem_scalar        , only: scalar_g
  !srf
  use io_params         , only: frqanl
  use mem_leaf          , only: leaf_g
  use micphys           , only: level

  !- use modules for Grell Parameterization
  use mem_grell_param   , only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp
  use mem_scratch1_grell, only: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,kstabi4d, &
       kstabm4d,xmb4d,edt4d,enup5d,endn5d,deup5d,dedn5d,zup5d,    &
       zdn5d,iruncon,zcup5d,pcup5d,prup5d,clwup5d,tup5d
  use mem_grell         , only: cuforc_g,cuforc_sh_g,grell_g,grell_g_sh
 

    
  use mem_radiate, only: ISWRTYP, ILWRTYP ! Intent(in)

!--(DMK-CCATT-INI)-----------------------------------------------------
  use ccatt_start, only:        &
       CCATT           ! intent(in)

!ML -- In case you want to output massflux
  use mem_stilt         , only: imassflx
!--(DMK-CCATT-OLD)-----------------------------------------------------
!  use catt_start, only:        &
!       CATT           ! intent(in)
!--(DMK-CCATT-FIM)-----------------------------------------------------
  
  implicit none
  include "i8.h"
  integer, intent(IN) :: iens
  integer,parameter :: CPTIME = 0. !orig: CPTIME = 7200.

  integer,parameter :: i_forcing = 1
  
  integer,parameter :: trigg = 1 ! trigg=1 aciona o gatilho original (kbcon) 
                                 ! trigg=2 aciona o gatilho do ecmwf

  integer,parameter :: autoconv = 1 ! =1, Kessler
                                    ! =2, Berry (must use trigg=2)
  
  integer,parameter :: do_cupar_mcphys_coupling = 1 ! direct link cupar-microphysics
  ! =0 , no coupling
  !------------------------ deep convection --------------------------------------------
  if(iens == 1) then 
     !
     !      Zero out tendencies initially
     if (TIME.eq.0.) then 
        cuparm_g(ngrid)%thsrc=0.
        cuparm_g(ngrid)%rtsrc=0.
		cuparm_g(ngrid)%clsrc=0.
		end if 

     if(INITIAL.eq.2.and.TIME.lt.CPTIME-dtlt) return

     if(mod(TIME,CONFRQ).lt.DTLT.or.time.lt. .01 .or.abs(time-cptime) .lt. 0.01) then !002

        iruncon=1
       
	IF( (ISWRTYP /= 4 .or. ILWRTYP /= 4) .and. AUTOCONV == 2 )  STOP ' berry formulation needs carma radiation'
        
        cuparm_g(ngrid)%thsrc=0.
        cuparm_g(ngrid)%rtsrc=0.
        cuparm_g(ngrid)%conprr=0.
        cuparm_g(ngrid)%clsrc=0.
       
        !srf - use the old way to define the cumulus forcing
	if(i_forcing /= 1) then
             call atob(mxp * myp * mzp,tend%THT(1),cuforc_g(ngrid)%lsfth(1,1,1))
             call atob(mxp * myp * mzp,tend%RTT(1),cuforc_g(ngrid)%lsfrt(1,1,1))
        endif	

        call cuparth_catt(CCATT,           &
             mynum,                        & !01
             mgmxp,                        & !02
             mgmyp,                        & !03
             mgmzp,                        & !04
             mzp,                          & !05
             mxp,                          & !06
             myp,                          & !07
             ia,                           & !08
             iz,                           & !09
             ja,                           & !10
             jz,                           & !11
             i0,                           & !12
             j0,                           & !13
             maxiens,                      & !15
             iens,                         & !16 
             ngrid,                        & !17
             ngrids_cp,                    & !18
             DTLT,                         & !19
             time,                         & !20
             basic_g(ngrid)%up   (1,1,1),  & !21
             basic_g(ngrid)%vp   (1,1,1),  & !22
             basic_g(ngrid)%wp   (1,1,1),  & !23
             basic_g(ngrid)%theta(1,1,1),  & !24
             basic_g(ngrid)%thp  (1,1,1),  & !24
             basic_g(ngrid)%pp   (1,1,1),  & !25
             basic_g(ngrid)%pi0  (1,1,1),  & !26
             basic_g(ngrid)%dn0  (1,1,1),  & !27
             basic_g(ngrid)%rv   (1,1,1),  & !28
             turb_g(ngrid)%tkep  (1,1,1),  & !29
             tkmin,                        & !30
	     micro_g(ngrid)%rcp(1,1,1),    &! liquid water
    	     micro_g(ngrid)%rrp(1,1,1),    &! pristine
    	     micro_g(ngrid)%rpp(1,1,1),    &
	     micro_g(ngrid)%rsp(1,1,1),    &
    	     micro_g(ngrid)%rap(1,1,1),    &
	     micro_g(ngrid)%rgp(1,1,1),    &
    	     micro_g(ngrid)%rhp(1,1,1),    &
!
	     grid_g(ngrid)%topt  (1,1),    & !29
             grid_g(ngrid)%RTGT  (1,1),    & !30
             !
             cuforc_g(ngrid)%lsfth(1,1,1)  ,& !33 
             cuforc_g(ngrid)%lsfrt(1,1,1)  ,& !34 
             tend%PT(1),                   & !35
             cuparm_g(ngrid)%THSRC (1,1,1),& !36 
             cuparm_g(ngrid)%RTSRC (1,1,1),& !37 
             cuparm_g(ngrid)%CLSRC (1,1,1),& !37 
             cuparm_g(ngrid)%CONPRR(1,1),  & !38      
!
!             extra3d(5,ngrid)%d3   (1,1,1),& !39 ! cloud/ice tendency
!             extra3d(1,ngrid)%d3   (1,1,1),& !39 ! ensemble output
!             extra3d(4,ngrid)%d3   (1,1,1),& !39 ! ensemble output
!             extra3d(9,ngrid)%d3   (1,1,1),& !39 ! perfil vertical de W do novo trigger (claudio Silva)
!             extra2d(3,ngrid)%d2   (1,1),  & !39 
!
             ierr4d,                       & !40
             jmin4d,                       & !41
             kdet4d,                       & !42
             k224d,                        & !43
             kbcon4d,                      & !44
             ktop4d,                       & !45
             kpbl4d,                       & !46
             kstabi4d,                     & !47
             kstabm4d,                     & !48
             xmb4d,                        & !49
             edt4d,                        & !50
             !
             zcup5d,                       & !51 
             pcup5d,                       & !52 
             enup5d,                       & !53
             endn5d,                       & !54
             deup5d,                       & !55 
             dedn5d,                       & !56
             zup5d,                        & !57
             zdn5d ,                       & !58
             prup5d,                       & !59
             clwup5d,                      & !60
             tup5d ,                       & !61
             !
             grell_g(ngrid)%upmf   (1,1),  & !62
             grell_g(ngrid)%dnmf   (1,1),  & !63
             grell_g(ngrid)%xierr  (1,1),  & !64
             grell_g(ngrid)%xktop  (1,1),  & !65
             grell_g(ngrid)%xkbcon (1,1),  & !66
             grell_g(ngrid)%xk22   (1,1),  & !67             
             grell_g(ngrid)%xjmin  (1,1),  & !68
             grell_g(ngrid)%xkdt   (1,1),  & !69
             grell_g(ngrid)%xiact_p(1,1),  & !70
             grell_g(ngrid)%xiact_c(1,1),  & !71
	     confrq,frqanl,                &
	     deltaxn(ngrid)*deltayn(ngrid),&
             leaf_g(ngrid)%patch_area(1,1,1), &
	     npatch,                       &
	     level,                        &
	     grid_g(ngrid)%glat     (1,1), & 
	     grid_g(ngrid)%glon     (1,1), & !
  	     turb_g(ngrid)%sflux_r  (1,1), & ! fluxos a serem usados em trigg_ecmwf
             turb_g(ngrid)%sflux_t  (1,1), &
	     trigg,autoconv )    	       

!--(DMK-CCATT-INI)------------------------------------------------------------
! [ML------------- Stilt - RAMS coupling  ------------------
      if (imassflx == 1) then
        call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz             &
           ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp                  &    
           ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d           &
           ,kstabi4d,kstabm4d,xmb4d,edt4d                              &
           ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d      & 
           ,iens)
      end if
! ------------- Stilt - RAMS coupling  ------------------ ML]
!--(DMK-CCATT-FIM)------------------------------------------------------------

     end if 

     call accum(int(mxp*myp*mzp,i8), tend%tht(1), cuparm_g(ngrid)%thsrc(1,1,1))
     call accum(int(mxp*myp*mzp,i8), tend%rtt(1), cuparm_g(ngrid)%rtsrc(1,1,1))

     call update(mxp*myp, cuparm_g(ngrid)%aconpr(1,1),cuparm_g(ngrid)%conprr(1,1),dtlt)

     if(do_cupar_mcphys_coupling == 1) then
       call cupar2mcphysics(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt,& 
                            cuparm_g(ngrid)%clsrc,    &
                            basic_g(ngrid)%theta,    & 
                            basic_g(ngrid)%pp,    & 
                            basic_g(ngrid)%pi0     )
     endif
     
  endif 

  !---------------   Shallow cumulus scheme -----------------------------------------
  if(iens == 2) then !006

	if(TIME.eq.0.) then !004
          cuparm_g_sh(ngrid)%thsrc=0.
          cuparm_g_sh(ngrid)%rtsrc=0.
        end if 

       if(INITIAL.eq.2.and.TIME.lt.CPTIME-dtlt) return
       if(mod(TIME,CONFRQ).lt.DTLT.or.time.lt. .01 .or. abs(time-cptime).lt. 0.01) then !005
        iruncon=1

        cuparm_g_sh(ngrid)%thsrc=0.
        cuparm_g_sh(ngrid)%rtsrc =0.           
        
        !srf - use the old way to define the cumulus forcing
	if(i_forcing /= 1) then
             call atob(mxp * myp * mzp,tend%THT(1),cuforc_sh_g(ngrid)%lsfth(1,1,1))
             call atob(mxp * myp * mzp,tend%RTT(1),cuforc_sh_g(ngrid)%lsfrt(1,1,1))
        endif	


        call cuparth_shal(CCATT,             &
             mynum,        		     &   
             mgmxp,        		     &   
             mgmyp,        		     &   
             mgmzp,        		     &   
             mzp,          		     &   
             mxp,          		     &   
             myp,          		     &   
             ia,           		     &   
             iz,           		     &   
             !
             ja,           		     &   
             jz,           		     &   
             i0,           		     &   
             j0,           		     &   
             maxiens,      		     &   
             iens,         		     &   
             ngrid,        		     &   
             ngrids_cp,    		     &   
             dtlt,         		     &   
             time,         		     &   
             !
             basic_g(ngrid)%up   (1,1,1),    &   
             basic_g(ngrid)%vp   (1,1,1),    &   
             basic_g(ngrid)%wp   (1,1,1),    &   
             basic_g(ngrid)%theta(1,1,1),    &   
             basic_g(ngrid)%thp  (1,1,1),    &   
             basic_g(ngrid)%pp   (1,1,1),    &   
             basic_g(ngrid)%pi0  (1,1,1),    &   
             basic_g(ngrid)%dn0  (1,1,1),    &   
             basic_g(ngrid)%rv   (1,1,1),    &   
             turb_g(ngrid)%tkep  (1,1,1),    &   
             !
             tkmin,                          &   
	     micro_g(ngrid)%rcp(1,1,1),      &! liquid water
!
             grid_g(ngrid)%topt     (1,1),   &   
             grid_g(ngrid)%RTGT     (1,1),   &   
             cuforc_sh_g(ngrid)%lsfth(1,1,1), &   
             cuforc_sh_g(ngrid)%lsfrt(1,1,1), &   
             tend%PT(1),		     &   
             cuparm_g_sh(ngrid)%thsrc(1,1,1),&   
             cuparm_g_sh(ngrid)%rtsrc(1,1,1),&   
!
!             extra3d(2,ngrid)%d3     (1,1,1),&   !39 !<< usando extra3d(2)
!             extra2d(2,ngrid)%d2     (1,1),  &   !39 !<< usando extra2d(2)
!
             !
             ierr4d,                         & 
             jmin4d,                         & 
             kdet4d,                         & 
             k224d,                          & 
             kbcon4d,                        & 
             ktop4d,                         & 
             kpbl4d,                         & 
             kstabi4d,                       & 
             kstabm4d,                       & 
             !
             xmb4d,                          & 
             edt4d,                          & 
             zcup5d,                         & 
             pcup5d,                         & 
             enup5d,                         & 
             endn5d,                         & 
             deup5d,                         & 
             dedn5d,                         & 
             zup5d,                          & 
             zdn5d,                          & 
             prup5d,                         & 
             clwup5d,                        & 
             tup5d,                          & 
             !
             grell_g_sh(ngrid)%upmf  (1,1),  & 
             grell_g_sh(ngrid)%xierr (1,1),  & 
             grell_g_sh(ngrid)%xktop (1,1),  & 
             grell_g_sh(ngrid)%xkbcon(1,1),  & 
             grell_g_sh(ngrid)%xk22  (1,1),  & 
!            grell_g   (ngrid)%xierr (1,1),  & !para uso futuro, inibir shallow se deep is ON
             grell_g_sh   (ngrid)%xierr (1,1),  & 
	     confrq,frqanl,                  &
	     deltaxn(ngrid)*deltayn(ngrid),  &
             leaf_g(ngrid)%patch_area(1,1,1),&
	     npatch,                         &
	     level,                          &
	     trigg,autoconv )    

!
! [ML------------- Stilt - RAMS coupling  ------------------
       if (imassflx == 1) then
         call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz           &
          ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp                  &    
          ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d           &
          ,kstabi4d,kstabm4d,xmb4d,edt4d                              &
          ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d      & 
          ,iens)
       end if
! ------------- Stilt - RAMS coupling  ------------------ ML]
!
     end if 

     call accum(int(mxp*myp*mzp,i8), tend%tht(1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     call accum(int(mxp*myp*mzp,i8), tend%rtt(1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))
	end if


  !--------- Convective Transport based on mass flux scheme ------------------------------------
  if(CCATT == 1 .and. iruncon == 1) then
     scratch%scr1(:)=0.
     call trans_conv_mflx(iens,scratch%scr1(1))
  end if

end subroutine CUPARM_GRELL_CATT

!------------------------------------------------------------
subroutine cupar2mcphysics(m1,m2,m3,ia,iz,ja,jz,ngrid,dtlt & 
                          ,clsrc                &
                          ,theta,pp,pi0)
  use micphys     ,only: level,mcphys_type
  use mem_micro   ,only: micro_g
  use mem_tend    ,only: tend
  use rconstants  ,only: cpi
  implicit none
  integer m1,m2,m3,ia,iz,ja,jz,k,i,j,ngrid
  real dtlt
  real, dimension(m1,m2,m3) :: theta, pp, pi0
  real, dimension(m1,m2,m3) :: clsrc! liquid/ice tendency from
                                    ! cumulus parameterization
   if(level < 2  .and. mcphys_type < 2 ) return
   
   if(level == 2 .and. mcphys_type < 2) then 
       call mcphysics2(m1,m2,m3,ia,iz,ja,jz,dtlt &
	             ,clsrc                  &
                     ,tend%rct(1)            &
		     ,tend%rtt(1)            )
   
   elseif(level == 3 .or. (mcphys_type==2 .or. mcphys_type==3))  then 

       call mcphysics3(m1,m2,m3,ia,iz,ja,jz,dtlt,cpi  &
            ,theta, pp, pi0 &
            ,clsrc          &! cumulus tendency            
	    ,tend%rct(1)    &! cloud water
            ,tend%rrt(1)    &! pristine
            ,tend%rtt(1)    )! total water tendency
            

    endif
return
end  subroutine cupar2mcphysics  

!------------------------------------------------------------------------
subroutine mcphysics2(m1,m2,m3,ia,iz,ja,jz,dtlt,clsrc,rct,rtt)
implicit none
integer m1,m2,m3,ia,iz,ja,jz,k,i,j
real dtlt
real, dimension(m1,m2,m3),intent(in   ) :: clsrc
real, dimension(m1,m2,m3),intent(inout) :: rct,rtt

do j = ja,jz
   do i = ia,iz
     do k = 1,m1
       rct(k,i,j)=rct(k,i,j)+clsrc(k,i,j)
       rtt(k,i,j)=rtt(k,i,j)+clsrc(k,i,j)
     enddo
   enddo
enddo
end subroutine mcphysics2

!------------------------------------------------------------------------
subroutine mcphysics3(m1,m2,m3,ia,iz,ja,jz,dtlt,cpi,theta,pp,pi0,clsrc &
                 ,rct,rrt,rtt)
implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,k,i,j
real dtlt,cpi
real, dimension(m1,m2,m3),intent(in)    :: clsrc
real, dimension(m1,m2,m3),intent(in)    :: theta,pp,pi0   
real, dimension(m1,m2,m3),intent(inout) :: rct,rrt,rtt
real ::tempk,tem1       
real, parameter :: tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf)

do j = ja,jz
   do i = ia,iz
     do k = 1,m1
            tempk = theta(k,i,j)*(pp(k,i,j)+pi0(k,i,j))*cpi ! air temp (Kelvin)
            tem1 = max(0.0, min(1.0, (tcr-tempk)*tcrf))
            !- splitting cumulus tendency into water and ice tendencies
            rct(k,i,j)=rct(k,i,j)+clsrc(k,i,j)*(1.0-tem1) ! cloud water
            rrt(k,i,j)=rrt(k,i,j)+clsrc(k,i,j)* tem1      ! pristine ice
            !- it must include also the ice/liq tendencies at rtt for
	    !- consistency, since rtt includes ice and liq mixing ratios
	    rtt(k,i,j)=rtt(k,i,j)+clsrc(k,i,j)
     enddo
   enddo
enddo

end subroutine mcphysics3
!------------------------------------------------------------------------

subroutine prepare_lsf(nnqparm,nnshcu,iwork)

!  use mem_grell   ,only: grell_g, grell_g_sh
  use mem_grell   ,only: cuforc_g,cuforc_sh_g
  use mem_tend    ,only: tend
  use mem_scratch ,only: scratch
  use mem_grid    ,only: time,ngrid,dtlt  
  use mem_cuparm  ,only: confrq ,cuparm_g_sh 
  use node_mod    ,only: mxp,myp,mzp ,ia,iz,ja,jz
  
  implicit none
  include "i8.h"
  character(len=3) :: forcing
  integer,intent(IN) :: nnqparm,nnshcu,iwork

! forcing='RAD' ! you must also coment the lines where lsfcupar is calculated at turb_diff routine
! forcing='PBL' 
! forcing='ADV' 
 forcing='ALL' 

 IF(mod(time,confrq).lt.dtlt .or. time .lt. .01) then  
    if(iwork.eq.1) then
        cuforc_sh_g(ngrid)%lsfth = 0. 
	cuforc_sh_g(ngrid)%lsfrt = 0.
        cuforc_g   (ngrid)%lsfth = 0. 
	cuforc_g   (ngrid)%lsfrt = 0.
        return
    endif

    IF(nnqparm >= 5)then! .and. nnshcu == 3) then !srf !<<<<<<<<333 !<<<<<<<< !<<<<<<<<
       !- for this specific combination:
       !  the forcing for shallow is only due to diffusion in PBL only (which is calculated in turb routines)
       !  the forcing for deep is due to radiation + 3dim advection
       if(iwork.eq.2) then

          call atob_long(int((mxp*myp*mzp),i8),tend%THT, cuforc_g(ngrid)%lsfth)! make lsfthsh = tht due radiation
          call atob_long(int((mxp*myp*mzp),i8),tend%RTT, cuforc_g(ngrid)%lsfrt)! make lsfrtsh = rtt due radiation
          return
       endif
       ! salva a tendencia total antes da adveccao nos arrais  vt3dm e vt3dn
       ! atencao: certifique -se que os arrays vt3dm e vt3dn nao sao usados(ou zerados) dentro da rotina
       ! de adveccao
       if(iwork.eq.3) then
         
          scratch%vt3dm = 0.
          scratch%vt3dn = 0.
          call atob_long(int((mxp*myp*mzp),i8), tend%THT, scratch%vt3dm)  ! make vt3dm = tht due all process before adv
          call atob_long(int((mxp*myp*mzp),i8), tend%RTT, scratch%vt3dn)  ! make vt3dn = rtt due all process before adv
          return
       endif
       if(iwork.eq.4) then
          ! salva a tendencia advectiva nos arrais vt3dm e vt3dn (ae1m1 => a = b - c)

          call ae1m1_l(int((mxp*myp*mzp),i8), scratch%vt3dm, tend%tht, scratch%vt3dm) ! salva tend advectiva thp em a(ivt3dm)
          call ae1m1_l(int((mxp*myp*mzp),i8), scratch%vt3dn, tend%rtt, scratch%vt3dn) ! salva tend advectiva rtp em a(ivt3dn)

          ! calcula o forcing para conveccao profunda = rad  + adv
          ! for deep convection    LSF =  radiation +advection (ae1p1=> A = B + C)
          call ae1p1_l(int((mxp*myp*mzp),i8), cuforc_g(ngrid)%lsfth,cuforc_g(ngrid)%lsfth, scratch%vt3dm(1))
          call ae1p1_l(int((mxp*myp*mzp),i8), cuforc_g(ngrid)%lsfrt,cuforc_g(ngrid)%lsfrt, scratch%vt3dn(1))
          return
       endif

    ELSE

     if(iwork.eq.2) then

        call atob_long(int((mxp*myp*mzp),i8),tend%THT, cuforc_sh_g(ngrid)%lsfth)! make lsfthsh = tht due radiation
        call atob_long(int((mxp*myp*mzp),i8),tend%RTT, cuforc_sh_g(ngrid)%lsfrt)! make lsfrtsh = rtt due radiation
        
        !-srf: for only rad forcing tests
        !-srf: you must also coment the lines where lsfcupar is calculated at turb_diff routine
        if    (forcing == 'RAD') then
          cuforc_sh_g(ngrid)%lsfth = cuforc_g(ngrid)%lsfth
          cuforc_sh_g(ngrid)%lsfrt = cuforc_g(ngrid)%lsfrt    
        endif
        !- end rad forcing only tests


        !-srf: for only turb forcing tests
        if    (forcing == 'PBL') then
         cuforc_sh_g(ngrid)%lsfth = 0.
	 cuforc_sh_g(ngrid)%lsfrt = 0.
         cuforc_g(ngrid)%lsfth = 0. 
	 cuforc_g(ngrid)%lsfrt = 0.
        endif
        !- end turb forcing only tests
        return
     endif
     
     !-srf: for rad forcing only tests
     !-srf: you must also coment the lines where lsfcupar is calculated at turb_diff routine
     if(forcing == 'RAD' .and. iwork >  2) then
        return
     endif
     !- end rad forcing only tests

     !-srf: for only turb forcing tests
     if(forcing == 'PBL' .and. iwork ==  3) then
          cuforc_sh_g(ngrid)%lsfth = cuforc_g(ngrid)%lsfth
          cuforc_sh_g(ngrid)%lsfrt = cuforc_g(ngrid)%lsfrt 
     endif    
     if(forcing == 'PBL' .and. iwork >=  3) then
        return
     endif
     !- end turb forcing only tests

     !-srf
     !at this point: PBL term is calculated in diffuse routines and this is the forcing for shallow cumulus


     ! salva a tendencia total antes da adveccao nos arrais  vt3dm e vt3dn
     ! atencao: certifique -se que os arrays vt3dm e vt3dn nao sao usados(ou zerados) dentro da rotina
     ! de adveccao
     if(iwork.eq.3) then
       
        scratch%vt3dm = 0.
	scratch%vt3dn = 0.
        call atob_long(int((mxp*myp*mzp),i8), tend%THT, scratch%vt3dm)  ! make vt3dm = tht due all process before adv
        call atob_long(int((mxp*myp*mzp),i8), tend%RTT, scratch%vt3dn)  ! make vt3dn = rtt due all process before adv

        return
     endif

     if(iwork.eq.4) then
        ! salva a tendencia advectiva nos arrais vt3dm e vt3dn (ae1m1 => a = b - c)

        call ae1m1_l(int((mxp*myp*mzp),i8), scratch%vt3dm, tend%tht, scratch%vt3dm) ! salva tend advectiva thp em a(ivt3dm)
        call ae1m1_l(int((mxp*myp*mzp),i8), scratch%vt3dn, tend%rtt, scratch%vt3dn) ! salva tend advectiva rtp em a(ivt3dn)

        !-srf: for only ADV forcing tests
        if(forcing == 'ADV') then
          call atob_long(int((mxp*myp*mzp),i8), scratch%vt3dm,  cuforc_g(ngrid)%lsfth)
          call atob_long(int((mxp*myp*mzp),i8), scratch%vt3dn, cuforc_g(ngrid)%lsfrt)
	  return  
	  !CALL teste(mzp,mxp,myp,ia,iz,ja,jz,tend%THT(1),scratch%vt3dm(1),3)
          !CALL teste(mzp,mxp,myp,ia,iz,ja,jz,tend%RTT(1),scratch%vt3dn(1),4)
        endif    
        !- end ADV forcing only tests

	!CALL teste(mzp,mxp,myp,ia,iz,ja,jz,tend%THT(1),scratch%vt3dm(1),3)
        !CALL teste(mzp,mxp,myp,ia,iz,ja,jz,tend%RTT(1),scratch%vt3dn(1),4)


        ! calcula o forcing para conveccao profunda = rad + pbl turb + adv
        ! for deep convection    LSF =  radiation + pbl_turb + advection (ae1p1=> A = B + C)
        call ae1p1_l(int((mxp*myp*mzp),i8), cuforc_g(ngrid)%lsfth,cuforc_sh_g(ngrid)%lsfth, scratch%vt3dm(1))
        call ae1p1_l(int((mxp*myp*mzp),i8), cuforc_g(ngrid)%lsfrt,cuforc_sh_g(ngrid)%lsfrt, scratch%vt3dn(1))
        
        !CALL teste(mzp,mxp,myp,ia,iz,ja,jz,cuforc_g(ngrid)%lsfth(1,1,1),cuforc_sh_g(ngrid)%lsfth(1,1,1),500)
	!CALL teste(mzp,mxp,myp,ia,iz,ja,jz,tend%THT(1),cuforc_g(ngrid)%lsfrt(1,1,1),6)

        return
     endif


     if(iwork.eq.5) then 

         call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfth(1,1,1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     	 call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfrt(1,1,1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))
	 
     endif
           
   ENDIF
ENDIF

end subroutine prepare_lsf

!---------------------------------------------------
subroutine RK_prepare_lsf(nnqparm,nnshcu,iwork)

  use mem_grell   ,only: cuforc_g,cuforc_sh_g
  use mem_tend    ,only: tend
  !use mem_scratch ,only: scratch
  use mem_grid    ,only: time,ngrid,dtlt  
  use mem_cuparm  ,only: confrq ,cuparm_g_sh 
  use node_mod    ,only: mxp,myp,mzp ,ia,iz,ja,jz
  use mem_radiate, only: ilwrtyp, iswrtyp, radiate_g
  use mem_grid, only:ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4,itopo
  use mem_basic, only: basic_g 
  
  implicit none
  include "i8.h"
  character(len=3) :: forcing
  integer,intent(IN) :: nnqparm,nnshcu,iwork
  !- scratchs (local arrays)
  real :: vt3da(mzp,mxp,myp)
  real :: vt3db(mzp,mxp,myp)
  real :: vt3dc(mzp,mxp,myp)
  real :: vt3dh(mzp,mxp,myp)
  real :: vt3dj(mzp,mxp,myp)
  real :: vt3dk(mzp,mxp,myp)
  real :: vt3di(mzp,mxp,myp)
  real :: vt3df(mzp,mxp,myp)
  real :: vt3dg(mzp,mxp,myp)
  real :: vt3de(mzp,mxp,myp)
  real :: vt3dd(mzp,mxp,myp)
  real :: vctr1(mzp)
  real :: vctr2(mzp)
  real :: scr1(mzp,mxp,myp)
  integer :: i,j,k
  !- parameter to define if include or not diffusion tendencies at forcing for deep convection
  logical,parameter :: forc_deep_pbl = .false.
 
 
 IF(mod(time,confrq).lt.dtlt .or. time .lt. .01) then  
        
       !- 
       !  the forcing for shallow is only due to diffusion in PBL only (which is calculated in turb routines)
       !  the forcing for deep is due to radiation + 3dim advection 
       if(iwork.eq.1) then
        
	 !----------- include radiation for theta
	 if(ilwrtyp + iswrtyp > 0) then
	      cuforc_g(ngrid)%lsfth(1:mzp,1:mxp,1:myp)= radiate_g(ngrid)%fthrd(1:mzp,1:mxp,1:myp)
	 endif

	 !----------- include advection for theta and rv (or should be rtp?)
         vt3dd=0.0
         vt3de=0.0
         vt3df=0.0
         vt3dg=0.0
         vt3dh=0.0
         vt3di=0.0
         vt3dj=0.0
         vt3dk=0.0
         vctr1=0.0
         vctr2=0.0
         do j = 1,myp
           do i = 1,mxp
             do k = 1,mzp
                vt3da(k,i,j) = basic_g(ngrid)%uc(k,i,j)*dtlt
                vt3db(k,i,j) = basic_g(ngrid)%vc(k,i,j)*dtlt
                vt3dc(k,i,j) = basic_g(ngrid)%wc(k,i,j)*dtlt
             end do
           end do
         end do     

         call fa_preptc(mzp,mxp,myp,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df,vt3dh,vt3di, &
                        vt3dj,vt3dk,vctr1,vctr2, &
	                basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v,&
	                grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv,grid_g(ngrid)%fmapt,&
			grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi, &
                        grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%dxu,&
			grid_g(ngrid)%dyv,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,hw4,dzm,dzt,zm,zt)
         !---- thc
	 scr1(1:mzp,1:mxp,1:myp) = basic_g(ngrid)%thc(1:mzp,1:mxp,1:myp) 
	 
         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,basic_g(ngrid)%thc,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
              call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim)

         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

         ! input:  thetac , lsfth,scr1, dtlt
         ! output: lsfth
         call advtndc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc,scr1,cuforc_g(ngrid)%lsfth,dtlt)

         !---- water vapor
	 scr1(1:mzp,1:mxp,1:myp) = basic_g(ngrid)%rv(1:mzp,1:mxp,1:myp) 
	 
         ! output: scr1,vt3dg
         call fa_xc(mzp,mxp,myp,ia,iz,1,myp,basic_g(ngrid)%rv,scr1,vt3da,vt3dd,vt3dg,vt3dh,vt3di)

         ! input: scalarp, scr1,vt3db,vt3de,vt3dj,vt3di
         ! output: scr1,vt3dg
         if (jdim .eq. 1)  &
              call fa_yc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,vt3db,vt3de,vt3dg,vt3dj,vt3di,jdim)

         ! input: scalarp, scr1,vt3dc,vt3df,vt3dk, vctr1,vctr2
         ! output: scr1,vt3dg
         call fa_zc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,vt3dc,vt3df,vt3dg,vt3dk,vctr1,vctr2)

         ! input: basic(ngrid)%rv, scalart,scr1, dtlt
         ! output: lsfrt = rad + adv
         call advtndc(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%rv,scr1,cuforc_g(ngrid)%lsfrt,dtlt)

	 !- here the forcings contain rad+adv for temp and adv for water vapor
         !-end of inclusion of the advection forcings
	 
         !- flag to include diffusion in pbl (only vertical) in deep convection forcing.
	 !- the pbl forcing (cuforc_sh_g(ngrid)%lsfth and %lsfrt) has been calculated
	 !- in the turbulence routines.
	 if(forc_deep_pbl) then
          ! calcula o forcing para conveccao profunda = rad + pbl turb + adv
          ! for deep convection    LSF =  radiation + pbl_turb + advection
          cuforc_g(ngrid)%lsfth(:,:,:) = cuforc_g(ngrid)%lsfth(:,:,:)+cuforc_sh_g(ngrid)%lsfth(:,:,:)
          cuforc_g(ngrid)%lsfrt(:,:,:) = cuforc_g(ngrid)%lsfrt(:,:,:)+cuforc_sh_g(ngrid)%lsfrt(:,:,:)
	 endif
	
     endif
     if(iwork.eq.2) then 

         call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfth(1,1,1), cuparm_g_sh(ngrid)%thsrc(1,1,1))
     	 call accum(int(mxp*myp*mzp,i8), cuforc_g(ngrid)%lsfrt(1,1,1), cuparm_g_sh(ngrid)%rtsrc(1,1,1))
	 
     endif
           
ENDIF

end subroutine RK_prepare_lsf

!-------------------------------------------------------------------
subroutine teste(m1,m2,m3,ia,iz,ja,jz,tht,lsfth,it)

  real, dimension(m1,m2,m3) :: tht,lsfth
  ix=0
!  return

	     print*,'----------------teste:',it
  do j=ja,jz
     do i=ia,iz
	do k=2,m1-1
	   if(j.eq.38.and.i.eq.49)then
	      !IF(tht(k,i,j).gt.1.e-14) THEN
	      !if(k.eq.2) print*,'----------------',it,j,i
	      if(it==1)print*,'TH RAD     =',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.
	      if(it==100)print*,'TH RAD+TURB=',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.
	      if(it==500)print*,'TH DP SH=',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.
	      if(it==2)print*,'RT RAD=',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.

	      if(it==30)print*,'TH ADV=',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.
	      if(it==40)print*,'RT ADV=',k,tht(k,i,j)*86400.,lsfth(k,i,j)*86400.

	      if(it==3)print*,'TH ADV=',k,lsfth(k,i,j)*86400.
	      if(it==4)print*,'RT ADV=',k,lsfth(k,i,j)*86400.
	      if(it==5)print*,'TH RAD+TURB+ADV=',k,lsfth(k,i,j)*86400.
	      if(it==6)print*,'RT RAD+TURB+ADV=',k,lsfth(k,i,j)*86400.
	   endif
	enddo
     enddo
  enddo
  !IF(ix.eq.10) stop
end subroutine teste
!-------------------------------------------------------------------
subroutine get_traning_grell(ngr, n2, n3, train)

  use mem_grid, only: &
       platn, plonn, xmn, ymn ! INTENT(IN)

  use node_mod, only: &
       mynum,         & ! INTENT(IN)
       nodei0,        & ! INTENT(IN)
       nodej0           ! INTENT(IN)

!!$  !tmp-srf para diferente difusao numerica >
!!$  use extras, only: extra2d, NA_EXTRA2D

  implicit none
  ! Arguments:
  integer, intent(IN) :: n2, n3, ngr
  real, intent(OUT)   :: train(n2,n3)
  ! Local Variables:
  real    :: rlat(n2,n3), rlon(n2,n3), rm
  integer :: i, j, np,ii,jj

!!$  print *, "DEBUG-ALF:get_train"

  !srf-define diferentes AKMINs para melhorar estabilidade 
  !srf-sobre os Andes
!!$  !testa numero de extras:
!!$  if (NA_EXTRA2D<5) call fatal_error('NA_EXTRA2d deve ser no minimo 5')
  !default 
!!$  extra2d(5,ngr)%d2(:,:) = 1.
  !train = 1.17 !valor do La Plata.
  !----

  !print*,n2,n3,platn(ngr),plonn(ngr),xmn(1,ngr),ymn(1,ngr)
  !-calculate lat, lon of each grid box T-points
  do j=1,n3
     do i=1,n2
        call xy_ll(rlat(i,j), rlon(i,j), platn(ngr), plonn(ngr), &
             xmn(i+nodei0(mynum,ngr),ngr), ymn(j+nodej0(mynum,ngr),ngr))
     enddo
  enddo

  do j=1,n3
     do i=1,n2

        !  train(i,j) = 1.17 !valor do La Plata.
        train(i,j) = 1.3 !valor do La Plata.! 30102011



        !regiao norte
        if (rlat(i,j) > -14.) then
           if (rlon(i,j)<-45.5) then 
              !             train(i,j) = 0.43! 30102011
              train(i,j) = 0.3
           endif
        endif

        !regiao NORDESTE
        if (rlat(i,j) > -18.5) then
           if (rlon(i,j) .ge. -45.5) then 
              !              train(i,j) = 0.25! 30102011
              train(i,j) = 0.3
           endif
        endif
        !regiao CENTRO OESTE
        if (rlat(i,j) > -25. .and. RLAT(I,J) < -12.) then
           if (rlon(i,j) < -45.5) then 
              train(i,j) = 0.75
           endif
        endif

        !regiao SULDESTE
        if (rlat(i,j) > -26.5 .and. RLAT(I,J) < -14.) then
           if (rlon(i,j) .ge. 45.5) then 
              train(i,j) = 0.32
           endif
        endif

        !REGIAO SUL
        if (rlat(i,j) > -34. .and. RLAT(I,J) .le. -26.5) then
           if (rlon(i,j) .ge. 45.5) then 
              !              train(i,j) = 0.93 ! 30102011
              train(i,j) = 1.3  
           endif
        endif

     enddo
  enddo

  !media movel ( 2 viz)
  do j=3,n3-2
     do i=3,n2-2
        rm=0.
        np=0
        do jj=j-1,j+1
           do ii=i-1,i+1
              np=np+1
              rm = rm + train(ii,jj)
           enddo
        enddo
        train(I,J)=rm/np
     enddo
  enddo
end subroutine get_traning_grell
