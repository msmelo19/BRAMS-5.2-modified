!--------------------------------------------------------------------------------!
! Cumulus Parameterization by G. Grell - versions 3d and GD-FIM                  !
! Implemented in BRAMS by S. Freitas @ Feb/2012                                  !
! Rafael Mello: included parallelization for spread/g3d_smoothh arrays               !
!--------------------------------------------------------------------------------!

MODULE CUPARM_GRELL3

  use ModNamelistFile, only: namelistFile
   
  use ModMessageSet, only: &
      PostRecvSendMsgs, &
      WaitRecvMsgs
       
  use mem_basic         , only: basic_g
  use mem_tend          , only: tend
  use mem_cuparm        , only: confrq,cuparm_g,cuparm_g_sh&
                               ,nnqparm
  
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
           		   	deltayn, &   ! INTENT(IN)
				npatch,  &   ! INTENT(IN)
				   ztn,  &   ! INTENT(IN)
				   zmn,  &   ! INTENT(IN)
				akminvar     ! INTENT(IN)
				

  use rconstants        , only: tkmin 
! use extras            , only: extra3d,extra2d,na_EXTRA3D ,na_EXTRA2D
  use mem_turb          , only: turb_g,akmin
  use mem_micro         , only: micro_g
! use mem_scratch       , only: scratch
!srf
  use io_params         , only: frqanl
  use mem_leaf          , only: leaf_g, isfcl
  use micphys           , only: level,mcphys_type

!- use modules for Grell Parameterization
  use mem_grell_param   , only: mgmxp,mgmyp,mgmzp,maxiens,ngrids_cp,&
       maxens_g3d ,                        & !INTENT(IN)
       maxens2_g3d,                        & !INTENT(IN)
       maxens3_g3d,                        & !INTENT(IN)
       ensdim_g3d ,                        & !INTENT(IN)
       maxens ,                        & !INTENT(IN)
       maxens2,                        & !INTENT(IN)
       maxens3,                        & !INTENT(IN)
       ensdim ,                        & !INTENT(IN)
       icoic                 

  use mem_scratch1_grell, only: ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d,   &
                                kstabi4d,kstabm4d,xmb4d,edt4d,pwav4d,               &
                                zup5d, zdn5d,iruncon, pcup5d, prup5d,prdn5d, &
				clwup5d,tup5d,enup5d,endn5d,deup5d,dedn5d,zcup5d,                                 &
                                up_massdetr5d, up_massentr5d,                       &
                                dd_massdetr5d, dd_massentr5d


  use mem_grell         , only: cuforc_g,cuforc_sh_g

  use mem_carma, only: carma

  use mem_radiate, only: ISWRTYP, ILWRTYP,radiate_g ! Intent(in)

!-----------Grell G3d - GD-FIM - GF
  use module_cu_g3    , only:  G3DRV
  use module_cu_gd_fim, only:  GRELLDRV_FIM
  use module_cu_gf    , only:  GFDRV
  use module_cu_gf2   , only:  GFDRV2

  USE Phys_const, only: cp, p00, tcrit, g, cpor , XL, rm,rgas
  
!----------- 

  use ccatt_start, only: ccatt 
  !newCode begin
  !use mem_jules, only: jules_g
  !newCode end
  
  implicit none

  TYPE g3d_ens_vars   
     REAL, POINTER, DIMENSION(:,:)  ::apr
     REAL, POINTER, DIMENSION(:,:)  ::accapr
     REAL, POINTER, DIMENSION(:,:)  ::weight
     !-----------
  END TYPE g3d_ens_vars
  TYPE (g3d_ens_vars)    , allocatable :: g3d_ens_g(:,:),g3d_ensm_g(:,:)
 
  TYPE g3d_vars   
     REAL, POINTER, DIMENSION(:,:  )  ::xmb_deep 
     REAL, POINTER, DIMENSION(:,:  )  ::err_deep 
     REAL, POINTER, DIMENSION(:,:  )  ::xmb_shallow
     REAL, POINTER, DIMENSION(:,:,:)  ::cugd_ttens
     REAL, POINTER, DIMENSION(:,:,:)  ::cugd_qvtens
     REAL, POINTER, DIMENSION(:,:,:)  ::thsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::rtsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::clsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::usrc
     REAL, POINTER, DIMENSION(:,:,:)  ::vsrc
     REAL, POINTER, DIMENSION(:,:,:)  ::mup
  END TYPE g3d_vars
  
  TYPE (g3d_vars)       , allocatable :: g3d_g(:),g3dm_g(:)
  
  integer ::    ids,ide, jds,jde, kds,kde            & 
               ,ims,ime, jms,jme, kms,kme            & 
               ,ips,ipe, jps,jpe, kps,kpe            & 
               ,its,ite, jts,jte, kts,kte  

  
  integer,parameter ::  imomentum=0 
  integer           ::  ishallow_g3 
  
  !- define if the training will be used or not 
  INTEGER,PARAMETER:: training=0
  character(len=255) :: g3d_training_file
    
  !- define if the lateral subsidence spread will be done or not 
  INTEGER :: g3d_spread ! 1=ON, 0=OFF
  INTEGER :: cugd_avedx
  
  !- define if the horizontal smoothing is to be done or not 
  INTEGER :: g3d_smoothh! 1=ON, 0=OFF
  
  !- define if the vertical smoothing is to be done or not
  INTEGER :: g3d_smoothv! 1=ON, 0=OFF

  !- number of members of prec ensemble 
  INTEGER,PARAMETER :: train_dim= 5

  CHARACTER(LEN=6),PARAMETER,DIMENSION(train_dim) :: pre_name=(/ &
      'apr_gr' & !
     ,'apr_w ' & !
     ,'apr_mc' & !
     ,'apr_st' & !
     ,'apr_as' & !
   /)
  
  INTEGER,PARAMETER :: apr_gr=001
  INTEGER,PARAMETER :: apr_w =002
  INTEGER,PARAMETER :: apr_mc=003
  INTEGER,PARAMETER :: apr_st=004
  INTEGER,PARAMETER :: apr_as=005
 

  integer,parameter :: CPTIME = 0. !orig: CPTIME = 7200.

  integer,parameter :: i_forcing = 1

  integer,parameter :: autoconv = 1!2 ! =1, Kessler
                                      ! =2, Berry 
  integer,parameter :: aerovap = 1!3  ! =1, orig
                                    ! =2, mix orig+new
				    ! =3, new 
  !- direct link cupar-microphysics
  LOGICAL,parameter :: do_cupar_mcphys_coupling = .true.
Contains
!-----------------------------------------
  subroutine nullify_grell3(g3d_ens,g3d,ndim_train)

    implicit none
    integer, intent(in) ::ndim_train
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer i
 
    do i=1,ndim_train
      if (associated(g3d_ens(i)%apr))    nullify (g3d_ens(i)%apr)
      if (associated(g3d_ens(i)%accapr)) nullify (g3d_ens(i)%accapr)
      if (associated(g3d_ens(i)%weight)) nullify (g3d_ens(i)%weight)
    enddo
    
    if (associated(g3d%xmb_deep))    nullify (g3d%xmb_deep)
    if (associated(g3d%err_deep))    nullify (g3d%err_deep)
    if (associated(g3d%xmb_shallow)) nullify (g3d%xmb_shallow)

    if (associated(g3d%cugd_ttens))  nullify (g3d%cugd_ttens)
    if (associated(g3d%cugd_qvtens)) nullify (g3d%cugd_qvtens)

    !print *,'LFR-DBG: nullify: ',associated(g3d%thsrc);call flush(6)
    if (associated(g3d%thsrc)) nullify (g3d%thsrc)
    if (associated(g3d%rtsrc)) nullify (g3d%rtsrc)
    if (associated(g3d%clsrc)) nullify (g3d%clsrc)
    if (associated(g3d%usrc )) nullify (g3d%usrc)
    if (associated(g3d%vsrc )) nullify (g3d%vsrc)
    if (associated(g3d%mup  )) nullify (g3d%mup)
 
  end subroutine nullify_grell3
!-----------------------------------------
   subroutine alloc_grell3(g3d_ens,g3d, m1, m2, m3, ng,ndim_train)
    implicit none
    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens
    type (g3d_vars) :: g3d
    integer, intent(in) :: m1, m2, m3, ng,ndim_train
    integer :: i

    
    do i=1,ndim_train
       allocate(g3d_ens(i)%apr   (m2,m3))   ; g3d_ens(i)%apr    =0.0
       allocate(g3d_ens(i)%accapr(m2,m3))   ; g3d_ens(i)%accapr =0.0
       allocate(g3d_ens(i)%weight(m2,m3))   ; g3d_ens(i)%weight =0.0
    enddo


    allocate (g3d%xmb_deep   (m2,m3))       ;g3d%xmb_deep   =0.0
    allocate (g3d%err_deep  (m2,m3))        ;g3d%err_deep   =0.0
    allocate (g3d%xmb_shallow(m2,m3))       ;g3d%xmb_shallow=0.0
    
    allocate (g3d%cugd_ttens (m1, m2, m3))  ;g3d%cugd_ttens =0.0
    allocate (g3d%cugd_qvtens(m1, m2, m3))  ;g3d%cugd_qvtens=0.0
    !print *,'LFR-DBG: Allocating ',m1,m2,m3
    allocate (g3d%thsrc(m1, m2, m3))  ;g3d%thsrc=0.0
    allocate (g3d%rtsrc(m1, m2, m3))  ;g3d%rtsrc=0.0
    allocate (g3d%clsrc(m1, m2, m3))  ;g3d%clsrc=0.0
    !print *,'LFR-DBG: Allocating ',associated(g3d%thsrc);call flush(6)
    if(imomentum==1 .and. nnqparm(ng) >= 4) then
     allocate (g3d%usrc(m1, m2, m3))  ;g3d%usrc=0.0
     allocate (g3d%vsrc(m1, m2, m3))  ;g3d%vsrc=0.0
    endif
    allocate (g3d%mup  (m1, m2, m3))  ;g3d%mup=0.0

  end subroutine alloc_grell3
  
!-----------------------------------------
  subroutine filltab_grell3(g3d_ens,g3d,g3d_ensm,g3dm,imean, m1, m2, m3, ng,ndim_train)
    use var_tables
    implicit none
    include "i8.h"

    type (g3d_ens_vars),dimension(ndim_train) :: g3d_ens,g3d_ensm
    type (g3d_vars) :: g3d,g3dm
    integer, intent(in) :: imean, m1, m2, m3, ng,ndim_train
    integer(kind=i8) :: npts 
    integer :: i
    character (len=4) :: arrprop
    ! Fill pointers to arrays into variable tables

    npts=m2*m3
    do i=1,ndim_train
     if (associated(g3d_ens(i)%apr))  &
         call InsertVTab (g3d_ens(i)%apr   ,g3d_ensm(i)%apr    &
         ,ng, npts, imean,  &
         trim(pre_name(i))//' :2:hist:mpti:mpt3')
     
     if (associated(g3d_ens(i)%accapr))  &
         call InsertVTab (g3d_ens(i)%accapr   ,g3d_ensm(i)%accapr    &
         ,ng, npts, imean,  &
         'acc'//trim(pre_name(i)(2:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3')
    enddo

    do i=1,ndim_train
     if (associated(g3d_ens(i)%weight))  &
         call InsertVTab (g3d_ens(i)%weight   ,g3d_ensm(i)%weight    &
         ,ng, npts, imean,  &
         'weight'//trim(pre_name(i)(4:len_trim(pre_name(i))))//' :2:hist:anal:mpti:mpt3')
     
    enddo
    
    if (associated(g3d%xmb_deep))  &
         call InsertVTab (g3d%xmb_deep   ,g3dm%xmb_deep    &
         ,ng, npts, imean,  &
         'UPMF :2:hist:anal:mpti:mpt3')
    if (associated(g3d%err_deep))  &
         call InsertVTab (g3d%err_deep   ,g3dm%err_deep    &
         ,ng, npts, imean,  &
         'XIERR :2:hist:anal:mpti:mpt3')
    if (associated(g3d%xmb_shallow))  &
         call InsertVTab (g3d%xmb_shallow   ,g3dm%xmb_shallow    &
         ,ng, npts, imean,  &
         'UPMFSH3 :2:hist:anal:mpti:mpt3')


    !- 3D Arrays
    npts=m1*m2*m3
     
    !- define if the arrays will exchange 1 row x 1 line (not in use anymore)
    arrprop=''
    
    if (associated(g3d%cugd_ttens))  &
         call InsertVTab (g3d%cugd_ttens     ,g3dm%cugd_ttens      &
         ,ng, npts, imean,  &
         'TTENS :3:hist:anal:mpti:mpt3'//trim(arrprop))

    if (associated(g3d%cugd_qvtens))  &
         call InsertVTab (g3d%cugd_qvtens    ,g3dm%cugd_qvtens     &
         ,ng, npts, imean,  &
         'QVTTENS :3:hist:anal:mpti:mpt3'//trim(arrprop))
 
    if (associated(g3d%thsrc))  &
         call InsertVTab (g3d%thsrc     ,g3dm%thsrc      &
         ,ng, npts, imean,  &
         'THSRC :3:hist:anal:mpti:mpt3'//trim(arrprop))

    if (associated(g3d%rtsrc))  &
         call InsertVTab (g3d%rtsrc     ,g3dm%rtsrc     &
         ,ng, npts, imean,  &
         'RTSRC :3:hist:anal:mpti:mpt3'//trim(arrprop))

    !- this array does not need to be parallelized (only column)
    if (associated(g3d%clsrc))  &
         call InsertVTab (g3d%clsrc     ,g3dm%clsrc     &
         ,ng, npts, imean,  &
         'CLSRC :3:hist:anal:mpti:mpt3')

    
    if(imomentum==1 .and. nnqparm(ng) >= 4) then
    !- these arrays does not need to be parallelized (only column)
      if (associated(g3d%usrc))  &
         call InsertVTab (g3d%usrc     ,g3dm%usrc     &
         ,ng, npts, imean,  &
         'USRC :3:hist:anal:mpti:mpt3')
      if (associated(g3d%vsrc))  &
         call InsertVTab (g3d%vsrc     ,g3dm%vsrc     &
         ,ng, npts, imean,  &
         'VSRC :3:hist:anal:mpti:mpt3')

    endif
    if (associated(g3d%mup))  &
         call InsertVTab (g3d%mup    ,g3dm%mup     &
         ,ng, npts, imean,  &
         'MUP :3:hist:anal:mpti:mpt3')

     
  end subroutine filltab_grell3
!-------------------------------------------------------------
 
subroutine CUPARM_GRELL3_CATT(OneGrid, iens,iinqparm,iinshcu)

    USE mem_radiate, ONLY: &
         ilwrtyp, iswrtyp        ! INTENT(IN)
    use ModGrid, only: &
         Grid
!ML -- In case you want to output massflux
    use mem_stilt         , only: imassflx
!

  implicit none
  
  include "i8.h"
  integer, intent(IN) :: iens,iinqparm,iinshcu
  type(Grid), pointer :: OneGrid ! intent(in)
  integer :: i
  real :: grid_length

  REAL, DIMENSION( mxp , myp ) :: aot500,temp2m 

 if(initial.eq.2.and.time.lt.cptime-dtlt) return
 if(mod(time,confrq).lt.dtlt.or.time.lt. .01 .or.abs(time-cptime) .lt. 0.01) then
        !-start convective transport of tracers
	iruncon=1
        g3d_g(ngrid)%thsrc = 0.0
        g3d_g(ngrid)%rtsrc = 0.0
        g3d_g(ngrid)%clsrc = 0.0
        g3d_g(ngrid)%cugd_ttens  = 0.0
        g3d_g(ngrid)%cugd_qvtens = 0.0
        g3d_g(ngrid)%mup         = 0.0
        cuparm_g(ngrid)%conprr   = 0.0
        if(imomentum == 1 .and. nnqparm(ngrid) >= 4) then
	  g3d_g(ngrid)%usrc = 0.0
          g3d_g(ngrid)%vsrc = 0.0
        endif
	ishallow_g3=0
	!srf - use the old way to define the cumulus forcing

	if(i_forcing /= 1) then
                !call check(mzp * myp * mzp,tend%THT ,cuforc_g(ngrid)%lsfth ,tend%RTT  ,cuforc_g(ngrid)%lsfrt)
        	call atob(mxp * myp * mzp,tend%THT  ,cuforc_g(ngrid)%lsfth     )
        	call atob(mxp * myp * mzp,tend%RTT  ,cuforc_g(ngrid)%lsfrt     )
        endif	
!srf ----- tmp
!       cuforc_g(ngrid)%lsfth=3.*cuforc_g(ngrid)%lsfth
!       cuforc_g(ngrid)%lsfrt=3.*cuforc_g(ngrid)%lsfrt
!srf ------tmp


	!- converting WRF setting to BRAMS
        ids=1   ;ide=mxp ;jds=1   ;jde=myp ;kds=1; kde=mzp	      
	ims=1   ;ime=mxp ;jms=1   ;jme=myp ;kms=1; kme=mzp			
	ips=ia+1;ipe=iz-2;jps=ja+1;jpe=jz-2;kps=1; kpe=mzp			  
	its=ia  ;ite=iz  ;jts=ja  ;jte=jz  ;kts=1; kte=mzp-1  

        grid_length=sqrt(deltaxn(ngrid)*deltayn(ngrid))

!print*  ,ids,ide, jds,jde, kds,kde                        & 
!        ,ims,ime, jms,jme, kms,kme			   & 
!        ,ips,ipe, jps,jpe, kps,kpe			   & 
!        ,its,ite, jts,jte, kts,kte			   & 
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------

 if(iinqparm==3) then  ! G3d scheme 
   !
   !- lateral spreading
   if(g3d_spread == 0 )cugd_avedx=1
   if(g3d_spread == 1 )cugd_avedx=3

   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   end if

   CALL G3DRV( mynum,i0,j0,time             &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
	      ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
	      ,basic_g(ngrid)%wp            & !
	      ,basic_g(ngrid)%rv            & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
	      ,XL			    & !
	      ,CP			    & !
	      ,G			    & !
	      ,rm                           &
	      ,p00                          &
	      ,cpor                         & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
!              
	      ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
!
              ,training &
!===
!             ,APR_GR(ims:ime,jms:jme)	       & !
!	      ,APR_W (ims:ime,jms:jme)	       & !
!	      ,APR_MC(ims:ime,jms:jme)	       & !
!	      ,APR_ST(ims:ime,jms:jme)	       & !
!	      ,APR_AS(ims:ime,jms:jme)         & !
!             ,APR_CAPMA(ims:ime,jms:jme)      & !
!	      ,APR_CAPME(ims:ime,jms:jme)      & !
!	      ,APR_CAPMI (ims:ime,jms:jme)     & !
!             ,xmb_deep(ims:ime,jms:jme)       & !
!             ,xmb_shallow (ims:ime,jms:jme)   & !
!	      ,XF_ENS			       & !
!	      ,PR_ENS			       & !
!===
	      ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
	      ,npatch                          &
              ,radiate_g(ngrid)%rshort         &
!===
!             ,edt_out                          		&
!             ,GDC						&
!	      ,GDC2						&
!	      ,kpbl						&
!	      ,k22_shallow					&
!	      ,kbcon_shallow      				& 
!             ,ktop_shallow					& 
!===
	      ,cugd_avedx					& 
	      ,imomentum          				& 
!srf          ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ensdim,maxiens,maxens,maxens2,maxens3,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
!==
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
!	      
              ,g3d_g(ngrid)%cugd_ttens    &
              ,g3d_g(ngrid)%cugd_qvtens   &
! forcings -  for deep/shallow
	      ,cuforc_g(ngrid)%	lsfth     & ! forcing for theta deep
	      ,cuforc_g(ngrid)%	lsfrt     & ! forcing for rv deep 
	      ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
	      ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow 
              ,level                      &
	      ,micro_g(ngrid)%rcp         & ! liquid water
    	      ,micro_g(ngrid)%rrp         & ! pristine
    	      ,micro_g(ngrid)%rpp         &
	      ,micro_g(ngrid)%rsp         &
    	      ,micro_g(ngrid)%rap         &
	      ,micro_g(ngrid)%rgp         &
    	      ,micro_g(ngrid)%rhp         &
	      ,aot500                     & ! aot at 500nm
                                     )
!
!- exchange border information for parallel run
   if( g3d_spread == 1 .or. g3d_smoothh == 1) then 
      call PostRecvSendMsgs(OneGrid%SendG3D, OneGrid%RecvG3D)
      call WaitRecvMsgs    (OneGrid%SendG3D, OneGrid%RecvG3D)
   endif
!
!
!- call routine to do the lateral spread, smooths and limiters/fixers 
   CALL conv_grell_spread3d_brams(mzp,mxp,myp,ia,iz,ja,jz,dtlt,level,cugd_avedx& 
	      ,XL			    &
	      ,CP			    &
	      ,G			    &
	      ,rm                           &
	      ,p00                          &
	      ,cpor                         & 
	      ,cuparm_g(ngrid)%CONPRR       &!preci rate
              ,basic_g(ngrid)%theta         &
              ,basic_g(ngrid)%thp           &
              ,basic_g(ngrid)%pp            &
              ,basic_g(ngrid)%pi0           &
	      ,basic_g(ngrid)%rv            &
              ,tend%PT                      &
	      ,micro_g(ngrid)%rcp           & ! liquid water
    	      ,micro_g(ngrid)%rrp           & ! pristine
    	      ,micro_g(ngrid)%rpp           &
	      ,micro_g(ngrid)%rsp           &
    	      ,micro_g(ngrid)%rap           &
	      ,micro_g(ngrid)%rgp           &
    	      ,micro_g(ngrid)%rhp           &
!	      ,
              ,g3d_g(ngrid)%THSRC           & ! temp tendency
              ,g3d_g(ngrid)%RTSRC           & ! rv tendency
              ,g3d_g(ngrid)%CLSRC           & ! cloud/ice tendency
              ,g3d_g(ngrid)%cugd_ttens      &
              ,g3d_g(ngrid)%cugd_qvtens     &
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  )


!-------------------------------------------------------------
 elseif(iinqparm==7) then  ! IFS version

   CALL cupar_ifs(mynum,i0,j0,time,mzp-1,mxp,myp  &                   
              ,ia,iz,ja,jz                  &
              ,dtlt         		    & !
              ,grid_length                  & !
	      ,confrq                       & !
              ,basic_g(ngrid)%dn0  (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%dn0  (1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%up   (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%up   (1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%vp   (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%vp   (1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%theta(2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%theta(1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%thp  (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%thp  (1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%pp   (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%pp   (1,1:mxp,1:myp) & !& !
              ,basic_g(ngrid)%pi0  (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%pi0  (1,1:mxp,1:myp) & !& !
	      ,basic_g(ngrid)%wp   (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%wp   (1,1:mxp,1:myp) & !& !
	      ,basic_g(ngrid)%rv   (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%rv   (1,1:mxp,1:myp) & !& !
	      ,basic_g(ngrid)%rtp  (2:mzp,1:mxp,1:myp)  ,basic_g(ngrid)%rtp  (1,1:mxp,1:myp) & !& !

!2-d
	      ,cuparm_g(ngrid)%conprr(1:mxp,1:myp)      & !
              ,grid_g(ngrid)%rtgt    (1:mxp,1:myp)      & !
	      ,zmn(2:mzp,ngrid)                         & !
	      ,ztn(2:mzp,ngrid)                         & !
              ,g3d_g(ngrid)%xmb_deep   (1:mxp,1:myp)    &
              ,g3d_g(ngrid)%xmb_shallow(1:mxp,1:myp)    &
	      ,grid_g(ngrid)%topt      (1:mxp,1:myp)    &
              ,leaf_g(ngrid)%patch_area                 &
	      ,npatch                                   &
              ,g3d_g(ngrid)%THSRC  (2:mzp,1:mxp,1:myp)  & ! temp tendency
              ,g3d_g(ngrid)%RTSRC  (2:mzp,1:mxp,1:myp)  & ! rv tendency
              ,g3d_g(ngrid)%CLSRC  (2:mzp,1:mxp,1:myp)  & ! cloud/ice tendency
              ,g3d_g(ngrid)%USRC   (2:mzp,1:mxp,1:myp)  & ! U tendency
              ,g3d_g(ngrid)%VSRC   (2:mzp,1:mxp,1:myp)  & ! V tendency
              ,g3d_g(ngrid)%MUP    (2:mzp,1:mxp,1:myp)  & ! updraft mass flux
! forcings -  for deep/shallow
	      ,cuforc_g   (ngrid)%lsfth(2:mzp,1:mxp,1:myp)  & ! forcing for theta deep
	      ,cuforc_g   (ngrid)%lsfrt(2:mzp,1:mxp,1:myp)  & ! forcing for rv deep 
	      ,cuforc_sh_g(ngrid)%lsfth(2:mzp,1:mxp,1:myp)  & ! forcing for theta shallow
	      ,cuforc_sh_g(ngrid)%lsfrt(2:mzp,1:mxp,1:myp)  & ! forcing for rv shallow 
!
  	      ,turb_g(ngrid)%sflux_r(1:mxp,1:myp)           &
              ,turb_g(ngrid)%sflux_t(1:mxp,1:myp)           &
              ,turb_g(ngrid)%tkep(2:mzp,1:mxp,1:myp) ,turb_g(ngrid)%tkep(1,1:mxp,1:myp)  &
              ,tend%PT                                      & !
              ,g3d_ens_g(apr_gr,ngrid)%apr(1:mxp,1:myp)  &
              ,g3d_ens_g(apr_w ,ngrid)%apr(1:mxp,1:myp)  &
              ,g3d_ens_g(apr_mc,ngrid)%apr(1:mxp,1:myp)  &
              ,g3d_ens_g(apr_st,ngrid)%apr(1:mxp,1:myp) &
              ,g3d_ens_g(apr_as,ngrid)%apr(1:mxp,1:myp) &
              )
!
!-------------------------------------------------------------
 elseif(iinqparm==4) then  ! GD FIM version
   !-no lateral spreading for this scheme
   cugd_avedx=1
   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   endif

   if(isfcl == 5) then 
     !newCode begin
     !temp2m(:,:) = jules_g(ngrid)%t2mj(:,:)
     !newCode end
   else 
     temp2m(:,:) =0.5*( basic_g(ngrid)%theta(1,:,:)* &
                       (basic_g(ngrid)%pp(1,:,:)+basic_g(ngrid)%pi0(1,:,:))/cp + &
                        basic_g(ngrid)%theta(2,:,:)*&
		       (basic_g(ngrid)%pp(2,:,:)+basic_g(ngrid)%pi0(2,:,:))/cp )
   endif

   if(iinshcu == 3) ishallow_g3=1
   
   CALL GRELLDRV_FIM(                        &
               mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
	      ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
	      ,basic_g(ngrid)%wp            & !
	      ,basic_g(ngrid)%rv            & !
	      ,basic_g(ngrid)%rtp           & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
!	      ,XL			    & !
!	      ,CP			    & !
!	      ,G			    & !
!	      ,rm                           &
	      ,p00                          &
	      ,cpor                         & !
	      ,rgas                         & !
	      ,zmn(:,ngrid)                 & !
	      ,ztn(:,ngrid)                 & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%err_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
	      ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
              ,training &
	      ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
	      ,npatch                          &
              ,radiate_g(ngrid)%rshort         &
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
              ,g3d_g(ngrid)%USRC       & ! U tendency
              ,g3d_g(ngrid)%VSRC       & ! V tendency
              ,g3d_g(ngrid)%MUP        & ! updraft mass flux
	      ,cuforc_g(ngrid)%	lsfth     & ! forcing for theta deep
	      ,cuforc_g(ngrid)%	lsfrt     & ! forcing for rv deep 
	      ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
	      ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow 
              ,level                      &
	      ,micro_g(ngrid)%rcp         & ! liquid water
	      ,aot500                     &! aot at 500nm
	      ,temp2m                     &! aot at 500nm
  	      ,turb_g(ngrid)%sflux_r      &
              ,turb_g(ngrid)%sflux_t      &
              ,turb_g(ngrid)%tkep         &
              ,TKMIN                      &
	      ,akmin(ngrid)               &
!- for convective transport-start
              ,ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	             & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
	      ,kstabi4d		     & 
	      ,kstabm4d		     & 
	      ,xmb4d		     & 
	      ,edt4d		     & 
	      ,pwav4d		     & 
	      ,pcup5d  		     & 
              ,up_massentr5d	     &        
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d   		     & 
!- for convective transport- end
          			     )
!-------------------------------------------------------------
 elseif(iinqparm==5) then  ! GF 2014 scheme 
   !
   !- no lateral spreading
   cugd_avedx=1
 
   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   endif

   if(iinshcu == 3) ishallow_g3=1
    
   CALL GFDRV( CCATT                        &                    
              ,mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
	      ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
	      ,basic_g(ngrid)%wp            & !
	      ,basic_g(ngrid)%rv            & !
	      ,basic_g(ngrid)%rtp           & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
	      ,XL			    & !
	      ,CP			    & !
	      ,G			    & !
	      ,rm                           &
	      ,p00                          &
	      ,cpor                         & !
	      ,rgas                         & !
	      ,zmn(:,ngrid)                 & !
	      ,ztn(:,ngrid)                 & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%err_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
	      ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
              ,training &
	      ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
	      ,npatch                          &
              ,radiate_g(ngrid)%rshort         &
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
              ,g3d_g(ngrid)%cugd_ttens    &
              ,g3d_g(ngrid)%cugd_qvtens   &
	      ,cuforc_g(ngrid)%	lsfth     & ! forcing for theta deep
	      ,cuforc_g(ngrid)%	lsfrt     & ! forcing for rv deep 
	      ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
	      ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow 
              ,level                      &
	      ,micro_g(ngrid)%rcp         & ! liquid water
    	      ,micro_g(ngrid)%rrp         & ! pristine
    	      ,micro_g(ngrid)%rpp         &
	      ,micro_g(ngrid)%rsp         &
    	      ,micro_g(ngrid)%rap         &
	      ,micro_g(ngrid)%rgp         &
    	      ,micro_g(ngrid)%rhp         &
	      ,aot500                     &! aot at 500nm
  	      ,turb_g(ngrid)%sflux_r      &
              ,turb_g(ngrid)%sflux_t      &
              ,turb_g(ngrid)%tkep         &
              ,TKMIN                      &
	      ,akmin(ngrid)               &
!- for convective transport-start
              ,ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	             & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
	      ,kstabi4d		     & 
	      ,kstabm4d		     & 
	      ,xmb4d		     & 
	      ,edt4d		     & 
	      ,pwav4d		     & 
	      ,pcup5d  		     & 
              ,up_massentr5d	     &        
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d   		     & 
!- for convective transport- end
          			     )
!-------------------------------------------------------------
 elseif(iinqparm==6) then  ! GF 2015 scheme 
   !
   !- no lateral spreading
   cugd_avedx=1
 
   if(ilwrtyp==4 .or. iswrtyp==4) THEN
   	aot500(:,:)=carma(ngrid)%aot(:,:,11)
   else
   	aot500(:,:)=0.0
   endif

   if(isfcl == 5) then 
     !newCode begin
     !temp2m(:,:) = jules_g(ngrid)%t2mj(:,:)
     !newCode end
   else 
     temp2m(:,:) =0.5*( basic_g(ngrid)%theta(1,:,:)* &
                       (basic_g(ngrid)%pp(1,:,:)+basic_g(ngrid)%pi0(1,:,:))/cp + &
                        basic_g(ngrid)%theta(2,:,:)*&
		       (basic_g(ngrid)%pp(2,:,:)+basic_g(ngrid)%pi0(2,:,:))/cp )
   endif

   if(iinshcu == 3) ishallow_g3=1
    
   CALL GFDRV2(mgmxp,mgmyp,mgmzp,ngrid,ngrids_cp,iens &
              ,mynum,i0,j0,time,mzp,mxp,myp &
              ,dtlt         		    & !
              ,grid_length                  & !
              ,autoconv                     & !
              ,aerovap                      & !
              ,basic_g(ngrid)%dn0           & !
	      ,cuparm_g(ngrid)%CONPRR       & !
              ,basic_g(ngrid)%up            & !
              ,basic_g(ngrid)%vp            & !
              ,basic_g(ngrid)%theta         & !
              ,basic_g(ngrid)%thp           & !
              ,basic_g(ngrid)%pp            & !
              ,basic_g(ngrid)%pi0           & !
	      ,basic_g(ngrid)%wp            & !
	      ,basic_g(ngrid)%rv            & !
	      ,basic_g(ngrid)%rtp           & !
              ,grid_g(ngrid)%RTGT           & !
              ,tend%PT                      & !
	      ,XL			    & !
	      ,CP			    & !
	      ,G			    & !
	      ,rm                           &
	      ,p00                          &
	      ,cpor                         & !
	      ,rgas                         & !
	      ,zmn(:,ngrid)                 & !
	      ,ztn(:,ngrid)                 & !
              ,g3d_ens_g(apr_gr,ngrid)%apr  &
              ,g3d_ens_g(apr_w ,ngrid)%apr  &
              ,g3d_ens_g(apr_mc,ngrid)%apr  &
              ,g3d_ens_g(apr_st,ngrid)%apr  &
              ,g3d_ens_g(apr_as,ngrid)%apr  &
              ,g3d_g(ngrid)%xmb_deep        &
              ,g3d_g(ngrid)%err_deep        &
              ,g3d_g(ngrid)%xmb_shallow     &
	      ,g3d_ens_g(apr_gr,ngrid)%weight &
              ,g3d_ens_g(apr_w ,ngrid)%weight &
              ,g3d_ens_g(apr_mc,ngrid)%weight &
              ,g3d_ens_g(apr_st,ngrid)%weight &
              ,g3d_ens_g(apr_as,ngrid)%weight &
              ,training &
	      ,grid_g(ngrid)%topt              &
              ,leaf_g(ngrid)%patch_area        &
	      ,npatch                          &
              ,radiate_g(ngrid)%rshort         &
	      ,cugd_avedx					& 
	      ,imomentum          				& 
              ,ensdim_g3d,maxiens,maxens_g3d,maxens2_g3d,maxens3_g3d,icoic      & 
              ,ishallow_g3                                      & 
	      ,ids,ide, jds,jde, kds,kde                        & 
              ,ims,ime, jms,jme, kms,kme                        & 
              ,ips,ipe, jps,jpe, kps,kpe                        & 
              ,its,ite, jts,jte, kts,kte                        & 
              ,g3d_g(ngrid)%THSRC      & ! temp tendency
              ,g3d_g(ngrid)%RTSRC      & ! rv tendency
              ,g3d_g(ngrid)%CLSRC      & ! cloud/ice tendency
              ,g3d_g(ngrid)%USRC       & ! U tendency
              ,g3d_g(ngrid)%VSRC       & ! V tendency
              ,g3d_g(ngrid)%MUP        & ! updraft mass flux
	      ,cuforc_g(ngrid)%	lsfth     & ! forcing for theta deep
	      ,cuforc_g(ngrid)%	lsfrt     & ! forcing for rv deep 
	      ,cuforc_sh_g(ngrid)%lsfth   & ! forcing for theta shallow
	      ,cuforc_sh_g(ngrid)%lsfrt   & ! forcing for rv shallow 
              ,level                      &
	      ,micro_g(ngrid)%rcp         & ! liquid water
	      ,aot500                     &! aot at 500nm
	      ,temp2m                     &! aot at 500nm
  	      ,turb_g(ngrid)%sflux_r      &
              ,turb_g(ngrid)%sflux_t      &
              ,turb_g(ngrid)%tkep         &
              ,TKMIN                      &
	      ,akmin(ngrid)               &
	      ,do_cupar_mcphys_coupling   &
!- for convective transport-start
              ,ierr4d  		     &  	     
	      ,jmin4d  		     & 
	      ,kdet4d  		     & 
	      ,k224d	             & 
	      ,kbcon4d 		     & 
	      ,ktop4d  		     & 
	      ,kpbl4d  		     & 
	      ,kstabi4d		     & 
	      ,kstabm4d		     & 
	      ,xmb4d		     & 
	      ,edt4d		     & 
	      ,pwav4d		     & 
	      ,pcup5d  		     & 
              ,up_massentr5d	     &        
	      ,up_massdetr5d	     &
	      ,dd_massentr5d	     &
	      ,dd_massdetr5d	     &
	      ,zup5d		     &
	      ,zdn5d   		     & 
	      ,prup5d  		     & 
	      ,prdn5d  		     & 
	      ,clwup5d 		     & 
	      ,tup5d   		     & 
!- for convective transport- end
          			     )
 endif 

!--- filling the output tendencies for the level k =1 and k=mzp
  g3d_g(ngrid)%THSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(2    ,1:mxp,1:myp)  
  g3d_g(ngrid)%RTSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(2    ,1:mxp,1:myp)  
  g3d_g(ngrid)%CLSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(2    ,1:mxp,1:myp)  
  g3d_g(ngrid)%THSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%THSRC(mzp-1,1:mxp,1:myp)  
  g3d_g(ngrid)%RTSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%RTSRC(mzp-1,1:mxp,1:myp)  
  g3d_g(ngrid)%CLSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%CLSRC(mzp-1,1:mxp,1:myp)  
  if(imomentum==1 .and. nnqparm(ngrid) >= 4) then
   g3d_g(ngrid)%USRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%USRC (2    ,1:mxp,1:myp)  
   g3d_g(ngrid)%VSRC(1  ,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (2    ,1:mxp,1:myp) 
   g3d_g(ngrid)%USRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%USRC (mzp-1,1:mxp,1:myp)  
   g3d_g(ngrid)%VSRC(mzp,1:mxp,1:myp)= g3d_g(ngrid)%VSRC (mzp-1,1:mxp,1:myp) 
  endif 
!-------------------------------------------------------------
 endif! 002
!-------------------------------------------------------------
! stores precipitation rate for each closure, only for output/training
 if (training > 0) then
   do i=1,train_dim
     call update(mxp*myp, g3d_ens_g(i,ngrid)%accapr,g3d_ens_g(i,ngrid)%apr,dtlt)
   enddo
 else
   ! temporary use to save cloud work functions for debug purposes
   do i=1,train_dim
    g3d_ens_g(i,ngrid)%accapr(:,:)=g3d_ens_g(i,ngrid)%apr(:,:)
   enddo
 endif
!----------------------------------------------------------

 
 call accum(int(mxp*myp*mzp,i8), tend%tht, g3d_g(ngrid)%thsrc)
 call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%rtsrc)

 if(imomentum == 1 .and. nnqparm(ngrid) >= 4) then 
  call accum(int(mxp*myp*mzp,i8), tend%ut, g3d_g(ngrid)%usrc)
  call accum(int(mxp*myp*mzp,i8), tend%vt, g3d_g(ngrid)%vsrc)
 endif

 call update(mxp*myp, cuparm_g(ngrid)%aconpr   ,cuparm_g(ngrid)%conprr   ,dtlt)

 if(do_cupar_mcphys_coupling) then
   call cupar2mcphysics(mzp,mxp,myp,ia,iz,ja,jz,ngrid,dtlt,& 
                        g3d_g  (ngrid)%clsrc   ,&
			basic_g(ngrid)%theta   ,& 
			basic_g(ngrid)%pp      ,&  
			basic_g(ngrid)%pi0      )
 else
  !if is not to have direct coupling include cloud/ice source at rtotal tendency
  call accum(int(mxp*myp*mzp,i8), tend%rtt, g3d_g(ngrid)%clsrc)
 			
 endif
!
!--------- Convective Transport based on mass flux scheme -
  if (CCATT == 1 .and. iruncon == 1 .and. (iinqparm==5 .or. iinqparm==6) ) then
     !print*,"convective transport turned off"
     !return
     !print*,"convective transport turned on"
     if(iinqparm==5 .and. iens .ne. 1 ) &
       stop 'conv transp with GF scheme version 2014 only for deep convection'
     
     !- this call convective transport for deep convection
     call trans_conv_mflx_GF(1,iinqparm)
     
     !- if shallow convection was solved by GF version 2015, call again
     !- the convective transport routine to include the transport
     !- by shallow convection scheme. 
     if(iinqparm==6 .and. iinshcu == 3)  call trans_conv_mflx_GF(2,iinqparm)
    
  endif 

! [ML------------- Stilt - BRAMS coupling  ------------------
  if (imassflx == 1 ) then 
       
       !-srf -  mass fluxes from deep convection
       if( iinqparm==5 .or. iinqparm==6 )                                 &
          call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
              ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp		  &    
              ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d 	  &
              ,kstabi4d,kstabm4d,xmb4d,edt4d				  &
              ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d	  & 
              ,1)! = iens 
       !-srf if shallow convection was solved by GF version 2015, call again
       !-    the convective transport routine to include the mass fluxes
       !-    from the shallow convection scheme. 
       if(  iinqparm==6 .and. iinshcu == 3 )                              &
          call prep_convflx_to_stilt(mzp,mxp,myp,ia,iz,ja,jz              &
              ,mgmxp,mgmyp,mgmzp,maxiens,ngrid,ngrids_cp		  &    
              ,ierr4d,jmin4d,kdet4d,k224d,kbcon4d,ktop4d,kpbl4d 	  &
              ,kstabi4d,kstabm4d,xmb4d,edt4d				  &
              ,zcup5d,pcup5d,enup5d,endn5d,deup5d,dedn5d,zup5d,zdn5d	  & 
              ,2)! = iens 
      
   endif
! ------------- Stilt - BRAMS coupling  ------------------ ML]
!
end subroutine CUPARM_GRELL3_CATT
!
!-------------------------------------------------------------------------------------------------
!
subroutine init_weights(ng,n2,n3,nnqparm)
implicit none
integer, intent(in)::ng,n2,n3,nnqparm
integer :: it,i,j
real sumx,hweight
!- ordem dos pesos
!apr_gr=001
!apr_w =002
!apr_mc=003
!apr_st=004
!apr_as=005

if(training == 0) return

!-- training on closures
if(training == 1) then
   if(nnqparm==3) hweight = 0.2
   if(nnqparm==5) hweight = 0.25
   do j=1,n3
     do i=1,n2
      do it=1,train_dim
    
       g3d_ens_g(it,ng)%weight(i,j)=hweight
       !print*,'weights=', it,i,j, g3d_ens_g(it,ng)%weight(i,j)
       
       !if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175 
       !if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.25
       !if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
       !if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
       if(nnqparm==5) then
       
        !-special treatment over the ocean       
        if(leaf_g(ng)%patch_area(i,j,1) .gt. 0.9 ) then ! water
           if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
	   if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.425
           if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.1667
           if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.1667
           if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.1667

        else ! land
       
           if(it==apr_as) g3d_ens_g(it,ng)%weight(i,j)=0.0
	   if(it==apr_st) g3d_ens_g(it,ng)%weight(i,j)=0.175
           if(it==apr_gr) g3d_ens_g(it,ng)%weight(i,j)=0.25
           if(it==apr_w ) g3d_ens_g(it,ng)%weight(i,j)=0.25
           if(it==apr_mc) g3d_ens_g(it,ng)%weight(i,j)=0.25
       endif
      endif 
      !g3d_ens_g(it,ng)%weight(i,j)=float(i+j)*exp(-(float(it-2))**2)*float(i*j)

    enddo;enddo;enddo

!-- training on CAPS
elseif(training == 2) then

    do j=1,n3; do i=1,n2
    
       g3d_ens_g(apr_gr,ng)%weight(i,j)=0.3333 
       g3d_ens_g(apr_w ,ng)%weight(i,j)=0.3333 
       g3d_ens_g(apr_mc,ng)%weight(i,j)=0.3333 
       g3d_ens_g(apr_st,ng)%weight(i,j)=0.0 
       g3d_ens_g(apr_as,ng)%weight(i,j)=0.0 

     enddo;enddo

endif



return! <<<<
if(training == 1) then
 !- normalize a 1 
 do j=1,n3
    do i=1,n2
     sumx=0.
     do it=1,train_dim
       sumx=sumx+g3d_ens_g(it,ng)%weight(i,j)
     enddo
      do it=1,train_dim
      g3d_ens_g(it,ng)%weight(i,j) = g3d_ens_g(it,ng)%weight(i,j)/sumx
     enddo
 enddo;enddo    
endif

end subroutine init_weights
!-------------------------------------------------------------
!-------------------------------------------------------------
SUBROUTINE conv_grell_spread3d_brams(m1,m2,m3,ia,iz,ja,jz,dt, &
               level,cugd_avedx,                              & 
	       XLV,CP,G,r_v,p00,cpor,                         & 					 
               conprr, theta,thetail,pp,pi0,                  &
	       rv,pt,rcp,rrp,rpp,rsp,rap,rgp,rhp,             &
	       RTHcuten,				      &
               RQVcuten,				      &
	       RQCcuten,				      &
               cugd_ttens,				      & 
	       cugd_qvtens,				      & 
               apr_gr,					      & 
	       apr_w,					      & 
	       apr_mc,					      & 
	       apr_st,					      & 
	       apr_as					      ) 
	       
IMPLICIT NONE

   INTEGER,      INTENT(IN   ) :: m1,m2,m3,ia,iz,ja,jz,level,cugd_avedx
   REAL,         INTENT(IN   ) :: dt
   REAL,         INTENT(IN   ) :: XLV, R_v
   REAL,         INTENT(IN   ) :: CP,G, cpor, p00
   
   REAL, DIMENSION(m1,m2,m3),INTENT(IN   ) ::     &
 	  	 theta   ,& 
          	 thetail ,& 
          	 pp	 ,& 
          	 pi0	 ,& 
        	 pt	 ,&
        	 rv      ,rcp,rrp,rpp,rsp,rap,rgp,rhp           

   
   REAL, DIMENSION(m2,m3),INTENT(INOUT) ::   &
               conprr,                       &
               apr_gr,			     & 
	       apr_w ,			     & 
	       apr_mc,			     & 
	       apr_st,			     & 
	       apr_as			      
   
   
   
   REAL, DIMENSION(m1,m2,m3),INTENT(INOUT) ::    &
                        RTHcuten,    &
                        RQVcuten,    &
                        RQCcuten,    &
                        cugd_ttens,  & 
	                cugd_qvtens			  
   

  ! local var 
   REAL  ::   exner,r_sol,r_liq,fxc,tempk,dfxcdt,outt
   INTEGER :: j,i,k,kk,jfs,jfe,ifs,ife,kts,kte,ii,jj
   INTEGER :: cugd_spread

   REAL, DIMENSION (m1,m2,m3) :: & ! orig (its-2:ite+2,kts:kte,jts-2:jte+2) ::     &
          RTHcuten_tmp, &  ! tmp RTHcuten
	  RQVcuten_tmp     ! tmp RQVcuten

   REAL, DIMENSION (m2,m3) :: & ! orig (its-2:ite+2,jts-2:jte+2) ::
          Qmem

   REAL   :: & ! orig (its-1:ite+1,jts-1:jte+1) :: 
          smTT,smTQ

   REAL, DIMENSION (m1) :: & ! orig (kts:kte) :: 
          conv_TRASHT,conv_TRASHQ

   REAL :: Qmem1,Qmem2,Qmemf,Thresh

  !-initial settings
  ! g3d_smoothh=0  ! 0 or 1: do horizontal smoothing
  ! g3d_smoothv=0  ! 0 or 1: do vertical smoothing
   cugd_spread=cugd_avedx/2 ! = 0/1 => no/do spreading

   RTHcuten_tmp  = 0.0
   RQVcuten_tmp  = 0.0
   Qmem       = 1.0
   smTT       = 0.0
   smTQ       = 0.0
   conv_TRASHT= 0.0
   conv_TRASHQ= 0.0
   jfs=ja
   jfe=jz
   ifs=ia
   ife=iz
   kts=2
   kte=m1-1 !check if this correct or should be kte=m1
   
   !if(g3d_smoothh ==1 .or. cugd_spread > 0) then
   !  jfs=1
   !  jfe=m3
   !  ifs=1
   !  ife=m2
   !endif
   
   !- store input tendencies
   ! *** jm note -- for smoothing this goes out one row/column beyond tile in i and j
   do j=1,m3
     do i=1,m2
         RTHcuten_tmp(:,i,j)=RTHcuten (:,i,j) 
         RQVcuten_tmp(:,i,j)=RQVcuten (:,i,j) 
     enddo
   enddo



! ---------------- spreading   section --------------
   do j=ja,jz
     do i=ia,iz
!
! for high res run, spread the subsidence
! this is tricky......only consider grid points where there was no rain,
! so cugd_tten and such are zero!
!
!      if do spreading
       if(cugd_spread > 0)then
         do k=kts,kte
	    do jj=j-1,j+1 ! only 3 neighboors
	      do ii=i-1,i+1 ! only 3 neighboors

               RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)     &
                                            +Qmem(ii,jj)*cugd_ttens(k,ii,jj)

               RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)     &
                                            +Qmem(ii,jj)*cugd_qvtens(k,ii,jj)
             enddo
           enddo
         enddo
!      end spreading
!
!      if not spreading
       elseif(cugd_spread == 0)then
         do k=kts,kte
           RTHcuten_tmp(k,i,j)=RTHcuten_tmp(k,i,j)+cugd_ttens (k,i,j)
           RQVcuten_tmp(k,i,j)=RQVcuten_tmp(k,i,j)+cugd_qvtens(k,i,j)
         enddo
       endif
!
     enddo  ! end i
   enddo  ! end j

! ----------------horizontal smoothing  section --------------

!- if not doing horizontal smoothing, get the final tendencies
  if(g3d_smoothh == 0)then 
      do j=ja,jz
         do i=ia,iz
            do k=kts,kte
               RTHcuten(k,i,j)=RTHcuten_tmp(k,i,j) 
               RQVcuten(k,i,j)=RQVcuten_tmp(k,i,j)
            enddo ! enf k
          enddo  ! end j
      enddo  ! end j
	  
!- if doing horizontal smoothing ...      
   else if(g3d_smoothh == 1)then   
      do k=kts,kte      
        do j=ja,jz
           do i=ia,iz

	    smTT = 0.0
            smTQ = 0.0
	    do jj=j-1,j+1 ! only 3 neighboors
	      do ii=i-1,i+1 ! only 3 neighboors
               smTT = smTT +RTHcuten_tmp(k,ii,jj)
	       smTQ = smTQ +RQVcuten_tmp(k,ii,jj)
            
              enddo  ! end ii
            enddo  ! end jj
	    
            RTHcuten(k,i,j)=(3.*RTHcuten_tmp(k,i,j) + smTT)/12.
            RQVcuten(k,i,j)=(3.*RQVcuten_tmp(k,i,j) + smTQ)/12.

            enddo  ! end i
          enddo  ! end j
       enddo  ! end k

   endif  ! g3d_smoothh
  !
  ! - checking and limiting moistening/heating rates
  !
   do j=ja,jz
      do i=ia,iz
        !--- moistening section ------
	Qmemf  = 1.0
        Thresh = 1.e-20
        do k=kts,kte              
         
	 if(RQVcuten(k,i,j) < 0.0) then
	    Qmem1 = rv(k,i,j)+RQVcuten(k,i,j)*dt
	    if(Qmem1 < Thresh)then
              Qmem1 = RQVcuten(k,i,j)
              Qmem2 = (Thresh-rv(k,i,j))/dt
              Qmemf = min(Qmemf,Qmem2/Qmem1)
              Qmemf = max(0.,Qmemf)
              Qmemf = min(1.,Qmemf)
            endif
          endif

         enddo  ! end k
         ! - limiting moistening 
         do k=kts,kte
          RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
          RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
          RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
	  cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf			
          cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf			  

         enddo ! end k
         
	 !- limiting precip for consistency
	 conprr (i,j) = conprr (i,j)*Qmemf
         apr_gr (i,j) = apr_gr (i,j)*Qmemf		      
	 apr_w  (i,j) = apr_w  (i,j)*Qmemf		      
	 apr_mc (i,j) = apr_mc (i,j)*Qmemf		      
	 apr_st (i,j) = apr_st (i,j)*Qmemf		      
	 apr_as (i,j) = apr_as (i,j)*Qmemf		      
         ! no futuro inclua tambem o limting para o fluxo de massa
	 ! xmb(i,j)=xmb(i,j)*qmemf
	 
	 !--- heating section ------
         Thresh=200. ! max heating/cooling rate allowed  K/day
!srf         Thresh=100. ! max heating/cooling rate allowed  K/day
         Qmemf=1.
         Qmem1=0.
         
	 do k=kts,kte
            Qmem1=abs(RTHcuten(k,i,j))*86400. 

            if(Qmem1 > Thresh)then
              Qmem2 = Thresh/Qmem1
              Qmemf = min(Qmemf,Qmem2)
              Qmemf = max(0.,Qmemf) 
            endif
         enddo
	 
         ! - limiting heating/cooling 
         do k=kts,kte
          RTHcuten   (k,i,j) = RTHcuten   (k,i,j)*Qmemf
          RQVcuten   (k,i,j) = RQVcuten   (k,i,j)*Qmemf
          RQCcuten   (k,i,j) = RQCcuten   (k,i,j)*Qmemf
	  cugd_ttens (k,i,j) = cugd_ttens (k,i,j)*Qmemf			
          cugd_qvtens(k,i,j) = cugd_qvtens(k,i,j)*Qmemf			  
         enddo ! end k
	 
	 !- limiting precip for consistency
	 conprr (i,j) = conprr (i,j)*Qmemf
         apr_gr (i,j) = apr_gr (i,j)*Qmemf		      
	 apr_w  (i,j) = apr_w  (i,j)*Qmemf		      
	 apr_mc (i,j) = apr_mc (i,j)*Qmemf		      
	 apr_st (i,j) = apr_st (i,j)*Qmemf		      
	 apr_as (i,j) = apr_as (i,j)*Qmemf		      
         ! no futuro inclua tambem o limting para o fluxo de massa
	 ! xmb(i,j)=xmb(i,j)*qmemf

     enddo  ! end i
   enddo  ! end j
  !
  ! ---  vertical smooth ------------
  ! 
   if (g3d_smoothv == 1)then
   
    do j=ja,jz
      do i=ia,iz

          do k=kts+2,kte-2
            conv_TRASHT(k)= .25*(RTHcuten(k-1,i,j)+2.*RTHcuten(k,i,j)+RTHcuten(k+1,i,j))
            conv_TRASHQ(k)= .25*(RQVcuten(k-1,i,j)+2.*RQVcuten(k,i,j)+RQVcuten(k+1,i,j))
          enddo
          do k=kts+2,kte-2
            RTHcuten(k,i,j)=conv_TRASHT(k)
            RQVcuten(k,i,j)=conv_TRASHQ(k)
          enddo
     enddo  ! end i
    enddo  ! end j

   endif
   
  ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
  ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
    
   !if(level <=2) then
   
    do j=ja,jz; do i=ia,iz; do k=kts,kte
    
	!if(RTHCUTEN (k,i,j) /= 0.0) then
	
           ! Converte tend da temperatura (OUTT) em tend de theta (OUTTEM)
           ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
           ! Exner's function = pp(k,i,j)+pi0(k,i,j)
            exner	   = pp(k,i,j) + pi0(k,i,j)
           ! tendencia do theta devida a conv profunda
           RTHcuten(k,i,J) = cp/exner * RTHCUTEN(k,i,J) !- theta(k,i,j)*pt(k,i,j)/exner
        !endif
	!RQVcuten(k,i,J) = RQVCUTEN(k,i,J)+ RQCCUTEN(k,i,J)
        !RQCcuten(k,i,J) = 0.
	
    enddo; enddo; enddo
    
   
   !elseif(level > 2) then 
   !
   ! do j=ja,jz; do i=ia,iz; do k=kts,kte
   !	    !
   !	 ! - tend na temperatura (para uso na conversão do thetail
   !	 outt=RTHCUTEN (k,i,j)
   !	 ! Exner's function = pp(k,i,j)+pi0(k,i,j)
   !	 exner= pp(k,i,j) + pi0(k,i,j)
   !	 if(outt /= 0.0 ) then
   !	   !
   !	   ! converte tend da temperatura (outt) em tend de theta (outtem)
   !	   ! cp*T=Pi*Theta => cp dT/dt = Theta*dPi/dt + Pi*dTheta/dt,
   !
   !	   ! tendencia do theta  devida a conv profunda
   !	   RTHCUTEN (k,i,j) = cp/exner * RTHCUTEN(k,i,j) - theta(k,i,j)*pt(k,i,j)/exner
   !
   !	 endif
   !
   !	 ! tendencia do theta_il devida a conv profunda
   !	 r_liq= max(0.,rcp(k,i,j) + rrp(k,i,j))
   !
   !	 r_sol= max(0.,rsp(k,i,j)+rpp(k,i,j)+ &
   !		       rap(k,i,j)+rgp(k,i,j)+  &
   !		       rhp(k,i,j))
   !	 
   !	 tempk = theta(k,i,j)*(exner)/cp ! air temp (Kelvin)
   !
   !	 if(tempk.le.253) then
   !	   fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.)) 
   !	   
   !	   dfxcdt = 2.83e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.))
   !	   
   !	  RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt ) 
   !     
   !     else
   !     
   !       fxc =   (2.5e6*r_liq+2.83e6*r_sol)/(cp*amax1(tempk,253.)) 
   !
   !!orig     dfxcdt = 2.5e6*OUTQC(I,K)*cuten(i)/(cp*amax1(tempk,253.)) - & 
   !!orig         fxc/(cp*amax1(tempk,253.)) * cp * OUTT(I,K)
!  !         
   !	  dfxcdt = 2.5e6*RQCCUTEN(k,i,J)/(cp*amax1(tempk,253.)) - & 
   !          	   fxc/(cp*amax1(tempk,253.)) * cp * OUTT
   !       
   !       RTHCUTEN (k,i,j) = (1./(1.+fxc))*( RTHCUTEN (k,i,j) - thetail(k,i,j)*dfxcdt ) 
   !     
   !     endif
   !
   ! enddo; enddo; enddo
   !endif
  !- tendencies at boundaries
   RTHcuten(1,ia:iz,ja:jz)=RTHcuten(2,ia:iz,ja:jz)
   RQVcuten(1,ia:iz,ja:jz)=RQVcuten(2,ia:iz,ja:jz)
   RQCcuten(1,ia:iz,ja:jz)=RQCcuten(2,ia:iz,ja:jz)
  
   RTHcuten(m1,ia:iz,ja:jz)=RTHcuten(m1-1,ia:iz,ja:jz)
   RQVcuten(m1,ia:iz,ja:jz)=RQVcuten(m1-1,ia:iz,ja:jz)
   RQCcuten(m1,ia:iz,ja:jz)=RQCcuten(m1-1,ia:iz,ja:jz)
  
END SUBROUTINE conv_grell_spread3d_brams

!-------------------------------------------------------------

  subroutine StoreNamelistFileAtCup_grell3(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    
    g3d_spread = oneNamelistFile%g3d_spread
    g3d_smoothh = oneNamelistFile%g3d_smoothh
    g3d_smoothv = oneNamelistFile%g3d_smoothv
    
  end subroutine StoreNamelistFileAtCup_grell3

END MODULE CUPARM_GRELL3

!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
!************************************************************************
! P. Bechtold              ECMWF/Reading
!                          last update 10/2006
! --------------------------------------------------------------------
! #### This program reads the GATE dataset (161 soundings) and
! #### calls the convection routine.
! 
! #### All thermodynamic and dynamic variables are supposed to be
! #### situated at the same vertical levels labeled from KLEV (first level
! #### above bottom) to 1 (top)
! #### Fluxes are supposed to be on half levels and have dimensions KLEV+1
! #### therefore, also half level pressure and geopotential height are given
! 
subroutine cupar_ifs( mynum,i0,j0,time,klev,m2,m3 &
                     ,ia,iz,ja,jz                 &
 		     ,dt			  &
		     ,dx			  &
		     ,confrq                      & !
!
		     ,rho	 ,rho_s 	  &
		     ,u 	 ,u_s		  & 
		     ,v 	 ,v_s		  & 
		     ,theta	 ,theta_s	  & 
		     ,thetail	 ,thetail_s	  & 
		     ,pp	 ,pp_s  	  & 
		     ,pi0	 ,pi0_s 	  & 
		     ,w 	 ,w_s		  & 
		     ,rv	 ,rv_s  	  & 
		     ,rtp	 ,rtp_s 	  & 
!-d2
!
		     ,pratec			  &
		     ,rtgt			  & 
		     ,zm			  & 
		     ,zt			  & 
		     ,xmb_deep 	                  &  
		     ,xmb_shallow		  & 
		     ,ht			  &
		     ,patch_area		  &
		     ,npat			  &
		     ,rthcuten  		  &
		     ,rqvcuten  		  &
		     ,rqccuten  		  &
		     ,rucuten  		          &
		     ,rvcuten  		          &
		     ,mup  		          &
		    ! forcings - for deep/shallow
		     ,rthften			  &
		     ,rqvften			  &
		     ,rthblten  		  &
		     ,rqvblten  		  &
		    ! 
		     ,sflux_r	          & 
		     ,sflux_t	          &
		     ,tke     ,tke_s	  &
		     ,pitend   &! )
                     ,apr_gr,apr_w,apr_mc,apr_st,apr_as)



!USE PARKIND1  ,ONLY : JPIM     ,JPRB
!
IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: mynum,i0,j0,klev,m2,m3,npat,ia,iz,ja,jz
  REAL   , INTENT(IN) :: dt, dx, time,confrq
  REAL   , INTENT(IN), DIMENSION(klev)  ::  zm,zt

  REAL   , INTENT(IN), DIMENSION(klev,m2,m3) ::U,    &
                                             V,    &
                                             W,    &
                                             rv,   &
                                             rtp,  &
	           			     theta   ,& 
                   			     thetail ,& 
                   			     pp      ,& 
                   			     pi0     ,& 
                                             rho     ,&
                                             tke     ,&   
                                  RTHFTEN,  RQVFTEN  ,&
                                  RTHBLTEN,RQVBLTEN

  REAL   , INTENT(IN), DIMENSION(klev+1,m2,m3) ::pitend !don´t change this

  REAL   , INTENT(IN), DIMENSION(m2,m3) ::   U_s,    &
                                             V_s,    &
                                             W_s,    &
                                             rv_s,   &
                                             rtp_s,  &
	           			     theta_s   ,& 
                   			     thetail_s ,& 
                   			     pp_s      ,& 
                   			     pi0_s     ,& 
                                             rho_s     ,&
                                             tke_s      

   REAL, DIMENSION(m2,m3,npat), INTENT(IN) :: patch_area
   REAL, DIMENSION(m2,m3     ), INTENT(IN) :: sflux_r,sflux_t,rtgt,ht 
   

   REAL, DIMENSION(m2,m3     ), INTENT(OUT) :: pratec,xmb_deep,xmb_shallow
   REAL, DIMENSION(klev,m2,m3), INTENT(OUT) :: RTHCUTEN,    &
                                               RQVCUTEN,    &
                                               RQCCUTEN,    & 
                                               RUCUTEN ,    & 
                                               RVCUTEN ,    &
					       MUP
 
   !REAL, DIMENSION(m2,m3 )            :: XLAND
   REAL, DIMENSION( m2,m3  ),  INTENT(OUT) ::                  &
                          apr_gr,apr_w,apr_mc,apr_st,apr_as
   INTEGER, PARAMETER :: JPRB = 8

!--- local arrays
! --------------------------------------------------------------------
!*Set Dimensions for convection call 
!for GATE soundings
!INTEGER, PARAMETER :: KLON = 161 ! number of soundings for GATE
!INTEGER, PARAMETER :: KLEV = 41  ! number of vertical levels
!for BRAMS coupling
INTEGER, PARAMETER :: KLON = 1 ! number of soundings for GATE
!
!
!
INTEGER, PARAMETER :: KTRAC= 2   ! number of chemical tracers
!               (put KTRAC to zero if you don't want/need tracers in your
!                host model)
 
! --------------------------------------------------------------------
!*Set horizontal + vertical loop bounds
!
INTEGER:: KIDIA=1        ! start index for  horizontal computations
INTEGER:: KFDIA=KLON     ! end index for    horizontal computations
INTEGER:: KTDIA=1        ! end index for vertical computations
                         ! defined as KLEV + 1 - KTDIA
! --------------------------------------------------------------------
!*Define convective in and output arrays
!
REAL(KIND=JPRB), DIMENSION(KLON,KLEV):: PGEO,PPRES,PT,PQ,PU,PV,PVERVEL, &
                               PRC,PRI,                              &
                               PTTEN, PQTEN, PRCTEN, PRITEN,        &
                               PUMF, PDMF, PMFUDE_RATE,PMFDDE_RATE,  &
                               PUTEN, PVTEN, PTU, PURV, PURCI,       &
                               zrvten,ztten,zq1,zq2,zqr,zadvt,zadvq
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,KTRAC):: PC, PCTEN
REAL(KIND=JPRB), DIMENSION(KTRAC):: PSCAV
REAL(KIND=JPRB), DIMENSION(KLON)    :: PCAPE, PCIN, PURAIN
INTEGER, DIMENSION(KLON)  :: KCLTOP, KCLBAS
!
REAL(KIND=JPRB), DIMENSION(KLON,KLEV+1):: PPRESH,PGEOH,PSHFLX,PLHFLX,PPRLFLX, PPRSFLX
INTEGER, DIMENSION(KLON)               :: KBOTSC, KTYPE
LOGICAL, DIMENSION(KLON)               :: LDLAND,LDCUM,LDSC
REAL(KIND=JPRB), DIMENSION(KLON)       :: PSSTRU, PSSTRV
REAL(KIND=JPRB), DIMENSION(KLEV)       :: PREF
! --------------------------------------------------------------------
!*Set switches for convection call 

INTEGER:: NSMAX     ! horizontal resolution/sepctral truncation
                         ! Nota: the convective adjustment time scale is 1 h for NSMAX <=319
                         !       and 1200 s for NSMAX > 319, so you can change this time scale either
                         !       by changing NSMAX or by changing it in sucumf.F90
!
REAL(KIND=JPRB)   :: PDTCONV  ! time step for call of convection scheme

! All other switches are defined in routine SUCUMF that you may modify
! -------------------------------------------------------------------------------------
!*Physical meaning of convective arrays

!INPUT:
                            ! LDLAND: Land Sea mask (TRUE if Land)
                            ! PGEO  : Geopotential on full levels (g*m)
                            ! PGEOH : Geopotential on half levels (g*m)
                            ! PPRES : full-level pressure (Pa)
                            ! PPRESH: half-level pressure (Pa)
                            ! PREF  : reference pressure of model levels or average pressure over domain 
                            ! PT    : temperature (K)
                            ! PQ    : specific humidity (kg/kg)
                            ! PVERVEL: vert. velocity (Pa/s)
                            ! PU    : u(m/s)
                            ! PV    : v(m/s)
                            ! PRC, PRI : specific cloud water and ice of environment (kg/kg)
                            ! PC    : chemical tracers ()
                            ! PSHFLX : sensible heat flux (W/m^2)
                            ! PLHFLX : latent heat flux (W/m^2)
                            !      NOTA: heat fluxes are  defined positive when directed upwards
                            !            so POSITIVE during day over land
                            !               - in Call of convection code sign is reversed
                            ! PTTEN, PQTEN  : T, r_v model tendency (K/s), (1/s)
                            ! PC    : chemical tracers ()
                            ! PCTEN          : tracer tendency (1/s)
                            ! PSCAV: tracer scavenging coefficient

!OUTPUT:
                            ! KCLTOP, KCLBAS : convective cloud top and base levels (model levels)
                            ! KTYPE : Type of convection (1=deep,2=shallow,3=mid-level
                            !                             0=None)
                            ! LDCUM : Logical, TRUE if cumulus convection
    !ignore this output     ! LDSC  : Logical, TRUE if stratocumulus convection
    !ignore this output     ! KBOTSC: Base of Stratocumulus (model level)
                            ! PTTEN, PQTEN  : T, r_v updated tendency (K/s), (1/s) = model +convective
                            ! PRCTEN, PRITEN : r_c, r_i convective tendency (1/s)
                            ! PPRLFLX, PPRSFLX: liquid/solid precipitation flux
                            !                      at each level (m/s)
                            ! PURAIN: updraft precipitation flux (no evaporation in downdraft)
                            ! PUMF, PDMF: updraft and downdraft mass flux (kg/(s m^2))
                            ! PMFUDE_RATE : updraft   detrainment rate (kg/(s m^3)) - for chemical modelling
                            ! PMFDDE_RATE : downdraft detrainment rate (kg/(s m^3))

                            ! PUTEN, PVTEN   : u, v  convective tendency (m/s^2) 
                            ! PCTEN          : convective tracer tendency (1/s) 
                            ! PTU   : updraft temperature (K)
                            ! PURV  : water vapor content in updraft (kg/kg)                                 
                            ! PURCI : cloud water+ice content in updraft (kg/kg)                                 
                            ! PCAPE, PCIN : Cape (J/kg), CIN (J/kg)

!DIAGNOSTC:                 zq1,zq2 : observed apparent heat/moisture sources (K/day)
 
!IMPORTANT !!!!!!!!!!!!!!!!!!
!NOTA: In prognostic applications you have to use as Input actual model tendencies
!      for T and   q_v
!      The convection scheme will then return updated tendencies -this is necesary
!      for shallow boundary-layer equilibrium closure and for organized entrainment.
!      If you do not have these tendencies then you better use a w* shallow closure instead
!      by putting the switch LMFSCL_WSTAR=.true. in sucumf.F90 (but you will still nedd the
!      surface fluxes as input)
!----------------- --------------------------------------------------------------------
!
!

!* some constants needed for comparison with data
!* conevction scheme contains its own set of constants
REAL(KIND=JPRB), PARAMETER:: XTJOUR = 86400. ! day in s
REAL(KIND=JPRB), PARAMETER:: XG = 9.80665 ! gravity constant
REAL(KIND=JPRB), PARAMETER:: XCP = 1005.46
REAL(KIND=JPRB), PARAMETER:: XLV= 2500800. , p00 = 1.e5  ,   rgas = 287.  ,cpor = xcp / rgas 
REAL(KIND=JPRB) :: ZEPS   
REAL(KIND=JPRB) :: zlon ,zdz,exner,cpdTdt
INTEGER :: JL, JK , I,J, kr! loop variables


!#include "cumastrn.intfb.h"
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!--srf array for grads output
!vargrads(nz,nlon,nvarmax)
INTEGER, PARAMETER :: NVARMAX = 100, gradsoutput=0
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,NVARMAX):: vargrads
INTEGER :: NREC,NV,NVAR,KLEVGRADS,NVARD2D,NVAR3D
CHARACTER (LEN=50) :: runname, runlabel
CHARACTER (LEN=40),DIMENSION(NVARMAX,2) :: gradsname
LOGICAL :: XLMFPEN
LOGICAL :: XLMFSCV
LOGICAL :: XLMFMID
LOGICAL :: XLMFDD
LOGICAL :: XLMFDUDV
LOGICAL :: XLMFUVDIS
LOGICAL :: XLMFSMOOTH
LOGICAL :: XLMFTRAC
LOGICAL :: XLMFWSTAR
INTEGER, SAVE :: read_namelist=0
NAMELIST /run/  runname, runlabel                        &
               ,XLMFPEN  ,XLMFSCV   ,XLMFMID   ,XLMFDD   &
               ,XLMFDUDV ,XLMFUVDIS ,XLMFTRAC 
!-------------------------------------------------------------------------------


#ifdef ifsoff 

    stop " CUPAR IFS not included, NNQPARM must be =< 6" 

#elif ifson

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

ZEPS=XLV/XCP*XTJOUR

  !- reads information in closure.run file - only at 1st time.
  if (read_namelist == 0) then
   print*,"Opening namelist closure.run at time=",time
   OPEN(15,FILE='closure.run',STATUS='old',FORM='formatted')
   READ(15,nml=run)
   CLOSE(15)
   print*,"running ifs with the options:"
   print*,"deep      -mid     -shallow=",  XLMFPEN,  XLMFMID  , XLMFSCV   
   print*,"downdrafts-momentum-tracers=",  XLMFDD,XLMFDUDV,XLMFTRAC
   call flush(6)  
   read_namelist = 1
  endif
  !
  NSMAX=319      ! horizontal resolution/sepctral truncation
                 ! Nota: the convective adjustment time scale is 1 h for NSMAX <=319
                 !       and 1200 s for NSMAX > 319, so you can change this time scale either
                 !       by changing NSMAX or by changing it in sucumf.F90
  !
  !PDTCONV = 900. ! time step for call of convection scheme
  PDTCONV = confrq
  
  CALL SUCST(54,20020211,0,0)
  CALL SUCUMF(NSMAX,KLEV,PREF &
           ,XLMFPEN  ,XLMFSCV   ,XLMFMID   ,XLMFDD   &
           ,XLMFDUDV ,XLMFUVDIS ,XLMFTRAC )
  CALL SU_YOETHF
  CALL SUPHLI
  CALL SUVDF
  CALL SUVDFS
  CALL SUCLDP

  !-loop at BRAMS grid
  DO j=ja,jz
   DO i=ia,iz


    do jl=1,klon ! 
       !  z(m)   P(hPa)   T(C)   Q(g/kg)  w(Pa/s) '
      do JK=1,KLEV
          kr=klev-jk+1
	  
          !PGEO(JL,JK)=PGEO(JL,JK)*XG
          PGEO(JL,JK)=XG*(zt(kr)*rtgt(i,j)+ht(i,j))
	  
          !PPRES(JL,JK)=PPRES(JL,JK)*1.E2
          PPRES(JL,JK)=((pp(kr,i,j)+pi0(kr,i,j))/xcp)**cpor*p00
	  
	  !PT(JL,JK)=PT(JL,JK)+273.16
	  PT(JL,JK)=theta(kr,i,j)*(pp(kr,i,j)+pi0(kr,i,j))/xcp
          
	  !PQ(JL,JK)=PQ(JL,JK)*1.E-3
	  PQ(JL,JK)=rv(kr,i,j)
          
	  PRC(JL,JK)=0.    ! initialize condensate
          PRI(JL,JK)=0.    !

          PVERVEL(JL,JK)= -XG*rho(kr,i,j) * w(kr,i  ,j  ) 
	  PU     (JL,JK)= 0.5*( u(kr,i,j) + u(kr,i-1,j  ) )
	  PV     (JL,JK)= 0.5*( v(kr,i,j) + v(kr,i  ,j-1) ) 
	  !print*,"X",i,j,kr,jk,PU (JL,JK), PT(JL,JK),PGEO(JL,JK);call flush(6)
  	  !print*,"Y",i,j,kr,jk,zt(kr),theta(kr,i,j) , rv(kr,i,j);call flush(6)
 	  !print*,"Z",i,j,kr,jk,PU (JL,JK), PT(JL,JK),PGEO(JL,JK);call flush(6)
     enddo

      !** 1b Define half-levels: In your model take them directly
      !                          from you Sigma coordinate definition if available

          PPRESH(jl,1) =PPRES(jl,1)
          PGEOH (jl,1) =pgeo (jl,1)
      
      do jk=2,klev
          PPRESH(jl,jk)   =0.5*(PPRES(jl,jk)+PPRES(jl,jk-1))
          PGEOH (jl,jk)   =0.5*(pgeo (jl,jk)+pgeo (jl,jk-1))
      enddo
          PPRESH(jl,klev+1) = 0.5*(PPRES(jl,klev) + ((pp_s(i,j)+pi0_s(i,j))/xcp)**cpor*p00 )!surface values
          PGEOH (jl,klev+1) = XG*ht(i,j)   !pgeo (jl,klev)-10.   ! normally PGEOH(jl,klev+1) is zero by def.
 
   end do ! jl

   zlon=real(klon)
   do jk=klev,1,-1
     pref(jk)=sum(ppres(:,jk))/zlon
   enddo
   !** Initialize fluxes etc for the present case

   !PSHFLX(:,KLEV+1)= 10. ! sensible heat flux W/m2  - surface value
   !PLHFLX(:,KLEV+1)=120. ! latent   heat flux W/m2
   PSHFLX(:,KLEV+1)= sflux_t(i,j)*rho(2,i,j)*XCP ! sensible heat flux W/m2  - surface value
   PLHFLX(:,KLEV+1)= sflux_r(i,j)*rho(2,i,j)*XLV ! latent   heat flux W/m2


   ! ### for this test case we just specify an arbitrary decrease of turb. fluxes with height ###
   !     in your model use your model values
   DO JK=KLEV,1,-1
      PLHFLX(:,JK)=.9*PLHFLX(:,JK+1)
      PSHFLX(:,JK)=.9*PSHFLX(:,JK+1)
   END DO
   ! ### 

   !LDLAND=.FALSE.
   !- check this later
   LDLAND=patch_area(i,j,1).lt.0.9 !<<<
   
   
   
   
   DO JK=1,KLEV
    DO JL=1,KLON 
      
      kr=klev-jk+1
      
      PTTEN(JL,JK)=0.  ! Recall, in prognostic applications
      PQTEN(JL,JK)=0.  ! PTTEN and PQTEN should be actual model tendencies 
                       ! for other variables set them to ZERO

      !* use GATE tendencies for humidity here
      !ptten(jl,jk)=-zq1(jl,jk)/xtjour
      !PQten(jl,jk)=zq2(jl,jk)/zeps
      !ptten(jl,jk)= 1./xtjour*(zqr(jl,jk)+zadvt(jl,jk))
      !PQten(jl,jk)= ZADVQ(jl,jk)
     
      exner= pp(kr,i,j)+pi0(kr,i,j)
      !-versao implementada em 2014 (only rad+adv as the forcing)
      PTten(jl,jk)= (exner*rthften(kr,i,j) + theta(kr,i,j)*pitend(kr+1,i,j))/xcp
      PQten(jl,jk)=  rqvften(kr,i,j)
      
      !------------      
      !-versao implementada em 2015 ( rad+adv+pbl as the forcing!)
      PTten(jl,jk)= PTten(jl,jk) + RTHBLTEN(kr,i,j)*exner/xcp
      PQten(jl,jk)= PQten(jl,jk) + RQVBLTEN(kr,i,j)
 
      !-versao implementada em 2015 (update T and rv wirh rad+adv+pbl tendencies)
      PT(JL,JK)=PT(JL,JK) + PTten(jl,jk)*dt
      PQ(JL,JK)=PQ(JL,JK) + PQten(jl,jk)*dt
      !------------      
 
           
      PUTEN (jl,jk)=0.  !u tend
      PVTEN (jl,jk)=0.  !v tend
      PRCTEN(jl,jk)=0. 
      PRITEN(jl,jk)=0.
    ENDDO
   ENDDO

   !** Initialize a boundary-layer and upper-tropospheric tracer
   !   scavenging coefficients should be zero by default.

      PCTEN(:,:,:)=0.
      PC(:,:,:)=0.
      PC(:,KLEV-6:KLEV,1)=1.
      PC(:,KLEV-30:KLEV-25,2)=1.
      PSCAV(1:KTRAC)=0.

    !!** 2b Call the Convection scheme
    !----------------------------------------------------------
    !
    !
    zrvten=pqten ! save tendencies to print out later only convective part
    ztten =ptten

    CALL CUMASTRN &
  (  KIDIA,    KFDIA,	 KLON,     KTDIA,    KLEV &
  ,  2,        2,	 LDLAND,   PDTCONV &
  ,  PT,       PQ,	PU,	  PV,	    PRC+PRI &
  ,  PVERVEL, -PLHFLX/XLV,	  -PSHFLX &
  ,  PPRES,    PPRESH,   PGEO,     PGEOH &
  ,  PTTEN,    PQTEN,	 PUTEN,    PVTEN &
  ,  PRCTEN,   PRITEN &
  ,  LDCUM,    KTYPE,	 KCLBAS,   KCLTOP &
  ,  KBOTSC,   LDSC &
  ,  PTU,      PURV,	 PURCI  &
  ,  PPRLFLX,  PPRSFLX,  PURAIN &
  ,  PUMF,     PDMF  &
  ,  PMFUDE_RATE,	 PMFDDE_RATE,	 PCAPE &
  ,  KTRAC,    PC,	 PCTEN,    PSCAV,apr_gr(i,j),apr_as(i,j) )

  ! only for CAPE diagnostic - this provides better/smoother CAPE estimation 
  !                            for forecaster than the PCAPE from CUMASTRN

    CALL CUANCAPE2(KIDIA,  KFDIA,   KLON,    KLEV,&
                   PPRES,  PPRESH,  PT,      PQ,    PCAPE,  PCIN)

  !-model feedback
  
   DO JL=1,KLON 
      !err_deep(i,j)=0.
      !if(LDCUM(jl)) then; err_deep(i,j)=1.
       RTHCUTEN(:,i,j)=0.0
       RQVCUTEN(:,i,j)=0.0
       RQCCUTEN(:,i,j)=0.0
       RUCUTEN (:,i,j)=0.0
       RVCUTEN (:,i,j)=0.0
       MUP     (:,i,j)=0.0
      !--rainfall output
      !--converting from m/s to mm/s
      pratec(i,j)= 1000.*( pprlflx(jl,klev+1)+pprsflx(jl,klev+1) )
      !
      if(pratec(i,j).gt.0.)then
         !-- mass flux (deep/shallow)
         xmb_deep(i,j)= pumf(jl,kclbas(jl))
         !print*,"ifs-rain=",i,j,pratec(i,j),xmb_deep(i,j),LDCUM(jl);call flush(6)
      
         DO JK=1,KLEV
       
           kr=klev-jk+1

           !- converting Dtemp/Dt to Dtheta/ Dt
    	   RTHCUTEN(kr,i,j)=(ptten(jl,jk)-ztten(jl,jk))*xcp/(pp(kr,i,j) + pi0(kr,i,j))![K/s]
	   !water vapor
	   RQVCUTEN(kr,i,j)= pqten(jl,jk)-zrvten(jl,jk)	 ! [kg/kg/s]
	   !ice + cloud
	   RQCCUTEN(kr,i,j)= prcten(jl,jk)+priten(jl,jk)  ! [kg/kg/s]	    
	   !
	   !- u,v tendencies
	   RUCUTEN(kr,i,j) = puten(jl,jk)  ! du/dt	  [m/s/s]
           RVCUTEN(kr,i,j) = pvten(jl,jk)  ! dv/dt	  [m/s/s]
	    
	   !print*,"tend",kr,pratec(i,j),RTHCUTEN(kr,i,j)*86400.,RQVCUTEN(kr,i,j)*86400.,&
	   !        RQCCUTEN(kr,i,j)*86400.;call flush
           MUP(kr,i,j)  =pumf(jl,jk)   !Mflxups [kg/s/m^2]
         ENDDO
      endif

   ENDDO
  
  !----------------------------------------------------------
  !-srf --- add outup for GrADS vizualization
  !- array output vargrads(nz,nlon,nvar)
  IF(gradsoutput==1)then! .or. pratec(i,j)>5.e-6 ) then 
     do jl=1,klon
       !if(jl==1) then
       !   do jk=klev,1,-1
       !    write(*,133) ppres(jl,jk)*1.e-2
       !  enddo
       ! endif
       do jk=1,klev
           zdz=xg/(ppresh(jl,jk+1)-ppresh(jl,jk))
           vargrads(jl,jk,1) =ppres(jl,jk)*1.e-2    !P  (hPa)
           vargrads(jl,jk,2) =pgeo(jl,jk)/xg        ! Z (m)
!- cumulus feedback
           vargrads(jl,jk,3) =(ptten(jl,jk)-ztten(jl,jk))*xtjour ! dT/dt     [K/day]
           vargrads(jl,jk,4) =(pqten(jl,jk)-zrvten(jl,jk))*zeps  ! dqv/dt    [K/day]
!- env forcing
!           vargrads(jl,jk,3) =(ztten (jl,jk)) ! dT/dt	 [K/s]
!           vargrads(jl,jk,4) =(zrvten(jl,jk)) ! dqv/dt	[kg/kg/s]
!
           vargrads(jl,jk,5) =(prcten(jl,jk)+priten(jl,jk))*zeps ! d(ql)/dt  [K/day]
           vargrads(jl,jk,6) =pumf(jl,jk)  !Mflxup   [kg/s/m^2]
           vargrads(jl,jk,7) =pdmf(jl,jk)  !Mflxdown [kg/s/m^2]
           vargrads(jl,jk,8) =pprlflx(jl,jk)*3.6e3  !Prflx_liq  [mm/h]
           vargrads(jl,jk,9) =pprsflx(jl,jk)*3.6e3  !Prflx_ice  [mm/h]
           vargrads(jl,jk,10) =purci(jl,jk)*1.e3    !rci_up     [g/kg]
           vargrads(jl,jk,11) =puten(jl,jk)*xtjour  ! du/dt     [m/s/day]
           vargrads(jl,jk,12) =pvten(jl,jk)*xtjour  ! dv/dt     [m/s/day]
           vargrads(jl,jk,13) =pu(jl,jk)!pc(jl,jk,1)          
           vargrads(jl,jk,14) =pv(jl,jk)!pc(jl,jk,1)+pcten(jl,jk,1)*pdtconv
           vargrads(jl,jk,15) =pt(jl,jk)!pc(jl,jk,2)
           vargrads(jl,jk,16) =1.e-3*pq(jl,jk)!pc(jl,jk,2)+pcten(jl,jk,2)*pdtconv 
      end do 
    enddo
    !- surface quantities
      vargrads(1:klon,1,17) = pprlflx(1:klon,klev+1)*3600.*1.e3
      vargrads(1:klon,1,18) = pprsflx(1:klon,klev+1)*3600.*1.e3
      do jl=1,klon
       vargrads(jl,1,19) =float(ktype (jl)) 
       vargrads(jl,1,20) = pgeo(jl,kcltop(jl))/xg 
       vargrads(jl,1,21) = pgeo(jl,kclbas(jl))/xg 
      enddo
      vargrads(1:klon,1,22) = (pcape(1:klon)) 
      vargrads(1:klon,1,23) = (pcin (1:klon)) 
    !===> number of total variables :
    nvar = 23
    if(nvar>nvarmax) stop 'nvarmax must be increased' 
    PRINT*,'Writing GrADS control file:',trim(runname)//'-out.ctl'
      gradsname(1,1)  ='P      '    ;gradsname(1,2)  =  'PRESSURE  (hPa)'
      gradsname(2,1)  ='Z      '    ;gradsname(2,2)  =  'height  (m)'
      gradsname(3,1)  ='dtdt   '    ;gradsname(3,2)  =  'temp tend (K/day)'
      gradsname(4,1)  ='dqdt  '    ;gradsname(4,2)  =  'qv tend (K/day)'
      gradsname(5,1)  ='dqldt  '    ;gradsname(5,2)  =  'qc tend (K/day)'
      gradsname(6,1)  ='mup    '    ;gradsname(6,2)  =  'm flux up (kg/s/m^2)'
      gradsname(7,1)  ='mdn    '    ;gradsname(7,2)  =  'qc tend  (kg/s/m^2)'
      gradsname(8,1)  ='pprl   '    ;gradsname(8,2)  =  'prec flux lig [mm/h]'
      gradsname(9,1)  ='ppri   '    ;gradsname(9,2)  =  'prec flux ice [mm/h]'
      gradsname(10,1) ='rci_ip '    ;gradsname(10,2) =  '??? [g/kg]'
      gradsname(11,1) ='dudt   '    ;gradsname(11,2) =  'u tend (m/s/day)'
      gradsname(12,1) ='dvdt   '    ;gradsname(12,2) =  'v tend m/s/day)'
      gradsname(13,1) ='us   '    ;gradsname(13,2) =  'env U (m/s)'
      gradsname(14,1) ='vs   '    ;gradsname(14,2) =  'env V (m/s)'
      gradsname(15,1) ='t   '    ;gradsname(15,2) =  'env T (K)'
      gradsname(16,1) ='q   '    ;gradsname(16,2) =  'env Q (g/kg)'

      gradsname(17,1) ='precip'     ;gradsname(17,2) =  'liquid surf prec [mm/h]'
      gradsname(18,1) ='prsol'     ;gradsname(18,2) =  'solid surf prec [mm/h]'
      gradsname(19,1) ='type  '    ;gradsname(19,2) =  'conv-type'
      gradsname(20,1) ='cltop '    ;gradsname(20,2) =  'cloud top index'
      gradsname(21,1) ='clbas '    ;gradsname(21,2) =  'cloud base indec'
      gradsname(22,1) ='pcape  '   ;gradsname(22,2) =  'pcape'
      gradsname(23,1) ='pcin   '   ;gradsname(23,2) =  'pcin'

       OPEN(20,file=trim(runname)//'-out.ctl',status='unknown')
       write(20,2001) '^'//trim(runname)//'-out.gra'
       write(20,2002) 'undef -9.99e33'
!       write(20,2002) 'options zrev'
       write(20,2002) 'title '//trim(runlabel)
       write(20,2003) 1,0.,1. ! units m/km
       write(20,2004) klon,1.,1.
       write(20,2005) klev,(ppres(1,jk)*1.e-2,jk=klev,1,-1)
       write(20,2006) 1,'00:00Z01JAN2000','1mn'
       write(20,2007) nvar
       do nv=1,nvar
  	   klevgrads=0
  	   if(nv<17) klevgrads=klev
  	   write(20,2008) gradsname(nv,1),klevgrads,gradsname(nv,2)
  	 enddo
  	 write(20,2002) 'endvars'
  	 
    2001 format('dset ',a)
    2002 format(a)
    2003 format('xdef ',i4,' linear ',2f15.3)
    2004 format('ydef ',i4,' linear ',2f15.3)
    2005 format('zdef ',i4,' levels ',60f6.0)
    2006 format('tdef ',i4,' linear ',2a15)
    2007 format('vars ',i4)
    2008 format(a10,i4,' 99 ',a40)!'[',a8,']')
  ! 2008 format(a10,i4,' 99' )
    2055 format(60f7.0)
     133 format (1x,F7.0)
    CLOSE(20)
    
     print*, 'opening GrADS file:',trim(runname)//'-out.gra'
     OPEN(19,FILE= trim(runname)//'-out.gra',  &
     !FORM='unformatted',ACCESS='direct',STATUS='unknown', RECL=(klon)) !intel
     FORM='unformatted',ACCESS='direct',STATUS='unknown', RECL=4*(klon))!pgi
     NREC=0
     do nv=1,nvar
  	klevgrads=1
  	if(nv<17) klevgrads=klev   ! nv = 17 -> last 2-dim variable
  	do jk=klevgrads,1,-1
  	  nrec=nrec+1
  	  WRITE(19,REC=nrec) real(vargrads(:,jk,nv),4)
  	enddo	 
     enddo
     close (19)
     !- text file using GATE format
     DO jl=1,klon
      write(4,*) klev,real(PDTCONV,4),NSMAX,LDLAND, real(PSHFLX(jl,KLEV+1),4)&
                     ,real(PLHFLX(jl,KLEV+1),4)    ,real(PPRESH(jl,klev+1),4)&
		     ,real(PGEOH (jl,klev+1),4)
      DO JK=1,KLEV
 
     
        write(4,*)PGEO(JL,JK),PPRES(JL,JK),PT(JL,JK),PQ(JL,JK),        &
                  PU(JL,JK),PV(JL,JK),PVERVEL(JL,JK),  &
		  ztten (jl,jk),zrvten(jl,jk), PPRESH(jl,jk),PGEOH(jl,jk) 
		  
      END DO
    END DO
    !- end of text file using GATE format
     
     
      !-srf-end--------------------------------------------------------------------------------
      ! print out tendencies in K/day in order to compare with observations
      
      do jl=1,klon
        
        write(8,*)'%Sounding ',jl
        
        if( ktrac>=1) then
        write(8,*)' %	  P	 Z	dT/dt	 dqv/dt d(ql)/dt Mflxup  Mflxdown&
        & Prflx_liq Prflx_ice	rci_up     du/dt     dv/dt   Trac1_ini   Trac1_new'
        write(8,*)' %  [hPa]	 [m]	    ---- [K/day] ----	     [kg/(sm^2]&
        &	 [mm/h] 	[g/kg]      [m/s/day]		 ---  [   ]  ---'
        
               do jk=1,klev
        	   zdz=xg/(ppresh(jl,jk+1)-ppresh(jl,jk))
        	   write(8,17)jk,ppres(jl,jk)*1.e-2,pgeo(jl,jk)/xg,&
        	 & (ptten(jl,jk)-ztten(jl,jk))*xtjour,(pqten(jl,jk)-zrvten(jl,jk))*zeps,(prcten(jl,jk)+priten(jl,jk))*zeps, &
               ! & pumf(jl,jk),pdmf(jl,jk),(pprlflx(jl,jk)+pprsflx(jl,jk))*3.6e3,&
        	 & pumf(jl,jk),pdmf(jl,jk),pprlflx(jl,jk)*3.6e3,pprsflx(jl,jk)*3.6e3,&
        	 & purci(jl,jk)*1.e3,&
        	 & puten(jl,jk)*xtjour,pvten(jl,jk)*xtjour, &
        	 & pc(jl,jk,1),pc(jl,jk,1)+pcten(jl,jk,1)*pdtconv ,&
        	 & pc(jl,jk,2),pc(jl,jk,2)+pcten(jl,jk,2)*pdtconv 
              end do 
       
        else
        write(8,*)' %	  P	 Z	dT/dt	 dqv/dt d(ql)/dt Mflxup  Mflxdown&
        &  Prflx   rci_up     du/dt    dv/dt'
        write(8,*)' %  [hPa]	 [m]	    ---- [K/day] ----	     [kg/(sm^2]&
        &    [mm/h]   [g/kg]	 [m/s/day]'
        
               do jk=1,klev
        	   write(8,18)jk,ppres(jl,jk)*1.e-2,pgeo(jl,jk)/xg,				    &
        	 & (ptten(jl,jk)-ztten(jl,jk))*xtjour,(pqten(jl,jk)-zrvten(jl,jk))*zeps,(prcten(jl,jk)+priten(jl,jk))*zeps, &
        	 & pumf(jl,jk),pdmf(jl,jk),(pprlflx(jl,jk)+pprsflx(jl,jk))*3.6e3,&
        	 & purci(jl,jk)*1.e3,&
        	 & puten(jl,jk)*xtjour,pvten(jl,jk)*xtjour
              end do 
        
        endif
        !
        ! print rainfall tend (mm/h) -> transform from m/s to mm/h
       
        write(8,'(a37,2f8.3)')'%liquid and solid surf precip [mm/h]:', &
        		       &pprlflx(jl,klev+1)*3600.*1.e3, pprsflx(jl,klev+1)*3600.*1.e3
        write(8,'(a38,3i4,2f9.0)')'%conv-type cloud-top cloud-base CAPE CIN:',ktype(jl),kcltop(jl),kclbas(jl),&
                                                                        &pcape(jl),pcin(jl)

    end do

   !17 format(i3,f7.0,f8.0,2f9.3,f8.3,4f9.4,2f9.2,4e11.3)
   17 format(i3,f7.0,f8.0,2f9.3,f8.3,5f9.4,2f9.2,4e11.3)
   18 format(i3,f7.0,f8.0,2f9.3,f8.3,4f9.4,2f9.2)
   19 format(2f9.1,6f12.1,2f8.4)

 
     ! Print mean profile over the whole period     
      write(9,*)'% P[hPa]    Z[m]    dT/dt_conv   dr_t/dt_conv    dT/dt_obs    dr_t/dt_obs    du/dt    dv/dt    Mfl    Mfl_obs'
      do jk=1,klev
           write(9,19)sum(ppres(:,jk))/zlon*1.e-2,sum(pgeo(:,jk))/(xg*zlon), &
           sum(ptten(:,jk)-ztten(:,jk))/zlon*xtjour,sum(pqten(:,jk)-zrvten(:,jk)+prcten(:,jk))/zlon*zeps,&         
           sum(zq1(:,jk))/zlon,-sum(zq2(:,jk))/zlon, &                                     
           sum(puten(:,jk))/zlon*xtjour,sum(pvten(:,jk))/zlon*xtjour,  &
           sum(pumf(:,jk)+pdmf(:,jk))/zlon,-sum(pvervel(:,jk))/zlon
      end do
   stop 222
  ENDIF

ENDDO;ENDDO ! loop i-j
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#endif

end subroutine cupar_ifs
!*************************************************************************************
subroutine check (m1,tht,ath,rtt,artt)
implicit none
integer, intent(in) :: m1
real, dimension(m1), intent(in) :: tht,ath,rtt,artt

integer k
do k=1,m1
print*,"check",k, tht(k),ath(k),rtt(k),artt(k)
enddo
print*,"mx1",maxval(tht),maxval(ath),minval(tht),minval(ath)
print*,"mx2",maxval(rtt),maxval(artt),minval(rtt),minval(artt)

end subroutine check
