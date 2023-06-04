!###########################################################################
!  Advection scheme for the Runge-Kutta dynamical core 
!    (after Wicker, Skamarock, 2002, MWR)
!  Dec/2015 by Saulo Freitas (INPE), Michael Baldauf (DWD)
!###########################################################################
!module ModAdvectc_rk
  !private
  !public :: advectc_rk
!contains

!-this routine is called by timestep_rk 
!call advectc_rk('V',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!call advectc_rk('THETAIL',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!call advectc_rk('PI',mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
!
!
  subroutine advectc_rk(varn,mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum)
    !use var_tables, only: scalar_table 
    use grid_dims, only: maxgrds
    use mem_tend, only: tend
    use var_tables, only: num_scalar, scalar_tab
    use mem_grid, only: ngrid, nzpmax, grid_g, dtlt, if_adap, jdim, time, &
       zt, zm, dzm, dzt, hw4,itopo

    use mem_basic, only: basic_g
    use mem_chem1, only: nspecies_transported
    use mem_stilt, only: stilt_g,iexev

    implicit none
    include "i8.h"
    integer, intent(in) :: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,mynum
    character(len=*), intent(in) :: varn
    !- the 2 parameters below should be moved to be defined in "RAMSIN"  file
    integer, parameter :: pd_or_mnt_constraint = 0
    integer, parameter :: order_h=3,order_v=3
    !- 
    !real, intent(in) :: dtlt
    !integer, dimension(ngrids), intent(in) :: num_scalar
    !integer, intent(in) :: itopo
    !real, dimension(mzp), intent(in) :: dzt,dzm,hw4,zm,zt
    !real, dimension(mzp,mxp,myp), intent(inout) :: ut,vt,wt
    !real, dimension(mxp,myp), intent(in) :: dxt,dxu,dxv,dyt,dyu,dyv,rtgu,rtgv,rtgt,f13t,f23t
    !real, dimension(mxp,myp), intent(in) :: fmapt,fmapu,fmapv,fmapui,fmapvi
    !real, dimension(mzp,mxp,myp), intent(in) :: up,vp,wp,uc,vc,wc,dn0,dn0u,dn0v
    !type(scalar_table), intent(inout) :: scalar_tab(maxsclr,ngrids)
    !-
    !
    ! Local Variables
    integer(kind=i8) :: mxyzp
    integer :: i,j,k,n
    !new    real, pointer :: scalarp(:,:,:)
    !new    real, pointer :: scalart(:)
    real, pointer :: scalarp, scalart
    integer :: i_scl,is,js,ks
    !- scratchs (local arrays)
    real :: vt3da(mzp,mxp,myp)
    real :: vt3db(mzp,mxp,myp)
    real :: vt3dc(mzp,mxp,myp)
    real :: vt3dh(mzp,mxp,myp)
    real :: vt3dj(mzp,mxp,myp)
    real :: vt3dk(mzp,mxp,myp)
    real :: vctr1(mzp)
    real :: vctr2(mzp)
     
    real :: mfx_wind(mzp,mxp,myp)
    real :: mfy_wind(mzp,mxp,myp)
    real :: mfz_wind(mzp,mxp,myp)

    mxyzp = mxp * myp * mzp

    vt3da=0.0
    vt3db=0.0
    vt3dc=0.0
    vt3dh=0.0
    vt3dj=0.0
    vt3dk=0.0
    vctr1=0.0
    vctr2=0.0
    mfx_wind=0.0 ;mfy_wind=0.0 ;mfz_wind =0.0      

    if (trim(varn) .eq. 'V') then

       !print*,"1adv wind";call flush(6)
       ! Advect  U, V, and W
       ! input: mzp,mxp,myp,ia,iz,ja,jz,izu,jzv
       !        basic_g%uc,%vc,%wc,%dn0,%dn0u,%dn0v
       !        grid_g%dxt,%dxu,%dxv,%dyt,%dyu,%dyv,%rtgt,%rtgu,%rtgv
       !        grid_g%f13t,%f23t,%fmapt,%fmapu,%fmapv,%fmapui,%fmapvi
       ! output: tend%ut,%vt,%wt
       !
       !--------------- U-advect 
       is=1
       js=0
       ks=0

       call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
            ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
            ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
            ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
            ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
            ! 
            ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
            ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
            ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
            ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
            ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)
      
       !- need to check which wind is to be advected 
       !  (currently: basic_g(ngrid)%uc)
       call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%uc &
                     ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*U 
                     ,vt3db    & ! similar for v
                     ,vt3dc    & ! similar for sigma_dot 
                     ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                     ,mfy_wind & ! similar for v
                     ,mfz_wind & ! similar for sigma_dot 
                     !
                     ,tend%ut_rk   &
		     ,is,js,ks     &
                     ,pd_or_mnt_constraint &
		     ,order_h,order_v      &
                     ,dtlt                 &
                     )
       
       !--------------- V-advect 
       is=0
       js=1
       ks=0
   
       call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
            ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
            ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
            ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
            ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
            !
            ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
            ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
            ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
            ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
            ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)

       !- need to check which wind is to be advected
       call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%vc &
                     ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*V
                     ,vt3db    & ! similar for v
                     ,vt3dc    & ! similar for sigma_dot 
                     ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                     ,mfy_wind & ! similar for v
                     ,mfz_wind & ! similar for sigma_dot 
                     !
                     ,tend%vt_rk   &
		     ,is,js,ks     &
                     ,pd_or_mnt_constraint &
		     ,order_h,order_v      &
                     ,dtlt                 &
                     )
       
       
       !--------------- W-advect 
       is=0
       js=0
       ks=1

       call mf_wind(mzp,mxp,myp,ia,iz,ja,jz,izu,jzv,itopo,hw4,jdim,dzt,dzm  &
            ,basic_g(ngrid)%uc,basic_g(ngrid)%vc,basic_g(ngrid)%wc              &
            ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v&
            ,grid_g(ngrid)%dxt,grid_g(ngrid)%dxu,grid_g(ngrid)%dxv     &
            ,grid_g(ngrid)%dyt,grid_g(ngrid)%dyu,grid_g(ngrid)%dyv     &
            !
            ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv  &
            ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t,grid_g(ngrid)%fmapt &
            ,grid_g(ngrid)%fmapu ,grid_g(ngrid)%fmapv                  &
            ,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi                 &
            ,vt3da,vt3db,vt3dc,mfx_wind,mfy_wind,mfz_wind,is,js,ks)
       
       
       !- need to check which wind is to be advected
       call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%wc &
                     ,vt3da    & ! uc*dn0u*fmapui*rtgu = rhou*W
                     ,vt3db    & ! similar for v
                     ,vt3dc    & ! similar for sigma_dot 
                     ,mfx_wind & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                     ,mfy_wind & ! similar for v
                     ,mfz_wind & ! similar for sigma_dot 
                     !
                     ,tend%wt_rk   &
		     ,is,js,ks     &
                     ,pd_or_mnt_constraint&
		     ,order_h,order_v &
                     ,dtlt                 &
                     )
       
    endif  ! endif of varn .eq. 'V'

    if (trim(varn) .eq. 'T'   .or. trim(varn) .eq. "PI"      .or. &
        trim(varn) .eq. "THA" .or. trim(varn) .eq. "THETAIL" .or. &
	trim(varn) .eq. 'SCALAR'  ) THEN

       !-------------- scalars advect
       is=0
       js=0
       ks=0
       
       ! input: basic_g%up,%uc,%vp,%vc,%wp,%wc,%dn0,%dn0u,%dn0v
       !        grid_g%rtgt,%rtgu,%rtgv,%fmapt,%fmapui,%fmapvi,%f13t,%f23t,%dxu,%dyv,%dxt,%dyt
       !        scalar_tab%var_p, %var_t
       ! output: scalar_tab%var_t
       
       if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR' ) then
        !- combine the 2-time levels wind fields for tracers
        do j = 1,myp
          do i = 1,mxp
             do k = 1,mzp
                vt3da(k,i,j) = (basic_g(ngrid)%up(k,i,j)  &
                              + basic_g(ngrid)%uc(k,i,j)) * 0.5
                vt3db(k,i,j) = (basic_g(ngrid)%vp(k,i,j)  &
                              + basic_g(ngrid)%vc(k,i,j)) * 0.5
                vt3dc(k,i,j) = (basic_g(ngrid)%wp(k,i,j)  &
                              + basic_g(ngrid)%wc(k,i,j)) * 0.5
             end do
          end do
        end do
       else
        do j = 1,myp
          do i = 1,mxp
             do k = 1,mzp
                vt3da(k,i,j) = basic_g(ngrid)%uc(k,i,j) 
                vt3db(k,i,j) = basic_g(ngrid)%vc(k,i,j)
                vt3dc(k,i,j) = basic_g(ngrid)%wc(k,i,j)
             end do
          end do
        end do
       endif

       ! input: vt3da,vt3db,vt3dc
       !        basic_g%dn0,%dn0u,%dn0v
       !        grid_g%rtgt,%rtgu,%rtgv,%fmapt,%fmapui,%fmapvi,%f13t,%f23t,%dxu,%dyv,%dxt,%dyt
       ! output:vt3da,vt3db,vt3dc,vt3dh,vt3dj,vt3dk
       !
       call fa_preptc_rk(mzp,mxp,myp    &
                   ,vt3da,vt3db,vt3dc,vt3dh,vt3dj,vt3dk,vctr1,vctr2              &
                   ,basic_g(ngrid)%dn0,basic_g(ngrid)%dn0u,basic_g(ngrid)%dn0v   &
                   ,grid_g(ngrid)%rtgt,grid_g(ngrid)%rtgu,grid_g(ngrid)%rtgv     &
                   ,grid_g(ngrid)%fmapt,grid_g(ngrid)%fmapui,grid_g(ngrid)%fmapvi&
                   ,grid_g(ngrid)%f13t,grid_g(ngrid)%f23t                        &
                   ,grid_g(ngrid)%dxu,grid_g(ngrid)%dyv                          &
                   ,grid_g(ngrid)%dxt,grid_g(ngrid)%dyt,hw4,dzm,dzt,zm,zt)
      
       if ( trim(varn) .eq. "THA" .and. iexev == 2) THEN
           !-get log(thetav)
           call prep_lnthetv(mzp,mxp,myp,ia,iz,ja,jz&
	                    ,basic_g(ngrid)%theta &
		            ,basic_g(ngrid)%rtp   &
		            ,basic_g(ngrid)%rv    &
		            ,stilt_g(ngrid)%lnthetav)
          
	  call advect_ws(mzp,mxp,myp,ia,iz,ja,jz &
	              ,stilt_g(ngrid)%lnthetav   &! advected field
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U 
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot 
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot 
                      !
                      ,stilt_g(ngrid)%lnthvadv & !tendency field including advection
		      ,is,js,ks,pd_or_mnt_constraint&
                      ,order_h,order_v               &
                      ,dtlt                &
                      )

          return
       endif !endif og varn .eq. 'THA'
       if ( trim(varn) .eq. "PI" .and. iexev == 2) THEN
       
          call advect_ws(mzp,mxp,myp,ia,iz,ja,jz&
	              ,basic_g(ngrid)%pc & !advected field
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U 
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot 
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot 
                      !
                      ,tend%pt_rk          & !tendency from advection
		      ,is,js,ks            &
                      ,pd_or_mnt_constraint&
                      ,order_h,order_v               &
                      ,dtlt                &
                      )


          return
       endif !endif og varn .eq. 'PI'

       if ( trim(varn) .eq. "THETAIL" ) THEN
       
          call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,basic_g(ngrid)%thc &
                      ,vt3da & ! uc*dn0u*fmapui*rtgu = rhou*U 
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot 
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot 
                      !
                      ,tend%tht_rk         &  
		      ,is,js,ks            &
                      ,pd_or_mnt_constraint&
                      ,order_h,order_v     &
                      ,dtlt                &
                      )


          return
       endif !endif of varn .eq. 'THETAIL'
       
       if (trim(varn) .eq. 'T' .or. trim(varn) .eq. 'SCALAR') THEN

         i_scl=num_scalar(ngrid)  !- all scalars 
        !i_scl= num_scalar(ngrid) - NSPECIES_TRANSPORTED !- only theta_il+water+tke

         do n=1,i_scl
  
          !- if RK or ABM3 schemes, THP/THC are not transported here
          if (scalar_tab(n,ngrid)%name == 'THC' .or. &
              scalar_tab(n,ngrid)%name == 'THP') cycle

          !new      scalarp => scalar_tab(n,ngrid)%var_p_3D
          !new      scalart => scalar_tab(n,ngrid)%var_t_1D
          scalarp => scalar_tab(n,ngrid)%var_p
          scalart => scalar_tab(n,ngrid)%var_t

          ! input: scalarp, scalart, dtlt
          ! output: scalart

          call advect_ws(mzp,mxp,myp,ia,iz,ja,jz,scalarp &
                      ,vt3da & ! 0.5(up+uc)*dn0u*fmapui*rtgu = rhou*U 
                      ,vt3db & ! similar for v
                      ,vt3dc & ! similar for sigma_dot 
                      ,vt3dh & ! fmapt*rtgti*dxt/dn0 = 1(rho dx)
                      ,vt3dj & ! similar for v
                      ,vt3dk & ! similar for sigma_dot 
                      !
                      ,scalart,is,js,ks    &
                      ,pd_or_mnt_constraint&
                      ,order_h,order_v               &
                      ,dtlt                &
                      )

         end do

       endif  !endif og varn .eq. 'T'

    endif !endif of varn .eq. 'T' .or. varn .eq. "PI"

  end subroutine advectc_rk

!---------------------------------------------------------------------

  subroutine fa_preptc_rk(m1,m2,m3,vt3da,vt3db,vt3dc, &
       vt3dh,vt3dj,vt3dk,vctr1, vctr2,dn0,dn0u,dn0v,  &
       rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t,  &
       dxu,dyv,dxt,dyt,hw4, dzm, dzt, zm, zt)


    implicit none
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    real, intent(in) :: hw4(m1), dzm(m1), dzt(m1), zm(m1), zt(m1)
    real, intent(inout) :: vt3da(m1,m2,m3)
    real, intent(inout) :: vt3db(m1,m2,m3)
    real, intent(inout) :: vt3dc(m1,m2,m3)
    real, intent(out) :: vt3dh(m1,m2,m3)
    real, intent(out) :: vt3dj(m1,m2,m3)
    real, intent(out) :: vt3dk(m1,m2,m3)
    real, intent(in) :: dn0(m1,m2,m3)
    real, intent(in) :: dn0u(m1,m2,m3)
    real, intent(in) :: dn0v(m1,m2,m3)
    real, intent(in) :: rtgt(m2,m3)
    real, intent(in) :: rtgu(m2,m3)
    real, intent(in) :: rtgv(m2,m3)
    real, intent(in) :: fmapt(m2,m3)
    real, intent(in) :: fmapui(m2,m3)
    real, intent(in) :: fmapvi(m2,m3)
    real, intent(in) :: f13t(m2,m3)
    real, intent(in) :: f23t(m2,m3)
    real, intent(in) :: dxu(m2,m3)
    real, intent(in) :: dyv(m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(out) :: vctr1(m1)
    real, intent(out) :: vctr2(m2)

    ! Local Variables
    integer :: j,i,k,im,ip,jm,jp
    real :: c1,c2,c3,c4,rtgti
    real :: vctr3(m1)


    ! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
    ! between past and current time levels) times dtlt.

    ! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
    !    and include 1/rtgt factor in VT3DC

    do j = 1,m3
       jm = max(1,j-1)
       jp = min(m3,j+1)
       do i = 1,m2
          im = max(1,i-1)
          ip = min(m2,i+1)
          rtgti = 1. / rtgt(i,j)
          c1 = .5 * dxu(i,j)
          c2 = .5 * dyv(i,j)
          c3 = dxt(i,j) * fmapt(i,j) * rtgti
          c4 = dyt(i,j) * fmapt(i,j) * rtgti

          do k = 1,m1-1
             vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
                  + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
                  + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
                  + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
                  + vt3dc(k,i,j) * rtgti
             !vt3dd(k,i,j) = c1 * vt3da(k,i,j)
             !vt3de(k,i,j) = c2 * vt3db(k,i,j)
             !vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
             vctr3(k) = 1. / dn0(k,i,j)
             !-- these are the metric factors for u,v and z directions
	     !-- divided by ( air_dens times grid spacing).
	     vt3dh(k,i,j) = c3 * vctr3(k)
             vt3dj(k,i,j) = c4 * vctr3(k)
             vt3dk(k,i,j) = dzt(k) * vctr3(k)
          end do

          !            vt3di(1,i,j) = dxu(i,j) / (dxu(i,j) + dxt(ip,j))
          !            vt3di(2,i,j) = dxu(i,j) / (dxu(i,j) + dxt(i,j))
          !            vt3di(3,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,jp))
          !            vt3di(4,i,j) = dyv(i,j) / (dyv(i,j) + dyt(i,j))
          ! temporary override
          !vt3di(1,i,j) = .5
          !vt3di(2,i,j) = .5
          !vt3di(3,i,j) = .5
          !vt3di(4,i,j) = .5
       end do
    end do

    do k = 1,m1-1
       vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
       vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
    end do

    ! Convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
    ! into mass fluxes times dtlt.

    do j = 1,m3
       do i = 1,m2
          c1 = fmapui(i,j) * rtgu(i,j)
          c2 = fmapvi(i,j) * rtgv(i,j)
          do k = 1,m1-1
             vt3da(k,i,j) = vt3da(k,i,j) * c1 * dn0u(k,i,j)
             vt3db(k,i,j) = vt3db(k,i,j) * c2 * dn0v(k,i,j)
             vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
                  * (dn0(k,i,j) + dn0(k+1,i,j))
          end do
       end do
    end do
  end subroutine fa_preptc_rk

  !---------------------------------------------------------------------
  subroutine advect_ws(mzp,mxp,myp,ia,iz,ja,jz,scp,ufx,vfx,wfx &
                      ,vt3dh,vt3dj,vt3dk,sct,is,js,ks          &
                      ,pd_or_mnt_constraint,order_h,order_v,dt)

    implicit none
    integer, intent(in) :: mzp !- z
    integer, intent(in) :: mxp !- x
    integer, intent(in) :: myp !- y
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: is,js,ks,pd_or_mnt_constraint,order_h,order_v
    real, intent(in)    :: dt
    real, intent(in)    :: scp(mzp,mxp,myp)
    real, intent(inout) :: sct(mzp,mxp,myp)
    real, dimension(mzp,mxp,myp), intent(in) :: ufx, vfx,wfx,vt3dh,vt3dj,vt3dk
    ! (u,v,w)fx = (u,v,w) * air_dens / length
    !---local vars
    !
    !- parameters for various advection                                        
    real, parameter :: f30 = 7./12.
    real, parameter :: f31 = 1./12.
    real, parameter :: f40 = 7./12.
    real, parameter :: f41 = 1./12.
    real, parameter :: f50 = 37./60.
    real, parameter :: f51 = 2./15.
    real, parameter :: f52 = 1./60.
    real, parameter :: f60 = 37./60.
    real, parameter :: f61 = 2./15.
    real, parameter :: f62 = 1./60.
    real, parameter :: eps = 1.e-8
    
    logical :: scalar
    integer :: n,i,j,k
    real :: ue,uw,vs,vn,wbot,wtop
    real,dimension(mzp+ks,  mxp+is,  myp+js) :: qx,qy,qz,qxl,qyl,qzl
    
    real,dimension(-2:mzp+3,-2:mxp+3,-2:myp+3) :: scr
!    real,allocatable, dimension(:,:,:) :: scale_in,scale_out

    real,dimension(0:mzp  ,0:mxp+1,0:myp  ) :: ufx_local
    real,dimension(0:mzp  ,0:mxp  ,0:myp+1) :: vfx_local
    real,dimension(0:mzp+1,0:mxp  ,0:myp  ) :: wfx_local


    real :: fifth_order
    real :: cr,dir,div_term,scl_low,flux_out,scale                                   
    real :: flux_out_x,flux_out_y,flux_out_z
    real :: qim3, qim2, qim1, qi, qip1, qip2,qip3
    real :: fq2, fq3, fq4, fq5, fq6, fq, flux_upwind

    
    !- polynomial interpolation operators
    ! for 5th and 6th orders
    fq(qim2,qim1,qi,qip1,qip2,qip3,dir) = &       
                    f50*(qip1 + qi) - f51*(qip2 + qim1) + f52*(qip3 + qim2) &
     - fifth_order* f52*(qip3-qim2-5.*(qip2-qim1)+10.*(qip1-qi))*dir

    !--- 4th order interpolation operator    
    fq4(qim1,qi,qip1,qip2)= f40*(qip1+qi)-f41*(qip2+qim1)

    !--- 3rd order interpolation operator
    fq3(qim1,qi,qip1,qip2,dir)= f40*(qip1+qi)-f41*(qip2+qim1) &
                              - f41*(3.*(qip1-qi)-(qip2-qim1))*dir

    !--- 2nd order interpolation operator
    fq2(qi,qip1)= 0.5*(qip1+qi)

    !--- 1st order or upwind scheme
    flux_upwind(qi, qip1, dir ) = 0.5*(1.+dir)*qi - 0.5*(dir-1.)*qip1
   !-----------------------------------------------------------------------------

   !- flag to determine if a scalar is being advected
   scalar = .false.
   IF(is==0 .and. js==0 .and. ks==0) scalar = .true.
   !
   !- copy input arrays
   do j=1,myp-1+js
    do i=1,mxp-1+is
     do k=1,mzp-1+ks
       scr(k,i,j)=scp(k,i,j)
   enddo;enddo;enddo
   do j=1,myp
    do i=1,mxp
     do k=1,mzp
       ufx_local(k,i,j)=ufx(k,i,j)
       vfx_local(k,i,j)=vfx(k,i,j)
       wfx_local(k,i,j)=wfx(k,i,j)
   enddo;enddo;enddo
    
   !-----------------------------------------------------------------------------
   ! Set x & y boundary values in halo zones
   ! for now ( this setting is neeeded only for serial runs)

   ! direction X - 
   do n=1,3
    scr(:,       1-n,:)=scr(:,       1,:)
    scr(:,mxp-1+n+is,:)=scr(:,mxp-1+is,:)
   enddo
   ufx_local(:,    0,:)=ufx_local(:,  1,:)
   ufx_local(:,mxp+1,:)=ufx_local(:,mxp,:)
   
   vfx_local(:,    0,:)=vfx_local(:,  1,:)
   !vfx_local(:,mxp+1,:)=vfx_local(:,mxp,:)

   wfx_local(:,    0,:)=wfx_local(:,  1,:)
  !wfx_local(:,mxp+1,:)=wfx_local(:,mxp,:)
   
   ! direction Y - 
   do n = 1,3
    scr(:,:,       1-n)=scr(:,:,       1)
    scr(:,:,myp-1+n+js)=scr(:,:,myp-1+js)
   enddo
   ufx_local(:,:,    0)=ufx_local(:,:,  1)
   !ufx_local(:,:,myp+1)=ufx_local(:,:,myp)
   
   vfx_local(:,:,    0)=vfx_local(:,:,  1)
   vfx_local(:,:,myp+1)=vfx_local(:,:,myp)
   
   wfx_local(:,:,    0)=wfx_local(:,:,  1)
   !wfx_local(:,:,myp+1)=wfx_local(:,:,myp)
   
  
   ! direction z- 
   do n = 1,3
    scr(1-n       ,:,:)=scr(1       ,:,:)
    scr(mzp-1+n+ks,:,:)=scr(mzp-1+ks,:,:)
   enddo
   ufx_local(0    ,:,:)=ufx_local(1,  :,:)
   !ufx_local(mzp+1,:,:)=ufx_local(mzp,:,:)
  
   vfx_local(0    ,:,:)=vfx_local(1,  :,:)
   !vfx_local(mzp+1,:,:)=vfx_local(mzp,:,:)
  
   wfx_local(0    ,:,:)=wfx_local(1  ,:,:)
   wfx_local(mzp+1,:,:)=wfx_local(mzp,:,:)
      
   !
   !-----------------------------------------------------------------------------
   !- compute interface values of scalars/wind
   IF(order_h == 1 ) then
      !- compute x-interface values upwind order
      do j = 1,myp
        do i = 1,mxp-1
         do k = 1,mzp

           dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+is,j+js))      
           qx(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
      enddo; enddo; enddo
      !- compute y-interface values upwind order
      do j = 1,myp-1
        do i = 1,mxp
          do k = 1,mzp

           dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+is,j+js))      
           qy(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
      enddo; enddo; enddo
    ENDIF
    IF(order_v == 1 ) then  !- compute z-interface values upwind order
      do j = 1,myp
       do i = 1,mxp
        do k = 1,mzp-1
	   
           dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+is,j+js))      
           qz(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
      enddo; enddo; enddo
    ENDIF
   !-----------------------------------------------------------------------------
   IF(order_h == 2 ) then
      !- compute x-interface values
      do j = 1,myp
        do i = 1,mxp-1
         do k = 1,mzp

            qx(k,i,j) = fq2(scr(k,i,j),scr(k,i+1,j))
      enddo; enddo; enddo
      !- compute y-interface values
      do j = 1,myp-1
         do i = 1,mxp
          do k = 1,mzp

             qy(k,i,j) = fq2(scr(k,i,j),scr(k,i,j+1))
      enddo; enddo; enddo
   ENDIF
   IF(order_v == 2 ) then  !- compute z-interface values upwind order
      !- compute z-interface values
      do j = 1,myp
        do i = 1,mxp
         do k = 1,mzp-1
            qz(k,i,j) = fq2(scr(k,i,j),scr(k+1,i,j))
      enddo; enddo; enddo
   ENDIF
   !-----------------------------------------------------------------------------
   IF(order_h == 3 ) then
      !- compute x-interface values
      do j = 1,myp
        do i = 1,mxp-1
         do k = 1,mzp

          dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+is,j+js))      
          qx(k,i,j) = fq3(scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j),dir)
	  
      enddo; enddo; enddo
      !- compute y-interface values
      do j = 1,myp-1
       do i = 1,mxp
        do k = 1,mzp
      
          dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+is,j+js))     
          qy(k,i,j) = fq3(scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1), scr(k,i,j+2),dir)
      enddo; enddo; enddo
   ENDIF
   IF(order_v == 3 ) then
      !- compute z-interface values 
      do j = 1,myp
       do i = 1,mxp
         do k = 1,mzp-1

          dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+is,j+js))     
          qz(k,i,j) = fq3(scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j), scr(k+2,i,j),dir)      
      enddo; enddo; enddo          
   ENDIF
   !-----------------------------------------------------------------------------
   IF(order_h == 4 ) then
     !- compute x-interface values 
     do j = 1,myp
       do i = 1,mxp-1
        do k = 1,mzp
	
           qx(k,i,j) = fq4(scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j))
     enddo; enddo; enddo          
     !- compute y-interface values
     do j = 1,myp-1
       do i = 1,mxp
        do k = 1,mzp

           qy(k,i,j) = fq4(scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1),scr(k,i,j+2))  
     enddo; enddo; enddo          
   ENDIF
   IF(order_v == 4 ) then
     !- compute z-interface values 
     do j = 1,myp
       do i = 1,mxp
        do k = 1,mzp-1
	   
           qz(k,i,j) = fq4(scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j),scr(k+2,i,j))  
     enddo; enddo; enddo
   ENDIF
   !-----------------------------------------------------------------------------
   IF(order_h == 5 .or. order_h == 6 ) then
     fifth_order = 1.0                ! 5th order
     if(order_h == 6)fifth_order = 0.0  ! 6th order
     !- compute x-interface values 
     do j = 1,myp
       do i = 1,mxp-1
        do k = 1,mzp
           dir = SIGN(1.0,ufx_local(k,i,j)+ufx_local(k+ks,i+is,j+js))      
           qx(k,i,j) = fq(scr(k,i-2,j),scr(k,i-1,j),scr(k,i,j),scr(k,i+1,j),scr(k,i+2,j),scr(k,i+3,j),dir)
     enddo; enddo; enddo
     !- compute y-interface values 
     do j = 1,myp-1
       do i = 1,mxp
        do k = 1,mzp
           dir = SIGN(1.0,vfx_local(k,i,j)+vfx_local(k+ks,i+is,j+js))      
           qy(k,i,j) = fq(scr(k,i,j-2),scr(k,i,j-1),scr(k,i,j),scr(k,i,j+1),scr(k,i,j+2),scr(k,i,j+3),dir)
     enddo; enddo; enddo
   ENDIF
   IF(order_v == 5 .or. order_v == 6 ) then
     fifth_order = 1.0                ! 5th order
     if(order_v == 6)fifth_order = 0.0  ! 6th order
     !- compute z-interface values 
     do j = 1,myp
       do i = 1,mxp
        do k = 1,mzp-1
           dir = SIGN(1.0,wfx_local(k,i,j)+wfx_local(k+ks,i+is,j+js))      
           qz(k,i,j) = fq(scr(k-2,i,j),scr(k-1,i,j),scr(k,i,j),scr(k+1,i,j),scr(k+2,i,j),scr(k+3,i,j),dir)
     enddo; enddo; enddo
   ENDIF   
   !-----------------------------------------------------------------------------
   !- positivity/monotonicity constraints
   !   
   IF(pd_or_mnt_constraint > 0 .and. scalar) then 
      
      !-compute x,y,z-interfaces values with upwind scheme
      do j = 1,myp
        do i = 1,mxp-1
         do k = 1,mzp

            dir = SIGN(1.0,ufx_local(k,i,j))
            !- upwind flux (pd and monotonic)
            qxl(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i+1,j),dir)
            !- difference of high order and upwind fluxes   
            qx (k,i,j)=qx (k,i,j)-qxl(k,i,j)
      enddo; enddo; enddo
      !
      do j = 1,myp-1
        do i = 1,mxp
         do k = 1,mzp

            dir = SIGN(1.0,vfx_local(k,i,j))
            !- upwind flux (pd and monotonic)
            qyl(k,i,j)=flux_upwind(scr(k,i,j),scr(k,i,j+1),dir)
            !- difference of high order and upwind fluxes
            qy (k,i,j)=qy (k,i,j)-qyl(k,i,j)
      enddo; enddo; enddo
      !
      do j = 1,myp
        do i = 1,mxp
         do k = 1,mzp-1

            dir = SIGN(1.0,wfx_local(k,i,j))      
            !- upwind flux (pd and monotonic)
            qzl(k,i,j)=flux_upwind(scr(k,i,j),scr(k+1,i,j),dir)
            !- difference of high order and upwind fluxes   
            qz (k,i,j)=qz (k,i,j)-qzl(k,i,j)
      enddo; enddo; enddo

      !-- this section only imposes positivity constraint on scalars
      IF(pd_or_mnt_constraint == 1) then 

       do j = ia,iz
        do i = ja,jz
         do k = 2,mzp
          ue  =ufx(k,i  ,j)
          uw  =ufx(k,i-1,j)
          vn  =vfx(k,i,j  ) 
          vs  =vfx(k,i,j-1) 
          wtop=wfx(k  ,i,j)
          wbot=wfx(k-1,i,j) 
         
          div_term = scp(k,i,j)*(vt3dh(k,i,j)*(ue   - uw) + &
                                 vt3dj(k,i,j)*(vn   - vs) + &
                                 vt3dk(k,i,j)*(wtop - wbot))
 
          !- 1st order update
          scl_low = scp(k,i,j) +   &
                   dt*(- vt3dh(k,i,j)*(ue  *qxl(k,i,j) - uw  *qxl(k,i-1,j)) &
                       - vt3dj(k,i,j)*(vn  *qyl(k,i,j) - vs  *qyl(k,i,j-1)) &
                       - vt3dk(k,i,j)*(wtop*qzl(k,i,j) - wbot*qzl(k-1,i,j)) &
                       + div_term &
                      )

          !- net flux out
          flux_out_x = vt3dh(k,i,j)*(max(0.,ue  *qx(k,i,j)) - min(0.,uw  *qx(k,i-1,j)))
          flux_out_y = vt3dj(k,i,j)*(max(0.,vn  *qy(k,i,j)) - min(0.,vs  *qy(k,i,j-1)))
          flux_out_z = vt3dk(k,i,j)*(max(0.,wtop*qz(k,i,j)) - min(0.,wbot*qz(k-1,i,j)))
         
          flux_out = flux_out_x + flux_out_y + flux_out_z
          
          !- include divergence term
          flux_out = (flux_out - div_term)*dt
         
          !- re-scale the fluxes to keep positivity
          !
          IF( flux_out > scl_low ) THEN
             !- scale factor (=< 1.)
            scale = max(0.,scl_low/(flux_out+eps))
            !
            !-faces - x
            if( ue*qx(k,i  ,j) > 0.  ) qx(k,i  ,j)=scale*qx(k,i  ,j)
            if( uw*qx(k,i-1,j) < 0.  ) qx(k,i-1,j)=scale*qx(k,i-1,j)
            !-faces - y
            if( vn*qy(k,i  ,j) > 0.  ) qy(k,i  ,j)=scale*qy(k,i  ,j)
            if( vs*qy(k,i,j-1) < 0.  ) qy(k,i,j-1)=scale*qy(k,i,j-1)
            !-faces - z
            if( wtop*qz(k  ,i,j) > 0.) qz(k  ,i,j)=scale*qz(k  ,i,j)
            if( wbot*qz(k-1,i,j) < 0.) qz(k-1,i,j)=scale*qz(k-1,i,j)
           ENDIF
       enddo;enddo;enddo
      ENDIF ! endif for pd_or_mnt_constraint == 1
      
      !-- this section imposes monotonicity constraint on scalars
      !-- flux renormalization for monotonicity constraint following Durran (1998)
      IF(pd_or_mnt_constraint == 2) then 
      !- reserved 
      ENDIF
      
      !-- this section imposes monotonicity constraint on scalars
      !-- flux renormalization of Blossey&Durran 2008 + Skamarock 2006
      IF(pd_or_mnt_constraint == 3) then 
      !- reserved 
      ENDIF ! endif for pd_or_mnt_constraint == 3
      
      
      
      !-- get back the total "flux" , but now renormalized to guarantee 
      !-- positivity and/or monotonicity
      do j = 1,myp
        do i = 1,mxp
          do k = 1,mzp-1
           qx(k,i,j) = qx(k,i,j) + qxl(k,i,j) 
           qy(k,i,j) = qy(k,i,j) + qyl(k,i,j)
           qz(k,i,j) = qz(k,i,j) + qzl(k,i,j)
      enddo;enddo;enddo

   ENDIF ! endif for pd_or_mnt_constraint > 1 and scalar .true.   
   !
   !--------------------------------------------------------
   !- create the tendency
    do j = ja,jz
       do i = ia,iz
          do k = 2,mzp-1
          
            ue  =0.5*( ufx(k,i  ,j) + ufx(k+ks,i  +is,j+js) )
            uw  =0.5*( ufx(k,i-1,j) + ufx(k+ks,i-1+is,j+js) )

            vn  =0.5*( vfx(k,i,j  ) + vfx(k+ks,i+is,j  +js) ) 
            vs  =0.5*( vfx(k,i,j-1) + vfx(k+ks,i+is,j-1+js) ) 

            wtop=0.5*( wfx(k  ,i,j) + wfx(k+ks  ,i+is,j+js) )
            wbot=0.5*( wfx(k-1,i,j) + wfx(k-1+ks,i+is,j+js) ) 
          
            div_term = scp(k,i,j)*(vt3dh(k,i,j)*(ue   - uw  ) + &
                                   vt3dj(k,i,j)*(vn   - vs  ) + &
                                   vt3dk(k,i,j)*(wtop - wbot))
           
            sct(k,i,j) = sct(k,i,j) &
                       - vt3dh(k,i,j)*(ue  *qx(k,i,j) - uw  *qx(k,i-1,j)) &
                       - vt3dj(k,i,j)*(vn  *qy(k,i,j) - vs  *qy(k,i,j-1)) &
                       - vt3dk(k,i,j)*(wtop*qz(k,i,j) - wbot*qz(k-1,i,j)) &
                       + div_term
          
!            div_term = scp(k,i,j)*(vt3dh(k,i,j)*(ufx(k,i,j) - ufx(k,i-1,j)) + &
!                                      vt3dj(k,i,j)*(vfx(k,i,j) - vfx(k,i,j-1)) + &
!                                      vt3dk(k,i,j)*(wfx(k,i,j) - wfx(k-1,i,j))   )
           
!            sct(k,i,j) = sct(k,i,j) &
!                       - vt3dh(k,i,j)*(ufx(k,i,j)*qx(k,i,j) - ufx(k,i-1,j)*qx(k,i-1,j)) &
!                       - vt3dj(k,i,j)*(vfx(k,i,j)*qy(k,i,j) - vfx(k,i,j-1)*qy(k,i,j-1)) &
!                       - vt3dk(k,i,j)*(wfx(k,i,j)*qz(k,i,j) - wfx(k-1,i,j)*qz(k-1,i,j)) &
!                        + div_term
          
          
          !if(k==10 .and. i==20 .and. j==20)&
          !print*,"2",scp  (k,i,j), sct(k,i,j),qx(k,i,j),qy(k,i,j),qz(k,i,j)
          
       
          end do
       end do
    end do
    !print*,"mx-mn=",maxval(sct),minval(sct)

  end subroutine advect_ws
  !---------------------------------------------------------------------

  subroutine mf_wind(m1,m2,m3,ia,iz,ja,jz,izu,jzv,itopo, hw4, jdim, dzt, dzm,&
                     uc,vc,wc,dn0,dn0u,dn0v,                                 &
		     dxt,dxu,dxv,dyt,dyu,dyv,rtgt,rtgu,rtgv,f13t,f23t,       &
                     fmapt,fmapu,fmapv,fmapui,fmapvi,                        &
                     !
                     flxu,flxv,flxw, mfx_wind,mfy_wind,mfz_wind,is,js,ks)

    implicit none
    integer, intent(in) :: m1
    integer, intent(in) :: m2
    integer, intent(in) :: m3
    integer, intent(in) :: ia
    integer, intent(in) :: iz
    integer, intent(in) :: ja
    integer, intent(in) :: jz
    integer, intent(in) :: izu
    integer, intent(in) :: jzv
    integer, intent(in) :: itopo, jdim,is,js,ks
    real, intent(in) :: dzt(m1), dzm(m1), hw4(m1)
    real, intent(in) :: uc  (m1,m2,m3)
    real, intent(in) :: vc  (m1,m2,m3)
    real, intent(in) :: wc  (m1,m2,m3)
    real, intent(in) :: dn0 (m1,m2,m3)
    real, intent(in) :: dn0u(m1,m2,m3)
    real, intent(in) :: dn0v(m1,m2,m3)
    real, intent(in) :: dxt(m2,m3)
    real, intent(in) :: dxu(m2,m3)
    real, intent(in) :: dxv(m2,m3)
    real, intent(in) :: dyt(m2,m3)
    real, intent(in) :: dyu(m2,m3)
    real, intent(in) :: dyv(m2,m3)
    real, intent(in) :: rtgt(m2,m3)
    real, intent(in) :: rtgu(m2,m3)
    real, intent(in) :: rtgv(m2,m3)
    real, intent(in) :: f13t(m2,m3)
    real, intent(in) :: f23t(m2,m3)
    real, intent(in) :: fmapt(m2,m3)
    real, intent(in) :: fmapu(m2,m3)
    real, intent(in) :: fmapv(m2,m3)
    real, intent(in) :: fmapui(m2,m3)
    real, intent(in) :: fmapvi(m2,m3)

    real, dimension(m1,m2,m3),intent(inout) :: mfx_wind,mfy_wind,mfz_wind,&
                                               flxu,flxv,flxw

    ! Local Variables
    integer :: j,i,k,jm,im
    real :: c1z,c1x,c1y

    ! Compute momentum fluxes flxu, flxv, flxw
    !- only necessary one time per timestep
    !
    if(is == 1 .and. js == 0 .and. ks == 0) then

     do j = 1,m3
       do i = 1,m2
          do k = 1,m1
             flxu(k,i,j) = uc(k,i,j) * dn0u(k,i,j) * rtgu(i,j)  &
                         * fmapui(i,j)
             flxv(k,i,j) = vc(k,i,j) * dn0v(k,i,j) * rtgv(i,j)  &
                         * fmapvi(i,j)
          end do
       end do
     end do

     if(itopo.eq.0) then
       do j = 1,m3
          do i = 1,m2
             do k = 1,m1-1
                flxw(k,i,j) = wc(k,i,j)  &
                            * .5 * (dn0(k,i,j) + dn0(k+1,i,j))
             end do
          end do
       end do
     else
       do j = 1,m3
          jm = max(j-1,1)
          do i = 1,m2
             im = max(i-1,1)
             do k = 1,m1-1
                flxw(k,i,j) = wc(k,i,j)  &
                     * .5 * (dn0(k,i,j) + dn0(k+1,i,j))  &
                     + hw4(k) * ((flxu(k,i,j) + flxu(k+1,i,j)  &
                     + flxu(k,im,j) + flxu(k+1,im,j)) * f13t(i,j)  &
                     + (flxv(k,i,j) + flxv(k+1,i,j)  &
                     + flxv(k,i,jm) + flxv(k+1,i,jm)) * f23t(i,j))
             end do
          end do
       end do
     end if
    endif
    
    if(is == 1 .and. js == 0 .and. ks == 0) then
    !- Compute metric factors to compute U tendency
     do j = ja,jz
       do i = ia,izu
          c1z = 1. / rtgu(i,j)
          c1x = c1z * fmapu(i,j) * dxu(i,j)
          c1y = c1z * fmapu(i,j) * dyu(i,j)

          do k = 2,m1-1
              mfx_wind(k,i,j)= c1x / dn0u(k,i,j) 
          end do

          do k = 2,m1-1
            mfy_wind(k,i,j)= c1y / dn0u(k,i,j) 
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j)= c1z * dzt(k) / dn0u(k,i,j) 
          end do
       end do
     end do


    elseif(is == 0 .and. js == 1 .and. ks == 0) then
     !- Compute metric factors to compute V tendency
     do j = ja,jzv
       do i = ia,iz
          c1z = 1. / rtgv(i,j)
          c1x = c1z * fmapv(i,j) * dxv(i,j)
          c1y = c1z * fmapv(i,j) * dyv(i,j)

          do k = 2,m1-1
             mfx_wind(k,i,j) = c1x / dn0v(k,i,j) 
          end do

          do k = 2,m1-1 
            mfy_wind(k,i,j) = c1y / dn0v(k,i,j) 
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j)= c1z * dzt(k) / dn0v(k,i,j)
          end do
       end do
     end do

    elseif(is == 0 .and. js == 0 .and. ks == 1) then
     !- Compute metric factors to compute W tendency
     do j = ja,jz
       do i = ia,iz
          c1z = 1. / rtgt(i,j)
          c1x = c1z * fmapt(i,j) * dxt(i,j)
          c1y = c1z * fmapt(i,j) * dyt(i,j)

          do k = 2,m1-1
            mfx_wind(k,i,j) = c1x / (dn0(k,i,j) + dn0(k+1,i,j))
          end do

          do k = 2,m1-1
            mfy_wind(k,i,j) = c1y / (dn0(k,i,j) + dn0(k+1,i,j))
          end do

          do k = 2,m1-1
            mfz_wind(k,i,j) = c1z * dzm(k) / (dn0(k,i,j) + dn0(k+1,i,j)) 
          end do
       end do
     end do
    endif
  end subroutine mf_wind

!     *********************************************************************
!end module ModAdvectc_rk
