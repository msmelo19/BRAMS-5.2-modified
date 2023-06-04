module dtset

  use grid_dims, only : &
       maxgrds           ! INTENT(IN)

  implicit none

  real    :: ssodx(maxgrds)
  integer :: idelx(maxgrds)
  real    :: sscourn(maxgrds)

contains

  subroutine dtset_new(mynum, nndtflg, dxtmax_local)

    use mem_grid, only : &
         ideltat,        & ! INTENT(IN)
         ngrids,         & ! INTENT(IN)
         deltaxn,        & ! INTENT(IN)
         nnxp,           & ! INTENT(IN)
         nnyp,           & ! INTENT(IN)
         nnzp,           & ! INTENT(IN)
         grid_g,         & ! INTENT(IN)
         sspct,          & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nxtnest,        & ! INTENT(IN)
         dtlongn,        & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         dtlong,         & ! INTENT(IN)    ! Falta passar para o escravo-ok
         nndtrat,        & ! INTENT(INOUT) ! Valor inicial passado-ok
         ! Calculado atualizacao nesta rot.
         ! Retirar da com.:Mest->escravo
         nnacoust,       & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nacoust,        & ! INTENT(IN)
         cflxy,          & ! INTENT(IN)    ! Calculado localmente em CFL
         cflz,           & ! INTENT(IN)    ! A ser calculado em modsched local
         iflag,          & ! INTENT(INOUT) ! definida localmente nesta rotina
         timmax,         & ! INTENT(IN)
         time,           & ! INTENT(IN)
	 dyncore_flag
    
    ! "cflxy" e "cflz" calculado localmente, precisa ser passado a todos os 
    ! escravos para se determinar o maior valor.
    
    use rconstants, only : &
         cp,               & ! INTENT(IN)
         cv,               & ! INTENT(IN)
         rgas                ! INTENT(IN)
    
    use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)
    ! Receber apos atualizacao-chamada de varf_read em RAMS_OUTPUT
    
    use grid_dims, only : &
         maxgrds,         & ! INTENT(IN)
         nzpmax             ! INTENT(IN)

    use io_params, only : &
         frqanl             ! INTENT(IN)
    
    use mem_stilt, only : &
         iexev             ! INTENT(IN)
    
    implicit none
    
    ! Arguments:
    integer, intent(in)  :: mynum
    integer, intent(out) :: nndtflg
    real, intent(in)     :: dxtmax_local(maxgrds)
    
    ! Local variables:
    integer, parameter  :: ndx=37,ndt=42
    integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
    real, dimension(maxgrds) :: dtlongn1
    real, dimension(nzpmax) :: vctr1
    real :: cflnumh, cflnumv, delx(ndx), delt(ndt)

    integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii
    real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft
    
    real :: dxta, dxtb, dxtc, dxtd
    character(len=8) :: c0, c1, c2, c3, c4
    character(len=*), parameter :: h="**(dtset_new)**"

    delx(:) = (/ &
         200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
         030000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
         001000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
         000150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
         000020.,    10.,     8.,    6.,    5.,    4.,    3.,  &
         000002.,     1. /)

    delt(:) = (/ &
         300.,  240.,   180.,  150.,  120.,   90.,   60.,  &
         050.,   40.,    30.,   20.,   15.,   12.,   10.,  &
         006.,    5.,     4.,    3.,   2.5,   2.0,   1.5,  &
         01.2,    1.0,     .8,    .6,   .5,    .4,    .3,  &
         00.2,     .1,     .08,   .06,  .05,   .04,   .03,  &
         00.02,    .01,    .008,  .006, .005,  .004,  .003/)

    idelx(1) = 0
    
    iabsdt = abs(ideltat)
    
    ! On the first call to this subroutine, initialize idelx, ssodx, dtlongn,
    ! nnacoust, and if required, nndtrat.
    
    if ( idelx(1)==0 ) then
       
       cflnumh = .90
       cflnumv = .90
       
       do ifm = 1,ngrids
          
          do id = ndx,1,-1
             if (delx(id)<=deltaxn(ifm))  idelx(ifm) = id
          enddo
          
          n2 = nnxp(ifm)
          n3 = nnyp(ifm)
          do k = 1,nnzp(ifm)
             vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
          enddo
          tmax = maxval(vctr1(1:nnzp(ifm)))
          ssmax = sqrt(cp / cv * rgas * tmax)
          
          nn2 = nnxp(ifm)
          nn3 = nnyp(ifm)

          if (mynum==0) then
             ! Master process
             dxtmax = max(grid_g(ifm)%dxt(1,1), &
                  grid_g(ifm)%dxt(nn2,1),       &
                  grid_g(ifm)%dxt(nn2,nn3),     &
                  grid_g(ifm)%dxt(1,nn3)        )
          else
             ! Slave (node) process
             dxtmax = dxtmax_local(ifm)
          endif
          
          ssodx(ifm) = ssmax * dxtmax
          
       enddo
       
       if ( ideltat==0 ) then
          
          sspct = 1.
          do ifm = 1,ngrids
             icm = nxtnest(ifm)
             if ( icm==0 ) then
                dtlongn(ifm) = dtlong
             else
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
             nnacoust(ifm) = nacoust
             sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
             sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          enddo
          
       else
          
          sscnmax = 0.
          do ifm = 1,ngrids
             icm = nxtnest(ifm)
             dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)
             
             ! For coarse grid(s) adjust dtlongn so that it is an integer
             ! divisor of FRQANL.  For nested grids, compute nndtrat(ifm) as the
             ! first integer greater than or equal to the timestep ratio between
             ! a nested grid's parent and the nested grid. Then adjust 
             ! dtlongn(ifm) for the nested grid to be the parent grid timestep 
             ! divided by nndtrat(ifm).
             
             if ( icm==0 ) then
                ntf = nint(frqanl / dtlongn(1))
                dtlongn(ifm) = frqanl / ntf
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
             
             ! Compute sst courant numbers (sstcourn(ifm)) for long timestep
             ! dtlongn.
             
             sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
             if ( sscourn(ifm)>sscnmax)  sscnmax = sscourn(ifm)
             
          enddo
          
          ! Define trial sspct0 in terms of sscnmax using nonlinear formula 
          ! intended to increase nnacoust as sspct decreases, but limit sspct0
          ! to a minimum of .2.
          
          sspct0 = min(1., (.95/sscnmax)**.5)
          
          if ( sspct0<.2 ) then
             print*, 'Sound speed percent is forced to be too low'
             stop 'low_sspct0'
          endif
          
          sspct = 1.
          do ifm = 1,ngrids
             nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
             sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          enddo
          
       endif
       
    endif
 
    ! check Courant numbers
    
    nndtflg = 0
    
    if ( ideltat>=0 ) then
       
      !if ( dyncore_flag == 0 .or. dyncore_flag == 1) then
       do ifm = 1,ngrids
          cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
          if ( cflxyz>1. ) then
             iflag = 1
	     !print*,'cfl z e xy= ',cflxy(ifm),cflz(ifm),mynum,iflag
             write(c0,"(f8.2)") cflxyz
             write(c3,"(f8.2)") cflxy(ifm)
             write(c4,"(f8.2)") cflz(ifm)
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             if ( dyncore_flag == 2 .or. dyncore_flag == 3) then
	       ! don't break here, there exist a better suited criterion in subroutine 'cfl1'
      
             else
               call fatal_error(h//" Model will stop because CFL limit exceeded:"//& 
                   &" clfxy="//trim(adjustl(c0))//" at proc "//trim(adjustl(c1))//&
                   &" on grid "//trim(adjustl(c2))//" cflxy="//trim(adjustl(c3))//&
	           &" cflz="//trim(adjustl(c4)))
             endif
	   end if
       enddo
      !endif
      !if ( dyncore_flag == 2 .or. dyncore_flag == 3) then
	! don't break here, there exist a better suited criterion in subroutine 'cfl1'
      !endif       
    else
       
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          cflxyz = max(cflxy(ifm)/cflnumh, cflz(ifm)/cflnumv)
          
          nndtrat1(ifm) = nndtrat(ifm)
          dtlongn1(ifm) = dtlongn(ifm)
          nnacoust1(ifm) = nnacoust(ifm)
          
          do id = ndt,1,-1
             if ( delt(id)*cflxyz<=dtlongn(ifm) ) idelt(ifm) = id
          enddo
          
          if ( idelt(ifm)>idelx(ifm)+4 ) then
             iflag = 1
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             call fatal_error(h//" Model will stop because adjustable timestep"//& 
                  &" forced to become too small at proc "//trim(adjustl(c1))//&
                  &" on grid "//trim(adjustl(c2)))
          else
             ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
             dtlongn(ifm) = delt(ii)
             
             ! For the coarse grid(s), adjust dtlongn(1) to not jump past an
             ! analysis write time or the end of the model run time.
             
             if ( icm==0 ) then
                timeleft =  &
                     min (timmax - time, frqanl - mod(time,frqanl))
                if ( dtlongn(1)>.95 * timeleft) then
                   dtlongn(ifm) = timeleft
                endif
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
          endif
          
          ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.
          
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)
       enddo
       
       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to 
       ! a minimum of .2.
       
       sspct0 = min(1., (.95/sscnmax) ** .5)
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif
       
       sspct = 1.
       do ifm = 1,ngrids
          !srf - forçando nacoust > 2 para iexev=2
          !nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==1)  nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==2)  nnacoust(ifm) = max(3, nint(sspct0 * sscourn(ifm) * 2. / .95))
          !srf

          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          
          ! If there are any updates to dtlongn(ifm), print out new values.
          
          if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
            write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),mynum,dtlongn(ifm),dtlongn1(ifm)
122          format('Timestep update1: ngrid, nndtrat, nnacoust, mynum,'  &
                  ,'new dtlongn, old dtlongn = ',4i4,2f13.6)
          endif
          
          ! If there are any updates to nndtrat(ifm) or others, set 
          ! nndtflg = 1 to flag new call to subroutine modsched and send new 
          ! stuff to nodes.
          
          if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
               nnacoust(ifm) /= nnacoust1(ifm) .or. &
               dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1
          
       enddo
       
    endif
    
    return
  end subroutine dtset_new

  ! ******************************************************************

  subroutine dump_courn()

    use mem_grid, only : &
         ngrids,         & ! INTENT(IN)
         dtlongn,        & ! INTENT(IN)
         nndtrat,        & ! INTENT(IN)
         nnacoust,       & ! INTENT(IN)
         sspct             ! INTENT(IN)

    implicit none

    integer :: ifm

    ! Print out initial values of dtlongn, nndtrat, nnacoust, sscourn, 
    ! and sspct
    write(*,"(a)") ' === Initial timestep info: ngrid, nndtrat, nnacoust,'//&
         '     dtlongn, sscourn,  sspct ==='
    do ifm = 1,ngrids
       write(*,"(29x,i3,1x,i8,1x,i9,f13.3,1x,2f8.3)") &
            ifm, nndtrat(ifm), nnacoust(ifm), dtlongn(ifm), sscourn(ifm), sspct
    enddo
    
  end subroutine dump_courn





!=================================================================




  subroutine SetDt(mynum, nndtflg, dxtmax_local)

    use mem_grid, only : &
         ideltat,        & ! INTENT(IN)
         ngrids,         & ! INTENT(IN)
         deltaxn,        & ! INTENT(IN)
         nnxp,           & ! INTENT(IN)
         nnyp,           & ! INTENT(IN)
         nnzp,           & ! INTENT(IN)
         grid_g,         & ! INTENT(IN)
         sspct,          & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nxtnest,        & ! INTENT(IN)
         dtlongn,        & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         dtlong,         & ! INTENT(IN)    ! Falta passar para o escravo-ok
         nndtrat,        & ! INTENT(INOUT) ! Valor inicial passado-ok
         ! Calculado atualizacao nesta rot.
         ! Retirar da com.:Mest->escravo
         nnacoust,       & ! INTENT(INOUT) ! Calculada localmente nesta rotina
         nacoust,        & ! INTENT(IN)
         cflxy,          & ! INTENT(IN)    ! Calculado localmente em CFL
         cflz,           & ! INTENT(IN)    ! A ser calculado em modsched local
         iflag,          & ! INTENT(INOUT) ! definida localmente nesta rotina
         timmax,         & ! INTENT(IN)
         time              ! INTENT(IN)
    
    ! "cflxy" e "cflz" calculado localmente, precisa ser passado a todos os 
    ! escravos para se determinar o maior valor.
    
    use rconstants, only : &
         cp,               & ! INTENT(IN)
         cv,               & ! INTENT(IN)
         rgas                ! INTENT(IN)
    
    use ref_sounding, only : &
         th01dn,             & ! INTENT(IN)
         pi01dn                ! INTENT(IN)
    ! Receber apos atualizacao-chamada de varf_read em RAMS_OUTPUT
    
    use grid_dims, only : &
         maxgrds,         & ! INTENT(IN)
         nzpmax             ! INTENT(IN)

    use io_params, only : &
         frqanl             ! INTENT(IN)
    

    use node_mod, only: &
         nodemxp, nodemyp

    use mem_stilt, only : &
         iexev             ! INTENT(IN)

    implicit none
    
    ! Arguments:
    integer, intent(in)  :: mynum
    integer, intent(out) :: nndtflg
    real, intent(in)     :: dxtmax_local(maxgrds)
    
    ! Local variables:
    integer, parameter  :: ndx=37,ndt=42
    integer, dimension(maxgrds) :: idelt,nndtrat1,nnacoust1
    real, dimension(maxgrds) :: dtlongn1
    real, dimension(nzpmax) :: vctr1
    real :: cflnumh, cflnumv, delx(ndx), delt(ndt)

    integer :: iabsdt,ifm,id,n2,n3,k,nn2,nn3,icm,ntf,ii
    real :: ssmax,tmax,dxtmax,sscnmax,sspct0,cflxyz,timeleft
    
    real :: dxta, dxtb, dxtc, dxtd
    character(len=8) :: c0, c1, c2
    character(len=*), parameter :: h="**(dtset_new)**"

    delx(:) = (/ &
         200000.,150000.,100000.,80000.,70000.,60000.,40000.,  &
         030000., 20000., 10000., 6000., 4000., 3000., 2000.,  &
         001000.,   800.,   600.,  500.,  400.,  300.,  200.,  &
         000150.,   100.,    80.,   60.,   50.,   40.,   30.,  &
         000020.,    10.,     8.,    6.,    5.,    4.,    3.,  &
         000002.,     1. /)

    delt(:) = (/ &
         300.,  240.,   180.,  150.,  120.,   90.,   60.,  &
         050.,   40.,    30.,   20.,   15.,   12.,   10.,  &
         006.,    5.,     4.,    3.,   2.5,   2.0,   1.5,  &
         01.2,    1.0,     .8,    .6,   .5,    .4,    .3,  &
         00.2,     .1,     .08,   .06,  .05,   .04,   .03,  &
         00.02,    .01,    .008,  .006, .005,  .004,  .003/)

    idelx(1) = 0
    
    iabsdt = abs(ideltat)
    
    ! initialize idelx, ssodx, dtlongn, nnacoust, and if required, nndtrat.
    
       
    cflnumh = .90
    cflnumv = .90
    
    do ifm = 1,ngrids
       
       do id = ndx,1,-1
          if (delx(id)<=deltaxn(ifm))  idelx(ifm) = id
       enddo
       
       n2 = nodemxp(mynum,ifm)
       n3 = nodemyp(mynum,ifm)
       do k = 1,nnzp(ifm)
          vctr1(k) = th01dn(k,1) * pi01dn(k,1) / cp
       enddo
       tmax = maxval(vctr1(1:nnzp(ifm)))
       ssmax = sqrt(cp / cv * rgas * tmax)
       
       nn2 = nodemxp(mynum,ifm)
       nn3 = nodemyp(mynum,ifm)
       
       dxtmax = max(grid_g(ifm)%dxt(1,1), &
            grid_g(ifm)%dxt(nn2,1),       &
            grid_g(ifm)%dxt(nn2,nn3),     &
            grid_g(ifm)%dxt(1,nn3)        )
       
       ssodx(ifm) = ssmax * dxtmax
       
    enddo
    
    if ( ideltat==0 ) then
       
       sspct = 1.
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          if ( icm==0 ) then
             dtlongn(ifm) = dtlong
          else
             dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
          endif
          nnacoust(ifm) = nacoust
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
       enddo
       
    else
       
       sscnmax = 0.
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          dtlongn(ifm) = delt(idelx(ifm)+iabsdt-1)
          
          ! For coarse grid(s) adjust dtlongn so that it is an integer
          ! divisor of FRQANL.  For nested grids, compute nndtrat(ifm) as the
          ! first integer greater than or equal to the timestep ratio between
          ! a nested grid's parent and the nested grid. Then adjust 
          ! dtlongn(ifm) for the nested grid to be the parent grid timestep 
          ! divided by nndtrat(ifm).
          
          if ( icm==0 ) then
             ntf = nint(frqanl / dtlongn(1))
             dtlongn(ifm) = frqanl / ntf
          else
             nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
             dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
          endif
          
          ! Compute sst courant numbers (sstcourn(ifm)) for long timestep
          ! dtlongn.
          
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if ( sscourn(ifm)>sscnmax)  sscnmax = sscourn(ifm)
          
       enddo
       
       ! Define trial sspct0 in terms of sscnmax using nonlinear formula 
       ! intended to increase nnacoust as sspct decreases, but limit sspct0
       ! to a minimum of .2.
       
       sspct0 = min(1., (.95/sscnmax)**.5)
       
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif
       
       sspct = 1.
       do ifm = 1,ngrids
          nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
       enddo
       
    endif
    
    
    ! check Courant numbers
    
    nndtflg = 0
    
    if ( ideltat>=0 ) then
       
       do ifm = 1,ngrids
          cflxyz = max(cflxy(ifm)/cflnumh,cflz(ifm)/cflnumv)
          if ( cflxyz>1. ) then
             iflag = 1
             write(c0,"(f8.2)") cflxyz
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             call fatal_error(h//" Model will stop because CFL limit exceeded:"//& 
                  &" clfxy="//trim(adjustl(c0))//" at proc "//trim(adjustl(c1))//&
                  &" on grid "//trim(adjustl(c2)))
          end if
       enddo
       
    else
       
       do ifm = 1,ngrids
          icm = nxtnest(ifm)
          cflxyz = max(cflxy(ifm)/cflnumh, cflz(ifm)/cflnumv)
          
          nndtrat1(ifm) = nndtrat(ifm)
          dtlongn1(ifm) = dtlongn(ifm)
          nnacoust1(ifm) = nnacoust(ifm)
          
          do id = ndt,1,-1
             if ( delt(id)*cflxyz<=dtlongn(ifm) ) idelt(ifm) = id
          enddo
          
          if ( idelt(ifm)>idelx(ifm)+4 ) then
             iflag = 1
             write(c1,"(i8)") mynum
             write(c2,"(i8)") ifm
             call fatal_error(h//" Model will stop because adjustable timestep"//& 
                  &" forced to become too small at proc "//trim(adjustl(c1))//&
                  &" on grid "//trim(adjustl(c2)))
          else
             ii = max(idelx(ifm)+iabsdt-1,idelt(ifm))
             dtlongn(ifm) = delt(ii)
             
             ! For the coarse grid(s), adjust dtlongn(1) to not jump past an
             ! analysis write time or the end of the model run time.
             
             if ( icm==0 ) then
                timeleft =  &
                     min (timmax - time, frqanl - mod(time,frqanl))
                if ( dtlongn(1)>.95 * timeleft) then
                   dtlongn(ifm) = timeleft
                endif
             else
                nndtrat(ifm) = min(10, nint(dtlongn(icm) / dtlongn(ifm)))
                dtlongn(ifm) = dtlongn(icm) / nndtrat(ifm)
             endif
          endif
          
          ! Compute sst courant numbers (sstcourn(ifm))for long timestep dtlongn.
          
          sscourn(ifm) = 2. * ssodx(ifm) * dtlongn(ifm)
          if (sscourn(ifm) .gt. sscnmax) sscnmax = sscourn(ifm)
       enddo
       
       ! Define trial sspct0 in terms of sscnmax using nonlinear formula
       ! intended to increase nnacoust as sspct decreases, but limit sspct0 to 
       ! a minimum of .2.
       
       sspct0 = min(1., (.95/sscnmax) ** .5)
       if ( sspct0<.2 ) then
          print*, 'Sound speed percent is forced to be too low'
          stop 'low_sspct0'
       endif
       
       sspct = 1.
       do ifm = 1,ngrids
          !srf - forçando nacoust > 2 para iexev=2
          !nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==1)  nnacoust(ifm) = max(2, nint(sspct0 * sscourn(ifm) * 2. / .95))
          if(iexev==2)  nnacoust(ifm) = max(3, nint(sspct0 * sscourn(ifm) * 2. / .95))
          !srf

          sspct = min(sspct, .95*float(nnacoust(ifm))/(2.*sscourn(ifm)))
          
          ! If there are any updates to dtlongn(ifm), print out new values.
          
          if (abs(dtlongn(ifm) - dtlongn1(ifm)) > 1.e-3) then
             write(6,122) ifm,nndtrat(ifm),nnacoust(ifm),dtlongn(ifm)
122          format('Timestep update: ngrid, nndtrat, nnacoust,'  &
                  ,' dtlongn = ',3i3,f10.3)
          endif
          
          ! If there are any updates to nndtrat(ifm) or others, set 
          ! nndtflg = 1 to flag new call to subroutine modsched and send new 
          ! stuff to nodes.
          
          if (nndtrat(ifm) /= nndtrat1(ifm) .or. &
               nnacoust(ifm) /= nnacoust1(ifm) .or. &
               dtlongn(ifm) /= dtlongn1(ifm) ) nndtflg = 1
          
       enddo
       
    endif
    
    return
  end subroutine SetDT
end module dtset
