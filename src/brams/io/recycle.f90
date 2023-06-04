!############################# Change Log ##################################
! 5.0.0
!
!###########################################################################
!  Copyright (C)  1990, 1995, 1999, 2000, 2003 - All Rights Reserved
!  Regional Atmospheric Modeling System - RAMS
!###########################################################################

subroutine recycle()

  use mem_grid, only: &
       ngrids, nnxp, nnyp, nzg, npatch, nzs, nnzp
  use var_tables, only : vtab_r, num_var
  use io_params, only: pastfn
  !srf -for carma AOT recycle
  use mem_aerad, only: nwave !INTENT(IN)
  use node_mod, only: &
       mzp, mxp, myp  ! INTENT(IN)
  use an_header,only: &
    anal_table,     &
    nvbtab
  use node_mod, only: &
       mzp, mxp, myp,  & ! INTENT(IN)
       izu, jzv,       & ! INTENT(IN)
       mynum,          & ! INTENT(IN)
       ibcon,          & ! INTENT(IN)
       nmachs,         &   ! INTENT(IN)
       mchnum, &
       master_num, &
       nodemxp, &
       nodemyp, &
       nodemzp, &
       nodei0, &
       nodej0
  USE ReadBcst, ONLY: &
       Broadcast
  use ParLib, only: &
       parf_bcast
 

  implicit none
  
  include "files.h"
  include "i8.h"
  
  integer, parameter :: i64 = selected_int_kind(14) !Kind for 64-bits Integer Numbers
  character(len=f_name_length) :: flnm
  !integer :: idtype,ng,nvars,lenf
  integer :: ng,nvars,lenf
  integer :: ierr
  real, allocatable :: scr1(:)
  real, allocatable :: scr2(:)
  character(len=f_name_length) :: flng
  integer(kind=i8)             :: npts
  integer(kind=i8)             :: fPosition
  character(len=*), parameter :: h="**(recycle)**"
  INTEGER :: npoints
  TYPE scridim
      real, allocatable :: scr(:,:,:,:)
  END TYPE scridim
  TYPE (scridim), DIMENSION(7) :: srcRead
  INTEGER :: ia,iz,ja,jz
  INTEGER :: m1,m2,m3
  INTEGER(kind=i64)  :: nzl,nxl,nyl,n4,j,k
 
  allocate (srcRead(2)%scr(1,nnxp(1),nnyp(1),1))
  allocate (srcRead(3)%scr(nnzp(1),nnxp(1),nnyp(1),1))
  allocate (srcRead(4)%scr(nzg,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(5)%scr(nzs,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(6)%scr(1,nnxp(1),nnyp(1),npatch))
  allocate (srcRead(7)%scr(1,nnxp(1),nnyp(1),nwave))

  flnm=pastfn(1:len_trim(pastfn)-9)
  
  !print*,'Reading assimilation fields from analysis file ',pastfn
  
  call rams_read_header(flnm(1:len_trim(flnm)))
  
  ! read the requested analysis file variables  
  !IF(mchnum==master_num) THEN
  !   DO nvars=1,nvbtab
  !      print*,'header:', anal_table(nvars)%string &
  !              , ' for grid:', anal_table(nvars)%ngrid   &
  !              , ' dim:', anal_table(nvars)%idim_type   &
  !              , ' npts:', anal_table(nvars)%nvalues 
  !    END DO
  ! END IF
 
  ia = nodei0(mynum,1)+1
  iz = nodei0(mynum,1)+nodemxp(mynum,1)
  ja = nodej0(mynum,1)+1
  jz = nodej0(mynum,1)+nodemyp(mynum,1)
  m1=nodemzp(mynum,1)
  m2=nodemxp(mynum,1)
  m3=nodemyp(mynum,1)

  do ng=1,ngrids
     
     DO nvars=1,num_var(ng)
        
        if(vtab_r(nvars,ng)%irecycle == 1) then

           !print*,'Reading assimilation field (from vtab):', vtab_r(nvars,ng)%name &
           !     , ' for grid:', ng                                     &
           !     , ' dim:', vtab_r(nvars,ng)%idim_type                  &
           !     , ' npts:', vtab_r(nvars,ng)%npts

           IF(mchnum==master_num) call FindFieldInAnalysisFile(vtab_r(nvars,ng)%name, ng,   &
                flnm(1:len_trim(flnm)), flng, npts, fPosition)
           npoints=npts
           CALL Broadcast(npoints, master_num, "npoints")  
           npts=npoints        
                
           allocate(scr1(npts), stat=ierr)
           if (ierr /= 0) call fatal_error(h//" allocating scr1")
           allocate(scr2(npts), stat=ierr)
           if (ierr /= 0) call fatal_error(h//" allocating scr2")

           IF(mchnum==master_num) call GetFieldInAnalysisFile(flng, npts, fPosition, scr1, scr2)
           
           CALL Broadcast(scr1, master_num, "scr1")

           if(vtab_r(nvars,ng)%idim_type == 4) then

              call unarrange_p(nnxp(ng),nnyp(ng),nzg,npatch  &
                   ,scr1(1),srcRead(4)%scr)
                   
           elseif(vtab_r(nvars,ng)%idim_type == 5) then
              call unarrange_p(nnxp(ng),nnyp(ng),nzs,npatch  &
                   ,scr1(1),srcRead(5)%scr)
              !srf
              !use this for 3d atmospheric fields
           elseif(vtab_r(nvars,ng)%idim_type == 3) then
              call unarrange(nnzp(ng),nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(3)%scr(:,:,:,1))

           elseif(vtab_r(nvars,ng)%idim_type == 6) then
              call rearrange_aot(npatch,nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(6)%scr(1,:,:,:))
              !srf
              !use this for 3d (NX,NY,NWAVE) CARMA AOT fields
           elseif(vtab_r(nvars,ng)%idim_type == 7) then
              call rearrange_aot(nwave,nnxp(ng),nnyp(ng)  &
                   ,scr1(1),srcRead(7)%scr(1,:,:,:))

              !use this for 2 dim:
           else
              call atob(vtab_r(nvars,ng)%npts,  &
                   scr1(1),srcRead(2)%scr(1,:,:,1) )
              
           endif
           
           SELECT CASE (vtab_r(nvars,ng)%idim_type)
           CASE (2)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = 1
               call mk_2_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,1), &
                              vtab_r(nvars,ng)%var_p, &
                              nnxp(ng), nnyp(ng), &
                              m2, m3, ia, iz, ja, jz)
            CASE (3)
               nzl = nnzp(1)
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = 1
               call mk_3_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(:,:,:,1), &
                              vtab_r(nvars,ng)%var_p, &
                              nnzp(ng),nnxp(ng), nnyp(ng), &
                              m1, m2, m3, ia, iz, ja, jz)
            CASE (4)               
               nzl = nzg
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch                              
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr, &
                              vtab_r(nvars,ng)%var_p, &
                              nzg,nnxp(ng), nnyp(ng),npatch, &
                              nzg, m2, m3,npatch, ia, iz, ja, jz)            
            CASE (5)
               nzl = nzs
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr, &
                              vtab_r(nvars,ng)%var_p, &
                              nzs,nnxp(ng), nnyp(ng),npatch, &
                              nzs, m2, m3,npatch, ia, iz, ja, jz)            
            CASE (6)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = npatch
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,:), &
                              vtab_r(nvars,ng)%var_p, &
                              1,nnxp(ng), nnyp(ng),npatch, &
                              1, m2, m3,npatch, ia, iz, ja, jz)
            CASE (7)
               nzl = 1
               nxl = nnxp(1)
               nyl = nnyp(1)
               n4 = nwave
               call mk_4_buff(srcRead(vtab_r(nvars,ng)%idim_type)%scr(1,:,:,:), &
                              vtab_r(nvars,ng)%var_p, &
                              1,nnxp(ng), nnyp(ng),nwave, &
                              1, m2, m3,nwave, ia, iz, ja, jz)
            CASE DEFAULT
               PRINT *, 'Wrong idim_type: ',vtab_r(nvars,ng)%idim_type
               STOP 'history_start'
            END SELECT
             
           deallocate(scr1)
           deallocate(scr2)

        endif

     enddo

  enddo
  !print *, 'end of recycle'
  return
end subroutine recycle

!******************************************************************

subroutine rearrange_aot(nwave, nx, ny, array_i, array_o)
   implicit none
   integer, intent(in) :: nx, ny, nwave
   real, intent(in)    :: array_i(nx,ny,nwave)
   real, intent(out)   :: array_o(nx,ny,nwave)
   
   integer :: i, j, w
   
   do w=1, nwave
      do j=1, ny
         do i=1, nx
	    array_o(i,j,w) = array_i(i,j,w)
	    !if(w==12 .and. array_o(i,j,w) > 0.01) print*,i,j,array_o(i,j,w)
	 enddo
      enddo
   enddo
!       open(19,file='aot.gra',         &
!            form='unformatted',access='direct',status='unknown',  &
!            recl=4*nx*ny)	  
!       nrec=1
!       write(19,rec=nrec) array_o(:,:,12)
!       close (19)
!       stop 333

end subroutine rearrange_aot

!*******************************************************************************

subroutine unarrange_p(n2,n3,n4,n5,a,b)
implicit none

integer :: n2,n3,n4,n5
real :: a(n2,n3,n4,n5),b(n4,n2,n3,n5)

integer :: i,j,k,ip

do ip = 1,n5
   do k = 1,n4
      do j = 1,n3
         do i = 1,n2
            b(k,i,j,ip) = a(i,j,k,ip)
         enddo
      enddo
   enddo
enddo
return
end
