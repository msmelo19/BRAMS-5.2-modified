 module mem_tuv
   USE ModTuv, ONLY: nw

  implicit none
  
  type mtuv
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_tauaer
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_wol
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_gol
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_taucld 
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_wcld
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: g_gcld
  end type
  
  type (mtuv), allocatable, dimension(:) :: carma_tuv
  
  INTEGER,DIMENSION(:),allocatable :: tuv2carma  

 contains

  subroutine alloc_carma_tuv()
    use mem_globrad, only : ntotal, nlayer

!    use mem_grid, only : ngrids,nnxp,nnyp
    use mem_grid, only: ngrids
    use node_mod, only: nodemxp,  nodemyp, mynum

    integer i
  
    integer :: ierr
    character(len=8) :: c0
    character(len=*), parameter :: h="**(alloc_carma_tuv)**"
    !print *,'LFR->test, nlayer: ',ntotal,nlayer; CALL flush(6)

    allocate (carma_tuv(ngrids),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(ngrids) fails with stat="//trim(adjustl(c0)))
     endif
    
    do i=1,ngrids
!     print*,'i=',i,ntotal,nlayer,nnxp(i),nnyp(i); call flush(6)
    
!     allocate (carma_tuv(i)%g_tauaer(ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)    
     allocate (carma_tuv(i)%g_tauaer(ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)    
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_tauaer fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv(i)%g_wol   (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv(i)%g_wol   (ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_wol fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv(i)%g_gol   (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv(i)%g_gol   (ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_gol fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv(i)%g_taucld(ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv(i)%g_taucld(ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_taucld fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv(i)%g_wcld  (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv(i)%g_wcld  (ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_wcld fails with stat="//trim(adjustl(c0)))
     endif

!     allocate (carma_tuv(i)%g_gcld  (ntotal,nlayer,nnxp(i),nnyp(i)),stat=ierr)
     allocate (carma_tuv(i)%g_gcld  (ntotal,nlayer,nodemxp(mynum,i),nodemyp(mynum,i)),stat=ierr)
     if (ierr/=0) then
         write(c0,"(i8)") ierr
         call fatal_error(h//"Allocating carma_tuv(i)%g_gcld fails with stat="//trim(adjustl(c0)))
     endif

!--(BRAMS-5.0-INI)---------------------------------------------------  
     carma_tuv(i)%g_tauaer = 0.
     carma_tuv(i)%g_wol    = 0.
     carma_tuv(i)%g_gol    = 0.
     carma_tuv(i)%g_taucld = 0.
     carma_tuv(i)%g_wcld   = 0.
     carma_tuv(i)%g_gcld   = 0.
!--(BRAMS-5.0-FIM)---------------------------------------------------  

    enddo  
 
    allocate(tuv2carma(nw),stat=ierr)
    if (ierr/=0) then
        write(c0,"(i8)") ierr
        call fatal_error(h//"Allocating tuv2carma fails with stat="//trim(adjustl(c0)))
    endif

!--(BRAMS-5.0-INI)---------------------------------------------------  
    tuv2carma = 0
!--(BRAMS-5.0-FIM)---------------------------------------------------  
    
    return
  end subroutine alloc_carma_tuv



end module mem_tuv
