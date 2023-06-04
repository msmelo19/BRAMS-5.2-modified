!###########################################################################
!  B - Regional Atmospheric Modeling System - RAMS
!###########################################################################


module mem_plume_chem1

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  USE ModNamelistFile, only: namelistFile

  include "i8.h"
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

  type plume_vars   
     real, pointer, dimension(:,:)  :: fct     
  end type plume_vars  
  type (plume_vars)    , allocatable :: plume_g (:,:,:), plumem_g (:,:,:)

  type plume_mean_vars   
     real, pointer, dimension(:,:)  :: fire_size     
     real, pointer, dimension(:,:)  :: flam_frac     
  end type plume_mean_vars  

  type (plume_mean_vars)    , allocatable :: plume_mean_g(:,:), plume_meanm_g(:,:)


  integer, parameter :: nveg_agreg      = 4
  integer, parameter :: tropical_forest = 1
  integer, parameter :: boreal_forest   = 2
  integer, parameter :: savannah        = 3
  integer, parameter :: grassland       = 4 ! must be equal to nveg_agreg
  character(len=20), parameter :: veg_name(nveg_agreg) = (/ &
                               'Tropical-Forest', &
                               'Boreal-Forest  ', &
                               'Savanna        ', &
                               'Grassland      ' /)
  character(len=20), parameter :: spc_suf(nveg_agreg) = (/ &
                               'agtf' , &  ! trop forest
                               'agef' , &  ! extratrop forest
                               'agsv' , &  ! savanna
                               'aggr'   /) ! grassland
  REAL                         :: prfrq 	    
  INTEGER                      :: plumerise 

contains
  !---------------------------------------------------------------

  subroutine alloc_plume_chem1(plume,plume_mean,n1,n2,n3,nspecies)

    use chem1_list, only : spc_alloc,spc_name,src,on
    implicit none

    integer,intent(in) :: n1,n2,n3,nspecies
    integer :: ispc,iv
    
    type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean

    integer:: imean_plume
    imean_plume = 1 !change this at emis_flam_smold routine also
         
    !print*,'----------------------------------------------------------------'
    !print*,' memory allocation for plumerise sources:'
    if(imean_plume /= on) then
    
    	  do ispc=1,nspecies
    	   !print*,'spc=',spc_name(ispc),'size=',n1,n2,n3

    	   if(spc_alloc(src,ispc) == on) then 

    	      do iv=1,nveg_agreg
    		  allocate (plume(iv,ispc)%fct(n2,n3))
    		  plume(iv,ispc)%fct(:,:) = 0.
    	      enddo
    	   endif
          enddo
     else 
 
          do iv=1,nveg_agreg
            allocate (plume_mean(iv)%flam_frac(n2,n3))
            plume_mean(iv)%flam_frac(:,:)=0.
          enddo
    endif


    do iv=1,nveg_agreg
       allocate (plume_mean(iv)%fire_size(n2,n3))
       plume_mean(iv)%fire_size(:,:)=0.
    enddo




    return
  end subroutine alloc_plume_chem1

  !---------------------------------------------------------------

  subroutine dealloc_plume_chem1(plume,plume_mean,nspecies)

   implicit none

   integer,intent(in) :: nspecies
   type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
   type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
   integer :: ispc,iv
  
    !  Deallocate arrays
    do ispc=1,nspecies
     do iv=1,nveg_agreg
       if (associated(plume(iv,ispc)%fct)) deallocate(plume(iv,ispc)%fct)
     enddo
    enddo
     do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) deallocate(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) deallocate(plume_mean(iv)%flam_frac)
     enddo
 end subroutine dealloc_plume_chem1

  !---------------------------------------------------------------

  subroutine nullify_plume_chem1(plume,plume_mean,nspecies)

    implicit none

    integer,intent(in) :: nspecies
    type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume
    type (plume_mean_vars) ,dimension(nveg_agreg         ) :: plume_mean
    
    integer :: ispc,iv
 
    do ispc=1,nspecies

      do iv=1,nveg_agreg
        if (associated(plume(iv,ispc)%fct)) nullify(plume(iv,ispc)%fct)
      enddo
    
    enddo
    do iv=1,nveg_agreg
       if (associated(plume_mean(iv)%fire_size)) nullify(plume_mean(iv)%fire_size)
       if (associated(plume_mean(iv)%flam_frac)) nullify(plume_mean(iv)%flam_frac)
    enddo
 
    return
  end subroutine nullify_plume_chem1

  !---------------------------------------------------------------

  subroutine filltab_plume_chem1(plume,plumem,plume_mean,plume_meanm&
                          ,imean,n1,n2,n3,nspecies,ng)

    use chem1_list, only: spc_alloc,spc_name,src,on 
    use mem_chem1, only: chem1_g

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    use var_tables, only: InsertVTab
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    implicit none

    integer, intent(in) :: imean,n1,n2,n3,nspecies,ng
    type (plume_vars)    ,dimension(nveg_agreg,nspecies) :: plume,plumem
    type (plume_mean_vars) ,dimension(nveg_agreg) ::plume_mean,plume_meanm

    integer :: ispc,iv  

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    integer(kind=i8) :: npts
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!   integer :: npts
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    integer:: imean_plume
    imean_plume = 1 !change this at emis_flam_smold routine also
         

    ! Fill pointers to arrays into variable tables
    ! 2d var
    npts=n2*n3

    if(imean_plume /= on) then

    	do ispc=1,nspecies
    	  if (associated(chem1_g(ispc,ng)%sc_p)) then  ! check this latter

!------- sources 
    	   
    	    if(spc_alloc(src,ispc) == on) then
    		do iv=1,nveg_agreg

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
    		    call InsertVTab(plume(iv,ispc)%fct,       &
                                    plumem(iv,ispc)%fct,      &
    				    ng, npts, imean,          &
                                    trim(spc_name(ispc))//'_'//trim(spc_suf(iv))//&
                                    ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!    		    call vtables2 (  plume(iv,ispc)%fct(1,1)&
!    				  , plumem(iv,ispc)%fct(1,1)&
!    				  ,ng, npts, imean, trim(spc_name(ispc))&
!    				  //'_'//trim(spc_suf(iv))//&
!    				  ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

    		enddo
    	    endif
    	   
    	  endif

    	enddo
!---  flam frac mean
    else
        do iv=1,nveg_agreg

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
           call InsertVTab(plume_mean(iv)%flam_frac,  &
             	           plume_meanm(iv)%flam_frac,  &
                           ng, npts, imean, 'flam_frac'&
                           //'_'//trim(spc_suf(iv))//&
                           ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!           call vtables2 (  plume_mean(iv)%flam_frac(1,1)&
!             	        , plume_meanm(iv)%flam_frac(1,1)&
!   		        ,ng, npts, imean, 'flam_frac'&
!        	        //'_'//trim(spc_suf(iv))//&
!        	        ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------

        enddo
    endif

!--- fire size   
    do iv=1,nveg_agreg

!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
       call InsertVTab(plume_mean(iv)%fire_size,  &
                       plume_meanm(iv)%fire_size, &
                       ng, npts, imean, 'firesize'&
                       //'_'//trim(spc_suf(iv))//&
                       ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-4-OLD)--------------------------------------------------------------------
!       call vtables2 (  plume_mean(iv)%fire_size(1,1)&
!        	     , plume_meanm(iv)%fire_size(1,1)&
!   		     ,ng, npts, imean, 'firesize'&
!        	     //'_'//trim(spc_suf(iv))//&
!        	     ' :2:hist:anal:mpti:mpt3:mpt1')
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------
   
    enddo
  end subroutine filltab_plume_chem1


!--(DMK-CCATT-BRAMS-5.0-INI)------------------------------------------------------------------
  subroutine StoreNamelistFileAtMem_plumeChem1(oneNamelistFile)
    type(namelistFile), pointer :: oneNamelistFile
    plumerise  = oneNamelistFile%plumerise
    prfrq = oneNamelistFile%prfrq
  end subroutine StoreNamelistFileAtMem_plumeChem1
!--(DMK-CCATT-BRAMS-5.0-FIM)------------------------------------------------------------------


  !---------------------------------------------------------------
end module mem_plume_chem1
