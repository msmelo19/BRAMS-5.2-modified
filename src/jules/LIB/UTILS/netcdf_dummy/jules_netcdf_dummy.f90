! "Dummy" netCDF code, to allow JULES to be built without netCDF libraries.
! If a netCDF procedure is called at run time, an error will be raised.

  module netcdf

  implicit none

!-------------------------------------------------------------------------------
! Scalar parameters.
!-------------------------------------------------------------------------------
! The following parameters are given arbitrary values, since they are not used.
! Declared as parameters rather than variables simply to agree with "proper" netCDF.  
  integer, parameter :: nf90_clobber = 0      
  integer, parameter :: nf90_float = 0     
  integer, parameter :: nf90_global = 0   
  integer, parameter :: nf90_noClobber = 0        
  integer, parameter :: nf90_noErr = 0 
  integer, parameter :: nf90_noWrite = 0
  integer, parameter :: nf90_unlimited = 0

!-------------------------------------------------------------------------------
! Interfaces.
!-------------------------------------------------------------------------------
  interface nf90_get_var
  module procedure nf90_get_var_real_1D,nf90_get_var_real_2D,nf90_get_var_real_3D
  end interface nf90_get_var

  interface nf90_put_att
  module procedure nf90_put_att_char,nf90_put_att_real
  end interface nf90_put_att

  interface nf90_put_var
  module procedure nf90_put_var_int_1d,nf90_put_var_real_1d,nf90_put_var_real_2d
  end interface nf90_put_var

  contains

!###############################################################################
!###############################################################################

  function nf90_close( ncID ) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_close

!###############################################################################
!###############################################################################

  function nf90_create( path,cmode,ncID,i1,i2 ) result(iresult)
  implicit none
  integer :: iresult
  character(len=*) :: path
  integer :: cmode,ncID
  integer, optional :: i1,i2
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_create

!###############################################################################
!###############################################################################

! Specific functions for the generic nf90_get_var.
! At present there are only versions for real variables, with 1 to 3 dimensions.

  function nf90_get_var_real_1D( ncID,varID,values,start,count,stride,map ) result(iresult)
! Module procedure with generic name nf90_get_var.
  implicit none
  integer :: iresult
  integer :: ncID,varID
  integer, optional :: start(:),count(:),stride(:),map(:)
  real :: values(:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_get_var_real_1D
!###############################################################################

  function nf90_get_var_real_2D( ncID,varID,values,start,count,stride,map ) result(iresult)
! Module procedure with generic name nf90_get_var.
  implicit none
  integer :: iresult
  integer :: ncID,varID
  integer, optional :: start(:),count(:),stride(:),map(:)
  real :: values(:,:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_get_var_real_2D
!###############################################################################

  function nf90_get_var_real_3D( ncID,varID,values,start,count,stride,map ) result(iresult)
! Module procedure with generic name nf90_get_var.
  implicit none
  integer :: iresult
  integer :: ncID,varID
  integer, optional :: start(:),count(:),stride(:),map(:)
  real :: values(:,:,:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_get_var_real_3D
!###############################################################################
!###############################################################################

  function nf90_def_dim(ncID,name,len,dimID) result(iresult)
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,len
  character(len=*), intent(in) :: name
  integer, intent(out) :: dimID
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_def_dim
!###############################################################################
!###############################################################################

  function nf90_def_var(ncID,name,xtype,dimIDs,varID) result(iresult)
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,xtype
  integer, intent(in), optional :: dimIDs(:)
  character(len=*), intent(in) :: name
  integer, intent(out) :: varID
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_def_var
!###############################################################################
!###############################################################################

  function nf90_enddef(ncID,i1,i2,i3,i4) result(iresult)
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID
  integer, intent(in), optional :: i1,i2,i3,i4
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_enddef
!###############################################################################
!###############################################################################

  function nf90_inquire(ncID) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_inquire
!###############################################################################
!###############################################################################

  function nf90_inquire_dimension(ncID,dimID,name,len) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID,dimID
  integer, optional :: len
  character(len=*),optional :: name
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_inquire_dimension
!###############################################################################

  function nf90_inq_dimid( ncID,name,dimID ) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID,dimID
  character(len=*) :: name
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_inq_dimid
!###############################################################################

  function nf90_inq_varid( ncID,name,varID ) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID,varID
  character(len=*) :: name
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_inq_varid
!###############################################################################

  function nf90_open( filename,mode,ncid,chunkSize ) result(iresult)
  implicit none
  integer :: iresult
  integer :: ncID,mode
  integer, optional :: chunkSize
  character(len=*) :: fileName
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_open
!###############################################################################
!###############################################################################

  function nf90_put_att_char(ncID,varID,name,values) result(iresult)
! Module procedure with generic name nf90_put_att.
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,varID
  character(len=*), intent(in) :: name,values
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_put_att_char

!###############################################################################
!###############################################################################
  function nf90_put_att_real(ncID,varID,name,values) result(iresult)
! Module procedure with generic name nf90_put_att.
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,varID
  character(len=*), intent(in) :: name
  real, intent(in) :: values
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_put_att_real
!###############################################################################
!###############################################################################

  function nf90_put_var_int_1d(ncID,varID,values,start,count,stride,map) result(iresult)
! Module procedure with generic name nf90_put_var
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,varID
  integer, intent(in) :: values(:)
  integer, intent(in), optional :: start(:),count(:),stride(:),map(:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_put_var_int_1d

!###############################################################################
!###############################################################################

  function nf90_put_var_real_1d(ncID,varID,values,start,count,stride,map) result(iresult)
! Module procedure with generic name nf90_put_var
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,varID
  real, intent(in) :: values(:)
  integer, intent(in), optional :: start(:),count(:),stride(:),map(:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_put_var_real_1d

!###############################################################################
!###############################################################################

  function nf90_put_var_real_2d(ncID,varID,values,start,count,stride,map) result(iresult)
! Module procedure with generic name nf90_put_var
  implicit none
  integer :: iresult
  integer, intent(in) :: ncID,varID
  real, intent(in) :: values(:,:)
  integer, intent(in), optional :: start(:),count(:),stride(:),map(:)
  iresult = 0
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_put_var_real_2d

!###############################################################################
!###############################################################################

  function nf90_strerror(ncerr) result(cresult)
  implicit none
  character(len=80) :: cresult
  integer :: ncerr
  cresult = ''
!DSM  write(*,*)'ERROR: attempting to use a netCDF procedure, but netCDF library not available.'
!DSM  write(*,*)'This is a ''dummy'' procedure.'
!DSM  write(*,*)'To use netCDF, remake model with the ''proper'' netCDF library.'
!DSM  stop
  end function nf90_strerror
!###############################################################################

  end module netcdf
