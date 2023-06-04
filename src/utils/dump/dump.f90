MODULE dump
  
   private
  
   public :: dumpvar,map2d,warning
  
   interface dumpvar
     module procedure dumpVar3D_real
     module procedure dumpVar2D_real
     module procedure dumpVar1D_real
     module procedure dumpVar3D_int
     module procedure dumpVar2D_int
     module procedure dumpVar1D_int
   end interface

  CONTAINS

      SUBROUTINE dumpVar3D_real(Var2BDump,filenameIn,fileExtension,xDim,yDim,zDim)
      !Subroutine for dump variables 3D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim,yDim,zDim !Dims of variable
      REAL, INTENT(IN) :: Var2BDump(xDim,yDim,zDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(3) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 3                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i,j,k
      
      REAL :: dset_data(xDim,yDim,zDim), data_out(xDim,yDim,zDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim,yDim,zDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
         DO j = 1, yDim
            DO k=1,zDim
               dset_data(i,j,k) = Var2BDump(i,j,k)
            END DO
         END DO
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      data_dims(2) = yDim
      data_dims(3) = zDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar3D_real
      
      SUBROUTINE dumpVar2D_real(Var2BDump,filenameIn,fileExtension,xDim,yDim)
      !Subroutine for dump variables 2D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim,yDim !Dims of variable
      REAL, INTENT(IN) :: Var2BDump(xDim,yDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 2                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i,j
      
      REAL :: dset_data(xDim,yDim), data_out(xDim,yDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim,yDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
         DO j = 1, yDim
               dset_data(i,j) = Var2BDump(i,j)
         END DO
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      data_dims(2) = yDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar2d_real
      
      SUBROUTINE dumpVar1D_real(Var2BDump,filenameIn,fileExtension,xDim)
      !Subroutine for dump variables 1D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim !Dims of variable
      REAL, INTENT(IN) :: Var2BDump(xDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(1) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 1                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i
      
      REAL :: dset_data(xDim), data_out(xDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
               dset_data(i) = Var2BDump(i)
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar1d_real
      
      SUBROUTINE dumpVar3D_int(Var2BDump,filenameIn,fileExtension,xDim,yDim,zDim)
      !Subroutine for dump variables 3D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim,yDim,zDim !Dims of variable
      INTEGER, INTENT(IN) :: Var2BDump(xDim,yDim,zDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(3) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 3                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i,j,k
      
      REAL :: dset_data(xDim,yDim,zDim), data_out(xDim,yDim,zDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(3) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim,yDim,zDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
         DO j = 1, yDim
            DO k=1,zDim
               dset_data(i,j,k) = Var2BDump(i,j,k)
            END DO
         END DO
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      data_dims(2) = yDim
      data_dims(3) = zDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar3D_int
      
      SUBROUTINE dumpVar2D_int(Var2BDump,filenameIn,fileExtension,xDim,yDim)
      !Subroutine for dump variables 2D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim,yDim !Dims of variable
      INTEGER, INTENT(IN) :: Var2BDump(xDim,yDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(2) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 2                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i,j
      
      REAL :: dset_data(xDim,yDim), data_out(xDim,yDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim,yDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
         DO j = 1, yDim
               dset_data(i,j) = Var2BDump(i,j)
         END DO
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      data_dims(2) = yDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar2d_int
      
      SUBROUTINE dumpVar1D_int(Var2BDump,filenameIn,fileExtension,xDim)
      !Subroutine for dump variables 1D in hdf5 format
      !  The Var2BDump must be real
      !  The size of filename+fileExtension must be less than 253 characters
      ! Author: Luiz Flavio
      USE HDF5 ! This module contains all necessary modules
      
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: xDim !Dims of variable
      INTEGER, INTENT(IN) :: Var2BDump(xDim) !Var to be dumped in a file filename above
      CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
      CHARACTER(LEN=*), INTENT(IN) :: fileExtension !Filenamen extension to be added
      
      CHARACTER(LEN=256) :: filename ! File name
      CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
      
      INTEGER(HID_T) :: file_id       ! File identifier
      INTEGER(HID_T) :: dset_id       ! Dataset identifier
      INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
      
      
      INTEGER(HSIZE_T), DIMENSION(1) :: dims  ! Dataset dimensions
      INTEGER     ::   rank = 1                        ! Dataset rank
      
      INTEGER     ::   error ! Error flag
      INTEGER     ::   i
      
      REAL :: dset_data(xDim), data_out(xDim) ! Data buffers
      INTEGER(HSIZE_T), DIMENSION(1) :: data_dims
      
      filename=trim(filenameIn)//trim(fileExtension)//'.h5'
      dims = (/xDim/)
      !
      ! Initialize the dset_data array.
      !
      DO i = 1, xDim
               dset_data(i) = Var2BDump(i)
      END DO
      
      !
      ! Initialize FORTRAN interface.
      !
      CALL h5open_f(error)
      
      !
      ! Create a new file using default properties.
      !
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      
      !
      ! Create the dataspace.
      !
      CALL h5screate_simple_f(rank, dims, dspace_id, error)
      
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, &
            dset_id, error)
      
         !
      ! Write the dataset.
      !
      data_dims(1) = xDim
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)
      
      !
      !
      ! Close the dataset.
      !
      CALL h5dclose_f(dset_id, error)
      
      !
      ! Close the file.
      !
      CALL h5fclose_f(file_id, error)
      
      !
      ! Close FORTRAN interface.
      !
      CALL h5close_f(error)
      
      END SUBROUTINE dumpvar1d_int
      
      subroutine map2d(var,filenameIn,lev,time,xDim,yDim,is2do,is2del)
            IMPLICIT NONE
      
         INTEGER, INTENT(IN) :: xDim,yDim !Dims of variable
         INTEGER, INTENT(IN) :: lev
         INTEGER,INTENT(IN) :: time
         REAL, INTENT(IN) :: Var(xDim,yDim) !Var to be dumped in a file filename above
         CHARACTER(LEN=*), INTENT(IN) :: filenameIn !Filename to be composed
         LOGICAL, INTENT(IN) :: is2do,is2del
         
         CHARACTER(LEN=3) :: clev
         CHARACTER(LEN=6) :: ctime
         CHARACTER(LEN=3) :: cy
         INTEGER :: i,j
         
         WRITE(clev,fmt='(I3.3)') lev
         WRITE(ctime,fmt='(I6.6)') time
         WRITE(cy,fmt='(I3.3)') yDim
         
         open(88,file=trim(filenameIn)//'.'//ctime//'.'//clev//'.dat')
         
         DO j=1,ydim
            write (88,fmt='('//cy//'(E16.6,1X))') (var(i,j),i=1,xdim)
         END DO
            
         CLOSE(88)
         
         open(88,file=trim(filenameIn)//'.'//ctime//'.'//clev//'.gp')
         write (88,fmt='(A)') 'set terminal png'
         write (88,fmt='(A)') 'set output "'//trim(filenameIn)//'.'//ctime//'.'//clev//'.png"'
         write (88,fmt='(A)') 'set pm3d map'
         write (88,fmt='(A)') 'splot "'//trim(filenameIn)//'.'//ctime//'.'//clev//'.dat'//'" matrix'
         CLOSE (88)
         
         IF(is2do) THEN
            CALL system('gnuplot '//trim(filenameIn)//'.'//ctime//'.'//clev//'.gp')
            IF(is2del) THEN
               CALL system('rm '//trim(filenameIn)//'.'//ctime//'.'//clev//'.dat')
               CALL system('rm '//trim(filenameIn)//'.'//ctime//'.'//clev//'.gp')
            END IF
         END IF
         

      end subroutine map2d
      
      subroutine warning(subpos,text,nlines)
         INTEGER,INTENT(IN) :: nlines
         Character(len=*), intent(in) :: text(nlines),subpos
         INTEGER :: sot,ntimes_bef

         sot=0
         DO n=1,nlines
            IF(len(text(n))>sot) sot=len(text(n))
         END DO
         
         ntimes_bef=sot/2-5
         if(ntimes_bef<0) ntimes_bef=0
         
         WRITE (*,fmt='(A)') ''
         write (*,fmt='(A)') repeat(' ',ntimes_bef)//'       ////'
         write (*,fmt='(A)') repeat(' ',ntimes_bef)//'      (O O)'
         write (*,fmt='(A)') repeat('-',ntimes_bef)//'oOO----(_)-'//repeat('-',ntimes_bef)
         write (*,fmt='(A)') 'From '//subpos//' :'
         DO n=1,nlines
              write (*,fmt='(A)') text(n)
         END DO
         write (*,fmt='(A)') repeat('-',ntimes_bef)//'-------------oOO'//repeat('-',ntimes_bef-5)
         
      END SUBROUTINE warning
      

END MODULE dump