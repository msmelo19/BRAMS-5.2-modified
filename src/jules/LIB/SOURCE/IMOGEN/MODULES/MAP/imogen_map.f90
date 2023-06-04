  MODULE imogen_map
  
    USE imogen_constants, ONLY : n_imogen_land
  
    IMPLICIT NONE

!------------------------------------------------------------------------------
! Module variables
!------------------------------------------------------------------------------
    INTEGER, DIMENSION(n_imogen_land) :: sgindinv

  CONTAINS

    SUBROUTINE get_imogen_map(imogenOrderFile)
!------------------------------------------------------------------------------
! Module imports
!------------------------------------------------------------------------------
      USE time_loc, ONLY :                                        &
!  imported scalars with intent (inout)
        latitude,longitude

      USE ancil_info, ONLY : land_pts,land_index,row_length
    
      USE file_utils, ONLY : closeFile,fileUnit,openFile
    
      USE inout, ONLY : formatAsc
      
      USE imogen_constants, ONLY : N_IMOGEN_LAND

      IMPLICIT NONE
      
      CHARACTER(len=*) ::                                         &
        imogenOrderFile ! Filename to read IMOGEN points order from
  
!------------------------------------------------------------------------------
! Local variable declarations
!------------------------------------------------------------------------------
      INTEGER :: indlat,indlon,i,j,l,map1(96,56),funit
      INTEGER :: sgjind(land_pts),sgind(land_pts)

      funit = fileUnit( formatAsc )
      CALL openFile(1,.FALSE.,funit,'read',formatAsc,imogenOrderFile,'old')
      READ(funit,*) map1
      CALL closeFile(funit,formatAsc)
  
      DO l=1,land_pts
        j = (land_index(l) - 1) / row_length + 1
        i = land_index(l) - (j-1) * row_length
      
        indlat = int((latitude(i,j) + 55.1) / 2.5) + 1
        indlon = int((longitude(i,j) + 180.1) / 3.75) + 1
        sgind(l) = map1(indlon, 57 - indlat)
        sgjind(l) = (indlat - 1) * 96 + indlon
      ENDDO

      sgindinv(:) = 0
      DO i=1,N_IMOGEN_LAND
        DO j=1,land_pts
          IF (i == sgind(j)) sgindinv(i) = j
        ENDDO
      ENDDO

    END SUBROUTINE get_imogen_map

  END MODULE imogen_map
