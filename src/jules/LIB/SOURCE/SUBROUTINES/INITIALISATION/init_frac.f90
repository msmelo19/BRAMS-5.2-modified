!###############################################################################
!###############################################################################
! Read fractional coverage for surface types.
! NB A check that sum(frac)=1 is performed in init_ic().
!###############################################################################
!###############################################################################
SUBROUTINE init_frac(nia,nja,npatch,patch_area,veg_fracarea,leaf_class)

  USE ancil_info, ONLY :  &
!  imported scalars with intent(in)
     land_pts  &
!  imported arrays with intent(out)
    ,frac, land_index

  USE file_utils, ONLY:  &
!  imported procedures
     closeFile,fileUnit,findTag,openFile  &
!  imported arrays with intent(out)
    ,irecPrev

  USE inout, ONLY :  &
!  imported scalar parameters
     formatLen,formatAsc,formatBin,formatNc,formatPP,jinUnit,tagAscBin,tagNc  &
!  imported scalars with intent(in)
    ,echo,nxIn,nyIn  &
!  imported scalars with intent(out)
    ,readFracIC  &
!  imported arrays with intent(in)
    ,mapInLand

  USE jules_netcdf, ONLY :  &
!  imported scalars with intent(in)
     ncType

  USE nstypes, ONLY :  &
!  imported scalars with intent(in)
     ntype

  USE readwrite_mod, ONLY :  &
!  imported procedures
     readvar3dComp

  USE switches, ONLY :  &
!  imported scalars with intent(in)
     l_veg_compete,routeOnly
!-------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER ::  &!  scalar parameters
   stashFrac = 216  !  STASH code for fractional coverage field. PP code = 1391.

  INTEGER ::  &!  local SCALARS
    fieldNum     &!  field number in file of first field in frac
   ,i            &!  loop counter
   ,inUnit       &!  unit used to connect to input file
   ,nfieldFile   &!  number of fields per time in a file
   ,nheaderField &!  number of header records before each field in file
   ,nheaderFile  &!  number of header records at start of file
   ,nheaderT     &!  number of headers at start of each time
   ,nlineField   &!  work
   ,readT        &!  time level to be read from file
   ,useIndex      !  index in irecPrev

  LOGICAL ::  &!  local SCALARS
    readFile   !  flag indicating if another file is to be read

  CHARACTER(len=formatLen) ::  &! local SCALARS
    fileFormat   !  format of file

  CHARACTER(len=150) ::  &! local SCALARS
    fileName  &!  the name of a file
   ,varName    !  name of variable in input file

   INTEGER, INTENT(IN)  :: nia,nja,npatch
   REAL,    INTENT(IN)  :: patch_area(nia,nja,npatch),veg_fracarea(nia,nja,npatch),leaf_class(nia,nja,npatch)
!------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! If data are not needed, nothing to do.
!-------------------------------------------------------------------------------
  IF ( routeOnly ) RETURN

  !WRITE(*,"(50('-'),/,a)") 'init_frac'

!-------------------------------------------------------------------------------
! Locate the start of this section in input file and establish whether
! fractional cover is to be read here or with initial condition.
!-------------------------------------------------------------------------------
  CALL findTag( jinUnit,'init_frac','>INIT_FRAC' )
  READ(jinUnit,*) readFracIC

!-------------------------------------------------------------------------------
! For runs with competing vegetation, insist that frac is read as initial condition.
!-------------------------------------------------------------------------------
  IF ( l_veg_compete .AND. .NOT.readFracIC ) THEN
    WRITE(*,*)'ERROR: init_frac: A run with competing vegetation must read frac with initial condition. Set readFracIC=T.'
    STOP
  ENDIF

  IF ( readFracIC ) THEN
    WRITE(*,*)'readFracIC=T, meaning fractional cover will be read with the initial condition.'
    WRITE(*,*)'frac will NOT be read by init_frac'
    RETURN
  ENDIF
!-------------------------------------------------------------------------------
! Establish where data are to be read from, and other details of format.  Only
! read parameters for the file format indicated.
!-------------------------------------------------------------------------------
  READ(jinUnit,*) readFile

  IF ( readFile ) THEN
!-------------------------------------------------------------------------------
!   An external file will be read.
!-------------------------------------------------------------------------------
    READ(jinUnit,*) fileFormat
    READ(jinUnit,*) fileName

    SELECT CASE ( fileFormat )

      CASE ( formatAsc,formatBin,formatPP )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_frac',tagAscBin )
        READ(jinUnit,*) nheaderFile,nheaderField
        READ(jinUnit,*) fieldNum

      CASE ( formatNc )
!       Locate the information in run control file.
        CALL findTag( jinUnit,'init_grid',tagNc,preInit=.TRUE. )
        READ(jinUnit,*) varName

      CASE default
        WRITE(*,*)'ERROR: init_frac: no code for fileFormat=',TRIM(fileFormat)
        STOP

    END SELECT

  ELSE   !  NOT readFile
!-------------------------------------------------------------------------------
!   Data will be read from run control file.
!   The first field will be read, no headers expected. Field number is
!   redundant for stdIn, but is used to set nfieldFile.
!-------------------------------------------------------------------------------
    fileFormat = formatAsc
    nheaderFile = 0
    nheaderField = 0
    fieldNum = 1

  ENDIF
!-------------------------------------------------------------------------------
! Open file.
!-------------------------------------------------------------------------------
  IF ( readFile ) THEN
    inUnit = fileUnit( fileFormat )  !  get unit
!   Use first arg to openFile to set recl (unformatted file) to be enough for a single value.
    !DSM CALL openFile( 1,.FALSE.,inUnit,'read',fileFormat,fileName,'old','init_frac',ncType )
  ELSE
    WRITE(*,*)'Reading fractional cover from the run control file.'
    WRITE(*,*)'Data must be the first field encountered.'
    inUnit = jinUnit
!   Locate the start of the data in the run control file.
    CALL findTag( inUnit,'init_frac','>DATA',preInit=.TRUE. )
  ENDIF

!-------------------------------------------------------------------------------
! Read data.
!-------------------------------------------------------------------------------
  IF ( inUnit==jinUnit .AND. nxIn*nyIn==1 ) THEN
!   If only need to read one space point from run control file, expect no new
!   line between fields (eg all on one line).  Calling readVar means we could
!   cope with headers in the run control file!  But there's no need since
!   annotation is already simple in this case.
    READ(jinUnit,*) frac(:,:)
  ELSE
!   Simplifying assumptions regarding input file. Otherwise have to read these in.
    readT      = 1                !  time level to read from file
    nfieldFile = fieldNum+ntype-1 !  # of fields in file. Set to last level of required field - OK while readT=1
    nheaderT   = 0                !  no headers at top of each time
    nlineField = 0                !  0 means will not attempt to read ASCII line-by-line

!   Set index to use for irecPrev with netCDF files - this irecPrev isn't changed,
!   but need to keep index within bounds.
    useIndex = inUnit
    IF ( fileFormat == formatNC ) useIndex = 1

    !DSM CALL readVar3dComp( readT,fieldNum,stashFrac,.FALSE.,irecPrev(useIndex)  &
    !DSM                    ,nfieldFile,nheaderFile  &
    !DSM                    ,nheaderT,nheaderField,nlineField,nxIn,nyIn,ntype,.FALSE. &
    !DSM                    ,inUnit,varName  &
    !DSM                    ,mapInLand(:),(/(i,i=1,land_pts)/)  &
    !DSM                    ,fileFormat,frac(:,:),'init_frac','init_frac',ncType )
  ENDIF

   !{DSM
      CALL frac_from_leaf( frac,land_pts,ntype,land_index,nia,nja,npatch,patch_area,veg_fracarea    &
                           ,leaf_class )
   !DSM}

!-------------------------------------------------------------------------------
! Close file if it is not the JULES control file
!-------------------------------------------------------------------------------
  IF ( inUnit /= jinUnit ) CALL closeFile( inUnit,fileFormat )

!------------------------------------------------------------------------------
! If requested, write to stdout the fractions for each gridbox.
!-------------------------------------------------------------------------------
  IF ( echo ) THEN
    DO i=1,ntype
      WRITE(*,"(a,i3,a,(10f6.3))") 'type ',i,' frac=',frac(:,i)
    ENDDO
  ENDIF

END SUBROUTINE init_frac

!{DSM
!------------------------------------------------------------------------------
!--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
!-------------------------------------------------------------------------------
SUBROUTINE brams2jules(veg,ntype)
   IMPLICIT NONE
   INTEGER, INTENT(in)              :: ntype
   CHARACTER (len=80), INTENT(out)  :: veg(ntype)
   
   IF (ntype /= 9 ) THEN
      PRINT*, "ntype=",ntype
      PRINT*, "ATENCAO... ntype <> 9, Deve-se ajustar a subrotina brams2jules em init_frac2brams.f90"
      STOP
   ENDIF
   veg(1)='06 07 20'       !tJ=1 => BT=broadleaf trees
   veg(2)='04 05 14'       !tJ=2 => NT=needleleaf trees
   veg(3)='15'             !tJ=3 => C3G=C3 (temperate) grass
   veg(4)='08 09'          !tJ=4 => C4G=C4 (tropical) grass
   veg(5)='11 12 13 18'    !tJ=5 => shrub
   veg(6)='19 21'          !tJ=6 => urban
   veg(7)='00 01 16 17'    !tJ=7 => lake=inland water
   veg(8)='03 10'          !tJ=8 => soil=bare soil
   veg(9)='02'             !tJ=9 => ice

END SUBROUTINE brams2jules



!--------------------------------------------------------------------------------------------------
!--- Encontra a fracao de vegetacao a partir do mapa da fracao de vegetacao definida pelo leaf3 ---
!--------------------------------------------------------------------------------------------------
SUBROUTINE frac_from_leaf( frac,land_pts,ntype,land_index,nia,nja,npatch,patch_area,veg_fracarea    &
                           ,leaf_class )
  USE ancil_info, ONLY :  &
!  imported scalars with intent(in) 
     row_length

   IMPLICIT NONE
   INTEGER, INTENT(IN)              :: land_pts,ntype,nia,nja,npatch
   INTEGER, INTENT(IN)              :: land_index(land_pts)
   REAL,    INTENT(OUT)             :: frac(land_pts,ntype)

   REAL,    INTENT(IN)  :: patch_area(nia,nja,npatch),veg_fracarea(nia,nja,npatch),leaf_class(nia,nja,npatch)
 
   INTEGER              :: j,i,l,tJ,n
   CHARACTER (LEN=80)   :: veg(ntype)
   CHARACTER (LEN=2)    :: tB_str          
   
   IF (ntype/=9) THEN
      print*, 'ntype=',ntype
      print*, 'ERRO!!! O modelo estah preparado apenas para 9 tipos de coberturas de solo'
      stop
   ENDIF
   
   !--- Convertendo o tipo de vegetacao do BRAMS para o JULES ---
   CALL brams2jules(veg,ntype)
  
   frac=0.
                      
   DO l=1,land_pts
      j = ( land_index(l)-1 ) / row_length + 1
      i = land_index(l) - ( j-1 ) * row_length

      IF (i > nia .or. j > nja) THEN
         PRINT*, 'ERRO!!! i > nia ou j > nja - i, nia, j, nja =',i, nia, j, nja
         STOP
      ENDIF

      DO n=2,npatch
         WRITE(tB_str,'(i2.2)') nint(leaf_class(i,j,n))
         
         !--- Encontrando o tipo correspondente ao JULES ---
         DO tJ=1,ntype+1    !--- o +1 eh apenas para checar se foi encontrado (condicao abaixo)
            IF (index(veg(tJ),tB_str)/=0) exit
         ENDDO
         
         !--- Checando se encontrou um indice valido para a vegetacao do JULES ---
         IF (tJ>ntype) THEN
            PRINT*, 'ERRO!!! Nao foi encontrado uma correspondencia entre BRAMS e JULES'
            STOP
         ENDIF
         
         !30Mai/2013 frac(l,tJ)=frac(l,tJ) + max(0.,patch_area(i,j,n)*veg_fracarea(i,j,n))
         frac(l,tJ)=frac(l,tJ) + max(0.,patch_area(i,j,n))

      ENDDO  !--- DO n=1,npatch
         
      !--- veg_fracarea(:,:,1) eh agua, assim deve-se incluir o patch_area(:,:,1) em tJ=7
      frac(l,7)=frac(l,7)+max(0.,patch_area(i,j,1))
      
      WHERE ( frac(l,:) <= 1.00004E-06 ) frac(l,:) = 1.00004E-06  ! valor minimo que o JULES trabalha

      frac(:,9)=0.0 !--- Mas ice fraction deve iniciar com zero.
      
      n=1
      DO WHILE (sum(frac(l,:)) > 1.0)
         IF (frac(l,n)>0.01) frac(l,n)=frac(l,n)-0.01                             !delta=0.01
         n=n+1
         if (n>8) n=1     
      ENDDO
      
      !--- Para garantir que a fracao total seja igual a 1 ---
      frac(l,1)=1.-( frac(l,2)+frac(l,3)+frac(l,4)     &
                     +frac(l,5)+frac(l,6)+frac(l,7)+frac(l,8)+frac(l,9) )
                     
      
   ENDDO  !l=1,land_pts
  
END SUBROUTINE frac_from_leaf

!DSM}
