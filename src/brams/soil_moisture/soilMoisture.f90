!========================================================================
! Soil Moisture Estimate for NWP Models
! Coded and implemented by Rodrigo Gevaerd and Saulo Freitas
! Ref.: Gevaerd, R. e S. R. Freitas, Estimativa operacional da umidade 
! do solo para iniciacao de modelos de previsao numerica da atmosfera. 
! Parte I: Descricao da metodologia e validacao. Rev. Bras. Meteo.,
! volume especial do LBA, 2007.
!========================================================================
MODULE soilMoisture

  use ModNamelistFile, only: namelistFile

  use memSoilMoisture, only: &
       soil_moist, soil_moist_fail, usdata_in, usmodel_in ! INTENT(IN)

  IMPLICIT NONE

  PUBLIC  :: soilMoistureInit
  PUBLIC  :: StoreNamelistFileAtSoilMoisture
!!$  PRIVATE :: gatherData
!!$  PRIVATE :: gatherData2D
!!$  PRIVATE :: gatherData4D
!  PUBLIC :: apiPrlatlon
!  PUBLIC :: interpolacao
!  PRIVATE :: changeDay

!!$  interface gatherData
!!$     module procedure gatherData2D, gatherData4D
!!$  end interface

contains

  SUBROUTINE soilMoistureInit(n1, n2, n3, mzg, mzs, npat, ifm, &
       theta, pi0, pp,                                           &
       soil_water, soil_energy, soil_text,                       &
       glat, glon,                                               &
       lpw_R                                                       &
       )

    USE mem_grid, ONLY : RUNTYPE,  &          ! INTENT(IN)
         iyear1, imonth1, idate1, itime1, &   ! INTENT(IN)
         nnxp, nnyp, nnzp,                &   ! INTENT(IN)
         GlobalSizes                          ! Subroutine

    USE io_params, ONLY : timstr     ! INTENT(IN)

    USE rconstants, ONLY : cpi       ! INTENT(IN)

    USE leaf_coms, ONLY : soilcp,  & ! INTENT(IN)
         slmsts,                   & ! INTENT(IN)
         slcpd                       ! INTENT(IN)

    USE mem_leaf, ONLY : stgoff,   & ! INTENT(IN)
         slmstr,                   & ! INTENT(IN)
         slz                         ! INTENT(IN)

    USE node_mod, ONLY: &
         nodei0, nodej0, & ! INTENT(IN)
         nodemxp, nodemyp, & ! INTENT(IN)
         nmachs,  & ! INTENT(IN)
         mynum,  &  ! INTENT(IN)
         mchnum, &  ! INTENT(IN)
         master_num ! INTENT(IN)

    USE ParLib, ONLY: &
         parf_bcast ! Subroutine

    USE mem_aerad, ONLY: &
         nwave ! INTENT(IN)

    USE ReadBcst, ONLY: &
         gatherData
!!$         PreProcAndGather, & ! Subroutine
!!$         LocalSizesAndDisp   ! Subroutine

    IMPLICIT NONE
    INCLUDE "i8.h"
    ! Arguments:
    INTEGER, INTENT(IN) :: n1, n2, n3, mzg, mzs, npat, ifm
!!$    REAL, INTENT(IN)    :: theta(n1,n2,n3)
!!$    REAL, INTENT(IN)    :: pi0(n1,n2,n3)
!!$    REAL, INTENT(IN)    :: pp(n1,n2,n3)
!!$    REAL, INTENT(IN)    :: glat(n2,n3)
!!$    REAL, INTENT(IN)    :: glon(n2,n3)
!!$    INTEGER, INTENT(IN) :: lpw(n2,n3)
!!$    REAL, INTENT(INOUT) :: soil_water(mzg,n2,n3,npat)
!!$    REAL, INTENT(INOUT) :: soil_energy(mzg,n2,n3,npat)
!!$    REAL, INTENT(IN)    :: soil_text(mzg,n2,n3,npat)
    REAL, INTENT(IN)    :: theta(:,:,:)         !(n1,n2,n3)
    REAL, INTENT(IN)    :: pi0(:,:,:)           !(n1,n2,n3)
    REAL, INTENT(IN)    :: pp(:,:,:)            !(n1,n2,n3)
    REAL, INTENT(IN)    :: glat(:,:)            !(n2,n3)
    REAL, INTENT(IN)    :: glon(:,:)            !(n2,n3)
    REAL, INTENT(IN) :: lpw_R(:,:)             !(n2,n3)
    REAL, INTENT(INOUT) :: soil_water(:,:,:,:)  !(mzg,n2,n3,npat)
    REAL, INTENT(INOUT) :: soil_energy(:,:,:,:) !(mzg,n2,n3,npat)
    REAL, INTENT(IN)    :: soil_text(:,:,:,:)   !(mzg,n2,n3,npat)

    ! Local Variables:

    include "files.h"
    
    INTEGER :: lpw(n2,n3)             !(n2,n3)
    INTEGER            :: i, j, k, ipat, nveg, nsoil
    REAL               :: c1, airtemp, pis
    REAL               :: globalSoilWater(mzg, nnxp(ifm), nnyp(ifm), npat)
    REAL               :: globalSoilText(mzg, nnxp(ifm), nnyp(ifm), npat)
    REAL               :: globalGlon(nnxp(ifm), nnyp(ifm))
    REAL               :: globalGlat(nnxp(ifm), nnyp(ifm))
    INTEGER            :: qi1, qi2, qj1, qj2, ncount,               &
         ii, jj, jc, ic, i1, j1, i2, j2, kk, ifname, k2, &
         ipref, ipref_start, icihourmin
    INTEGER            :: n4us
    INTEGER            :: nlat, nlon
    REAL, ALLOCATABLE  :: slz_us(:)  
    REAL               :: latni, latnf, lonni, lonnf, ilatn, ilonn, &
         ilats, ilons, latn, lonn, lats, lons, dlatr, dlonr
    LOGICAL            :: there, theref
    CHARACTER(len=f_name_length) :: usdata, usmodel
    CHARACTER(len=50)  :: pref
    CHARACTER(len=2)   :: cidate, cimon
    CHARACTER(len=1)   :: cgrid
    CHARACTER(len=4)   :: ciyear
    CHARACTER(len=4)   :: cihourmin
    REAL, ALLOCATABLE  :: api_us(:,:,:), prlat(:,:), prlon(:,:), usdum(:),api_temp(:,:,:)
    REAL               :: DIF_TIME, DUMMY
    INTEGER            :: INT_DIF_TIME, idate2, imonth2, iyear2, hourmin, &
         DA, SAIR
    INTEGER            :: ierr
!!$    INTEGER, PARAMETER :: idim_type_min = 2
!!$    INTEGER, PARAMETER :: idim_type_max = 7
    INTEGER, PARAMETER :: idim_type     = 4
!!$    INTEGER            :: il1(nmachs)
!!$    INTEGER            :: ir2(nmachs)
!!$    INTEGER            :: jb1(nmachs)
!!$    INTEGER            :: jt2(nmachs)
!!$    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: maxLocalSize
!!$    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxSizeGathered
!!$    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxsizeFullField
!!$    INTEGER            :: globalSize(idim_type_min:idim_type_max)
!!$    REAL, ALLOCATABLE  :: localChunk(:)
!!$    REAL, ALLOCATABLE  :: gathered(:)
!!$    REAL, ALLOCATABLE  :: fullField(:)
    CHARACTER(len=16)  :: varn
    INTEGER            :: i0, j0, ia, iz, ja, jz,irec,nrec
    LOGICAL            :: general
    NAMELIST /gradeumso/ latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon

    lpw=int(lpw_R)
    ! Initial grid definitions for grid ifm
    i0 = nodei0(mynum,ifm)
    j0 = nodej0(mynum,ifm)
    ia = nodei0(mynum,ifm) + 1
    iz = nodei0(mynum,ifm) + nodemxp(mynum,ifm)
    ja = nodej0(mynum,ifm) + 1
    jz = nodej0(mynum,ifm) + nodemyp(mynum,ifm)

    iyear2  = iyear1
    imonth2 = imonth1
    idate2  = idate1

    ! Determinacao do tipo de produto de umidade
    DO i=256,1,-1
       IF (usdata_in(I:I)=='/') THEN
          ipref_start = i + 1
          EXIT
       ENDIF

       !   print*,'i=',i
       !   print*,usdata(i+1:ifname)
       !   stop
       !  endif
    ENDDO


    ! DEFINICAO DA ESPESSURA DAS CAMADAS
    ipref = LEN_TRIM(usdata_in)
    pref = usdata_in(ipref_start:ipref)
    !print*,'soil data=', pref
    
    IF (pref=='SM_v2.') THEN
       n4us = 6                    ! modelo V2 com 6 camadas
       ALLOCATE(slz_us(0:n4us))
       slz_us = (/-3.0, -2.0, -1.0, -0.5, -0.25, -0.1, 0. /)    
    ELSEIF (pref=='GL_SM.GPCP.' .OR.  pref=='GL_SM.GPNR.' ) THEN
       n4us = 8                    ! modelo GLSM V2 com 8 camadas
       ALLOCATE(slz_us(0:n4us))
       slz_us = (/-4.5, -2.5, -1.75, -1.0, -0.5, -0.25, -0.13, -0.05, 0./)    
    ELSEIF (pref=='gfs025gr.soilm.') THEN
       n4us = 4                    ! modelo do GFS 4 camadas 
       ALLOCATE(slz_us(0:n4us))
       slz_us = (/-2.0, -1.0,  -0.40, -0.10, 0./) 
    ELSE
       n4us = 4                    ! modelo original com 4 camadas
       ALLOCATE(slz_us(n4us))
       slz_us = (/-2.4, -0.4, -0.1, 0./)  
    ENDIF

    ! COMPOSICAO DO NOME DOS ARQUIVOS DE ENTRADA E SAIDA

    IF ( (runtype(1:7)=='HISTORY') .AND.                 &
         ((SOIL_MOIST=='h') .OR. (SOIL_MOIST=='H') .OR.  &
         (SOIL_MOIST=='a') .OR. (SOIL_MOIST=='A'))) THEN

       dif_time = timstr

       !IF (TIMEUNIT == 'h') THEN
       !INT_DIF_TIME = FLOOR(DIF_TIME/24.)
       !ELSEIF (TIMEUNIT == 'm') THEN
       !INT_DIF_TIME = FLOOR(DIF_TIME/1440.)
       !ELSEIF (TIMEUNIT == 's') THEN
       INT_DIF_TIME = FLOOR(DIF_TIME/86400.)
       !ENDIF

       CALL changeDay(idate1, imonth1, iyear1, INT_DIF_TIME,   &
            idate2, imonth2, iyear2)

    ELSE

       INT_DIF_TIME = 0

    ENDIF

    IF ((SOIL_MOIST_FAIL=='l') .OR. (SOIL_MOIST_FAIL=='L')) THEN
       ! Looking for until 5 days old file
       DA = 5
    ELSE
       DA = 1
    ENDIF

    SAIR = 0

    DO I=1,DA

       WRITE (cidate, '(I2.2)') idate2
       WRITE (cimon,  '(I2.2)') imonth2
       WRITE (ciyear, '(I4)')   iyear2
       WRITE (cgrid,  '(I1)')   ifm

       ! Calculating the hour of simulation
       IF ((ITIME1>=0000) .AND. (ITIME1<1200)) THEN
          hourmin = 0000
          IF (pref=='GL_SM.GPCP.' .OR. pref=='GL_SM.GPNR.' .OR. pref=='gfs025gr.soilm.') hourmin = 00
       ELSE
          hourmin = 1200
          IF (pref=='GL_SM.GPCP.' .OR. pref=='GL_SM.GPNR.' .OR. pref=='gfs025gr.soilm.') hourmin = 12
       ENDIF

       IF (pref=='GL_SM.GPCP.' .OR. pref=='GL_SM.GPNR.' .OR. pref=='gfs025gr.soilm.') THEN
          WRITE (cihourmin, '(I2.2)') hourmin
          icihourmin = 2
       ELSE
          WRITE (cihourmin, '(I4.4)') hourmin
          icihourmin = 4
       ENDIF
       !ipref=INDEX(usdata_in,' ')
       !pref=usdata_in(ipref-2:ipref-1)
       !ipref = len_trim(usdata_in)

       !     if (usdata_in(ipref_start:ipref) == 'us'.OR.usdata_in(ipref_start:ipref) == 'SM') then
       !        pref = usdata_in(ipref_start:ipref)
       !     endif
       !     if (usdata_in(ipref_start:ipref) == 'SM_v2.' .or. usdata_in(ipref_start:ipref) == 'GL_SM.GPCP.') then
       !        pref = usdata_in(ipref_start:ipref)
       !     endif

       IF (pref=='us') THEN
          cihourmin  = ''
          icihourmin = 0
       ENDIF

       usdata = TRIM(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.vfm'

       ifname = LEN_TRIM(usdata)

       IF (mchnum==master_num) INQUIRE(file=usdata(1:ifname),exist=theref)
       CALL parf_bcast(theref, master_num)

       IF (.NOT.theref) &
            usdata = TRIM(usdata_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'.gra'

       usmodel = TRIM(usmodel_in)//ciyear//cimon//cidate//cihourmin(1:icihourmin)//'_g'//cgrid//'.mod'

       IF (mchnum==master_num) THEN
          PRINT *, 'Looking up Soil Moisture Date Files: '
          PRINT *, '  usdata : ', usdata(1:LEN_TRIM(usdata))
          PRINT *, '  usmodel: ', usmodel(1:LEN_TRIM(usmodel))

          ifname = LEN_TRIM(usmodel)
          INQUIRE(file=usmodel(1:ifname), exist=there)
          IF(.NOT.there) THEN

             ifname = LEN_TRIM(usdata)
             INQUIRE(file=usdata(1:ifname), exist=there)
             IF (there) THEN
                SAIR = 1
             ELSE
                PRINT *, 'Files:', usdata(1:ifname), &
                     ' and: ', usmodel(1:ifname), ' not Found!'
             ENDIF
          ELSE
             SAIR = 1
          ENDIF
       ENDIF !(mchnum==master_num)

       CALL parf_bcast(SAIR, master_num)

       IF (SAIR==1) EXIT

       CALL changeDay(idate1, imonth1, iyear1, (INT_DIF_TIME - I),   &
            idate2, imonth2, iyear2)

    ENDDO

    ifname = LEN_TRIM(usmodel)

    IF (mchnum==master_num) INQUIRE(file=usmodel(1:ifname), exist=there)
    CALL parf_bcast(there, master_num)

!- 8/9/2015 srf
!    because the the soil moisture dataset interpolated to the model grid
!    is never checked again, we will allways require the interpolation to
!    be done every model initialization. THis will prevents soil moisture 
!    interpolated for given grid specification being erroneously be used 
!    for a different model configuration
!
!    IF (.NOT.there) THEN
!-srf end
       ifname = LEN_TRIM(usdata)
       IF (mchnum==master_num) INQUIRE(file=usdata(1:ifname), exist=there)
       CALL parf_bcast(there, master_num)

       IF (.NOT.there) THEN
          IF (mchnum==master_num) THEN
             PRINT *, '  usdata : ', usdata(1:ifname)
             PRINT *, '  usmodel: ', usmodel(1:LEN_TRIM(usmodel))
          END IF
          IF ((SOIL_MOIST_FAIL=='s') .OR. (SOIL_MOIST_FAIL=='S')) THEN
             CALL fatal_error('* Heterogenous Soil Moisture Init. ERROR! '//&
                  '  Program will Stop!')
             !!STOP
          ELSE
             IF (mchnum==master_num) THEN
                PRINT *, '  Homogeneous Soil Moisture Initialization.'
             END IF
             RETURN
          ENDIF
       ENDIF

       !srf-2015
       !IF (mchnum==master_num) THEN
          !PRINT*,'----------------------------------------------'
          !PRINT*,'  Homogenous Soil Moisture initialization in'
          !PRINT*,'     point outside the South America'
          !PRINT*,'----------------------------------------------'
       !END IF

       c1 = 0.5*cpi

       DO j=1,n3
          DO i=1,n2

             k2  = lpw(i,j)
             pis = c1*(pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
             airtemp = theta(k2,i,j)*pis

             DO ipat=2,npat
                DO k=1,mzg
                   nsoil = NINT(soil_text(k,i,j,ipat))
                   soil_water(k,i,j,ipat) = &
                        MAX(soilcp(nsoil), slmstr(k)*slmsts(nsoil))
                   soil_energy(k,i,j,ipat) = &
                        (airtemp - 273.15 + stgoff(k))*  &
                        (slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6) + &
                        soil_water(k,i,j,ipat)*3.34e8
                ENDDO
             ENDDO
          ENDDO
       ENDDO

       !--------------------------------------------------------- HETEROGENEOUS
       IF (mchnum==master_num) THEN
          PRINT*,'----------------------------------------------'
          PRINT*,'       Soil Moisture Initialization'
          PRINT*,'----------------------------------------------'
       ENDIF

       ! DADOS DA GRADE DE PRECIPITACAO

       IF (mchnum==master_num) &
            INQUIRE(file=TRIM(usdata_in)//'_ENT', exist=general)
       CALL parf_bcast(general, master_num)

       IF (general) THEN
          IF (mchnum==master_num) THEN
             OPEN (unit=93,file=usdata_in(1:len_trim(usdata_in))//'_ENT',status='old')
             READ (unit=93,nml=gradeumso) !latni, latnf, lonni, lonnf, ilatn, ilonn, nlat, nlon
             CLOSE (unit=93, status='keep')
          ENDIF
          CALL parf_bcast(latni, master_num)
          CALL parf_bcast(latnf, master_num)
          CALL parf_bcast(lonni, master_num)
          CALL parf_bcast(lonnf, master_num)
          CALL parf_bcast(ilatn, master_num)
          CALL parf_bcast(ilonn, master_num)
          CALL parf_bcast(nlat,  master_num)
          CALL parf_bcast(nlon,  master_num)

       ELSEIF (pref=='us') THEN
          latni =  -45.
          latnf =   12.9477
          lonni =  -82.
          lonnf =  -30.055
          ilatn =    0.0359477
          ilonn =    0.0382513
          nlat  = 1613
          nlon  = 1359
       ELSEIF (pref=='SM') THEN
          latni =  -50.125
          latnf =   40.125
          lonni = -120.125
          lonnf =   60.125
          ilatn =    0.250
          ilonn =    0.250
          nlat  =  362
          nlon  =  722
       ELSEIF (pref=='SM_v2.') THEN
          latni =  -59.875              !SM v2  TRMM global
          latnf =   59.875
          lonni = -179.875
          lonnf =  179.875
          ilatn =    0.250
          ilonn =    0.250
          nlat  =  480
          nlon  = 1440
       ELSEIF (pref=='GL_SM.GPCP.') THEN
          latni =  -89.5            !  GPCP global
          latnf =   89.5
          lonni = -179.5
          lonnf =  179.5
          ilatn =    1.
          ilonn =    1.
          nlat  =  180
          nlon  =  360
       ELSEIF (pref=='GL_SM.GPNR.') THEN
          latni =  -89.875             !  TRMM/NAVY + GPCP global
          latnf =   89.875
          lonni = -179.875
          lonnf =  179.875
          ilatn =    0.250
          ilonn =    0.250
          nlat  =  720
          nlon  = 1440
       ELSEIF (pref=='gfs025gr.soilm.') THEN
          latni =  -89.875             !  GFS
          latnf =   89.875
          lonni =  -179.875 ! 0.125
          lonnf =   179.875 !359.875
          ilatn =    0.250
          ilonn =    0.250
          nlat  =  721
          nlon  = 1441
       ELSE
          !!print *, 'Unexpected prefix (',pref,') for soil moisture'
          !!stop 'Program will STOP!'
          CALL fatal_error('Unexpected prefix for soil moisture')
       ENDIF

       ALLOCATE(prlat(nlon,nlat), STAT=ierr)
       IF (ierr/=0) CALL fatal_error("ERROR allocating prlats (soilMoistureInit)")
       ALLOCATE(prlon(nlon,nlat), STAT=ierr)
       IF (ierr/=0) CALL fatal_error("ERROR allocating prlon (soilMoistureInit)")

       CALL apiPrlatlon(nlon, nlat, prlat, prlon, ilatn, ilonn, latni, lonni)

       ALLOCATE(api_us(nlon,nlat,n4us), STAT=ierr)
       IF (ierr/=0) CALL fatal_error("ERROR allocating api_us (soilMoistureInit)")
       ALLOCATE(api_temp(nlon,nlat,n4us), STAT=ierr)
       IF (ierr/=0) CALL fatal_error("ERROR allocating api_temp (soilMoistureInit)")
       ALLOCATE(usdum(n4us), STAT=ierr)
       IF (ierr/=0) CALL fatal_error("ERROR allocating usdum (soilMoistureInit)")

       IF (mchnum==master_num) THEN
          PRINT *, '-------------------------------------- Grid=', ifm
          PRINT *, 'Opening soil moisture data= ', TRIM(usdata)
       ENDIF

       IF (.NOT.theref) THEN ! arquivo .gra

          IF (pref=='us' .OR. pref=='SM' .OR. pref=='gfs025gr.soilm.') THEN   ! arquivo .gra acesso direto
             IF (mchnum==master_num) THEN
                
	       IF (pref=='us' .OR. pref=='SM' )then
		OPEN(2, status='OLD', form='unformatted', access='direct', &
                     recl=4*nlat*nlon*n4us, file=usdata(1:len_trim(usdata)))
                READ(UNIT=2,REC=1) api_us        ! water content
                
	       ELSEIF(pref=='gfs025gr.soilm.') THEN

		OPEN(2, status='OLD', form='unformatted', access='direct', &
                     recl=4*nlat*nlon, file=usdata(1:len_trim(usdata)))
                NREC=0
		DO IREC=1,n4us
		 NREC=NREC+1
		 READ(UNIT=2,REC=NREC) api_us(:,:,n4us-IREC+1)       ! water content
		ENDDO
		!print*, 'soilw0_10cm', MAXVAL(api_us(:,:,4)), MINVAL(api_us(:,:,4))
                !print*, 'soilw10_40cm', MAXVAL(api_us(:,:,3)), MINVAL(api_us(:,:,3))
                !print*, 'soilw40_100cm', MAXVAL(api_us(:,:,2)), MINVAL(api_us(:,:,2))
                !print*, 'soilw100_200cm', MAXVAL(api_us(:,:,1)), MINVAL(api_us(:,:,1))
		
		DO IREC=1,n4us
		 NREC=NREC+1
		 READ(UNIT=2,REC=NREC) api_temp(:,:,n4us-IREC+1)       ! water content
		ENDDO
                !print*, 'tsoil0_10cm', MAXVAL(api_temp(:,:,4)), MINVAL(api_temp(:,:,4))
                !print*, 'tsoil10_40cm', MAXVAL(api_temp(:,:,3)), MINVAL(api_temp(:,:,3))
                !print*, 'tsoil40_100cm', MAXVAL(api_temp(:,:,2)), MINVAL(api_temp(:,:,2))
                !print*, 'tsoil100_200cm', MAXVAL(api_temp(:,:,1)), MINVAL(api_temp(:,:,1))
               	
		!- shifting data from 0->360 to -180->+180       
	        DO K=1,n4us
		 DO J=1,nlat
		  DO i=1,nlon/2
		   DUMMY                 =api_us(i,j,k)
		   api_us(i,j,k)         =api_us(i+nlon/2-1,j,k)
	       	   api_us(i+nlon/2-1,j,k)=DUMMY
		   !
		   DUMMY                   =api_temp(i,j,k)
		   api_temp(i,j,k)         =api_temp(i+nlon/2-1,j,k)
	       	   api_temp(i+nlon/2-1,j,k)=DUMMY
		ENDDO;ENDDO ;ENDDO
	       
	       
	       ENDIF
               CLOSE(2)
             ENDIF
	     
             CALL parf_bcast(api_us, INT(nlon,i8), INT(nlat,i8), INT(n4us,i8), &
                  master_num)
             
	     IF(pref/='gfs025gr.soilm.') THEN
              CALL swap32(api_us, nlat*nlon*n4us) ! Verify before call swap32 - Demerval
              IF ( (api_us(nlat/2,nlon/2,1)<0) .OR. &
                   (api_us(nlat/2,nlon/2,1)>1))     &
                    CALL swap32(api_us, nlat*nlon*n4us)
             ENDIF
	     
             IF (mchnum==master_num) &
                  PRINT*,'--------------------------------------'      

          ELSE     ! arquivo .gra acesso sequencial
             IF (mchnum==master_num) THEN
                OPEN(2, status='OLD', form='unformatted', file=usdata(1:len_trim(usdata)))  
                DO k=1,n4us
                   READ(2) ((api_us(i,j,k),i=1,nlon), j=1,nlat) ! wetness
                ENDDO
                CLOSE(2)
             ENDIF
             CALL parf_bcast(api_us, INT(nlon,i8), INT(nlat,i8), INT(n4us,i8), &
                  master_num)
          ENDIF

       ELSE  ! arquivo .vfm
          IF (mchnum==master_num) THEN
             OPEN(UNIT=2, FILE=usdata(1:len_trim(usdata)), FORM='formatted', STATUS='old')
             CALL vfirec(2, api_us, (nlat*nlon*n4us), 'LIN')
             CLOSE(2)
          ENDIF
          CALL parf_bcast(api_us, INT(nlon,i8), INT(nlat,i8), INT(n4us,i8), &
               master_num)

       ENDIF
       ! geva 20.10.04

!!$       print '(A,5I4)', "DEBUG-ALF:soilMoistureInit:mynum,nlat,nlon=", &
!!$            mynum, nlat,nlon
!!$       print '(A,5I4)', "DEBUG-ALF:soilMoistureInit:mynum,n2,n3,n4us,mzg=", &
!!$            mynum, n2,n3,n4us,mzg
!!$       call flush(6)

       ! Gathering Data
       varn = 'SOIL_WATER'
       call gatherData(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,                    &
            soil_water, GlobalSoilWater)
       varn = 'SOIL_TEXT'
       call gatherData(idim_type, varn, ifm, mzg, nnxp(ifm), nnyp(ifm), &
            npat, nmachs, mchnum, mynum, master_num,                    &
            soil_text, GlobalSoilText)
       varn = 'GLON'
       call gatherData(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
            glon, globalGlon)
       varn = 'GLAT'
       call gatherData(2, varn, ifm, nnxp(ifm), nnyp(ifm), &
            nmachs, mchnum, mynum, master_num,             &
            glat, globalGlat)

       ! loop no dominio global do modelo
       DO j=1,nnyp(ifm)
          DO i=1,nnxp(ifm)                  
             ! evitando pontos fora do dominio da grade de umidade de solo
             IF ( globalGlat(i,j)<latni .OR. globalGlat(i,j)>latnf .OR. &
                  globalGlon(i,j)<lonni .OR. &
		  globalGlon(i,j)>lonnf) CYCLE 
             !print*,'ij=',i,j
             CALL interpolacao(globalGlon(i,j), globalGlat(i,j), nlon, nlat, &
                  prlat, prlon, i1, i2, ic, j1, j2, jc)

             IF (ic>=0 .AND. jc>=0) THEN
                !print*,ic,jc,i,j,ifm
                dlonr = 0.5*(globalGlon(n2,j) - globalGlon(1,j))/&
                     float(nnxp(ifm)-1)
                dlatr = 0.5*(globalGlat(i,n3) - globalGlat(i,1))/&
                     float(nnyp(ifm)-1)
                qi1   = INT(dlonr/ilonn+0.5)
                qi2   = INT(dlonr/ilonn+0.5)
                qj1   = INT(dlatr/ilatn+0.5)
                qj2   = INT(dlatr/ilatn+0.5)

                DO k=1,n4us
                   ncount = 0
                   usdum(k)=0.
                   DO jj=MAX(1,jc-qj1),MIN(nlat,jc+qj2)
                      DO ii=MAX(1,ic-qi1),MIN(nlon,ic+qi2)
                         IF (api_us(ii,jj,k)>1.E-5) THEN
                            DO ipat=2,npat
                               ncount = ncount + 1
                               usdum(k) = usdum(k) + api_us(ii,jj,k)
!!$                               IF (pref=='us' .OR. pref=='SM') THEN
!!$                                  usdum(k) = usdum(k) + api_us(ii,jj,k)   ! umidade lida em m3/m3 (vol.) - v1 (us e SM)
!!$                               ENDIF
!!$                               IF ( pref=='SM_v2.'      .OR. &
!!$                                    pref=='GL_SM.GPCP.' .OR. &
!!$                                    pref=='GL_SM.GPNR.') THEN
!!$                                  !!nsoil = NINT(globalSoilText(k,i,j,ipat))
!!$                                  usdum(k) = usdum(k) + &
!!$                                       api_us(ii,jj,k) !*slmsts(nsoil)  ! umidade lida em % (armazenamento) - v2 (SM_v2.)                            
!!$                               ENDIF
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
                   usdum(k) = usdum(k)/(float(ncount) + 1.E-10)
                ENDDO

                DO k=mzg,1,-1
                   DO kk=n4us,1,-1
		  
                      IF (slz(k)>=slz_us(kk)) THEN
                         DO ipat=2,npat
                            nsoil = NINT(globalSoilText(k,i,j,ipat))
                            
			    !- umidade em mm3/mm3 (GFS case)
			    globalSoilWater(k,i,j,ipat) = usdum(kk+1)
			    !if(ipat==2)print*,"smoist",k,i,j,globalSoilWater(k,i,j,ipat)
                            
			    ! umidade lida em % (armazenamento) 
                            IF ( pref=='SM_v2.'      .OR. &
                                 pref=='GL_SM.GPCP.' .OR. &
                                 pref=='GL_SM.GPNR.'  ) THEN
                               globalSoilWater(k,i,j,ipat) = &
                                    globalSoilWater(k,i,j,ipat)*slmsts(nsoil)
                            ENDIF
                            IF (usdum(kk+1)<1.e-5) &
                                 globalSoilWater(k,i,j,ipat) = slmstr(k)*slmsts(nsoil) !oceano
                                 globalSoilWater(k,i,j,ipat) = MAX(soilcp(nsoil), &
                                                               MIN(globalSoilWater(k,i,j,ipat), &
                                                               slmsts(nsoil)))
                         ENDDO
                         EXIT !GOTO 222
                      ELSEIF (slz(k)<slz_us(1)) THEN
                         DO ipat=2,npat
                            nsoil = NINT(globalSoilText(k,i,j,ipat))
			   
			    !- umidade em mm3/mm3 (GFS case)
                            globalSoilWater(k,i,j,ipat) = usdum(1)
                            !if(ipat==2)print*,"smoist",k,i,j,globalSoilWater(k,i,j,ipat)
			    ! umidade lida em % (armazenamento) - v2 (SM_v2.)
                            IF ( pref=='SM_v2.'      .OR. &
                                 pref=='GL_SM.GPCP.' .OR. &
                                 pref=='GL_SM.GPNR.'  ) THEN
                                 globalSoilWater(k,i,j,ipat) = &
                                    globalSoilWater(k,i,j,ipat)*slmsts(nsoil)
                            ENDIF
                            IF (usdum(1)<1.e-5)   &
                                 globalSoilWater(k,i,j,ipat) = slmstr(k)*slmsts(nsoil)!oceano
                                 globalSoilWater(k,i,j,ipat) = MAX(soilcp(nsoil), &
                                                               MIN(globalSoilWater(k,i,j,ipat), &
                                                               slmsts(nsoil)))
                         ENDDO
                         EXIT !GOTO 222
                      ENDIF
                   ENDDO
                   !!222                CONTINUE
                ENDDO
             ENDIF
          ENDDO
       ENDDO

       DEALLOCATE(api_us, usdum, prlat, prlon,api_temp)

       IF (mchnum==master_num) THEN
!!$          Global_soil_water = RESHAPE(fullField, &
!!$               (/mzg, nnxp(ifm), nnyp(ifm), npat/))
          OPEN(2, status='UNKNOWN', form='unformatted', access='direct', &
               recl=4*nnxp(ifm)*nnyp(ifm)*mzg*npat, file=usmodel(1:len_trim(usmodel)))
          WRITE(UNIT=2,REC=1) globalSoilWater
          CLOSE(2)
       ENDIF
!!$       DEALLOCATE(fullField)
!!$       DEALLOCATE(gathered)
!!$       DEALLOCATE(localChunk)

       ! Scattering Local Data
      CALL mk_4_buff(globalSoilWater(:,:,:,:), soil_water(:,:,:,:), &
            mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)

!- 8/9/2015 srf
!    because the the soil moisture dataset interpolated to the model grid
!    is never checked again, we will allways require the interpolation to
!    be done every model initialization. THis will prevents soil moisture 
!    interpolated for given grid specification being erroneously be used 
!    for a different model configuration
!
!    ELSE
!       IF (mchnum==master_num) THEN
!          PRINT*,'-------------------------------------- Grid=', ifm
!          PRINT*,'Opening soil moisture file= ', TRIM(usmodel)
!          OPEN(2, status='OLD', form='unformatted', access='direct', &
!               recl=4*nnxp(ifm)*nnyp(ifm)*mzg*npat, file=usmodel(1:len_trim(usmodel)))
!          READ(UNIT=2,REC=1) globalSoilWater
!          CLOSE(2)
!          PRINT*,'--------------------------------------'
!       ENDIF
!
!       CALL parf_bcast(globalSoilWater, INT(mzg,i8), INT(nnxp(ifm),i8), &
!            INT(nnyp(ifm),i8), INT(npat,i8), master_num)
!
!       ! Distributing local information about Soil Water
!       CALL mk_4_buff(globalSoilWater(:,:,:,:), soil_water(:,:,:,:), &
!            mzg, nnxp(ifm), nnyp(ifm), npat, mzg, n2, n3, npat, ia, iz, ja, jz)
!
!    ENDIF
!-srf end


    !----- recalculate soil_energy
    c1 = 0.5*cpi

    DO j=1,n3
       DO i=1,n2

          k2      = lpw(i,j)
          pis     = c1*(pi0(k2-1,i,j) + pi0(k2,i,j) + pp(k2-1,i,j) + pp(k2,i,j))
          airtemp = theta(k2,i,j)*pis

          DO ipat=2,npat
             DO k=1,mzg
               nsoil = NINT(soil_text(k,i,j,ipat))
                soil_energy(k,i,j,ipat) = (airtemp - 273.15 + stgoff(k))*   &
                     (slcpd(nsoil) + soil_water(k,i,j,ipat)*4.186e6)      + &
                     soil_water(k,i,j,ipat)*3.34e8
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE soilMoistureInit


!!$  ! Recreating Global Information (Gathering data)
!!$  SUBROUTINE gatherData2D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, nwave, &
!!$       nmachs, mchnum, mynum, master_num,                                   &
!!$       localData2D, globalData2D)
!!$
!!$    USE ReadBcst, ONLY: &
!!$         PreProcAndGather, & ! Subroutine
!!$         LocalSizesAndDisp   ! Subroutine
!!$    USE mem_grid, ONLY : &
!!$         GlobalSizes         ! Subroutine
!!$    USE ParLib, ONLY: &
!!$         parf_bcast ! Subroutine
!!$
!!$    IMPLICIT NONE
!!$    INCLUDE "i8.h"
!!$    ! Arguments:
!!$    INTEGER, INTENT(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
!!$         nwave, nmachs, mchnum, mynum, master_num
!!$    CHARACTER(LEN=16), INTENT(IN) :: varn
!!$    REAL, INTENT(IN)              :: localData2D(:,:)
!!$    REAL, INTENT(OUT)             :: globalData2D(:,:)
!!$    ! Local Variables:
!!$    CHARACTER(LEN=16)  :: localVarn
!!$    INTEGER            :: ierr
!!$    INTEGER, PARAMETER :: idim_type_min = 2
!!$    INTEGER, PARAMETER :: idim_type_max = 7
!!$    INTEGER            :: il1(nmachs)
!!$    INTEGER            :: ir2(nmachs)
!!$    INTEGER            :: jb1(nmachs)
!!$    INTEGER            :: jt2(nmachs)
!!$    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: maxLocalSize
!!$    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxSizeGathered
!!$    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxsizeFullField
!!$    INTEGER            :: globalSize(idim_type_min:idim_type_max)
!!$    REAL, ALLOCATABLE  :: localChunk(:)
!!$    REAL, ALLOCATABLE  :: gathered(:)
!!$    REAL, ALLOCATABLE  :: fullField(:)
!!$
!!$    ! Recreating Global information about Soil Water
!!$    ! grid dependent, field independent constants for gather and unpacking
!!$    ! as a function of idim_type
!!$    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
!!$    maxLocalSize = MAXVAL(localSize(mynum,:))
!!$    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating localChunk (gatherData)")
!!$    ENDIF
!!$    CALL CopyLocalChunk(localData2D(1,1), localChunk, &
!!$         LocalSize(mynum,idim_type))
!!$    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
!!$    maxSizeGathered = MAXVAL(sizeGathered)
!!$    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating gathered (gatherData)")
!!$    ENDIF
!!$    ! grid dependent field sizes as a function of idim_type
!!$    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
!!$    IF (mchnum==master_num) THEN
!!$       sizeFullField(:) = globalSize(:)
!!$    ELSE
!!$       sizeFullField(:) = 1
!!$    END IF
!!$    maxSizeFullField = MAXVAL(sizeFullField)
!!$    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating fullField (gatherData)")
!!$    ENDIF
!!$    localVarn = trim(varn)
!!$    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
!!$         il1, ir2, jb1, jt2, localSize, disp,                 &
!!$         localSize(mynum,idim_type), LocalChunk,              &
!!$         sizeGathered(idim_type), gathered,                   &
!!$         sizeFullField(idim_type), fullField                  )
!!$    
!!$    IF (mchnum==master_num) THEN
!!$       globalData2D = RESHAPE(fullField, (/nnxp, nnyp/))
!!$    ENDIF
!!$
!!$    call parf_bcast(globalData2D, int(nnxp,i8), int(nnyp,i8), master_num)
!!$
!!$    DEALLOCATE(fullField)
!!$    DEALLOCATE(gathered)
!!$    DEALLOCATE(localChunk)
!!$    
!!$  END SUBROUTINE gatherData2D
!!$
!!$  ! Recreating Global Information (Gathering data)
!!$  SUBROUTINE gatherData4D(idim_type, varn, ifm, mzg, nnxp, nnyp, npat, nwave, &
!!$       nmachs, mchnum, mynum, master_num,                                   &
!!$       localData4D, globalData4D)
!!$
!!$    USE ReadBcst, ONLY: &
!!$         PreProcAndGather, & ! Subroutine
!!$         LocalSizesAndDisp   ! Subroutine
!!$    USE mem_grid, ONLY : &
!!$         GlobalSizes         ! Subroutine
!!$    USE ParLib, ONLY: &
!!$         parf_bcast ! Subroutine
!!$
!!$    IMPLICIT NONE
!!$    INCLUDE "i8.h"
!!$    ! Arguments:
!!$    INTEGER, INTENT(IN)           :: idim_type, ifm, mzg, nnxp, nnyp, npat, &
!!$         nwave, nmachs, mchnum, mynum, master_num
!!$    CHARACTER(LEN=16), INTENT(IN) :: varn
!!$    REAL, INTENT(IN)              :: localData4D(:,:,:,:)
!!$    REAL, INTENT(OUT)             :: globalData4D(:,:,:,:)
!!$    ! Local Variables:
!!$    CHARACTER(LEN=16)  :: localVarn
!!$    INTEGER            :: ierr
!!$    INTEGER, PARAMETER :: idim_type_min = 2
!!$    INTEGER, PARAMETER :: idim_type_max = 7
!!$    INTEGER            :: il1(nmachs)
!!$    INTEGER            :: ir2(nmachs)
!!$    INTEGER            :: jb1(nmachs)
!!$    INTEGER            :: jt2(nmachs)
!!$    INTEGER            :: localSize(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: disp(nmachs,idim_type_min:idim_type_max)
!!$    INTEGER            :: maxLocalSize
!!$    INTEGER            :: sizeGathered(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxSizeGathered
!!$    INTEGER            :: sizeFullField(idim_type_min:idim_type_max)
!!$    INTEGER            :: maxsizeFullField
!!$    INTEGER            :: globalSize(idim_type_min:idim_type_max)
!!$    REAL, ALLOCATABLE  :: localChunk(:)
!!$    REAL, ALLOCATABLE  :: gathered(:)
!!$    REAL, ALLOCATABLE  :: fullField(:)
!!$
!!$    ! Recreating Global information about Soil Water
!!$    ! grid dependent, field independent constants for gather and unpacking
!!$    ! as a function of idim_type
!!$    CALL LocalSizesAndDisp(ifm, il1, ir2, jb1, jt2, localSize, disp)
!!$    maxLocalSize = MAXVAL(localSize(mynum,:))
!!$    ALLOCATE(localChunk(maxLocalSize), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating localChunk (gatherData)")
!!$    ENDIF
!!$    CALL CopyLocalChunk(localData4D(1,1,1,1), localChunk, &
!!$         LocalSize(mynum,idim_type))
!!$    sizeGathered(:) = disp(nmachs,:) + localSize(nmachs,:)
!!$    maxSizeGathered = MAXVAL(sizeGathered)
!!$    ALLOCATE(gathered(maxSizeGathered), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating gathered (gatherData)")
!!$    ENDIF
!!$    ! grid dependent field sizes as a function of idim_type
!!$    CALL GlobalSizes(ifm, nmachs, nwave, globalSize)
!!$    IF (mchnum==master_num) THEN
!!$       sizeFullField(:) = globalSize(:)
!!$    ELSE
!!$       sizeFullField(:) = 1
!!$    END IF
!!$    maxSizeFullField = MAXVAL(sizeFullField)
!!$    ALLOCATE(fullField(sizeFullField(idim_type)), stat=ierr)
!!$    IF (ierr/=0) THEN
!!$       CALL fatal_error("Error allocating fullField (gatherData)")
!!$    ENDIF
!!$    localVarn = trim(varn)
!!$    CALL PreProcAndGather(.FALSE., ifm, idim_type, localVarn, &
!!$         il1, ir2, jb1, jt2, localSize, disp,                 &
!!$         localSize(mynum,idim_type), LocalChunk,              &
!!$         sizeGathered(idim_type), gathered,                   &
!!$         sizeFullField(idim_type), fullField                  )
!!$    
!!$    IF (mchnum==master_num) THEN
!!$       globalData4D = RESHAPE(fullField, (/mzg, nnxp, nnyp, npat/))
!!$    ENDIF
!!$
!!$    call parf_bcast(globalData4D, &
!!$         int(mzg,i8), int(nnxp,i8), int(nnyp,i8), int(npat,i8), master_num)
!!$
!!$    DEALLOCATE(fullField)
!!$    DEALLOCATE(gathered)
!!$    DEALLOCATE(localChunk)
!!$    
!!$  END SUBROUTINE gatherData4D
  subroutine StoreNamelistFileAtSoilMoisture(oneNamelistFile)
    implicit none
    type(namelistFile), pointer :: oneNamelistFile
    soil_moist = oneNamelistFile%soil_moist
    soil_moist_fail = oneNamelistFile%soil_moist_fail
    usdata_in = oneNamelistFile%usdata_in
    usmodel_in = oneNamelistFile%usmodel_in
  end subroutine StoreNamelistFileAtSoilMoisture
END MODULE soilMoisture

  !
  ! prlatlon
  !----------------------------------------------------------------
  ! SUB-ROTINA QUE ESTABELECE LATITUDES E LONGITUDES DOS PONTOS DE  
  ! GRADE DO CAMPO DE PRECIPITACAO
  SUBROUTINE apiPrlatlon(nlon, nlat, prlat, prlon, ilatn, ilonn, latni, lonni)
    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN) :: nlon, nlat
    REAL, INTENT(OUT)   :: prlat(nlon,nlat) !(nlon,nlat)
    REAL, INTENT(OUT)   :: prlon(nlon,nlat) !(nlon,nlat)
    REAL, INTENT(IN)    ::  ilatn, ilonn, latni, lonni
    ! Local Variables:
    INTEGER :: i, j

    DO j=1,nlat
       DO i=1,nlon
          prlon(i,j) = lonni + (i-1)*ilonn
          prlat(i,j) = latni + (j-1)*ilatn
       ENDDO
    ENDDO

  END SUBROUTINE apiPrlatlon

  !----------------------------------------------------------------
  ! interpolacao
  !----------------------------------------------------------------
  ! SUB-ROTINA QUE REALIZA INTERPOLACAO ENTRE GRADES (RAMS E UMIDADE DO SOLO)  
  SUBROUTINE interpolacao(glon, glat, nlon, nlat, prlat, prlon, &
       i1, i2, ic, j1, j2, jc)

    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(OUT) :: i1, i2, ic, j1, j2, jc
    REAL, INTENT(IN)     :: glat, glon
    INTEGER, INTENT(IN)  :: nlon, nlat
    REAL, INTENT(IN)     :: prlat(nlon,nlat)
    REAL, INTENT(IN)     :: prlon(nlon,nlat)
    ! Local Variables
    REAL    :: diffx1, diffx2, diffy1, diffy2
    INTEGER :: i, j

    DO i=1,nlon
       IF (glon<=prlon(i,1)) EXIT
    ENDDO
    i2 = i
    i1 = i-1

    DO j=1,nlat
       IF (glat<=prlat(1,j)) EXIT
    ENDDO
    j2 = j
    j1 = j-1

    diffx1 =   glon - prlon(i1,j1)
    diffx2 = -(glon - prlon(i1,j2))
    diffy1 =   glat - prlat(i1,j1)
    diffy2=  -(glat - prlat(i2,j1))

    jc = j1
    ic = i1
    IF (diffx1>diffx2) ic = i2
    IF (diffy1>diffy2) jc = j2

    IF (i1<1 .OR. i1>nlon .OR. j1<1 .OR. j1>nlat) THEN
       ic = -9999
       jc = -9999
    ENDIF

  END SUBROUTINE interpolacao

  !----------------------------------------------------------------

  SUBROUTINE changeDay(idate1, imonth1, iyear1, INT_DIF_TIME, &
       idate2, imonth2, iyear2)
    IMPLICIT NONE
    ! Arguments:
    INTEGER, INTENT(IN)  :: idate1, imonth1, iyear1, INT_DIF_TIME
    INTEGER, INTENT(OUT) :: idate2, imonth2, iyear2
    ! Local Variables:
    INTEGER :: i, increm, DMES(12)

    ! Initiate DMES
    DMES = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)

    iyear2  = iyear1
    imonth2 = imonth1
    idate2  = idate1

    increm  = 1

    IF (INT_DIF_TIME<1) increm = increm*(-1)

    DO i=1,ABS(INT_DIF_TIME)
       idate2 = idate2 + increm
       IF (idate2<1) THEN
          imonth2 = imonth2 + increm
          IF (imonth2<1) THEN
             imonth2 = 12
             iyear2  = iyear2 - 1
          ENDIF
          idate2 = DMES(imonth2)
       ELSEIF (idate2>DMES(imonth2)) THEN
          imonth2 = imonth2 + increm
          IF (imonth2>12) THEN
             imonth2 = 1
             iyear2  = iyear2 + 1
          ENDIF
          idate2 = 1
       ENDIF
    ENDDO

  END SUBROUTINE changeDay


!-srf END MODULE soilMoisture

!============================================================================

SUBROUTINE Swap32(A, N)
  !
  !      REVERSE ORDER OF BYTES IN INTEGER*4 WORD, or REAL*4
  !
  IMPLICIT NONE
  ! Arguments:
  INTEGER, INTENT(IN)            :: n
  INTEGER(kind=4), INTENT(INOUT) :: a(n)
  ! Local Varaibles:
  CHARACTER(LEN=1) :: jtemp(4)
  CHARACTER(LEN=1) :: ktemp
  !
  ! Local variables
  INTEGER :: i, itemp
  
  EQUIVALENCE (jtemp(1), itemp)
  !
  SAVE
  !
  DO i=1,n
     itemp    = a(i)
     ktemp    = jtemp(4)
     jtemp(4) = jtemp(1)
     jtemp(1) = ktemp
     ktemp    = jtemp(3)
     jtemp(3) = jtemp(2)
     jtemp(2) = ktemp
     a(i)     = itemp
  ENDDO
  
END SUBROUTINE Swap32
