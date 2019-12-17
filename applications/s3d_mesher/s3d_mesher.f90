Program e3d_mesher
USE Sigma_User, ONLY : mySetup,mySigma,myProcess
USE iniparser, ONLY : inip_makeDirectory
implicit none
REAL*8 DZi,DZo,DL,Dzz,dx,DA
INTEGER nR,nT,nZ,nN,nM,nP,nT1,nT2,nDEl,deltaNR
INTEGER nX,nY
REAL*8 dbW,dbH,dbL,dSize
REAL*8 xCut, yCut, dalphaCut1,dalphaCut2 
REAL*8,  ALLOCATABLE ::  dCoor(:,:)
INTEGER, ALLOCATABLE ::  kVert(:,:),knpr(:)
INTEGER, ALLOCATABLE :: knpr_inP(:),knpr_inM(:),knpr_outP(:),knpr_outM(:)
INTEGER, ALLOCATABLE :: knpr_zMin(:),knpr_zMax(:)
REAL*8 dBOX(3,2),daux,dnX,dnY,dnz
CHARACTER*(50) CaseFile,command
INTEGER NEL,NVT
character(8)  :: cdate
character(10) :: ctime
character(5)  :: czone
integer,dimension(8) :: values
integer i,idir

LOGICAL bExist
CHARACTER cExtrud3DFile*120

call date_and_time(cdate,ctime,czone,values)

!WRITE(CaseFile,'(5A)') 'Case_',cdate(3:4),cdate(5:6),cdate(7:8),ctime(1:6)
#if defined( WIN32 ) || defined( _WIN32 ) || defined( __WIN32__ )
 CaseFile='_data\meshDir'
WRITE(*,*) 'Windows'
#else
 CaseFile='_data/meshDir'
WRITE(*,*) 'Not Windows'
#endif
WRITE(*,*) CaseFile

CALL ReadPar()


IF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."HOLLOWCYLINDER") THEN
!  WRITE(command,'(2A)') 'mkdir ',TRIM(ADJUSTL(CaseFile))
!  CALL system(TRIM(ADJUSTL(command)))
 call inip_makeDirectory(TRIM(ADJUSTL(CaseFile)))

 nR = mySetup%m_nR
 nT = mySetup%m_nT
 nZ = mySetup%m_nZ
 Dzo = mySigma%Dz_out
 Dzi = mySigma%Dz_in
 DL  = mySigma%L

 IF (mySetup%MeshResolution.eq.1) THEN
  if (mySetup%m_nR.eq.0)  nR  = 2
 END IF
 IF (mySetup%MeshResolution.eq.2) THEN
  if (mySetup%m_nR.eq.0)  nR  =  3
 END IF
 IF (mySetup%MeshResolution.eq.3)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  4
 END IF
 IF (mySetup%MeshResolution.eq.4)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  5
 END IF
 IF (mySetup%MeshResolution.eq.5)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  6
 END IF
 
 IF (mySetup%MeshResolution.ne.0)  THEN
  nT = INT(1.0d0*DBLE(nR)*(3.14*(Dzo+Dzi)*0.5)/(Dzo-Dzi))
  nZ = INT(1.0d0*DBLE(nR)*mySigma%L/(DZo-DZi))
 END IF

 IF (MOD(nZ,2).EQ.1) THEN
  WRITE(*,*) "number of elements in Z is uneven and is to be set to: ", nZ+1
  nZ = nZ + 1
 END IF
 
 WRITE(*,'(A,3I0)') "Resolution of the hollow cylinder nR,nT,nZ: ", nR,nT,nZ
 
 CALL SetUpMesh_HC()

 CALL Parametrization_HC()

 CALL WriteMesh_YXZ()
 
ELSEIF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."FULLCYLINDER") THEN

!  WRITE(command,'(2A)') 'mkdir ',TRIM(ADJUSTL(CaseFile))
!  CALL system(TRIM(ADJUSTL(command)))
 call inip_makeDirectory(TRIM(ADJUSTL(CaseFile)))

 nN = mySetup%m_nT;     nM = mySetup%m_nR;    nZ = mySetup%m_nZ;  nP = mySetup%m_nP  
 Dzo = mySigma%Dz_out;  DL  = mySigma%L
 
 IF (MOD(nZ,2).EQ.1) THEN
  WRITE(*,*) "number of elements in Z is uneven and is to be set to: ", nZ+1
  nZ = nZ + 1
 END IF
 
 CALL SetUpMesh_FC()
 
 CALL Parametrization_FC()

 CALL WriteMesh_YXZ()

ELSEIF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."TWINSCREW") THEN

!  WRITE(command,'(2A)') 'mkdir ',TRIM(ADJUSTL(CaseFile))
!  CALL system(TRIM(ADJUSTL(command)))
 call inip_makeDirectory(TRIM(ADJUSTL(CaseFile)))

 DZo = mySigma%Dz_out;  DZi = mySigma%Dz_in; DA  = mySigma%a
 DL  = mySigma%L;       dx  = mySigma%W;     DZz = mySigma%Dzz

 nR  = mySetup%m_nR
 nZ  = mySetup%m_nZ
 nT1 = mySetup%m_nT1
 nT2 = mySetup%m_nT2

 !!! calculate the relevant sizes
 IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."ROUND") THEN
  xCut = (DZz**2d0 - DZo**2d0 + DA**2d0)/(4d0*DA)
  yCut = SQRT((0.5*DZo)**2d0 - (xCut-0.5*DA)**2d0) 
 ENDIF
 IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."STRAIGHT") THEN
  yCut = SQRT((0.5*DZo)**2d0 - (0.5*DA)**2d0) + dx
  xCut = -SQRT((0.5*DZo)**2d0 - yCut**2d0) + 0.5*DA
 END IF
 dalphaCut1 = +(-2d0*DATAN(1d0) + DATAN((xCut-0.5*DA)/yCut))
 dalphaCut2 = -(-2d0*DATAN(1d0) + DATAN((xCut-0.5*DA)/yCut))
 
 !barrel thickness
 dbW = (DZo-DZi)*0.5d0
 !barrel height through the V cut
 dbH = 2d0*yCut
 !outer barrel circumference
 dbL = ABS(dalphaCut1)*DZo

 deltaNR = 0
 

 IF (mySetup%MeshResolution.ne.0) THEN
  nR = 2

 1 CONTINUE
  dSize = dbW / DBLE(nR)
  nT1 = MAX(1,INT(0.80d0*dbH/dSize))
  nT2 = MAX(1,INT(0.30d0*dbL/dSize))

  nDEL=2*NR*(NT1+NT2) + NT1
  if (nDEl.lt.72) THEN
   nR = nR + 1
   GOTO 1
  end if
  
  !!! shift the resolution to this factor higher (because in 2D there were too few elements)
  deltaNR = nR - 2
 END IF

 IF (mySetup%MeshResolution.eq.1) THEN
  if (mySetup%m_nR.eq.0)  nR  = 2
!  if (mySetup%m_nT1.eq.0) nT1 = 5
!  if (mySetup%m_nT2.eq.0) nT2 = 12
 END IF
 IF (mySetup%MeshResolution.eq.2) THEN
  if (mySetup%m_nR.eq.0)  nR  = 3
!  if (mySetup%m_nT1.eq.0) nT1 = 6
!  if (mySetup%m_nT2.eq.0) nT2 = 17
 END IF
 IF (mySetup%MeshResolution.eq.3)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  4
!  if (mySetup%m_nT1.eq.0) nT1 =  7
!  if (mySetup%m_nT2.eq.0) nT2 = 22 
 END IF
 IF (mySetup%MeshResolution.eq.4)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  5
!  if (mySetup%m_nT1.eq.0) nT1 =  8
!  if (mySetup%m_nT2.eq.0) nT2 = 28 
 END IF
 IF (mySetup%MeshResolution.eq.5)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  6
!  if (mySetup%m_nT1.eq.0) nT1 =  10
!  if (mySetup%m_nT2.eq.0) nT2 =  35 
 END IF
 IF (mySetup%MeshResolution.eq.6)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  7
!  if (mySetup%m_nT1.eq.0) nT1 =  10
!  if (mySetup%m_nT2.eq.0) nT2 =  35 
 END IF
 IF (mySetup%MeshResolution.ge.7)  THEN
  if (mySetup%m_nR.eq.0)  nR  =  8
!  if (mySetup%m_nT1.eq.0) nT1 =  10
!  if (mySetup%m_nT2.eq.0) nT2 =  35 
 END IF
 
 IF (mySetup%MeshResolution.ne.0) THEN
  nR = nR + deltaNR
  dSize = dbW / DBLE(nR)
  nT1 = MAX(1,INT(0.80d0*dbH/dSize))
  nT2 = MAX(1,INT(0.30d0*dbL/dSize))
 END IF
 
 !2DMesh resolution
 nDEL=2*NR*(NT1+NT2) + NT1
!  if (nDEl.lt.60) THEN
!   mySetup%MeshResolution = mySetup%MeshResolution + 1
!   mySetup%m_nR = 0
!   GOTO 1
!  end if
 
 if (mySetup%m_nZ.eq.0) nZ = INT(1.0d0*DBLE(nR)*mySigma%L/(DZo-DZi))
 
 IF (MOD(nZ,2).EQ.1) THEN
  WRITE(*,*) "number of elements in Z is uneven and is to be set to: ", nZ+1
  nZ = nZ + 1
 END IF
 write(*,*) dbW,dbH,dbL,dSize
 WRITE(*,'(A,I0,6(","I0))') "resolution parameters are [nR,nT1,nT2,nZ] + n2D_EL + ndRes: ",nR,nT1,nT2,nZ,nDEl,deltaNR,nDEl*nZ
 
 CALL SetUpMesh_TSE()

 CALL Parametrization_TSE()

 CALL WriteMesh_YXZ()

ELSEIF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."BOX") THEN

 dBOX = mySetup%m_BOX
 if (mySetup%nBoxElem.le.0) then
  nX = mySetup%m_nX
  nY = mySetup%m_nY
  nZ = mySetup%m_nZ
 else
  iDir = -1
  daux = 1d30
  do i=1,3
   if ((mySetup%m_BOX(i,2)-mySetup%m_BOX(i,1)).lt.daux) then
    iDir = i
    daux = mySetup%m_BOX(i,2)-mySetup%m_BOX(i,1)
   end if
   write(*,*) mySetup%m_BOX(i,2)-mySetup%m_BOX(i,1)
  end do

  dnX = (mySetup%m_BOX(1,2)-mySetup%m_BOX(1,1))/(mySetup%m_BOX(iDir,2)-mySetup%m_BOX(iDir,1))
  dnY = (mySetup%m_BOX(2,2)-mySetup%m_BOX(2,1))/(mySetup%m_BOX(iDir,2)-mySetup%m_BOX(iDir,1))
  dnZ = (mySetup%m_BOX(3,2)-mySetup%m_BOX(3,1))/(mySetup%m_BOX(iDir,2)-mySetup%m_BOX(iDir,1))
 write(*,*) 'nX,nY,nZ,nEl'
 write(*,*) dnX,dnY,dnz
  
  iDir = nint((mySetup%nBoxElem/(dnx*dny*dnz))**0.3333d0)
  nX = iDir*nint(dnx)
  nY = iDir*nint(dny)
  nZ = iDir*nint(dnz)
 end if

 write(*,*) 'nX,nY,nZ,nEl'
 write(*,*) nX,nY,nZ,nX*nY*nZ

 CALL SetUpMesh_BOX()
 
 CALL Parametrization_BOX()
 
 CALL WriteMesh_XYZ()
 
ELSEIF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."OFF") THEN

 WRITE(*,*) "No mesh is being created ... (Hexmesher is turned OFF) "

ELSE

 WRITE(*,*) "No mesh is being created ... (unknown meshing keyword:'",ADJUSTL(TRIM(mySetup%cMesher)),"')"
 
END IF

 CONTAINS 
!-----------------------------------------------
SUBROUTINE ReadPar()

 cExtrud3DFile = '_data/Extrud3D_0.dat'
INQUIRE(file=cExtrud3DFile,Exist=bExist)
if (bExist) then
 CALL ReadS3Dfile(cExtrud3DFile)
ELSE
 write(*,*) 'file: "',ADJUSTL(TRIM(cExtrud3DFile)),'" does not exist!'
END IF

END SUBROUTINE ReadPar
!-----------------------------------------------
SUBROUTINE SetUpMesh_FC()
REAL*8 dAlpha,daux,PX1,PX2,PY,PZ,dR
INTEGER i,j,k,kk,nhalf
REAL*8 PX_max,PX_min,PX_mid,PY_max,PY_min,PY_mid,ratio,dZOO
INTEGER i1,i2,i3,i4,i5,i6,i7,i8

REAL*8 P(2),dRadius
integer ni,nnz,ll,mm1,mm2,iel

!NVT = (1*(NR+1)*(nT))*(NZ+1)
NVT = (NP*(NN*NN + NN) + NM*2*NN*NP +1 )*(NZ+1)
NEL = NZ*NP*(NN*NN + NM*2*NN)
!NEL = (1*NR*(nT) )*NZ

write(*,*) 'NN,NM,NZ,NP ',NN,NM,NZ,NP
write(*,*) 'NEL,NVT ',NEL,NVT

ALLOCATE(dcoor(3,NVT))
ALLOCATE(kVert(8,NEL))
ALLOCATE(knpr(NVT))
knpr = 0

nnz=nz
if (nz.eq.0) nnz=1

kk = 1

DO k=0,NZ

 dCoor(:,kk) = [0d0, 0d0 ,DL*DBLE(k)/DBLE(nnZ)]
 kk = kk + 1
 
 DO j=1,NN
  
  ni = 2*NP*j

  DO  i = 1,2*NP*j
  
   P = [0d0,(Dzo/2d0)*(DBLE(NN)/DBLE(NM+NN))*DBLE(j)/DBLE(NN)]
   
   dAlpha = 8d0*DATAN(1d0) * dble(i-1)/dble(ni)
   dRadius = 1d0 !+ 0.25d0*DBLE(ABS(j-mod(i,2*j)))
!    write(*,*) i,dRadius
   P = dRadius * P

   
   dCoor(:,kk) = [P(1)*cos(dAlpha) - P(2)*sin(dAlpha),&
                  P(1)*sin(dAlpha) + P(2)*cos(dAlpha),&
                  DL*DBLE(k)/DBLE(nnZ)]
   kk = kk + 1

  END DO
 END DO
 
 DO j=1,NM

  ni = 2*NP*NN

  DO  i = 1,2*NP*NN
   dAlpha = 8d0*DATAN(1d0) * dble(i-1)/dble(ni)
   dRadius = 1d0 !+ 0.25d0*DBLE(ABS(j-mod(i,2*j)))
   P = [0d0,((Dzo/2d0)*(DBLE(NN)/DBLE(NM+NN)))+(Dzo/2d0)*(DBLE(NM)/DBLE(NM+NN))*DBLE(j)/DBLE(NM)]
   P = dRadius * P
   
   dCoor(:,kk) = [P(1)*cos(dAlpha) - P(2)*sin(dAlpha),&
                  P(1)*sin(dAlpha) + P(2)*cos(dAlpha),&
                  DL*DBLE(k)/DBLE(nnZ)]
   kk = kk + 1
  END DO
 END DO
 
END DO


iel=0

do k=1,NZ
! k=1

mm1 = (k-1)*(NP*(NN*NN + NN) + NM*2*NN*NP +1 )
mm2 =    k *(NP*(NN*NN + NN) + NM*2*NN*NP +1 )
kk =  NP*2+1
ll =  1


DO  j = 1,NP-1
 i=2*(j-1)+1
 iel = iel + 1
 KVERT(:,iel) = [mm1+i+1,mm1+i+2,mm1+i+3,mm1+1,i+1+mm2,i+2+mm2,i+3+mm2,1+mm2]
END DO
j=NP
i=2*(j-1)+1
iel = iel + 1
KVERT(:,iel) = [mm1+i+1,mm1+i+2,mm1+2,mm1+1,i+1+mm2,i+2+mm2,2+mm2,1+mm2]

DO j=2,NN
 DO  i = 1,2*NP*j
  if (mod(i,2*j).ne.j.and.mod(i,2*j).ne.j+1) then
   if (i.lt.2*NP*j) then
    kk = kk + 1 
    ll = ll + 1 
    iel = iel + 1
    KVERT(:,iel) = [mm1+kk,mm1+kk+1,mm1+ll+1,mm1+ll, kk+mm2,kk+1+mm2,ll+1+mm2,ll+mm2]
   else
    kk = kk + 1 
    ll = ll + 1 
    iel = iel + 1
    KVERT(:,iel) = [mm1+kk,mm1+ll+1,mm1+ll+1-2*NP*(j-1),mm1+ll, kk+mm2,ll+1+mm2,ll+1-2*NP*(j-1)+mm2,ll+mm2]
   end if
  else
   if (mod(i,2*j).eq.j) then
    kk = kk + 1
    iel = iel + 1
    KVERT(:,iel) = [mm1+kk,mm1+kk+1,mm1+kk+2,mm1+ll+1, kk+mm2,kk+1+mm2,kk+2+mm2,ll+1+mm2]
   else
    kk = kk + 1
   end if
  end if
 END DO 
END DO 

DO j=1,NM
 DO  i = 1,2*NP*NN-1
    kk = kk + 1
    ll = ll + 1 
    iel = iel + 1
    KVERT(:,iel) = [mm1+kk,mm1+kk+1,mm1+ll+1,mm1+ll, kk+mm2,kk+1+mm2,ll+1+mm2,ll+mm2]
 END DO
 i=2*NP*NN
 kk = kk + 1
 ll = ll + 1 
 iel = iel + 1
 KVERT(:,iel) = [mm1+kk,mm1+ll+1,mm1+ll+1-i,mm1+ll, kk+mm2,ll+1+mm2,ll+1-i+mm2,ll+mm2]
END DO

END DO

END SUBROUTINE SetUpMesh_FC
!-----------------------------------------------
SUBROUTINE SetUpMesh_BOX()
implicit none
! REAL*8 dAlpha,daux,PX1,PX2,PY,PZ,dR
! INTEGER i,j,k,kk,nhalf
! REAL*8 dalphaCut1,dalphaCut2 
! REAL*8 PX_max,PX_min,PX_mid,PY_max,PY_min,PY_mid,ratio,dZOO
INTEGER i,j,k,kk
INTEGER i1,i2,i3,i4,i5,i6,i7,i8
REAL*8 PX,PY,PZ

nel = 0

NVT = (nX+1)*(nY+1)*(NZ+1)
NEL =  nx*nY*NZ

ALLOCATE(dcoor(3,NVT))
ALLOCATE(kVert(8,NEL))
ALLOCATE(knpr(NVT))
knpr = 0

kk = 1

DO k=0,NZ
 DO j=0,nY
  DO i=0,nX

   PX = dBOX(1,1) + (DBLE(i)/DBLE(NX))*(dBOX(1,2)-dBOX(1,1))
   PY = dBOX(2,1) + (DBLE(j)/DBLE(NY))*(dBOX(2,2)-dBOX(2,1))
   PZ = dBOX(3,1) + (DBLE(k)/DBLE(NZ))*(dBOX(3,2)-dBOX(3,1))

   dCoor(:,kk) = [PX,PY,PZ]
   kk = kk + 1

  END DO

 END DO
END DO

kk = 1
DO k=1,nz
 DO j=1,nY
  DO i=1,nX

   i1 = (k-1)*((nY+1)*(nX+1)) + (j-1)*(nX+1) + i
   i2 = i1 + 1
   IF (MOD(i1,nX+1).eq.0) i2=i1-(nX+1)+1
   i3 = i1 + (nX+1) + 1
   IF (MOD(i1,nX+1).eq.0) i3=i3-(nX+1)
   i4 = i1 + (nX+1)
   i5 = i1 + (nY+1)*(nX+1)
   i6 = i5 + 1
   IF (MOD(i1,nX+1).eq.0) i6=i6-(nX+1)
   i7 = i5 + (nX+1) + 1
   IF (MOD(i1,nX+1).eq.0) i7=i7-(nX+1)
   i8 = i5 + (nX+1)

   KVERT(:,kk) = [i1,i2,i3,i4,i5,i6,i7,i8]
!    write(*,*) [i1,i2,i3,i4,i5,i6,i7,i8]

   kk = kk + 1
  END DO
 END DO
END DO

END SUBROUTINE SetUpMesh_BOX
!-----------------------------------------------
SUBROUTINE SetUpMesh_HC()
REAL*8 dAlpha,daux,PX1,PX2,PY,PZ,dR
INTEGER i,j,k,kk,nhalf
REAL*8 dalphaCut1,dalphaCut2 
REAL*8 PX_max,PX_min,PX_mid,PY_max,PY_min,PY_mid,ratio,dZOO
INTEGER i1,i2,i3,i4,i5,i6,i7,i8

nel = 0

NVT = (1*(NR+1)*(nT))*(NZ+1)
NEL = (1*NR*(nT) )*NZ

ALLOCATE(dcoor(3,NVT))
ALLOCATE(kVert(8,NEL))
ALLOCATE(knpr(NVT))
knpr = 0

dalphaCut1 = 0d0
dalphaCut2 = 8d0*DATAN(1d0)

kk = 1

DO k=0,NZ
 DO j=0,NR

  DO i=1,nT
   dAlpha = dalphaCut1 + DBLE(i)*(dalphaCut2-dalphaCut1)/(nT)

   daux = DBLE(j)/DBLE(NR)
   dR = 0.5d0*(DZi +  daux*(DZo-DZi))

   PX1 = 0d0 + dR*dCOS(dAlpha)
   PY = dR*dSIN(dAlpha)
   PZ = DL*DBLE(k)/DBLE(NZ)

   dCoor(:,kk) = [PX1,PY,PZ]
   kk = kk + 1

  END DO

 END DO
END DO

kk = 1
DO k=1,nz
 DO j=1,NR
  DO i=1,nT

   i1 = (k-1)*((NR+1)*(nT)) + (j-1)*(nT) + i
   i2 = i1 + 1
   IF (MOD(i1,nT).eq.0) i2=i1-(nT)+1
   i3 = i1 + (nT) + 1
   IF (MOD(i1,nT).eq.0) i3=i3-(nT)
   i4 = i1 + (nT)
   i5 = i1 + (NR+1)*(nT)
   i6 = i5 + 1
   IF (MOD(i1,nT).eq.0) i6=i6-(nT)
   i7 = i5 + (nT) + 1
   IF (MOD(i1,nT).eq.0) i7=i7-(nT)
   i8 = i5 + (nT)

   KVERT(:,kk) = [i1,i2,i3,i4,i5,i6,i7,i8]

   kk = kk + 1
  END DO
 END DO
END DO

END SUBROUTINE SetUpMesh_HC
!-----------------------------------------------
SUBROUTINE SetUpMesh_TSE()
REAL*8 dAlpha,daux,PX1,PX2,PY,PZ,dR
INTEGER i,j,k,kk,nhalf
REAL*8 dalphaCut1,dalphaCut2 
REAL*8 PX_max,PX_min,PX_mid,PY_max,PY_min,PY_mid,ratio,dZOO
INTEGER i1,i2,i3,i4,i5,i6,i7,i8

nel = 0

NVT = (2*(NR+1)*(NT1+NT2))*(NZ+1)
NEL = (2*NR*(NT1+NT2) + NT1)*NZ
!NEL = 2*NR*(NT1+NT2)*NZ

ALLOCATE(dcoor(3,NVT))
ALLOCATE(kVert(8,NEL))
ALLOCATE(knpr(NVT))
knpr = 0

 !!! calculate the relevant sizes
 IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."ROUND") THEN
  xCut = (DZz**2d0 - DZo**2d0 + DA**2d0)/(4d0*DA)
  yCut = SQRT((0.5*DZo)**2d0 - (xCut-0.5*DA)**2d0) 
 ENDIF
 IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."STRAIGHT") THEN
  yCut = SQRT((0.5*DZo)**2d0 - (0.5*DA)**2d0) + dx
  xCut = -SQRT((0.5*DZo)**2d0 - yCut**2d0) + 0.5*DA
 END IF
 dalphaCut1 = +(-2d0*DATAN(1d0) + DATAN((xCut-0.5*DA)/yCut))
 dalphaCut2 = -(-2d0*DATAN(1d0) + DATAN((xCut-0.5*DA)/yCut))
! write(*,*) DZO,DZI,DA,DL,dx,DZZ,nR,NT1,nt2,nz
! WRITE(*,*) ADJUSTL(TRIM(mySigma%cZwickel))
! WRITE(*,*) xCut,yCut,dalphaCut1,dalphaCut2
! pause
! pause
kk = 1

nhalf = nvt/2

DO k=0,NZ
 DO j=0,NR

  DO i=1,NT2
   dAlpha = dalphaCut1 + DBLE(i)*(dalphaCut2-dalphaCut1)/(NT2)

   daux = DBLE(j)/DBLE(NR)
   dR = 0.5d0*(DZi +  daux*(DZo-DZi))

   PX1 =  0.5*DA + dR*dCOS(dAlpha)
   PX2 = -0.5*DA - dR*dCOS(dAlpha)
   PY = dR*dSIN(dAlpha)
   PZ = DL*DBLE(k)/DBLE(NZ)

   dCoor(:,kk) = [PX1,PY,PZ]
   dCoor(:,nhalf+kk) = [PX2,PY,PZ]
!    WRITE(*,'(3ES11.4)') [PX,PY,PZ]
   kk = kk + 1

  END DO

  DO i=1,NT1
   dAlpha = dalphaCut2 + DBLE(i)*(8d0*DATAN(1d0)-(dalphaCut2-dalphaCut1))/(NT1)

   PX_max = 0.5*DA + 0.5*DZo*dCOS(dAlpha)
   PY_max =          0.5*DZo*dSIN(dAlpha)
   PX_min = 0.5*DA
   PY_min = 0.0      

   ratio  = 1.0-(PX_max-xCut)/(PX_Max-PX_Min)
   PY_mid = PY_min + ratio*(PY_max-PY_min)
   dZOO   = DZo*PY_mid/PY_max
!    write(*,*) dzoo,ratio
!    pause

   daux = DBLE(j)/DBLE(NR)
   dR = 0.5d0*(DZi +  daux*(DZOO-DZi))

   PX1 =  0.5*DA + dR*dCOS(dAlpha)
   PX2 = -0.5*DA - dR*dCOS(dAlpha)
  

   PY = dR*dSIN(dAlpha)
   PZ = DL*DBLE(k)/DBLE(NZ)

   dCoor(:,kk) = [PX1,PY,PZ]
   dCoor(:,nhalf+kk) = [PX2,PY,PZ]
!    WRITE(*,'(3ES11.4)') [PX,PY,PZ]
   kk = kk + 1

  END DO

 END DO
END DO

nhalf = NR*(NT1+NT2)*NZ

kk = 1
DO k=1,nz
 DO j=1,NR
  DO i=1,NT1+NT2

   i1 = (k-1)*((NR+1)*(NT1+NT2)) + (j-1)*(NT1+NT2) + i
   i2 = i1 + 1
   IF (MOD(i1,NT1+NT2).eq.0) i2=i1-(NT1+NT2)+1
   i3 = i1 + (NT1+NT2) + 1
   IF (MOD(i1,NT1+NT2).eq.0) i3=i3-(NT1+NT2)
   i4 = i1 + (NT1+NT2)
   i5 = i1 + (NR+1)*(NT1+NT2)
   i6 = i5 + 1
   IF (MOD(i1,NT1+NT2).eq.0) i6=i6-(NT1+NT2)
   i7 = i5 + (NT1+NT2) + 1
   IF (MOD(i1,NT1+NT2).eq.0) i7=i7-(NT1+NT2)
   i8 = i5 + (NT1+NT2)

   KVERT(:,kk) = [i1,i2,i3,i4,i5,i6,i7,i8]
   KVERT(:,nhalf+kk) = [nvt/2+i1,nvt/2+i2,nvt/2+i3,nvt/2+i4,nvt/2+i5,nvt/2+i6,nvt/2+i7,nvt/2+i8]
!   WRITE(*,'(8I8)') 
   kk = kk + 1
  END DO
 END DO
END DO

kk = 1
 DO k=1,NZ
  DO i=1,NT1

   i1 = (k-1)*((NR+1)*(NT1+NT2)) + NR*(NT1+NT2) + NT2 + i -1
   i2 = i1 + 1
   i3 = i2 + nvt/2
   i4 = i3 -1
   i5 = (k-0)*((NR+1)*(NT1+NT2)) + NR*(NT1+NT2) + NT2 + i - 1
   i6 = i5 + 1
   i7 = i6 + nvt/2
   i8 = i7 -1

   KVERT(:,2*nhalf+kk) = [i1,i2,i3,i4,i5,i6,i7,i8]
!   write(*,'(8I8)') [i1,i2,i3,i4,i5,i6,i7,i8]
   kk = kk + 1
  END DO
 END DO

END SUBROUTINE SetUpMesh_TSE
!-----------------------------------------------
SUBROUTINE Parametrization_FC()
INTEGER i,j,k
INTEGER nHalf,ivt

!! Outer cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/out.par',UNIT=1)
WRITE(1,*) 2*NN*NP*(NZ+1), ' Wall'
WRITE(1,'(A2,4E14.6,A)') "'7", 0e0,  0e0, 0e0, 0.5*DZo, " 1.0 1.0 0.0'"
k = 0
DO j=1,nZ+1
 DO i=1,2*NN*NP
  ivt = (NP*(NN*NN + NN) + NM*2*NN*NP +1 )*(j)- 2*NN*NP + i
!   ivt = (NR+1)*(nT)*(j-1) + NR*(nT) + i
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
!   WRITE(1,'(I6)') ivt
  k = k + 1
END DO
END DO
CLOSE(1)

!! Outflow
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z+.par',UNIT=1)
WRITE(1,*) NP*(NN*NN + NN) + NM*2*NN*NP +1, ' Outflow'
WRITE(1,'(A,E14.6,A)') "'4 0.0 0.0 1.0", -DL, "'"
k = 0
DO j=1,NP*(NN*NN + NN) + NM*2*NN*NP +1
  ivt = (NP*(NN*NN + NN) + NM*2*NN*NP +1)*NZ + j
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
END DO
CLOSE(1)
! 
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z-.par',UNIT=1)

IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
 WRITE(1,*) NP*(NN*NN + NN) + NM*2*NN*NP +1, ' Outflow'
ELSE
 WRITE(1,*) NP*(NN*NN + NN) + NM*2*NN*NP +1, ' Inflow10'
END IF


WRITE(1,'(A,E14.6,A)') "'4 0.0 0.0 1.0", 0d0, "'"
k = 0
DO j=1,NP*(NN*NN + NN) + NM*2*NN*NP +1
  ivt = j
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
END DO
CLOSE(1)

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/file.prj',UNIT=2)
WRITE(2,'(A)') 'Mesh.tri'
WRITE(2,'(A)') 'out.par'
WRITE(2,'(A)') 'z-.par'
WRITE(2,'(A)') 'z+.par'
CLOSE(2)

END SUBROUTINE Parametrization_FC
!-----------------------------------------------
SUBROUTINE Parametrization_HC()
INTEGER i,j,k
INTEGER nHalf,ivt
character*200 cFormat
! INTEGER, ALLOCATABLE :: knpr_inP(:),knpr_inM(:),knpr_outP(:),knpr_outM(:)
! INTEGER, ALLOCATABLE :: knpr_zMin(:),knpr_zMax(:)

nHalf = nvt/2

!! Inner cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/in.par',UNIT=1)
WRITE(1,*) (nT)*(NZ+1), ' Inflow13'

if (mySigma%InnerDiamNParam.eq.0) then
 WRITE(1,'(A2,4E14.6,A)') "'7", 0e0,  0e0, 0e0, 0.5*DZi, " 1.0 1.0 0.0'"
else
 cFormat = '(A,I2,A,xxE14.6,A)'
 write(cFormat(9:10),'(I2.2)') mySigma%InnerDiamNParam*2
 WRITE(1,cFormat) "'",mySigma%InnerDiamNParam*2+6," 0.0 0.0 0.0", mySigma%InnerDiamDParam,mySigma%InnerDiamZParam, " 1.0 1.0 0.0'"
end if

k = 0
DO j=1,nZ+1
 DO i=1,nT
  ivt = (NR+1)*(nT)*(j-1) + i
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
END DO
END DO
CLOSE(1)


!! Outer cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/out.par',UNIT=1)
WRITE(1,*) (nT)*(NZ+1), ' Wall'
WRITE(1,'(A2,4E14.6,A)') "'7", 0e0,  0e0, 0e0, 0.5*DZo, " 1.0 1.0 0.0'"
k = 0
DO j=1,nZ+1
 DO i=1,nT
  ivt = (NR+1)*(nT)*(j-1) + NR*(nT) + i
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
END DO
END DO
CLOSE(1)


!! Outflow
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z+.par',UNIT=1)

IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
 WRITE(1,*) (nT)*(NR+1), ' Outflow'
ELSE
 WRITE(1,*) (nT)*(NR+1), ' Symmetry110'
END IF
WRITE(1,'(A,E14.6,A)') "'4 0.0 0.0 1.0", -DL, "'"
k = 0
DO j=1,nR+1
 DO i=1,nT
  ivt = (NR+1)*(nT)*(nZ) + (j-1)*(nT) + i
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
 END DO
END DO
CLOSE(1)

!! Inflow
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z-.par',UNIT=1)
IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
 WRITE(1,*) (nT)*(NR+1), ' Outflow'
ELSE
 WRITE(1,*) (nT)*(NR+1), ' Inflow10'
END IF
WRITE(1,'(A,E14.6,A)') "'4 0.0 0.0 1.0", 0.0, "'"
k = 0
DO j=1,nR+1
 DO i=1,nT
  ivt = (j-1)*(nT) + i
  knpr(ivt) = 1
  WRITE(1,'(I6)') ivt
  k = k + 1
 END DO
END DO
CLOSE(1)

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/file.prj',UNIT=1)
WRITE(1,'(A)') 'Mesh.tri'
WRITE(1,'(A)') 'in.par'
WRITE(1,'(A)') 'out.par'
WRITE(1,'(A)') 'z-.par'
WRITE(1,'(A)') 'z+.par'
CLOSE(1)

END SUBROUTINE Parametrization_HC
!-----------------------------------------------
SUBROUTINE Parametrization_BOX()
implicit none
INTEGER j

!! Inner cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/x-.par',UNIT=1)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/x+.par',UNIT=2)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/y-.par',UNIT=3)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/y+.par',UNIT=4)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z-.par',UNIT=5)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z+.par',UNIT=6)

WRITE(1,*) (nY+1)*(NZ+1), ' Wall'
WRITE(1,'(A2,4E14.6,A)') "'4", 1e0, 0e0,  0e0, -dBox(1,1),"'"
WRITE(2,*) (nY+1)*(NZ+1), ' Wall'
WRITE(2,'(A2,4E14.6,A)') "'4", 1e0, 0e0,  0e0, -dBox(1,2),"'"
WRITE(3,*) (nX+1)*(NZ+1), ' Wall'
WRITE(3,'(A2,4E14.6,A)') "'4", 0e0, 1e0,  0e0, -dBox(2,1),"'"
WRITE(4,*) (nX+1)*(NZ+1), ' Wall'
WRITE(4,'(A2,4E14.6,A)') "'4", 0e0, 1e0,  0e0, -dBox(2,2),"'"
WRITE(5,*) (nX+1)*(NY+1), ' Wall'
WRITE(5,'(A2,4E14.6,A)') "'4", 0e0, 0e0,  1e0, -dBox(3,1),"'"
WRITE(6,*) (nX+1)*(NY+1), ' Wall'
WRITE(6,'(A2,4E14.6,A)') "'4", 0e0, 0e0,  1e0, -dBox(3,2),"'"

DO j=1,nvt
!  write(*,*) ABS(dcoor(1,j)-dBox(1,1))
 if (ABS(dcoor(1,j)-dBox(1,1)).lt.1d-8) write(1,'(i6)') j
 if (ABS(dcoor(1,j)-dBox(1,2)).lt.1d-8) write(2,'(i6)') j
 if (ABS(dcoor(2,j)-dBox(2,1)).lt.1d-8) write(3,'(i6)') j
 if (ABS(dcoor(2,j)-dBox(2,2)).lt.1d-8) write(4,'(i6)') j
 if (ABS(dcoor(3,j)-dBox(3,1)).lt.1d-8) write(5,'(i6)') j
 if (ABS(dcoor(3,j)-dBox(3,2)).lt.1d-8) write(6,'(i6)') j
END DO

CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)
CLOSE(6)

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/file.prj',UNIT=1)
write(1,'(A)') 'Mesh.tri'
write(1,'(A)') 'x-.par'
write(1,'(A)') 'x+.par'
write(1,'(A)') 'y-.par'
write(1,'(A)') 'y+.par'
write(1,'(A)') 'z-.par'
write(1,'(A)') 'z+.par'
close(1)

END SUBROUTINE Parametrization_BOX
!-----------------------------------------------
SUBROUTINE Parametrization_TSE()
INTEGER i,j,k
INTEGER nHalf,ivt
! INTEGER, ALLOCATABLE :: knpr_inP(:),knpr_inM(:),knpr_outP(:),knpr_outM(:)
! INTEGER, ALLOCATABLE :: knpr_zMin(:),knpr_zMax(:)

nHalf = nvt/2

!! Inner cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/in+.par',UNIT=1)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/in-.par',UNIT=2)
WRITE(1,*) (NT1+NT2)*(NZ+1), ' Inflow11'
WRITE(1,'(A2,4E15.7,A)') "'7", 0e0,  0.5*DA, 0e0, 0.49*DZi, " 1.0 1.0 0.0'"
WRITE(2,*) (NT1+NT2)*(NZ+1), ' Inflow12'
WRITE(2,'(A2,4E15.7,A)') "'7", 0e0, -0.5*DA, 0e0, 0.49*DZi, " 1.0 1.0 0.0'"
k = 0
DO j=1,nZ+1
 DO i=1,NT1+NT2
  ivt = (NR+1)*(NT1+NT2)*(j-1) + i
  knpr(ivt) = 1
  knpr(ivt + nHalf) = 1
  WRITE(1,'(I6)') ivt
  WRITE(2,'(I6)') ivt + nHalf
  k = k + 1
END DO
END DO
CLOSE(1)
CLOSE(2)


!! Outer cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/out+.par',UNIT=1)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/out-.par',UNIT=2)
WRITE(1,*) (NT2+1)*(NZ+1), ' Wall'
WRITE(1,'(A2,4E15.7,A)') "'7", 0e0,  0.5*DA, 0e0, 0.5*DZo, " 1.0 1.0 0.0'"
WRITE(2,*) (NT2+1)*(NZ+1), ' Wall'
WRITE(2,'(A2,4E15.7,A)') "'7", 0e0, -0.5*DA, 0e0, 0.5*DZo, " 1.0 1.0 0.0'"
k = 0
DO j=1,nZ+1
 DO i=1,NT2
  ivt = (NR+1)*(NT1+NT2)*(j-1) + NR*(NT1+NT2) + i
  knpr(ivt) = 1
  knpr(ivt + nHalf) = 1
  WRITE(1,'(I6)') ivt
  WRITE(2,'(I6)') ivt + nHalf
  k = k + 1
END DO
ivt = (NR+1)*(NT1+NT2)*(j)
knpr(ivt) = 1
knpr(ivt + nHalf) = 1
WRITE(1,'(I6)') ivt
WRITE(2,'(I6)') ivt + nHalf
k = k + 1
END DO
CLOSE(1)
CLOSE(2)

!! Lines and straight
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/line1.par',UNIT=1)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/line2.par',UNIT=2)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/line3.par',UNIT=3)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/line4.par',UNIT=4)
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/straight.par',UNIT=5)
WRITE(1,*) (NZ+1), ' Wall'
WRITE(1,'(A)') "'1 0.0'"
WRITE(2,*) (NZ+1), ' Wall'
WRITE(2,'(A)') "'1 0.0'"
WRITE(3,*) (NZ+1), ' Wall'
WRITE(3,'(A)') "'1 0.0'"
WRITE(4,*) (NZ+1), ' Wall'
WRITE(4,'(A)') "'1 0.0'"
WRITE(5,*) 4*(NZ+1), ' Wall'
IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."STRAIGHT") THEN
 WRITE(5,'(A,ES12.4,A)') "'7 0.0 0.0 0.0 ", yCut, " 1.0 0.0 0.0'"
END IF
IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."ROUND") THEN
 WRITE(5,'(A,ES12.4,A)') "'7 0.0 0.0 0.0 ", DZz*0.5, " 1.0 1.0 0.0'"
END IF
k = 0
DO j=1,nZ+1
  ivt = (NR+1)*(NT1+NT2)*(j-1) + NR*(NT1+NT2) + NT2
  WRITE(1,'(I6)') ivt
  WRITE(5,'(I6)') ivt
  WRITE(2,'(I6)') ivt + nHalf
  WRITE(5,'(I6)') ivt + nHalf
  ivt = (NR+1)*(NT1+NT2)*(j)
  WRITE(3,'(I6)') ivt
  WRITE(5,'(I6)') ivt
  WRITE(4,'(I6)') ivt + nHalf
  WRITE(5,'(I6)') ivt + nHalf
  k = k + 1
END DO
CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)

!! Outflow
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z+.par',UNIT=1)
IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
 WRITE(1,*) 2*(NT2+NT1)*(NR+1), ' Outflow'
ELSE
 WRITE(1,*) 2*(NT2+NT1)*(NR+1), ' Symmetry110'
END IF
WRITE(1,'(A,E15.7,A)') "'4 0.0 0.0 1.0", -DL, "'"
k = 0
DO j=1,nR+1
 DO i=1,NT1+NT2
  ivt = (NR+1)*(NT1+NT2)*(nZ) + (j-1)*(NT1+NT2) + i
  knpr(ivt) = 1
  knpr(ivt + nHalf) = 1
  WRITE(1,'(I6)') ivt
  WRITE(1,'(I6)') ivt + nHalf
  k = k + 1
 END DO
END DO
CLOSE(1)

!! Inflow
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/z-.par',UNIT=2)
IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
 WRITE(2,*) 2*(NT2+NT1)*(NR+1), ' Outflow'
ELSE
 WRITE(2,*) 2*(NT2+NT1)*(NR+1), ' Inflow10'
END IF
WRITE(2,'(A,E15.7,A)') "'4 0.0 0.0 1.0", 0.0, "'"
k = 0
DO j=1,nR+1
 DO i=1,NT1+NT2
  ivt = (j-1)*(NT1+NT2) + i
  knpr(ivt) = 1
  knpr(ivt + nHalf) = 1
  WRITE(2,'(I6)') ivt
  WRITE(2,'(I6)') ivt + nHalf
  k = k + 1
 END DO
END DO
CLOSE(2)

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/file.prj',UNIT=2)
WRITE(2,'(A)') 'Mesh.tri'
WRITE(2,'(A)') 'in+.par'
WRITE(2,'(A)') 'in-.par'
WRITE(2,'(A)') 'out+.par'
WRITE(2,'(A)') 'out-.par'
WRITE(2,'(A)') 'z-.par'
WRITE(2,'(A)') 'z+.par'
WRITE(2,'(A)') 'line2.par'
WRITE(2,'(A)') 'line2.par'
WRITE(2,'(A)') 'line3.par'
WRITE(2,'(A)') 'line4.par'
WRITE(2,'(A)') 'straight.par'
CLOSE(2)

END SUBROUTINE Parametrization_TSE
!-----------------------------------------------
SUBROUTINE WriteMesh_YXZ()
INTEGER k

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/Mesh.tri',UNIT=2)
WRITE(2,*) "Coarse mesh exported by Partitioner"
WRITE(2,*) "Parametrisierung PARXC, PARYC, TMAXC"
WRITE(2,'(2(I6),A30)') nel,nvt," 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"
WRITE(2,*) "DCORVG"
DO k=1,nvt
 WRITE(2,*)dCoor(2,k),dCoor(1,k),dCoor(3,k)
END DO

WRITE(2,*) 'KVERT'
DO k=1,nel
 WRITE(2,'(8I8)') KVERT(:,k)
END DO

WRITE(2,*) 'KNPR'
DO k=1,nvt
 WRITE(2,'(I3)') knpr(k)
END DO

CLOSE(2)

END SUBROUTINE WriteMesh_YXZ
!-----------------------------------------------
SUBROUTINE WriteMesh_XYZ()
INTEGER k

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/Mesh.tri',UNIT=2)
WRITE(2,*) "Coarse mesh exported by Partitioner"
WRITE(2,*) "Parametrisierung PARXC, PARYC, TMAXC"
WRITE(2,'(2(I6),A30)') nel,nvt," 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"
WRITE(2,*) "DCORVG"
DO k=1,nvt
 WRITE(2,*)dCoor(1,k),dCoor(2,k),dCoor(3,k)
END DO

WRITE(2,*) 'KVERT'
DO k=1,nel
 WRITE(2,'(8I8)') KVERT(:,k)
END DO

WRITE(2,*) 'KNPR'
DO k=1,nvt
 WRITE(2,'(I3)') knpr(k)
END DO

CLOSE(2)

END SUBROUTINE WriteMesh_XYZ

END Program e3d_mesher
