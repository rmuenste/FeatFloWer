Program e3d_mesher
USE Sigma_User, ONLY : mySetup,mySigma,myProcess
implicit none
REAL*8 DZi,DZo,DL
INTEGER nR,nT,nZ
REAL*8,  ALLOCATABLE ::  dCoor(:,:)
INTEGER, ALLOCATABLE ::  kVert(:,:),knpr(:)
INTEGER, ALLOCATABLE :: knpr_inP(:),knpr_inM(:),knpr_outP(:),knpr_outM(:)
INTEGER, ALLOCATABLE :: knpr_zMin(:),knpr_zMax(:)
CHARACTER*(50) CaseFile,command
INTEGER NEL,NVT
character(8)  :: cdate
character(10) :: ctime
character(5)  :: czone
integer,dimension(8) :: values

LOGICAL bExist
CHARACTER cExtrud3DFile*120

call date_and_time(cdate,ctime,czone,values)

!WRITE(CaseFile,'(5A)') 'Case_',cdate(3:4),cdate(5:6),cdate(7:8),ctime(1:6)
 CaseFile='_data/meshDir'
WRITE(*,*) CaseFile

CALL ReadPar()

IF (ADJUSTL(TRIM(mySetup%cMesher)).EQ."HOLLOWCYLINDER") THEN
 WRITE(command,'(2A)') 'mkdir ',TRIM(ADJUSTL(CaseFile))
 CALL system(TRIM(ADJUSTL(command)))

 nR = mySetup%m_nR;     nT = mySetup%m_nT;    nZ = mySetup%m_nZ
 Dzo = mySigma%Dz_out;  Dzi = mySigma%Dz_in;  DL  = mySigma%L

 IF (MOD(nZ,2).EQ.1) THEN
  WRITE(*,*) "number of elements in Z is uneven and is to be set to: ", nZ+1
  nZ = nZ + 1
 END IF
 
 CALL SetUpMesh()

 CALL Parametrization()

 CALL WriteMesh()
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
SUBROUTINE SetUpMesh()
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

END SUBROUTINE SetUpMesh
!-----------------------------------------------
SUBROUTINE Parametrization()
INTEGER i,j,k
INTEGER nHalf,ivt
! INTEGER, ALLOCATABLE :: knpr_inP(:),knpr_inM(:),knpr_outP(:),knpr_outM(:)
! INTEGER, ALLOCATABLE :: knpr_zMin(:),knpr_zMax(:)

nHalf = nvt/2

!! Inner cylinder
OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/in.par',UNIT=1)
WRITE(1,*) (nT)*(NZ+1), ' Inflow11'
WRITE(1,'(A2,4E14.6,A)') "'7", 0e0,  0e0, 0e0, 0.5*DZi, " 1.0 1.0 0.0'"

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
WRITE(1,*) (nT)*(NR+1), ' Outflow'
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

END SUBROUTINE Parametrization
!-----------------------------------------------
SUBROUTINE WriteMesh()
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

OPEN(FILE=ADJUSTL(TRIM(CaseFile))//'/file.prj',UNIT=2)
WRITE(2,'(A)') 'Mesh.tri'
WRITE(2,'(A)') 'in.par'
WRITE(2,'(A)') 'out.par'
WRITE(2,'(A)') 'z-.par'
WRITE(2,'(A)') 'z+.par'
CLOSE(2)

END SUBROUTINE WriteMesh

END Program e3d_mesher
