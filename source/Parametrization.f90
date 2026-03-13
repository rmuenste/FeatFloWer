MODULE Parametrization

USE PP3D_MPI, ONLY:myid,showid,myMPI_Barrier,master,sort1d
USE var_QuadScalar

IMPLICIT NONE

TYPE tBndr
 INTEGER nVerts,nFaces
 LOGICAL, ALLOCATABLE :: Vert(:)
 INTEGER, ALLOCATABLE :: Face(:,:)
 REAL*8, ALLOCATABLE :: CoorList(:,:)
END TYPE tBndr

TYPE tParBndr
 TYPE(tBndr), ALLOCATABLE :: Bndr(:)
 CHARACTER :: Names*200,Types*200,Parameters*200
 INTEGER :: nBndrPar,CGAL_ID,Dimens
 REAL*8, ALLOCATABLE :: dBndrPar(:)
END TYPE tParBndr

INTEGER nBnds,iBnds
TYPE (tParBndr), ALLOCATABLE :: myParBndr(:)
integer, parameter :: JSON_ENTRY_MASTER = 0
integer, parameter :: JSON_ENTRY_SUB_COARSE = 1
integer, parameter :: JSON_ENTRY_SUB_PART = 2

CONTAINS
!----------------------------------------------------------------------------------------
SUBROUTINE InitBoundaryStructure(kvert,kedge)
implicit none
INTEGER kvert(8,*),kedge(12,*)
INTEGER ndof,i,j,k,ivt1,ivt2,iLoc
INTEGER iSym(3)
INTEGER Neigh(2,12)
DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
CHARACTER cAux*50,cType*50
INTEGER   iType, iPhase, jPhase, iInterface,iTemp
logical :: bOK
 
 ndof = KNVT(ilev) + KNET(ilev) + KNAT(ilev) + KNEL(ilev)
 ALLOCATE(myBoundary%bWall(ndof))
 ALLOCATE(myBoundary%bSlip(ndof))
 ALLOCATE(myBoundary%iInflow(ndof))
 ALLOCATE(myBoundary%iTemperature(ndof))
 ALLOCATE(myBoundary%iPhase(ndof))
 ALLOCATE(myBoundary%bOutflow(ndof))
 ALLOCATE(myBoundary%bSymmetry(3,ndof))
 ALLOCATE(BndrForce(ndof))
 ALLOCATE(myBoundary%LS_zero(ndof))
 ALLOCATE(myBoundary%bDisp_DBC(ndof)) 

 myBoundary%bWall     = .FALSE.
 myBoundary%bSlip     = .FALSE.
 myBoundary%iInflow   = 0
 myBoundary%iTemperature   = 0
 myBoundary%iPhase    = 0
 myBoundary%bOutflow  = .FALSE.
 myBoundary%bSymmetry = .FALSE.
 myBoundary%LS_zero   = 0
 BndrForce            = .FALSE.
 myBoundary%bDisp_DBC = .FALSE. 

 jPhase = 0
 DO iBnds=1,nBnds

  cAux = ADJUSTL(TRIM(myParBndr(iBnds)%Types))
  iType = 0
  IF (cAux(1:6).EQ.'Inflow') THEN
   READ(cAux(7:),'(I3)') iType
  END IF
  iTemp = 0
  IF (cAux(1:11).EQ.'Temperature') THEN
   READ(cAux(12:),'(I3)') iTemp
  END IF
  iPhase = 0
  IF (cAux(1:5).EQ.'Phase') THEN
   READ(cAux(6:),'(I2)') iPhase
  END IF
  iInterface = 0
  IF (cAux(1:6).EQ.'Bubble') THEN
   READ(cAux(7:),'(I2)') iInterface
  END IF  
  
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(i) = .TRUE.     
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Slip') myBoundary%bSlip(i) = .TRUE.     
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(i) = .TRUE.
      BndrForce(i) = .TRUE.
    END IF
    IF (iType.NE.0) myBoundary%iInflow(i) = iType
    IF (iTemp.GT.0) myBoundary%iTemperature(i) = iTemp
    IF (iPhase.GT.0) myBoundary%iPhase(i) = iPhase
    IF (iInterface.GT.0) myBoundary%LS_zero(i) = iInterface
!     IF (iPhase.GT.0) jPhase = jPhase + 1
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(i) = .TRUE.
    iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
    IF (ILOC.EQ.1) THEN
     READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
     IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,i) = .TRUE.
     IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,i) = .TRUE.
     IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,i) = .TRUE.
    END IF
   END IF
  END DO

  k=1
  DO i=1,NEL
   DO j=1,12
    IF (k.eq.kedge(j,i)) THEN
     ivt1 = kvert(Neigh(1,j),i)
     ivt2 = kvert(Neigh(2,j),i)
     IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt1).AND.myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt2)) THEN
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+k) = .TRUE.
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Slip') myBoundary%bSlip(nvt+k) = .TRUE.
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+k) = .TRUE.      
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
        myBoundary%bWall(nvt+k) = .TRUE.
        BndrForce(nvt+k) = .TRUE.
      END IF
      IF (iType.NE.0) myBoundary%iInflow(nvt+k) = iType
      IF (iTemp.GT.0) myBoundary%iTemperature(nvt+k) = iTemp
      IF (iPhase.GT.0) myBoundary%iPhase(nvt+k) = iPhase
      IF (iInterface.GT.0) myBoundary%LS_zero(nvt+k) = iInterface
!       IF (iPhase.GT.0) jPhase = jPhase + 1
      IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+k) = .TRUE.
      iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
      IF (ILOC.EQ.1) THEN
       READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
       IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+k) = .TRUE.
       IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+k) = .TRUE.
       IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+k) = .TRUE.
      END IF
     END IF
     k = k + 1
    END IF
   END DO
  END DO

  DO i=1,NAT
   IF (myParBndr(iBnds)%Bndr(ILEV)%Face(1,i).NE.0) THEN
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+net+i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Slip') myBoundary%bSlip(nvt+net+i) = .TRUE.
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+net+i) = .TRUE.    
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
      myBoundary%bWall(nvt+net+i) = .TRUE.
      BndrForce(nvt+net+i) = .TRUE.
    END IF
    IF (iType.NE.0) myBoundary%iInflow(nvt+net+i) = iType
    IF (iTemp.GT.0) myBoundary%iTemperature(nvt+net+i) = iTemp
    IF (iPhase.GT.0) myBoundary%iPhase(nvt+net+i) = iPhase
    IF (iInterface.GT.0) myBoundary%LS_zero(nvt+net+i) = iInterface
!     IF (iPhase.GT.0) jPhase = jPhase + 1
    IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+net+i) = .TRUE.
    iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
    IF (ILOC.EQ.1) THEN
     READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
     IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+net+i) = .TRUE.
     IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+net+i) = .TRUE.
     IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+net+i) = .TRUE.
    END IF
   END IF   
  END DO

  DO i=1,NEL
   bOK = myParBndr(iBnds)%Bndr(ILEV)%Vert(kvert(1,i))
   DO j=2,8
     bOK = (bOK.AND.myParBndr(iBnds)%Bndr(ILEV)%Vert(kvert(j,i)))
   END DO
   IF (bOK) THEN
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') myBoundary%bWall(nvt+net+nat+i) = .TRUE.
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Disp_DBC') myBoundary%bDisp_DBC(nvt+net+nat+i) = .TRUE.      
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'WallF') THEN
       myBoundary%bWall(nvt+net+nat+i) = .TRUE.
       BndrForce(nvt+net+nat+i) = .TRUE.
     END IF
     IF (iType.NE.0) myBoundary%iInflow(nvt+net+nat+i) = iType
     IF (iTemp.GT.0) myBoundary%iTemperature(nvt+net+nat+i) = iTemp
     IF (iPhase.GT.0) myBoundary%iPhase(nvt+net+nat+i) = iPhase
     IF (iPhase.GT.0) jPhase = jPhase + 1
     IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Outflow') myBoundary%bOutFlow(nvt+net+nat+i) = .TRUE.
     iLoc = INDEX(myParBndr(iBnds)%Types,'Symmetry')
     IF (ILOC.EQ.1) THEN
      READ(myParBndr(iBnds)%Types(9:11),'(3I1)') iSym
      IF (iSym(1).EQ.1) myBoundary%bSymmetry(1,nvt+net+nat+i) = .TRUE.
      IF (iSym(2).EQ.1) myBoundary%bSymmetry(2,nvt+net+nat+i) = .TRUE.
      IF (iSym(3).EQ.1) myBoundary%bSymmetry(3,nvt+net+nat+i) = .TRUE.
     END IF
   END IF
  END DO
 
 END DO

END SUBROUTINE InitBoundaryStructure
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ReviseWallBC(mgMesh,ilevel)
integer :: ilevel
type(tMultiMesh) :: mgMesh
!integer ivt,ndof
integer inode,i,j,k,iel,iat,ivt1,ivt2,ivt3,ivt4,nn
integer nn_000,nn_nvt,nn_net,nn_nat,nn_nel
REAL*8, ALLOCATABLE :: daux(:),out(:,:)
real*8 iPlus
INTEGER NeighA(4,6),NeighU(4,6)
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
DATA NeighU/1,2,3,4,1,6,9,5,2,7,10,6,3,8,11,7,4,5,12,8,9,10,11,12/


! 
!  ndof = mgMesh%level(ilevel)%nvt + &
!         mgMesh%level(ilevel)%net + &
!         mgMesh%level(ilevel)%nat + &
!         mgMesh%level(ilevel)%nel
!         
!  DO iBnds=1,nBnds
!   IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') THEN
!    DO ivt=1,ndof
!      IF (.not.mgMesh%BndryNodes(ivt)%bOuterPoint) THEN
!       myBoundary%bWall(ivt) = .FALSE.
!      END IF
!    END DO
!   END IF
!  END DO
nn_000 = 0 
nn_NVT = mgMesh%level(ilevel)%nvt 
nn_NET = mgMesh%level(ilevel)%nvt + mgMesh%level(ilevel)%net
nn_NAT = mgMesh%level(ilevel)%nvt + mgMesh%level(ilevel)%net + mgMesh%level(ilevel)%nat
nn_NEL = mgMesh%level(ilevel)%nvt + mgMesh%level(ilevel)%net + mgMesh%level(ilevel)%nat + mgMesh%level(ilevel)%nel

allocate(out(2,nn_NEL))
out = 0d0

ALLOCATE(daux(nn_NEL))
daux = 0d0

DO iBnds=1,nBnds
 IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Wall') THEN
  
  myBoundary%bWall(nn_000+1:nn_NET) = .FALSE.
  myBoundary%bWall(nn_nat+1:nn_NEL) = .FALSE.
  
  k = 1
  DO iel=1,mgMesh%level(ilevel)%nel
    DO j=1,6
      IF (k.eq.mgMesh%level(ilevel)%karea(j,iel)) THEN
        iat = nn_NET + k
        IF (myBoundary%bWall(iat).and.mgMesh%BndryNodes(iat)%bOuterPoint) THEN
!         IF (myBoundary%bWall(iat).and.mgMesh%BndryNodes(iat)%bOuterPoint.and.(.not.myBoundary%bOutFlow(iat))) THEN
          
          ! CornerPoints
          ivt1 = mgMesh%level(ilevel)%kvert(NeighA(1,j),iel)
          ivt2 = mgMesh%level(ilevel)%kvert(NeighA(2,j),iel)
          ivt3 = mgMesh%level(ilevel)%kvert(NeighA(3,j),iel)
          ivt4 = mgMesh%level(ilevel)%kvert(NeighA(4,j),iel)
          myBoundary%bWall(ivt1) = .TRUE.
          myBoundary%bWall(ivt2) = .TRUE.
          myBoundary%bWall(ivt3) = .TRUE.
          myBoundary%bWall(ivt4) = .TRUE.
          
          ! EdgePoints
          ivt1 = mgMesh%level(ilevel)%kedge(NeighU(1,j),iel)
          ivt2 = mgMesh%level(ilevel)%kedge(NeighU(2,j),iel)
          ivt3 = mgMesh%level(ilevel)%kedge(NeighU(3,j),iel)
          ivt4 = mgMesh%level(ilevel)%kedge(NeighU(4,j),iel)
          myBoundary%bWall(nn_nvt+ivt1) = .TRUE.
          myBoundary%bWall(nn_nvt+ivt2) = .TRUE.
          myBoundary%bWall(nn_nvt+ivt3) = .TRUE.
          myBoundary%bWall(nn_nvt+ivt4) = .TRUE.
          
          ! FacePoint
          myBoundary%bWall(iat) = .TRUE.
        ELSE
          myBoundary%bWall(iat) = .FALSE.
        END IF
        k = k + 1
      END IF
    END DO
  END DO
 END IF
END DO

DO iel=1,mgMesh%level(ilevel)%nel
  DO j=1,6
     k = mgMesh%level(ilevel)%karea(j,iel)
     iat = nn_NET + k
     
     iPlus = 0d0
     IF (mgMesh%BndryNodes(iat)%bOuterPoint) THEN
      IF (myBoundary%bOutflow(iat)) iPlus = 1d0
      
      ! CornerPoints
      ivt1 = mgMesh%level(ilevel)%kvert(NeighA(1,j),iel)
      ivt2 = mgMesh%level(ilevel)%kvert(NeighA(2,j),iel)
      ivt3 = mgMesh%level(ilevel)%kvert(NeighA(3,j),iel)
      ivt4 = mgMesh%level(ilevel)%kvert(NeighA(4,j),iel)
      out(:,ivt1) = out(:,ivt1) + [1d0,iPlus]
      out(:,ivt2) = out(:,ivt2) + [1d0,iPlus]
      out(:,ivt3) = out(:,ivt3) + [1d0,iPlus]
      out(:,ivt4) = out(:,ivt4) + [1d0,iPlus]
      
      ! EdgePoints
      ivt1 = mgMesh%level(ilevel)%kedge(NeighU(1,j),iel)
      ivt2 = mgMesh%level(ilevel)%kedge(NeighU(2,j),iel)
      ivt3 = mgMesh%level(ilevel)%kedge(NeighU(3,j),iel)
      ivt4 = mgMesh%level(ilevel)%kedge(NeighU(4,j),iel)
      out(:,nn_nvt+ivt1) = out(:,nn_nvt+ivt1) + [1d0,iPlus]
      out(:,nn_nvt+ivt2) = out(:,nn_nvt+ivt2) + [1d0,iPlus]
      out(:,nn_nvt+ivt3) = out(:,nn_nvt+ivt3) + [1d0,iPlus]
      out(:,nn_nvt+ivt4) = out(:,nn_nvt+ivt4) + [1d0,iPlus]
      
      out(:,iat) = out(:,iat) + [1d0,iPlus]
     END IF
  END DO
END DO

if (bParallel) then
 if (myid.ne.0) THEN
  CALL E013SUM(out(1,:))
  CALL E013SUM(out(2,:))

  daux = 0d0
  do i=1,nn_NEL
   if (myBoundary%bWall(i)) then
    daux(i) = 1d0
   end if
  end do
  CALL E013SUM(daux)
  do i=1,NN_NEL
   if (daux(i).ge.1d0) then
    myBoundary%bWall(i) = .true.
   end if
  end do
 end if
end if

do i=1,nn_NEL
 if (myBoundary%bWall(i).and.myBoundary%bOutflow(i).and.mgMesh%BndryNodes(i)%bOuterPoint) then
   if (out(1,i).eq.out(2,i).and.(out(1,i).gt.0)) THEN
!     write(*,*) myid, i,' / ', nn_NEL
    myBoundary%bWall(i) = .false.
   end if
 end if
end do

DEALLOCATE(daux,out)

END SUBROUTINE ReviseWallBC
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeQ2Nodes(dCoor)
REAL*8 dCoor(3,*)
INTEGER I1,I2

 ILEV=NLMAX
 CALL SETLEV(2)
 i1 = NVT+1
 i2 = NVT+NET+NAT

 NLMAX=NLMAX+1
 ILEV=NLMAX
 CALL SETLEV(2)

 DO iBnds = 1, nBnds
  
   CALL Parametrize(dCoor,i1,i2,ilev)

 END DO

 NLMAX=NLMAX-1
 ILEV=NLMAX
 CALL SETLEV(2)

END SUBROUTINE ParametrizeQ2Nodes
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ProlongateParametrization_STRCT(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

INTEGER I1,I2

 DO iBnds = 1, nBnds

  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%Vert)) THEN
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Vert(mgMesh%level(ilevel)%NVT))
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Face(2,mgMesh%level(ilevel)%NAT))

  myParBndr(iBnds)%Bndr(ilevel)%nVerts =0
  myParBndr(iBnds)%Bndr(ilevel)%nFaces =0
  myParBndr(iBnds)%Bndr(ilevel)%Vert =.FALSE.
  myParBndr(iBnds)%Bndr(ilevel)%Face = 0

  I1 = ilevel-1
  I2 = ilevel

  CALL GetHigherStructures(mgMesh%level(i1)%kvert,&
                           mgMesh%level(i1)%kedge,&
                           mgMesh%level(i2)%kvert,&
                           mgMesh%level(i2)%karea,&
                           mgMesh%level(i2)%kadj,&
                           mgMesh%level(i1)%nel,&
                           mgMesh%level(i1)%nvt,&
                           mgMesh%level(i1)%nat,&
                           mgMesh%level(i1)%net)
  END IF

 END DO

END SUBROUTINE ProlongateParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndr(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

INTEGER I1,I2

 DO iBnds = 1, nBnds

  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%Vert)) THEN
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Vert(mgMesh%level(ilevel)%NVT))
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%Face(2,mgMesh%level(ilevel)%NAT))

  myParBndr(iBnds)%Bndr(ilevel)%nVerts =0
  myParBndr(iBnds)%Bndr(ilevel)%nFaces =0
  myParBndr(iBnds)%Bndr(ilevel)%Vert =.FALSE.
  myParBndr(iBnds)%Bndr(ilevel)%Face = 0

  I1 = ilevel-1
  I2 = ilevel

  CALL GetHigherStructures(mgMesh%level(i1)%kvert,&
                           mgMesh%level(i1)%kedge,&
                           mgMesh%level(i2)%kvert,&
                           mgMesh%level(i2)%karea,&
                           mgMesh%level(i2)%kadj,&
                           mgMesh%level(i1)%nel,&
                           mgMesh%level(i1)%nvt,&
                           mgMesh%level(i1)%nat,&
                           mgMesh%level(i1)%net)
  END IF

  CALL Parametrize(mgmesh%level(ilevel)%dcorvg,1,mgmesh%level(ilevel)%nvt,ilevel)

 END DO

END SUBROUTINE ParametrizeBndr
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndryPoints(mgMesh,ilevel)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ilevel
type(tMultiMesh) :: mgMesh

 DO iBnds = 1, nBnds

  CALL Parametrize(mgmesh%level(ilevel)%dcorvg,1,mgmesh%level(ilevel)%nvt,ilevel)

 END DO

END SUBROUTINE ParametrizeBndryPoints
!===============================================================================================
!                             Subroutine Parametrize
!===============================================================================================
SUBROUTINE Parametrize(DCORVG,NVT1,NVT2,ilevel)
 implicit none
 real*8 DCORVG(3,*)
 integer, intent(in) :: NVT1,NVT2
 integer, intent(inout) :: ilevel

 integer :: i,j
 real*8 :: dx,dy,dz,dist,dFact
 real*8 :: px,py,pz,dScale
 real*8 :: RX,RY,RZ,RAD,DFX,DFY,DFZ
 real*8 :: DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2

 !===============================================================
 !         Handle the Bubble parametrization type (1-6)
 !===============================================================
 IF (myParBndr(iBnds)%Types(1:6).EQ.'Bubble') THEN
  IF (itns.ge.1) THEN
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i).and.myid.ne.0) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     CALL ParametrizePointOnTheBubble(dx,dy,dz,i)
     DCORVG(1,i) = dx
     DCORVG(2,i) = dy
     DCORVG(3,i) = dz
    END IF
   END DO
   RETURN
  END IF
 END IF

 !===============================================================
 !         Handle parametrization type (0)
 !===============================================================
 IF (myParBndr(iBnds)%nBndrPar.EQ.0) RETURN

 !===============================================================
 !         Handle parametrization type (1)
 !===============================================================
 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ILEV)%CoorList)) THEN
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   ALLOCATE (myParBndr(iBnds)%Bndr(ILEV)%CoorList(3,j))
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     myParBndr(iBnds)%Bndr(ILEV)%CoorList(:,j) = [dx,dy,dz]
    END IF
   END DO
  ELSE
   j = 0
   DO i=NVT1,NVT2
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
     px = myParBndr(iBnds)%Bndr(ILEV)%CoorList(1,j)
     py = myParBndr(iBnds)%Bndr(ILEV)%CoorList(2,j)
     pz = myParBndr(iBnds)%Bndr(ILEV)%CoorList(3,j)
     DCORVG(1,i) = px
     DCORVG(2,i) = py
     DCORVG(3,i) = pz
    END IF
   END DO
  END IF
 END IF

 !===============================================================
 !         Handle Plane parametrization type (4)
 !===============================================================
 ! dA, dB, dC, dD are the parameters of the general plane equation:
 ! dA * x + dB * y + dC * z + dD = 0
 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
  ! Compute the squared length of the plane normal
  ! Note: if the length of the normal is 1 then
  ! the is no difference between the DSquare and the actual length
  ! of the vector
  DSquare = dA*dA + dB*dB + dC*dC

  j = 0
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)
    dy = DCORVG(2,i)
    dz = DCORVG(3,i)

    ! First project the point (dx, dy, dz) onto
    ! the plane normal  
    dSum = (dx*dA + dy*dB + dz*dC + dD)

    ! Project the point onto the plane 
    ! by moving the point 'distance' amount
    ! opposite to the normal direction which 
    ! takes the point directly onto the plane
    px = dx - dA*dSum/dSquare
    py = dy - dB*dSum/dSquare
    pz = dz - dC*dSum/dSquare

    ! Write back the projected coordinates
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF

 !===============================================================
 !         Handle parametrization type (7)
 !===============================================================
 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)

  j = 0
  dScale = 1d0 !+ 1.0e-6
  IF (NVT1.ne.1) dScale = 1d0 !+ 0.2e-6
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
!     IF (nvt1.ne.1) write(*,*) myid,i,DCORVG(:,i),SIZE(myParBndr(iBnds)%Bndr(ILEV)%Vert)
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    dFact = dScale*RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF
 
 !===============================================================
 !  Handle varying radius cylinder parametrization type (10)
 !===============================================================
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.EQ.10) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD1 = myParBndr(iBnds)%dBndrPar(4)
  Z1 = myParBndr(iBnds)%dBndrPar(5)
  RAD2 = myParBndr(iBnds)%dBndrPar(6)
  Z2 = myParBndr(iBnds)%dBndrPar(7)
  DFX = myParBndr(iBnds)%dBndrPar(8)
  DFY = myParBndr(iBnds)%dBndrPar(9)
  DFZ = myParBndr(iBnds)%dBndrPar(10)

  j = 0
  DO i=NVT1,NVT2
   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(DCORVG(3,i)-Z1)/(Z2-Z1)
!     WRITE(*,*) dz
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF


END SUBROUTINE Parametrize
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE Parametrize_old(DCORVG,NVT,ilevel)
 REAL*8 DCORVG(3,*)
 INTEGER NVT,NVT1
 integer :: ilevel
 INTEGER i,j
 REAL*8 dx,dy,dz,dist,dFact
 REAL*8 px,py,pz
 REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
 REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,Z1,Z2

 IF (myParBndr(iBnds)%nBndrPar.EQ.0) RETURN

 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  IF (.NOT.ALLOCATED(myParBndr(iBnds)%Bndr(ilevel)%CoorList)) THEN
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   ALLOCATE (myParBndr(iBnds)%Bndr(ilevel)%CoorList(3,j))
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
     dx = DCORVG(1,i)
     dy = DCORVG(2,i)
     dz = DCORVG(3,i)
     myParBndr(iBnds)%Bndr(ilevel)%CoorList(:,j) = [dx,dy,dz]
    END IF
   END DO
  ELSE
   j = 0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
     j = j + 1
     px = myParBndr(iBnds)%Bndr(ilevel)%CoorList(1,j)
     py = myParBndr(iBnds)%Bndr(ilevel)%CoorList(2,j)
     pz = myParBndr(iBnds)%Bndr(ilevel)%CoorList(3,j)
     DCORVG(1,i) = px
     DCORVG(2,i) = py
     DCORVG(3,i) = pz
    END IF
   END DO
  END IF
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
  DSquare = dA*dA + dB*dB + dC*dC

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)
    dy = DCORVG(2,i)
    dz = DCORVG(3,i)

    dSum = (dx*dA + dy*dB + dz*dC + dD)
    px = dx - dA*dSum/dSquare
    py = dy - dB*dSum/dSquare
    pz = dz - dC*dSum/dSquare
    
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN

 IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Bubble'.AND.TIMENS.LT.1d-5) THEN
  RETURN
 END IF
 
! IF (ADJUSTL(TRIM(myParBndr(iBnds)%Types)).EQ.'Bubble'.AND.TIMENS.GT.1d-5) THEN
!  DO i=NVT1,NVT2
!   IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i).and.myid.ne.0) THEN
!    j = j + 1
!    dx = DCORVG(1,i)
!    dy = DCORVG(2,i)
!    dz = DCORVG(3,i)
!    CALL ParametrizePointOnTheBubble(dx,dy,dz,i)
!    DCORVG(1,i) = dx
!    DCORVG(2,i) = dy
!    DCORVG(3,i) = dz
!   END IF
!  END DO
!
!  RETURN
! END IF
 
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF
 
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.EQ.10) THEN

  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD1 = myParBndr(iBnds)%dBndrPar(4)
  Z1 = myParBndr(iBnds)%dBndrPar(5)
  RAD2 = myParBndr(iBnds)%dBndrPar(6)
  Z2 = myParBndr(iBnds)%dBndrPar(7)
  DFX = myParBndr(iBnds)%dBndrPar(8)
  DFY = myParBndr(iBnds)%dBndrPar(9)
  DFZ = myParBndr(iBnds)%dBndrPar(10)

  j = 0
  DO i=1,NVT
   IF (myParBndr(iBnds)%Bndr(ilevel)%Vert(i)) THEN
    j = j + 1
    dx = DCORVG(1,i)-RX
    dy = DCORVG(2,i)-RY
    dz = DCORVG(3,i)-RZ
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(DCORVG(3,i)-Z1)/(Z2-Z1)
!     WRITE(*,*) dz
    dFact = RAD/dist
    px = (1d0-DFX)*DCORVG(1,i) + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*DCORVG(2,i) + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*DCORVG(3,i) + DFZ*(RZ + dFact*dz)
    DCORVG(1,i) = px
    DCORVG(2,i) = py
    DCORVG(3,i) = pz
   END IF
  END DO
 END IF


END SUBROUTINE Parametrize_old
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE GetHigherStructures(KVERT1,KEDGE1,KVERT2,KAREA2,KADJ2,NEL1,NVT1,NAT1,NET1)
 INTEGER NAT1,NVT1,NEL1,NET1
 INTEGER KADJ2(6,*),KAREA2(6,*),KVERT2(8,*)
 INTEGER KEDGE1(12,*),KVERT1(8,*)
 INTEGER i,j,k,iel,iat,kat
 INTEGER kel(8)
 INTEGER jel(4),jat(4),jarea(4),jvt(4)
 INTEGER map(4,6),Neigh(2,12)
 DATA Neigh/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
 DATA map/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

 DO i=1,NAT1
  iel = myParBndr(iBnds)%Bndr(ILEV-1)%Face(1,i)
  iat = myParBndr(iBnds)%Bndr(ILEV-1)%Face(2,i)
  IF (iel.NE.0) THEN

   kel(1)=iel
   kel(2)=KADJ2(3,kel(1))
   kel(3)=KADJ2(3,kel(2))
   kel(4)=KADJ2(3,kel(3))
   kel(5)=KADJ2(6,kel(1))
   kel(6)=KADJ2(3,kel(5))
   kel(7)=KADJ2(3,kel(6))
   kel(8)=KADJ2(3,kel(7))

   IF (iat.EQ.1) THEN
    jel(1) = kel(1); jat(1) = 1
    jel(2) = kel(2); jat(2) = 1
    jel(3) = kel(3); jat(3) = 1
    jel(4) = kel(4); jat(4) = 1
   END IF

   IF (iat.EQ.2) THEN
    jel(1) = kel(1); jat(1) = 2
    jel(2) = kel(2); jat(2) = 5
    jel(3) = kel(6); jat(3) = 5
    jel(4) = kel(5); jat(4) = 2
   END IF

   IF (iat.EQ.3) THEN
    jel(1) = kel(2); jat(1) = 2
    jel(2) = kel(3); jat(2) = 5
    jel(3) = kel(7); jat(3) = 5
    jel(4) = kel(6); jat(4) = 2
   END IF

   IF (iat.EQ.4) THEN
    jel(1) = kel(3); jat(1) = 2
    jel(2) = kel(4); jat(2) = 5
    jel(3) = kel(8); jat(3) = 5
    jel(4) = kel(7); jat(4) = 2
   END IF

   IF (iat.EQ.5) THEN
    jel(1) = kel(4); jat(1) = 2
    jel(2) = kel(1); jat(2) = 5
    jel(3) = kel(5); jat(3) = 5
    jel(4) = kel(8); jat(4) = 2
   END IF

   IF (iat.EQ.6) THEN
    jel(1) = kel(5); jat(1) = 1
    jel(2) = kel(6); jat(2) = 1
    jel(3) = kel(7); jat(3) = 1
    jel(4) = kel(8); jat(4) = 1
   END IF

   jarea(1) = KAREA2(jat(1),jel(1))
   jarea(2) = KAREA2(jat(2),jel(2))
   jarea(3) = KAREA2(jat(3),jel(3))
   jarea(4) = KAREA2(jat(4),jel(4))

   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(1)) = [jel(1),jat(1)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(2)) = [jel(2),jat(2)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(3)) = [jel(3),jat(3)]
   myParBndr(iBnds)%Bndr(ILEV)%Face(:,jarea(4)) = [jel(4),jat(4)]

   myParBndr(iBnds)%Bndr(ILEV)%nFaces = myParBndr(iBnds)%Bndr(ILEV)%nFaces + 4

   DO kat = 1,4
    DO j=1,4
     jvt(j) = KVERT2(map(j,jat(kat)),jel(kat))
     IF (.NOT.myParBndr(iBnds)%Bndr(ILEV)%Vert(jvt(j))) THEN
      myParBndr(iBnds)%Bndr(ILEV)%nVerts = myParBndr(iBnds)%Bndr(ILEV)%nVerts + 1
      myParBndr(iBnds)%Bndr(ILEV)%Vert(jvt(j)) = .TRUE.
     END IF
    END DO
   END DO

  END IF
 END DO

 DO i=1,NVT1
  IF (myParBndr(iBnds)%Bndr(ILEV-1)%Vert(i)) myParBndr(iBnds)%Bndr(ILEV)%Vert(i) = .TRUE.
 END DO

 k=1
 DO i=1,NEL1
  DO j=1,12
   IF (k.eq.KEDGE1(j,i)) THEN
    jvt(1) = KVERT1(Neigh(1,j),i)
    jvt(2) = KVERT1(Neigh(2,j),i)
    IF (myParBndr(iBnds)%Bndr(ILEV-1)%Vert(jvt(1)).AND.myParBndr(iBnds)%Bndr(ILEV-1)%Vert(jvt(2))) THEN
     myParBndr(iBnds)%Bndr(ILEV)%Vert(NVT1+k) = .TRUE.
    END IF
    k = k + 1
   END IF
  END DO
 END DO

END SUBROUTINE GetHigherStructures
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE InitParametrization_STRCT(mesh,ilevel)
 type(tMesh) :: mesh
 integer :: ilevel
 INTEGER i,iVert,iLong,iAux,iError
 CHARACTER cFile*200,string*200
 integer :: iunit = 333
 integer :: istat,LenStr
 CHARACTER cExt*3
 integer, allocatable :: my_ParList(:,:)
 logical :: bUsedJSON

 CALL GetFileList()
 
 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Vert(&
    mesh%NVT))
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Face(2,&
    mesh%NAT))
  cFile = ADJUSTL(TRIM(cProjectFolder))
  iLong = LEN(ADJUSTL(TRIM(cFile)))+1

  LenStr = LEN(ADJUSTL(TRIM(myParBndr(iBnds)%Names)))
  READ(myParBndr(iBnds)%Names(1:LenStr-4),"(A)") string
  
  IF (INDEX(myParBndr(iBnds)%Names,'.par').ne.0) THEN
   cExt = 'par'
   IF (myid.lt.1)     THEN 
    WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".par"
   ELSE
    WRITE(cFile(iLong:),"(A,A1,A4,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".par"
   END IF
  END IF
  
  IF (INDEX(myParBndr(iBnds)%Names,'.pls').ne.0) THEN
   cExt = 'pls'
   IF (myid.lt.1)     THEN 
    WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".pls"
   ELSE
    WRITE(cFile(iLong:),"(A,A1,A4,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".pls"
   END IF
  END IF

  myParBndr(iBnds)%Bndr(1)%Vert = .FALSE.
  myParBndr(iBnds)%Bndr(1)%Face = 0
  bUsedJSON = .false.

  if (cExt.eq.'par' .and. partition_uses_json_mesh()) then
    bUsedJSON = load_boundary_from_json(trim(cFile), myParBndr(iBnds)%Bndr(ilevel)%nVerts, &
      myParBndr(iBnds)%Types, myParBndr(iBnds)%Parameters, myParBndr(iBnds)%Bndr(1)%Vert)
    if (bUsedJSON) then
      myParBndr(iBnds)%Bndr(1)%nVerts = myParBndr(iBnds)%Bndr(ilevel)%nVerts
    end if
  end if

  if (.not. bUsedJSON) then
    OPEN(UNIT = iunit, FILE = TRIM(ADJUSTL(cFile)), action='read',iostat=istat)
    if(istat .ne. 0)then
      write(*,*)"Could not open file for reading. ",TRIM(ADJUSTL(cFile))
      stop          
    end if

    READ(iunit,*)  myParBndr(iBnds)%Bndr(ilevel)%nVerts, myParBndr(iBnds)%Types
    READ(iunit,*)  myParBndr(iBnds)%Parameters

    if (cExt.eq.'par') then
     DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
      READ(iunit,*) iVert
      myParBndr(iBnds)%Bndr(1)%Vert(iVert) = .TRUE.
     END DO
    ELSE ! 'pls'
     allocate(my_ParList(4,myParBndr(iBnds)%Bndr(1)%nVerts))
     DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
      READ(iunit,*) my_ParList(1,i),my_ParList(2,i),my_ParList(3,i),my_ParList(4,i)
     END DO
    END IF

    CLOSE(iunit)
  END IF

!  WRITE(*,'(I,3A)') myid,"|",myParBndr(iBnds)%Parameters,"|"
  READ(myParBndr(iBnds)%Parameters,*,IOSTAT=iError) myParBndr(iBnds)%nBndrPar
  myParBndr(iBnds)%Dimens = 0
  IF (iError.NE.0.OR.myParBndr(iBnds)%nBndrPar.EQ.0) THEN
    myParBndr(iBnds)%nBndrPar = 0
  ELSE
    IF (myParBndr(iBnds)%nBndrPar.GT.0) THEN
      ALLOCATE(myParBndr(iBnds)%dBndrPar(myParBndr(iBnds)%nBndrPar))
      READ(myParBndr(iBnds)%Parameters,*) iAux,(myParBndr(iBnds)%dBndrPar(i),i=1,myParBndr(iBnds)%nBndrPar)
      myParBndr(iBnds)%Dimens = 3 ! ==> default:surface | will need to be redifened if we introduce specific line/point/volume parametrization ... 
    ELSE
      READ(myParBndr(iBnds)%Parameters,*) iAux,myParBndr(iBnds)%CGAL_ID
      myParBndr(iBnds)%Dimens = ABS(iAux)
      myParBndr(iBnds)%nBndrPar = 0
!       if (myid.eq.1) WRITE(*,*) myParBndr(iBnds)%Dimens,myParBndr(iBnds)%CGAL_ID
    END IF
  END IF
!  WRITE(*,'(2I8,<myParBndr(iBnds)%nBndrPar>D12.4)') myid,myParBndr(iBnds)%nBndrPar,myParBndr(iBnds)%dBndrPar

  if (cExt.eq.'par') then
   CALL GetParFaces(mesh%kvert,mesh%karea,mesh%NEL)
  else
   CALL GetPlsFaces(mesh%kvert,mesh%karea,mesh%NEL,my_ParList,myParBndr(iBnds)%Bndr(1)%nVerts)
   deallocate(my_ParList)
  end if
  
  CALL GetEdges()

 END DO

!   PAUSE

END SUBROUTINE InitParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE DeterminePointParametrization_STRCT(mgMesh,ilevel)
 integer :: ilevel
 type(tMultiMesh) :: mgMesh
 integer inode,i,j,k,iel,iat,ivt1,ivt2,ivt3,ivt4,nn
 REAL*8, ALLOCATABLE :: daux(:)
 INTEGER NeighA(4,6),NeighU(4,6)
 DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
 DATA NeighU/1,2,3,4,1,6,9,5,2,7,10,6,3,8,11,7,4,5,12,8,9,10,11,12/

 ALLOCATE(mgMesh%BndryNodes(mgMesh%level(ilevel+1)%nvt))
 ALLOCATE(daux(mgMesh%level(ilevel+1)%nvt))
 daux = 0d0
 
 DO iel=1,mgMesh%level(ilevel)%nel
  DO j=1,6
   k = mgMesh%level(ilevel)%karea(j,iel)
   daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k) = &
   daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k) + 1d0
  END DO
 END DO
 
 ILEV = ilevel
 if (bParallel) CALL E013SUM(daux)

 k = 1
 DO iel=1,mgMesh%level(ilevel)%nel
   DO j=1,6
     IF (k.eq.mgMesh%level(ilevel)%karea(j,iel)) THEN
       IF (daux(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k).eq.1d0) THEN
         ! CornerPoint
         ivt1 = mgMesh%level(ilevel)%kvert(NeighA(1,j),iel)
         ivt2 = mgMesh%level(ilevel)%kvert(NeighA(2,j),iel)
         ivt3 = mgMesh%level(ilevel)%kvert(NeighA(3,j),iel)
         ivt4 = mgMesh%level(ilevel)%kvert(NeighA(4,j),iel)
         mgMesh%BndryNodes(ivt1)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt2)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt3)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(ivt4)%bOuterPoint = .TRUE.
         ! EdgePoint
         ivt1 = mgMesh%level(ilevel)%kedge(NeighU(1,j),iel)
         ivt2 = mgMesh%level(ilevel)%kedge(NeighU(2,j),iel)
         ivt3 = mgMesh%level(ilevel)%kedge(NeighU(3,j),iel)
         ivt4 = mgMesh%level(ilevel)%kedge(NeighU(4,j),iel)
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt1)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt2)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt3)%bOuterPoint = .TRUE.
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+ivt4)%bOuterPoint = .TRUE.
         ! FacePoint
         mgMesh%BndryNodes(mgMesh%level(ilevel)%nvt+mgMesh%level(ilevel)%net+k)%bOuterPoint = .TRUE.
       ELSE
!          WRITE(*,*)
       END IF
       k = k + 1
     END IF
   END DO
 END DO

 IF (myid.ne.MASTER) then
   DO iNode=1,mgMesh%level(ilevel+1)%nvt
       mgMesh%BndryNodes(iNode)%nPoint   = 0
       mgMesh%BndryNodes(iNode)%nLine    = 0
       mgMesh%BndryNodes(iNode)%nSurf    = 0
       mgMesh%BndryNodes(iNode)%nVolume  = 0
   END DO
   DO iBnds = 1, nBnds
    DO iNode=1,mgMesh%level(ilevel+1)%nvt
      IF (myParBndr(iBnds)%Bndr(ilevel+1)%Vert(iNode)) THEN
        IF (myParBndr(iBnds)%Dimens.eq.1) mgMesh%BndryNodes(iNode)%nPoint   = mgMesh%BndryNodes(iNode)%nPoint + 1
        IF (myParBndr(iBnds)%Dimens.eq.2) mgMesh%BndryNodes(iNode)%nLine    = mgMesh%BndryNodes(iNode)%nLine + 1
        IF (myParBndr(iBnds)%Dimens.eq.3) mgMesh%BndryNodes(iNode)%nSurf    = mgMesh%BndryNodes(iNode)%nSurf + 1
        IF (myParBndr(iBnds)%Dimens.eq.4) mgMesh%BndryNodes(iNode)%nVolume  = mgMesh%BndryNodes(iNode)%nVolume + 1
      END IF
    END DO
   END DO

   DO iNode=1,mgMesh%level(ilevel+1)%nvt
     IF (mgMesh%BndryNodes(iNode)%nPoint .ne.0) THEN
      ALLOCATE(mgMesh%BndryNodes(iNode)%P(mgMesh%BndryNodes(iNode)%nPoint))
      mgMesh%BndryNodes(iNode)%P = 0
     END IF
     IF (mgMesh%BndryNodes(iNode)%nLine  .ne.0) THEN
      ALLOCATE(mgMesh%BndryNodes(iNode)%L(mgMesh%BndryNodes(iNode)%nLine))
      mgMesh%BndryNodes(iNode)%L = 0
     END IF
     IF (mgMesh%BndryNodes(iNode)%nSurf  .ne.0) THEN
      ALLOCATE(mgMesh%BndryNodes(iNode)%S(mgMesh%BndryNodes(iNode)%nSurf))
      mgMesh%BndryNodes(iNode)%S = 0
     END IF
     IF (mgMesh%BndryNodes(iNode)%nVolume.ne.0) THEN
      ALLOCATE(mgMesh%BndryNodes(iNode)%V(mgMesh%BndryNodes(iNode)%nVolume))
      mgMesh%BndryNodes(iNode)%V = 0
     END IF
   END DO

   DO iNode=1,mgMesh%level(ilevel+1)%nvt
     mgMesh%BndryNodes(iNode)%ParamTypes(1:4)=.FALSE.
     IF (mgMesh%BndryNodes(iNode)%nPoint .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(1)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nLine  .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(2)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nSurf  .ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(3)=.TRUE.
     IF (mgMesh%BndryNodes(iNode)%nVolume.ne.0) mgMesh%BndryNodes(iNode)%ParamTypes(4)=.TRUE.
   END DO
   
   mgMesh%BndryNodes(:)%nPoint  = 0
   mgMesh%BndryNodes(:)%nLine   = 0
   mgMesh%BndryNodes(:)%nSurf   = 0
   mgMesh%BndryNodes(:)%nVolume = 0
   
   DO iBnds = 1, nBnds
     DO iNode=1,mgMesh%level(ilevel+1)%nvt
       IF (myParBndr(iBnds)%Bndr(ilevel+1)%Vert(iNode)) THEN
         IF (myParBndr(iBnds)%Dimens.eq.1) THEN
           mgMesh%BndryNodes(iNode)%nPoint   = mgMesh%BndryNodes(iNode)%nPoint + 1
           nn = mgMesh%BndryNodes(iNode)%nPoint
           mgMesh%BndryNodes(iNode)%P(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.2) THEN
           mgMesh%BndryNodes(iNode)%nLine    = mgMesh%BndryNodes(iNode)%nLine + 1
           nn = mgMesh%BndryNodes(iNode)%nLine
           mgMesh%BndryNodes(iNode)%L(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.3) THEN 
           mgMesh%BndryNodes(iNode)%nSurf    = mgMesh%BndryNodes(iNode)%nSurf + 1
           nn = mgMesh%BndryNodes(iNode)%nSurf
           mgMesh%BndryNodes(iNode)%S(nn) = iBnds
         END IF
         IF (myParBndr(iBnds)%Dimens.eq.4) THEN
           mgMesh%BndryNodes(iNode)%nVolume  = mgMesh%BndryNodes(iNode)%nVolume + 1
           nn = mgMesh%BndryNodes(iNode)%nVolume
           mgMesh%BndryNodes(iNode)%V(nn) = iBnds
         END IF
       END IF
     END DO
   END DO
  
   IF (bBoundaryCheck) THEN
    DO iNode=1,mgMesh%level(ilevel+1)%nvt
      IF (mgMesh%BndryNodes(iNode)%bOuterPoint) THEN
       IF (.NOT.(mgMesh%BndryNodes(iNode)%ParamTypes(1).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(2).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(3).OR.&
                 mgMesh%BndryNodes(iNode)%ParamTypes(4))) then
           WRITE(*,'(I4,A,I8,A)') myid,'BNDRY POINT:',iNode,' is not assigned to any parametrizations ...!'
       END IF
      END IF
    END DO
   END IF

 END IF
 
 if (bParallel) then
  daux = 0d0
  do i=1,mgMesh%level(ilevel+1)%nvt
   if (mgMesh%BndryNodes(i)%bOuterPoint) then
    daux(i) = 1d0
   end if
  end do
  CALL E013SUM(daux)
  do i=1,mgMesh%level(ilevel+1)%nvt
   if (daux(i).ge.1d0) then
    mgMesh%BndryNodes(i)%bOuterPoint = .true.
   end if
  end do
 end if
 
 deallocate(daux)
 
END SUBROUTINE DeterminePointParametrization_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizeBndryPoints_STRCT(mgMesh,ii)
USE PP3D_MPI, ONLY: myid,coarse,myMPI_Barrier
IMPLICIT NONE
integer :: ii
type(tMultiMesh) :: mgMesh
REAL*8 x,y,z,cpx,cpy,cpz,d_temp,xiS(3)
INTEGER iNode,ipc,iBnds,iS,W
integer :: iRepeat,nRepeat=8

 IF (myid.eq.MASTER) return
 
 DO iNode = 1, mgmesh%level(ii)%nvt
   IF ((mgMesh%BndryNodes(iNode)%ParamTypes(1).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(2).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(3).OR.&
        mgMesh%BndryNodes(iNode)%ParamTypes(4))) then
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(1)) then
     
       GOTO 1
     END IF
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(2)) then
        iBnds = mgMesh%BndryNodes(iNode)%L(1)
        ipc   = myParBndr(iBnds)%CGAL_ID
        x = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(1,iNode)
        y = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(2,iNode)
        z = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(3,iNode)
        call projectonboundaryid(x,y,z,cpx,cpy,cpz,d_temp,ipc)
        mg_mesh%level(ii)%dcorvg(1,iNode) = cpx/dCGALtoRealFactor
        mg_mesh%level(ii)%dcorvg(2,iNode) = cpy/dCGALtoRealFactor
        mg_mesh%level(ii)%dcorvg(3,iNode) = cpz/dCGALtoRealFactor
       GOTO 1
     END IF
     IF (mgMesh%BndryNodes(iNode)%ParamTypes(3)) then
       xiS = 0d0
       W = DBLE(mgMesh%BndryNodes(iNode)%nSurf)
       DO iRepeat=1,nRepeat
       DO iS=1,mgMesh%BndryNodes(iNode)%nSurf
         iBnds = mgMesh%BndryNodes(iNode)%S(iS)
         IF (myParBndr(iBnds)%nBndrPar.NE.0) THEN !analytical
           x = mg_mesh%level(ii)%dcorvg(1,iNode)
           y = mg_mesh%level(ii)%dcorvg(2,iNode)
           z = mg_mesh%level(ii)%dcorvg(3,iNode)
           if (myParBndr(iBnds)%nBndrPar.EQ.1) then
            if (.not.mgMesh%BndryNodes(iNode)%bx) then
             mgMesh%BndryNodes(iNode)%x = [x,y,z]
             mgMesh%BndryNodes(iNode)%bx = .true.
            else 
             x = mgMesh%BndryNodes(iNode)%x(1)
             y = mgMesh%BndryNodes(iNode)%x(2)
             z = mgMesh%BndryNodes(iNode)%x(3)
            end if
           end if
           CALL projectonanalyticplane(x,y,z,cpx,cpy,cpz,iBnds)
           mg_mesh%level(ii)%dcorvg(1,iNode) = cpx
           mg_mesh%level(ii)%dcorvg(2,iNode) = cpy
           mg_mesh%level(ii)%dcorvg(3,iNode) = cpz
         ELSE  !CGAL
           ipc   = myParBndr(iBnds)%CGAL_ID
           x = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(1,iNode)
           y = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(2,iNode)
           z = dCGALtoRealFactor*mg_mesh%level(ii)%dcorvg(3,iNode)
           call projectonboundaryid(x,y,z,cpx,cpy,cpz,d_temp,ipc)
           !write(*,*)x,y,z,cpx,cpy,cpz
           mg_mesh%level(ii)%dcorvg(1,iNode) = cpx/dCGALtoRealFactor
           mg_mesh%level(ii)%dcorvg(2,iNode) = cpy/dCGALtoRealFactor
           mg_mesh%level(ii)%dcorvg(3,iNode) = cpz/dCGALtoRealFactor
         END IF
       ENDDO
       ENDDO
       GOTO 1
     END IF

   END IF
1  CONTINUE 
 END DO

END SUBROUTINE ParametrizeBndryPoints_STRCT
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE projectonanalyticplane(x1,y1,z1,x2,y2,z2,iBnds)
REAL*8 x1,y1,z1,x2,y2,z2
REAL*8 dx,dy,dz,dist,dFact
REAL*8 px,py,pz
REAL*8 RX,RY,RZ,RAD,DFX,DFY,DFZ
REAL*8 DA,DB,DC,DD,DSQUARE,DSUM,RAD1,RAD2,DZ1,DZ2
INTEGER iBnds,iSec,nSec
REAL*8, ALLOCATABLE :: Sec(:,:)

 IF (myParBndr(iBnds)%nBndrPar.EQ.0) THEN
  x2 = x1
  y2 = y1
  z2 = z1
  RETURN
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.1) THEN
  x2 = x1
  y2 = y1
  z2 = z1
  RETURN
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.4) THEN
  dA = myParBndr(iBnds)%dBndrPar(1)
  dB = myParBndr(iBnds)%dBndrPar(2)
  dC = myParBndr(iBnds)%dBndrPar(3)
  dD = myParBndr(iBnds)%dBndrPar(4)
!------------------------------------------------
  DSquare = dA*dA + dB*dB + dC*dC
  dx = x1
  dy = y1
  dz = z1
  dSum = (dx*dA + dy*dB + dz*dC + dD)
  px = dx - dA*dSum/dSquare
  py = dy - dB*dSum/dSquare
  pz = dz - dC*dSum/dSquare
  x2 = px
  y2 = py
  z2 = pz
  RETURN
 END IF

 IF (myParBndr(iBnds)%nBndrPar.EQ.7) THEN
 
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  RAD = myParBndr(iBnds)%dBndrPar(4)
  DFX = myParBndr(iBnds)%dBndrPar(5)
  DFY = myParBndr(iBnds)%dBndrPar(6)
  DFZ = myParBndr(iBnds)%dBndrPar(7)
!------------------------------------------------
  dx = x1-RX
  dy = y1-RY
  dz = z1-RZ
  dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
  dFact = RAD/dist
  px = (1d0-DFX)*x1 + DFX*(RX + dFact*dx)
  py = (1d0-DFY)*y1 + DFY*(RY + dFact*dy)
  pz = (1d0-DFZ)*z1 + DFZ*(RZ + dFact*dz)
  x2 = px
  y2 = py
  z2 = pz
  RETURN
 END IF
 
 ! Parametrization with respect to a varying radius cylinder
 IF (myParBndr(iBnds)%nBndrPar.GE.10) THEN
 
  nSec = (myParBndr(iBnds)%nBndrPar - 6)/2
  allocate(Sec(2,nSec))
  RX  = myParBndr(iBnds)%dBndrPar(1)
  RY  = myParBndr(iBnds)%dBndrPar(2)
  RZ  = myParBndr(iBnds)%dBndrPar(3)
  do iSec=1,nSec
   Sec(2,iSec) = myParBndr(iBnds)%dBndrPar(3+iSec)
  end do
  do iSec=1,nSec
   Sec(1,iSec) = myParBndr(iBnds)%dBndrPar(3+nSec+iSec)
  end do
  DFX = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+3)
  DFY = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+4)
  DFZ = myParBndr(iBnds)%dBndrPar(3+2*(nSec-1)+5)
  
!------------------------------------------------

  dx = x1-RX
  dy = y1-RY
  dz = z1-RZ
  do iSec=1,nSec-1
   RAD1= Sec(1,iSec  )
   RAD2= Sec(1,iSec+1)
   DZ1= Sec(2,iSec  )
   DZ2= Sec(2,iSec+1)
   if (dz.ge.DZ1.and.dz.le.DZ2) then
    dist = DSQRT(DFX*dx**2d0 + DFY*dy**2d0 + DFZ*dz**2d0)
    RAD = RAD1 + (RAD2-RAD1)*(z1-DZ1)/(DZ2-DZ1)
    dFact = RAD/dist
    px = (1d0-DFX)*x1 + DFX*(RX + dFact*dx)
    py = (1d0-DFY)*y1 + DFY*(RY + dFact*dy)
    pz = (1d0-DFZ)*z1 + DFZ*(RZ + dFact*dz)
   end if
  end do
  x2 = px
  y2 = py
  z2 = pz
  deallocate(Sec)
  RETURN
 END IF

END SUBROUTINE projectonanalyticplane
!
!----------------------------------------------------------------------------------------
!
 SUBROUTINE InitParametrization(mesh,ilevel)
 type(tMesh) :: mesh
 integer :: ilevel
 INTEGER i,iVert,iLong,iAux,iError
 CHARACTER cFile*200,string*200
 integer :: iunit = 333
 integer :: istat,LenStr
 CHARACTER cExt*3
 integer, allocatable :: my_ParList(:,:)
 logical :: bUsedJSON

 CALL GetFileList()

 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Vert(&
    mesh%NVT))
  ALLOCATE (myParBndr(iBnds)%Bndr(1)%Face(2,&
    mesh%NAT))
  cFile = ADJUSTL(TRIM(cProjectFolder))
  iLong = LEN(ADJUSTL(TRIM(cFile)))+1

  LenStr = LEN(ADJUSTL(TRIM(myParBndr(iBnds)%Names)))
  READ(myParBndr(iBnds)%Names(1:LenStr-4),"(A)") string
 
  IF (INDEX(myParBndr(iBnds)%Names,'.par').ne.0) THEN
   cExt = 'par'
   IF (myid.lt.1)     THEN 
    WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".par"
   ELSE
    WRITE(cFile(iLong:),"(A,A1,A4,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".par"
   END IF
  END IF
  
  IF (INDEX(myParBndr(iBnds)%Names,'.pls').ne.0) THEN
   cExt = 'pls'
   IF (myid.lt.1)     THEN 
    WRITE(cFile(iLong:),"(2A)") ADJUSTL(TRIM(string)),".pls"
   ELSE
    WRITE(cFile(iLong:),"(A,A1,A4,A)") ADJUSTL(TRIM(string)),"_",cProjectNumber,".pls"
   END IF
  END IF

!   write(*,*) ADJUSTL(TRIM(cFile))
!   pause
  myParBndr(iBnds)%Bndr(1)%Vert = .FALSE.
  myParBndr(iBnds)%Bndr(1)%Face = 0
  bUsedJSON = .false.

  if (cExt.eq.'par' .and. partition_uses_json_mesh()) then
    bUsedJSON = load_boundary_from_json(trim(cFile), myParBndr(iBnds)%Bndr(ilevel)%nVerts, &
      myParBndr(iBnds)%Types, myParBndr(iBnds)%Parameters, myParBndr(iBnds)%Bndr(1)%Vert)
    if (bUsedJSON) then
      myParBndr(iBnds)%Bndr(1)%nVerts = myParBndr(iBnds)%Bndr(ilevel)%nVerts
    end if
  end if

  if (.not. bUsedJSON) then
    OPEN(UNIT = iunit, FILE = TRIM(ADJUSTL(cFile)), action='read',iostat=istat)
    if(istat .ne. 0)then
      write(*,*)"Could not open file for reading. ",TRIM(ADJUSTL(cFile))
      stop          
    end if

    READ(iunit,*)  myParBndr(iBnds)%Bndr(ilevel)%nVerts, myParBndr(iBnds)%Types
    READ(iunit,*)  myParBndr(iBnds)%Parameters

    if (cExt.eq.'par') then
     DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
      READ(iunit,*) iVert
      myParBndr(iBnds)%Bndr(1)%Vert(iVert) = .TRUE.
     END DO
    ELSE ! 'pls'
     allocate(my_ParList(4,myParBndr(iBnds)%Bndr(1)%nVerts))
     DO i=1, myParBndr(iBnds)%Bndr(1)%nVerts
      READ(iunit,*) my_ParList(1,i),my_ParList(2,i),my_ParList(3,i),my_ParList(4,i)
     END DO
    END IF

    CLOSE(iunit)
  END IF

!  WRITE(*,'(I,3A)') myid,"|",myParBndr(iBnds)%Parameters,"|"
  READ(myParBndr(iBnds)%Parameters,*,IOSTAT=iError) myParBndr(iBnds)%nBndrPar
  IF (iError.NE.0.OR.myParBndr(iBnds)%nBndrPar.EQ.0) THEN
    myParBndr(iBnds)%nBndrPar = 0
  ELSE
    IF (myParBndr(iBnds)%nBndrPar.GT.0) THEN
      ALLOCATE(myParBndr(iBnds)%dBndrPar(myParBndr(iBnds)%nBndrPar))
      READ(myParBndr(iBnds)%Parameters,*) iAux,(myParBndr(iBnds)%dBndrPar(i),i=1,myParBndr(iBnds)%nBndrPar)
      myParBndr(iBnds)%Dimens = 3 ! ==> default:surface | will need to be redifened if we introduce specific line/point/volume parametrization ... 
    ELSE
      READ(myParBndr(iBnds)%Parameters,*) iAux,myParBndr(iBnds)%CGAL_ID
      myParBndr(iBnds)%Dimens = ABS(iAux)
    END IF
  END IF
!  WRITE(*,'(2I8,<myParBndr(iBnds)%nBndrPar>D12.4)') myid,myParBndr(iBnds)%nBndrPar,myParBndr(iBnds)%dBndrPar
!  PAUSE

  if (cExt.eq.'par') then
   CALL GetParFaces(mesh%kvert,mesh%karea,mesh%NEL)
  else
   CALL GetPlsFaces(mesh%kvert,mesh%karea,mesh%NEL,my_ParList,myParBndr(iBnds)%Bndr(1)%nVerts)
   deallocate(my_ParList)
  end if
  
  CALL GetEdges()

  CALL Parametrize(mesh%dcorvg,1,mesh%NVT,ilevel)

 END DO


END SUBROUTINE InitParametrization
!
!----------------------------------------------------------------------------------------
!
logical function load_boundary_from_json(cfile, nVerts, type_string, param_string, vertMask)
  implicit none
  character(len=*), intent(in) :: cfile
  integer, intent(out) :: nVerts
  character(len=*), intent(inout) :: type_string
  character(len=*), intent(inout) :: param_string
  logical, intent(inout) :: vertMask(:)
  character(len=:), allocatable :: json_file, sub_key, entry_key, source_name, json_text
  integer :: entry_kind
  integer :: entry_start, entry_end
  logical :: success
  integer, allocatable :: node_values(:)
  integer :: iNode, node_id

  load_boundary_from_json = .false.
  nVerts = 0

  call derive_json_lookup(trim(cfile), json_file, sub_key, entry_key, entry_kind, &
       source_name, success)
  if (.not. success) return

  call read_file_to_string(trim(json_file), json_text, success)
  if (.not. success) then
    write(*,*) "JSON boundary file not found: ", trim(json_file)
    return
  end if

  select case (entry_kind)
  case (JSON_ENTRY_MASTER)
    success = find_json_value_span(json_text, "master", entry_start, entry_end)
  case (JSON_ENTRY_SUB_COARSE)
    success = locate_sub_entry(json_text, trim(sub_key), "coarse", entry_start, entry_end)
  case (JSON_ENTRY_SUB_PART)
    success = locate_sub_part(json_text, trim(sub_key), trim(entry_key), entry_start, entry_end)
  case default
    success = .false.
  end select

  if (.not. success) then
    write(*,*) "JSON boundary entry not found: ", trim(entry_key)
    if (allocated(json_text)) deallocate(json_text)
    return
  end if

  if (.not. parse_string_field(json_text(entry_start:entry_end), "type", type_string)) then
    type_string = ""
  end if
  if (.not. parse_string_field(json_text(entry_start:entry_end), "parameter", param_string)) then
    param_string = ""
  end if

  if (.not. parse_int_list_field(json_text(entry_start:entry_end), "nodes", node_values, nVerts)) then
    if (allocated(json_text)) deallocate(json_text)
    return
  end if

  vertMask = .FALSE.
  if (nVerts > 0) then
    do iNode = 1, nVerts
      node_id = node_values(iNode)
      if (node_id >= 1 .and. node_id <= size(vertMask)) then
        vertMask(node_id) = .TRUE.
      end if
    end do
  end if

  if (allocated(node_values)) deallocate(node_values)
  if (allocated(json_text)) deallocate(json_text)
  load_boundary_from_json = .true.
end function load_boundary_from_json

logical function partition_uses_json_mesh()
  use var_QuadScalar, only: cPartitionFormat
  implicit none
  character(len=16) :: fmt
  integer :: i

  fmt = cPartitionFormat
  do i = 1, len(fmt)
    select case (fmt(i:i))
    case ("A":"Z")
      fmt(i:i) = char(iachar(fmt(i:i)) + 32)
    end select
  end do
  partition_uses_json_mesh = (trim(fmt) == "json")
end function

subroutine derive_json_lookup(cfile, json_file, sub_key, entry_key, entry_kind, source_name, success)
  implicit none
  character(len=*), intent(in) :: cfile
  character(len=:), allocatable, intent(out) :: json_file
  character(len=:), allocatable, intent(out) :: sub_key
  character(len=:), allocatable, intent(out) :: entry_key
  character(len=:), allocatable, intent(out) :: source_name
  integer, intent(out) :: entry_kind
  logical, intent(out) :: success
  character(len=:), allocatable :: dir_part, file_name
  character(len=:), allocatable :: base_dir, base_name, digits
  logical :: ok, has_digits

  json_file = ""
  sub_key = ""
  entry_key = ""
  source_name = ""
  entry_kind = JSON_ENTRY_MASTER
  success = .false.

  call split_last_component(trim(cfile), dir_part, file_name, ok)
  if (.not. ok) return

  call split_base_directory(dir_part, base_dir, sub_key)
  if (len_trim(base_dir) == 0) return

  call split_file_components(file_name, base_name, digits, has_digits, ok)
  if (.not. ok) return
  source_name = trim(base_name)

  if (len_trim(base_dir) > 0) then
    if (base_dir(len_trim(base_dir):len_trim(base_dir)) == '/' .or. &
        base_dir(len_trim(base_dir):len_trim(base_dir)) == '\') then
      json_file = trim(base_dir)//trim(source_name)//".json"
    else
      json_file = trim(base_dir)//"/"//trim(source_name)//".json"
    end if
  else
    json_file = trim(source_name)//".json"
  end if

  if (len_trim(sub_key) == 0) then
    entry_kind = JSON_ENTRY_MASTER
    entry_key = trim(source_name)
  else
    if (has_digits) then
      entry_kind = JSON_ENTRY_SUB_PART
      call build_part_entry_name(source_name, file_name, digits, entry_key)
    else
      entry_kind = JSON_ENTRY_SUB_COARSE
      entry_key = trim(source_name)
    end if
  end if

  success = (len_trim(json_file) > 0)
end subroutine

subroutine split_last_component(full_path, directory, filename, success)
  implicit none
  character(len=*), intent(in) :: full_path
  character(len=:), allocatable, intent(out) :: directory
  character(len=:), allocatable, intent(out) :: filename
  logical, intent(out) :: success
  integer :: len_path, idx
  character(len=:), allocatable :: work

  directory = ""
  filename = ""
  success = .false.
  work = trim(full_path)
  len_path = len_trim(work)
  if (len_path <= 0) return

  idx = len_path
  do while (idx >= 1)
    if (work(idx:idx) == '/' .or. work(idx:idx) == '\') exit
    idx = idx - 1
  end do
  if (idx <= 0 .or. idx == len_path) return

  directory = trim(work(:idx-1))
  filename = trim(work(idx+1:len_path))
  success = (len_trim(filename) > 0)
end subroutine

subroutine split_base_directory(path_dir, base_dir, sub_key)
  implicit none
  character(len=*), intent(in) :: path_dir
  character(len=:), allocatable, intent(out) :: base_dir
  character(len=:), allocatable, intent(out) :: sub_key
  integer :: len_dir, idx
  character(len=:), allocatable :: work, candidate

  work = trim(path_dir)
  len_dir = len_trim(work)
  if (len_dir <= 0) then
    base_dir = ""
    sub_key = ""
    return
  end if

  idx = len_dir
  do while (idx >= 1)
    if (work(idx:idx) == '/' .or. work(idx:idx) == '\') exit
    idx = idx - 1
  end do

  if (idx <= 0) then
    base_dir = trim(work)
    sub_key = ""
    return
  end if

  candidate = trim(work(idx+1:len_dir))
  if (is_subfolder_token(candidate)) then
    base_dir = trim(work(:idx-1))
    sub_key = candidate
  else
    base_dir = trim(work)
    sub_key = ""
  end if
end subroutine

logical function is_subfolder_token(token)
  implicit none
  character(len=*), intent(in) :: token
  character(len=:), allocatable :: lower
  integer :: len_tok, i

  lower = token
  call to_lower_string(lower)
  len_tok = len_trim(lower)
  if (len_tok < 4) then
    is_subfolder_token = .false.
    return
  end if
  if (lower(1:3) /= "sub") then
    is_subfolder_token = .false.
    return
  end if
  do i = 4, len_tok
    if (lower(i:i) < "0" .or. lower(i:i) > "9") then
      is_subfolder_token = .false.
      return
    end if
  end do
  is_subfolder_token = .true.
end function

subroutine split_file_components(file_name, base_name, digits, has_digits, success)
  implicit none
  character(len=*), intent(in) :: file_name
  character(len=:), allocatable, intent(out) :: base_name
  character(len=:), allocatable, intent(out) :: digits
  logical, intent(out) :: has_digits
  logical, intent(out) :: success
  character(len=:), allocatable :: work, extension, root, cleaned_root
  integer :: len_name, dot_pos, idx

  base_name = ""
  digits = ""
  has_digits = .false.
  success = .false.

  work = trim(file_name)
  len_name = len_trim(work)
  if (len_name <= 0) return

  dot_pos = len_name
  do while (dot_pos >= 1)
    if (work(dot_pos:dot_pos) == '.') exit
    dot_pos = dot_pos - 1
  end do
  if (dot_pos <= 1) return

  extension = work(dot_pos:len_name)
  root = work(:dot_pos-1)
  idx = len_trim(root)
  do while (idx >= 1)
    if (root(idx:idx) >= '0' .and. root(idx:idx) <= '9') then
      idx = idx - 1
    else
      exit
    end if
  end do

  if (idx < len_trim(root)) then
    if (len_trim(root) - idx >= 4) then
      has_digits = .true.
      digits = root(idx+1:)
    else
      has_digits = .false.
      digits = ""
    end if
  else
    has_digits = .false.
    digits = ""
  end if

  if (has_digits) then
    if (idx <= 0) then
      cleaned_root = ""
    else
      cleaned_root = root(:idx)
    end if
    cleaned_root = strip_trailing_char(cleaned_root, '_')
    if (len_trim(cleaned_root) == 0) then
      if (idx > 0) then
        cleaned_root = root(:idx)
      else
        cleaned_root = root
      end if
    end if
  else
    cleaned_root = root
  end if

  if (len_trim(cleaned_root) == 0) return
  base_name = trim(cleaned_root)//trim(extension)
  success = .true.
end subroutine

subroutine build_part_entry_name(base_name, file_name, digits, entry_name)
  implicit none
  character(len=*), intent(in) :: base_name
  character(len=*), intent(in) :: file_name
  character(len=*), intent(in) :: digits
  character(len=:), allocatable, intent(out) :: entry_name
  character(len=:), allocatable :: base_root, base_ext
  character(len=:), allocatable :: file_root, file_ext
  integer :: len_digits

  len_digits = len_trim(digits)
  if (len_digits == 0) then
    entry_name = trim(base_name)
    return
  end if

  call split_name_and_ext(base_name, base_root, base_ext)
  call split_name_and_ext(file_name, file_root, file_ext)

  if (.not. strings_equal_ignore_case(base_ext, file_ext)) then
    entry_name = trim(file_name)
    return
  end if

  if (len_trim(file_root) <= len_trim(base_root)) then
    entry_name = trim(file_name)
    return
  end if

  file_root = file_root(:len_trim(file_root)-len_digits)
  file_root = strip_trailing_char(file_root, '_')

  if (.not. strings_equal_ignore_case(file_root, base_root)) then
    entry_name = trim(file_name)
    return
  end if

  entry_name = trim(base_root)//"."//trim(digits)//trim(base_ext)
end subroutine

subroutine split_name_and_ext(full_name, name_root, name_ext)
  implicit none
  character(len=*), intent(in) :: full_name
  character(len=:), allocatable, intent(out) :: name_root
  character(len=:), allocatable, intent(out) :: name_ext
  integer :: len_name, dot_pos
  character(len=:), allocatable :: work

  work = trim(full_name)
  len_name = len_trim(work)
  dot_pos = len_name
  do while (dot_pos >= 1)
    if (work(dot_pos:dot_pos) == '.') exit
    dot_pos = dot_pos - 1
  end do
  if (dot_pos <= 0) then
    name_root = work
    name_ext = ""
    return
  end if
  name_root = work(:dot_pos-1)
  name_ext = work(dot_pos:len_name)
end subroutine

function strip_trailing_char(text, target) result(output)
  implicit none
  character(len=*), intent(in) :: text
  character, intent(in) :: target
  character(len=:), allocatable :: output
  integer :: last_idx

  output = trim(text)
  do
    last_idx = len_trim(output)
    if (last_idx <= 0) exit
    if (output(last_idx:last_idx) /= target) exit
    if (last_idx == 1) then
      output = ""
    else
      output = output(:last_idx-1)
    end if
  end do
end function

subroutine to_lower_string(text)
  implicit none
  character(len=*), intent(inout) :: text
  integer :: i

  do i = 1, len(text)
    select case (text(i:i))
    case ("A":"Z")
      text(i:i) = char(iachar(text(i:i)) + 32)
    end select
  end do
end subroutine

logical function strings_equal_ignore_case(a, b)
  implicit none
  character(len=*), intent(in) :: a
  character(len=*), intent(in) :: b
  character(len=:), allocatable :: aa, bb
  integer :: len_a, len_b, i

  aa = trim(a)
  bb = trim(b)
  len_a = len_trim(aa)
  len_b = len_trim(bb)
  if (len_a /= len_b) then
    strings_equal_ignore_case = .false.
    return
  end if
  call to_lower_string(aa)
  call to_lower_string(bb)
  do i = 1, len_a
    if (aa(i:i) /= bb(i:i)) then
      strings_equal_ignore_case = .false.
      return
    end if
  end do
  strings_equal_ignore_case = .true.
end function

subroutine read_file_to_string(filename, content, success)
  implicit none
  character(len=*), intent(in) :: filename
  character(len=:), allocatable, intent(out) :: content
  logical, intent(out) :: success
  integer(kind=8) :: file_size
  logical :: exists
  integer :: ios
  integer, parameter :: JSON_UNIT = 947
  integer :: buffer_len

  content = ""
  success = .false.
  inquire(file=trim(filename), exist=exists, size=file_size)
  if (.not. exists) return
  if (file_size <= 0) return

  buffer_len = int(file_size)
  if (buffer_len <= 0) return

  if (allocated(content)) deallocate(content)
  allocate(character(len=buffer_len) :: content)
  open(unit=JSON_UNIT, file=trim(filename), access="stream", &
       form="unformatted", status="old", action="read", iostat=ios)
  if (ios /= 0) then
    deallocate(content)
    return
  end if
  read(JSON_UNIT, iostat=ios) content
  close(JSON_UNIT)
  if (ios /= 0) then
    deallocate(content)
    return
  end if
  success = .true.
end subroutine

logical function locate_sub_entry(json_text, sub_key, target_key, start_pos, end_pos)
  implicit none
  character(len=*), intent(in) :: json_text
  character(len=*), intent(in) :: sub_key
  character(len=*), intent(in) :: target_key
  integer, intent(out) :: start_pos, end_pos
  integer :: subs_start, subs_end
  integer :: local_start, local_end
  integer :: sub_start, sub_end

  locate_sub_entry = .false.
  if (.not. find_json_value_span(json_text, "subs", subs_start, subs_end)) return
  if (.not. find_json_value_span(json_text(subs_start:subs_end), sub_key, local_start, local_end)) return
  sub_start = subs_start + local_start - 1
  sub_end = subs_start + local_end - 1
  if (.not. find_json_value_span(json_text(sub_start:sub_end), target_key, local_start, local_end)) return
  start_pos = sub_start + local_start - 1
  end_pos = sub_start + local_end - 1
  locate_sub_entry = .true.
end function

logical function locate_sub_part(json_text, sub_key, part_key, start_pos, end_pos)
  implicit none
  character(len=*), intent(in) :: json_text
  character(len=*), intent(in) :: sub_key
  character(len=*), intent(in) :: part_key
  integer, intent(out) :: start_pos, end_pos
  integer :: subs_start, subs_end
  integer :: local_start, local_end
  integer :: sub_start, sub_end
  integer :: parts_start, parts_end

  locate_sub_part = .false.
  if (.not. find_json_value_span(json_text, "subs", subs_start, subs_end)) return
  if (.not. find_json_value_span(json_text(subs_start:subs_end), sub_key, local_start, local_end)) return
  sub_start = subs_start + local_start - 1
  sub_end = subs_start + local_end - 1
  if (.not. find_json_value_span(json_text(sub_start:sub_end), "parts", local_start, local_end)) return
  parts_start = sub_start + local_start - 1
  parts_end = sub_start + local_end - 1
  if (.not. find_json_value_span(json_text(parts_start:parts_end), part_key, local_start, local_end)) return
  start_pos = parts_start + local_start - 1
  end_pos = parts_start + local_end - 1
  locate_sub_part = .true.
end function

logical function find_json_value_span(text, key, start_pos, end_pos)
  implicit none
  character(len=*), intent(in) :: text
  character(len=*), intent(in) :: key
  integer, intent(out) :: start_pos, end_pos
  character(len=:), allocatable :: token
  integer :: location, idx, len_text
  logical :: ok
  character :: ch

  find_json_value_span = .false.
  start_pos = 1
  end_pos = 0
  len_text = len(text)
  token = '"'//trim(key)//'"'
  location = index(text, token)
  if (location <= 0) return
  idx = location + len_trim(token)
  do while (idx <= len_text)
    if (text(idx:idx) == ':') exit
    idx = idx + 1
  end do
  if (idx > len_text) return
  idx = idx + 1
  idx = idx + skip_whitespace(text, idx)
  if (idx > len_text) return
  ch = text(idx:idx)
  select case (ch)
  case ('{')
    call match_block(text, idx, '{', '}', end_pos, ok)
    if (.not. ok) return
    start_pos = idx
  case ('[')
    call match_block(text, idx, '[', ']', end_pos, ok)
    if (.not. ok) return
    start_pos = idx
  case default
    call find_simple_value_end(text, idx, end_pos)
    start_pos = idx
  end select
  find_json_value_span = (end_pos >= start_pos)
end function

integer function skip_whitespace(text, start_idx)
  implicit none
  character(len=*), intent(in) :: text
  integer, intent(in) :: start_idx
  integer :: pos, len_text

  skip_whitespace = 0
  len_text = len(text)
  pos = start_idx
  do while (pos <= len_text)
    select case (text(pos:pos))
    case (' ', char(9), char(10), char(13))
      pos = pos + 1
    case default
      exit
    end select
  end do
  skip_whitespace = pos - start_idx
end function

subroutine match_block(text, start_idx, open_ch, close_ch, end_idx, success)
  implicit none
  character(len=*), intent(in) :: text
  integer, intent(in) :: start_idx
  character, intent(in) :: open_ch
  character, intent(in) :: close_ch
  integer, intent(out) :: end_idx
  logical, intent(out) :: success
  integer :: depth, i, len_text
  logical :: in_string

  success = .false.
  depth = 0
  in_string = .false.
  len_text = len(text)
  do i = start_idx, len_text
    if (text(i:i) == '"' .and. .not. is_escaped(text, i)) then
      in_string = .not. in_string
    else if (.not. in_string) then
      if (text(i:i) == open_ch) then
        depth = depth + 1
      else if (text(i:i) == close_ch) then
        depth = depth - 1
        if (depth == 0) then
          end_idx = i
          success = .true.
          return
        end if
      end if
    end if
  end do
end subroutine

logical function is_escaped(text, pos)
  implicit none
  character(len=*), intent(in) :: text
  integer, intent(in) :: pos

  if (pos <= 1) then
    is_escaped = .false.
  else
    is_escaped = (text(pos-1:pos-1) == '\')
  end if
end function

subroutine find_simple_value_end(text, start_idx, end_idx)
  implicit none
  character(len=*), intent(in) :: text
  integer, intent(in) :: start_idx
  integer, intent(out) :: end_idx
  integer :: i, len_text
  logical :: in_string

  len_text = len(text)
  in_string = .false.
  end_idx = len_text
  do i = start_idx, len_text
    if (text(i:i) == '"' .and. .not. is_escaped(text, i)) then
      in_string = .not. in_string
    else if (.not. in_string) then
      if (text(i:i) == ',' .or. text(i:i) == '}' .or. text(i:i) == ']') then
        end_idx = i - 1
        exit
      end if
    end if
  end do
  do while (end_idx >= start_idx)
    if (text(end_idx:end_idx) == ' ' .or. text(end_idx:end_idx) == char(9) .or. &
        text(end_idx:end_idx) == char(10) .or. text(end_idx:end_idx) == char(13)) then
      end_idx = end_idx - 1
    else
      exit
    end if
  end do
end subroutine

logical function parse_string_field(text, key, value)
  implicit none
  character(len=*), intent(in) :: text
  character(len=*), intent(in) :: key
  character(len=*), intent(out) :: value
  integer :: s, e
  character(len=:), allocatable :: buffer

  parse_string_field = .false.
  value = ""
  if (.not. find_json_value_span(text, key, s, e)) return
  buffer = text(s:e)
  if (.not. decode_json_string(buffer, value)) return
  parse_string_field = .true.
end function

logical function decode_json_string(raw_text, dest)
  implicit none
  character(len=*), intent(in) :: raw_text
  character(len=*), intent(out) :: dest
  integer :: start_idx, end_idx, len_raw, i, out_idx
  character :: ch

  decode_json_string = .false.
  dest = ""
  len_raw = len_trim(raw_text)
  start_idx = 1
  do while (start_idx <= len_raw .and. is_json_whitespace(raw_text(start_idx:start_idx)))
    start_idx = start_idx + 1
  end do
  end_idx = len_raw
  do while (end_idx >= start_idx .and. is_json_whitespace(raw_text(end_idx:end_idx)))
    end_idx = end_idx - 1
  end do
  if (start_idx > end_idx) return
  if (raw_text(start_idx:start_idx) /= '"' .or. raw_text(end_idx:end_idx) /= '"') return
  out_idx = 1
  i = start_idx + 1
  do while (i <= end_idx - 1)
    ch = raw_text(i:i)
    if (ch == '\' .and. i < end_idx - 1) then
      i = i + 1
      ch = raw_text(i:i)
    end if
    if (out_idx <= len(dest)) dest(out_idx:out_idx) = ch
    out_idx = out_idx + 1
    i = i + 1
  end do
  if (out_idx <= len(dest)) dest(out_idx:) = ' '
  decode_json_string = .true.
end function

logical function is_json_whitespace(ch)
  implicit none
  character, intent(in) :: ch

  select case (ch)
  case (' ', char(9), char(10), char(13))
    is_json_whitespace = .true.
  case default
    is_json_whitespace = .false.
  end select
end function

integer function count_numeric_tokens(buffer)
  implicit none
  character(len=*), intent(in) :: buffer
  integer :: i
  logical :: in_token
  character :: ch

  count_numeric_tokens = 0
  in_token = .false.
  do i = 1, len(buffer)
    ch = buffer(i:i)
    select case (ch)
    case (' ', char(9), char(10), char(13))
      in_token = .false.
    case default
      if (.not. in_token) then
        count_numeric_tokens = count_numeric_tokens + 1
        in_token = .true.
      end if
    end select
  end do
end function

logical function parse_int_list_field(text, key, values, count)
  implicit none
  character(len=*), intent(in) :: text
  character(len=*), intent(in) :: key
  integer, allocatable, intent(out) :: values(:)
  integer, intent(out) :: count
  integer :: s, e, ios
  character(len=:), allocatable :: buffer

  parse_int_list_field = .false.
  if (.not. find_json_value_span(text, key, s, e)) return
  buffer = text(s:e)
  call sanitize_number_buffer(buffer)
  count = count_numeric_tokens(buffer)
  if (count < 0) count = 0
  if (allocated(values)) deallocate(values)
  if (count == 0) then
    allocate(values(0))
    parse_int_list_field = .true.
    return
  end if
  allocate(values(count))
  read(buffer, *, iostat=ios) values
  if (ios /= 0) then
    deallocate(values)
    return
  end if
  parse_int_list_field = .true.
end function

subroutine sanitize_number_buffer(buffer)
  implicit none
  character(len=*), intent(inout) :: buffer
  integer :: i

  do i = 1, len(buffer)
    select case (buffer(i:i))
    case ('[',']',',')
      buffer(i:i) = ' '
    case default
      if (iachar(buffer(i:i)) < 32) buffer(i:i) = ' '
    end select
  end do
end subroutine
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE GetFileList()
 CHARACTER(LEN=200) :: string,cFile
 CHARACTER cWD*200,cSub*200
 INTEGER lenCommand,i,iPos,LenStr,iEnd

 if (myid.eq.1.or.(.not.bParallel))  Write(*,*)'Project File: "'//ADJUSTL(TRIM(cProjectFile))//'"'
 
 nBnds = 0
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cProjectFile)))
 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") THEN
    nBnds = nBnds + 1
   END IF
   IF (cFile(LenStr-3:LenStr).EQ.".pls") THEN
    nBnds = nBnds + 1
   END IF
  END IF
 END DO

 if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Number of boundary parametrizations: ',nBnds
 
 ALLOCATE (myParBndr(nBnds))

 if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Warning: Allocated an additional boundary level'
 DO iBnds = 1, nBnds
  ALLOCATE (myParBndr(iBnds)%Bndr(NLMIN:NLMAX+1))
 END DO

 nBnds = 0
 REWIND(1)

 DO
  READ(1,FMT='(200A)',IOSTAT=iEnd) string
  IF (iEnd.EQ.-1) EXIT
  LenStr = LEN(ADJUSTL(TRIM(string)))
  IF (LenStr.gt.4) THEN
   cFile = ADJUSTL(TRIM(string))
   IF (cFile(LenStr-3:LenStr).EQ.".par") THEN
    if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Boundary parametrization: "'//ADJUSTL(TRIM(cFile))//'"'
    nBnds = nBnds + 1
    READ(cFile(1:LenStr),"(A)") myParBndr(nBnds)%Names
   END IF
   IF (cFile(LenStr-3:LenStr).EQ.".pls") THEN
    if (myid.eq.1.or.(.not.bParallel)) Write(*,*)'Boundary parametrization: "'//ADJUSTL(TRIM(cFile))//'"'
    nBnds = nBnds + 1
    READ(cFile(1:LenStr),"(A)") myParBndr(nBnds)%Names
   END IF
  END IF
 END DO
 CLOSE(1)

END SUBROUTINE GetFileList
!----------------------------------------------------------------------------------------
SUBROUTINE GetEdges()

END SUBROUTINE GetEdges
!----------------------------------------------------------------------------------------
SUBROUTINE GetParFaces(KVERT,KAREA,NEL)
 INTEGER KVERT(8,*),KAREA(6,*),NEL
 INTEGER i,iel,iat,iarea
 INTEGER IVT1,IVT2,IVT3,IVT4
 INTEGER map(4,6)
 DATA map/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

 myParBndr(iBnds)%Bndr(ILEV)%nFaces = 0

 DO iel=1,NEL
  DO iat=1,6
   
   ivt1 = KVERT(map(1,iat),iel)
   ivt2 = KVERT(map(2,iat),iel)
   ivt3 = KVERT(map(3,iat),iel)
   ivt4 = KVERT(map(4,iat),iel)

   IF (myParBndr(iBnds)%Bndr(1)%Vert(ivt1).AND.myParBndr(iBnds)%Bndr(1)%Vert(ivt2).AND.&
       myParBndr(iBnds)%Bndr(1)%Vert(ivt3).AND.myParBndr(iBnds)%Bndr(1)%Vert(ivt4)) THEN
       myParBndr(iBnds)%Bndr(1)%nFaces = myParBndr(iBnds)%Bndr(1)%nFaces+ 1
       iarea = KAREA(iat,iel)
       myParBndr(iBnds)%Bndr(1)%Face(1,iarea) = iel
       myParBndr(iBnds)%Bndr(1)%Face(2,iarea) = iat
   END IF

  END DO
 END DO
 
END SUBROUTINE GetParFaces
!----------------------------------------------------------------------------------------
SUBROUTINE GetPlsFaces(KVERT,KAREA,NEL,l,ln)
 INTEGER KVERT(8,*),KAREA(6,*),NEL,l(4,*),ln
 INTEGER i,iel,iat,iarea
 INTEGER iFace(4),jFace(4)
 INTEGER map(4,6),lnn,j,k
 DATA map/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

 myParBndr(iBnds)%Bndr(ILEV)%nFaces = 0
 lnn = 0
 
 DO j=1,ln
 
  jFace = l(:,j)
  call sort1d(jFace,4)
  
  DO iel=1,NEL
   DO iat=1,6
    
    iFace = KVERT(map(:,iat),iel)
    call sort1d(iFace,4)

    if (iFace(1).eq.jFace(1).and.iFace(2).eq.jFace(2).and.&
        iFace(3).eq.jFace(3).and.iFace(4).eq.jFace(4)) then
        myParBndr(iBnds)%Bndr(1)%nFaces = myParBndr(iBnds)%Bndr(1)%nFaces+ 1
        iarea = KAREA(iat,iel)
        myParBndr(iBnds)%Bndr(1)%Face(1,iarea) = iel
        myParBndr(iBnds)%Bndr(1)%Face(2,iarea) = iat
        do k=1,4
         if (.not.myParBndr(iBnds)%Bndr(1)%Vert(iFace(k))) then
          myParBndr(iBnds)%Bndr(1)%Vert(iFace(k)) = .True.
          lnn = lnn + 1
         end if
        end do
    END IF

   END DO
  END DO
 END DO
 
 ln = lnn

END SUBROUTINE GetPlsFaces
!
!----------------------------------------------------------------------------------------
!
SUBROUTINE ParametrizePointOnTheBubble(px,py,pz,iPoint)
REAL*8  :: px,py,pz

INTEGER i,j,k,kk,jClose,iTriang,iDef,iPoint
INTEGER nXX,iAux,iMinDist
REAL*8 cAux(3),dAux,P(3),ProjP(3),dTriang(3,3),dTestVector(4),AX2,AY2,AZ2,AX3,AY3,AZ3,ProjPdist,dist
REAL*8 tx,ty,tz
PARAMETER (nXX = 8)
REAL*8  :: dMonitor(nXX+1),cMonitor(3,nXX+1)
INTEGER :: iMonitor(nXX+1),jSet(nXX)
!!!!!!!!!!!!!!!!!!!!!!!!!!1
REAL*8 PP(3,9),dBAS(9),T(3),S(3),R(3),PPP(3,3), dMin_dL,dL,distP,distE
INTEGER lineFace(3,4)
LOGICAL bFound,bSFound
DATA LineFace/1,3,2, 3,5,4, 5,7,6, 7,1,8/
!!!1
INTEGER iIF

iIF = myBoundary%LS_zero(iPoint)

P = [px,py,pz]

dMonitor = 1d30
iMonitor = 0

DO i=1,myTSurf(iIF)%nT
 dist = (myTSurf(iIF)%T(i)%C(1,9)-P(1))**2d0 + (myTSurf(iIF)%T(i)%C(2,9)-P(2))**2d0 + (myTSurf(iIF)%T(i)%C(3,9)-P(3))**2d0
 IF (dist.lt.dMonitor(nXX)) THEN
  dMonitor(nXX) = dist
  iMonitor(nXX) = i
  DO kk=nXX-1,1,-1
   IF (dMonitor(kk+1).lt.dMonitor(kk)) THEN
    iAux = iMonitor(kk)
    dAux = dMonitor(kk)
    iMonitor(kk) = iMonitor(kk+1)
    dMonitor(kk) = dMonitor(kk+1)
    iMonitor(kk+1) = iAux
    dMonitor(kk+1) = dAux
   END IF
  END DO
 END IF
END DO

bFound = .FALSE.

DO i = 1,nXX

iMinDist = iMonitor(i)

PP(:,1) = myTSurf(iIF)%T(iMinDist)%C(:,1)
PP(:,2) = myTSurf(iIF)%T(iMinDist)%C(:,3)
PP(:,3) = myTSurf(iIF)%T(iMinDist)%C(:,5)
PP(:,4) = myTSurf(iIF)%T(iMinDist)%C(:,7)
PP(:,5) = myTSurf(iIF)%T(iMinDist)%C(:,2)
PP(:,6) = myTSurf(iIF)%T(iMinDist)%C(:,4)
PP(:,7) = myTSurf(iIF)%T(iMinDist)%C(:,6)
PP(:,8) = myTSurf(iIF)%T(iMinDist)%C(:,8)
PP(:,9) = myTSurf(iIF)%T(iMinDist)%C(:,9)

CALL ProjectPointOntoInterphaseSurf(P,T,PP)

IF (bFound) THEN
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  A ----'
 EXIT
END IF

END DO

IF (.NOT.bFound) THEN

 dMin_dL = 1d30

 DO i = 1,nXX

  iMinDist = iMonitor(i)

  DO j=1,4
   PPP(:,1) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(1,j))
   PPP(:,2) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(2,j))
   PPP(:,3) = myTSurf(iIF)%T(iMinDist)%C(:,LineFace(3,j))
   bSFound = .FALSE.
   CALL ProjectPointOntoInterphaseLine(P,T,PPP,dL)	
   IF (bSFound) THEN 
    IF (dMin_dL.gt.dL) THEN
      bFound = .TRUE.
      dMin_dL = dL
      S = T
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  E ----'
    END IF
   END IF
  END DO
 END DO

 IF (bFound) THEN
  T = S
  distE = dMin_dL
 END IF
END IF

! IF (.NOT.bFound) THEN
!  write(*,*) 'funny things are happenning .. '
!  dMin_dL = 1d30
 bSFound = .FALSE.
 DO i = 1,nXX
  iMinDist = iMonitor(i)
  DO j=1,9
   R = myTSurf(iIF)%T(iMinDist)%C(:,j)
   dL = (R(1)-P(1))**2d0+(R(2)-P(2))**2d0+(R(3)-P(3))**2d0
   IF (dMin_dL.gt.dL) THEN
     dMin_dL = dL
     S = R
     bSFound = .TRUE.
     distP = dMin_dL
!      DIST = SQRT((px-t(1))**2d0+(py-t(2))**2d0+(pz-t(3))**2d0)
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
!      if (myid.eq.10.and.iPoint.eq.394) write(*,*) '----  P ----'
   END IF
  END DO
 END DO

 IF (distP.lt.distE.and.bSFound) THEN
  T = S
 END IF
! END IF

tx = T(1)
ty = T(2)
tz = T(3)

DIST = SQRT((px-tx)**2d0+(py-ty)**2d0+(pz-tz)**2d0)

px = Tx
py = Ty
pz = Tz

! if (myid.eq.10.and.iPoint.eq.394) write(*,*) P
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) T
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) DIST
! if (myid.eq.10.and.iPoint.eq.394) write(*,*) '==== ===== ===== ====== ====== ====== ===== ===== '

 CONTAINS

 SUBROUTINE ProjectPointOntoInterphaseSurf(P,T,C)
 REAL*8 pP(2),P(3),T(3),pPx,pPy
 REAL*8 C(3,9),Der(2),Val(3),ddd
 REAL*8 XXX1,XXX2,XXX3,YYY1,YYY2,YYY3,df1,df2,dg1,dg2
 REAL*8 updateX,updateY
 INTEGER :: nW=10,iW
 REAL*8 dEPS,dBoxEPS
 
 dEPS = 1d-5
 dBoxEPS = 1d-5
 pP = 0d0
 
 CALL Get2DQ2basis(pP(1),pP(2),P,T,C,Val(1))
!  WRITE(*,'(3ES12.4)') pP,Val(1) 
 DO iW = 1,nW

  XXX1 = pP(1) - dEPS
  XXX2 = pP(1) 
  XXX3 = pP(1) + dEPS
  YYY1 = pP(2) - dEPS
  YYY2 = pP(2) 
  YYY3 = pP(2) + dEPS
  
!   CALL Get2DQ2basis(pPx,pP(2),P,T,C,Val(2))
  CALL Get2DQ2basis(XXX1,YYY2,P,T,C,Val(1))
  CALL Get2DQ2basis(XXX2,YYY2,P,T,C,Val(2))
  CALL Get2DQ2basis(XXX3,YYY2,P,T,C,Val(3))
  df1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  df2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  CALL Get2DQ2basis(XXX2,YYY1,P,T,C,Val(1))
  CALL Get2DQ2basis(XXX2,YYY3,P,T,C,Val(3))
  dg1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  dg2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  
  updateX = 1.0d0*df1/df2
  updateY = 1.0d0*dg1/dg2
  IF (df2.eq.0d0) updateX = 0d0
  IF (dg2.eq.0d0) updateY = 0d0

  ddd = SQRT(updateX**2d0 +updateY**2d0)
  
!   IF (myid.eq.1.and.iPoint.eq.32) WRITE(*,'(I,10ES12.4)') iW,Val(2),updateX,updateY,pP,dg1,dg2

  if (ddd.lt.1d-3) DEPS = 1d-6
  if (ddd.lt.1d-4) DEPS = 1d-7
  if (ddd.lt.1d-5) DEPS = 1d-8
  if (ddd.lt.1d-6) DEPS = 1d-9
    
  pP(1) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(1) - updateX))
  pP(2) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(2) - updateY))
  
  if (ddd.lt.1d-5) THEN
   bFound = .TRUE.
   EXIT
  END IF

 END DO

 IF (bFound) THEN
  XXX2 = pP(1) 
  YYY2 = pP(2) 
  CALL Get2DQ2basis(XXX2,YYY2,P,T,C,Val(2))
 END IF
 
END SUBROUTINE ProjectPointOntoInterphaseSurf

SUBROUTINE Get2DQ2basis(X1,X2,P,T,C,DIST)
 REAL*8 DHELP(9),X1,X2,P(3),C(3,9),T(3),DIST
 REAL*8 Q4,Q2
 PARAMETER (Q4=.25D0,Q2=.5D0)
 integer j

 DHELP(1)= Q4*(1D0-X1)*(1D0-X2)*X1*X2
 DHELP(2)=-Q4*(1D0+X1)*(1D0-X2)*X1*X2
 DHELP(3)= Q4*(1D0+X1)*(1D0+X2)*X1*X2
 DHELP(4)=-Q4*(1D0-X1)*(1D0+X2)*X1*X2
 DHELP(5)=-Q2*(1D0-X1*X1)*(1D0-X2)*X2
 DHELP(6)= Q2*(1D0+X1)*(1D0-X2*X2)*X1
 DHELP(7)= Q2*(1D0-X1*X1)*(1D0+X2)*X2
 DHELP(8)=-Q2*(1D0-X1)*(1D0-X2*X2)*X1
 DHELP(9)= (1D0-X1*X1)*(1D0-X2*X2)

 T = 0d0
 DO j=1,9
  T(1) = T(1) + DHELP(j)*C(1,j)
  T(2) = T(2) + DHELP(j)*C(2,j)
  T(3) = T(3) + DHELP(j)*C(3,j)
 END DO

 DIST = (T(1)-P(1))**2d0+(T(2)-P(2))**2d0+(T(3)-P(3))**2d0
 
END SUBROUTINE Get2DQ2basis

SUBROUTINE ProjectPointOntoInterphaseLine(P,T,C,dLineDIST)
 REAL*8 pP(1),P(3),T(3),pPx,pPy,dLineDIST
 REAL*8 C(3,3),Der(1),Val(3),ddd,daux
 REAL*8 XXX1,XXX2,XXX3,df1,df2
 REAL*8 updateX
 INTEGER :: nW=10,iW
 REAL*8 dEPS,dBoxEPS
 
 dEPS = 1d-5
 dBoxEPS = 1d-5
 
 pP = 0d0
 
 CALL Get1DQ2basis(pP(1),P,T,C,Val(1))

 DO iW = 1,nW

  XXX1 = pP(1) - dEPS
  XXX2 = pP(1) 
  XXX3 = pP(1) + dEPS
  
  CALL Get1DQ2basis(XXX1,P,T,C,Val(1))
  CALL Get1DQ2basis(XXX2,P,T,C,Val(2))
  CALL Get1DQ2basis(XXX3,P,T,C,Val(3))
  df1 = (Val(3)              - Val(1)) / (2d0*DEPS)
  df2 = (Val(3) - 2d0*Val(2) + Val(1)) / (DEPS*DEPS)
  
  updateX = 1.0d0*df1/df2
  IF (df2.eq.0d0) updateX = 0d0

  ddd = ABS(updateX)
  
!   IF (myid.eq.1.and.iPoint.eq.32) WRITE(*,'(I,10ES12.4)') iW,Val(2),updateX,updateY,pP,dg1,dg2

  if (ddd.lt.1d-3) DEPS = 1d-6
  if (ddd.lt.1d-4) DEPS = 1d-7
  if (ddd.lt.1d-5) DEPS = 1d-8
  if (ddd.lt.1d-6) DEPS = 1d-9
!   if (ddd.lt.1d-3) DEPS = 0.1*deps
    
  pP(1) =   MAX(-1d0-dBoxEPS,MIN(1d0+dBoxEPS,pP(1) - updateX))
  
  if (ddd.lt.1d-5) THEN
!    WRITE(*,*) pP,iW,i
   bSFound = .TRUE.
   EXIT
  END IF

 END DO
 
 IF (bSFound) THEN
  XXX2 = pP(1)
  CALL Get1DQ2basis(XXX2,P,T,C,dLineDIST)
 END IF

END SUBROUTINE ProjectPointOntoInterphaseLine

SUBROUTINE Get1DQ2basis(X1,P,T,C,DIST)
 REAL*8 DHELP(3),X1,P(3),C(3,3),T(3),DIST
 REAL*8 Q4,Q2
 PARAMETER (Q4=.25D0,Q2=.5D0)
 integer j

 DHELP(1)= -Q2*(1D0-X1)*X1
 DHELP(2)=  Q2*(1D0+X1)*X1
 DHELP(3)= (1D0-X1)*(1d0+X1)

 T = 0d0
 DO j=1,3
  T(1) = T(1) + DHELP(j)*C(1,j)
  T(2) = T(2) + DHELP(j)*C(2,j)
  T(3) = T(3) + DHELP(j)*C(3,j)
 END DO

 DIST = (T(1)-P(1))**2d0+(T(2)-P(2))**2d0+(T(3)-P(3))**2d0
 
END SUBROUTINE Get1DQ2basis

END SUBROUTINE ParametrizePointOnTheBubble

END MODULE Parametrization
