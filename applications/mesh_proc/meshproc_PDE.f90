MODULE MeshProcPDE

USE def_LinScalar
USE types, ONLY: tMultiMesh,tMesh
USE var_QuadScalar, ONLY : UMF_CMat,UMF_lMat
USE PP3D_MPI, ONLY:myid,master
USE UMFPackSolver, ONLY : myUmfPack_Factorize,myUmfPack_Solve
USE Parametrization, ONLY: ParametrizeBndryPoints_STRCT

TYPE(lScalar3) MeshDef
REAL*8, ALLOCATABLE :: OrigCoor(:,:)
Logical :: bDefTensor = .false.
contains

SUBROUTINE MeshDefPDE()

ILEV = mg_Mesh%nlmax

do i=1,MeshDef%ndof
 if (.not.mg_Mesh%BndryNodes(i)%bOuterPoint) then
  mg_mesh%level(ilev)%dcorvg(1,i) = OrigCoor(1,i)
  mg_mesh%level(ilev)%dcorvg(2,i) = OrigCoor(2,i)
  mg_mesh%level(ilev)%dcorvg(3,i) = OrigCoor(3,i)
 end if
END DO

CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)

MeshDef%valX = mg_mesh%level(ilev)%dcorvg(1,1:MeshDef%ndof)
MeshDef%valY = mg_mesh%level(ilev)%dcorvg(2,1:MeshDef%ndof)
MeshDef%valZ = mg_mesh%level(ilev)%dcorvg(3,1:MeshDef%ndof)

! Impose Dirichlet BC for all boundary points
WRITE(*,*) 'Set Dirichlet BC'
CALL LinSc_Knpr_MeshDef()

if (bDefTensor) then
 CALL Boundary_MeshDef_MatS(S11Mat,S22Mat,S33Mat,&
                            S12Mat,S13Mat,S23Mat,&
                            S21Mat,S31Mat,S32Mat,&
                            lMat%LdA)
ELSE
 CALL Boundary_MeshDef_Mat(DMat,lMat%LdA)
END IF

MeshDef%defX = 0d0
MeshDef%defY = 0d0
MeshDef%defZ = 0d0


if (bDefTensor) then
 CALL MatDef_MeshDefS()
else
 CALL MatDef_MeshDef()
end if

CALL Boundary_MeshDef_Def()

CALL GetRes_MeshDef()

! Initialize Umfack / factorization of the matrix
WRITE(*,*) 'Init_UmfPack'
if (bDefTensor) then
 CALL Init_UmfPackS()
ELSE
 CALL Init_UmfPack()
END IF

if (bDefTensor) then
 CALL Solve_UmfPackS()
ELSE
 CALL Solve_UmfPack()
END IF

mg_mesh%level(ilev)%dcorvg(1,1:MeshDef%ndof) = MeshDef%valX
mg_mesh%level(ilev)%dcorvg(2,1:MeshDef%ndof) = MeshDef%valY
mg_mesh%level(ilev)%dcorvg(3,1:MeshDef%ndof) = MeshDef%valZ
                                            
! CALL MatDef_MeshDefS()
! 
! CALL Boundary_MeshDef_Def()
! 
! CALL GetRes_MeshDef()

END SUBROUTINE MeshDefPDE          
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr_MeshDef()
INTEGER i

DO i=1,MeshDef%ndof

 if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
  MeshDef%knprX(i) = 1
  MeshDef%knprY(i) = 1
  MeshDef%knprZ(i) = 1
 end if
 
END DO

END SUBROUTINE LinSc_Knpr_MeshDef
!
! ----------------------------------------------
!
SUBROUTINE Init_MeshDefField()
 
MeshDef%ndof = NVT
ALLOCATE(MeshDef%knprX(nvt),MeshDef%knprY(nvt),MeshDef%knprZ(nvt))
ALLOCATE(MeshDef%valX_old(nvt),MeshDef%valY_old(nvt),MeshDef%valZ_old(nvt))
ALLOCATE(MeshDef%valX(nvt),MeshDef%valY(nvt),MeshDef%valZ(nvt))
ALLOCATE(MeshDef%defX(nvt),MeshDef%defY(nvt),MeshDef%defZ(nvt))

END SUBROUTINE Init_MeshDefField
!
! ----------------------------------------------
!
SUBROUTINE Init_UmfPackS()

IF (.not.ALLOCATED(UMF_CMat)) ALLOCATE (UMF_CMat(9*lMat%na))
IF (.not.ALLOCATED(UMF_lMat%ColA)) ALLOCATE (UMF_lMat%ColA(9*lMat%na))
IF (.not.ALLOCATED(UMF_lMat%LdA)) ALLOCATE (UMF_lMat%LdA(3*lMat%nu+1))
UMF_lMat%nu   = 3*lMat%nu
UMF_lMat%na   = 9*lMat%na

k = 0
UMF_lMat%LdA(1) = 1

do i=1,lMat%nu
 kk = 0
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S11Mat(j)
  UMF_lMat%ColA(k) = 0*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k) = S12Mat(j)
  UMF_lMat%ColA(k) = 1*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k) = S13Mat(j)
  UMF_lMat%ColA(k) = 2*lMat%nu + lMat%ColA(j)
 end do
 UMF_lMat%LdA(0*lMat%nu + i+1) = UMF_lMat%LdA(0*lMat%nu + i) + kk 
end do

do i=1,lMat%nu
 kk = 0
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S21Mat(j)
  UMF_lMat%ColA(k) = 0*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S22Mat(j)
  UMF_lMat%ColA(k) = 1*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S23Mat(j)
  UMF_lMat%ColA(k) = 2*lMat%nu + lMat%ColA(j)
 end do
 UMF_lMat%LdA(1*lMat%nu + i+1) = UMF_lMat%LdA(1*lMat%nu + i) + kk 
end do

do i=1,lMat%nu
 kk = 0
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S31Mat(j)
  UMF_lMat%ColA(k) = 0*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S32Mat(j)
  UMF_lMat%ColA(k) = 1*lMat%nu + lMat%ColA(j)
 end do
 do j=lMat%LdA(i),lMat%LdA(i+1)-1
  k  = k  + 1
  kk = kk + 1
  UMF_CMat(k)      = S33Mat(j)
  UMF_lMat%ColA(k) = 2*lMat%nu + lMat%ColA(j)
 end do
 UMF_lMat%LdA(2*lMat%nu + i+1) = UMF_lMat%LdA(2*lMat%nu + i) + kk 
end do

CALL myUmfPack_Factorize(UMF_CMat,UMF_lMat)

pause
END SUBROUTINE Init_UmfPackS
!
! ----------------------------------------------
!
SUBROUTINE Init_UmfPack()

IF (.not.ALLOCATED(UMF_CMat)) ALLOCATE (UMF_CMat(lMat%na))
UMF_CMat = DMat
IF (.not.ALLOCATED(UMF_lMat%ColA)) ALLOCATE (UMF_lMat%ColA(lMat%na))
IF (.not.ALLOCATED(UMF_lMat%LdA)) ALLOCATE (UMF_lMat%LdA(lMat%nu+1))
UMF_lMat%ColA = lMat%ColA
UMF_lMat%LdA  = lMat%LdA
UMF_lMat%nu   = lMat%nu
UMF_lMat%na   = lMat%na
CALL myUmfPack_Factorize(UMF_CMat,UMF_lMat)

END SUBROUTINE Init_UmfPack
!
! ----------------------------------------------
!
SUBROUTINE Solve_UmfPackS()
real*8, allocatable :: sol(:),rhs(:)

MeshDef%valX_old = MeshDef%valX
MeshDef%valY_old = MeshDef%valY
MeshDef%valZ_old = MeshDef%valZ

allocate(sol(3*MeshDef%ndof))
allocate(rhs(3*MeshDef%ndof))

sol(0*MeshDef%ndof+1:) = MeshDef%valX
sol(1*MeshDef%ndof+1:) = MeshDef%valY
sol(2*MeshDef%ndof+1:) = MeshDef%valZ

rhs(0*MeshDef%ndof+1:) = MeshDef%defX
rhs(1*MeshDef%ndof+1:) = MeshDef%defY
rhs(2*MeshDef%ndof+1:) = MeshDef%defZ

CALL myUmfPack_Solve(sol,rhs,UMF_CMat,UMF_lMat,1)

DO i=1,MeshDef%ndof
 MeshDef%valX(i) = MeshDef%valX_old(i) - sol(0*MeshDef%ndof+i)
 MeshDef%valY(i) = MeshDef%valY_old(i) - sol(1*MeshDef%ndof+i)
 MeshDef%valZ(i) = MeshDef%valZ_old(i) - sol(2*MeshDef%ndof+i)
END DO

END SUBROUTINE Solve_UmfPackS
!
! ----------------------------------------------
!
SUBROUTINE Solve_UmfPack()

MeshDef%valX_old = MeshDef%valX
MeshDef%valY_old = MeshDef%valY
MeshDef%valZ_old = MeshDef%valZ

CALL myUmfPack_Solve(MeshDef%valX,MeshDef%defX,UMF_CMat,UMF_lMat,1)
CALL myUmfPack_Solve(MeshDef%valY,MeshDef%defY,UMF_CMat,UMF_lMat,1)
CALL myUmfPack_Solve(MeshDef%valZ,MeshDef%defZ,UMF_CMat,UMF_lMat,1)

DO i=1,MeshDef%ndof
 MeshDef%valX(i) = MeshDef%valX_old(i) - MeshDef%valX(i)
 MeshDef%valY(i) = MeshDef%valY_old(i) - MeshDef%valY(i)
 MeshDef%valZ(i) = MeshDef%valZ_old(i) - MeshDef%valZ(i)
END DO

END SUBROUTINE Solve_UmfPack
!
! ----------------------------------------------
!
SUBROUTINE GetRes_MeshDef()
REAL*8 resX,resY,resZ

CALL LL21(MeshDef%defX,MeshDef%ndof,resX)
CALL LL21(MeshDef%defY,MeshDef%ndof,resY)
CALL LL21(MeshDef%defZ,MeshDef%ndof,resZ)

WRITE(*,'(A,3ES12.4)') 'residuals:',resX,resY,resZ

END SUBROUTINE GetRes_MeshDef
!
! ----------------------------------------------
!
SUBROUTINE MatDef_MeshDef()
REAL*8 resX,resY,resZ

CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defX,1d0,0d0)
CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defY,1d0,0d0)
CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defZ,1d0,0d0)

END SUBROUTINE MatDef_MeshDef
!
! ----------------------------------------------
!
SUBROUTINE MatDef_MeshDefS()
REAL*8 resX,resY,resZ

CALL LAX17(S11Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defX,1d0,0d0)
CALL LAX17(S12Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defX,1d0,1d0)
CALL LAX17(S13Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defX,1d0,1d0)

CALL LAX17(S21Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defY,1d0,0d0)
CALL LAX17(S22Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defY,1d0,1d0)
CALL LAX17(S23Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defY,1d0,1d0)

CALL LAX17(S31Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defZ,1d0,0d0)
CALL LAX17(S32Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defZ,1d0,1d0)
CALL LAX17(S33Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defZ,1d0,1d0)

END SUBROUTINE MatDef_MeshDefS
!
! ----------------------------------------------
!
SUBROUTINE Boundary_MeshDef_Def()
INTEGER i

DO i=1,MeshDef%ndof
 IF (MeshDef%knprX(i).eq.1) THEN
  MeshDef%defX(i) = 0d0
 END IF
 IF (MeshDef%knprY(i).eq.1) THEN
  MeshDef%defY(i) = 0d0
 END IF
 IF (MeshDef%knprZ(i).eq.1) THEN
  MeshDef%defZ(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_MeshDef_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_MeshDef_Mat(DA,KLD)
REAL*8  DA(*)
INTEGER KLD(*),ICOL,I

DO I=1,MeshDef%ndof
 if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
   DA(KLD(I))=1d0
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA(ICOL)=0d0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_MeshDef_Mat
!
! ----------------------------------------------
!
SUBROUTINE Boundary_MeshDef_MatS(DS11,DS22,DS33,&
           DS12,DS13,DS23,DS21,DS31,DS32,KLD)
REAL*8  DS11(*),DS22(*),DS33(*)
REAL*8  DS12(*),DS13(*),DS23(*),DS21(*),DS31(*),DS32(*)
INTEGER KLD(*),ICOL,I

DO I=1,MeshDef%ndof
 if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
  ICOL = KLD(I)
  DS12(ICOL) = 0d0
  DS13(ICOL) = 0d0
  DO ICOL=KLD(I)+1,KLD(I+1)-1
   DS11(ICOL) = 0d0
   DS12(ICOL) = 0d0
   DS13(ICOL) = 0d0
  END DO
  
  ICOL = KLD(I)
  DS23(ICOL) = 0d0
  DS21(ICOL) = 0d0
  DO ICOL=KLD(I)+1,KLD(I+1)-1
   DS22(ICOL) = 0d0
   DS23(ICOL) = 0d0
   DS21(ICOL) = 0d0
  END DO

  ICOL = KLD(I)
  DS31(ICOL) = 0d0
  DS32(ICOL) = 0d0
  DO ICOL=KLD(I)+1,KLD(I+1)-1
   DS33(ICOL) = 0d0
   DS31(ICOL) = 0d0
   DS32(ICOL) = 0d0
  END DO
 END IF
END DO

END SUBROUTINE Boundary_MeshDef_MatS
!
! ----------------------------------------------
!
SUBROUTINE InitMeshDef()

ILEV = mg_Mesh%nlmax

! Building up the matrix strucrures
CALL Create_MatStruct()

if (bDefTensor) then
 ! Building up the Deformation tensor operator
 CALL Create_ConstDeformationMat()
else
 ! Building up the Diffusion operator
 CALL Create_ConstDiffMat()
end if

! Initialize the MeshDef field structures
WRITE(*,*) 'Init MeshDef structures'
CALL Init_MeshDefField()

if (.not.allocated(OrigCoor)) allocate(OrigCoor(3,MeshDef%ndof))
OrigCoor(1:3,:) = mg_mesh%level(ilev)%dcorvg(1:3,1:MeshDef%ndof)

END SUBROUTINE InitMeshDef

END MODULE MeshProcPDE
