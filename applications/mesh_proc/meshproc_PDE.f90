MODULE MeshProcPDE

USE def_LinScalar
USE types, ONLY: tMultiMesh,tMesh
USE var_QuadScalar, ONLY : UMF_CMat,UMF_lMat
USE PP3D_MPI, ONLY:myid,master
USE UMFPackSolver, ONLY : myUmfPack_Factorize,myUmfPack_Solve
USE Parametrization, ONLY: ParametrizeBndryPoints_STRCT

TYPE(lScalar3) MeshDef

contains

SUBROUTINE BuildUpMatrixStruct()

ILEV = mg_Mesh%nlmax

! Building up the matrix strucrures
CALL Create_MatStruct()

! Building up the Diffusion operator
CALL Create_ConstDiffMat()

! Initialize the MeshDef field structures
WRITE(*,*) 'Init MeshDef structures'
CALL Init_MeshDefField()

CALL ParametrizeBndryPoints_STRCT(mg_mesh,ilev)

MeshDef%valX = mg_mesh%level(ilev)%dcorvg(1,1:MeshDef%ndof)
MeshDef%valY = mg_mesh%level(ilev)%dcorvg(2,1:MeshDef%ndof)
MeshDef%valZ = mg_mesh%level(ilev)%dcorvg(3,1:MeshDef%ndof)

! Impose Dirichlet BC for all boundary points
WRITE(*,*) 'Set Dirichlet BC'
CALL LinSc_Knpr_MeshDef()

CALL Boundary_MeshDef_Mat(DMat,lMat%LdA)

MeshDef%defX = 0d0
MeshDef%defY = 0d0
MeshDef%defZ = 0d0

CALL MatDef_MeshDef()

CALL Boundary_MeshDef_Def()

CALL GetRes_MeshDef()

! Initialize Umfack / factorization of the matrix
WRITE(*,*) 'Init_UmfPack'
CALL Init_UmfPack()

CALL Solve_UmfPack()

mg_mesh%level(ilev)%dcorvg(1,1:MeshDef%ndof) = MeshDef%valX
mg_mesh%level(ilev)%dcorvg(2,1:MeshDef%ndof) = MeshDef%valY
mg_mesh%level(ilev)%dcorvg(3,1:MeshDef%ndof) = MeshDef%valZ
                                            
CALL MatDef_MeshDef()

CALL Boundary_MeshDef_Def()

CALL GetRes_MeshDef()

END SUBROUTINE BuildUpMatrixStruct          
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
SUBROUTINE Solve_UmfPack()

! DO i=1,MeshDef%ndof
!  write(*,*) MeshDef%defX(i),MeshDef%defY(i),MeshDef%defZ(i)
! END DO
! pause
MeshDef%valX_old = MeshDef%valX
MeshDef%valY_old = MeshDef%valY
MeshDef%valZ_old = MeshDef%valZ

CALL myUmfPack_Solve(MeshDef%valX,MeshDef%defX,UMF_CMat,UMF_lMat,1)
CALL myUmfPack_Solve(MeshDef%valY,MeshDef%defY,UMF_CMat,UMF_lMat,1)
CALL myUmfPack_Solve(MeshDef%valZ,MeshDef%defZ,UMF_CMat,UMF_lMat,1)

! CALL LLC1(MeshDef%valX_old,MeshDef%valX,MeshDef%ndof,1D0,1D0)
! CALL LLC1(MeshDef%valY_old,MeshDef%valY,MeshDef%ndof,1D0,1D0)
! CALL LLC1(MeshDef%valZ_old,MeshDef%valZ,MeshDef%ndof,1D0,1D0)

DO i=1,MeshDef%ndof
 MeshDef%valX(i) = MeshDef%valX_old(i) - MeshDef%valX(i)
 MeshDef%valY(i) = MeshDef%valY_old(i) - MeshDef%valY(i)
 MeshDef%valZ(i) = MeshDef%valZ_old(i) - MeshDef%valZ(i)
!  write(*,*) MeshDef%valX(i),MeshDef%valY(i),MeshDef%valZ(i)
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
END MODULE MeshProcPDE
