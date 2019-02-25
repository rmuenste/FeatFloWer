MODULE MeshProcPDE

USE def_LinScalar
USE types, ONLY: tMultiMesh,tMesh
USE var_QuadScalar, ONLY : UMF_CMat,UMF_lMat,myBoundary
USE PP3D_MPI, ONLY:myid,master
USE UMFPackSolver, ONLY : myUmfPack_Factorize,myUmfPack_Solve
USE Parametrization, ONLY: ParametrizeBndryPoints_STRCT
USE MeshProcDef, ONLY : norm_u,norm_v,norm_w,norm_d,cOutputFolder


TYPE(lScalar3) MeshDef
REAL*8, ALLOCATABLE :: OrigCoor(:,:)
contains

SUBROUTINE MeshDefPDE(bDefTensor)
Logical :: bDefTensor

ILEV = mg_Mesh%nlmax

!mg_Mesh%BndryNodes(:)%bOuterPoint = .false.

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
 CALL MatDef_MeshDefS()
else
 CALL MatDef_MeshDef()
end if

CALL Boundary_MeshDef_Def()

if (bDefTensor) then
 CALL Boundary_MeshDef_MatS(S11Mat,S22Mat,S33Mat,&
                            S12Mat,S13Mat,S23Mat,&
                            S21Mat,S31Mat,S32Mat,&
                            lMat%LdA)
ELSE
 CALL Boundary_MeshDef_Mat(DMat,lMat%LdA)
END IF

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

do i=1,MeshDef%ndof
 if (.not.mg_Mesh%BndryNodes(i)%bOuterPoint) then
  mg_mesh%level(ilev)%dcorvg(1,i) = MeshDef%valX(i)
  mg_mesh%level(ilev)%dcorvg(2,i) = MeshDef%valY(i)
  mg_mesh%level(ilev)%dcorvg(3,i) = MeshDef%valZ(i)
 else
!   write(*,'(I0,6ES12.4)') i,MeshDef%valX(i),MeshDef%defX(i),&
!   MeshDef%valY(i),MeshDef%defY(i),MeshDef%valZ(i),MeshDef%defZ(i)
 end if
END DO

! mg_mesh%level(ilev)%dcorvg(1,1:MeshDef%ndof) = MeshDef%valX
! mg_mesh%level(ilev)%dcorvg(2,1:MeshDef%ndof) = MeshDef%valY
! mg_mesh%level(ilev)%dcorvg(3,1:MeshDef%ndof) = MeshDef%valZ
                                            
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
integer iC,kk

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


! open(unit=3,file=adjustl(trim(cOutputFolder))//'/mat.dat')
! write(3,'(A,I0,A,I0)') 'nu= 3*',lMat%nu,', na= ',UMF_lMat%na
! write(3,'(A10," ",A10," ",A12)') 'i','j','a_ij'
! 
! kk = 0
! iC = 0
! do i=1,UMF_lMat%nu
!  kk = kk + 1
!  if (kk.gt.lMat%nu) then
!   kk = 0
!   iC= iC + 1
!  end if
!  
!  do j=UMF_lMat%LdA(i),UMF_lMat%LdA(i+1)-1
!   k  = UMF_lMat%ColA(j)
!   if (iC.eq.0) write(3,'(A1," "I10," ",I10," ",ES12.4)') 'x',i,k,UMF_CMat(j)
!   if (iC.eq.1) write(3,'(A1," "I10," ",I10," ",ES12.4)') 'y',i,k,UMF_CMat(j)
!   if (iC.eq.2) write(3,'(A1," "I10," ",I10," ",ES12.4)') 'z',i,k,UMF_CMat(j)
!   if (iC.eq.3) write(3,'(A1," "I10," ",I10," ",ES12.4)') '*',i,k,UMF_CMat(j)
!  end do
! end do
! close(3)
! pause

CALL myUmfPack_Factorize(UMF_CMat,UMF_lMat)

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

sol(0*MeshDef%ndof+1:) = 0d0!MeshDef%valX
sol(1*MeshDef%ndof+1:) = 0d0!MeshDef%valY
sol(2*MeshDef%ndof+1:) = 0d0!MeshDef%valZ

rhs(0*MeshDef%ndof+1:) = MeshDef%defX
rhs(1*MeshDef%ndof+1:) = MeshDef%defY
rhs(2*MeshDef%ndof+1:) = MeshDef%defZ

CALL myUmfPack_Solve(sol,rhs,UMF_CMat,UMF_lMat,1)

DO i=1,MeshDef%ndof
 MeshDef%valX(i) = MeshDef%valX_old(i) + sol(0*MeshDef%ndof+i)
 MeshDef%valY(i) = MeshDef%valY_old(i) + sol(1*MeshDef%ndof+i)
 MeshDef%valZ(i) = MeshDef%valZ_old(i) + sol(2*MeshDef%ndof+i)
END DO

deallocate(sol)
deallocate(rhs)

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
 MeshDef%valX(i) = MeshDef%valX_old(i) + MeshDef%valX(i)
 MeshDef%valY(i) = MeshDef%valY_old(i) + MeshDef%valY(i)
 MeshDef%valZ(i) = MeshDef%valZ_old(i) + MeshDef%valZ(i)
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
MeshDef%valX,MeshDef%defX,-1d0,0d0)
CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defY,-1d0,0d0)
CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defZ,-1d0,0d0)

END SUBROUTINE MatDef_MeshDef
!
! ----------------------------------------------
!
SUBROUTINE MatDef_MeshDefS()
REAL*8 resX,resY,resZ

CALL LAX17(S11Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defX,-1d0,0d0)
CALL LAX17(S12Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defX,-1d0,1d0)
CALL LAX17(S13Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defX,-1d0,1d0)

CALL LAX17(S21Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defY,-1d0,0d0)
CALL LAX17(S22Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defY,-1d0,1d0)
CALL LAX17(S23Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defY,-1d0,1d0)

CALL LAX17(S31Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valX,MeshDef%defZ,-1d0,0d0)
CALL LAX17(S32Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valY,MeshDef%defZ,-1d0,1d0)
CALL LAX17(S33Mat,lMat%ColA,lMat%LdA,lMat%nu,&
MeshDef%valZ,MeshDef%defZ,-1d0,1d0)

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
SUBROUTINE InitMeshDef(bDefTensor)
Logical :: bDefTensor

ILEV = mg_Mesh%nlmax

! Building up the matrix strucrures
CALL Create_MatStruct()

if (bDefTensor) then
 ! Building up the Deformation tensor operator
 CALL Create_ConstDeformationMat()
else
 ! Building up the Diffusion operator
 CALL Create_ConstDiffMat(1d0)
end if

! Initialize the MeshDef field structures
WRITE(*,*) 'Init MeshDef structures'
CALL Init_MeshDefField()

if (.not.allocated(OrigCoor)) allocate(OrigCoor(3,MeshDef%ndof))
OrigCoor(1:3,:) = mg_mesh%level(ilev)%dcorvg(1:3,1:MeshDef%ndof)

END SUBROUTINE InitMeshDef
!----------------------------------------------------------
SUBROUTINE InitBoundaryFields()
EXTERNAL E013
integer ndof

ILEV = mg_Mesh%nlmax
ndof = mg_Mesh%level(ilev)%nvt + mg_Mesh%level(ilev)%net + mg_Mesh%level(ilev)%nat + mg_Mesh%level(ilev)%nel
allocate(norm_u(ndof),norm_v(ndof),norm_w(ndof),norm_d(ndof))
allocate(myBoundary%LS_zero(ndof))
allocate(myBoundary%iPhase(ndof))
myBoundary%LS_zero = 0
myBoundary%iPhase = 1

do i=1,ndof
 if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
  myBoundary%LS_zero(i) = 1
 end if
END DO

END SUBROUTINE InitBoundaryFields
!----------------------------------------------------------
SUBROUTINE RecoverSurfaceNormals()
EXTERNAL E013
integer ndof,k

ILEV = mg_Mesh%nlmax
ndof = mg_Mesh%level(ilev)%nvt + mg_Mesh%level(ilev)%net + mg_Mesh%level(ilev)%nat + mg_Mesh%level(ilev)%nel

norm_u = 0d0
norm_v = 0d0
norm_w = 0d0
norm_d = 0d0

CALL SETLEV(2)
CALL RecoverNormals(norm_u,norm_v,norm_w,norm_d,&
 mg_mesh%level(ilev)%kvert,&
 mg_mesh%level(ilev)%karea,&
 mg_mesh%level(ilev)%kedge,&
 mg_mesh%level(ilev)%dcorvg,&
 E013)

k = 0
do i=1,ndof
 if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
  k = k + 1
  norm_u(i) = norm_u(i)/norm_d(i)
  norm_v(i) = norm_v(i)/norm_d(i)
  norm_w(i) = norm_w(i)/norm_d(i)
 end if
END DO

! ndof = mg_Mesh%level(ilev)%nvt
! open(unit=3,file=adjustl(trim(cOutputFolder))//'/norm.csv')
! !write(3,'(A,I0,A,I0)') 'nu= ',ndof, ', n_bc= ',k 
! write(3,'(A10,6(A13))') 'id, ','n_x, ','n_y, ','n_z, ','bc_x, ','bc_y, ','bc_z'
! do i=1,ndof
!  if (mg_Mesh%BndryNodes(i)%bOuterPoint) then
!   write(3,'(I10,6(A1,ES12.4))') i,', ',norm_u(i),', ',norm_v(i),', ',norm_w(i),', ',&
!   mg_mesh%level(ilev)%dcorvg(1,i),', ',mg_mesh%level(ilev)%dcorvg(2,i),', ',mg_mesh%level(ilev)%dcorvg(3,i)
!  end if
! end do
! close(3)
! 
! open(unit=3,file=adjustl(trim(cOutputFolder))//'/end_coor.csv')
! !write(3,'(A,I0,A,I0)') 'nV= ',ndof 
! write(3,'(7(A10," "))') 'id, ','x, ','y, ','z, '
! do i=1,ndof
!  write(3,'(I10,6(A1,ES12.4))') i,', ',mg_mesh%level(ilev)%dcorvg(1,i),', ',mg_mesh%level(ilev)%dcorvg(2,i),', ',mg_mesh%level(ilev)%dcorvg(3,i)
! end do
! close(3)
! 
! open(unit=3,file=adjustl(trim(cOutputFolder))//'/start_coor.csv')
! ! write(3,'(A,I0,A,I0)') 'nV= ',ndof
! write(3,'(7(A10," "))') 'id, ','x, ','y, ','z, '
! do i=1,ndof
!  write(3,'(I10,6(A1,ES12.4))') i,', ',OrigCoor(1,i),', ',OrigCoor(2,i),', ',OrigCoor(3,i)
! end do
! close(3)
! pause

END SUBROUTINE RecoverSurfaceNormals
!----------------------------------------------------------
SUBROUTINE Correct_myQ2Coor()
implicit none
REAL*8 PX,PY,PZ,dl1,dl2,dl3
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,ivt5,ivt6,ivt7,ivt8
INTEGER NeighE(2,12),NeighA(4,6),NeighAE(4,6)
DATA NeighE /1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA /1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
DATA NeighAE/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/

ILEV = mg_Mesh%nlmax

k=1
DO i=1,mg_Mesh%level(ilev)%nel
 DO j=1,6
  IF (k.eq.mg_Mesh%level(ilev)%karea(j,i)) THEN
   ivt1 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%kedge(NeighAE(1,j),i)
   ivt2 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%kedge(NeighAE(2,j),i)
   ivt3 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%kedge(NeighAE(3,j),i)
   ivt4 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%kedge(NeighAE(4,j),i)
   
   IF (myBoundary%LS_zero(mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+k).eq.0) THEN
   IF (myBoundary%LS_zero(ivt1).ne.0.OR.myBoundary%LS_zero(ivt2).ne.0.OR.&
       myBoundary%LS_zero(ivt3).ne.0.OR.myBoundary%LS_zero(ivt4).ne.0) THEN
   
    DL1 = (mg_Mesh%level(ilev)%dcorvg(1,ivt1)-mg_Mesh%level(ilev)%dcorvg(1,ivt3))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(2,ivt1)-mg_Mesh%level(ilev)%dcorvg(2,ivt3))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(3,ivt1)-mg_Mesh%level(ilev)%dcorvg(3,ivt3))**2d0 
    DL2 = (mg_Mesh%level(ilev)%dcorvg(1,ivt2)-mg_Mesh%level(ilev)%dcorvg(1,ivt4))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(2,ivt2)-mg_Mesh%level(ilev)%dcorvg(2,ivt4))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(3,ivt2)-mg_Mesh%level(ilev)%dcorvg(3,ivt4))**2d0 

    IF (DL1.lt.DL2) THEN
     PX = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(1,ivt1)+mg_Mesh%level(ilev)%dcorvg(1,ivt3))
     PY = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(2,ivt1)+mg_Mesh%level(ilev)%dcorvg(2,ivt3))
     PZ = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(3,ivt1)+mg_Mesh%level(ilev)%dcorvg(3,ivt3))
    ELSE
     PX = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(1,ivt2)+mg_Mesh%level(ilev)%dcorvg(1,ivt4))
     PY = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(2,ivt2)+mg_Mesh%level(ilev)%dcorvg(2,ivt4))
     PZ = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(3,ivt2)+mg_Mesh%level(ilev)%dcorvg(3,ivt4))
    END IF
   
    mg_Mesh%level(ilev)%dcorvg(:,mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+k)=[PX,PY,PZ]

   END IF
   END IF
   
   k = k + 1
  END IF
 END DO
END DO

DO i=1,mg_Mesh%level(ilev)%nel
 PX = 0d0
 PY = 0d0
 PZ = 0d0
 ivt1 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(1,i)
 ivt2 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(2,i)
 ivt3 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(3,i)
 ivt4 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(4,i)
 ivt5 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(5,i)
 ivt6 = mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%karea(6,i)

 IF (myBoundary%LS_zero(ivt1).ne.0.OR.myBoundary%LS_zero(ivt2).ne.0.OR.&
     myBoundary%LS_zero(ivt3).ne.0.OR.myBoundary%LS_zero(ivt4).ne.0.OR.&
     myBoundary%LS_zero(ivt5).ne.0.OR.myBoundary%LS_zero(ivt6).ne.0) THEN
   
    DL1 = (mg_Mesh%level(ilev)%dcorvg(1,ivt1)-mg_Mesh%level(ilev)%dcorvg(1,ivt6))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(2,ivt1)-mg_Mesh%level(ilev)%dcorvg(2,ivt6))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(3,ivt1)-mg_Mesh%level(ilev)%dcorvg(3,ivt6))**2d0 
    DL2 = (mg_Mesh%level(ilev)%dcorvg(1,ivt2)-mg_Mesh%level(ilev)%dcorvg(1,ivt4))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(2,ivt2)-mg_Mesh%level(ilev)%dcorvg(2,ivt4))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(3,ivt2)-mg_Mesh%level(ilev)%dcorvg(3,ivt4))**2d0 
    DL3 = (mg_Mesh%level(ilev)%dcorvg(1,ivt3)-mg_Mesh%level(ilev)%dcorvg(1,ivt5))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(2,ivt3)-mg_Mesh%level(ilev)%dcorvg(2,ivt5))**2d0 +\
          (mg_Mesh%level(ilev)%dcorvg(3,ivt3)-mg_Mesh%level(ilev)%dcorvg(3,ivt5))**2d0 

    IF (DL1.lt.DL2.and.DL1.lt.DL3) THEN
     PX = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(1,ivt1)+mg_Mesh%level(ilev)%dcorvg(1,ivt6))
     PY = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(2,ivt1)+mg_Mesh%level(ilev)%dcorvg(2,ivt6))
     PZ = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(3,ivt1)+mg_Mesh%level(ilev)%dcorvg(3,ivt6))
    END IF
    IF (DL2.lt.DL1.and.DL2.lt.DL3) THEN
     PX = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(1,ivt2)+mg_Mesh%level(ilev)%dcorvg(1,ivt4))
     PY = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(2,ivt2)+mg_Mesh%level(ilev)%dcorvg(2,ivt4))
     PZ = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(3,ivt2)+mg_Mesh%level(ilev)%dcorvg(3,ivt4))
    END IF
    IF (DL3.lt.DL1.and.DL3.lt.DL2) THEN
     PX = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(1,ivt3)+mg_Mesh%level(ilev)%dcorvg(1,ivt5))
     PY = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(2,ivt3)+mg_Mesh%level(ilev)%dcorvg(2,ivt5))
     PZ = 0.5d0*(mg_Mesh%level(ilev)%dcorvg(3,ivt3)+mg_Mesh%level(ilev)%dcorvg(3,ivt5))
    END IF
   
    mg_Mesh%level(ilev)%dcorvg(:,mg_Mesh%level(ilev)%nvt+mg_Mesh%level(ilev)%net+mg_Mesh%level(ilev)%nat+i)=[PX,PY,PZ]
    
 END IF
END DO

END SUBROUTINE Correct_myQ2Coor
!----------------------------------------------------------
END MODULE MeshProcPDE
