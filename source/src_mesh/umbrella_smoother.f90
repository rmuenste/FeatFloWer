module umbrella_smoother
!-------------------------------------------------------------------------------------------------
! A module that contains several variants of the Laplacian 
! smoother 'Umbrella'. Besides the standard version of this 
! smoother type, versions with user-defined weighting functions
! are available
!-------------------------------------------------------------------------------------------------
  use var_QuadScalar, only : tMultiMesh, tQuadScalar
  use geometry_processing, only: calcDistanceFunction

  ! No implicit variables in this module
  implicit none

  contains

!
!------------------------------------------------------------------------------------------------
! Default version of the umbrella type smoother 
!------------------------------------------------------------------------------------------------
! @param myTime Simulation time
! @param mgMesh The input mesh that calculations will be performed on
! @param qscStruct Simulation values that are required for weight calculation
  subroutine us_UmbrellaSmoother(myTime, umb_steps, mgMesh, qscStruct)
  use var_QuadScalar, only: knel,knvt,knet,knat
  use PP3D_MPI, only: myid,coarse
  use def_FEAT
  implicit none

  ! Simulation time
  real*8 :: myTime

  integer, dimension(:), intent(in) :: umb_steps

  type(tMultiMesh), intent(inout) :: mgMesh

  type(tQuadScalar), intent(inout) :: qscStruct

  ! local variables
  integer :: nUsedSteps,iSwitch

  character(len=60) :: cFile

  real*8 , ALLOCATABLE :: a1(:),a2(:),a3(:),a4(:),a5(:),a6(:)

  integer :: jlev, ndof
  
  IF (myid.ne.0)then
  
    ndof = KNEL(mgMesh%nlmax) + KNVT(mgMesh%nlmax) + KNET(mgMesh%nlmax) + KNAT(mgMesh%nlmax)

    ALLOCATE(a1(ndof))
    ALLOCATE(a2(ndof))
    ALLOCATE(a3(ndof))
    ALLOCATE(a4(ndof))
    ALLOCATE(a5(ndof))
    ALLOCATE(a6(ndof))
    a1 = 0d0
    a2 = 0d0
    a3 = 0d0
    a4 = 0d0
    a5 = 0d0
    a6 = 0d0
    
    DO ILEV=NLMIN,nlmax
     
     JLEV = ILEV-1
     CALL SETLEV(2)

     nUsedSteps = umb_steps(ILEV)

     IF (nUsedSteps.gt.0) THEN

      CALL us_EdgeRunner(a1,a2,a3,a4,a5,a6,mgMesh, qscStruct, ilev,&
                         mgMesh%level(ilev)%dcorvg,&
                         mgMesh%level(ilev)%kvert,&
                         mgMesh%level(ilev)%kedge,&
                         mgMesh%level(ilev)%karea,&
                         mgMesh%level(ilev)%nel,&
                         mgMesh%level(ilev)%nvt,&
                         mgMesh%level(ilev)%net,&
                         mgMesh%level(ilev)%nat,&
                         mgMesh%level(jlev)%kvert,&
                         mgMesh%level(jlev)%kedge,&
                         mgMesh%level(jlev)%karea,&
                         mgMesh%level(jlev)%nel,&
                         mgMesh%level(jlev)%nvt,&
                         mgMesh%level(jlev)%net,&
                         mgMesh%level(jlev)%nat,&
                         nUsedSteps,myTime,1)


     end IF

     CALL us_ProlongateCoordinates(&
                                   mgMesh%level(ilev)%dcorvg,&
                                   mgMesh%level(ilev)%karea,&
                                   mgMesh%level(ilev)%kvert,&
                                   mgMesh%level(ilev)%kedge,&
                                   mgMesh%level(ilev)%nel,&
                                   mgMesh%level(ilev)%nvt,&
                                   mgMesh%level(ilev)%net,&
                                   mgMesh%level(ilev)%nat)
    
    end DO
  
  end if
  
  ILEV = mgMesh%nlmin
  CALL SETLEV(2)
  
  CALL us_ExchangeNodeValuesOnCoarseLevel(&
                                         mgMesh%level(ilev)%dcorvg,&
                                         mgMesh%level(ilev)%kvert,&
                                         mgMesh%level(ilev)%nvt,&
                                         mgMesh%level(ilev)%nel)
  
  end subroutine us_UmbrellaSmoother
  !
  ! --------------------------------------------------------------
  !
  subroutine us_EdgeRunner(f,x,y,z,w,v,mgMesh, qscStruct, ilevel, dcorvg,&
                           kvert,kedge,karea,nel,nvt,net,nat,&
                           kvert2,kedge2,karea2,nel2,nvt2,net2,nat2,&
                           nProjStep,myTime,iSwitch)

  use PP3D_MPI
  !use PP3D_MPI, only: myid,coarse
  use Parametrization, only: ParametrizeBndr
  use Sigma_User, only: mySigma
  use var_QuadScalar, only: tMultiMesh,FictKNPR
  use geometry_processing, ONLY : calcDistanceFunction, dEpsDist
  
  !USE STL_Processing, ONLY : dEpsDist
  implicit none
  
  ! Parameters
  real*8, dimension(:) :: f,x,y,z,w,v

  type(tMultiMesh), intent(inout) :: mgMesh

  type(tQuadScalar), intent(inout) :: qscStruct

  integer, intent(in) :: ilevel

  real*8, dimension(:,:) :: dcorvg
  
  integer, dimension(:,:), intent(inout) :: kedge, kvert, karea

  integer :: nel,nvt,nat,net

  integer, dimension(:,:), intent(inout) :: kedge2, kvert2, karea2

  integer :: nel2,nvt2,nat2,net2
  
  integer :: nProjStep

  real*8 :: myTime

  integer :: iSwitch
  
  ! local variables
  integer NeighE(2,12)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  integer i,j,k,ivt1,ivt2,iProjStep,iaux,iel,iSTL
  real*8 WeightE,P1(3),P2(3),daux2,daux1,PX,PY,PZ,dScale1,dScale2
  real*8 :: dOmega = 0.2d0
  real*8 DIST,dIII
  real*8 :: dCrit1,dCrit2 
  real*8 dFactor,dKernel,dPower
  REAL*4, ALLOCATABLE :: myVol(:)
  real*8, ALLOCATABLE :: DXXX(:)
  
  DO k=nvt+1,nvt+net
   v(k) = 1d0
  end DO
  
  CALL E013SUM(v)
  
  ALLOCATE(myVol(nel))
  ALLOCATE(DXXX(nvt))
  
  if (myid.eq.1) write(*,'(A$)') 'MeshSmoothening:{'
  DO iProjStep=1,nProjStep
  
   if (myid.eq.1) write(*,'(I0,A1$)') iProjStep," "

   myVol = 0e0
   CALL  SETARE(myVol,nel,kvert,dcorvg)
  
  f(1:nvt) = 0d0
  w(1:nvt) = 0d0
  
   DO iel=1,nel
    DO i=1,8
     j = kvert(i,iel)
     f(j) = f(j) + abs(myVol(iel))
     w(j) = w(j) + 1d0
    end DO
   end DO
   CALL E013SUM(f)
   CALL E013SUM(w)
  
   DO i=1,nvt
    f(i) = f(i)/w(i)
   end DO
  
   qscStruct%AuxU = dEpsDist
   qscStruct%AuxV = dEpsDist

   CALL calcDistanceFunction(dcorvg,kvert2,kedge2,karea2,&
                             nel2,nvt2,nat2,net2,&
                             qscStruct%AuxU,qscStruct%AuxV,qscStruct%AuxW)
  
   DO i=1,nvt
     PX = dcorvg(1,i)
     PY = dcorvg(2,i)
     PZ = dcorvg(3,i)
  
     CALL GetWeight(PX,PY,PZ,myTime,dFactor)
  
     f(i) = dFactor*f(i)
   end DO
  
   x(1:nvt) = 0d0
   y(1:nvt) = 0d0
   z(1:nvt) = 0d0
   w(1:nvt) = 0d0
  
   k=1
   DO i=1,nel
    DO j=1,12
     IF (k.eq.kedge(j,i)) THEN
      ivt1 = kvert(NeighE(1,j),i)
      ivt2 = kvert(NeighE(2,j),i)
      P1(:) = dcorvg(:,ivt1)
      P2(:) = dcorvg(:,ivt2)
  
      daux1 = ABS(f(ivt1))
      daux2 = ABS(f(ivt2))
      WeightE = 1d0/(v(nvt + k))
  
      x(ivt1) = x(ivt1) + WeightE*P2(1)*daux2
      y(ivt1) = y(ivt1) + WeightE*P2(2)*daux2
      z(ivt1) = z(ivt1) + WeightE*P2(3)*daux2
      w(ivt1) = w(ivt1) + WeightE*daux2
  
      x(ivt2) = x(ivt2) + WeightE*P1(1)*daux1
      y(ivt2) = y(ivt2) + WeightE*P1(2)*daux1
      z(ivt2) = z(ivt2) + WeightE*P1(3)*daux1
      w(ivt2) = w(ivt2) + WeightE*daux1
  
      k = k + 1
     end IF
    end DO
   end DO
  
   CALL E013SUM(x)
   CALL E013SUM(y)
   CALL E013SUM(z)
   CALL E013SUM(w)
  
   DO i=1,nvt
     PX = x(i)/w(i)
     PY = y(i)/w(i)
     PZ = z(i)/w(i)
     dcorvg(1,i) = MAX(0d0,(1d0-dOmega))*dcorvg(1,i) + dOmega*PX
     dcorvg(2,i) = MAX(0d0,(1d0-dOmega))*dcorvg(2,i) + dOmega*PY
     dcorvg(3,i) = MAX(0d0,(1d0-dOmega))*dcorvg(3,i) + dOmega*PZ
   end DO
  
   CALL ParametrizeBndr(mgMesh,ilevel)
  
  end DO
  
  if (myid.eq.1) write(*,'(A1)') "}"
  
  DEALLOCATE(myVol,DXXX)
  
   CONTAINS
  
  subroutine GetWeight(x,y,z,t,f)
  IMPLICIT NONE
  real*8 x,y,z,t,f,YY,dBase
  real*8 d1,d2,d3,d33,d4,d5
  real*8 f1,f2,f3,f4,f5
  integer iaux
  
  
  IF (iSwitch.EQ.0) THEN
   f = 1d0
  ELSE
  
   d1 = 0.5d0*mySigma%Dz_Out - SQRT(X*X + Y*Y)
   d1 = 5d0*6d0*(6d0/mySigma%Dz_Out)*d1
   d2 = 5d0*6d0*(6d0/mySigma%Dz_Out)*qscStruct%AuxU(i)
   d3 = 5d0*6d0*(6d0/mySigma%Dz_Out)*qscStruct%AuxV(i)
!   d2 = 1.0d0
!   d3 = 1.0d0

! if(myid.eq.1)then
!   write(*,*)'mySigma%Dz',mySigma%Dz_Out
!   write(*,*)'i',i
!   write(*,*)'AuxU',qscStruct%AuxU(i)
!   write(*,*)'d1, d2, d3',d1, d2, d3
! end if
  
   CALL KernelFunction(d1,f1)
   CALL KernelFunction(d2,f2)
   CALL KernelFunction(d3,f3)
! if(myid.eq.1)then
!   write(*,*)'f1, f2, f3',f1, f2, f3
! end if

  
   f = MIN(f1*f2*f3,25d0)
   f = f**2.3d0
!   if(myid.eq.1)then
!     write(*,*)'f',f
!   end if
!   call myMPI_Barrier()
!   pause
  
  end IF
  
  end subroutine GetWeight
  
  subroutine KernelFunction(d,f)
  real*8 d,daux,f
  
  if (d.lt.0d0) d = 2.5d0*abs(d)
  
  IF (d.lt.3.0d0) THEN
   daux = 2d0 + 1.0d0*(3.0d0-d)
  ELSE
   daux = 2d0 - 0.2d0*(d-3.0d0)
  end IF
  
  
  f = max(daux,0.8d0)
  
  end subroutine KernelFunction
  
  end
  !
  ! --------------------------------------------------------------
  !
  subroutine us_ProlongateCoordinates(dcorvg,karea,kvert,kedge,nel,nvt,net,nat)
  IMPLICIT NONE
  real*8  dcorvg(3,*)
  integer kvert(8,*),kedge(12,*),karea(6,*),nel,nvt,net,nat
  real*8 PX,PY,PZ
  integer i,j,k,ivt1,ivt2,ivt3,ivt4
  integer NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
  
  k=1
  DO i=1,nel
   DO j=1,12
    IF (k.eq.kedge(j,i)) THEN
     ivt1 = kvert(NeighE(1,j),i)
     ivt2 = kvert(NeighE(2,j),i)
     PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
     PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
     PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
     dcorvg(:,nvt+k)=[PX,PY,PZ]
     k = k + 1
    end IF
   end DO
  end DO
  
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.karea(j,i)) THEN
     ivt1 = kvert(NeighA(1,j),i)
     ivt2 = kvert(NeighA(2,j),i)
     ivt3 = kvert(NeighA(3,j),i)
     ivt4 = kvert(NeighA(4,j),i)
     PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
     PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
     PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
     dcorvg(:,nvt+net+k)=[PX,PY,PZ]
     k = k + 1
    end IF
   end DO
  end DO
  
  DO i=1,nel
   PX = 0d0
   PY = 0d0
   PZ = 0d0
   DO j=1,8
    PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
    PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
    PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
   end DO
   dcorvg(:,nvt+net+nat+i)=[PX,PY,PZ]
  end DO
  
  end subroutine us_ProlongateCoordinates
  !
  ! ----------------------------------------------
  !
  subroutine us_ExchangeNodeValuesOnCoarseLevel(dcorvg,kvert,nvt,nel)
  USE PP3D_MPI
  
  IMPLICIT NONE
  real*8 dcorvg(3,*)
  integer kvert(8,*),nvt,nel
  integer pN,i,j,pID
  integer iAUX(nvt)
  real*8 dAUX(3,nvt)
  
   IF (myid.NE.0) THEN
    pN = NVT
    CALL SendI_myMPI(pN,0)
    CALL SendD_myMPI(dcorvg,3*pN,0)
   ELSE
  
    DO i=1,nvt
     dcorvg(:,i) = 0d0
    end DO
    iAux = 0
  
    DO pID=1,subnodes
     CALL RECVI_myMPI(pN ,pID)
     CALL RECVD_myMPI(dAux,3*pN,pID)
  
     DO I=1,pN
      J = coarse%pVERTLINK(pID,I)
      iAux(J)     = iAux(J)     + 1
      dcorvg(:,j) = dcorvg(:,j) + dAux(:,i)
     end DO
  
    end DO
  
    DO I=1,nvt
     dcorvg(1,i) = dcorvg(1,i)/DBLE(iAux(i))
     dcorvg(2,i) = dcorvg(2,i)/DBLE(iAux(i))
     dcorvg(3,i) = dcorvg(3,i)/DBLE(iAux(i))
    end DO
   
  end IF
  
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  
  
  end subroutine us_ExchangeNodeValuesOnCoarseLevel

end module umbrella_smoother
