MODULE Transport_Q1

USE def_LinScalar
USE PP3D_MPI, ONLY:E011Sum,E011Knpr,Comm_NLComplete,&
    Comm_Maximum,Comm_Summ,myid,master,CommSum,Comm_SummN,myMPI_barrier
USE Transport_Q2P1, ONLY: QuadSc,ParKNPR,mgDiffCoeff,&
    myBoundary,myQ2Coor,&
    MoveInterfacePoints,myALE,Properties,getmeshvelocity,Temperature
USE var_QuadScalar, ONLY: myMG,myHeatObjects,Properties,dIntegralHeat
USE mg_LinScalar, ONLY : mgProlRestInit
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,MyMaterials

IMPLICIT NONE

TYPE(lScalar) Tracer
TYPE(lScalar3) Tracer3

CHARACTER*25 :: CInitFile="#data/LS02"

REAL*8 dArea1,dFlux1,dArea2,dFlux2,dHeatSource

include 'LinSc_user_include.h'

! The handler function for the initial condition
procedure(LinSc_IC), pointer :: LinSc_IC_ptr => null()
! The handler function for the boundary condition
procedure(LinSc_BC), pointer :: LinSc_BC_ptr => null()
! The handler function for the source term
procedure(LinSc_SRC), pointer :: LinSc_SRC_ptr => null()

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE SetPointersToCase()
implicit none

! LinSc_IC_ptr  => LinSc_InitCond
! LinSc_BC_ptr  => Boundary_LinSc_Val
! LinSc_SRC_ptr => AddSource

END SUBROUTINE SetPointersToCase
!
! ----------------------------------------------
!
SUBROUTINE Transport_LinScalar(mfile,INL)
INTEGER mfile,INL
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

NLMAX = NLMAX + 1

thstep = 0.5d0*tstep

! advect the scalar field
IF (myid.ne.0) THEN

 CALL Build_LinSc_Convection()
 IF (Tracer%prm%AFC) CALL InitAFC_General_LinScalar()

! Assemble the right hand side
 CALL LCL1(Tracer%def,Tracer%ndof)
 CALL Matdef_LinScalar(Tracer,1,0)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def()

! Store the constant right hand side
 Tracer%rhs = Tracer%def

! Assemble the defect vector and fine level matrix
 CALL Matdef_LinScalar(Tracer,-1,1)
 CALL E011Sum(Tracer%def)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def()

! Save the old solution
 CALL LCP1(Tracer%val(NLMAX)%x,Tracer%val_old,Tracer%ndof)

! Compute the defect
 CALL Resdfk_General_LinScalar(Tracer,ResTemp,DefTemp,RhsTemp)

END IF

! WRITE(*,'(I5,3D12.4)') myid,ResTemp,DefTemp,RhsTemp

CALL COMM_Maximum(RhsTemp)
DefTempCrit=MAX(RhsTemp*Tracer%prm%defCrit,Tracer%prm%MinDef)

CALL Protocol_linScalar(mfile,Tracer,0,&
     ResTemp,DefTemp,DefTempCrit," Scalar advection ")

DO INL=1,Tracer%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_LinScalar(Tracer,ParKNPR,&
Boundary_LinSc_Val,Boundary_LinSc_Mat)

!CALL myMPI_barrier()
!write(*,*) 'sdf sdf  s  fsdf sd fs  f',myid
!pause

IF (myid.ne.0) THEN

!!!!          Checking the quality of the result           !!!!
! Restore the constant right hand side
 Tracer%def = Tracer%rhs

! Assemble the defect vector and fine level matrix
 CALL Matdef_LinScalar(Tracer,-1,0)
 CALL E011Sum(Tracer%def)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def()

! Save the old solution
 CALL LCP1(Tracer%val(NLMAX)%x,Tracer%val_old,Tracer%ndof)

! Compute the defect
 CALL Resdfk_General_LinScalar(Tracer,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
CALL COMM_Maximum(RhsTemp)
CALL Protocol_linScalar(mfile,Tracer,INL,&
     ResTemp,DefTemp,RhsTemp)

IF ((DefTemp.LE.DefTempCrit).AND.&
    (INL.GE.Tracer%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

NLMAX = NLMAX - 1

END SUBROUTINE Transport_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Init_LinScalar(mfile)
INTEGER mfile

NLMAX = NLMAX + 1

IF (myid.ne.master) THEN

 ! Building up the matrix strucrures
 CALL Create_MatStruct()

 ! Building up the matrix strucrures
 CALL Create_AFCStruct()

 ! Building up linear operators
 ! Mass matrix
 CALL Create_MassMat()

! Mass matrix
 CALL Create_LMassMat()

! Convection matrix (only allocation)
 CALL Create_LKonvMat()

! Diffusion matrix 
 CALL Create_ConstDiffMat(Properties%DiffCoeff(1))

! Iteration matrix (only allocation)
 CALL Create_AMat()

! Initialize the scalar quantity
 CALL Initialize(Tracer)

! ! Set the types of boundary conditions (set up knpr)
!  CALL Create_Knpr(LinSc_Knpr)

END IF


Tracer%cName        = "Tracer"
Tracer%prm%SolvIter = 8
Tracer%prm%AFC      = .TRUE.
Tracer%prm%NLmin    = 4
Tracer%prm%NLmax    = 8
Tracer%prm%defCrit  =1d-4
Tracer%prm%epsCrit  =1d-3
Tracer%prm%MinDef   =1d-10
Tracer%prm%SolvType =1

NLMAX = NLMAX - 1

END SUBROUTINE Init_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE InitHeatObjects()

ALLOCATE(myHeatObjects%Block(Tracer%ndof))
ALLOCATE(myHeatObjects%Wire(Tracer%ndof))
ALLOCATE(myHeatObjects%Channel(Tracer%ndof))
ALLOCATE(myHeatObjects%Segment(Tracer%ndof))
myHeatObjects%Block = 0d0
myHeatObjects%Wire  = 0d0
myHeatObjects%Channel  = 0d0
myHeatObjects%Segment  = 0 ! Air

END SUBROUTINE InitHeatObjects
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar()
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

CALL LinSc_InitCond(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL Boundary_LinSc_Val()

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE SetTracerToLoadedTemperatue(sub_BC)
USE PP3D_MPI, ONLY : myid
EXTERNAL sub_BC

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

! Set the types of boundary conditions (set up knpr)
CALL Create_Knpr(LinSc_Knpr)

if (myid.ne.0) then
 Tracer%Val(NLMAX)%x = Temperature
end if

! Set boundary conditions
CALL sub_BC()

NLMAX = NLMAX - 1

END SUBROUTINE SetTracerToLoadedTemperatue
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar_EWIKON(sub_IC,sub_BC)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
EXTERNAL sub_IC,sub_BC

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

! Set the types of boundary conditions (set up knpr)
CALL Create_Knpr(LinSc_Knpr)

! Set The initial Conditions
CALL sub_IC(mg_mesh%level(ilev)%dcorvg)

if (myid.ne.0) then
 Tracer%Val(NLMAX)%x = Temperature
end if

! Set boundary conditions
CALL sub_BC()

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_EWIKON
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar_XSE(sub_IC,sub_BC)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

EXTERNAL sub_IC,sub_BC
integer i

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

! Set the types of boundary conditions (set up knpr)
CALL Create_Knpr(LinSc_Knpr_XSE)

! Set The initial Conditions
CALL sub_IC(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL sub_BC()

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_XSE
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar_General(sub_IC,sub_BC)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

EXTERNAL sub_IC,sub_BC
integer i

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

! Set the types of boundary conditions (set up knpr)
CALL Create_Knpr(LinSc_Knpr)

! Set The initial Conditions
CALL sub_IC(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL sub_BC()

if (myid.ne.0) then
 Temperature = Tracer%Val(NLMAX)%x 
end if

! do i=1,size(Temperature)
!  if (Temperature(i).gt.201d0) then
!   write(*,*) "asdfadad:" ,Temperature(i)
!  end if
! end do

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_General
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr(dcorvg)
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.2d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i,iSeg,jSeg,k

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer%knpr(I) = 0

 IF (myBoundary%iTemperature(i).ne.0) THEN
  Tracer%knpr(I) = myBoundary%iTemperature(i)
 END IF
 
END DO

if (myid.ne.0) then
 DO jSeg=1,mySigma%NumberOfSeg
  DO i=1,Tracer%ndof
   iSeg = myHeatObjects%Segment(i)
   IF (iSeg.eq.jSeg) THEN
    IF (mySigma%mySegment(iSeg)%TemperatureBC.eq.'CONSTANT') THEN
     Tracer%knpr(I) = 3
    END IF
   END IF
  END DO
 END DO
end if

END SUBROUTINE LinSc_Knpr
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr_XSE(dcorvg)
use, intrinsic :: ieee_arithmetic

REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.2d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i,iSeg,jSeg,k
real*8 :: myInf

if(ieee_support_inf(myInf))then
  myInf = ieee_value(myInf, ieee_negative_inf)
endif

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer%knpr(I) = 0

 IF (myBoundary%iInflow(i).eq.10) Tracer%knpr(I) = 1
 
 IF (myBoundary%bWall(i).and.myProcess%Ta.ne.myInf) THEN
  Tracer%knpr(I) = 2
 END IF
 
 IF (myBoundary%iInflow(i).gt.10.and.myProcess%Ti.ne.myInf) THEN
  Tracer%knpr(I) = 3
 END IF

END DO

END SUBROUTINE LinSc_Knpr_XSE
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr_Weber(dcorvg)
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.2d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer%knpr(I) = 0

 DIST = SQRT((X-PX)**2d0+(Y-PY)**2d0) - RAD
 IF (abs(Z-74.68d0).lt.1d-4) THEN
  Tracer%knpr(I) = 1
 END IF

END DO

END SUBROUTINE LinSc_Knpr_Weber
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Def()
INTEGER i

DO i=1,Tracer%ndof
 IF (Tracer%knpr(i).ne.0) THEN
  Tracer%def(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_LinSc_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Mat(VA,KLD,KNPR)
REAL*4  VA(*)
INTEGER KLD(*),KNPR(*),ICOL,I

DO I=1,Tracer%ndof
 IF (KNPR(I).ne.0) THEN
   VA(KLD(I))=1E-16
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    VA(ICOL)=0E0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_LinSc_Mat
!
! ----------------------------------------------
!
SUBROUTINE Build_LinSc_Convection()
INTEGER I,J,jLEV
EXTERNAL E011
integer invt, inet, inat, inel

ILEV=NLMAX
JLEV = ILEV-1
CALL SETLEV(2)

Kmat = 0d0

CALL Conv_LinSc2(QuadSc%valU,QuadSc%valV,QuadSc%valW,Kmat,&
lMat%nu,lMat%ColA,lMat%LdA,&
mg_mesh%level(ilev)%kvert,&
mg_mesh%level(ilev)%karea,&
mg_mesh%level(ilev)%kedge,&
mg_mesh%level(ilev)%dcorvg,&
mg_mesh%level(ilev)%kadj,&
mg_mesh%level(jlev)%kvert,&
mg_mesh%level(jlev)%karea,&
mg_mesh%level(jlev)%kedge,&
mg_mesh%level(jlev)%nel,&
mg_mesh%level(jlev)%nvt,&
mg_mesh%level(jlev)%net,&
mg_mesh%level(jlev)%nat,E011)

!ILEV=NLMAX
!CALL SETLEV(2)

!CALL Conv_LinSc1(QuadSc%valU,QuadSc%valV,QuadSc%valW,Kmat,&
!lMat%nu,lMat%ColA,lMat%LdA,KWORK(L(LVERT)),KWORK(L(LAREA)),&
!KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),E011)

END SUBROUTINE Build_LinSc_Convection
!
! ----------------------------------------------
!
SUBROUTINE dump_LinScalar_out()
USE var_QuadScalar, only:knvt
INTEGER ifilen,itwx,i
CHARACTER COFile*15
DATA ifilen/0/
integer :: iend 

 IF (myid.EQ.0) RETURN

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='#ns/LS        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1,A1)') 0,itwx,'_'
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2,A1)') itwx,'_'
 IF (myid.lt.10) WRITE(COFile(10:11),'(I1,I1)') 0,myid
 IF (myid.ge.10) WRITE(COFile(10:11),'(I2)') myid
 OPEN (UNIT=2,FILE=COFile,FORM="FORMATTED")

 iend = KNVT(NLMAX)
 DO I=1,iend
  WRITE(2,'(G19.12)') Tracer%val(NLMAX)%x(i)
 END DO

 CLOSE(2)

END SUBROUTINE dump_LinScalar_out
!
! ----------------------------------------------
!
SUBROUTINE dump_LinScalar_in(cdump)
USE var_QuadScalar, only:knvt
implicit none
INTEGER IFl,itwx,i
CHARACTER CIFile*15,cdump*(*)
integer :: iend 

 IF (myid.EQ.0) RETURN

 IFl=LEN(TRIM(ADJUSTL(cdump)))
 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.lt.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I1,I1)') '_',0,myid
 IF (myid.ge.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I2)') '_',myid
 OPEN (UNIT=1,FILE=CIFile,STATUS="OLD",FORM="FORMATTED")

 iend = KNVT(NLMAX)
 DO I=1,iend
  READ(1,*) Tracer%val(NLMAX)%x(i)
 END DO

 CLOSE(1)

END SUBROUTINE dump_LinScalar_in
!
! ----------------------------------------------
!
subroutine InitLinearOperators(mfile, mgMesh)
use PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
use var_QuadScalar, only : tMultiMesh

implicit none

integer, intent(in) :: mfile

type(tMultiMesh), intent(inout) :: mgMesh

myHeatObjects%Block = 0d0
myHeatObjects%Wire  = 0d0
myHeatObjects%Channel  = 0d0
myHeatObjects%Segment  = 0 ! Air
call updateHeatGeometry(mfile)

! call OperatorRegenaration(1)
! call OperatorRegenaration(2)
! call OperatorRegenaration(3)

end subroutine InitLinearOperators
!
! ----------------------------------------------
!
SUBROUTINE updateHeatGeometry(mfile)
use geometry_processing, only : calcDistanceFunction_heat, dEpsDist

integer, intent(in) :: mfile

integer :: i,j,ivt,iSeg,jSeg,iMat,ndof
REAL*8 dAvgDiffCoeff
real*8, allocatable :: volume(:)

REAL :: tttt0,tttt1

CALL myMPI_Barrier()
CALL ZTIME(tttt0)

if (myid.ne.0) then

 ILEV=NLMAX
 CALL SETLEV(2)
 QuadSc%AuxU = dEpsDist
 QuadSc%AuxV = dEpsDist
 QuadSc%AuxW = dEpsDist

 !MixerKNPR(:) = 0

 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."HEAT") THEN
  CALL calcDistanceFunction_heat(mg_mesh%level(ilev)%dcorvg,&
                            mg_mesh%level(ilev)%kvert,&
                            mg_mesh%level(ilev)%kedge,&
                            mg_mesh%level(ilev)%karea,&
                            mg_mesh%level(ilev)%nel,&
                            mg_mesh%level(ilev)%nvt,&
                            mg_mesh%level(ilev)%nat,&
                            mg_mesh%level(ilev)%net,&
                            QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%defU)
 END IF

 CALL myMPI_Barrier()
 CALL ZTIME(tttt1)
 IF (myid.eq.1) WRITE(mterm,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"
 IF (myid.eq.1) WRITE(mfile,"(A,F6.1,A)") "Time used for FINE mesh distance estimation was: ", tttt1-tttt0, "s!"

 NLMAX=NLMAX+1
 ILEV=NLMAX
 ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net + mg_mesh%level(ilev)%nat

 DO I=1,mg_mesh%level(ilev)%nel
  
  dAvgDiffCoeff = 0d0
  do j=1,8
   ivt = mg_mesh%level(ilev)%kvert(j,i)
   iSeg = myHeatObjects%Segment(ivt)
   if (iSeg.eq.0) then
    iMat= 0
   else
    iMat = mySigma%mySegment(iSeg)%MatInd
   end if
   dAvgDiffCoeff = dAvgDiffCoeff + 0.125d0*myMaterials(iMat)%Alpha
  end do
  
  mgDiffCoeff(ilev)%x(i) = dAvgDiffCoeff
 END DO

! Mass matrix
 CALL Create_MassMat()

! Mass Lumped matrix
 CALL Create_LMassMat()

 NLMAX = NLMAX - 1
end if

allocate(volume(0:mySigma%NumberOfSeg))

if (myid.ne.0) then
 DO jSeg=0,mySigma%NumberOfSeg
  Volume(jSeg) = 0d0
  DO i=1,Tracer%ndof
   iSeg = myHeatObjects%Segment(i)
   IF (iSeg.eq.jSeg) THEN
    Volume(jSeg) = Volume(jSeg) + MLmat(i)
   END IF
  END DO
 END DO
end if

CALL Comm_SummN(Volume,mySigma%NumberOfSeg+1)
mySigma%mySegment(:)%Volume = Volume(1:mySigma%NumberOfSeg)
if (myid.eq.1) write(*,'(A,10ES12.4)') "Volume of segments: ",Volume
deallocate(volume)
! pause

if (myid.ne.0) then

 NLMAX=NLMAX+1
 
 ILEV=NLMAX
 CALL SETLEV(2)

 ! Diffusion matrix 
 CALL Create_LambdaDiffMat()

 ! Mass matrix
 CALL Create_CpRhoMassMat()

 ! Convection matrix
 CALL Create_RhoCpConvMat(QuadSc%valU,QuadSc%valV,QuadSc%valW)

 ! Mass Lumped matrix
 CALL Create_LRhoCpMassMat()

 NLMAX = NLMAX - 1
end if

END SUBROUTINE  updateHeatGeometry
!
! ----------------------------------------------
!
include 'LinSc_transport_extensions.f90'
include 'LinSc_user.f90'
!
! ----------------------------------------------
!
END MODULE Transport_Q1
