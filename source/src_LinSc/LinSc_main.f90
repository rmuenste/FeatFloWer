MODULE Transport_Q1

USE def_LinScalar
USE PP3D_MPI, ONLY:E011Sum,E011Knpr,Comm_NLComplete,&
    Comm_Maximum,Comm_Summ,myid,master,CommSum
USE Transport_Q2P1, ONLY: QuadSc,ParKNPR,mgDiffCoeff,&
    myBoundary,myQ2Coor,&
    MoveInterfacePoints,myALE,Properties,getmeshvelocity

IMPLICIT NONE

TYPE(lScalar) Tracer
TYPE(lScalar3) Tracer3
CHARACTER*25 :: CInitFile="#data/LS02"

CONTAINS
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
 CALL Matdef_General_LinScalar(Tracer,1,0)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def()

! Store the constant right hand side
 Tracer%rhs = Tracer%def

! Set dirichlet boundary conditions on the solution
! CALL Boundary_LinSc_Val(DWORK(L(LCORVG)))

! Assemble the defect vector and fine level matrix
 CALL Matdef_General_LinScalar(Tracer,-1,1)
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

IF (myid.ne.0) THEN

!!!!          Checking the quality of the result           !!!!
! Restore the constant right hand side
 Tracer%def = Tracer%rhs

! Assemble the defect vector and fine level matrix
 CALL Matdef_General_LinScalar(Tracer,-1,0)
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
SUBROUTINE Init_LinScalar

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
 CALL Create_DiffMat(mgDiffCoeff(NLMAX)%x)

! Iteration matrix (only allocation)
 CALL Create_AMat()

! Initialize the scalar quantity
 CALL Initialize(Tracer)

! Set the types of boundary conditions (set up knpr)
 CALL Create_Knpr(LinSc_Knpr)

END IF


Tracer%cName = "Tracer"
Tracer%prm%SolvIter = 1
Tracer%prm%AFC = .TRUE.
IF (Tracer%prm%AFC) THEN
 Tracer%prm%NLmin = 2
ELSE
 Tracer%prm%NLmin = 2
END IF
Tracer%prm%NLmax   =20
Tracer%prm%defCrit =1d-6
Tracer%prm%epsCrit =1d-3
Tracer%prm%MinDef  =1d-16
Tracer%prm%SolvType=1
Tracer%prm%iMass=1

NLMAX = NLMAX - 1

END SUBROUTINE Init_LinScalar
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
CALL Boundary_LinSc_Val(mg_mesh%level(ilev)%dcorvg)

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr(dcorvg)
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.5d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer%knpr(I) = 0

 DIST = SQRT((X-PX)**2d0+(Y-PY)**2d0) - RAD
 IF (DIST.LT.1d-4) THEN
  Tracer%knpr(I) = 1
 END IF
 IF (X.LT.1d-4) THEN
  Tracer%knpr(I) = 1
 END IF

!  IF (Y.LT.-0.0249d0) THEN
!   xx = X - 0.015d0
!   IF (ABS(xx).LT.0.0149d0) THEN
!    Tracer%knpr(I) = 1
!   END IF
!   xx = X - 0.485d0
!   IF (ABS(xx).LT.0.0149d0) THEN
!    Tracer%knpr(I) = 1
!   END IF
!  END IF

END DO

END SUBROUTINE LinSc_Knpr
!
! ----------------------------------------------
!
SUBROUTINE LinSc_InitCond(dcorvg)
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0, RY = 0.0d0, RZ = 2.4d0
REAL*8 :: RADx = 0.20d0,RADs=0.040
REAL*8 DIST
INTEGER i

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 Tracer%val(NLMAX)%x(i) = 0.0d0

END DO

END SUBROUTINE LinSc_InitCond
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Def()
INTEGER i

DO i=1,Tracer%ndof
 IF (Tracer%knpr(i).eq.1) THEN
  Tracer%def(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_LinSc_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val(dcorvg)
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0,RY = 0.0d0
REAL*8 :: RADx = 0.2d0
REAL*8 DIST
INTEGER i

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 IF (Tracer%knpr(i).eq.1) THEN
  IF (X.LT.1d-4) THEN
   Tracer%val(NLMAX)%x(i) = 0d0
  ELSE
   Tracer%val(NLMAX)%x(i) = 1d0
  END IF
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Mat(VA,KLD,KNPR)
REAL*4  VA(*)
INTEGER KLD(*),KNPR(*),ICOL,I

DO I=1,Tracer%ndof
 IF (KNPR(I).eq.1) THEN
   VA(KLD(I))=1E0
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

write(*,*)'evil subroutine'
stop
Kmat = 0d0
!CALL Conv_LinSc2(QuadSc%valU,QuadSc%valV,QuadSc%valW,Kmat,&
!lMat%nu,lMat%ColA,lMat%LdA,KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),&
!DWORK(L(LCORVG)),KWORK(L(LADJ)),KWORK(L(KLVERT(JLEV))),&
!KWORK(L(KLAREA(JLEV))),KWORK(L(KLEDGE(JLEV))),KNEL(JLEV),&
!KNVT(JLEV),KNET(JLEV),KNAT(JLEV),E011)

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
include 'LinSc_transport_extensions.f90'
!
! ----------------------------------------------
!
END MODULE Transport_Q1
