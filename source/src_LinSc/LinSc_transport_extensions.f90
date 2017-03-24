!
! ----------------------------------------------
!
SUBROUTINE Transport_Q1_displacement(mfile,INL)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
INTEGER mfile,INL
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

NLMAX = NLMAX + 1

thstep = 0d0*tstep

IF (myid.ne.0) THEN

! advect the scalar field
CALL Create_NewDiffMat_Q1(myALE%Q2Coor_old,Properties%DiffCoeff(1))

! Assemble the right hand side
 CALL LCL1(Tracer3%defX,Tracer3%ndof)
 CALL LCL1(Tracer3%defY,Tracer3%ndof)
 CALL LCL1(Tracer3%defZ,Tracer3%ndof)
 CALL Matdef_General_LinScalar_Q1(Tracer3,1,0)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def_Q1()

! Store the constant right hand side
 Tracer3%rhs = 0d0 !Tracer%def

 thstep = tstep

! Set dirichlet boundary conditions on the solution
 CALL Boundary_LinSc_Val_Q1(mg_mesh%level(NLMAX)%dcorvg)

! Assemble the defect vector and fine level matrix
 CALL Matdef_General_LinScalar_Q1(Tracer3,-1,1)
 CALL E011Sum(Tracer3%defX)
 CALL E011Sum(Tracer3%defY)
 CALL E011Sum(Tracer3%defZ)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def_Q1()

! Save the old solution
 CALL LCP1(Tracer3%valX(NLMAX)%x,Tracer3%valX_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valY(NLMAX)%x,Tracer3%valY_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valZ(NLMAX)%x,Tracer3%valZ_old,Tracer3%ndof)

! Compute the defect
 CALL Resdfk_General_LinScalar_Q1(Tracer3,ResTemp,DefTemp,RhsTemp)

END IF

CALL COMM_Maximum(RhsTemp)
DefTempCrit=MAX(RhsTemp*Tracer3%prm%defCrit,Tracer3%prm%MinDef)

CALL Protocol_linScalarQ1(mfile,Tracer3,0,&
     ResTemp,DefTemp,DefTempCrit," Scalar advection ")

DO INL=1,Tracer3%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_LinScalar_Q1(Tracer3,ParKNPR,&
Boundary_LinSc_Val,Boundary_LinSc_Mat)

IF (myid.ne.0) THEN

! Restore the constant right hand side
 Tracer3%defX = 0d0 !Tracer%rhs
 Tracer3%defY = 0d0 !Tracer%rhs
 Tracer3%defZ = 0d0 !Tracer%rhs

! Assemble the defect vector and fine level matrix
 CALL Matdef_General_LinScalar_Q1(Tracer3,-1,0)
 CALL E011Sum(Tracer3%defX)
 CALL E011Sum(Tracer3%defY)
 CALL E011Sum(Tracer3%defZ)

! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def_Q1()

! Save the old solution
 CALL LCP1(Tracer3%valX(NLMAX)%x,Tracer3%valX_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valY(NLMAX)%x,Tracer3%valY_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valZ(NLMAX)%x,Tracer3%valZ_old,Tracer3%ndof)

! Compute the defect
 CALL Resdfk_General_LinScalar_Q1(Tracer3,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
CALL COMM_Maximum(RhsTemp)

!CALL Protocol_linScalar(mfile,Tracer3,INL,&
!     ResTemp,DefTemp,RhsTemp)
!CALL Protocol_linScalarQ1(mfile,Tracer3,0,&
!     ResTemp,DefTemp,DefTempCrit," Scalar advection ")

IF ((DefTemp.LE.DefTempCrit).AND.&
    (INL.GE.Tracer%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

myQ2coor(1,:) = Tracer3%valX(NLMAX)%x
myQ2coor(2,:) = Tracer3%valY(NLMAX)%x
myQ2coor(3,:) = Tracer3%valZ(NLMAX)%x

NLMAX = NLMAX - 1

ILEV=NLMAX
CALL SETLEV(2)

!CALL CoommunicateCoareGrid()
!CALL GetMeshVelocity()

END SUBROUTINE Transport_Q1_displacement
!
! ----------------------------------------------
!
SUBROUTINE Protocol_linScalarQ1(mfile,myScalar,nINL,&
           ResScalar,DefScalar,RhsScalar,cTitle)
TYPE(lscalar3), INTENT(INOUT) :: myScalar
INTEGER nINL,mfile
INTEGER i,length
REAL*8 ResScalar,DefScalar,RhsScalar
CHARACTER C1*14,C2*14,C3*14
CHARACTER, OPTIONAL:: cTitle*(*)
integer :: mterm = 6



IF (myid.eq.showID) THEN
length =  LEN(myScalar%cName)

C1='              '
C2='              '
C3='              '
WRITE(C1(1:3+length),'(A3,A7)') 'Res',myScalar%cName
WRITE(C2(1:3+length),'(A3,A7)') 'Def',myScalar%cName
WRITE(C3(1:7+length),'(A7,A7)') 'GlobDef',myScalar%cName

IF (PRESENT(cTitle)) THEN
 length = LEN(cTitle)
 IF (MOD(length,2).eq.1) length = length + 1
 length = (80-length)/2
END IF

IF (nINL.EQ.0) THEN

 IF (PRESENT(cTitle)) THEN
  WRITE(*,*) cTitle
 ELSE
  WRITE(*,*)
 END IF

 WRITE(*,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(MFILE,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(*,5)
 WRITE(MFILE,5)
 WRITE(*,'(A8,6XA10,5(6X,D10.4))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(MFILE,'(A8,6XA10,5(6X,D10.4))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(*,5)
 WRITE(MFILE,5)
 WRITE(*,'(I8,5(6X,D10.4))') 0,ResScalar,DefScalar
 WRITE(MFILE,'(I8,5(6X,D10.4))') 0,ResScalar,DefScalar
ELSE
 WRITE(*,'(I8,5(6X,D10.4))') nINL,ResScalar,DefScalar,RhsScalar
 WRITE(MFILE,'(I8,5(6X,D10.4))') nINL,ResScalar,DefScalar,RhsScalar
END IF

END IF

5  FORMAT(80('-'))
4  FORMAT(80('-'))

END SUBROUTINE Protocol_linScalarQ1
!
! ----------------------------------------------
!
!
! ----------------------------------------------
!
SUBROUTINE Init_Disp_Q1()

NLMAX = NLMAX + 1

IF (myid.ne.master) THEN

 ! Building up the matrix strucrures
 CALL Create_MatStruct_Q1()

 ! Building up the matrix strucrures
 CALL Create_AFCStruct_Q1()

 ! Building up linear operators
 ! Mass matrix
 CALL Create_MassMat_Q1()

! Mass matrix
 CALL Create_LMassMat_Q1()

! Convection matrix (only allocation)
 CALL Create_LKonvMat_Q1()

! Diffusion matrix 
 CALL Create_DiffMat_Q1(mgDiffCoeff(NLMAX)%x)

! Iteration matrix (only allocation)
 CALL Create_AMat_Q1()

END IF

 CALL Initialize_Q1(Tracer3)

! Set the types of boundary conditions (set up knpr)
 CALL Create_Knpr_Q1(LinSc_Knpr)

Tracer3%cName = "Tracer"
Tracer3%prm%SolvIter = 50
Tracer3%prm%AFC = .TRUE.
IF (Tracer3%prm%AFC) THEN
 Tracer3%prm%NLmin = 1
ELSE
 Tracer3%prm%NLmin = 1
END IF
Tracer3%prm%NLmax   =10
Tracer3%prm%defCrit =1d-6
Tracer3%prm%epsCrit =1d-3
Tracer3%prm%MinDef  =1d-16
Tracer3%prm%SolvType=1
Tracer3%prm%iMass=1

NLMAX = NLMAX - 1

END SUBROUTINE Init_Disp_Q1
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar_Q1()
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

CALL LinSc_InitCond_Q1(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL Boundary_LinSc_Val_Q1(mg_mesh%level(ilev)%dcorvg)

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr_Q1(dcorvg)
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.5d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i

DO i=1,Tracer3%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer3%knprX(I) = 0
 Tracer3%knprY(I) = 0
 Tracer3%knprZ(I) = 0

!  SPP1506
 IF (X.gt.+0.499d0) Tracer3%knprX(I) = 1
 IF (X.lt.-0.499d0) Tracer3%knprX(I) = 1

 IF (Y.gt.+3.749d0) Tracer3%knprY(I) = 1
 IF (Y.lt.-2.749d0) Tracer3%knprY(I) = 1
 
 IF (Z.gt.+0.099d0) Tracer3%knprZ(I) = 1
 IF (Z.lt.+0.001d0) Tracer3%knprZ(I) = 1
 
 IF (myBoundary%LS_zero(i)) Tracer3%knprX(I) = 1
 IF (myBoundary%LS_zero(i)) Tracer3%knprY(I) = 1
 IF (myBoundary%LS_zero(i)) Tracer3%knprZ(I) = 1
  
END DO

END SUBROUTINE LinSc_Knpr_Q1
!
! ----------------------------------------------
!
SUBROUTINE LinSc_InitCond_Q1(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0, RY = 0.0d0, RZ = 2.4d0
REAL*8 :: RADx = 0.20d0,RADs=0.040
REAL*8 DIST
INTEGER i

DO i=1,Tracer3%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 Tracer3%valX(NLMAX)%x(i) = X
 Tracer3%valY(NLMAX)%x(i) = Y
 Tracer3%valZ(NLMAX)%x(i) = Z

END DO

END SUBROUTINE LinSc_InitCond_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Def_Q1()
INTEGER i

DO i=1,Tracer3%ndof
 IF (Tracer3%knprX(i).eq.1) THEN
  Tracer3%defX(i) = 0d0
 END IF
 IF (Tracer3%knprY(i).eq.1) THEN
  Tracer3%defY(i) = 0d0
 END IF
 IF (Tracer3%knprZ(i).eq.1) THEN
  Tracer3%defZ(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_LinSc_Def_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val_Q1(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0,RY = 0.0d0
REAL*8 :: RADx = 0.2d0
REAL*8 DIST
INTEGER i

DO i=1,Tracer3%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 IF (Tracer3%knprX(i).eq.1) THEN
  Tracer3%valX(NLMAX)%x(i) = myQ2Coor(1,i)
 END IF

 IF (Tracer3%knprY(i).eq.1) THEN
  Tracer3%valY(NLMAX)%x(i) = myQ2Coor(2,i)
 END IF
 
 IF (Tracer3%knprZ(i).eq.1) THEN
  Tracer3%valZ(NLMAX)%x(i) = myQ2Coor(3,i)
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Mat_Q1(VA,KLD,KNPR)
REAL*4  VA(*)
INTEGER KLD(*),KNPR(*),ICOL,I

DO I=1,Tracer3%ndof
 IF (KNPR(I).eq.1) THEN
   VA(KLD(I))=1E0
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    VA(ICOL)=0E0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_LinSc_Mat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Build_LinSc_Convection_Q1()
INTEGER I,J,jLEV
EXTERNAL E011

ILEV=NLMAX
JLEV = ILEV-1
CALL SETLEV(2)

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

END SUBROUTINE Build_LinSc_Convection_Q1
!
! ----------------------------------------------
!
SUBROUTINE dump_LinScalar_out_Q1()
USE var_QuadScalar, only:knvt
INTEGER ifilen,itwx,i
CHARACTER COFile*15
DATA ifilen/0/

 IF (myid.EQ.0) RETURN

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='#ns/LS        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1,A1)') 0,itwx,'_'
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2,A1)') itwx,'_'
 IF (myid.lt.10) WRITE(COFile(10:11),'(I1,I1)') 0,myid
 IF (myid.ge.10) WRITE(COFile(10:11),'(I2)') myid
 OPEN (UNIT=2,FILE=COFile,FORM="FORMATTED")

 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.12)') Tracer3%valX(NLMAX)%x(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.12)') Tracer3%valY(NLMAX)%x(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.12)') Tracer3%valZ(NLMAX)%x(i)
 END DO

 CLOSE(2)

END SUBROUTINE dump_LinScalar_out_Q1
!
! ----------------------------------------------
!
SUBROUTINE dump_LinScalar_in_Q1(cdump)
USE var_QuadScalar, only:knvt
INTEGER IFl,itwx,i
CHARACTER CIFile*15,cdump*(*)

 IF (myid.EQ.0) RETURN

 IFl=LEN(TRIM(ADJUSTL(cdump)))
 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.lt.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I1,I1)') '_',0,myid
 IF (myid.ge.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I2)') '_',myid
 OPEN (UNIT=1,FILE=CIFile,STATUS="OLD",FORM="FORMATTED")

 DO I=1,KNVT(NLMAX)
  READ(1,*) Tracer3%valX(NLMAX)%x(i)
 END DO
 DO I=1,KNVT(NLMAX)
  READ(1,*) Tracer3%valY(NLMAX)%x(i)
 END DO
 DO I=1,KNVT(NLMAX)
  READ(1,*) Tracer3%valZ(NLMAX)%x(i)
 END DO

 CLOSE(1)

END SUBROUTINE dump_LinScalar_in_Q1
!
