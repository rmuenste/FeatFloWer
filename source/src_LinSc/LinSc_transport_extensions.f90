!
! ----------------------------------------------
!
SUBROUTINE Transport_Q1_displacement(mfile,INL)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
INTEGER mfile,INL
REAL*8  ResTemp(3),DefTemp(3),DefTempCrit(3),RhsTemp(3)
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

NLMAX = NLMAX + 1

thstep = 0d0*tstep

IF (myid.ne.0) THEN

! advect the scalar field
CALL Create_NewDiffMat_Q1(myALE%Q2Coor_old,Properties%DiffCoeff(1))

! Set dirichlet boundary conditions on the solution
CALL Boundary_LinSc_Val_Q1()

! Assemble the defect vector and fine level matrix
IF (ADJUSTL(TRIM(Tracer3%prm%cEquation)).eq.'Deformation') THEN
 CALL Matdef_DeformationTensor_LinScalar(myALE%Q2Coor_old,Tracer3,-1,1,Properties%DiffCoeff(1))
ELSE
 CALL Matdef_Laplace_LinScalar(myALE%Q2Coor_old,Tracer3,-1,1,Properties%DiffCoeff(1))
END IF

CALL E011Sum(Tracer3%defX)
CALL E011Sum(Tracer3%defY)
CALL E011Sum(Tracer3%defZ)

! Set dirichlet boundary conditions on the defect
CALL Boundary_LinSc_Def_Q1()

!! Save the old solution
CALL LCP1(Tracer3%valX,Tracer3%valX_old,Tracer3%ndof)
CALL LCP1(Tracer3%valY,Tracer3%valY_old,Tracer3%ndof)
CALL LCP1(Tracer3%valZ,Tracer3%valZ_old,Tracer3%ndof)

!! Compute the defect
CALL Resdfk_General_LinScalar_Q1(Tracer3,ResTemp,DefTemp,RhsTemp)

end if

CALL COMM_Maximum(RhsTemp(1))
CALL COMM_Maximum(RhsTemp(2))
CALL COMM_Maximum(RhsTemp(3))

DefTempCrit(1)=MAX((RhsTemp(1))*Tracer3%prm%defCrit,Tracer3%prm%MinDef)
DefTempCrit(2)=MAX((RhsTemp(2))*Tracer3%prm%defCrit,Tracer3%prm%MinDef)
DefTempCrit(3)=MAX((RhsTemp(3))*Tracer3%prm%defCrit,Tracer3%prm%MinDef)


do INL=1,Tracer3%prm%NLmax
INLComplete = 0

! Calling the solver
CALL Solve_General_MGLinScalar(Tracer3,&
                               Boundary_LinSc_Val_Q1,&
                               Boundary_LinSc_XYZMat,&
                               mfile)

IF (myid.ne.0) THEN

! Restore the constant right hand side
 Tracer3%defX = 0d0 !Tracer%rhs
 Tracer3%defY = 0d0 !Tracer%rhs
 Tracer3%defZ = 0d0 !Tracer%rhs

! Assemble the defect vector and fine level matrix
 IF (ADJUSTL(TRIM(Tracer3%prm%cEquation)).eq.'Deformation') THEN
  CALL Matdef_DeformationTensor_LinScalar(myALE%Q2Coor_old,Tracer3,-1,0,Properties%DiffCoeff(1))
 ELSE
  CALL Matdef_Laplace_LinScalar(myALE%Q2Coor_old,Tracer3,-1,0,Properties%DiffCoeff(1))
 END IF

 CALL E011Sum(Tracer3%defX)
 CALL E011Sum(Tracer3%defY)
 CALL E011Sum(Tracer3%defZ)

 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def_Q1()

! Save the old solution
 CALL LCP1(Tracer3%valX,Tracer3%valX_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valY,Tracer3%valY_old,Tracer3%ndof)
 CALL LCP1(Tracer3%valZ,Tracer3%valZ_old,Tracer3%ndof)

! Compute the defect
 CALL Resdfk_General_LinScalar_Q1(Tracer3,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
CALL COMM_Maximum(RhsTemp(1))
CALL COMM_Maximum(RhsTemp(2))
CALL COMM_Maximum(RhsTemp(3))

!CALL Protocol_linScalar(mfile,Tracer,INL,&
!     ResTemp,DefTemp,RhsTemp)

IF ((DefTemp(1).LE.DefTempCrit(1)).AND.&
    (DefTemp(2).LE.DefTempCrit(2)).AND.&
    (DefTemp(3).LE.DefTempCrit(3)).AND.&
    (INL.GE.Tracer%prm%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO

1 CONTINUE

myQ2coor(1,:) = Tracer3%valX
myQ2coor(2,:) = Tracer3%valY
myQ2coor(3,:) = Tracer3%valZ

2 CONTINUE

NLMAX = NLMAX - 1

CALL CoommunicateCoarseGrid()

CALL GetMeshVelocity()

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
SUBROUTINE Init_Disp_Q1(log_unit)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE var_QuadScalar
implicit none
integer, intent(in) :: log_unit
integer :: n, ndof


NLMAX = NLMAX + 1

 ! Building up the matrix strucrures
 CALL Create_MatStruct_Q1()

! Iteration matrix (only allocation)
 CALL Create_AMat_Q1()

 CALL Initialize_Q1(Tracer3)

 ILEV=NLMAX
 CALL SETLEV(2)

! Set the types of boundary conditions (set up knpr)
 call LinSc_Knpr_Q1(mg_mesh%level(ilev)%dcorvg)

 ! Prolongation Restriction Matrices
 ALLOCATE(mg_E011Prol(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011Rest(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011ProlM(NLMIN:NLMAX-1))
 ALLOCATE(mg_E011RestM(NLMIN:NLMAX-1))

 DO ILEV=NLMIN,NLMAX-1
  N = KNVT(ILEV+1) + 2*KNET(ILEV+1) + 4*KNAT(ILEV+1) + 8*KNEL(ILEV+1)
  ALLOCATE(mg_E011Prol(ILEV)%a(N))
  ALLOCATE(mg_E011Rest(ILEV)%a(N))
  ALLOCATE(mg_E011ProlM(ILEV)%ColA(N))
  ALLOCATE(mg_E011RestM(ILEV)%ColA(N))
  mg_E011ProlM(ILEV)%na = N
  mg_E011RestM(ILEV)%na = N
  NDOF = KNVT(ILEV+1)!+KNET(ILEV+1)+KNAT(ILEV+1)+KNEL(ILEV+1)
  mg_E011ProlM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011ProlM(ILEV)%LdA(NDOF+1))
  NDOF = KNVT(ILEV)!+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)
  mg_E011RestM(ILEV)%nu = NDOF
  ALLOCATE(mg_E011RestM(ILEV)%LdA(NDOF+1))
 END DO

CALL InitializeProlRest(Tracer3)


Tracer3%cName = "Disp"

CALL GetDispParameters(Tracer3%prm,Tracer3%cName,log_unit)

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

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr_Q1(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
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

 IF (myBoundary%LS_zero(i).ne.0.or.myBoundary%bDisp_DBC(i)) THEN
   Tracer3%knprX(I) = 1
   Tracer3%knprY(I) = 1
   Tracer3%knprZ(I) = 1
 END IF

 IF (myBoundary%bSymmetry(1,i)) THEN
   Tracer3%knprX(I) = 1
 END IF

 IF (myBoundary%bSymmetry(2,i)) THEN
   Tracer3%knprY(I) = 1
 END IF

 IF (myBoundary%bSymmetry(3,i)) THEN
   Tracer3%knprZ(I) = 1
 END IF

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

 Tracer3%valX(i) = X
 Tracer3%valY(i) = Y
 Tracer3%valZ(i) = Z
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
SUBROUTINE Boundary_LinSc_Val_Q1()
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0,RY = 0.0d0
REAL*8 :: RADx = 0.2d0
REAL*8 DIST
INTEGER i

DO i=1,Tracer3%ndof

 IF (Tracer3%knprX(i).eq.1) THEN
  Tracer3%valX(i) = myQ2Coor(1,i)
 END IF

 IF (Tracer3%knprY(i).eq.1) THEN
  Tracer3%valY(i) = myQ2Coor(2,i)
 END IF
 
 IF (Tracer3%knprZ(i).eq.1) THEN
  Tracer3%valZ(i) = myQ2Coor(3,i)
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
SUBROUTINE InitializeProlRest(lSc)
implicit none
TYPE(lScalar3), INTENT(INOUT), TARGET :: lSc

 MyMG%MinLev  = NLMIN
 myMG%MedLev  = NLMIN
 myMG%MaxLev  = NLMAX

 IF(myid.eq.showid) WRITE(*,*) "Initialization of displacement prolongation matrix"
 myMG%bProlRest => lSc%bProlRest
 MyMG%cVariable = "Displacement"
 CALL mgProlRestInit()

END SUBROUTINE InitializeProlRest  
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
  WRITE(2,'(G18.12)') Tracer3%valX(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.12)') Tracer3%valY(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.12)') Tracer3%valZ(i)
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
  READ(1,*) Tracer3%valX(i)
 END DO
 DO I=1,KNVT(NLMAX)
  READ(1,*) Tracer3%valY(i)
 END DO
 DO I=1,KNVT(NLMAX)
  READ(1,*) Tracer3%valZ(i)
 END DO

 CLOSE(1)

END SUBROUTINE dump_LinScalar_in_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_XYZMat(DA11,DA22,DA33,KLD,&
           KNPRU,KNPRV,KNPRW,NDOF)
REAL*8  DA11(*),DA22(*),DA33(*)
INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF

DO I=1,NDOF
 IF (KNPRU(I).EQ.1) THEN
   ICOL = KLD(I)
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
   END DO
 END IF
 IF (KNPRV(I).EQ.1) THEN
   ICOL = KLD(I)
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA22(ICOL) = 0d0
   END DO
 END IF
 IF (KNPRW(I).EQ.1) THEN
   ICOL = KLD(I)
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA33(ICOL) = 0d0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_LinSc_XYZMat

