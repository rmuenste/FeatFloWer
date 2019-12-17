SUBROUTINE Transport_LinScalar_General(sub_BC,sub_SRC,mfile,INL)
INTEGER mfile,INL
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

EXTERNAL sub_BC,sub_SRC

NLMAX = NLMAX + 1

thstep = 0.5d0*tstep

! advect the scalar field
IF (myid.ne.0) THEN

 Tracer%oldSol = Tracer%val(NLMAX)%x
 
 IF (Tracer%prm%AFC) CALL InitAFC_General_LinScalar()

! Assemble the right hand side
 CALL LCL1(Tracer%def,Tracer%ndof)
 CALL Matdef_LinScalar(Tracer,1,0)

 ! Add the source term to the RHS
 CALL sub_SRC()

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
CALL Solve_General_LinScalar(Tracer,ParKNPR,sub_BC,Boundary_LinSc_Mat)

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

if (myid.ne.master) then
 Temperature = Tracer%val(NLMAX)%x
end if

NLMAX = NLMAX - 1


END SUBROUTINE Transport_LinScalar_General
!
! ----------------------------------------------
!
SUBROUTINE Transport_LinScalar_XSE(sub_BC,sub_SRC,mfile,INL)
INTEGER mfile,INL
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

EXTERNAL sub_BC,sub_SRC

NLMAX = NLMAX + 1

thstep = 0.5d0*tstep

! advect the scalar field
IF (myid.ne.0) THEN

 Tracer%oldSol = Tracer%val(NLMAX)%x
 
 IF (Tracer%prm%AFC) CALL InitAFC_General_LinScalar()

! Assemble the right hand side
 CALL LCL1(Tracer%def,Tracer%ndof)
 CALL Matdef_LinScalar(Tracer,1,0)

 ! Add the source term to the RHS
 CALL sub_SRC()

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
CALL Solve_General_LinScalar(Tracer,ParKNPR,sub_BC,Boundary_LinSc_Mat)

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

if (myid.ne.master) then
 Temperature = Tracer%val(NLMAX)%x
end if

NLMAX = NLMAX - 1

CALL COMM_SUMM(dIntegralHeat)
IF (myid.eq.1) WRITE(*,*) 'IntegralHeatRate_[W]: ',dIntegralHeat*(myThermodyn%density*myThermodyn%cp)
IF (myid.eq.1) WRITE(mfile,*) 'IntegralHeatRate_[W]: ',dIntegralHeat*(myThermodyn%density*myThermodyn%cp)

ilev=nlmax
call setlev(2)
CALL IntegrateOutflowTemp(mg_mesh%level(ilev)%dcorvg,&
                          mg_mesh%level(ilev)%karea,&
                          mg_mesh%level(ilev)%kvert,&
                          mg_mesh%level(ilev)%nel,&
                          mfile)

END SUBROUTINE Transport_LinScalar_XSE
!
! ----------------------------------------------
!
SUBROUTINE Transport_LinScalar_EWIKON(sub_BC,sub_SRC,mfile,INL)
INTEGER mfile,INL
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
REAL*8 tstep_old,thstep_old
INTEGER INLComplete,I,J

EXTERNAL sub_BC,sub_SRC

NLMAX = NLMAX + 1

thstep = 0.5d0*tstep

! advect the scalar field
IF (myid.ne.0) THEN

 Tracer%oldSol = Tracer%val(NLMAX)%x
 
 ! Convection matrix
 CALL Create_RhoCpConvMat(QuadSc%valU,QuadSc%valV,QuadSc%valW)
! CALL Build_LinSc_Convection()
 IF (Tracer%prm%AFC) CALL InitAFC_General_LinScalar()

! Assemble the right hand side
 CALL LCL1(Tracer%def,Tracer%ndof)
 CALL Matdef_LinScalar_EWIKON(Tracer,1,0)

! Add the source term to the RHS
 CALL sub_SRC()

! Add the boundary heat flux (explicit part)
 CALL AddBoundaryHeatFlux(1)
 
 ! Set dirichlet boundary conditions on the defect
 CALL Boundary_LinSc_Def()

 ! Store the constant right hand side
 Tracer%rhs = Tracer%def

! Assemble the defect vector and fine level matrix
 CALL Matdef_LinScalar_EWIKON(Tracer,-1,1)

! Add the boundary heat flux (implicit part)
 CALL AddBoundaryHeatFlux(2)
   
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
sub_BC,Boundary_LinSc_Mat)

!CALL myMPI_barrier()
!write(*,*) 'sdf sdf  s  fsdf sd fs  f',myid
!pause

IF (myid.ne.0) THEN

!!!!          Checking the quality of the result           !!!!
! Restore the constant right hand side
 Tracer%def = Tracer%rhs

! Assemble the defect vector and fine level matrix
 CALL Matdef_LinScalar_EWIKON(Tracer,-1,0)

 ! Add the boundary heat flux
!  CALL  AddHeatFlux(dArea,dFlux)

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

CALL COMM_SUMM(dArea1)
CALL COMM_SUMM(dArea2)
CALL COMM_SUMM(dFlux1)
CALL COMM_SUMM(dFlux2)

if (myid.eq.showid) then
 WRITE(mterm,'(A,8ES12.4)') 'Time[s]_Area[cm2]_HeatFlux[kW]_HeatSource[kW]: ',timens, dArea1,1e-3*dFlux1, dArea2,1e-3*dFlux2,dHeatSource
 WRITE(mfile,'(A,8ES12.4)') 'Time[s]_Area[cm2]_HeatFlux[kW]_HeatSource[kW]: ',timens, dArea1,1e-3*dFlux1, dArea2,1e-3*dFlux2,dHeatSource
end if

CALL IntegrateOutputQuantities(mfile)

if (myid.ne.master) then
 Temperature = Tracer%val(NLMAX)%x
end if

NLMAX = NLMAX - 1


END SUBROUTINE Transport_LinScalar_EWIKON
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

myALE%Q2Coor_old = myQ2Coor

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

CALL Protocol_linScalar_Disp_Q1(mfile,Tracer3,0,&
       ResTemp,DefTemp,DefTempCrit," Laplace equation ")

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
CALL Protocol_linScalar_Disp_Q1(mfile,Tracer3,INL,&
       ResTemp,DefTemp,DefTempCrit," Laplace equation ")

IF ((DefTemp(1).LE.DefTempCrit(1)).AND.&
    (DefTemp(2).LE.DefTempCrit(2)).AND.&
    (DefTemp(3).LE.DefTempCrit(3)).AND.&
    (INL.GE.Tracer3%prm%NLmin)) INLComplete = 1

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
 WRITE(*,'(A8,6XA10,5(6X,ES10.3))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(MFILE,'(A8,6XA10,5(6X,ES10.3))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(*,5)
 WRITE(MFILE,5)
 WRITE(*,'(I8,5(6X,ES10.3))') 0,ResScalar,DefScalar
 WRITE(MFILE,'(I8,5(6X,ES10.3))') 0,ResScalar,DefScalar
ELSE
 WRITE(*,'(I8,5(6X,ES10.3))') nINL,ResScalar,DefScalar,RhsScalar
 WRITE(MFILE,'(I8,5(6X,ES10.3))') nINL,ResScalar,DefScalar,RhsScalar
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

 IF (myBoundary%LS_zero(i).ne.0) THEN
  Tracer3%valY(i) = myQ2Coor(2,i)+0.5
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
SUBROUTINE Build_LinSc_Convection_Q1(U,V,W)
REAL*8 U(*),V(*),W(*)
INTEGER I,J,jLEV
EXTERNAL E011

ILEV=NLMAX
JLEV = ILEV-1
CALL SETLEV(2)

Kmat = 0d0
IF (myid.eq.showid) write(*,*) 'Regenerating K Matrix for Q1'

CALL Conv_LinSc2(U,V,W,Kmat,&
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
  WRITE(2,'(G18.11)') Tracer3%valX(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.11)') Tracer3%valY(i)
 END DO
 DO I=1,KNVT(NLMAX)
  WRITE(2,'(G18.11)') Tracer3%valZ(i)
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
!
! ----------------------------------------------
!
SUBROUTINE AddBoundaryHeatFlux(iSwitch)
INTEGER iSwitch
EXTERNAL E011


if (myid.ne.master) then
 ilev = NLMAX
 call setlev(2)
!  CALL AddBoundaryHeatFluxSub(Amat,lMat%LdA,lMat%ColA,&
!                              Tracer%def,Tracer%oldSol,&
!                              mg_mesh%level(ilev)%kvert,&
!                              mg_mesh%level(ilev)%karea,&
!                              mg_mesh%level(ilev)%kedge,&
!                              mg_mesh%level(ilev)%dcorvg,&
!                              E011,dArea1,dFlux1,iSwitch)
 CALL AddConductiveHeatFluxSub(Amat,lMat%LdA,lMat%ColA,&
                             Tracer%def,Tracer%oldSol,&
                             mg_mesh%level(ilev)%kvert,&
                             mg_mesh%level(ilev)%karea,&
                             mg_mesh%level(ilev)%kedge,&
                             mg_mesh%level(ilev)%dcorvg,&
                             E011,dArea2,dFlux2,iSwitch)
end if

END SUBROUTINE AddBoundaryHeatFlux
!
! ----------------------------------------------
!
SUBROUTINE IntegrateOutputQuantities(mfile)
EXTERNAL E011
REAL*8 dQuant(12),dSensorTemperature(2),P(3),Q(3),dist
REAL*8 :: dTotalEnthalpy(2)=0d0
integer mfile,iS,iSeg,i

Q = myProcess%TemperatureSensorCoor

if (myid.ne.master) then
 ilev = NLMAX
 call setlev(2)
 dSensorTemperature = 0d0
 do i=1,mg_mesh%level(ilev)%nvt
  P = mg_mesh%level(ilev)%dcorvg(:,i)
  dist = SQRT((P(1)-Q(1))**2d0 + (P(2)-Q(2))**2d0 + (P(3)-Q(3))**2d0) 
  IF (dist.le.myProcess%TemperatureSensorRadius) then
   dSensorTemperature(1) = dSensorTemperature(1) + MLmat(i)*Tracer%oldSol(i)
   dSensorTemperature(2) = dSensorTemperature(2) + MLmat(i)
  END IF
 end do

 dTotalEnthalpy(1) = 0d0
 do i=1,mg_mesh%level(ilev)%nvt
   dTotalEnthalpy(1) = dTotalEnthalpy(1) + MLRhoCpMat(i)*Tracer%oldSol(i)*1e-10
 end do

 ilev = NLMAX
 call setlev(2)
 dQuant = 0d0
 CALL IntegrateOutputQuantitiesSub(Tracer%oldSol,&
                                   mg_mesh%level(ilev)%kvert,&
                                   mg_mesh%level(ilev)%karea,&
                                   mg_mesh%level(ilev)%kedge,&
                                   mg_mesh%level(ilev)%dcorvg,&
                                   E011,dQuant)
end if

CALL Comm_SummN(dTotalEnthalpy,1)
CALL Comm_SummN(dSensorTemperature,2)
CALL Comm_SummN(dQuant,12)

if (itns.eq.1) then
 dTotalEnthalpy(2) = dTotalEnthalpy(1)
end if

IF (myid.eq.1) then
 write(MTERM,'(A,20ES12.4)') 'SensorTemperature_t[s]_Tsens[C]_Vsens[cm3]_dH[kJ]: ',timens,dSensorTemperature(1)/dSensorTemperature(2),dSensorTemperature(2),dTotalEnthalpy(2)-dTotalEnthalpy(1)
 write(MFILE,'(A,20ES12.4)') 'SensorTemperature_t[s]_Tsens[C]_Vsens[cm3]_dH[kJ]: ',timens,dSensorTemperature(1)/dSensorTemperature(2),dSensorTemperature(2),dTotalEnthalpy(2)-dTotalEnthalpy(1)
 iS = 0
 write(MTERM,'(A,20ES12.4)') 'IntQuantBLOCK_t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
 write(MFILE,'(A,20ES12.4)') 'IntQuantBLOCK_t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
 iS = 4
 write(MTERM,'(A,20ES12.4)') 'IntQuantWIRE__t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
 write(MFILE,'(A,20ES12.4)') 'IntQuantWIRE__t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
 iS = 8
 write(MTERM,'(A,20ES12.4)') 'IntQuantMELT__t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
 write(MFILE,'(A,20ES12.4)') 'IntQuantMELT__t[s]_A[cm2]_V[cm2]_Ta[C]_Tv[C]: ',timens,dQuant(iS+1),dQuant(iS+3),dQuant(iS+2)/dQuant(iS+1),dQuant(iS+4)/dQuant(iS+3)
end if

iS = 0
iSeg = 0
dHeatSource = 0d0

do 
 iSeg = iSeg + 1
 
 mySigma%mySegment(iSeg)%UseHeatSource = 0d0
 
 if (iSeg.gt.mySigma%NumberOfSeg) exit 
 IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
!   IF (dQuant(iS+4)/dQuant(iS+3).gt.291d0) then
  IF (dSensorTemperature(1)/dSensorTemperature(2).gt.291d0) then
   mySigma%mySegment(iSeg)%UseHeatSource  =  mySigma%mySegment(iSeg)%HeatSourceMin
  END IF
  IF (dSensorTemperature(1)/dSensorTemperature(2).lt.289d0) then
!  IF (dQuant(iS+4)/dQuant(iS+3).lt.289d0) then
   mySigma%mySegment(iSeg)%UseHeatSource  =  mySigma%mySegment(iSeg)%HeatSourceMax
  END IF
  dHeatSource = mySigma%mySegment(iSeg)%UseHeatSource
  GOTO 1
 END IF
end do
1 continue

return

END SUBROUTINE IntegrateOutputQuantities
!
! ----------------------------------------------
!
SUBROUTINE Assemble_LinScOperators_XSE(mfile)
USE var_QuadScalar, ONLY : Screw,Viscosity,Shearrate

integer mfile
REAL*8 myDiffCoeff_melt,myDiffCoeff_steel,daux,dHeat
REAL*8, ALLOCATABLE ::  AlphaDiff(:)
INTEGER iel,ivt,i

dHeat = 0d0

if (myid.ne.0) then 

  NLMAX = NLMAX + 1
  
  ! Convection matrix
  CALL Build_LinSc_Convection_Q1(QuadSc%valU,QuadSc%valV,QuadSc%valW)

  ! Diffusion matrix 
!  myDiffCoeff = Properties%DiffCoeff(1)
  !      W         cm3        (k)g . K           J         cm3   1       1        cm3     1     cm2 . cm      1    cm2
  ! ------------* ----- * ------------ =  ------------* ----- * ---- = --------* ----- * --- =  ---------- = ----* ----
  !    m . K        g          (k)J           s . m         1    J      s . m      1      1       100cm .s    100    s
 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!   -     ---   Diffusion Matrix Setup --- -    !!!
  myDiffCoeff_melt = 0.01d0*myThermodyn%lambda/(myThermodyn%density*myThermodyn%cp)
  myDiffCoeff_steel = 1E-1 ! cm2/s
  if (myid.eq.1)  write(*,'(A,2ES12.4)') ' ThermalDiffCoeff_steel_/_melt_[cm2/s]: ',myDiffCoeff_steel,myDiffCoeff_melt
  
  ALLOCATE(AlphaDiff(mg_mesh%level(NLMAX)%nel))
  AlphaDiff = 0d0
  
  DO iel = 1,mg_mesh%level(NLMAX)%nel
   do ivt=1,8
    i = mg_mesh%level(NLMAX)%kvert(ivt,iel)
    IF (Screw(i).ge.0d0) THEN
     AlphaDiff(iel) = AlphaDiff(iel) + 0.125d0*myDiffCoeff_melt
    ELSE
     AlphaDiff(iel) = AlphaDiff(iel) + 0.125d0*myDiffCoeff_steel
    END IF
   END DO
  END DO
  
  CALL Create_XSE_DiffMat(AlphaDiff)
  
  DEALLOCATE(AlphaDiff)
  !!!   -     ---   Diffusion Matrix Setup --- -    !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Mass matrix
  CALL Create_MassMat()

 ! Lumped Mass matrix
  CALL Create_LMassMat()

  NLMAX = NLMAX - 1

end if

!   pause
END SUBROUTINE Assemble_LinScOperators_XSE
!
! ----------------------------------------------
!
SUBROUTINE Assemble_LinScOperators_General()
USE var_QuadScalar, ONLY : distance

REAL*8 myDiffCoeff_melt,myDiffCoeff_steel
REAL*8, ALLOCATABLE ::  AlphaDiff(:)
INTEGER iel,ivt,i

if (myid.ne.0) then 

  NLMAX = NLMAX + 1
  
  ! Convection matrix
  CALL Build_LinSc_Convection_Q1(QuadSc%valU,QuadSc%valV,QuadSc%valW)

  ! Diffusion matrix 
!  myDiffCoeff = Properties%DiffCoeff(1)
  !      W         cm3        (k)g . K           J         cm3   1       1        cm3     1     cm2 . cm      1    cm2
  ! ------------* ----- * ------------ =  ------------* ----- * ---- = --------* ----- * --- =  ---------- = ----* ----
  !    m . K        g          (k)J           s . m         1    J      s . m      1      1       100cm .s    100    s
 
  myDiffCoeff_melt = 0.01d0*myThermodyn%lambda/(myThermodyn%density*myThermodyn%cp)
  myDiffCoeff_steel = 5E-2 ! cm2/s
  if (myid.eq.1)  write(*,'(A,2ES12.4)') ' ThermalDiffCoeff_steel_/_melt_[cm2/s]: ',myDiffCoeff_steel,myDiffCoeff_melt
  
  ALLOCATE(AlphaDiff(mg_mesh%level(NLMAX)%nel))
  AlphaDiff = 0d0
  
  DO iel = 1,mg_mesh%level(NLMAX)%nel
   do ivt=1,8
    i = mg_mesh%level(NLMAX)%kvert(ivt,iel)
    IF (distance(i).ge.0d0) THEN
     AlphaDiff(iel) = AlphaDiff(iel) + 0.125d0*myDiffCoeff_melt
    ELSE
     AlphaDiff(iel) = AlphaDiff(iel) + 0.125d0*myDiffCoeff_steel
    END IF
   END DO
  END DO
  
!   write(*,*) 'size(AlphaDiff)=',size(AlphaDiff)
!  CALL Create_ConstDiffMat(myDiffCoeff_melt)  
  CALL Create_DIE_DiffMat(AlphaDiff)

  DEALLOCATE(AlphaDiff)
  
  ! Mass matrix
  CALL Create_MassMat()

 ! Lumped Mass matrix
  CALL Create_LMassMat()
  
  NLMAX = NLMAX - 1

end if

END SUBROUTINE Assemble_LinScOperators_General
!
! ----------------------------------------------
!
SUBROUTINE IntegrateOutflowTemp(dcorvg,karea,kvert,nel,mfile)
INTEGER mfile
REAL*8 dcorvg(3,*),T_Avg
INTEGER karea(6,*),kvert(8,*),nel
!---------------------------------
INTEGER NeighA(4,6)
REAL*8 P(3),dA,dVolFlow,dVolFlowT,dArea
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
LOGICAL BB(4)

if (myid.ne.0) then
dVolFlow = 0d0
dVolFlowT = 0d0
dArea = 0d0

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   
   BB(1) = (myBoundary%bSymmetry(1,ivt1).and.myBoundary%bSymmetry(2,ivt1))
   BB(2) = (myBoundary%bSymmetry(1,ivt2).and.myBoundary%bSymmetry(2,ivt2))
   BB(3) = (myBoundary%bSymmetry(1,ivt3).and.myBoundary%bSymmetry(2,ivt3))
   BB(4) = (myBoundary%bSymmetry(1,ivt4).and.myBoundary%bSymmetry(2,ivt4))
   
   IF (BB(1).and.BB(2).and.BB(3).and.BB(4)) THEN
!    IF (abs(dcorvg(3,ivt1)-1.0d0*mySigma%L).lt.1d-3.and.&
!        abs(dcorvg(3,ivt2)-1.0d0*mySigma%L).lt.1d-3.and.&
!        abs(dcorvg(3,ivt3)-1.0d0*mySigma%L).lt.1d-3.and.&
!        abs(dcorvg(3,ivt4)-1.0d0*mySigma%L).lt.1d-3) THEN

       CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
       dArea = dArea + DABS(dA)
       dVolFlow  = dVolFlow  + dA*QuadSc%ValW(nvt+net+k)
       dVolFlowT = dVolFlowT + dA*QuadSc%ValW(nvt+net+k)*Tracer%val(NLMAX+1)%x(nvt+net+k)
   END IF
   k = k + 1
  END IF
 END DO
END DO
end if

CALL COMM_SUMM(dVolFlow)
CALL COMM_SUMM(dVolFlowT)
CALL COMM_SUMM(dArea)

if (myid.ne.0) then
 T_Avg = dVolFlowT/dVolFlow
end if

IF (myid.eq.1) WRITE(MTERM,'(A,3ES14.6)') 'OutflowTemp_[C]_Area_[mm2]_VolFlow_[l/h]: ',T_Avg,dArea,3.6d0*dVolFlow
IF (myid.eq.1) WRITE(MFILE,'(A,3ES14.6)') 'OutflowTemp_[C]_Area_[mm2]_VolFlow_[l/h]: ',T_Avg,dArea,3.6d0*dVolFlow


END SUBROUTINE IntegrateOutflowTemp
!
! ----------------------------------------------
!

!  CALL AddLumpedHeatFlux(mg_mesh%level(nlmax)%dcorvg,&
!                         mg_mesh%level(nlmax)%karea,&
!                         mg_mesh%level(nlmax)%kvert,&
!                         mg_mesh%level(nlmax)%nvt,&
!                         mg_mesh%level(nlmax)%nel,&
!                         mg_mesh%level(nlmax)%net,1/2)
! SUBROUTINE AddLumpedHeatFlux(dcorvg,karea,kvert,nvt,nel,net,iSwitch)
! REAL*8 dcorvg(3,*)
! INTEGER karea(6,*),kvert(8,*),nel,net,nvt
! INTEGER iSwitch
! !---------------------------------
! INTEGER NeighA(4,6)
! REAL*8 P(3),dA,tLocal
! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
! INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
! 
! if (iSwitch.eq.1) then
!  dFlux = 0d0
!  dArea = 0d0
! end if
! 
! k=1
! DO i=1,nel
!  DO j=1,6
!   IF (k.eq.karea(j,i)) THEN
!    ivt1 = (NeighA(1,j),i)
!    ivt2 = kvert(NeighA(2,j),i)
!    ivt3 = kvert(NeighA(3,j),i)
!    ivt4 = kvert(NeighA(4,j),i)
!    IF (myBoundary%iTemperature(ivt1).eq.2.and. &
!        myBoundary%iTemperature(ivt2).eq.2.and. &
!        myBoundary%iTemperature(ivt3).eq.2.and. &
!        myBoundary%iTemperature(ivt4).eq.2) THEN
!        CALL GET_area(dcorvg(1:3,ivt1),dcorvg(1:3,ivt2),dcorvg(1:3,ivt3),dcorvg(1:3,ivt4),dA)
!        if (iSwitch.eq.1) then
!         TLocal = 0.25d0*(Tracer%val(NLMAX)%x(ivt1) + Tracer%val(NLMAX)%x(ivt2) + &
!                          Tracer%val(NLMAX)%x(ivt3) + Tracer%val(NLMAX)%x(ivt4))
!         Tracer%def(ivt1) = Tracer%def(ivt1) + 1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP*myProcess%AirTemperature
!         Tracer%def(ivt2) = Tracer%def(ivt2) + 1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP*myProcess%AirTemperature
!         Tracer%def(ivt3) = Tracer%def(ivt3) + 1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP*myProcess%AirTemperature
!         Tracer%def(ivt4) = Tracer%def(ivt4) + 1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP*myProcess%AirTemperature
!         dFlux = dFlux     + myProcess%HeatTransferCoeff*(dA*1e-4)*(myProcess%AirTemperature-tLocal)
!         dArea = dArea     + dA
!        end if
!        if (iSwitch.eq.2) then
!         Amat(lMat%LdA(ivt1)) = Amat(lMat%LdA(ivt1)) + REAL(1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP)
!         Amat(lMat%LdA(ivt2)) = Amat(lMat%LdA(ivt2)) + REAL(1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP)
!         Amat(lMat%LdA(ivt3)) = Amat(lMat%LdA(ivt3)) + REAL(1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP)
!         Amat(lMat%LdA(ivt4)) = Amat(lMat%LdA(ivt4)) + REAL(1d3*myProcess%HeatTransferCoeff*(0.25d0*dA)*TSTEP)
!        end if
!    END IF
!    k = k + 1
!   END IF
!  END DO
! END DO
! 
! END SUBROUTINE AddLumpedHeatFlux
