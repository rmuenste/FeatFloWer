MODULE Transport_Q1

USE def_LinScalar
USE PP3D_MPI, ONLY:E011Sum,E011Knpr,Comm_NLComplete,&
    Comm_Maximum,Comm_Summ,myid,master,CommSum,Comm_SummN,myMPI_barrier
USE Transport_Q2P1, ONLY: QuadSc,ParKNPR,mgDiffCoeff,&
    myBoundary,myQ2Coor,&
    MoveInterfacePoints,myALE,Properties,getmeshvelocity
USE var_QuadScalar, ONLY: myMG,myHeatObjects
USE mg_LinScalar, ONLY : mgProlRestInit
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,MyMaterials

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

!CALL myMPI_barrier()
!write(*,*) 'sdf sdf  s  fsdf sd fs  f',myid
!pause

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
 CALL Create_DiffMat(mgDiffCoeff(NLMAX)%x)

! Iteration matrix (only allocation)
 CALL Create_AMat()

! Initialize the scalar quantity
 CALL Initialize(Tracer)

! Set the types of boundary conditions (set up knpr)
 CALL Create_Knpr(LinSc_Knpr)

END IF


Tracer%cName        = "Tracer"
Tracer%prm%SolvIter = 1
Tracer%prm%AFC      = .TRUE.
Tracer%prm%NLmin    = 2
Tracer%prm%NLmax    =20
Tracer%prm%defCrit  =1d-6
Tracer%prm%epsCrit  =1d-3
Tracer%prm%MinDef   =1d-16
Tracer%prm%SolvType =1

NLMAX = NLMAX - 1

END SUBROUTINE Init_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE InitHeatObjects()

ALLOCATE(myHeatObjects%Block(Tracer%ndof))
ALLOCATE(myHeatObjects%Wire(Tracer%ndof))
ALLOCATE(myHeatObjects%Segment(Tracer%ndof))
myHeatObjects%Block = 0d0
myHeatObjects%Wire  = 0d0
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
SUBROUTINE InitCond_GeneralLinScalar(sub_IC,sub_BC)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
EXTERNAL sub_IC,sub_BC

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

CALL sub_IC(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL sub_BC()

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_GeneralLinScalar
!
! ----------------------------------------------
!
SUBROUTINE LinSc_Knpr(dcorvg)
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.2d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 Tracer%knpr(I) = 0

 IF (myBoundary%iTemperature(i).GT.0) THEN
  Tracer%knpr(I) = myBoundary%iTemperature(i)
 END IF

END DO

END SUBROUTINE LinSc_Knpr
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
SUBROUTINE LinSc_InitCond(dcorvg)
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0, RY = 0.0d0, RZ = 2.4d0
REAL*8 :: RADx = 0.20d0,RADs=0.040
REAL*8 DIST
INTEGER i,iSeg,iMat

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 Tracer%val(NLMAX)%x(i) = 0d0

END DO

END SUBROUTINE LinSc_InitCond
!
! ----------------------------------------------
!
SUBROUTINE LinSc_InitCond_EWIKON(dcorvg)
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0, RY = 0.0d0, RZ = 2.4d0
REAL*8 :: RADx = 0.20d0,RADs=0.040
REAL*8 DIST
INTEGER i,iSeg,iMat

DO i=1,Tracer%ndof
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 iSeg = myHeatObjects%Segment(i)
 
 IF (iSeg.eq.0) Tracer%val(NLMAX)%x(i) = myProcess%AirTemperature
 IF (iSeg.ne.0) then
  Tracer%val(NLMAX)%x(i) = mySigma%mySegment(iSeg)%InitTemp
 END IF

END DO

END SUBROUTINE LinSc_InitCond_EWIKON
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
SUBROUTINE Boundary_LinSc_Val_Weber()
REAL*8 X,Y,Z
REAL*8 :: x0 = -0d0, y0 = -22.1d0, r0=3.0d0
REAL*8 XT,YT,ZT,xR,yR,Tx,Ty
REAL*8 :: PI=4d0*DATAN(1d0),dAlpha
REAL*8 :: ProfXT(11,2),ProfYT(11,2)
REAL*8 dist,frac
INTEGER i,l
DATA ProfXT/0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,&
            203d0,203d0,193d0,206d0,208d0,209d0,208d0,206d0,195d0,200d0,200d0/
DATA ProfYT/0.0d0,0.2d0,0.4d0,0.6d0,0.8d0,1.0d0,1.2d0,1.4d0,1.6d0,1.8d0,2.0d0,&
            203d0,205d0,207d0,208d0,209d0,209d0,209d0,209d0,209d0,209d0,209d0/

DO i=1,Tracer%ndof
 X = mg_mesh%level(ilev)%dcorvg(1,i)
 Y = mg_mesh%level(ilev)%dcorvg(2,i)
 Z = mg_mesh%level(ilev)%dcorvg(3,i)

 dAlpha = (22.5d0)*PI/180d0
 XT = X*cos(dAlpha) - Y*sin(dAlpha)
 YT = X*sin(dAlpha) + Y*cos(dAlpha)
 ZT = Z
 
 IF (Tracer%knpr(i).eq.1) THEN
  
  dist = SQRT((XT-x0)**2d0 + (YT-y0)**2d0)
  if (dist.lt.r0) then
   yR = 1d0 + (xT-x0)/r0
   xR = 2d0-(1d0 + (yT-y0)/r0)
   do l=1,10
    if (xR.gt.ProfXT(l,1).and.xR.le.ProfXT(l+1,1)) THEN
     frac = (xR - ProfXT(l,1))/(ProfXT(l+1,1) - ProfXT(l,1))
     Tx = ProfXT(l,2) + frac*(ProfXT(l+1,2)-ProfXT(l,2))
    end if
   end do
   do l=1,10
    if (yR.gt.ProfYT(l,1).and.yR.le.ProfYT(l+1,1)) THEN
     frac = (yR - ProfYT(l,1))/(ProfYT(l+1,1) - ProfYT(l,1))
     Ty = ProfYT(l,2) + frac*(ProfYT(l+1,2)-ProfYT(l,2)) - 209d0
    end if
   end do
   Tracer%val(NLMAX)%x(i)= Tx + Ty
  else
   Tracer%val(NLMAX)%x(i)= 200d0
  end if
    
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val_Weber
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val()
REAL*8 X,Y,Z
REAL*8 :: PI=4d0*DATAN(1d0)
INTEGER i,l

DO i=1,Tracer%ndof
 X = mg_mesh%level(ilev)%dcorvg(1,i)
 Y = mg_mesh%level(ilev)%dcorvg(2,i)
 Z = mg_mesh%level(ilev)%dcorvg(3,i)

 IF (Tracer%knpr(i).eq.1) THEN
  
   Tracer%val(NLMAX)%x(i)= 1d0
    
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val_EWIKON()
REAL*8 X,Y,Z
REAL*8 :: PI=4d0*DATAN(1d0)
INTEGER i,l

DO i=1,Tracer%ndof
 X = mg_mesh%level(ilev)%dcorvg(1,i)
 Y = mg_mesh%level(ilev)%dcorvg(2,i)
 Z = mg_mesh%level(ilev)%dcorvg(3,i)

 IF (Tracer%knpr(i).eq.1) THEN
  
   Tracer%val(NLMAX)%x(i)= myProcess%AirTemperature
    
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val_EWIKON
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
subroutine AddSource()

return

end subroutine AddSource
!
! ----------------------------------------------
!
subroutine AddSource_EWIKON()
integer i,iSeg,iMat
real*8 dSource

!   ! Q_dot = [kW]
!   ! rho   = [g/cm3]
!   ! Cp    = [kJ/(kg*K)] = [J/(g*K)]
!   ! V     = [cm3]
!   
!   ! q_dot         Q_dot                    kW                       1000 * J/s           1000 K
!   !__________ =  ____________   =   [_______________________] = [ _______________ ] = [__________ ]
!   ! rho * Cp     V * rho * Cp        cm3 * (g/cm3) * J/(g*K)           J/K                 s
! 
DO i=1,Tracer%ndof
 iSeg = myHeatObjects%Segment(i)
 if (iSeg.eq.0) then
  iMat= 0
 else
  iMat = mySigma%mySegment(iSeg)%MatInd
 end if
 
 IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE'.OR.TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'MELT') THEN
  dSource = 1e3*mySigma%mySegment(iSeg)%HeatSource/(mySigma%mySegment(iSeg)%Volume*myMaterials(iMat)%cp*myMaterials(iMat)%Density)
!   if (dSource.ne.0d0) write(*,*) dSource
  Tracer%def(i) = Tracer%def(i) + MLmat(i)*dSource*tstep
 END IF
END DO

! DO i=1,Tracer%ndof
!  iSeg = myHeatObjects%Segment(i)
!  IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
!   Tracer%def(i) = Tracer%def(i) + MLmat(i)*mySigma%mySegment(iSeg)%HeatSource*tstep
!  END IF
! END DO


END subroutine AddSource_EWIKON
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
                            QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW)
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
 !   
 ! if (myid.eq.1) write(*,*) mgDiffCoeff(NLMAX)%x
 ! Diffusion matrix 
 CALL Create_DiffMat(mgDiffCoeff(NLMAX)%x)

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

END SUBROUTINE  updateHeatGeometry
!
! ----------------------------------------------
!
include 'LinSc_transport_extensions.f90'
!
! ----------------------------------------------
!
END MODULE Transport_Q1
