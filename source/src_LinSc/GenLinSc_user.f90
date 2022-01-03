!
! ----------------------------------------------
!
SUBROUTINE Boundary_GenLinSc_Q1_Val(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
REAL*8 dcorvg(3,*)
REAL*8 X,Y,Z,dFrac,dTemp
REAL*8 :: RX = 0.0d0,RY = 0.0d0
REAL*8 :: RADx = 0.2d0
REAL*8 DIST
INTEGER i,ifld

do iFld=1,GenLinScalar%nOfFields
 DO i=1,GenLinScalar%ndof
 
  X = dcorvg(1,i)
  Y = dcorvg(2,i)
  Z = dcorvg(3,i)
  
  DIST = sqrt((X-5.5d0)**2d0 + (Y+6.5d0)**2d0 + (Z+0.5)**2d0)
  dFrac = DIST/1.5d0
  dTemp = 240d0 - 20d0*dFrac
  if (Dist.gt.1.75) dTemp = 260d0
!   dFrac = (5d0+X)/10d0
!   dTemp = 220d0 + 30d0*dFrac
  
  
  
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'k') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = 0d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 1d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.3) GenLinScalar%Fld(iFld)%val(i) = 12d0
  END IF
 
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'epsilon') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = 4d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 5d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.3) GenLinScalar%Fld(iFld)%val(i) = 6d0
  END IF

  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'nut') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = 7d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 8d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.3) GenLinScalar%Fld(iFld)%val(i) = 9d0
  END IF

  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'FAC') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = -1d0
  END IF
  
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'FAC') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 1d0
  END IF

  END DO
END DO

END SUBROUTINE Boundary_GenLinSc_Q1_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_GenLinSc_Q1_Mat(DA,KLD,&
           KNPR,NDOF)
REAL*8  DA(*)
INTEGER KLD(*),KNPR(*),ICOL,I,NDOF

DO I=1,NDOF
 IF (KNPR(I).ne.0) THEN
   ICOL = KLD(I)
   DA(ICOL) = 1d0
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA(ICOL) = 0d0
   END DO
 END IF
END DO

END SUBROUTINE Boundary_GenLinSc_Q1_Mat
!
! ----------------------------------------------
!
SUBROUTINE Knpr_GenLinSc_HEATALPHA_Q1(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,myMultiMat
IMPLICIT NONE
REAL*8 dcorvg(3,*)
real*8 X,Y,Z,DIST
INTEGER i,j,ifld
integer iInflow,iSubInflow,iMat,mySubinflow
REAL*8 dInnerRadius,dOuterRadius,dMassFlow,dVolFlow,daux,dInnerInflowRadius,dDensity
REAL*8 dCenter(3),dNormal(3),dProfil(3),dScale
real*8, dimension(11) :: x_arr, y_arr, CC, DD, MM
REAL*8 :: PI=dATAN(1d0)*4d0,myTwoPI=2d0*dATAN(1d0)*4d0
REAL*8 :: U_bar, h, normalizedTime, val,dFact
logical bBC,iBC

DO i=1,GenLinScalar%ndof

  bBC = .false.
  
  X = dcorvg(1,i)
  Y = dcorvg(2,i)
  Z = dcorvg(3,i)
  
  iInflow = ABS(myBoundary%iInflow(i))
  
  if (iInflow.gt.0) then
   IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
    mySubinflow = 1
    IF (myProcess%myInflow(iInflow)%iBCType.eq.1) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     ddensity      = myMultiMat%Mat(iMat)%Thermodyn%density
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = RotParabolicVelo3D(dMassFlow,dDensity,dOuterRadius)
     daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
     if (daux.ne.0d0) then
      bBC=.true.
     END IF
    END IF

    IF (myProcess%myInflow(iInflow)%iBCType.eq.2) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     ddensity      = myMultiMat%Mat(iMat)%Thermodyn%density
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = RotDoubleParabolicVelo3D(dMassFlow,dDensity,dInnerRadius,dOuterRadius)
     daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
     if (daux.ne.0d0) then
      bBC=.true.
     END IF
    END IF
    
    IF (myProcess%myInflow(iInflow)%iBCType.eq.3) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     ddensity      = myMultiMat%Mat(iMat)%Thermodyn%density
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dProfil = FlatVelo3D(dMassFlow,dDensity,dOuterRadius)
     daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
     if (daux.ne.0d0) then
      bBC=.true.
     END IF
    END IF
    
    IF (myProcess%myInflow(iInflow)%iBCType.eq.4) then
     dCenter       = myProcess%myInflow(iInflow)%center
     dNormal       = myProcess%myInflow(iInflow)%normal
     dMassFlow     = myProcess%myInflow(iInflow)%massflowrate
     iMat          = myProcess%myInflow(iInflow)%Material
     ddensity      = myMultiMat%Mat(iMat)%Thermodyn%density
     douterradius  = myProcess%myInflow(iInflow)%outerradius
     dinnerradius  = myProcess%myInflow(iInflow)%innerradius
     dProfil = CurvedFlatVelo3D(dMassFlow,dDensity,dInnerRadius,dOuterRadius)
     daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
     if (daux.ne.0d0) then
      bBC=.true.
     END IF
    END IF
    
   ELSE
   
    DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
     IF (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCType.eq.1) then
      dCenter       = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%center
      dNormal       = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%normal
      dMassFlow     = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate
      iMat          = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material
      ddensity      = myMultiMat%Mat(iMat)%Thermodyn%density
      douterradius  = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius
      dinnerradius  = myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius
      dProfil = RotParabolicVelo3D(dMassFlow,dDensity,dOuterRadius)
      daux = dProfil(1)**2d0 + dProfil(2)**2d0 + dProfil(3)**2d0
      if (daux.ne.0d0) then
       mySubinflow = iSubInflow
       bBC=.true.
      END IF
     END IF
    END DO
    
   END IF
   
   do iFld=1,GenLinScalar%nOfFields
    IF (bBC) THEN
!      write(*,*) 'mySubInflow',mySubInflow,'iInflow',iInflow
     GenLinScalar%fld(iFld)%knpr(i) = mySubInflow
    ELSE
     GenLinScalar%fld(iFld)%knpr(i) = 0
    END IF
   END DO
 END IF
END DO !i

return

 CONTAINS
include '../include/ProfileFunctions.f90'

END SUBROUTINE Knpr_GenLinSc_HEATALPHA_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_GenLinSc_HEATALPHA_Q1_Val(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
USE Sigma_User, ONLY: mySigma,myThermodyn,myProcess,myMultiMat
REAL*8 dcorvg(3,*)
REAL*8 X,Y,Z,dFrac,dTemp
REAL*8 DIST,TempBC
INTEGER i,ifld,mySubinflow,iInflow,iMat

if (myid.eq.0) return

DO i=1,GenLinScalar%ndof
 
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)
 
 do iFld=1,GenLinScalar%nOfFields
 
  iInflow=abs(myBoundary%iInflow(i))
  mySubinflow = GenLinScalar%Fld(iFld)%knpr(i)
  iMat = 0
  TempBC = 0d0
  
  IF (mySubinflow.ne.0) then
!    write(*,*) iInflow,mySubinflow
   IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
    iMat   = myProcess%myInflow(iInflow)%Material
    TempBC = myProcess%myInflow(iInflow)%Temperature
   ELSE
    iMat   = myProcess%myInflow(iInflow)%mySubInflow(mySubinflow)%Material
    TempBC = myProcess%myInflow(iInflow)%mySubInflow(mySubinflow)%Temperature
   END IF
   
   IF (iFld.eq.1) then
    GenLinScalar%Fld(iFld)%val(i) = TempBC
   ELSE
    IF (iMat.eq.iFld-1) THEN
     GenLinScalar%Fld(iFld)%val(i) = 1d0
    else
     GenLinScalar%Fld(iFld)%val(i) = 0d0
    end if
   END IF
!    WRITE(*,*) TRIM(GenLinScalar%Fld(iFld)%cName), iInflow,mySubinflow, iMAt, GenLinScalar%Fld(iFld)%val(i), TempBC
  END IF
 END DO
END DO

! pause

END SUBROUTINE Boundary_GenLinSc_HEATALPHA_Q1_Val
!
! ----------------------------------------------
!
SUBROUTINE Knpr_GenLinSc_Q1(dcorvg)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
IMPLICIT NONE
REAL*8 dcorvg(3,*),X,Y,Z,DIST,xx
REAL*8 :: PX=0.5d0,PY=0.2d0,PZ=0.2d0,RAD=0.050d0
INTEGER i,ifld

do iFld=1,GenLinScalar%nOfFields
 DO i=1,GenLinScalar%ndof
  X = dcorvg(1,i)
  Y = dcorvg(2,i)
  Z = dcorvg(3,i)
  GenLinScalar%fld(iFld)%knpr(i) = 0

  IF (myBoundary%iInflow(i).eq.171) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 1
  END IF
  
  IF (myBoundary%iInflow(i).eq.172) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 2
  END IF
  
  IF (myBoundary%iInflow(i).eq.2.and.(TRIM(GenLinScalar%Fld(iFld)%cName).eq.'FAC')) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 1
  END IF
  
  IF (myBoundary%bWall(i).and.(TRIM(GenLinScalar%Fld(iFld)%cName).eq.'FAC')) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 2
  END IF

 END DO
END DO


END SUBROUTINE Knpr_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Boundary_GenLinSc_Q1_Def()
INTEGER i,iFld

DO iFld=1,GenLinScalar%nOfFields
 DO i=1,GenLinScalar%ndof
  IF (GenLinScalar%Fld(iFld)%knpr(i).ne.0) THEN
   GenLinScalar%Fld(iFld)%def(i) = 0d0
  END IF
 END DO
END DO

END SUBROUTINE Boundary_GenLinSc_Q1_Def
