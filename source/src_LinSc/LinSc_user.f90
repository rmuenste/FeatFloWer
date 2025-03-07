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
SUBROUTINE LinSc_InitCond_XSE(dcorvg)
!use Sigma_User, only : myProcess
REAL*8, dimension(:,:), pointer :: dcorvg
REAL*8 X,Y,Z
REAL*8 :: RX = 0.0d0, RY = 0.0d0, RZ = 2.4d0
REAL*8 :: RADx = 0.20d0,RADs=0.040
REAL*8 DIST
INTEGER i,iSeg,iMat

RETURN

END SUBROUTINE LinSc_InitCond_XSE
!
! ----------------------------------------------
!
SUBROUTINE LinSc_InitCond_General(dcorvg)
!use Sigma_User, only : myProcess
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

 Tracer%val(NLMAX)%x(i) = myProcess%T0
  
END DO

END SUBROUTINE LinSc_InitCond_General
!
! ----------------------------------------------
!
SUBROUTINE LinSc_InitCond_Weber(dcorvg)
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

 Tracer%val(NLMAX)%x(i) = 205d0
  
END DO

END SUBROUTINE LinSc_InitCond_Weber
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
 
 IF (iSeg.eq.0) Tracer%val(NLMAX)%x(i) = myProcess%AmbientTemperature
 IF (iSeg.ne.0) then
  Tracer%val(NLMAX)%x(i) = mySigma%mySegment(iSeg)%InitTemp
 END IF

END DO

END SUBROUTINE LinSc_InitCond_EWIKON
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val_XSE()

REAL*8 X,Y,Z
INTEGER i,iInflow
REAL*8 :: T1=35d0,T2=45d0
real*8 dc(3),dR,dRR,dT,TempBC

DO i=1,Tracer%ndof

 X = mg_mesh%level(nlmax)%dcorvg(1,i)
 Y = mg_mesh%level(nlmax)%dcorvg(2,i)
 Z = mg_mesh%level(nlmax)%dcorvg(3,i)
 
 IF (Tracer%knpr(i).ge.100.and.Tracer%knpr(i).lt.1000) THEN
  Tracer%val(NLMAX)%x(i)= myProcess%SegThermoPhysProp(Tracer%knpr(i)-100)%T_Const
 END IF

 IF (Tracer%knpr(i).eq.1) THEN
  Tracer%val(NLMAX)%x(i)= myProcess%T0
    
!     iInflow = 1
!     
!     dC = myProcess%myInflow(iInflow)%center
!     
!     IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.0) THEN
!      TempBC = myProcess%myInflow(iInflow)%Temperature
!     END IF
!     
!     IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.1) THEN
!      dRR = myProcess%myInflow(iInflow)%outerradius
!      dR = SQRT((dC(1)-X)**2d0 + (dC(2)-Y)**2d0 + (dC(2)-Z)**2d0)
!      dR = Min(dR,dRR)
!      dT = myProcess%myInflow(iInflow)%TemperatureRange
!      
!      TempBC = myProcess%myInflow(iInflow)%Temperature + dT*(dRR-dR)/dRR
!    END IF
!    IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.2) THEN
!      dRR = myProcess%myInflow(iInflow)%outerradius
!      dR = SQRT((dC(1)-X)**2d0 + (dC(2)-Y)**2d0 + (dC(2)-Z)**2d0)
!      dR = Min(dR,dRR)
!      dT = myProcess%myInflow(iInflow)%TemperatureRange
!      
!      TempBC = myProcess%myInflow(iInflow)%Temperature + dT*(dRR-dR)*(dRR+dR)/(dRR*dRR)
!     END IF
!   
!   Tracer%val(NLMAX)%x(i)= TempBC
 END IF
 
 IF (Tracer%knpr(i).eq.2) THEN
!   if (z.lt.11.25d0) THEN
!    Tracer%val(NLMAX)%x(i)= T1
!   else
!    Tracer%val(NLMAX)%x(i)= T2
!   END IF
  Tracer%val(NLMAX)%x(i)= myProcess%Ta
 END IF
 
 IF (Tracer%knpr(i).eq.3) THEN
  Tracer%val(NLMAX)%x(i)= myProcess%Ti
 END IF

 IF (Tracer%knpr(i).ge.1000.and.Tracer%knpr(i).lt.2000) THEN
  iInflow = Tracer%knpr(i) - 1000

  dC = myProcess%myInflow(iInflow)%center
  
  IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.0) THEN
   TempBC = myProcess%myInflow(iInflow)%Temperature
  END IF
  
  IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.1) THEN
   dRR = myProcess%myInflow(iInflow)%outerradius
   dR = SQRT((dC(1)-X)**2d0 + (dC(2)-Y)**2d0 + (dC(3)-Z)**2d0)
   dT = myProcess%myInflow(iInflow)%TemperatureRange
   dR = Min(dR,dRR)
   
   TempBC = myProcess%myInflow(iInflow)%Temperature + dT*(1d0-(dRR-dR)/dRR)
!    write(*,'(A,18ES12.4)') 'L',X,Y,Z,dR,dRR,dT,dC,TempBC
 END IF
 IF (myProcess%myInflow(iInflow)%Temperaturetype.eq.2) THEN
   dRR = myProcess%myInflow(iInflow)%outerradius
   dR = SQRT((dC(1)-X)**2d0 + (dC(2)-Y)**2d0 + (dC(3)-Z)**2d0)
   dT = myProcess%myInflow(iInflow)%TemperatureRange
   dR = Min(dR,dRR)
   
   TempBC = myProcess%myInflow(iInflow)%Temperature + dT*(1d0-(dRR-dR)*(dRR+dR)/(dRR*dRR))
!   write(*,'(A,8ES12.4)') 'Q',dRR,dR,dT,dC,TempBC
  END IF
  Tracer%val(NLMAX)%x(i)= TempBC
  
 END IF
 
 IF (Tracer%knpr(i).ge.2000) THEN
  TempBC = DBLE(Tracer%knpr(i) - 2000)
  Tracer%val(NLMAX)%x(i)= TempBC
 END IF

END DO

END SUBROUTINE Boundary_LinSc_Val_XSE
!
! ----------------------------------------------
!
SUBROUTINE Boundary_LinSc_Val_General()
use, intrinsic :: ieee_arithmetic

REAL*8 X,Y,Z
INTEGER i
real*8 :: myInf

if(ieee_support_inf(myInf))then
  myInf = ieee_value(myInf, ieee_negative_inf)
endif

DO i=1,Tracer%ndof

 X = mg_mesh%level(ilev)%dcorvg(1,i)
 Y = mg_mesh%level(ilev)%dcorvg(2,i)
 Z = mg_mesh%level(ilev)%dcorvg(3,i)
 
 IF (Tracer%knpr(i).eq.1) THEN
  Tracer%val(NLMAX)%x(i)= myProcess%T0
 END IF
 
 IF (Tracer%knpr(i).eq.2) THEN
  Tracer%val(NLMAX)%x(i)= myProcess%Ta
 END IF
 
END DO

END SUBROUTINE Boundary_LinSc_Val_General
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
INTEGER i,l,iSeg

DO i=1,Tracer%ndof
 X = mg_mesh%level(ilev)%dcorvg(1,i)
 Y = mg_mesh%level(ilev)%dcorvg(2,i)
 Z = mg_mesh%level(ilev)%dcorvg(3,i)

 IF (Tracer%knpr(i).eq.1) THEN
  
   Tracer%val(NLMAX)%x(i)= 45d0
    
 END IF

 IF (Tracer%knpr(i).eq.2) THEN
  
   Tracer%val(NLMAX)%x(i)= 80d0
    
 END IF

 IF (Tracer%knpr(i).eq.3) THEN

   Tracer%val(NLMAX)%x(i)= myProcess%MeltInflowTemperature

 END IF

 IF (Tracer%knpr(i).eq.4) THEN
  
   iSeg = myHeatObjects%Segment(i)
   Tracer%val(NLMAX)%x(i)= mySigma%mySegment(iSeg)%InitTemp
    
 END IF
 
 IF (Tracer%knpr(i).eq.5) THEN

   Tracer%val(NLMAX)%x(i)= 310d0

 END IF

!  IF (Tracer%knpr(i).eq.6) THEN
!   
!    Tracer%val(NLMAX)%x(i)= 25d0
!     
!  END IF
! 
!  IF (Tracer%knpr(i).eq.7) THEN
!   
!    Tracer%val(NLMAX)%x(i)= 200d0
!     
!  END IF

END DO

END SUBROUTINE Boundary_LinSc_Val_EWIKON
!
! ----------------------------------------------
!
subroutine AddSource()

return

end subroutine AddSource
!
! ----------------------------------------------
!
subroutine AddSource_XSE()
USE var_QuadScalar, ONLY : Screw,Shell,Viscosity,Shearrate,dIntegralHeat,mySegmentIndicator

integer i,iSeg
real*8 daux

!    g      1   1         g     
! --------*---*--- = ----------  
!  cm . s   s   s     cm .s.s.s 

dIntegralHeat = 0d0

DO i = 1,Tracer%ndof

 
  daux = Viscosity(i)*(Shearrate(i)**2d0)
  
!  IF (Screw(i).ge.0d0) THEN
!  IF (Screw(i).ge.0d0.and.Shell(i).ge.0d0) THEN
    Tracer%def(i) = Tracer%def(i) + daux * MlMat(i)*tstep
    dIntegralHeat = dIntegralHeat + daux * MlMat(i)
!  END IF
  
END DO

end subroutine AddSource_XSE
!
! ----------------------------------------------
!
subroutine AddSource_General()

return

end subroutine AddSource_General
!
! ----------------------------------------------
!
subroutine AddSource_EWIKON()
integer i,iSeg,iMat
real*8 dSource

!BASIC units L:[cm], M[g], T[s], t[K]
!Q_dot = [kW]
!V     = [cm3]
!
!Q_dot    =  1000 [W] = 1000 * [kg*m2/s3] = 1e3* [1e3g * 1e4cm2 / s3] = 1e10 g*cm2/s3
! 
DO i=1,Tracer%ndof
 iSeg = myHeatObjects%Segment(i)
 if (iSeg.eq.0) then
  iMat= 0
 else
  iMat = mySigma%mySegment(iSeg)%MatInd
 end if
 
 IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE'.OR.TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'MELT') THEN
  dSource = 1e10*mySigma%mySegment(iSeg)%UseHeatSource/(mySigma%mySegment(iSeg)%Volume)
!  dSource = 1e3*mySigma%mySegment(iSeg)%HeatSource/(mySigma%mySegment(iSeg)%Volume*myMaterials(iMat)%cp*myMaterials(iMat)%Density)
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
