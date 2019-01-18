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
 
 IF (iSeg.eq.0) Tracer%val(NLMAX)%x(i) = myProcess%AirTemperature
 IF (iSeg.ne.0) then
  Tracer%val(NLMAX)%x(i) = mySigma%mySegment(iSeg)%InitTemp
 END IF

END DO

END SUBROUTINE LinSc_InitCond_EWIKON
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
  
   Tracer%val(NLMAX)%x(i)= 45d0
    
 END IF

 IF (Tracer%knpr(i).eq.2) THEN
  
   Tracer%val(NLMAX)%x(i)= 60d0
    
 END IF

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
