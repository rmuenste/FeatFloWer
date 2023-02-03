
!
! ----------------------------------------------
!
SUBROUTINE CheckFBMElem8(P8,F8,I8,X8,V8,iE)
USE particle_step, ONLY : GetPointFromElement
real*8 p8(3,8),x8(8),v8(3,8)
integer f8(8),i8(8),ie

real*8 IntersectionPoint(3),LineEquation(3),daux,relC(3),dxBas
real*8 :: ParticleRadius=0.4d0
integer iC,jC,nC,i,iX,j,iVert
REAL*8 PX,PY,PZ,RX,RY,RZ,VX,VY,VZ,dS,dVal,dist,dd
logical bFound
REAL*8, allocatable :: Centers(:,:),Velocities(:,:),UniqueCenters(:,:),UniqueVelocities(:,:)
logical, allocatable :: Unique(:)
integer nUniqueCenters,iUniqueCenter
real*8 mindist

 nC = 0
 iC = 0
 do i=1,8
  if (F8(i).ne.0) nC = nC + 1
 end do

! if (nC.ne.0.and.nC.ne.8) THEN
 if (nC.ne.0) THEN

  allocate(Velocities(3,nC),centers(3,nC),unique(nC))
  unique = .true.
  
  iC = 0
  DO i=1,8
   PX = P8(1,i)
   PY = P8(2,i)
   PZ = P8(3,i)
   IF (F8(i).ne.0) then
    CALL fbm_Recover_Weigths(PX,PY,PZ,I8(i),RX,RY,RZ,VX,VY,VZ)
    iC = iC + 1
    centers(:,iC) = [RX,RY,RZ]
    Velocities(:,iC) = [VX,VY,VZ]
   end if
  end do
 
  DO iC=1,nC-1
   PX = centers(1,iC)
   PY = centers(2,iC)
   PZ = centers(3,iC)
   DO jC=iC+1,nC
    if (unique(jC)) then
     RX = centers(1,jC)
     RY = centers(2,jC)
     RZ = centers(3,jC)
     
     dist = SQRT((PX-RX)**2d0 + (PY-RY)**2d0 + (PZ-RZ)**2d0)
     
     if (dist.lt.1d-2*ParticleRadius) unique(jC) = .false.
    end if
   END DO
  END DO
  
  nUniqueCenters = 0
  DO iC=1,nC
   if (unique(iC)) nUniqueCenters = nUniqueCenters + 1
  end do
  
  ALLOCATE(UniqueCenters(3,nUniqueCenters))
  ALLOCATE(UniqueVelocities(3,nUniqueCenters))
   
  nUniqueCenters = 0
  DO iC=1,nC
   if (unique(iC)) THEN
    nUniqueCenters = nUniqueCenters +1
    UniqueCenters(:,nUniqueCenters)    = centers(:,iC)
    UniqueVelocities(:,nUniqueCenters) = Velocities(:,iC)
   end if
  end do
  
  if (nUniqueCenters.gt.1) then
!    write(*,*) 'nUniqueCenters: ',nUniqueCenters
  end if
  
  DO i=1,8
   PX = P8(1,i)
   PY = P8(2,i)
   PZ = P8(3,i)
   
   mindist = 1d8
   DO iC = 1,nUniqueCenters
    dist = (PX-UniqueCenters(1,iC))**2d0 + (PY-UniqueCenters(2,iC))**2d0 + (PZ-UniqueCenters(3,iC))**2d0
    if (dist.lt.mindist) then
     mindist = dist
     iUniqueCenter = iC
     RX = UniqueCenters(1,iC)
     RY = UniqueCenters(2,iC)
     RZ = UniqueCenters(3,iC)
     VX = UniqueVelocities(1,iC)
     VY = UniqueVelocities(2,iC)
     VZ = UniqueVelocities(3,iC)
    end if
   END DO
   
   LineEquation = [PX-RX,PY-RY,PZ-RZ]
   
   daux = sqrt(LineEquation(1)**2d0 + LineEquation(2)**2d0 + LineEquation(3)**2d0)

   if (daux.gt.ParticleRadius) then

!    dVal = 1d0-min(1d0,max(0d0,(daux - ParticleRadius)/0.15d0))
      
    daux = ParticleRadius/daux
  
    IntersectionPoint = [RX + daux*LineEquation(1),&
                         RY + daux*LineEquation(2),&
                         RZ + daux*LineEquation(3)]
    
    dist = 1d8
    do j=1,8                     
     dd= sqrt((IntersectionPoint(1)-P8(1,j))**2d0 + &
              (IntersectionPoint(2)-P8(2,j))**2d0 + &
              (IntersectionPoint(3)-P8(3,j))**2d0) 
     if (dd.lt.dist) then
      dist = dd
      iVert = j
     end if
    end do

    iX = 0
    bFound = .false.
    CALL GetPointFromElement(P8(:,1:8),IntersectionPoint,iVert,bFound,relC,iX)
    
    if (bFound) then
     CALL ReturnE011Bas(i,dVal,relC(1),relC(2),relC(3))
    else
      dVal = 0d0
    end if
   else
    dVal = 1d0
   end if
   
   IF (dVal.gt.x8(i)) then
!     x8(i) = dVal
    if (dVal.eq.1d0) THEN 
     x8(i) = 1d0 !dVal
    else
     x8(i) = 0d0 !dVal
    end if
    v8(:,i) = [VX,VY,VZ]
   end if
   
  end do
  
  DEALLOCATE(UniqueCenters,Centers,Velocities,UniqueVelocities,Unique)
 end if
 
 if (iC.eq.0) then
  DO i=1,8
   dval = 0d0
   x8(i) = max(x8(i),dVal)
  end do
 end if

 DO i=1,8
  if (F8(i).ne.0) THEN
   x8(i) = 1d0
  end if
 end do
  
 contains 
 
 subroutine ReturnE011Bas(iBas,dBas,x1,x2,x3)
 real*8 x1,x2,x3
 integer ibas
 real*8 dbas
 real*8 :: Q8=0.125d0
 
 if (iBas.eq.1) dbas = Q8*(1D0-X1)*(1D0-X2)*(1D0-X3)
 if (iBas.eq.2) dbas = Q8*(1D0+X1)*(1D0-X2)*(1D0-X3)
 if (iBas.eq.3) dbas = Q8*(1D0+X1)*(1D0+X2)*(1D0-X3)
 if (iBas.eq.4) dbas = Q8*(1D0-X1)*(1D0+X2)*(1D0-X3)
 if (iBas.eq.5) dbas = Q8*(1D0-X1)*(1D0-X2)*(1D0+X3)
 if (iBas.eq.6) dbas = Q8*(1D0+X1)*(1D0-X2)*(1D0+X3)
 if (iBas.eq.7) dbas = Q8*(1D0+X1)*(1D0+X2)*(1D0+X3)
 if (iBas.eq.8) dbas = Q8*(1D0-X1)*(1D0+X2)*(1D0+X3)
 
 end subroutine ReturnE011Bas
                       
END SUBROUTINE CheckFBMElem8
!
! ----------------------------------------------
!
SUBROUTINE CheckFBMElem27(P27,F27,I27,X27,V27,iE)
USE particle_step, ONLY : GetPointFromElement
real*8 p27(3,27),x27(27),v27(3,27)
integer f27(27),i27(27),ie

real*8 IntersectionPoint(3),LineEquation(3),daux,relC(3),dxBas
real*8 :: ParticleRadius=0.4d0
integer iC,jC,nC,i,iX,j,iVert
REAL*8 PX,PY,PZ,RX,RY,RZ,VX,VY,VZ,dS,dVal,dist,dd
logical bFound
REAL*8, allocatable :: Centers(:,:),Velocities(:,:),UniqueCenters(:,:),UniqueVelocities(:,:)
logical, allocatable :: Unique(:)
integer nUniqueCenters,iUniqueCenter
real*8 mindist

 nC = 0
 do i=1,27
  if (F27(i).ne.0) nC = nC + 1
 end do

! if (nC.ne.0.and.nC.ne.27) THEN
 if (nC.ne.0) THEN

  allocate(Velocities(3,nC),centers(3,nC),unique(nC))
  unique = .true.
  
  iC = 0
  DO i=1,27
   PX = P27(1,i)
   PY = P27(2,i)
   PZ = P27(3,i)
   IF (F27(i).ne.0) then
    CALL fbm_Recover_Weigths(PX,PY,PZ,I27(i),RX,RY,RZ,VX,VY,VZ)
    iC = iC + 1
    centers(:,iC) = [RX,RY,RZ]
    Velocities(:,iC) = [VX,VY,VZ]
   end if
  end do
 
  DO iC=1,nC-1
   PX = centers(1,iC)
   PY = centers(2,iC)
   PZ = centers(3,iC)
   DO jC=iC+1,nC
    if (unique(jC)) then
     RX = centers(1,jC)
     RY = centers(2,jC)
     RZ = centers(3,jC)
     
     dist = SQRT((PX-RX)**2d0 + (PY-RY)**2d0 + (PZ-RZ)**2d0)
     
     if (dist.lt.1d-2*ParticleRadius) unique(jC) = .false.
    end if
   END DO
  END DO
  
  nUniqueCenters = 0
  DO iC=1,nC
   if (unique(iC)) nUniqueCenters = nUniqueCenters + 1
  end do
  
  ALLOCATE(UniqueCenters(3,nUniqueCenters))
  ALLOCATE(UniqueVelocities(3,nUniqueCenters))
   
  nUniqueCenters = 0
  DO iC=1,nC
   if (unique(iC)) THEN
    nUniqueCenters = nUniqueCenters +1
    UniqueCenters(:,nUniqueCenters)    = centers(:,iC)
    UniqueVelocities(:,nUniqueCenters) = Velocities(:,iC)
   end if
  end do
  
  if (nUniqueCenters.gt.1) then
!    write(*,*) 'nUniqueCenters: ',nUniqueCenters
  end if
  
  DO i=1,27
   PX = P27(1,i)
   PY = P27(2,i)
   PZ = P27(3,i)
   
   mindist = 1d8
   DO iC = 1,nUniqueCenters
    dist = (PX-UniqueCenters(1,iC))**2d0 + (PY-UniqueCenters(2,iC))**2d0 + (PZ-UniqueCenters(3,iC))**2d0
    if (dist.lt.mindist) then
     mindist = dist
     iUniqueCenter = iC
     RX = UniqueCenters(1,iC)
     RY = UniqueCenters(2,iC)
     RZ = UniqueCenters(3,iC)
     VX = UniqueVelocities(1,iC)
     VY = UniqueVelocities(2,iC)
     VZ = UniqueVelocities(3,iC)
    end if
   END DO
   
   LineEquation = [PX-RX,PY-RY,PZ-RZ]
   
   daux = sqrt(LineEquation(1)**2d0 + LineEquation(2)**2d0 + LineEquation(3)**2d0)

   if (daux.gt.ParticleRadius) then

   dVal = 1d0-min(1d0,max(0d0,(daux - ParticleRadius)/0.15d0))
      
!     daux = ParticleRadius/daux
!   
!     IntersectionPoint = [RX + daux*LineEquation(1),&
!                          RY + daux*LineEquation(2),&
!                          RZ + daux*LineEquation(3)]
!     
!     dist = 1d8
!     do j=1,8                     
!      dd= sqrt((IntersectionPoint(1)-P27(1,j))**2d0 + &
!               (IntersectionPoint(2)-P27(2,j))**2d0 + &
!               (IntersectionPoint(3)-P27(3,j))**2d0) 
!      if (dd.lt.dist) then
!       dist = dd
!       iVert = j
!      end if
!     end do
! 
!     iX = 0
!     bFound = .false.
!     CALL GetPointFromElement(P27(:,1:8),IntersectionPoint,iVert,bFound,relC,iX)
!     
!     if (bFound) then
!      CALL ReturnE013Bas(i,dVal,relC(1),relC(2),relC(3))
!     else
!       dVal = 0d0
!     end if
   else
    dVal = 1d0
   end if
   
   IF (dVal.gt.x27(i)) then
    x27(i) = dVal
!     if (dVal.eq.1d0) THEN 
!      x27(i) = 1d0 !dVal
!     else
!      x27(i) = 0d0 !dVal
!     end if
    v27(:,i) = [VX,VY,VZ]
   end if
   
  end do
  
  DEALLOCATE(UniqueCenters,Centers,Velocities,UniqueVelocities,Unique)
 end if
 
 if (iC.eq.0) then
  DO i=1,27
   dval = 0d0
   x27(i) = max(x27(i),dVal)
  end do
 end if

 DO i=1,27
  if (F27(i).ne.0) THEN
   x27(i) = 1d0
  end if
 end do
  
 contains 
 
 subroutine ReturnE013Bas(iBas,dBas,x1,x2,x3)
 real*8 x1,x2,x3
 integer ibas
 real*8 dbas
 real*8 :: Q8=0.125d0
 
 if (iBas.eq.1) dbas = -Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0-X3)
 if (iBas.eq.2) dbas =  Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0-X3)
 if (iBas.eq.3) dbas = -Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0-X3)
 if (iBas.eq.4) dbas = Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0-X3)
 if (iBas.eq.5) dbas = Q8*X1*(1D0-X1)*X2*(1D0-X2)*X3*(1D0+X3)
 if (iBas.eq.6) dbas = -Q8*X1*(1D0+X1)*X2*(1D0-X2)*X3*(1D0+X3)
 if (iBas.eq.7) dbas = Q8*X1*(1D0+X1)*X2*(1D0+X2)*X3*(1D0+X3)
 if (iBas.eq.8) dbas = -Q8*X1*(1D0-X1)*X2*(1D0+X2)*X3*(1D0+X3)
 if (iBas.eq.9) dbas = Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0-X3)
 if (iBas.eq.10) dbas = -Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
 if (iBas.eq.11) dbas = -Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0-X3)
 if (iBas.eq.12) dbas = Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
 if (iBas.eq.13) dbas = Q8*X1*(1D0-X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.14) dbas = -Q8*X1*(1D0+X1)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.15) dbas = Q8*X1*(1D0+X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.16) dbas = -Q8*X1*(1D0-X1)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.17) dbas = -Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*X3*(1D0+X3)
 if (iBas.eq.18) dbas = Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
 if (iBas.eq.19) dbas = Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*X3*(1D0+X3)
 if (iBas.eq.20) dbas = -Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
 if (iBas.eq.21) dbas = -Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0-X3)
 if (iBas.eq.22) dbas = -Q8*(2D0-2D0*X1**2D0)*X2*(1D0-X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.23) dbas = Q8*X1*(1D0+X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
 if (iBas.eq.24) dbas = Q8*(2D0-2D0*X1**2D0)*X2*(1D0+X2)*(2D0-2D0*X3**2D0)
 if (iBas.eq.25) dbas = -Q8*X1*(1D0-X1)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
 if (iBas.eq.26) dbas = Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*X3*(1D0+X3)
 if (iBas.eq.27) dbas = Q8*(2D0-2D0*X1**2D0)*(2D0-2D0*X2**2D0)*(2D0-2D0*X3**2D0)
 
 end subroutine ReturnE013Bas
                       
END SUBROUTINE CheckFBMElem27
!
! ----------------------------------------------
!
SUBROUTINE Set_RHSBC_QuadScalar()
  INTEGER i
  REAL*8 daux,Alpha_U,Alpha_V,Alpha_W

  ilev = NLMAX
  call setlev(2)
  
  DO i=1,QuadSc%ndof
  
    Alpha_U = FBMWeight(i)
    Alpha_V = FBMWeight(i)
    Alpha_W = FBMWeight(i)
    
    IF (QuadSc%knprU(nlmax)%x(i).NE.0) Alpha_U = 1d0
    IF (QuadSc%knprV(nlmax)%x(i).NE.0) Alpha_V = 1d0
    IF (QuadSc%knprW(nlmax)%x(i).NE.0) Alpha_W = 1d0

    
    QuadSc%defU(i) = (1d0-Alpha_U)*QuadSc%defU(i) + Alpha_U*FBMVelocity(1,i)
    QuadSc%defV(i) = (1d0-Alpha_V)*QuadSc%defV(i) + Alpha_V*FBMVelocity(2,i)
    QuadSc%defW(i) = (1d0-Alpha_W)*QuadSc%defW(i) + Alpha_W*FBMVelocity(3,i)

  END DO

END SUBROUTINE Set_RHSBC_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Set_MATBC_QuadScalar_sub(DA11,DA22,DA33,KLD,&
    KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8  DA11(*),DA22(*),DA33(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF
  REAL*8 daux,Alpha_U,Alpha_V,Alpha_W

  DO I=1,NDOF
  
   Alpha_U = FBMWeight(i)
   Alpha_V = FBMWeight(i)
   Alpha_W = FBMWeight(i)
    
   IF (QuadSc%knprU(nlmax)%x(i).NE.0) Alpha_U = 1d0
   IF (QuadSc%knprV(nlmax)%x(i).NE.0) Alpha_V = 1d0
   IF (QuadSc%knprW(nlmax)%x(i).NE.0) Alpha_W = 1d0
    
   ICOL = KLD(I)
   DA11(ICOL) = (1d0-Alpha_U)*DA11(ICOL) + Alpha_U*1d0
   DA22(ICOL) = (1d0-Alpha_V)*DA22(ICOL) + Alpha_V*1d0
   DA33(ICOL) = (1d0-Alpha_W)*DA33(ICOL) + Alpha_W*1d0
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = (1d0-Alpha_U)*DA11(ICOL)
    DA22(ICOL) = (1d0-Alpha_V)*DA22(ICOL)
    DA33(ICOL) = (1d0-Alpha_W)*DA33(ICOL)
   END DO

  END DO
  
END SUBROUTINE Set_MATBC_QuadScalar_sub
!
! ----------------------------------------------
!
SUBROUTINE Set_MATBC_QuadScalar_9_sub(DA11,DA22,DA33,DA12,DA13,DA23,DA21,DA31,DA32,&
    KLD,KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8 DA11(*),DA22(*),DA33(*),DA12(*),DA13(*),DA23(*),DA21(*),DA31(*),DA32(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF
  REAL*8 daux,Alpha_U,Alpha_V,Alpha_W

  DO I=1,NDOF
   Alpha_U = FBMWeight(i)
   Alpha_V = FBMWeight(i)
   Alpha_W = FBMWeight(i)
  
   IF (QuadSc%knprU(nlmax)%x(i).NE.0) Alpha_U = 1d0
   IF (QuadSc%knprV(nlmax)%x(i).NE.0) Alpha_V = 1d0
   IF (QuadSc%knprW(nlmax)%x(i).NE.0) Alpha_W = 1d0
   
   ICOL = KLD(I)
   DA11(ICOL) = (1d0-Alpha_U)*DA11(ICOL) + Alpha_U*1d0
   DA22(ICOL) = (1d0-Alpha_V)*DA22(ICOL) + Alpha_V*1d0
   DA33(ICOL) = (1d0-Alpha_W)*DA33(ICOL) + Alpha_W*1d0
   
   DA12(ICOL) = (1d0-Alpha_U)*DA12(ICOL)
   DA13(ICOL) = (1d0-Alpha_U)*DA13(ICOL)
   DA21(ICOL) = (1d0-Alpha_V)*DA21(ICOL)
   DA23(ICOL) = (1d0-Alpha_V)*DA23(ICOL)
   DA31(ICOL) = (1d0-Alpha_W)*DA31(ICOL)
   DA32(ICOL) = (1d0-Alpha_W)*DA32(ICOL)
   
   DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = (1d0-Alpha_U)*DA11(ICOL)
    DA22(ICOL) = (1d0-Alpha_V)*DA22(ICOL)
    DA33(ICOL) = (1d0-Alpha_W)*DA33(ICOL)
    
    DA12(ICOL) = (1d0-Alpha_U)*DA12(ICOL)
    DA13(ICOL) = (1d0-Alpha_U)*DA13(ICOL)
    DA21(ICOL) = (1d0-Alpha_V)*DA21(ICOL)
    DA23(ICOL) = (1d0-Alpha_V)*DA23(ICOL)
    DA31(ICOL) = (1d0-Alpha_W)*DA31(ICOL)
    DA32(ICOL) = (1d0-Alpha_W)*DA32(ICOL)
   END DO
   
  END DO

END SUBROUTINE Set_MATBC_QuadScalar_9_sub
