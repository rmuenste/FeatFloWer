MODULE module_levelset

USE def_SCALAR
USE PP3D_MPI, ONLY:CommSum,CommSumHalf,Comm_NLComplete,&
              Comm_Maximum,Comm_Summ,knprmpi,myid,master


IMPLICIT NONE

TYPE(mg_scalar) :: levelset
TYPE(mg_vector) :: Density,SurfT_mat,normal
REAL*8 :: Density_Secondary=1d0,Density_Primary=1d0
REAL*8 :: SurfaceTension=1d-2,smth_dist=3.15d-3,ls_eps
REAL*8 :: DGX=-9.81d0,DGY=0d0,DGZ=0d0 ! Gravity field
INTEGER, ALLOCATABLE :: interphase(:)
REAL*8 , ALLOCATABLE :: dnterphase(:)
INTEGER n_interphase

! -------------- Subroutines -------------------
CONTAINS
! ----------------------------------------------
SUBROUTINE Init_levelset

CALL Init_General_Scalar(levelset,Knpr_levelset)

ILEV=NLMAX
CALL SETLEV(2)

! Set dirichlet boundary conditions on the solution
CALL InitialBC_levelset(DWORK(L(KLCAG(ILEV))),&
levelset%l(ILEV)%val,levelset%l(ILEV)%ndof)

CALL Boundary_levelset_Val(DWORK(L(KLCAG(ILEV))))

CALL SETUP_LS(.true.)

! Reconstruct the normals after advecting the distance function
CALL GetNormal()
! Do the renormalization of the normals
CALL RenormalizeNormal(normal%l(ILEV)%val(1),&
     normal%l(ILEV)%val(NAT+1),normal%l(ILEV)%val(2*NAT+1))

levelset%cName = "LevSet"
levelset%NLmin = 2
levelset%NLmax = 20
levelset%defCrit =2d-4
levelset%epsCrit =1d-2
levelset%MinDef  =1d-16
levelset%SolvIter=25
levelset%SolvType=1
ls_eps=(2d0*smth_dist)**2d0

END SUBROUTINE Init_levelset
!
! ----------------------------------------------
!
SUBROUTINE Transport_levelset(mfile,INL)
INTEGER mfile
REAL*8  ResTemp,DefTemp,DefTempCrit,RhsTemp
INTEGER INL,INLComplete

!-------------------------------------------------
! Advecting the levelset function 
!-------------------------------------------------
! GOTO 1

IF (myid.NE.master) THEN

ILEV=NLMAX
CALL SETLEV(2)

CALL InitializeAFC_General_Scalar()

! Assemble the right hand side
CALL LCL1(levelset%l(ILEV)%def,levelset%l(ILEV)%ndof)
CALL Matdef_General_Scalar(levelset,1,0)

! Set dirichlet boundary conditions on the defect
CALL Boundary_levelset_Def()

! Store the constant right hand side
levelset%l(ILEV)%rhs = levelset%l(ILEV)%def

! Set dirichlet boundary conditions on the solution
CALL Boundary_levelset_Val(DWORK(L(KLCAG(ILEV))))

! Assemble the defect vector and fine level matrix
CALL Matdef_General_Scalar(levelset,-1,1)
CALL CommSum(levelset%l(NLEV)%def,NLEV)          ! PARALLEL

! Set dirichlet boundary conditions on the defect
CALL Boundary_levelset_Def()

! Save the old solution
CALL LCP1(levelset%l(NLEV)%val,levelset%l(NLEV)%val_old,&
     levelset%l(NLEV)%ndof)

! Compute the defect
CALL Resdfk_General_Scalar(levelset,ResTemp,DefTemp,RhsTemp)

END IF

CALL COMM_Maximum(RhsTemp)
DefTempCrit=MAX(RhsTemp*levelset%defCrit,levelset%MinDef)

CALL Protocol_General_Scalar(mfile,levelset,0,&
     ResTemp,DefTemp,DefTempCrit)

DO INL=1,levelset%NLmax
INLComplete = 0

! Solve the linear system
IF (myid.NE.master) THEN

! Calling the solver
CALL Solve_General_Scalar(levelset,Boundary_levelset_Val)

!!!!          Checking the quality of the result           !!!!
! Restore the constant right hand side
levelset%l(ILEV)%def = levelset%l(ILEV)%rhs

! Assemble the defect vector and fine level matrix
CALL Matdef_General_Scalar(levelset,-1,0)
CALL CommSum(levelset%l(NLEV)%def,NLEV)

! Set dirichlet boundary conditions on the defect
CALL Boundary_levelset_Def()

! Save the old solution
CALL LCP1(levelset%l(NLEV)%val,levelset%l(NLEV)%val_old,&
     levelset%l(NLEV)%ndof)

! Compute the defect
CALL Resdfk_General_Scalar(levelset,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
CALL COMM_Maximum(RhsTemp)
CALL Protocol_General_Scalar(mfile,levelset,INL,&
     ResTemp,DefTemp,RhsTemp)

IF ((DefTemp.LE.DefTempCrit).AND.&
    (INL.GE.levelset%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 1

END DO


1 CONTINUE

! Compute the deviation of the volume 
! fraction of the secondary phase
CALL Correct_levelset()
! GOTO 2

!-------------------------------------------------
! Levelset Reinitialization by solution of the
! hyperbolic PDE
!-------------------------------------------------

IF (myid.NE.master) THEN

ILEV=NLMAX
CALL SETLEV(2)

! levelset%l(ILEV)%val_old = levelset%l(ILEV)%val
! dnterphase = levelset%l(ilev)%val
! Reconstruct the normals after advecting the distance function
CALL GetNormal()
! Store the 'exact' interphase values
CALL GetInterphase(KWORK(L(LAREA)),KWORK(L(LADJ)))
! Do the renormalization of the normals
CALL RenormalizeNormal(normal%l(ILEV)%val(1),&
     normal%l(ILEV)%val(NAT+1),normal%l(ILEV)%val(2*NAT+1))

CALL InitializeAFC_Reinitialization()

! Assemble the right hand side
CALL LCL1(levelset%l(ILEV)%def,levelset%l(ILEV)%ndof)
CALL Matdef_General_Scalar(levelset,1,0)

! Assemble source/sink
CALL Source_reinitialization()

! Set dirichlet boundary conditions on the defect
CALL Boundary_reinitialization_Def()

! Store the constant right hand side
levelset%l(ILEV)%rhs = levelset%l(ILEV)%def

! Set dirichlet boundary conditions on the solution
CALL Boundary_reinitialization_Val(DWORK(L(KLCAG(ILEV))))

! Assemble the defect vector and fine level matrix
CALL Matdef_General_Scalar(levelset,-1,1)
CALL CommSum(levelset%l(NLEV)%def,NLEV)          ! PARALLEL

! Set dirichlet boundary conditions on the defect
CALL Boundary_reinitialization_Def()
CALL Boundary_reinitialization_mat(KWORK(L(KLLDA(NLEV))))

! Save the old solution
CALL LCP1(levelset%l(NLEV)%val,levelset%l(NLEV)%val_old,&
     levelset%l(NLEV)%ndof)

! Compute the defect
CALL Resdfk_General_Scalar(levelset,ResTemp,DefTemp,RhsTemp)

END IF

CALL COMM_Maximum(RhsTemp)
DefTempCrit=MAX(RhsTemp*levelset%defCrit,levelset%MinDef)

CALL Protocol_General_Scalar(mfile,levelset,0,&
     ResTemp,DefTemp,DefTempCrit)

DO INL=1,levelset%NLmax
INLComplete = 0

! Solve the linear system
IF (myid.NE.master) THEN

! Calling the solver
CALL Solve_General_Scalar(levelset,Boundary_reinitialization_Val)

!!!!          Checking the quality of the result           !!!!
! Restore the constant right hand side
levelset%l(ILEV)%def = levelset%l(ILEV)%rhs

! Assemble the defect vector and fine level matrix
CALL Matdef_General_Scalar(levelset,-1,0)
CALL CommSum(levelset%l(NLEV)%def,NLEV)

! Set dirichlet boundary conditions on the defect
CALL Boundary_reinitialization_Def()
CALL Boundary_reinitialization_mat(KWORK(L(KLLDA(NLEV))))

! Save the old solution
CALL LCP1(levelset%l(NLEV)%val,levelset%l(NLEV)%val_old,&
     levelset%l(NLEV)%ndof)

! Compute the defect
CALL Resdfk_General_Scalar(levelset,ResTemp,DefTemp,RhsTemp)

END IF

! Checking convergence rates against criterions
RhsTemp=DefTemp
CALL COMM_Maximum(RhsTemp)
CALL Protocol_General_Scalar(mfile,levelset,INL,&
     ResTemp,DefTemp,RhsTemp)

IF ((DefTemp.LE.DefTempCrit).AND.&
    (INL.GE.levelset%NLmin)) INLComplete = 1

CALL COMM_NLComplete(INLComplete)
IF (INLComplete.eq.1) GOTO 2

END DO

2 CONTINUE

! ! Reconstruct the normals after advecting the distance function
! CALL GetNormal()

!-------------------------------------------------
! Redefine the density function and rebuild the 
! pressure-poisson operator
!-------------------------------------------------
! CALL SETUP_LS(.false.)

END SUBROUTINE Transport_levelset
!
! ----------------------------------------------
!
SUBROUTINE Knpr_levelset(knpr,dcorag,nat)
INTEGER knpr(*),nat
REAL*8  dcorag(3,*),PX,PY,PZ

DO I=1,NAT
 PX = dcorag(1,I)
 PY = dcorag(2,I)
 PZ = dcorag(3,I)
 KNPR(I)=0
END DO

END SUBROUTINE Knpr_levelset
!
! ----------------------------------------------
!
SUBROUTINE Boundary_levelset_Val(dcorag)
REAL*8  dcorag(3,*),PX,PY,PZ
REAL*8 :: RX = 0.5d0,RY = 0.2d0, DIST

DO i=1,levelset%l(ILEV)%ndof
 PX = dcorag(1,I)
 PY = dcorag(2,I)
 PZ = dcorag(3,I)
 IF (levelset%l(ILEV)%knpr(i).eq.1) THEN
  CONTINUE
 END IF
END DO

END SUBROUTINE Boundary_levelset_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_levelset_Def

DO i=1,levelset%l(ILEV)%ndof
 IF (levelset%l(ILEV)%knpr(i).eq.1) THEN
  levelset%l(ILEV)%def(i) = 0d0
 END IF
END DO

END SUBROUTINE Boundary_levelset_Def
!
! ----------------------------------------------
!
SUBROUTINE GetInterphase(KAREA,KADJ)
REAL*8 DDD
INTEGER KAREA(6,*),KADJ(6,*)
INTEGER iel,iat,NJALFA,NIALFA

j = 0
DO iel=1,NEL
 NJALFA=0
 NIALFA=0
 DO iat=1,6
  i = KAREA(iat,iel)
  IF (levelset%l(ilev)%val(i).LE.0d0) THEN
   NJALFA=NJALFA+1
  ELSE
   NIALFA=NIALFA+1
  ENDIF
 END DO
 IF (.NOT.((NIALFA.EQ.6).OR.(NJALFA.EQ.6))) THEN
  DO iat=1,6
   IF (KADJ(iat,iel).LT.iel) THEN
    j = j + 1
    i = KAREA(iat,iel)
    interphase(j)=i
!     DDD = DSQRT(normal%l(ilev)%val(0*NAT+i)**2d0 + &
!                 normal%l(ilev)%val(1*NAT+i)**2d0 + &
!                 normal%l(ilev)%val(2*NAT+i)**2d0)
    dnterphase(j) = levelset%l(ilev)%val(i)!/DDD
   END IF
  END DO
 END IF
END DO

n_interphase = j

END SUBROUTINE GetInterphase
!
! ----------------------------------------------
!
SUBROUTINE Boundary_reinitialization_Val(dcorag)
REAL*8  dcorag(3,*),PX,PY,PZ

DO j=1,n_interphase
 i  = interphase(j)
 levelset%l(ILEV)%val(i) = dnterphase(j)
END DO

END SUBROUTINE Boundary_reinitialization_Val
!
! ----------------------------------------------
!
SUBROUTINE Boundary_reinitialization_Def

DO j=1,n_interphase
 i = interphase(j)
 levelset%l(ILEV)%def(i) = 0d0
END DO

END SUBROUTINE Boundary_reinitialization_Def
!
! ----------------------------------------------
!
SUBROUTINE Boundary_reinitialization_mat(KLD)
INTEGER KLD(*),ICOL

DO j=1,n_interphase
 i = interphase(j)
 levelset%l(ILEV)%Amat(KLD(I))=1E0
 DO ICOL=KLD(I)+1,KLD(I+1)-1
  levelset%l(ILEV)%Amat(ICOL)=0E0
 END DO
END DO

END SUBROUTINE Boundary_reinitialization_mat
!
! ----------------------------------------------
!
SUBROUTINE InitialBC_levelset(dcorvg,values,nu)
REAL*8 dcorvg(3,*),values(*)
REAL*8 X,Y,Z
REAL*8 :: RX = 0.5d0,RY = 0.205d0
REAL*8 :: RY1 = 0.155d0,RY2 = 0.255d0
REAL*8 :: RADx = 0.05d0, RADy = 0.1d0
REAL*8 DIST,dDist
INTEGER nu

DO i=1,nu
 X = dcorvg(1,i)
 Y = dcorvg(2,i)
 Z = dcorvg(3,i)

 IF (Y.GE.RY1.AND.Y.LE.RY2) THEN
  DIST = ABS(X-RX)
  values(i) = DIST - RADx  !+1d0
 ELSE
  IF (Y.LT.RY1) THEN
   DIST = DSQRT((RX-X)**2d0 + (RY1-Y)**2d0)
   values(i) = DIST - RADx  !+1d0
  ELSE
   DIST = DSQRT((RX-X)**2d0 + (RY2-Y)**2d0)
   values(i) = DIST - RADx  !+1d0
  END IF
 END IF


!  DIST = DSQRT((RX-X)**2d0 + (RY-Y)**2d0)
!  values(i) = DIST - RADx  !+1d0

END DO


END SUBROUTINE InitialBC_levelset
!
! ----------------------------------------------
!
SUBROUTINE GetNormal()

! Recovery of the normal direction
ILEV=NLMAX
CALL  SETLEV (2)
IF (IELT.EQ.0) THEN
!   CALL GETNORMALDG(levelset%l(ILEV)%val,Normal,KWORK(L(LVERT)),&
!        KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),&
!        VWORK(L(KLVOL(NLMAX))),ilev,E031)
END IF
IF (IELT.EQ.1) THEN
!   CALL GETNORMALDG(levelset%l(ILEV)%val,Normal,KWORK(L(LVERT)),&
!        KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),&
!        VWORK(L(KLVOL(NLMAX))),ilev,E030)
END IF
IF (IELT.EQ.2) THEN
  CALL GETNORMALNP(levelset%l(ILEV)%val,Normal%l(ILEV)%val,&
       KWORK(L(LVERT)),&
       KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),&
       VWORK(L(KLVOL(NLMAX))),ilev,EM31)
END IF
IF (IELT.EQ.3) THEN
  CALL GETNORMALNP(levelset%l(ILEV)%val,Normal%l(ILEV)%val,&
       KWORK(L(LVERT)),&
       KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),&
       VWORK(L(KLVOL(NLMAX))),ilev,EM30)
END IF

END SUBROUTINE GetNormal
!
! ----------------------------------------------
!
SUBROUTINE Source_reinitialization()

! IF (IELT.EQ.2) CALL RI_MASSDG(levelset%l(ILEV)%val,  &
!     levelset%l(ILEV)%def,NU,KWORK(KCOLA),KWORK(KLDA),&
!     KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)), &
!     KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),TSTEP,ls_eps,E031)
! IF (IELT.EQ.3) CALL RI_MASSDG(levelset%l(ILEV)%val,  &
!     levelset%l(ILEV)%def,NU,KWORK(KCOLA),KWORK(KLDA),&
!     KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)), &
!     KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),TSTEP,ls_eps,E030)
IF (IELT.EQ.2) CALL RI_MASSNP(levelset%l(ILEV)%val,  &
    levelset%l(ILEV)%def,NU,KWORK(KCOLA),KWORK(KLDA),&
    KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)), &
    KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),TSTEP,ls_eps,EM31)
IF (IELT.EQ.3) CALL RI_MASSNP(levelset%l(ILEV)%val,  &
    levelset%l(ILEV)%def,NU,KWORK(KCOLA),KWORK(KLDA),&
    KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)), &
    KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),TSTEP,ls_eps,EM30)

END SUBROUTINE Source_reinitialization
!
! ----------------------------------------------
!
SUBROUTINE GravForce(DML,RHSx,RHSy,RHSz)
REAL*8 DML(*),RHSx(*),RHSy(*),RHSz(*)

IF (Density_Primary.NE.Density_Secondary) THEN

 IF (DGX.NE.0d0) THEN
  DO I=1,NU
   IF (levelset%l(ILEV)%val(I).LT.0d0) THEN
    RHSx(I) = RHSx(I) - &
    TSTEP*DGX*DML(I)*(Density_Primary-Density_Secondary)
   END IF
  END DO
 END IF

 IF (DGY.NE.0d0) THEN
  DO I=1,NU
   IF (levelset%l(ILEV)%val(I).LT.0d0) THEN
    RHSy(I) = RHSy(I) - &
    TSTEP*DGY*DML(I)*(Density_Primary-Density_Secondary)
   END IF
  END DO
 END IF

 IF (DGZ.NE.0d0) THEN
  DO I=1,NU
   IF (levelset%l(ILEV)%val(I).LT.0d0) THEN
    RHSz(I) = RHSz(I) - &
    TSTEP*DGZ*DML(I)*(Density_Primary-Density_Secondary)
   END IF
  END DO
 END IF

END IF

END SUBROUTINE GravForce
!
! ----------------------------------------------
!
SUBROUTINE XSurfTens()
INTEGER IDEF

DO ILEV=NLMIN,NLMAX

 CALL SETLEV(2)

 IF (ILEV.EQ.NLMAX) THEN
  IDEF=1
 ELSE
  IDEF=0
 END IF

 IF (IELT.EQ.0) THEN
!   CALL SURFTENSDG(SurfT_Mat%l(ILEV)%val,DWORK(KF1),DWORK(KF2),&
!   DWORK(KF3),normal,levelset%l(ILEV)%val,smth_dist,&
!   SurfaceTension,KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),&
!   KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
!   DWORK(L(LCORVG)),KWORK(L(LADJ)),TSTEP,E031,IDEF)
 END IF
 IF (IELT.EQ.1) THEN
!   CALL SURFTENSDG(SurfT_Mat%l(ILEV)%val,DWORK(KF1),DWORK(KF2),&
!   DWORK(KF3),normal,levelset%l(ILEV)%val,smth_dist,&
!   SurfaceTension,KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),&
!   KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
!   DWORK(L(LCORVG)),KWORK(L(LADJ)),TSTEP,E030,IDEF)
 END IF
 IF (IELT.EQ.2) THEN
  CALL SURFTENSNP(SurfT_Mat%l(ILEV)%val,DWORK(KF1),DWORK(KF2),&
  DWORK(KF3),normal%l(ILEV)%val,levelset%l(ILEV)%val,&
  2d0*smth_dist,&
  SurfaceTension,KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),&
  KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
  DWORK(L(LCORVG)),KWORK(L(LADJ)),TSTEP,EM31,IDEF)
 END IF
 IF (IELT.EQ.3) THEN
  CALL SURFTENSNP(SurfT_Mat%l(ILEV)%val,DWORK(KF1),DWORK(KF2),&
  DWORK(KF3),normal%l(ILEV)%val,levelset%l(ILEV)%val,&
  2d0*smth_dist,&
  SurfaceTension,KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),&
  KWORK(L(LAREA)),KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),&
  DWORK(L(LCORVG)),KWORK(L(LADJ)),TSTEP,EM30,IDEF)
 END IF

END DO

ILEV=NLEV
CALL SETLEV(2)

END SUBROUTINE XSurfTens
!
! ----------------------------------------------
!
SUBROUTINE SETUP_LS(bInit)
LOGICAL bInit

IF (bInit) ALLOCATE (Density%l(NLMIN:NLMAX))
IF (bInit) ALLOCATE (SurfT_Mat%l(NLMIN:NLMAX))
IF (bInit) ALLOCATE (Normal%l(NLMIN:NLMAX))
IF (bInit) ALLOCATE (interphase(KNAT(NLMAX)))
IF (bInit) ALLOCATE (dnterphase(KNAT(NLMAX)))

DO ILEV=NLMAX,NLMIN,-1

 CALL SETLEV(2)

 IF (bInit) THEN
  ALLOCATE (Density%l(ILEV)%val(1:KNU(ILEV)))
  ALLOCATE (SurfT_Mat%l(ILEV)%val(1:KNA(ILEV)))
  ALLOCATE (Normal%l(ILEV)%val(1:3*KNAT(ILEV)))
 END IF

 IF (ILEV.NE.NLMIN) THEN
  IF(myid.ne.0) CALL RESTR_general_scalar(levelset)
  IF(myid.ne.0) CALL CommSumHalf(levelset%l(ILEV-1)%val,ILEV-1)
 END IF

 CALL GetDensity(Density%l(ILEV)%val,levelset%l(ILEV)%val,&
      levelset%l(ILEV)%ndof)

 CALL KNPRMPI(KWORK(L(KLNPR(ILEV))+NVT),2,ILEV) ! PARALLEL
 CALL LCL1  (DWORK(L(KLC(ILEV))),KNC(ILEV))
 CALL PROJMA(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),          &
             KWORK(L(KLLDC(ILEV))),KWORK(L(LNPR)),                &
             KWORK(L(LAREA)),KWORK(L(LADJ)),DWORK(L(KLM(ILEV))),  &
             Density%l(ILEV)%val,DWORK(KB1),DWORK(KB2),DWORK(KB3),&
             KWORK(KCOLB),KWORK(KLDB),KWORK(L(KLABD(ILEV))),ILEV)
 CALL KNPRMPI(KWORK(L(KLNPR(ILEV))+NVT),0,ILEV) ! PARALLEL

END DO

CONTAINS
! ----------------------------------------------
SUBROUTINE GetDensity(RHO,LS,nn)
REAL*8 RHO(*),LS(*)
INTEGER in,nn

DO in=1,nn
 IF (LS(in).LT.0d0) THEN
   RHO(in) = Density_Secondary
 ELSE
   RHO(in) = Density_Primary
 END IF
END DO

END SUBROUTINE GetDensity

END SUBROUTINE SETUP_LS
!
! ----------------------------------------------
!
SUBROUTINE InitializeAFC_Reinitialization()

! Construct the convection operator
CALL GetConvection_reinitialization(normal%l(ILEV)%val(1),&
     normal%l(ILEV)%val(NAT+1),normal%l(ILEV)%val(2*NAT+1),&
     DWORK(L(KLDK(ILEV))),levelset%l(ilev)%val,&
     KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LAREA)),&
     KWORK(L(LEDGE)),KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)))

CALL GetAFCStuctures_Scalar(DWORK(L(KLDK(ILEV))),KWORK(KCOLA),&
     KWORK(KLDA),NU,KWORK(L(LISEP)),KWORK(L(LIAUX)),&
     KWORK(L(LINOD)),KWORK(L(LJNOD)),DWORK(L(LAEDGE)),NEDGE)

END SUBROUTINE InitializeAFC_Reinitialization
!
! ----------------------------------------------
!
SUBROUTINE GetConvection_reinitialization(U1,U2,U3,DK,DLS,&
           KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG)
REAL*8  U1(*),U2(*),U3(*),DK(*),DCORVG(*),DLS(*)
INTEGER KCOLA(*),KLDA(*),KVERT(*),KAREA(*),KEDGE(*),KINT(*)
INTEGER NA

NA=KLDA(NU+1)-1
IF (IELT.EQ.0) CALL RI_CONVDG(U1,U2,U3,DLS,DK,NA,&
    KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,ls_eps,E031)
IF (IELT.EQ.1) CALL RI_CONVDG(U1,U2,U3,DLS,DK,NA,&
    KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,ls_eps,E030)
IF (IELT.EQ.2) CALL RI_CONVNP(U1,U2,U3,DLS,DK,NA,&
    KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,ls_eps,EM31)
IF (IELT.EQ.3) CALL RI_CONVNP(U1,U2,U3,DLS,DK,NA,&
    KCOLA,KLDA,KVERT,KAREA,KEDGE,KINT,DCORVG,ls_eps,EM30)

END SUBROUTINE GetConvection_reinitialization
!
! ----------------------------------------------
!
SUBROUTINE Correct_levelset()
REAL*8 DVOL,DSURF,DCORR

IF (myid.ne.master) THEN
 ILEV = NLMAX
 CALL SETLEV(2)

 IF (IELT.EQ.2) CALL LS_CORRNP(levelset%l(ilev)%val,&
   DVOL,DSURF,DCORR,KWORK(L(LVERT)),KWORK(L(LAREA)),&
   KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),&
   VWORK(L(KLVOL(ILEV))),EM31)
 IF (IELT.EQ.3) CALL LS_CORRNP(levelset%l(ilev)%val,&
   DVOL,DSURF,DCORR,KWORK(L(LVERT)),KWORK(L(LAREA)),&
   KWORK(L(LEDGE)),KWORK(L(KLINT(ILEV))),DWORK(L(LCORVG)),&
   VWORK(L(KLVOL(ILEV))),EM30)
END IF

CALL COMM_SUMM(DVOL)

IF (myid.eq.ShowID) WRITE(*,4)

IF (myid.eq.ShowID) WRITE(*,5) &
"Volume of secondary phase:", DVOL,"Fraction wrt initial:", DVOL/1.7853d-4

IF (myid.eq.ShowID) THEN
 OPEN(987,FILE='#data/Volume.txt',ACCESS = 'APPEND')
 WRITE(987,*) timens,dvol,DVOL/1.7853d-4
 CLOSE(987)
END IF

5  FORMAT(A26,1X,G9.3,1X,A21,1X,G11.5)
4  FORMAT(80('-'))

END SUBROUTINE Correct_levelset 
!
! ----------------------------------------------
!
SUBROUTINE RenormalizeNormal(n1,n2,n3)
REAL*8 n1(*),n2(*),n3(*),DDD

DO I=1,NAT
 DDD = DSQRT(n1(I)**2d0 + n2(I)**2d0 + n3(I)**2d0)
 n1(I) = n1(I)/DDD
 n2(I) = n2(I)/DDD
 n3(I) = n3(I)/DDD
END DO

END SUBROUTINE RenormalizeNormal

END MODULE
