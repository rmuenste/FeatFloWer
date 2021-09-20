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
  
  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'temp') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = dTemp
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 260d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.3) GenLinScalar%Fld(iFld)%val(i) = 250d0
  END IF

  IF (TRIM(GenLinScalar%Fld(iFld)%cName).eq.'alpha') then
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.1) GenLinScalar%Fld(iFld)%val(i) = 0d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.2) GenLinScalar%Fld(iFld)%val(i) = 1d0
   IF (GenLinScalar%Fld(iFld)%knpr(i).eq.3) GenLinScalar%Fld(iFld)%val(i) = 0.5d0
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

  IF (myBoundary%iInflow(i).eq.-1) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 1
  END IF
  
  IF (myBoundary%iInflow(i).eq.-2) THEN
    GenLinScalar%fld(iFld)%knpr(i) = 2
  END IF

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

  IF (myBoundary%LS_zero(i).ne.0) THEN
!     GenLinScalar%fld(iFld)%knpr(i) = 3
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
