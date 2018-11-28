MODULE def_QuadScalar

USE PP3D_MPI, ONLY:E011Sum,E011DMat,myid,showID,MGE013,COMM_SUMMN,&
                   COMM_Maximum,COMM_SUMM,COMM_NLComplete,&
                   myMPI_Barrier
USE var_QuadScalar
USE mg_QuadScalar, ONLY : MG_Solver,mgProlRestInit,mgProlongation,myMG,mgLev
USE UMFPackSolver, ONLY : myUmfPack_Factorize

use, intrinsic :: ieee_arithmetic

IMPLICIT NONE

CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Create_QuadMatStruct()
INTEGER iSymm,nERow
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA
INTEGER I,J,MatSize,NDOF
EXTERNAL E013,coefst

 ALLOCATE(mg_qMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX
   CALL SETLEV(2)

   ndof = mg_mesh%level(ilev)%nvt+&
          mg_mesh%level(ilev)%net+&
          mg_mesh%level(ilev)%nat+&
          mg_mesh%level(ilev)%nel

   MatSize = 300*NDOF

   ALLOCATE(TempColA(MatSize))
   ALLOCATE(mg_qMat(ILEV)%LdA(NDOF+1))
   mg_qMat(ILEV)%nu = NDOF
   mg_qMat(ILEV)%na = MatSize
   iSymm =   0
   nERow =   300

   CALL AP7(TempColA,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%na,&
            mg_qMat(ILEV)%nu,E013,iSymm,nERow,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea)

   IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
      "M,K,D,A matrix structure created",mg_qMat(ILEV)%nu,mg_qMat(ILEV)%na

   ALLOCATE(mg_qMat(ILEV)%ColA(mg_qMat(ILEV)%na))
   mg_qMat(ILEV)%ColA(:) = TempColA(1:mg_qMat(ILEV)%na)
   DEALLOCATE(TempColA)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)
 qMat => mg_qMat(NLMAX)

END SUBROUTINE Create_QuadMatStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_QuadLinMatStruct()
INTEGER iSymm,nERow
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA,TempLdA
INTEGER I,J,MatSize,NDOF
CHARACTER*10 myFile
EXTERNAL E011,E013,E010,coefst

 ALLOCATE (mg_qlMat(NLMIN:NLMAX))
 ALLOCATE (mg_lqMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX
  CALL SETLEV(2)

  NDOF = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

  MatSize = 16*27*mg_mesh%level(ilev)%nel

  ALLOCATE(TempColA(MatSize))
  ALLOCATE(mg_qlMat(ILEV)%LdA(NDOF+1))
  mg_qlMat(ILEV)%nu = NDOF
  mg_qlMat(ILEV)%na = MatSize
  iSymm =   0
  nERow =   16

  CALL AP9(TempColA,mg_qlMat(ILEV)%LdA,mg_qlMat(ILEV)%na,&
           mg_qlMat(ILEV)%nu,E013,E010,nERow,&
           mg_mesh%level(ILEV)%kvert,&
           mg_mesh%level(ILEV)%kedge,&
           mg_mesh%level(ILEV)%karea)

  mg_qlMat(ILEV)%na = 4*mg_qlMat(ILEV)%na

  ALLOCATE(mg_qlMat(ILEV)%ColA(mg_qlMat(ILEV)%na))

  CALL MatStructQ2P1(TempColA,mg_qlMat(ILEV)%ColA,mg_qlMat(ILEV)%LdA,&
       mg_qlMat(ILEV)%na,mg_qlMat(ILEV)%nu)

  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
   "B matrix structure created",mg_qlMat(ILEV)%nu,mg_qlMat(ILEV)%na

!   CALL OutputMatrixStuct("MatB",mg_qlMat(ILEV))

  MatSize = 4*27*mg_mesh%level(ilev)%nel
  ALLOCATE (TempLdA(4*mg_mesh%level(ilev)%nel+1))
  mg_lqMat(ILEV)%nu = mg_mesh%level(ilev)%nel
  mg_lqMat(ILEV)%na = MatSize
  iSymm =   0
  nERow =   27

  CALL AP9(TempColA,TempLdA,mg_lqMat(ILEV)%na,mg_lqMat(ILEV)%nu,E010,E013,nERow,&
           mg_mesh%level(ILEV)%kvert,&
           mg_mesh%level(ILEV)%kedge,&
           mg_mesh%level(ILEV)%karea)

  mg_lqMat(ILEV)%nu = 4*mg_lqMat(ILEV)%nu
  mg_lqMat(ILEV)%na = 4*mg_lqMat(ILEV)%na

  ALLOCATE(mg_lqMat(ILEV)%LdA(mg_lqMat(ILEV)%nu+1))
  ALLOCATE(mg_lqMat(ILEV)%ColA(mg_lqMat(ILEV)%na))

  CALL MatStructP1Q2(TempLdA,mg_lqMat(ILEV)%LdA,&
                     TempColA,&
                     mg_lqMat(ILEV)%ColA,&
                     MatSize,mg_mesh%level(ilev)%nel)

  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "BT matrix structure created",mg_lqMat(ILEV)%nu,mg_lqMat(ILEV)%na

  DEALLOCATE(TempColA,TempLdA)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat => mg_qlMat(NLMAX)
 lqMat => mg_lqMat(NLMAX)

END SUBROUTINE Create_QuadLinMatStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_LinMatStruct()

 ALLOCATE (mg_lMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlMat => mg_qlMat(ILEV)
  qlPMat => mg_qlPMat(ILEV)

  CALL Get_CMatLen(qlMat%LdA,qlMat%ColA,&
                   mg_mesh%level(ILEV)%kvert,&
                   mg_mesh%level(ilev)%nel,&
                   mg_lMat(ILEV)%na,1)



  mg_lMat(ILEV)%nu = 4*mg_mesh%level(ilev)%nel
  ALLOCATE (mg_lMat(ILEV)%LdA(mg_lMat(ILEV)%nu+1),mg_lMat(ILEV)%ColA(mg_lMat(ILEV)%na))

  mg_lMat(ILEV)%LdA=0
  mg_lMat(ILEV)%ColA=0

  CALL Get_CMatStruct(mg_lMat(ILEV)%LdA,mg_lMat(ILEV)%ColA,qlMat%LdA,qlMat%ColA,&
       mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%nel,mg_lMat(ILEV)%na,1)

  ! CALL OutputMatrixStuct("MatC",mg_lMat(ILEV))
  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "C matrix structure created",mg_lMat(ILEV)%nu,mg_lMat(ILEV)%na

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat => mg_qlMat(NLMAX)
 lMat  => mg_lMat(NLMAX)

END SUBROUTINE Create_LinMatStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_ParLinMatStruct()

 ALLOCATE (mg_lPMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlPMat => mg_qlPMat(ILEV)
  



  CALL Get_CMatLen(qlPMat%LdA,qlPMat%ColA,&
                   mg_mesh%level(ilev)%kvert,&
                   mg_mesh%level(ilev)%nel,&
                   mg_lPMat(ILEV)%na,2)

  mg_lPMat(ILEV)%nu = 4*NEL
  ALLOCATE (mg_lPMat(ILEV)%LdA(mg_lPMat(ILEV)%nu+1),mg_lPMat(ILEV)%ColA(mg_lPMat(ILEV)%na))
  mg_lPMat(ILEV)%LdA=0
  mg_lPMat(ILEV)%ColA=0

  CALL Get_CMatStruct(mg_lPMat(ILEV)%LdA,mg_lPMat(ILEV)%ColA,qlPMat%LdA,qlPMat%ColA,&
       mg_mesh%level(ilev)%kvert,mg_mesh%level(ilev)%nel,&
       mg_lPMat(ILEV)%na,2)

  ! CALL OutputMatrixStuct("MaPC",mg_lPMat(ILEV))
  IF (myid.eq.showID) WRITE(MTERM,'(A40,2I10)') &
  "Parallel C matrix structure created",mg_lPMat(ILEV)%nu,mg_lPMat(ILEV)%na
  ! pause
 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlPMat => mg_lqMat(NLMAX)
 lPMat => mg_lPMat(NLMAX)

END SUBROUTINE Create_ParLinMatStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_MRhoMat()
EXTERNAL E013
REAL*8  DML
INTEGER I,J

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Mmat))      ALLOCATE(mg_Mmat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_MlRhomat))  ALLOCATE(mg_MlRhomat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_MlRhoPmat)) ALLOCATE(mg_MlRhoPmat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Mmat(ILEV)%a)) THEN
   ALLOCATE(mg_Mmat(ILEV)%a(qMat%na))
  END IF

  mg_Mmat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [MRho] & [MlRho]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  CALL BuildMRhoMat(mgDensity(ILEV)%x,mg_Mmat(ILEV)%a,qMat%na,qMat%ColA,qMat%LdA,&
  mg_mesh%level(ILEV)%kvert,&
  mg_mesh%level(ILEV)%karea,&
  mg_mesh%level(ILEV)%kedge,&
  mg_mesh%level(ILEV)%dcorvg,&
  E013)

  if(bSteadyState)then
    mg_MMat(ILEV)%a = 0d0
  end if

  IF (.not.ALLOCATED(mg_MlRhomat(ILEV)%a)) ALLOCATE(mg_MlRhomat(ILEV)%a(qMat%nu))

!  IF (myid.eq.showID) WRITE(MTERM,*) "Assembling MLRho Matrix on Level [", ILEV,"]"
  DO I=1,qMat%nu
   DML = 0d0
   DO J=qMat%LdA(I),qMat%LdA(I+1)-1
    DML = DML + mg_Mmat(ILEV)%a(J)
   END DO
   mg_MlRhomat(ILEV)%a(I) = DML
  END DO

  IF (.not.ALLOCATED(mg_MlRhoPmat(ILEV)%a)) ALLOCATE(mg_MlRhoPmat(ILEV)%a(qMat%nu))
  mg_MlRhoPmat(ILEV)%a = mg_MlRhomat(ILEV)%a
  CALL E013SUM(mg_MlRhoPmat(ILEV)%a)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat      => mg_qMat(NLMAX)
 Mmat      => mg_Mmat(NLMAX)%a
 MlRhomat  => mg_MlRhomat(NLMAX)%a
 MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tMMat = myStat%tMMat + (myStat%t1-myStat%t0)
 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

END SUBROUTINE Create_MRhoMat
!
! ----------------------------------------------
!
SUBROUTINE Create_MMat()
EXTERNAL E013
REAL*8  DML
INTEGER I,J

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Mmat))  ALLOCATE(mg_Mmat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_MlMat))  ALLOCATE(mg_MlMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Mmat(ILEV)%a)) THEN
   ALLOCATE(mg_Mmat(ILEV)%a(qMat%na))
  ELSE
   mg_Mmat(ILEV)%a=0d0
  END IF


  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [MRho] & [MlRho]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL BuildMMat(mg_Mmat(ILEV)%a,qMat%na,qMat%ColA,qMat%LdA,&
                 mg_mesh%level(ILEV)%kvert,&
                 mg_mesh%level(ILEV)%karea,&
                 mg_mesh%level(ILEV)%kedge,&
                 mg_mesh%level(ILEV)%dcorvg,&
                 E013)


  IF (.not.ALLOCATED(mg_MlMat(ILEV)%a)) ALLOCATE(mg_MlMat(ILEV)%a(qMat%nu))

  DO I=1,qMat%nu
   DML = 0d0
   DO J=qMat%LdA(I),qMat%LdA(I+1)-1
    DML = DML + mg_Mmat(ILEV)%a(J)
   END DO
   mg_MlMat(ILEV)%a(I) = DML
  END DO

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat      => mg_qMat(NLMAX)
 Mmat      => mg_Mmat(NLMAX)%a
 MlMat     => mg_MlMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tMMat = myStat%tMMat + (myStat%t1-myStat%t0)
 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

END SUBROUTINE Create_MMat
!
! ----------------------------------------------
!
SUBROUTINE Create_CMat(knprU,knprV,knprW,knprP,coarse_lev,coarse_solver) !(C)
INTEGER coarse_lev,coarse_solver
TYPE(mg_kVector) :: knprU(*),knprV(*),knprW(*),knprP(*)
INTEGER i,j,iEntry,jCol

 CALL ZTIME(myStat%t0)

 IF (.NOT.ALLOCATED(mg_CMat)) ALLOCATE(mg_CMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat      => mg_qMat(ILEV)
  lMat      => mg_lMat(ILEV)
  qlMat     => mg_qlMat(ILEV)
  lqMat     => mg_lqMat(ILEV)
  BXMat     => mg_BXMat(ILEV)%a
  BYMat     => mg_BYMat(ILEV)%a
  BZMat     => mg_BZMat(ILEV)%a
  BTXMat    => mg_BTXMat(ILEV)%a
  BTYMat    => mg_BTYMat(ILEV)%a
  BTZMat    => mg_BTZMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a

  IF (.NOT.ALLOCATED(mg_CMat(ILEV)%a)) ALLOCATE(mg_CMat(ILEV)%a(lMat%na))
  mg_CMat(ILEV)%a=0d0
  CALL Get_CMat(MlRhoPmat,mg_CMat(ILEV)%a,lMat%LdA,lMat%ColA,&
       BXMat,BYMat,BZMat,qlMat%LdA,qlMat%ColA,&
       BTXMat,BTYMat,BTZMat,lqMat%LdA,lqMat%ColA, &
       knprU(ILEV)%x,knprV(ILEV)%x,knprW(ILEV)%x,lMat%nu,qMat%nu)

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B{T} MRho{-1} B]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
   !CALL OutputMatrix("MATC",lMat,mg_CMat(ILEV)%a,ILEV)
 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

!  WRITE(*,*) myid,mg_qMat(2)%nu
!  WRITE(outfile(8:9),'(A1,I1)') "K",myid
!  OPEN(FILE=outfile,UNIT=987)
!  WRITE(987,*) mg_lPMat(2)%ColA
!  CLOSE(987)
!  STOP

!   DO I=1,mg_qlMat(2)%na
!   WRITE(987,'(G12.4)') mg_qlMAt(2)%ColA(I)
!   WRITE(987,'(G12.4)') mg_lqMAt(2)%ColA(I)
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BXMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BYMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BZMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTXMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTYMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTZMat(2)%a(I)))
!   END DO
!   DO I=1,mg_lMat(2)%na
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_CMat(2)%a(I)))
!   END DO
!  CLOSE(987)
!  STOP

!  CALL OutputMatrix("cMAP",mg_lPMat(2),mg_CPMat(2)%a,2)
! !  CALL OutputMatrix("CMAT",mg_lMat(2),mg_CMat(2)%a,2)
!  STOP

 qMat      => mg_qMat(NLMAX)
 lMat      => mg_lMat(NLMAX)
 qlMat     => mg_qlMat(NLMAX)
 lqMat     => mg_lqMat(NLMAX)
 BXMat     => mg_BXMat(NLMAX)%a
 BYMat     => mg_BYMat(NLMAX)%a
 BZMat     => mg_BZMat(NLMAX)%a
 BTXMat    => mg_BTXMat(NLMAX)%a
 BTYMat    => mg_BTYMat(NLMAX)%a
 BTZMat    => mg_BTZMat(NLMAX)%a
 CMat      => mg_CMat(NLMAX)%a
 MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

 ! LU factorization of C-matrix needed for UMFPACK
 IF (myid.eq.0) THEN

  ILEV = coarse_lev
  CALL SETLEV(2)
  lMat      => mg_lMat(ILEV)
  CMat      => mg_CMat(ILEV)%a

  IF (bNoOutFlow) THEN
   DO I=lMat%LdA(1)+1,lMat%LdA(2)-1
    CMat(I) = 0d0
   END DO
  END IF

  IF (coarse_solver.EQ.2) THEN
   IF (.not.ALLOCATED(UMF_CMat)) ALLOCATE (UMF_CMat(lMat%na))
   UMF_CMat = CMat
   IF (.not.ALLOCATED(UMF_lMat%ColA)) ALLOCATE (UMF_lMat%ColA(lMat%na))
   IF (.not.ALLOCATED(UMF_lMat%LdA)) ALLOCATE (UMF_lMat%LdA(lMat%nu+1))
   UMF_lMat%ColA = lMat%ColA
   UMF_lMat%LdA  = lMat%LdA
   UMF_lMat%nu   = lMat%nu
   UMF_lMat%na   = lMat%na
   CALL myUmfPack_Factorize(UMF_CMat,UMF_lMat)
  END IF


  IF (coarse_solver.EQ.4.OR.coarse_solver.EQ.3) THEN

     crsSTR%A%nu = lMat%nu/4
     crsSTR%A%na = lMat%na/16   !!!! /16????????????????

     IF (.not.ALLOCATED(crsSTR%A%LdA)) ALLOCATE(crsSTR%A%Lda(crsSTR%A%nu+1))
     IF (.not.ALLOCATED(crsSTR%A%ColA)) ALLOCATE(crsSTR%A%ColA(crsSTR%A%na))
     IF (.not.ALLOCATED(crsSTR%A_MAT)) ALLOCATE(crsSTR%A_MAT(crsSTR%A%na))
     IF (.not.ALLOCATED(crsSTR%A_RHS)) ALLOCATE(crsSTR%A_RHS(crsSTR%A%nu))
     IF (.not.ALLOCATED(crsSTR%A_SOL)) ALLOCATE(crsSTR%A_SOL(crsSTR%A%nu))

     crsSTR%A%LdA(1) = 1
     DO i=1,crsSTR%A%nu
      j = 4*i-3
      crsSTR%A%LdA(i+1) = crsSTR%A%LdA(i) + (lMat%LdA(j+1)-lMat%LdA(j))/4
     END DO

     iEntry = 0
     DO i=1,crsSTR%A%nu
      j = 4*(i-1) + 1
      DO jCol = lMat%LdA(j),lMat%LdA(j+1)-1,4
       iEntry = iEntry + 1
       crsSTR%A%ColA(iEntry) = (lMat%ColA(jCol)-1)/4+1
       crsSTR%A_MAT(iEntry) = CMat(jCol)
      END DO
     END DO

     DO i=1,crsSTR%A%nu
      j = 4*(i-1) + 1
      if (KNPRP(ilev)%x(i).gt.0) crsSTR%A_Mat(crsSTR%A%LDA(i)) = 1d-8
     END DO
  
     IF (coarse_solver.EQ.4) THEN
      CALL myUmfPack_Factorize(crsSTR%A_MAT,crsSTR%A)
     END IF

  END IF

 END IF

 CALL ZTIME(myStat%t1)
 myStat%tCMat = myStat%tCMat + (myStat%t1-myStat%t0)
 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"


END SUBROUTINE Create_CMat!
! ----------------------------------------------
!
SUBROUTINE Create_ParCMat(knprU,knprV,knprW) !(C)
TYPE(mg_kVector) :: knprU(*),knprV(*),knprW(*)
INTEGER i

 CALL ZTIME(myStat%t0)

 IF (.NOT.ALLOCATED(mg_CPMat)) ALLOCATE(mg_CPMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat      => mg_qMat(ILEV)
  lPMat     => mg_lPMat(ILEV)
  qlPMat    => mg_qlPMat(ILEV)
  lqMat     => mg_lqMat(ILEV)
  BTXMat    => mg_BTXMat(ILEV)%a
  BTYMat    => mg_BTYMat(ILEV)%a
  BTZMat    => mg_BTZMat(ILEV)%a
  BXPMat    => mg_BXPMat(ILEV)%a
  BYPMat    => mg_BYPMat(ILEV)%a
  BZPMat    => mg_BZPMat(ILEV)%a
  MlRhoPmat => mg_MlRhoPmat(ILEV)%a

 IF (.NOT.ALLOCATED(mg_CPMat(ILEV)%a)) ALLOCATE(mg_CPMat(ILEV)%a(lPMat%na))
  mg_CPMat(ILEV)%a=0d0
  CALL Get_CMat(MlRhoPmat,mg_CPMat(ILEV)%a,lPMat%LdA,lPMat%ColA,&
       BXPMat,BYPMat,BZPMat,qlPMat%LdA,qlPMat%ColA,&
       BTXMat,BTYMat,BTZMat,lqMat%LdA,lqMat%ColA, &
       knprU(ILEV)%x,knprV(ILEV)%x,knprW(ILEV)%x,lPMat%nu,qMat%nu)

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " Par [B{T} MRho{-1} B] [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
   !CALL OutputMatrix("MAPC",lPMat,mg_CPMat(ILEV)%a,ILEV)
 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

!  WRITE(*,*) myid,mg_qMat(2)%nu
!  WRITE(outfile(8:9),'(A1,I1)') "K",myid
!  OPEN(FILE=outfile,UNIT=987)
!  WRITE(987,*) mg_lPMat(2)%ColA
!  CLOSE(987)
!  STOP

!   DO I=1,mg_qlMat(2)%na
!   WRITE(987,'(G12.4)') mg_qlMAt(2)%ColA(I)
!   WRITE(987,'(G12.4)') mg_lqMAt(2)%ColA(I)
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BXMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BYMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BZMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTXMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTYMat(2)%a(I)))
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_BTZMat(2)%a(I)))
!   END DO
!   DO I=1,mg_lMat(2)%na
!   WRITE(987,'(G12.4)') MAX(1d-18,ABS(mg_CMat(2)%a(I)))
!   END DO
!  CLOSE(987)
!  STOP

!  CALL OutputMatrix("cMAP",mg_lPMat(2),mg_CPMat(2)%a,2)
! !  CALL OutputMatrix("CMAT",mg_lMat(2),mg_CMat(2)%a,2)
!  STOP

 qMat      => mg_qMat(NLMAX)
 lPMat     => mg_lPMat(NLMAX)
 qlPMat    => mg_qlPMat(NLMAX)
 lqMat     => mg_lqMat(NLMAX)
 BTXMat    => mg_BTXMat(NLMAX)%a
 BTYMat    => mg_BTYMat(NLMAX)%a
 BTZMat    => mg_BTZMat(NLMAX)%a
 BXPMat    => mg_BXPMat(NLMAX)%a
 BYPMat    => mg_BYPMat(NLMAX)%a
 BZPMat    => mg_BZPMat(NLMAX)%a
 CPMat     => mg_CPMat(NLMAX)%a
 MlRhoPmat => mg_MlRhoPmat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tCMat = myStat%tCMat + (myStat%t1-myStat%t0)
 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

END SUBROUTINE Create_ParCMat
!
! ----------------------------------------------
!
SUBROUTINE Create_BMat() !(B)
INTEGER nERow,pNEL
INTEGER I,J
real*8 ddx,ddy,ddz
CHARACTER*10 myFile
EXTERNAL E011,E013

 IF (.NOT.ALLOCATED(mg_BXMat)) ALLOCATE(mg_BXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BYMat)) ALLOCATE(mg_BYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BZMat)) ALLOCATE(mg_BZMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  qlMat => mg_qlMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BXMat(ILEV)%a)) ALLOCATE(mg_BXMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BYMat(ILEV)%a)) ALLOCATE(mg_BYMat(ILEV)%a(qlMat%na))
  IF (.NOT.ALLOCATED(mg_BZMat(ILEV)%a)) ALLOCATE(mg_BZMat(ILEV)%a(qlMat%na))
  mg_BXMat(ILEV)%a=0d0
  mg_BYMat(ILEV)%a=0d0
  mg_BZMat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BMatP1(mg_BXMat(ILEV)%a,mg_BYMat(ILEV)%a,&
       mg_BZMat(ILEV)%a,qlMat%LdA,qlMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       qlMat%na,E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 IF (.NOT.ALLOCATED(mg_BTXMat)) ALLOCATE(mg_BTXMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTYMat)) ALLOCATE(mg_BTYMat(NLMIN:NLMAX))
 IF (.NOT.ALLOCATED(mg_BTZMat)) ALLOCATE(mg_BTZMat(NLMIN:NLMAX))

DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)

  lqMat => mg_lqMat(ILEV)
  IF (.NOT.ALLOCATED(mg_BTXMat(ILEV)%a)) ALLOCATE(mg_BTXMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTYMat(ILEV)%a)) ALLOCATE(mg_BTYMat(ILEV)%a(lqMat%na))
  IF (.NOT.ALLOCATED(mg_BTZMat(ILEV)%a)) ALLOCATE(mg_BTZMat(ILEV)%a(lqMat%na))
  mg_BTXMat(ILEV)%a=0d0
  mg_BTYMat(ILEV)%a=0d0
  mg_BTZMat(ILEV)%a=0d0


  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [B{T}]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
  CALL Build_BTMatP1(mg_BTXMat(ILEV)%a,mg_BTYMat(ILEV)%a,&
       mg_BTZMat(ILEV)%a,lqMat%LdA,lqMat%ColA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       lqMat%na,E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

!  CALL OutputMatrix("bXMM",mg_qlMat(2),mg_BXMat(2)%a,2)
!  CALL OutputMatrix("bYMM",mg_qlMat(2),mg_BYMat(2)%a,2)
!  CALL OutputMatrix("bZMM",mg_qlMat(2),mg_BZMat(2)%a,2)
!  CALL OutputMatrix("bXTM",mg_lqMat(2),mg_BTXMat(2)%a,2)
!  CALL OutputMatrix("bYTM",mg_lqMat(2),mg_BTYMat(2)%a,2)
!  CALL OutputMatrix("bZTM",mg_lqMat(2),mg_BTZMat(2)%a,2)
!  STOP

55 CONTINUE

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat  => mg_qlMat(NLMAX)
 BXMat => mg_BXMat(NLMAX)%a
 BYMat => mg_BYMat(NLMAX)%a
 BZMat => mg_BZMat(NLMAX)%a

 lqMat  => mg_lqMat(NLMAX)
 BTXMat => mg_BTXMat(NLMAX)%a
 BTYMat => mg_BTYMat(NLMAX)%a
 BTZMat => mg_BTZMat(NLMAX)%a

END SUBROUTINE Create_BMat
!
! ----------------------------------------------
!
SUBROUTINE BuildUpPressureCommunicator(lScalar,lPScalar)
TYPE(TLinScalar), INTENT(INOUT), TARGET :: lScalar
TYPE(TParLinScalar), INTENT(INOUT), TARGET ::  lPScalar

 DO ILEV=NLMIN,NLMAX
  IF (myid.ne.0) THEN
   CALL Create_GlobalP1CommNumbering(lScalar,lPScalar,&
                      mg_mesh%level(ILEV)%dcorvg,&
                      mg_mesh%level(ILEV)%kvert,&
                      mg_mesh%level(ILEV)%kedge,&
                      mg_mesh%level(ILEV)%karea,&
                      mg_mesh%level(ILEV)%nvt,&
                      mg_mesh%level(ILEV)%net,&
                      mg_mesh%level(ILEV)%nat,&
                      mg_mesh%level(ILEV)%nel)   
  end if
  end do

END SUBROUTINE BuildUpPressureCommunicator
!
! ----------------------------------------------
!
SUBROUTINE Create_QuadLinParMatStruct(myPLinSc) !(B)
TYPE(TParLinScalar) myPLinSc
INTEGER pNEL,MatSize
INTEGER, ALLOCATABLE :: TempLdB(:)

 ALLOCATE (mg_BXPMat(NLMIN:NLMAX))
 ALLOCATE (mg_BYPMat(NLMIN:NLMAX))
 ALLOCATE (mg_BZPMat(NLMIN:NLMAX))
 ALLOCATE (mg_qlPMat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlMat => mg_qlMat(ILEV)
  BXMat => mg_BXMat(ILEV)%a
  BYMat => mg_BYMat(ILEV)%a
  BZMat => mg_BZMat(ILEV)%a

  CALL ParPresComm_Init(qlMat%ColA,qlMat%LdA,qlMat%nu,mg_mesh%level(ILEV)%nel,ILEV)

  mg_qlPMat(ILEV)%nu = qlMat%nu
  ALLOCATE(mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1))
  mg_qlPMat(ILEV)%LdA = 0
  ! Get the Par_LD_B structure
  CALL Create_ParB_LD(mg_qlPMat(ILEV)%LdA,qlMat%LdA,qlMat%nu,ILEV)
  MatSize = mg_qlPMat(ILEV)%LdA(mg_qlPMat(ILEV)%nu+1)
  IF (myid.eq.showID) WRITE(MTERM,*) "Parallel B Matrix and Structure created", &
  mg_qlPMat(ILEV)%nu,MatSize
  ALLOCATE(mg_qlPMat(ILEV)%ColA(MatSize))
  mg_qlPMat(ILEV)%ColA = 0

  ALLOCATE(mg_BXPMat(ILEV)%a(MatSize))
  ALLOCATE(mg_BYPMat(ILEV)%a(MatSize))
  ALLOCATE(mg_BZPMat(ILEV)%a(MatSize))

  ALLOCATE(TempLdB(mg_qlPMat(ILEV)%nu+1))
  TempLdB = mg_qlPMat(ILEV)%LdA

  ! Get the Par_COL_B structure
  CALL Create_ParB_COLMAT(BXMat,BYMat,BZMat,&
       mg_BXPMat(ILEV)%a,mg_BYPMat(ILEV)%a,mg_BZPMat(ILEV)%a,&
       mg_qlPMat(ILEV)%LdA,mg_qlPMat(ILEV)%ColA,TempLdB,&
       qlMat%LdA,qlMat%ColA,qlMat%nu,&
       mg_mesh%level(ILEV)%nel,&
       pNEL,ILEV)

  DEALLOCATE(TempLdB)

 END DO

 myPLinSc%ndof=pNEL
 IF (ALLOCATED(myPLinSc%Val)) DEALLOCATE(myPLinSc%Val)
 ALLOCATE (myPLinSc%Val(4*pNEL))

! CALL OutputMatrixStuct("BPm3",mg_qlPMat(3))
!  CALL OutputMatrix("bXPM",mg_qlPMat(2),mg_BXPMat(2)%a,2)
!  CALL OutputMatrix("bYPM",mg_qlPMat(2),mg_BYPMat(2)%a,2)
!  CALL OutputMatrix("bZPM",mg_qlPMat(2),mg_BZPMat(2)%a,2)
!  STOP

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat  => mg_qlMat(NLMAX)
 qlPMat => mg_qlPMat(NLMAX)
 BXPMat => mg_BXPMat(NLMAX)%a
 BYPMat => mg_BYPMat(NLMAX)%a
 BZPMat => mg_BZPMat(NLMAX)%a
 BXMat  => mg_BXMat(NLMAX)%a
 BYMat  => mg_BYMat(NLMAX)%a
 BZMat  => mg_BZMat(NLMAX)%a

END SUBROUTINE Create_QuadLinParMatStruct
!
! ----------------------------------------------
!
SUBROUTINE Fill_QuadLinParMat() !(B)
INTEGER pNEL
INTEGER, ALLOCATABLE :: TempLdB(:)

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qlMat => mg_qlMat(ILEV)
  BXMat => mg_BXMat(ILEV)%a
  BYMat => mg_BYMat(ILEV)%a
  BZMat => mg_BZMat(ILEV)%a

  ALLOCATE(TempLdB(mg_qlPMat(ILEV)%nu+1))
  TempLdB = mg_qlPMat(ILEV)%LdA

  CALL Create_ParB_COLMAT(BXMat,BYMat,BZMat,&
       mg_BXPMat(ILEV)%a,mg_BYPMat(ILEV)%a,mg_BZPMat(ILEV)%a,&
       mg_qlPMat(ILEV)%LdA,mg_qlPMat(ILEV)%ColA,TempLdB,&
       qlMat%LdA,qlMat%ColA,qlMat%nu,NEL,pNEL,ILEV)

  DEALLOCATE(TempLdB)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qlMat  => mg_qlMat(NLMAX)
 qlPMat => mg_qlPMat(NLMAX)
 BXPMat => mg_BXPMat(NLMAX)%a
 BYPMat => mg_BYPMat(NLMAX)%a
 BZPMat => mg_BZPMat(NLMAX)%a
 BXMat  => mg_BXMat(NLMAX)%a
 BYMat  => mg_BYMat(NLMAX)%a
 BZMat  => mg_BZMat(NLMAX)%a

END SUBROUTINE Fill_QuadLinParMat
!
! ----------------------------------------------
!
SUBROUTINE Create_DiffMat(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Dmat)) ALLOCATE(mg_Dmat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Dmat(ILEV)%a)) THEN
   ALLOCATE(mg_Dmat(ILEV)%a(qMat%na))
  END IF

  mg_Dmat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [D]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  if(bNonNewtonian) THEN
    CALL DIFFQ2_NNEWT(myScalar%valU, myScalar%valV,myScalar%valW, &
         Temperature,&
         mg_Dmat(ILEV)%a,qMat%na,qMat%ColA,qMat%LdA,&
         mg_mesh%level(ILEV)%kvert,&
         mg_mesh%level(ILEV)%karea,&
         mg_mesh%level(ILEV)%kedge,&
         mg_mesh%level(ILEV)%dcorvg,&
         E013)
  else 
    CALL DIFFQ2_NEWT(mg_Dmat(ILEV)%a,qMat%na,qMat%ColA,&
         qMat%LdA,&
         mg_mesh%level(ILEV)%kvert,&
         mg_mesh%level(ILEV)%karea,&
         mg_mesh%level(ILEV)%kedge,&
         mg_mesh%level(ILEV)%dcorvg,&
         E013)
  end if

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 DMat => mg_DMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tDMat = myStat%tDMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_DiffMat
!
! ----------------------------------------------
!
SUBROUTINE Create_SMat(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013
INTEGER i

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_S11mat)) ALLOCATE(mg_S11mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S22mat)) ALLOCATE(mg_S22mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S33mat)) ALLOCATE(mg_S33mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S12mat)) ALLOCATE(mg_S12mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S13mat)) ALLOCATE(mg_S13mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S23mat)) ALLOCATE(mg_S23mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S21mat)) ALLOCATE(mg_S21mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S31mat)) ALLOCATE(mg_S31mat(NLMIN:NLMAX))
 IF (.not.ALLOCATED(mg_S32mat)) ALLOCATE(mg_S32mat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_S11mat(ILEV)%a)) ALLOCATE(mg_S11mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S22mat(ILEV)%a)) ALLOCATE(mg_S22mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S33mat(ILEV)%a)) ALLOCATE(mg_S33mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S12mat(ILEV)%a)) ALLOCATE(mg_S12mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S13mat(ILEV)%a)) ALLOCATE(mg_S13mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S23mat(ILEV)%a)) ALLOCATE(mg_S23mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S21mat(ILEV)%a)) ALLOCATE(mg_S21mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S31mat(ILEV)%a)) ALLOCATE(mg_S31mat(ILEV)%a(qMat%na))
  IF (.not.ALLOCATED(mg_S32mat(ILEV)%a)) ALLOCATE(mg_S32mat(ILEV)%a(qMat%na))

  mg_S11mat(ILEV)%a=0d0
  mg_S22mat(ILEV)%a=0d0
  mg_S33mat(ILEV)%a=0d0
  mg_S12mat(ILEV)%a=0d0
  mg_S13mat(ILEV)%a=0d0
  mg_S23mat(ILEV)%a=0d0
  mg_S21mat(ILEV)%a=0d0
  mg_S31mat(ILEV)%a=0d0
  mg_S32mat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [S]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF
!  IF (myid.eq.showID) WRITE(MTERM,*) "Assembling S Matrix on Level", ILEV

  CALL CUBATURESTRESS(myScalar%valU, myScalar%valV,myScalar%valW, &
       Temperature,&
       mg_S11mat(ILEV)%a,mg_S22mat(ILEV)%a,mg_S33mat(ILEV)%a,&
       mg_S12mat(ILEV)%a,mg_S13mat(ILEV)%a,mg_S23mat(ILEV)%a,&
       mg_S21mat(ILEV)%a,mg_S31mat(ILEV)%a,mg_S32mat(ILEV)%a,&
       qMat%na,qMat%ColA,qMat%LdA,&
       mg_mesh%level(ILEV)%kvert,&
       mg_mesh%level(ILEV)%karea,&
       mg_mesh%level(ILEV)%kedge,&
       mg_mesh%level(ILEV)%dcorvg,&
       E013)

 END DO

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 S11Mat => mg_S11Mat(NLMAX)%a
 S22Mat => mg_S22Mat(NLMAX)%a
 S33Mat => mg_S33Mat(NLMAX)%a
 S12Mat => mg_S12Mat(NLMAX)%a
 S13Mat => mg_S13Mat(NLMAX)%a
 S23Mat => mg_S23Mat(NLMAX)%a
 S21Mat => mg_S21Mat(NLMAX)%a
 S31Mat => mg_S31Mat(NLMAX)%a
 S32Mat => mg_S32Mat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tSMat = myStat%tSMat + (myStat%t1-myStat%t0)

END SUBROUTINE Create_SMat
!
! ----------------------------------------------
!
SUBROUTINE Create_KMat(myScalar)
TYPE(TQuadScalar) myScalar
INTEGER LINT
EXTERNAL E013
! Assembly for Convection Kmat

 CALL ZTIME(myStat%t0)

 IF (.not.ALLOCATED(mg_Kmat)) ALLOCATE(mg_Kmat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

  CALL SETLEV(2)
  qMat => mg_qMat(ILEV)

  IF (.not.ALLOCATED(mg_Kmat(ILEV)%a)) THEN
   ALLOCATE(mg_Kmat(ILEV)%a(qMat%na))
  END IF

  mg_Kmat(ILEV)%a=0d0

  IF (myid.eq.showID) THEN
   IF (ILEV.EQ.NLMIN) THEN
    WRITE(MTERM,'(A,I1,A)', advance='no') " [KRho]: [", ILEV,"]"
   ELSE
    WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
   END IF
  END IF

  LINT = KLINT(ILEV)

  CALL CONVQ2(mgDensity(ILEV)%x,myScalar%valU,myScalar%valV,myScalar%valW,&
  myALE%MeshVelo,&
  mg_Kmat(ILEV)%a,qMat%nu,qMat%ColA,qMat%LdA,&
  mg_mesh%level(ILEV)%kvert,&
  mg_mesh%level(ILEV)%karea,&
  mg_mesh%level(ILEV)%kedge,&
  mg_mesh%level(ILEV)%dcorvg,&
  E013)

 END DO

 IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='no') " |"

 ILEV=NLMAX
 CALL SETLEV(2)

 qMat  => mg_qMat(NLMAX)
 KMat => mg_KMat(NLMAX)%a

 CALL ZTIME(myStat%t1)
 myStat%tKMat = myStat%tKMat + (myStat%t1-myStat%t0)
!  WRITE(*,*) myid,myStat%tKMat,myStat%t1-myStat%t0

END SUBROUTINE Create_KMat
!
! ----------------------------------------------
!
SUBROUTINE AddStressToRHS(myScalar)
TYPE(TQuadScalar) myScalar
EXTERNAL E013

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)

 CALL STRESS(myScalar%valU,myScalar%valV,myScalar%valW,&
 Temperature,&
 myScalar%defU, myScalar%defV, myScalar%defW,&
 Viscosity,&
 mg_mesh%level(ILEV)%kvert,&
 mg_mesh%level(ILEV)%karea,&
 mg_mesh%level(ILEV)%kedge,&
 mg_mesh%level(ILEV)%dcorvg,&
 E013 ) 

END IF

END SUBROUTINE AddStressToRHS
!
! ----------------------------------------------
!
SUBROUTINE GetGradVelo_rhs(myScalar,x)
TYPE(TQuadScalar) myScalar
REAL*8 x(*)
EXTERNAL E013

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)

 myScalar%defU=0.0d0
 myScalar%defV=0.0d0
 myScalar%defW=0.0d0

 CALL GetGradVelo_rhs_sub(x,myScalar%defU,myScalar%defV,myScalar%defW,&
                          mg_mesh%level(ilev)%kvert,&
                          mg_mesh%level(ilev)%karea,&
                          mg_mesh%level(ilev)%kedge,&
                          mg_mesh%level(ilev)%dcorvg,&
                          E013)

END IF

END SUBROUTINE GetGradVelo_rhs
!
! ----------------------------------------------
!
SUBROUTINE GetGradVelo_val(myScalar,iComp,dRho)
TYPE(TQuadScalar) myScalar
REAL*8 dRho
INTEGER iComp,i

IF (myid.ne.0) THEN

 ILEV=NLMAX
 CALL SETLEV(2)
 MlRhoPMat => mg_MlRhoPmat(ILEV)%a

 IF (iComp.EQ.1) THEN
  DO i=1,myScalar%ndof
   myScalar%valUx(i) = dRho*myScalar%defU(i)/MlRhoPMat(i)
   myScalar%valUy(i) = dRho*myScalar%defV(i)/MlRhoPMat(i)
   myScalar%valUz(i) = dRho*myScalar%defW(i)/MlRhoPMat(i)
  END DO
 END IF

 IF (iComp.EQ.2) THEN
  DO i=1,myScalar%ndof
   myScalar%valVx(i) = dRho*myScalar%defU(i)/MlRhoPMat(i)
   myScalar%valVy(i) = dRho*myScalar%defV(i)/MlRhoPMat(i)
   myScalar%valVz(i) = dRho*myScalar%defW(i)/MlRhoPMat(i)
  END DO
 END IF

 IF (iComp.EQ.3) THEN
  DO i=1,myScalar%ndof
   myScalar%valWx(i) = dRho*myScalar%defU(i)/MlRhoPMat(i)
   myScalar%valWy(i) = dRho*myScalar%defV(i)/MlRhoPMat(i)
   myScalar%valWz(i) = dRho*myScalar%defW(i)/MlRhoPMat(i)
  END DO
 END IF

END IF

END SUBROUTINE GetGradVelo_val
!
!----------------------------------------------
!
SUBROUTINE Create_AMat()
INTEGER NA,complete

IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  ALLOCATE(mg_A11mat(NLMIN:NLMAX))
  ALLOCATE(mg_A22mat(NLMIN:NLMAX))
  ALLOCATE(mg_A33mat(NLMIN:NLMAX))
  ALLOCATE(mg_A12mat(NLMIN:NLMAX))
  ALLOCATE(mg_A13mat(NLMIN:NLMAX))
  ALLOCATE(mg_A23mat(NLMIN:NLMAX))
  ALLOCATE(mg_A21mat(NLMIN:NLMAX))
  ALLOCATE(mg_A31mat(NLMIN:NLMAX))
  ALLOCATE(mg_A32mat(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') "Allocation of Aii Matrix on Level [", ILEV,"]"
    END IF
    IF (ILEV.EQ.NLMAX) THEN
     WRITE(MTERM,'(A,I1,A)', advance='yes') ", [",ILEV,"]"
    END IF
    IF (ILEV.NE.NLMAX.AND.ILEV.NE.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF

   NA = mg_qMat(ILEV)%na
   ALLOCATE(mg_A11mat(ILEV)%a(NA))
   ALLOCATE(mg_A22mat(ILEV)%a(NA))
   ALLOCATE(mg_A33mat(ILEV)%a(NA))
   ALLOCATE(mg_A12mat(ILEV)%a(NA))
   ALLOCATE(mg_A13mat(ILEV)%a(NA))
   ALLOCATE(mg_A23mat(ILEV)%a(NA))
   ALLOCATE(mg_A21mat(ILEV)%a(NA))
   ALLOCATE(mg_A31mat(ILEV)%a(NA))
   ALLOCATE(mg_A32mat(ILEV)%a(NA))
  END DO
  A11Mat => mg_A11mat(NLMAX)%a
  A22Mat => mg_A22mat(NLMAX)%a
  A33Mat => mg_A33mat(NLMAX)%a
  A12Mat => mg_A12mat(NLMAX)%a
  A13Mat => mg_A13mat(NLMAX)%a
  A23Mat => mg_A23mat(NLMAX)%a
  A21Mat => mg_A21mat(NLMAX)%a
  A31Mat => mg_A31mat(NLMAX)%a
  A32Mat => mg_A32mat(NLMAX)%a
 ELSE
  ALLOCATE(mg_A11mat(NLMIN:NLMAX))
  ALLOCATE(mg_A22mat(NLMIN:NLMAX))
  ALLOCATE(mg_A33mat(NLMIN:NLMAX))
  DO ILEV=NLMIN,NLMAX

   IF (myid.eq.showID) THEN
    IF (ILEV.EQ.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') "Allocation of A Matrix on Level [", ILEV,"]"
    END IF
    IF (ILEV.EQ.NLMAX) THEN
     WRITE(MTERM,'(A,I1,A)', advance='yes') ", [",ILEV,"]"
    END IF
    IF (ILEV.NE.NLMAX.AND.ILEV.NE.NLMIN) THEN
     WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
    END IF
   END IF

   NA = mg_qMat(ILEV)%na
   ALLOCATE(mg_A11mat(ILEV)%a(NA))
   ALLOCATE(mg_A22mat(ILEV)%a(NA))
   ALLOCATE(mg_A33mat(ILEV)%a(NA))
  END DO
  A11Mat => mg_A11mat(NLMAX)%a
  A22Mat => mg_A22mat(NLMAX)%a
  A33Mat => mg_A33mat(NLMAX)%a
 END IF

END SUBROUTINE Create_AMat
!
! ----------------------------------------------
!
SUBROUTINE InitializeQuadScalar(myScalar)
TYPE(TQuadScalar) myScalar
INTEGER NDOF,N

 ALLOCATE(myScalar%knprU(NLMIN:NLMAX))
 ALLOCATE(myScalar%knprV(NLMIN:NLMAX))
 ALLOCATE(myScalar%knprW(NLMIN:NLMAX))

 DO ILEV=NLMIN,NLMAX
  !NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)

  NDOF = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

  ALLOCATE(myScalar%knprU(ILEV)%x(NDOF))
  ALLOCATE(myScalar%knprV(ILEV)%x(NDOF))
  ALLOCATE(myScalar%knprW(ILEV)%x(NDOF))
 END DO

 ILEV=NLMAX

 !NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)
 NDOF = mg_mesh%level(ilev)%nvt+&
        mg_mesh%level(ilev)%net+&
        mg_mesh%level(ilev)%nat+&
        mg_mesh%level(ilev)%nel

 myScalar%ndof = NDOF
 IF (ALLOCATED(mg_qMat)) myScalar%na = mg_qMat(ILEV)%na

 ! U component
 ALLOCATE(myScalar%auxU(myScalar%ndof))
 ALLOCATE(myScalar%rhsU(myScalar%ndof))
 ALLOCATE(myScalar%defU(myScalar%ndof))
 ALLOCATE(myScalar%valU_old(myScalar%ndof))
 ALLOCATE(myScalar%valU(myScalar%ndof))

 ! V component
 ALLOCATE(myScalar%auxV(myScalar%ndof))
 ALLOCATE(myScalar%rhsV(myScalar%ndof))
 ALLOCATE(myScalar%defV(myScalar%ndof))
 ALLOCATE(myScalar%valV_old(myScalar%ndof))
 ALLOCATE(myScalar%valV(myScalar%ndof))

 ! W component
 ALLOCATE(myScalar%auxW(myScalar%ndof))
 ALLOCATE(myScalar%rhsW(myScalar%ndof))
 ALLOCATE(myScalar%defW(myScalar%ndof))
 ALLOCATE(myScalar%valW_old(myScalar%ndof))
 ALLOCATE(myScalar%valW(myScalar%ndof))

 ! Derivatives
 ALLOCATE(myScalar%valUx(myScalar%ndof))
 ALLOCATE(myScalar%valUy(myScalar%ndof))
 ALLOCATE(myScalar%valUz(myScalar%ndof))
 ALLOCATE(myScalar%valVx(myScalar%ndof))
 ALLOCATE(myScalar%valVy(myScalar%ndof))
 ALLOCATE(myScalar%valVz(myScalar%ndof))
 ALLOCATE(myScalar%valWx(myScalar%ndof))
 ALLOCATE(myScalar%valWy(myScalar%ndof))
 ALLOCATE(myScalar%valWz(myScalar%ndof))

 ! Prolongation Restriction Matrices
 ALLOCATE(mg_E013Prol(NLMIN:NLMAX-1))
 ALLOCATE(mg_E013Rest(NLMIN:NLMAX-1))
 ALLOCATE(mg_E013ProlM(NLMIN:NLMAX-1))
 ALLOCATE(mg_E013RestM(NLMIN:NLMAX-1))

 DO ILEV=NLMIN,NLMAX-1

  !N = KNVT(ILEV+1) + 3*KNET(ILEV+1) + 9*KNAT(ILEV+1) + 27*KNEL(ILEV+1)
  !N = KNVT(ILEV+1) + 3*KNET(ILEV+1) + 9*KNAT(ILEV+1) + 27*KNEL(ILEV+1)
  N =    mg_mesh%level(ilev+1)%nvt+&
      3* mg_mesh%level(ilev+1)%net+&
      9* mg_mesh%level(ilev+1)%nat+&
      27*mg_mesh%level(ilev+1)%nel

  ALLOCATE(mg_E013Prol(ILEV)%a(N))
  ALLOCATE(mg_E013Rest(ILEV)%a(N))
  ALLOCATE(mg_E013ProlM(ILEV)%ColA(N))
  ALLOCATE(mg_E013RestM(ILEV)%ColA(N))
  mg_E013ProlM(ILEV)%na = N
  mg_E013RestM(ILEV)%na = N

! NDOF = KNVT(ILEV+1)+KNET(ILEV+1)+KNAT(ILEV+1)+KNEL(ILEV+1)

  NDOF = mg_mesh%level(ilev+1)%nvt+&
         mg_mesh%level(ilev+1)%net+&
         mg_mesh%level(ilev+1)%nat+&
         mg_mesh%level(ilev+1)%nel

  mg_E013ProlM(ILEV)%nu = NDOF
  ALLOCATE(mg_E013ProlM(ILEV)%LdA(NDOF+1))

  !NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)
  NDOF = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

  mg_E013RestM(ILEV)%nu = NDOF
  ALLOCATE(mg_E013RestM(ILEV)%LdA(NDOF+1))


 END DO

 ! Other Multilevel Vectors
 ALLOCATE(myScalar%def(NLMIN:NLMAX))
 ALLOCATE(myScalar%rhs(NLMIN:NLMAX))
 ALLOCATE(myScalar%aux(NLMIN:NLMAX))
 ALLOCATE(myScalar%sol(NLMIN:NLMAX))

 ! For CC
 ALLOCATE(myScalar%dsol(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX

  !NDOF = KNVT(ILEV)+KNET(ILEV)+KNAT(ILEV)+KNEL(ILEV)
  NDOF = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

  ALLOCATE(myScalar%def(ILEV)%x(3*NDOF))
  ALLOCATE(myScalar%rhs(ILEV)%x(3*NDOF))
  ALLOCATE(myScalar%aux(ILEV)%x(3*NDOF))
  ALLOCATE(myScalar%sol(ILEV)%x(3*NDOF))

  ! For CC
  ALLOCATE(myScalar%dsol(ILEV)%x(3*NDOF))
 END DO

 ILEV=NLMAX

END SUBROUTINE InitializeQuadScalar
!
! ----------------------------------------------
!
SUBROUTINE InitializeLinScalar(myScalar)
TYPE(TLinScalar) myScalar
INTEGER NDOF,NDOF2,N

 ALLOCATE(myScalar%knprP(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX
  NDOF = 4*mg_mesh%level(ilev)%nel
  ALLOCATE(myScalar%knprP(ILEV)%x(NDOF))
  myScalar%knprP(ilev)%x = 0
 END DO

 ILEV=NLMAX

 NDOF = 4*mg_mesh%level(ilev)%nel

 NDOF2 = mg_mesh%level(ilev)%nvt+&
         mg_mesh%level(ilev)%net+&
         mg_mesh%level(ilev)%nat+&
         mg_mesh%level(ilev)%nel

 NVT  = mg_mesh%level(ilev)%nvt

 myScalar%ndof = NDOF
 IF (ALLOCATED(mg_lMat)) myScalar%na = mg_lMat(ILEV)%na

 ALLOCATE(myScalar%valP_old(myScalar%ndof))
 ! For CC
 ALLOCATE(myScalar%P_old(myScalar%ndof))
 ALLOCATE(myScalar%P_new(myScalar%ndof))

 ALLOCATE(myScalar%ST_P(myScalar%ndof))
 ALLOCATE(myScalar%Q2(NDOF2))
 ALLOCATE(myScalar%valP_GMV(myScalar%ndof*8))

 ALLOCATE(myScalar%valP(NLMIN:NLMAX))

 ! For CC
 ALLOCATE(myScalar%dvalP(NLMIN:NLMAX))

 ALLOCATE(myScalar%rhsP(NLMIN:NLMAX))
 ALLOCATE(myScalar%defP(NLMIN:NLMAX))
 ALLOCATE(myScalar%auxP(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX
  NDOF = 4*mg_mesh%level(ilev)%nel
  ALLOCATE(myScalar%valP(ILEV)%x(NDOF))

  ! For CC
  ALLOCATE(myScalar%dvalP(ILEV)%x(NDOF))

  ALLOCATE(myScalar%defP(ILEV)%x(NDOF))
  ALLOCATE(myScalar%auxP(ILEV)%x(NDOF))
  ALLOCATE(myScalar%rhsP(ILEV)%x(NDOF))
 END DO

 ALLOCATE(mg_E012Prol(NLMIN:NLMAX-1))
 DO ILEV=NLMIN,NLMAX-1
  NDOF = 4*mg_mesh%level(ilev)%nel
  N = NDOF*32
  ALLOCATE(mg_E012Prol(ILEV)%a(N))
 END DO

 ILEV=NLMAX

END SUBROUTINE InitializeLinScalar
!
! ----------------------------------------------
!
SUBROUTINE InitializeProlRest(qScalar,lScalar)
TYPE(TLinScalar), INTENT(INOUT), TARGET :: lScalar
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar

 MyMG%MinLev             = qScalar%prm%MGprmIn%MinLev
 MyMG%MedLev             = qScalar%prm%MGprmIn%MedLev
 myMG%MaxLev=NLMAX
!  myMG%MinLev=NLMIN
!  myMG%MaxLev=NLMAX
 IF(myid.eq.showid) WRITE(*,*) "Initialization of velocity prolongation matrix"
 myMG%bProlRest => lScalar%bProlRest
 MyMG%cVariable = "Velocity"
 CALL mgProlRestInit()

 MyMG%MinLev = lScalar%prm%MGprmIn%MinLev
 MyMG%MedLev = lScalar%prm%MGprmIn%MedLev
 myMG%MaxLev = NLMAX

 IF(myid.eq.showid) WRITE(*,*) "Initialization of pressure prolongation matrix"
 myMG%bProlRest => qScalar%bProlRest
 MyMG%cVariable = "Pressure"
 CALL mgProlRestInit()

END SUBROUTINE InitializeProlRest
!
! ----------------------------------------------
!
SUBROUTINE Matdef_General_LinScalar(lScalar,qScalar,lPScalar,idef)
INTEGER :: i,j,jj,idef
TYPE(TLinScalar) lScalar
TYPE(TParLinScalar) lPScalar
TYPE(TQuadScalar) qScalar
REAL*8 def1,def2

IF (idef.eq.1) THEN
!  qScalar%valU(NLMAX)%x = 1d0
!  qScalar%valV(NLMAX)%x = 1d0
!  qScalar%valW(NLMAX)%x = 1d0
 CALL BT_Mul_U_mod(lqMat%ColA,lqMat%LdA,BTXMat,BTYMat,BTZMat,&
      qScalar%valU,qScalar%valV,qScalar%valW,&
      lScalar%defP(NLMAX)%x,qScalar%ndof,lScalar%ndof,TSTEP)
!  CALL LL21(lScalar%defP(NLMAX)%x,lScalar%ndof,def2)
!  WRITE(*,*) "Before",myid,def2
 ELSE

  CALL GetParPressure(lScalar%valP(NLMAX)%x,lPScalar%Val)
  CALL C_Mul_q (CMat,lMat%LdA,lMat%ColA,&
       CPMat,lPMat%LdA,lPMat%ColA,&
       lScalar%ValP(NLMAX)%x,lPScalar%Val,lScalar%defP(NLMAX)%x,lMat%nu)
!  CALL LL21(lScalar%defP(NLMAX)%x,lScalar%ndof,def2)
!  WRITE(*,*) "After",myid,def2
END IF

END SUBROUTINE Matdef_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_General_LinScalar(myScalar,&
           resScalarP,defScalar,rhsScalar)
TYPE(TLinScalar), INTENT(INOUT) :: myScalar
REAL*8  resScalarP,defScalar,rhsScalar,RESF,RESU
REAL*8  RESF_P,RESU_P

CALL LL21 (myScalar%rhsP(NLMAX)%x,myScalar%ndof,RESF_P)
RESF=MAX(1D-16,RESF_P)

CALL LL21 (myScalar%defP(NLMAX)%x,myScalar%ndof,RESU_P)
RESU=MAX(1D-16,RESU_P)

resScalarP = RESU_P/RESF
defScalar = RESU
rhsScalar = RESF

END SUBROUTINE Resdfk_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Solve_General_LinScalar(lScalar,lPScalar,qScalar,Bndry_Mat,Bndry_Def,mfile)
INTEGER mfile
TYPE(TLinScalar), INTENT(INOUT), TARGET :: lScalar
TYPE(TParLinScalar), INTENT(INOUT), TARGET ::  lPScalar
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar
REAL*8 resP,resDP,res,ddd
REAL*8 dnormu1,dnormu2,dnormu3,dnormu
INTEGER :: I,J
EXTERNAL Bndry_Mat,Bndry_Def

! REAL*8, ALLOCATABLE, DIMENSION(:) :: d1,D2,d3

 CALL ZTIME(myStat%t0)

 CALL LL21 (qScalar%ValU,qScalar%ndof,dnormu1)
 CALL LL21 (qScalar%ValV,qScalar%ndof,dnormu2)
 CALL LL21 (qScalar%ValW,qScalar%ndof,dnormu3)
 dnormu=sqrt(dnormu1*dnormu1+dnormu2*dnormu2+dnormu3*dnormu3)
 IF (ABS(dnormu).LT.1D-8) dnormu=1D-8
 CALL COMM_Maximum(dnormu)

  DO ILEV = NLMIN,NLMAX
   CALL SETLEV(2)
   CALL Bndry_Mat(mg_CMat(ILEV)%a,mg_lMat(ILEV)%LdA,lScalar%knprP(ILEV)%x,nel,1)
   if (myid.ne.0) CALL Bndry_Mat(mg_CPMat(ILEV)%a,mg_lPMat(ILEV)%LdA,lScalar%knprP(ILEV)%x,nel,0)
  END DO
  if (myid.ne.0) then 
   ILEV = NLMAX
   CALL SETLEV(2)
   CALL Bndry_Def(lScalar%rhsP(ILEV)%x,lScalar%knprP(ILEV)%x,nel)
  end if
 
 ! Set up the pointers for the Pressure MG driver 
 MyMG%A   => mg_CMat
 MyMG%AP  => mg_CPMat
 MyMG%L   => mg_lMat
 MyMG%LP  => mg_lPMat
 MyMG%X   => lScalar%valP
 MyMG%D   => lScalar%defP
 MyMG%B   => lScalar%rhsP
 MyMG%AUX => lScalar%auxP
 MyMG%XP  => lPScalar%Val
 MyMG%KNPRP => lScalar%knprP
 
 ! Set up the parameters for the Pressure MG driver 
 MyMG%cVariable     = "Pressure"
 MyMG%MinLev        = lScalar%prm%MGprmIn%MinLev
 MyMG%MedLev        = lScalar%prm%MGprmIn%MedLev
 MyMG%CrsSolverType = lScalar%prm%MGprmIn%CrsSolverType
 MyMG%CrsRelaxPrm   = lScalar%prm%MGprmIn%CrsRelaxPrm
 
 dnormu1 = DBLE(NLMAX)
 CALL COMM_Maximum(dnormu1)
 MyMG%MaxLev = NINT(dnormu1)
 MyMG%bProlRest => lScalar%bProlRest
 MyMG%MinIterCycle = lScalar%prm%MGprmIn%MinIterCycle
 MyMG%MaxIterCycle = lScalar%prm%MGprmIn%MaxIterCycle
 MyMG%nIterCoarse = lScalar%prm%MGprmIn%nIterCoarse
 MyMG%DefImprCoarse = lScalar%prm%MGprmIn%DefImprCoarse
 MyMG%nSmootherSteps = lScalar%prm%MGprmIn%nSmootherSteps
 MyMG%RLX            = lScalar%prm%MGprmIn%RLX
 MyMG%CycleType = lScalar%prm%MGprmIn%CycleType
 MyMG%Criterion1 = lScalar%prm%MGprmIn%Criterion1
 MyMG%Criterion2 = dnormu*lScalar%prm%MGprmIn%Criterion2/TSTEP

!-------------------  P - Component -------------------!
 CALL MG_Solver(mfile,mterm)
 lScalar%prm%MGprmOut%UsedIterCycle = myMG%UsedIterCycle
 lScalar%prm%MGprmOut%nIterCoarse   = myMG%nIterCoarse
 lScalar%prm%MGprmOut%DefInitial    = myMG%DefInitial
 lScalar%prm%MGprmOut%DefFinal      = myMG%DefFinal
 lScalar%prm%MGprmOut%RhoMG1        = myMG%RhoMG1
 lScalar%prm%MGprmOut%RhoMG2        = myMG%RhoMG2

 CALL ZTIME(myStat%t1)
 myStat%tMGP = myStat%tMGP + (myStat%t1-myStat%t0)
 myStat%iLinP = myStat%iLinP + lScalar%prm%MGprmOut%UsedIterCycle
!  IF (myid.eq.showid) write(*,*) myStat%t1-myStat%t0, "time needed ..."

END SUBROUTINE Solve_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Matdef_general_QuadScalar(myScalar,idef)
EXTERNAL E013
INTEGER :: idef
INTEGER I,J
TYPE(TQuadScalar) myScalar
REAL*8 daux
REAL tttx1,tttx0
! the formula for crank-nicholson may be wrong

! Build up the matrix
 IF (idef.eq.-1) THEN
  DO ILEV=NLMIN,NLMAX

    !!-------------------    POINTER Setup  -------------------!!
    IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
     A11Mat     => mg_A11Mat(ILEV)%a
     A22Mat     => mg_A22Mat(ILEV)%a
     A33Mat     => mg_A33Mat(ILEV)%a
     A12Mat     => mg_A12Mat(ILEV)%a
     A13Mat     => mg_A13Mat(ILEV)%a
     A23Mat     => mg_A23Mat(ILEV)%a
     A21Mat     => mg_A21Mat(ILEV)%a
     A31Mat     => mg_A31Mat(ILEV)%a
     A32Mat     => mg_A32Mat(ILEV)%a
    ELSE
     A11Mat     => mg_A11Mat(ILEV)%a
     A22Mat     => mg_A22Mat(ILEV)%a
     A33Mat     => mg_A33Mat(ILEV)%a
    END IF

    IF (myMatrixRenewal%S.GE.1) THEN
     S11Mat   => mg_S11Mat(ILEV)%a
     S22Mat   => mg_S22Mat(ILEV)%a
     S33Mat   => mg_S33Mat(ILEV)%a
     S12Mat   => mg_S12Mat(ILEV)%a
     S13Mat   => mg_S13Mat(ILEV)%a
     S23Mat   => mg_S23Mat(ILEV)%a
     S21Mat   => mg_S21Mat(ILEV)%a
     S31Mat   => mg_S31Mat(ILEV)%a
     S32Mat   => mg_S32Mat(ILEV)%a
    END IF

    IF (myMatrixRenewal%D.GE.1) THEN
     DMat     => mg_DMat(ILEV)%a
    END IF

    IF (myMatrixRenewal%K.GE.1) THEN
     KMat     => mg_KMat(ILEV)%a
    END IF

    IF (myMatrixRenewal%M.GE.1) THEN
     MMat     => mg_MMat(ILEV)%a
     MlRhoMat => mg_MlRhoMat(ILEV)%a
    END IF

    qMat     => mg_qMat(ILEV)
    !!-------------------    POINTER Setup  -------------------!!

    !!-------------------  MATRIX Assembly -------------------!!
    IF (bNonNewtonian) THEN
     IF (myMatrixRenewal%S.EQ.0) THEN!K-Alternative fehlt
      IF (myMatrixRenewal%K.GE.1) THEN
       DO I=1,qMat%nu
        J = qMat%LdA(I)
        daux = MlRhoMat(I) + thstep*(2d0*DMat(J)+KMat(J))
        A11mat(J) = daux
        A22mat(J) = daux
        A33mat(J) = daux
        DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
         daux = thstep*(2d0*DMat(J)+KMat(J))! 
         A11mat(J) =  daux
         A22mat(J) =  daux
         A33mat(J) =  daux
        END DO
       END DO
      ELSE
       DO I=1,qMat%nu
        J = qMat%LdA(I)
        daux = MlRhoMat(I) + thstep*(2d0*DMat(J))!  
        A11mat(J) = daux
        A22mat(J) = daux
        A33mat(J) = daux
        DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
         daux =  thstep*(2d0*DMat(J))!
         A11mat(J) =  daux
         A22mat(J) =  daux
         A33mat(J) =  daux
        END DO
       END DO 
      END IF     
     ELSE
      IF (myMatrixRenewal%K.GE.1) THEN
       DO I=1,qMat%nu
        J = qMat%LdA(I)
        daux = KMat(J)
        A11Mat(J) = MlRhoMat(I) + thstep*(S11Mat(J) + daux)
        A22Mat(J) = MlRhoMat(I) + thstep*(S22Mat(J) + daux)
        A33Mat(J) = MlRhoMat(I) + thstep*(S33Mat(J) + daux)
        A12Mat(J) = thstep*S12Mat(J)
        A13Mat(J) = thstep*S13Mat(J)
        A23Mat(J) = thstep*S23Mat(J)
        A21Mat(J) = thstep*S21Mat(J)
        A31Mat(J) = thstep*S31Mat(J)
        A32Mat(J) = thstep*S32Mat(J)
        DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
         daux = KMat(J)
         A11mat(J) = thstep*(S11Mat(J) + daux)
         A22mat(J) = thstep*(S22Mat(J) + daux)
         A33mat(J) = thstep*(S33Mat(J) + daux)
         A12Mat(J) = thstep*S12Mat(J)
         A13Mat(J) = thstep*S13Mat(J)
         A23Mat(J) = thstep*S23Mat(J)
         A21Mat(J) = thstep*S21Mat(J)
         A31Mat(J) = thstep*S31Mat(J)
         A32Mat(J) = thstep*S32Mat(J)
        END DO
       END DO

      ELSE

       DO I=1,qMat%nu
        J = qMat%LdA(I)
        A11Mat(J) = MlRhoMat(I) + thstep*S11Mat(J)
        A22Mat(J) = MlRhoMat(I) + thstep*S22Mat(J)
        A33Mat(J) = MlRhoMat(I) + thstep*S33Mat(J)
        A12Mat(J) = thstep*S12Mat(J)
        A13Mat(J) = thstep*S13Mat(J)
        A23Mat(J) = thstep*S23Mat(J)
        A21Mat(J) = thstep*S21Mat(J)
        A31Mat(J) = thstep*S31Mat(J)
        A32Mat(J) = thstep*S32Mat(J)
        DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
         A11mat(J) = thstep*S11Mat(J)
         A22mat(J) = thstep*S22Mat(J)
         A33mat(J) = thstep*S33Mat(J)
         A12Mat(J) = thstep*S12Mat(J)
         A13Mat(J) = thstep*S13Mat(J)
         A23Mat(J) = thstep*S23Mat(J)
         A21Mat(J) = thstep*S21Mat(J)
         A31Mat(J) = thstep*S31Mat(J)
         A32Mat(J) = thstep*S32Mat(J)
        END DO
       END DO
      END IF
     END IF
    ELSE

     IF (myMatrixRenewal%K.GE.1) THEN

      DO I=1,qMat%nu
       J = qMat%LdA(I)
       daux = MlRhoMat(I) + thstep*(DMat(J)+KMat(J))
       A11mat(J) = daux
       A22mat(J) = daux
       A33mat(J) = daux
       DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
        daux =  +thstep*(DMat(J)+KMat(J))
        A11mat(J) =  daux
        A22mat(J) =  daux
        A33mat(J) =  daux
       END DO
      END DO

     ELSE

      DO I=1,qMat%nu
       J = qMat%LdA(I)
       daux = MlRhoMat(I) + thstep*(DMat(J))
       A11mat(J) = daux
       A22mat(J) = daux
       A33mat(J) = daux
       DO J=qMat%LdA(I)+1,qMat%LdA(I+1)-1
        daux =  +thstep*(DMat(J))
        A11mat(J) =  daux
        A22mat(J) =  daux
        A33mat(J) =  daux
       END DO
      END DO

     END IF

    END IF

  END DO
 END IF
    !!-------------------  MATRIX Assembly -------------------!!

    !!-------------------    POINTER Setup  -------------------!!
 ILEV=NLMAX
 CALL SETLEV(2)

 IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
  A12Mat     => mg_A12Mat(ILEV)%a
  A13Mat     => mg_A13Mat(ILEV)%a
  A23Mat     => mg_A23Mat(ILEV)%a
  A21Mat     => mg_A21Mat(ILEV)%a
  A31Mat     => mg_A31Mat(ILEV)%a
  A32Mat     => mg_A32Mat(ILEV)%a
 ELSE
  A11Mat     => mg_A11Mat(ILEV)%a
  A22Mat     => mg_A22Mat(ILEV)%a
  A33Mat     => mg_A33Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%S.GE.1) THEN
  S11Mat   => mg_S11Mat(ILEV)%a
  S22Mat   => mg_S22Mat(ILEV)%a
  S33Mat   => mg_S33Mat(ILEV)%a
  S12Mat   => mg_S12Mat(ILEV)%a
  S13Mat   => mg_S13Mat(ILEV)%a
  S23Mat   => mg_S23Mat(ILEV)%a
  S21Mat   => mg_S21Mat(ILEV)%a
  S31Mat   => mg_S31Mat(ILEV)%a
  S32Mat   => mg_S32Mat(ILEV)%a
 END IF

 IF (myMatrixRenewal%D.GE.1) THEN
  DMat     => mg_DMat(ILEV)%a
 END IF

 IF (myMatrixRenewal%K.GE.1) THEN
  KMat     => mg_KMat(ILEV)%a
 END IF

 IF (myMatrixRenewal%M.GE.1) THEN
  MMat     => mg_MMat(ILEV)%a
  MlRhoMat => mg_MlRhoMat(ILEV)%a
 END IF

 qMat     => mg_qMat(ILEV)
    !!-------------------    POINTER Setup  -------------------!!

 ! Build up the defect
 IF (idef.eq. 1) THEN
   myScalar%defU = 0d0
   myScalar%defV = 0d0
   myScalar%defW = 0d0

   IF (myMatrixRenewal%M.GE.1) THEN
    myScalar%defU = myScalar%defU + MlRhoMat*myScalar%valU
    myScalar%defV = myScalar%defV + MlRhoMat*myScalar%valV
    myScalar%defW = myScalar%defW + MlRhoMat*myScalar%valW
   END IF

   IF (myMatrixRenewal%D.GE.1.AND.(.NOT.bNonNewtonian)) THEN
    CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defU,-thstep,1d0)
    CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defV,-thstep,1d0)
    CALL LAX17(DMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defW,-thstep,1d0)
   END IF

   IF (myMatrixRenewal%K.GE.1) THEN
    CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defU,-thstep,1d0)
    CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defV,-thstep,1d0)
    CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defW,-thstep,1d0)
   END IF

   IF (myMatrixRenewal%S.GE.1) THEN
    CALL LAX17(S11Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defU,-thstep,1d0)
    CALL LAX17(S12Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defU,-thstep,1d0)
    CALL LAX17(S13Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defU,-thstep,1d0)

    CALL LAX17(S21Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defV,-thstep,1d0)
    CALL LAX17(S22Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defV,-thstep,1d0)
    CALL LAX17(S23Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defV,-thstep,1d0)

    CALL LAX17(S31Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defW,-thstep,1d0)
    CALL LAX17(S32Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defW,-thstep,1d0)
    CALL LAX17(S33Mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defW,-thstep,1d0)
   END IF

   IF(bNonNewtonian.AND.myMatrixRenewal%S.EQ.0) THEN
    CALL ZTIME(tttx0)
    
    CALL STRESS(myScalar%valU,myScalar%valV,myScalar%valW,&
    Temperature,&
    myScalar%defU, myScalar%defV, myScalar%defW,&
    Viscosity,&
    mg_mesh%level(ILEV)%kvert,&
    mg_mesh%level(ILEV)%karea,&
    mg_mesh%level(ILEV)%kedge,&
    mg_mesh%level(ILEV)%dcorvg,&
    E013 ) ! S*u

    CALL ZTIME(tttx1)
    myStat%tSMat = myStat%tSMat + (tttx1-tttx0)
   END IF

 ELSE

   IF (bNonNewtonian) THEN

    IF(myMatrixRenewal%S.GE.1) THEN
     CALL LAX17(A11mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defU,-1d0,1d0)
     CALL LAX17(A12Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defU,-1d0,1d0)
     CALL LAX17(A13Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defU,-1d0,1d0)

     CALL LAX17(A21Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defV,-1d0,1d0)
     CALL LAX17(A22mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defV,-1d0,1d0)
     CALL LAX17(A23Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defV,-1d0,1d0)

     CALL LAX17(A31Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valU,myScalar%defW,-1d0,1d0)
     CALL LAX17(A32Mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valV,myScalar%defW,-1d0,1d0)
     CALL LAX17(A33mat,qMat%ColA,qMat%LdA,qMat%nu,&
     myScalar%valW,myScalar%defW,-1d0,1d0)

    ELSE
     IF (myMatrixRenewal%M.GE.1) THEN
      myScalar%defU = myScalar%defU - MlRhoMat*myScalar%valU
      myScalar%defV = myScalar%defV - MlRhoMat*myScalar%valV
      myScalar%defW = myScalar%defW - MlRhoMat*myScalar%valW
     END IF

     IF (myMatrixRenewal%K.GE.1) THEN
      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
      myScalar%valU,myScalar%defU,-thstep,1d0)
      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
      myScalar%valV,myScalar%defV,-thstep,1d0)
      CALL LAX17(KMat,qMat%ColA,qMat%LdA,qMat%nu,&
      myScalar%valW,myScalar%defW,-thstep,1d0)
     END IF

     CALL ZTIME(tttx0)

     CALL STRESS(myScalar%valU,myScalar%valV,myScalar%valW,&
     Temperature,&
     myScalar%defU, myScalar%defV, myScalar%defW,&
     Viscosity,&
     mg_mesh%level(ILEV)%kvert,&
     mg_mesh%level(ILEV)%karea,&
     mg_mesh%level(ILEV)%kedge,&
     mg_mesh%level(ILEV)%dcorvg,&
     E013 ) ! S*u

     CALL ZTIME(tttx1)
     myStat%tSMat = myStat%tSMat + (tttx1-tttx0)
    END IF   

   ELSE

    CALL LAX17(A11mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valU,myScalar%defU,-1d0,1d0)
    CALL LAX17(A22mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valV,myScalar%defV,-1d0,1d0)
    CALL LAX17(A33mat,qMat%ColA,qMat%LdA,qMat%nu,&
    myScalar%valW,myScalar%defW,-1d0,1d0)

   END IF

 END IF

END SUBROUTINE Matdef_general_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_General_QuadScalar(myScalar,&
           resScalarU,resScalarV,resScalarW,defScalar,rhsScalar)
TYPE(TQuadScalar), INTENT(INOUT) :: myScalar
REAL*8  resScalarU,resScalarV,resScalarW,defScalar,rhsScalar,RESF,RESU
REAL*8  RESF_U, RESF_V, RESF_W, RESU_U, RESU_V, RESU_W

CALL LL21 (myScalar%rhsU,myScalar%ndof,RESF_U)
CALL LL21 (myScalar%rhsV,myScalar%ndof,RESF_V)
CALL LL21 (myScalar%rhsW,myScalar%ndof,RESF_W)
RESF=MAX(1D-16,RESF_U,RESF_V,RESF_W)

CALL LL21 (myScalar%auxU,myScalar%ndof,RESU_U)
CALL LL21 (myScalar%auxV,myScalar%ndof,RESU_V)
CALL LL21 (myScalar%auxW,myScalar%ndof,RESU_W)
RESU=MAX(1D-16,RESU_U,RESU_V,RESU_W)

!WRITE(*,'(I1,6D12.4)') myid,RESU_U,RESU_V,RESU_W,RESF
resScalarU = RESU_U/RESF
resScalarV = RESU_V/RESF
resScalarW = RESU_W/RESF
defScalar = RESU
rhsScalar = RESF
! write(*,*) RESU,RESF

END SUBROUTINE Resdfk_General_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Solve_General_QuadScalar(myScalar,Bndry_Val,Bndry_Mat,Bndry_Mat_9,mfile)
INTEGER mfile
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: myScalar
REAL*8 daux,nrm_U,nrm_V,nrm_W
INTEGER ndof
EXTERNAL Bndry_Val,Bndry_Mat,Bndry_Mat_9

nrm_U = 0d0
nrm_V = 0d0
nrm_W = 0d0

 CALL ZTIME(myStat%t0)

 IF (myid.ne.0) THEN
  DO ILEV = NLMIN,NLMAX
   IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
    CALL Bndry_Mat_9(mg_A11mat(ILEV)%a,mg_A22mat(ILEV)%a,mg_A33mat(ILEV)%a,&
         mg_A12mat(ILEV)%a,mg_A13mat(ILEV)%a,mg_A23mat(ILEV)%a,&
         mg_A21mat(ILEV)%a,mg_A31mat(ILEV)%a,mg_A32mat(ILEV)%a,&
         mg_qMat(ILEV)%LdA,myScalar%knprU(NLMAX)%x,myScalar%knprV(NLMAX)%x,&
         myScalar%knprW(NLMAX)%x,mg_qMat(ILEV)%nu)
    CALL E013UVWMAT(mg_A11mat(ILEV)%a,mg_A22mat(ILEV)%a,mg_A33mat(ILEV)%a,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%nu)
   ELSE
    CALL Bndry_Mat(mg_A11mat(ILEV)%a,mg_A22mat(ILEV)%a,mg_A33mat(ILEV)%a,&
         mg_qMat(ILEV)%LdA,myScalar%knprU(NLMAX)%x,myScalar%knprV(NLMAX)%x,&
         myScalar%knprW(NLMAX)%x,mg_qMat(ILEV)%nu)
    CALL E013UVWMAT(mg_A11mat(ILEV)%a,mg_A22mat(ILEV)%a,mg_A33mat(ILEV)%a,mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%nu)
!    CALL E013Mat(mg_A11mat(ILEV)%a,mg_A22mat(ILEV)%a,mg_A33mat(ILEV)%a,&
!         mg_qMat(ILEV)%LdA,mg_qMat(ILEV)%nu)
   END IF
  END DO
  CALL LL21 (myScalar%defU,myScalar%ndof,nrm_U)
  CALL LL21 (myScalar%defV,myScalar%ndof,nrm_V)
  CALL LL21 (myScalar%defW,myScalar%ndof,nrm_W)
  CALL LCL1 (myScalar%valU,myScalar%ndof)
  CALL LCL1 (myScalar%valV,myScalar%ndof)
  CALL LCL1 (myScalar%valW,myScalar%ndof)
 END IF

 daux = MAX(nrm_U,nrm_V,nrm_W)
 CALL COMM_Maximum(daux)

!--------------- Set up the MG driver -----------------!
 IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
  MyMG%A11  => mg_A11Mat
  MyMG%A22  => mg_A22Mat
  MyMG%A33  => mg_A33Mat
  MyMG%A12  => mg_A12Mat
  MyMG%A13  => mg_A13Mat
  MyMG%A23  => mg_A23Mat
  MyMG%A21  => mg_A21Mat
  MyMG%A31  => mg_A31Mat
  MyMG%A32  => mg_A32Mat
 ELSE
  MyMG%A11    => mg_A11Mat
  MyMG%A22    => mg_A22Mat
  MyMG%A33    => mg_A33Mat
 END IF
 MyMG%L    => mg_qMat
 MyMG%D    => myScalar%def
 MyMG%AUX  => myScalar%aux
 MyMG%KNPRU => myScalar%knprU(NLMAX)%x
 MyMG%KNPRV => myScalar%knprV(NLMAX)%x
 MyMG%KNPRW => myScalar%knprW(NLMAX)%x
 MyMG%cVariable = "Velocity"
 MyMG%bProlRest => myScalar%bProlRest
! Variable specific settings 
 MyMG%MinIterCycle       = myScalar%prm%MGprmIn%MinIterCycle
 MyMG%MaxIterCycle       = myScalar%prm%MGprmIn%MaxIterCycle
 MyMG%nIterCoarse        = myScalar%prm%MGprmIn%nIterCoarse
 MyMG%DefImprCoarse      = myScalar%prm%MGprmIn%DefImprCoarse
 MyMG%nSmootherSteps     = myScalar%prm%MGprmIn%nSmootherSteps
 MyMG%CycleType          = myScalar%prm%MGprmIn%CycleType
 MyMG%Criterion1         = myScalar%prm%MGprmIn%Criterion1
 MyMG%Criterion2         = myScalar%prm%MGprmIn%Criterion2*daux
 MyMG%RLX                = myScalar%prm%MGprmIn%RLX
 MyMG%MinLev             = myScalar%prm%MGprmIn%MinLev
 MyMG%MedLev             = myScalar%prm%MGprmIn%MedLev
 MyMG%CrsSolverType      = myScalar%prm%MGprmIn%CrsSolverType
 MyMG%CrsRelaxPrm        = myScalar%prm%MGprmIn%CrsRelaxPrm
 MyMG%CrsRelaxParPrm     = myScalar%prm%MGprmIn%CrsRelaxParPrm

 daux = DBLE(NLMAX)
 CALL COMM_Maximum(daux)
 MyMG%MaxLev             = NINT(daux)

!-------------------  U - Component -------------------!
 ndof = SIZE(myScalar%ValU)

 myScalar%sol(NLMAX)%x(0*ndof+1:1*ndof) = myScalar%ValU
 myScalar%sol(NLMAX)%x(1*ndof+1:2*ndof) = myScalar%ValV
 myScalar%sol(NLMAX)%x(2*ndof+1:3*ndof) = myScalar%ValW
 MyMG%X    => myScalar%sol

 myScalar%rhs(NLMAX)%x(0*ndof+1:1*ndof) = myScalar%defU
 myScalar%rhs(NLMAX)%x(1*ndof+1:2*ndof) = myScalar%defV
 myScalar%rhs(NLMAX)%x(2*ndof+1:3*ndof) = myScalar%defW
 MyMG%B    => myScalar%rhs

 CALL MG_Solver(mfile,mterm)

 myScalar%ValU = myScalar%sol(NLMAX)%x(0*ndof+1:1*ndof)
 myScalar%ValV = myScalar%sol(NLMAX)%x(1*ndof+1:2*ndof)
 myScalar%ValW = myScalar%sol(NLMAX)%x(2*ndof+1:3*ndof)

 myScalar%prm%MGprmOut(1)%UsedIterCycle = myMG%UsedIterCycle
 myScalar%prm%MGprmOut(1)%nIterCoarse   = myMG%nIterCoarse
 myScalar%prm%MGprmOut(1)%DefInitial    = myMG%DefInitial
 myScalar%prm%MGprmOut(1)%DefFinal      = myMG%DefFinal
 myScalar%prm%MGprmOut(1)%RhoMG1        = myMG%RhoMG1
 myScalar%prm%MGprmOut(1)%RhoMG2        = myMG%RhoMG2
! !-------------------  V - Component -------------------!
!  myScalar%sol(NLMAX)%x = myScalar%ValV
!  MyMG%X    => myScalar%sol
!  myScalar%rhs(NLMAX)%x = myScalar%defV
!  MyMG%B    => myScalar%rhs
!  CALL MG_Solver(mfile,mterm)
!  myScalar%ValV = myScalar%sol(NLMAX)%x
!  myScalar%prm%MGprmOut(2)%UsedIterCycle = myMG%UsedIterCycle
!  myScalar%prm%MGprmOut(2)%nIterCoarse   = myMG%nIterCoarse
!  myScalar%prm%MGprmOut(2)%DefInitial    = myMG%DefInitial
!  myScalar%prm%MGprmOut(2)%DefFinal      = myMG%DefFinal
!  myScalar%prm%MGprmOut(2)%RhoMG1        = myMG%RhoMG1
!  myScalar%prm%MGprmOut(2)%RhoMG2        = myMG%RhoMG2
!  !-------------------  W - Component -------------------!
!  myScalar%sol(NLMAX)%x = myScalar%ValW
!  MyMG%X    => myScalar%sol
!  myScalar%rhs(NLMAX)%x = myScalar%defW
!  MyMG%B    => myScalar%rhs
!  CALL MG_Solver(mfile,mterm)
!  myScalar%ValW = myScalar%sol(NLMAX)%x
!  myScalar%prm%MGprmOut(3)%UsedIterCycle = myMG%UsedIterCycle
!  myScalar%prm%MGprmOut(3)%nIterCoarse   = myMG%nIterCoarse
!  myScalar%prm%MGprmOut(3)%DefInitial    = myMG%DefInitial
!  myScalar%prm%MGprmOut(3)%DefFinal      = myMG%DefFinal
!  myScalar%prm%MGprmOut(3)%RhoMG1        = myMG%RhoMG1
!  myScalar%prm%MGprmOut(3)%RhoMG2        = myMG%RhoMG2
!------------------------------------------------------!

 ! Update the solution
 IF (myid.ne.0) THEN
  CALL LLC1(myScalar%valU_old,myScalar%valU,&
       myScalar%ndof,1D0,1D0)
  CALL LLC1(myScalar%valV_old,myScalar%valV,&
       myScalar%ndof,1D0,1D0)
  CALL LLC1(myScalar%valW_old,myScalar%valW,&
       myScalar%ndof,1D0,1D0)

  ! Set dirichlet boundary conditions on the solution
  CALL Bndry_Val()
 END IF

 CALL ZTIME(myStat%t1)
 myStat%tMGUVW = myStat%tMGUVW + (myStat%t1-myStat%t0)
 myStat%iLinUVW = myStat%iLinUVW + myScalar%prm%MGprmOut(1)%UsedIterCycle &
                                 + myScalar%prm%MGprmOut(2)%UsedIterCycle &
                                 + myScalar%prm%MGprmOut(3)%UsedIterCycle

END SUBROUTINE Solve_General_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE ProlongateSolutionSub(qScalar,lScalar,Bndry_Val)
TYPE(TLinScalar), INTENT(INOUT), TARGET :: lScalar
TYPE(TQuadScalar), INTENT(INOUT), TARGET :: qScalar
INTEGER J,K,ndof
REAL*8 dnorm
EXTERNAL Bndry_Val

IF (myid.ne.0) THEN

 ILEV=NLMAX
 CALL SETLEV(2)

 J = NLMAX
 K = NLMAX-1
 mgLev = J

 IF(myid.eq.showid) WRITE(*,*) "Prolongation of velocity solution to a higher level"

 MyMG%cVariable = "Velocity"
 MyMG%KNPRU => qScalar%knprU(J)%x
 MyMG%KNPRV => qScalar%knprV(J)%x
 MyMG%KNPRW => qScalar%knprW(J)%x

 ndof = KNVT(K) + KNAT(K) + KNET(K) + KNEL(K)
 qScalar%sol(K)%x(0*ndof+1:1*ndof) = qScalar%ValU(1:ndof)
 qScalar%sol(K)%x(1*ndof+1:2*ndof) = qScalar%ValV(1:ndof)
 qScalar%sol(K)%x(2*ndof+1:3*ndof) = qScalar%ValW(1:ndof)
 MyMG%X    => qScalar%sol
 MyMG%AUX  => qScalar%aux

 CALL mgProlongation()

 ndof = KNVT(J) + KNAT(J) + KNET(J) + KNEL(J)
 qScalar%ValU = qScalar%aux(J)%x(0*ndof+1:1*ndof)
 qScalar%ValV = qScalar%aux(J)%x(1*ndof+1:2*ndof)
 qScalar%ValW = qScalar%aux(J)%x(2*ndof+1:3*ndof)

 ! Set dirichlet boundary conditions on the solution
 CALL Bndry_Val()

 IF(myid.eq.showid) WRITE(*,*) "Prolongation of pressure solution to a higher level"
 MyMG%cVariable = "Pressure"

 ndof = 4*KNEL(K)
 lScalar%valP(K)%x = lScalar%valP(J)%x(1:ndof)
!  CALL ll21(lScalar%valP(K)%x,ndof,dNorm)
!  WRITE(*,*) myid, "dnorm1 is: ", dnorm
 MyMG%X   => lScalar%valP
 MyMG%AUX => lScalar%auxP

 MyMG%KNPRP => lScalar%knprp

 CALL mgProlongation()

 lScalar%valP(J)%x = MyMG%AUX(J)%x
!  CALL ll21(lScalar%valP(J)%x,ndof*8,dNorm)
!  WRITE(*,*) myid, "dnorm2 is: ", dnorm

END IF

END SUBROUTINE ProlongateSolutionSub
!
! ----------------------------------------------
!
SUBROUTINE Protocol_QuadScalar(mfile,myScalar,nINL,&
           ResU,ResV,ResW,DefScalar,RhsScalar,cTitle)
TYPE(TQuadScalar), INTENT(INOUT) :: myScalar
INTEGER nINL,mfile
INTEGER i,length
REAL*8 ResU,ResV,ResW,DefScalar,RhsScalar
CHARACTER C1*20,C2*20,C3*20,C4*20,C5*20,C6*20,C7*20,C8*20
CHARACTER, OPTIONAL:: cTitle*(*)

IF (myid.eq.showID) THEN
length =  LEN(myScalar%cName)

 C1='              '
 C4='              '
 C5='              '
 C6='              '
WRITE(C1(1:3+length),'(A3,A7)') 'Res',myScalar%cName
 C2=C1; C3=C1;
WRITE(C1(length+1:length+2),'(A2)') '_U'
WRITE(C2(length+1:length+2),'(A2)') '_V'
WRITE(C3(length+1:length+2),'(A2)') '_W'
WRITE(C4(1:3+length),'(A3,A7)') 'Def',myScalar%cName
WRITE(C5(1:7+length),'(A7,A7)') 'GlobDef',myScalar%cName
WRITE(C6(1:8+length),'(A8,A7)') 'nMGcycIt',myScalar%cName
WRITE(C7(1:8+length),'(A8,A7)') 'nMGcrsIt',myScalar%cName
WRITE(C8(1:7+length),'(A7,A7)') 'MGRates',myScalar%cName

IF (PRESENT(cTitle)) THEN
 length = LEN(cTitle)
 IF (MOD(length,2).eq.1) length = length + 1
 length = (104-length)/2
END IF

IF (nINL.EQ.0) THEN
 IF (PRESENT(cTitle)) THEN
  WRITE(MTERM,4) !cTitle
  WRITE(MFILE,4) !cTitle
 ELSE
  WRITE(MTERM,5)
  WRITE(MFILE,5)
 END IF

 WRITE(MTERM,'(A3,5(1X,A11),(1X,A12),(1X,A12),(2X,A12))')&
       "INL",TRIM(C1),TRIM(C2),TRIM(C3),TRIM(C4),TRIM(C5),TRIM(C6),TRIM(C7),TRIM(C8)
 WRITE(MFILE,'(A3,5(1X,A11),(1X,A12),(1X,A12),(2X,A12))')&
       "INL",TRIM(C1),TRIM(C2),TRIM(C3),TRIM(C4),TRIM(C5),TRIM(C6),TRIM(C7),TRIM(C8)
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(I3,5(1X,D11.4),A1)') &
       0,ResU,ResV,ResW,DefScalar,RhsScalar,"*"
 WRITE(MFILE,'(I3,5(1X,D11.4),A1)') &
       0,ResU,ResV,ResW,DefScalar,RhsScalar,"*"
ELSE
 WRITE(MTERM,'(I3,5(1X,D11.4),2(4X,I5,4X),2X,D11.4)') nINL,ResU,ResV,ResW,DefScalar,&
 RhsScalar, myScalar%prm%MGprmOut(1)%UsedIterCycle, myScalar%prm%MGprmOut(1)%nIterCoarse,&
 myScalar%prm%MGprmOut(1)%RhoMG1
 WRITE(MFILE,'(I3,5(1X,D11.4),2(4X,I5,4X),2X,D11.4)') nINL,ResU,ResV,ResW,DefScalar,&
 RhsScalar, myScalar%prm%MGprmOut(1)%UsedIterCycle, myScalar%prm%MGprmOut(1)%nIterCoarse,&
 myScalar%prm%MGprmOut(1)%RhoMG1
END IF

END IF

5  FORMAT(104('-'))
4  FORMAT(104('-'))

END SUBROUTINE Protocol_QuadScalar
!
! ----------------------------------------------
!
SUBROUTINE Protocol_LinScalar(mfile,myScalar,cTitle)
TYPE(TLinScalar), INTENT(INOUT) :: myScalar
INTEGER nINL,mfile
INTEGER i,length
character(len=32) :: C0,C1,C2,C3,C4,C5
CHARACTER, OPTIONAL:: cTitle*(*)

IF (myid.eq.showID) THEN
length =  LEN(myScalar%cName)

 C0='              '
 C1='              '
 C2='              '
 C3='              '
 C4='              '
 C5='              '
WRITE(C0(1:6+length),'(A6,A7)') 'nMGcyc',myScalar%cName
WRITE(C1(1:7+length),'(A7,A7)') 'DefInit',myScalar%cName
WRITE(C2(1:9+length),'(A9,A7)') 'CoarseIte',myScalar%cName
WRITE(C5(1:8+length),'(A8,A7)') 'DefFinal',myScalar%cName
WRITE(C3(1:6+length),'(A6,A7)') 'RhoMGN',myScalar%cName
WRITE(C4(1:6+length),'(A6,A7)') 'RhoMGA',myScalar%cName

IF (PRESENT(cTitle)) THEN
 length = LEN(cTitle)
 IF (MOD(length,2).eq.1) length = length + 1
 length = (104-length)/2
END IF

! IF (PRESENT(cTitle)) THEN
!  WRITE(MTERM,4) cTitle
!  WRITE(MFILE,4) cTitle
! ELSE
!  WRITE(MTERM,5)
!  WRITE(MFILE,5)
! END IF
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(8(1X,A13))')&
 TRIM(C0),TRIM(C1),TRIM(C5),TRIM(C2),TRIM(C3),TRIM(C4)
 WRITE(MFILE,'(8(1X,A13))')&
 TRIM(C0),TRIM(C1),TRIM(C5),TRIM(C2),TRIM(C3),TRIM(C4)
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(4XI10,2(2X,D11.4),4XI10,2(2X,D11.4))') myScalar%prm%MGprmOut%UsedIterCycle,&
       myScalar%prm%MGprmOut%DefInitial,myScalar%prm%MGprmOut%DefFinal,&
       myScalar%prm%MGprmOut%nIterCoarse,myScalar%prm%MGprmOut%RhoMG1,myScalar%prm%MGprmOut%RhoMG2
 WRITE(MFILE,'(4XI10,2(2X,D11.4),3XI10,2(2X,D11.4))') myScalar%prm%MGprmOut%UsedIterCycle,&
       myScalar%prm%MGprmOut%DefInitial,myScalar%prm%MGprmOut%DefFinal,&
       myScalar%prm%MGprmOut%nIterCoarse,myScalar%prm%MGprmOut%RhoMG1,myScalar%prm%MGprmOut%RhoMG2

END IF

5  FORMAT(104('-'))
4  FORMAT(104('-'))


END SUBROUTINE Protocol_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE QuadScP1toQ2(lSc,qSc)
TYPE(TLinScalar) lSc
TYPE(TQuadScalar) qSc

ILEV=NLMAX
CALL SETLEV(2)

lSc%Q2   = 0d0
qSc%auxU = 0d0

CALL IntP1toQ2(myQ2Coor,&
               mg_mesh%level(ilev)%kvert,&
               mg_mesh%level(ilev)%kedge,&
               mg_mesh%level(ilev)%karea,&
               mg_MlRhomat(ILEV)%a,&
               lSc%valP(ILEV)%x,&
               lSc%Q2,qSc%auxU,&
               mg_mesh%level(ilev)%nel,&
               mg_mesh%level(ilev)%nvt,&
               mg_mesh%level(ilev)%net,&
               mg_mesh%level(ilev)%nat)

END SUBROUTINE QuadScP1toQ2
!
! ----------------------------------------------
!
SUBROUTINE IntP1toQ2(dcorvg,kvert,kedge,karea,Ml,P1,Q2,Aux,nel,nvt,net,nat)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
INTEGER nel,nvt,net,nat
REAL*8  dcorvg(3,*),Ml(*),P1(4,*),Q2(*),Aux(*)
REAL*8 dC(3),dP(3),Val
INTEGER i,iel,ivt,iet,iat

DO iel=1,nel

dC(1) = 0d0
dC(2) = 0d0
dC(3) = 0d0
DO i = 1,8
 ivt = kvert(i,iel)
 dC(1) = dC(1) +  dcorvg(1,ivt)
 dC(2) = dC(2) +  dcorvg(2,ivt)
 dC(3) = dC(3) +  dcorvg(3,ivt)
END DO
dC(1) = 0.125d0*dC(1)
dC(2) = 0.125d0*dC(2)
dC(3) = 0.125d0*dC(3)

DO i = 1,8
 ivt   = kvert(i,iel)
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

DO i=1,12
 iet = kedge(i,iel)
 ivt = nvt + iet
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

DO i = 1,6
 iat   = karea(i,iel)
 ivt = nvt + net + iat
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

ivt = nvt + net + nat + iel
Val   = P1(1,iel)
Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
Aux(ivt) = Aux(ivt) + Ml(ivt)

END DO

CALL E013Sum(Q2)
CALL E013Sum(Aux)

DO i=1,nvt + net + nat + nel
Q2(i) = Q2(i)/Aux(i)
END DO

END SUBROUTINE IntP1toQ2
!
! ----------------------------------------------
!
SUBROUTINE QuadScP1toQ2Periodic(lSc,qSc)
TYPE(TLinScalar) lSc
TYPE(TQuadScalar) qSc

ILEV=NLMAX
CALL SETLEV(2)

lSc%Q2   = 0d0
qSc%auxU = 0d0

CALL IntP1toQ2Periodic(myQ2Coor,&
               mg_mesh%level(ilev)%kvert,&
               mg_mesh%level(ilev)%kedge,&
               mg_mesh%level(ilev)%karea,&
               mg_MlRhomat(ILEV)%a,&
               lSc%valP(ILEV)%x,&
               lSc%Q2,qSc%auxU,qSc%auxV,qSc%auxW,&
               mg_mesh%level(ilev)%nel,&
               mg_mesh%level(ilev)%nvt,&
               mg_mesh%level(ilev)%net,&
               mg_mesh%level(ilev)%nat)

END SUBROUTINE QuadScP1toQ2Periodic
!
! ----------------------------------------------
!
SUBROUTINE IntP1toQ2Periodic(dcorvg,kvert,kedge,karea,Ml,P1,Q2,Aux,PP,MM,nel,nvt,net,nat)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
INTEGER nel,nvt,net,nat
REAL*8  dcorvg(3,*),Ml(*),P1(4,*),Q2(*),Aux(*),PP(*),MM(*)
REAL*8 dC(3),dP(3),Val,ddsum(2)
INTEGER i,iel,ivt,iet,iat,ndof

IF (myid.ne.0) THEN
DO iel=1,nel

dC(1) = 0d0
dC(2) = 0d0
dC(3) = 0d0
DO i = 1,8
 ivt = kvert(i,iel)
 dC(1) = dC(1) +  dcorvg(1,ivt)
 dC(2) = dC(2) +  dcorvg(2,ivt)
 dC(3) = dC(3) +  dcorvg(3,ivt)
END DO
dC(1) = 0.125d0*dC(1)
dC(2) = 0.125d0*dC(2)
dC(3) = 0.125d0*dC(3)

DO i = 1,8
 ivt   = kvert(i,iel)
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

DO i=1,12
 iet = kedge(i,iel)
 ivt = nvt + iet
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

DO i = 1,6
 iat   = karea(i,iel)
 ivt = nvt + net + iat
 dP(:) = dcorvg(:,ivt)
 Val   = P1(1,iel) + (P1(2,iel)*(dP(1)-dC(1))) + (P1(3,iel)*(dP(2)-dC(2))) &
                   + (P1(4,iel)*(dP(3)-dC(3)))
 Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
 Aux(ivt) = Aux(ivt) + Ml(ivt)
END DO

ivt = nvt + net + nat + iel
Val   = P1(1,iel)
Q2(ivt)  = Q2(ivt)  + Val*Ml(ivt)
Aux(ivt) = Aux(ivt) + Ml(ivt)

END DO

ddsum(1:2) = 0d0
DO i = 1,nvt+net+nat+nel
 ddSum(1) = Q2(i) + ddSum(1)
 ddSum(2) = Aux(i) + ddSum(2)
END DO

END IF

CALL COMM_SUMMN(ddSum,2)

IF (myid.ne.0) THEN

ndof = nvt + net + nat + nel
PP(1:ndof) = dcorvg(3,1:ndof)
CALL E013PerSum(Q2,Aux,PP)

DO i=1,nvt + net + nat + nel
  Q2(i) = Q2(i)/Aux(i) - ddsum(1)/ddsum(2)
END DO

END IF

END SUBROUTINE IntP1toQ2Periodic
!
! ----------------------------------------------
!
SUBROUTINE PressureToGMV(myScalar)
TYPE(TLinScalar) myScalar

ILEV = NLMAX

CALL ProlPressure(myScalar%ValP(ILEV)%x,myScalar%ValP_GMV,&
                  KWORK(L(KLADJ(ILEV+1))),KWORK(L(KLVERT(ILEV+1))),&
                  DWORK(L(KLCVG(ILEV+1))),KNEL(ILEV))

END SUBROUTINE PressureToGMV
!
! ----------------------------------------------
!
SUBROUTINE MatStructQ2P1(ColAo,ColAn,LdA,na,nu)
IMPLICIT NONE
INTEGER ColAo(*),ColAn(*),LdA(*),na,nu
INTEGER i,j,k,nn,ijEntry

nn = 0
DO i=1,nu
 DO j=LdA(i),LdA(i+1)-1
  ijEntry = ColAo(j)
  DO k=1,4
   nn = nn + 1
   ColAn(nn) = 4*(ijEntry-1) + k
  END DO
 END DO
END DO

DO i=2,nu+1
 LdA(i) = 4*(LdA(i)-1)+1
END DO

END SUBROUTINE MatStructQ2P1
!
! ----------------------------------------------
!
SUBROUTINE MatStructP1Q2(LdAo,LdAn,ColAo,ColAn,na,nu)
IMPLICIT NONE
INTEGER ColAo(*),ColAn(*),LdAo(*),LdAn(*),na,nu
INTEGER i,j,k,nn,iPos

nn = 0
DO i=1,nu
 DO k=1,4
  DO j=LdAo(i),LdAo(i+1)-1
   nn = nn + 1
   ColAn(nn) = ColAo(j)
  END DO
 END DO
END DO

nn = 0
iPos = -26
DO i=1,nu
 DO k=1,4
  nn = nn + 1
  iPos = iPos + 27
  LdAn(nn) = iPos
 END DO
END DO
nn = nn + 1
LdAn(nn) = LdAn(nn-1) + 27

END SUBROUTINE MatStructP1Q2
!
! ----------------------------------------------
!
SUBROUTINE OutputMatrixStuct(cFile,myMat)
TYPE(TMatrix) myMat
CHARACTER*4 cFile
CHARACTER*10 myFile
INTEGER I,J


IF (myid.NE.0) THEN

 IF (myid.LT.10) THEN
  WRITE(myFile(1:10),'(A4,A1,I1,A4)') cFile,"0",myid,".txt"
 ELSE
  WRITE(myFile(1:10),'(A4,I2,A4)') cFile,myid,".txt"
 END IF
!  WRITE(*,*) myid,myFile
 OPEN(987,FILE=myFile)

 DO I=1,myMat%nu
  WRITE(987,'(I10,I10,A1,1000I12)') I,myMat%LdA(I+1)-myMat%LdA(I),&
!   myMat%LdA(I)
  "|",(myMat%ColA(J),J=myMat%LdA(I),myMat%LdA(I+1)-1)
 END DO
 CLOSE(987)

END IF

END SUBROUTINE OutputMatrixStuct
!
! ----------------------------------------------
!
SUBROUTINE OutputMatrix(cFile,myMat,Mat,II)
!  CALL OutputMatrix("MMAT",mg_qMat(2),mg_MMat(2)%a,2)
TYPE(TMatrix) myMat
REAL*8 Mat(*),DD
CHARACTER*4 cFile
CHARACTER*12 myFile
INTEGER I,J,II


! IF (myid.NE.0) THEN

 IF (myid.LT.10) THEN
  WRITE(myFile(1:12),'(A4,A1,I1,A1,I1,A4)') cFile,"L",II,"0",myid,".txt"
 ELSE
  WRITE(myFile(1:12),'(A4,A1,I1,I2,A4)') cFile,"L",II,myid,".txt"
 END IF
!  WRITE(*,*) myid,myFile
 OPEN(987,FILE=myFile)

 DO I=1,myMat%nu
!   WRITE(987,'(I12,A1,I12)') I,"|",myMat%LdA(I+1)-myMat%LdA(I)
  DD = 0d0
  DO J=myMat%LdA(I),myMat%LdA(I+1)-1
   DD = DD + Mat(J)
  END DO
  WRITE(987,'(I10,I10,A1,1000D12.4)') I,myMat%LdA(I+1)-myMat%LdA(I),&
  "|",MAX(1d-18,ABS(DD)),(MAX(1d-18,ABS(Mat(J))),J=myMat%LdA(I),myMat%LdA(I+1)-1)
 END DO
 CLOSE(987)

! END IF

END SUBROUTINE OutputMatrix
!
! ----------------------------------------------
!
SUBROUTINE EvaluateDragLift9_old(U,V,W,P,CYL,dfV,dfP)
REAL*8 U(*),V(*),W(*),P(*)
LOGICAL CYL(*)
REAL*8 dfV(3),dfP(3)
INTEGER I,J,K,JJ

IF (myid.ne.0) THEN

ILEV=NLMAX
CALL SETLEV(2)

S11Mat   => mg_S11Mat(ILEV)%a
S22Mat   => mg_S22Mat(ILEV)%a
S33Mat   => mg_S33Mat(ILEV)%a
S12Mat   => mg_S12Mat(ILEV)%a
S13Mat   => mg_S13Mat(ILEV)%a
S23Mat   => mg_S23Mat(ILEV)%a
S21Mat   => mg_S21Mat(ILEV)%a
S31Mat   => mg_S31Mat(ILEV)%a
S32Mat   => mg_S32Mat(ILEV)%a
BXMat    => mg_BXMat(ILEV)%a
BYMat    => mg_BYMat(ILEV)%a
BZMat    => mg_BZMat(ILEV)%a
qMat     => mg_qMat(ILEV)

dfV = 0d0; dfP = 0d0  
DO I=1,qMat%nu
 IF (CYL(I)) THEN
  DO J=qMat%LdA(I),qMat%LdA(I+1)-1
   K = qMat%ColA(J)
   dfV(1) = dfV(1) - (S11Mat(J)*U(K) + S12Mat(J)*V(K) + S13Mat(J)*W(K))
   dfV(2) = dfV(2) - (S21Mat(J)*U(K) + S22Mat(J)*V(K) + S23Mat(J)*W(K))
   dfV(3) = dfV(3) - (S31Mat(J)*U(K) + S32Mat(J)*V(K) + S33Mat(J)*W(K))
  END DO
  DO J=qlMat%LdA(I),qlMat%LdA(I+1)-1
   K = qlMat%ColA(J)
   dfP(1) = dfP(1) + BXMat(J)*P(K)
   dfP(2) = dfP(2) + BYMat(J)*P(K)
   dfP(3) = dfP(3) + BZMat(J)*P(K)
  END DO
 END IF
END DO

END IF

CALL COMM_SUMM(dfV(1))
CALL COMM_SUMM(dfV(2))
CALL COMM_SUMM(dfV(3))
CALL COMM_SUMM(dfP(1))
CALL COMM_SUMM(dfP(2))
CALL COMM_SUMM(dfP(3))

END SUBROUTINE EvaluateDragLift9_old
!
! ----------------------------------------------
!
SUBROUTINE EvaluateDragLift_old(U,V,W,P,CYL,dfV,dfP)
REAL*8 U(*),V(*),W(*),P(*)
LOGICAL CYL(*)
REAL*8 dfV(3),dfP(3)
INTEGER I,J,K,JJ

IF (myid.ne.0) THEN

ILEV=NLMAX
CALL SETLEV(2)

DMat     => mg_DMat(ILEV)%a
BXMat    => mg_BXMat(ILEV)%a
BYMat    => mg_BYMat(ILEV)%a
BZMat    => mg_BZMat(ILEV)%a
qMat     => mg_qMat(ILEV)

! CALL GetParPressure(P,PP)

dfV = 0d0; dfP = 0d0  
DO I=1,qMat%nu
 IF (CYL(I)) THEN
  DO J=qMat%LdA(I),qMat%LdA(I+1)-1
   K = qMat%ColA(J)
   dfV(1) = dfV(1) - DMat(J)*U(K)
   dfV(2) = dfV(2) - DMat(J)*V(K)
   dfV(3) = dfV(3) - DMat(J)*W(K)
  END DO
  DO J=qlMat%LdA(I),qlMat%LdA(I+1)-1
   K = qlMat%ColA(J)
   dfP(1) = dfP(1) + BXMat(J)*P(K)
   dfP(2) = dfP(2) + BYMat(J)*P(K)
   dfP(3) = dfP(3) + BZMat(J)*P(K)
  END DO
 END IF
END DO

! WRITE(*,'(2I8,3D12.4)') myid,JJ,df

END IF

CALL COMM_SUMM(dfV(1))
CALL COMM_SUMM(dfV(2))
CALL COMM_SUMM(dfV(3))
CALL COMM_SUMM(dfP(1))
CALL COMM_SUMM(dfP(2))
CALL COMM_SUMM(dfP(3))

END SUBROUTINE EvaluateDragLift_old
!
! ----------------------------------------------
!
SUBROUTINE EvaluateDragLift(U,V,W,P,Nu,FBM)
REAL*8 U(*),V(*),W(*),P(*),Nu(*)
INTEGER FBM(*)
INTEGER I,J,K,JJ
EXTERNAL E013

ILEV=NLMAX
CALL SETLEV(2)
CALL GetForces(Properties%ForceScale,U,V,W,P,FBM,Nu,&
               mg_mesh%level(ILEV)%kvert,&
               mg_mesh%level(ILEV)%karea,&
               mg_mesh%level(ILEV)%kedge,&
               mg_mesh%level(ILEV)%dcorvg,&
               E013)

END SUBROUTINE EvaluateDragLift
!
! ----------------------------------------------
!
SUBROUTINE GetVeloParameters(myParam,cName,mfile)
IMPLICIT NONE
TYPE(tParamV) myParam
CHARACTER*7 cName
INTEGER mfile
INTEGER iEnd,iAt,iEq,istat
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile=90909090
CHARACTER(20) :: out_string
  
OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action='read',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open data file: ",myDataFile  
  stop          
end if

! IF (myid.eq.showid) WRITE(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
! IF (myid.eq.showid) WRITE(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
myParam%MGprmIn%RLX = 0.66d0
myParam%MGprmIn%CrsSolverType = 1

DO
 READ (UNIT=myFile,FMT='(A100)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 CALL StrStuct()
!   IF (myid.eq.showid) write(*,*) myid,iAt,iEq,bOK
  IF (bOK) THEN

  READ(string(1:iAt-1),*) cVar

!   IF (myid.eq.showid) write(*,*) myid,cVar,iAt,iEq,string
  IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN

   READ(string(iAt+1:iEq-1),*) cPar
!    IF (myid.eq.showid) write(*,*) myid,cPar

   SELECT CASE (TRIM(ADJUSTL(cPar)))

    ! Problems with GCC: 
    ! The gfortran compiler does not permit a statement like write(*,'(A,I)') where
    ! the length of the integer I is not specified.
    ! A way to solve this problem is to write the integer to a string:
    ! write(string,'(i20)')myParam%iMass
    ! that is long enough to hold all reasonable integer values
    ! then adjust the length of the string write it out:
    ! write(*,'(A)')adjustl(string)
    CASE ("iMass")
    READ(string(iEq+1:),*) myParam%iMass
    write(out_string,'(i20)')myParam%iMass
    call write_param_int(mfile,cVar,cPar,out_string,myParam%iMass)

    CASE ("SolvType")
    READ(string(iEq+1:),*) param
    IF (TRIM(ADJUSTL(param)).EQ."Jacobi") myParam%SolvType = 1
    IF (TRIM(ADJUSTL(param)).EQ."BiCGSTAB")myParam%SolvType = 2
    call write_param_int(mfile,cVar,cPar,out_string,myParam%SolvType)

    CASE ("defCrit")
    READ(string(iEq+1:),*) myParam%defCrit
    call write_param_real(mfile,cVar,cPar,out_string,myParam%defCrit)

    CASE ("Alpha")
    READ(string(iEq+1:),*) myParam%Alpha
    call write_param_real(mfile,cVar,cPar,out_string,myParam%Alpha)

    CASE ("MinDef")
    READ(string(iEq+1:),*) myParam%MinDef
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MinDef)

    CASE ("NLmin")
    READ(string(iEq+1:),*) myParam%NLmin
    call write_param_int(mfile,cVar,cPar,out_string,myParam%NLmin)

    CASE ("NLmax")
    READ(string(iEq+1:),*) myParam%NLmax
    call write_param_int(mfile,cVar,cPar,out_string,myParam%NLmax)

    CASE ("MGMinLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinLev
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MinLev)

    CASE ("MGMedLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MedLev
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MedLev)

    CASE ("MGMinIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MinIterCycle)

    CASE ("MGMaxIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MaxIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MaxIterCycle)

    CASE ("MGSmoothSteps")
    READ(string(iEq+1:),*) myParam%MGprmIn%nSmootherSteps
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%nSmootherSteps)

    CASE ("MGIterCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%nIterCoarse
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%nIterCoarse)

    CASE ("MG_VANKA")
    READ(string(iEq+1:),*) myParam%MGprmIn%VANKA
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%VANKA)

    CASE ("MGCrsRelaxPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%CrsRelaxPrm
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%CrsRelaxPrm)

    CASE ("MGCrsRelaxParPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%CrsRelaxParPrm
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%CrsRelaxParPrm)

    CASE ("MGCrsSolverType")
    READ(string(iEq+1:),*) myParam%MGprmIn%CrsSolverType
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%CrsSolverType)

    CASE ("MGDefImprCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%DefImprCoarse
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%DefImprCoarse)

    CASE ("MGCriterion1")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion1
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%Criterion1)

    CASE ("MGCriterion2")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion2
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%Criterion2)

    CASE ("MGCycType")
    READ(string(iEq+1:),*) param
    myParam%MGprmIn%CycleType = TRIM(ADJUSTL(param))
    IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType
    IF (myid.eq.showid) write(mfile,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType

    CASE ("MGRelaxPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%RLX
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%RLX)
  END SELECT

  END IF
 END IF
END DO

! IF (myid.eq.showid) write(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
! IF (myid.eq.showid) write(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

CLOSE (myFile)

myParam%MGprmIn%MaxLev = NLMAX
out_string = ""
write(out_string,'(i20)')myParam%iMass
call write_param_int(mfile,"Velo","MGMaxLev",out_string,myParam%MGprmIn%MaxLev)

CONTAINS

SUBROUTINE StrStuct()
IMPLICIT NONE
INTEGER i,n

n = len(string)
iAt = 0
iEq = 0
DO i=1,n
 IF (string(i:i).EQ. '@') iAt = i
 IF (string(i:i).EQ. '=') iEq = i
END DO

bOk=.FALSE.
IF (iAt.ne.0.AND.iEq.ne.0) bOk=.TRUE.

END SUBROUTINE StrStuct

END SUBROUTINE GetVeloParameters
!
! ----------------------------------------------
!
subroutine write_param_real(unitid,cVar,cPar,out_string,val)
  implicit none
  integer :: unitid
  character(len=*) :: cVar
  character(len=*) :: cPar
  character(len=*) :: out_string
  Real*8 :: val

  out_string = " "
  write(out_string,'(E16.8)')val
  IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",adjustl(out_string)
  IF (myid.eq.showid) write(unitid,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",adjustl(out_string)

end subroutine write_param_real
!
! ----------------------------------------------
!
subroutine write_param_int(unitid,cVar,cPar,out_string,val)
  implicit none
  integer :: unitid
  character(len=*) :: cVar
  character(len=*) :: cPar
  character(len=*) :: out_string
  integer :: val

  out_string = " "
  write(out_string,'(I20)')val
  IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",adjustl(out_string)
  IF (myid.eq.showid) write(unitid,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",adjustl(out_string)

end subroutine write_param_int
!
! ----------------------------------------------
!
SUBROUTINE GetPresParameters(myParam,cName,mfile)
IMPLICIT NONE
TYPE(tParamP) myParam
INTEGER mfile
CHARACTER*7 cName
INTEGER iEnd,iAt,iEq,istat
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile=90909090
character(len=20) :: out_string
  
OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action='read',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open data file: ",myDataFile  
  stop          
end if

!IF (myid.eq.showid) write(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
!IF (myid.eq.showid) write(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
myParam%MGprmIn%RLX = 0.66d0

DO
 READ (UNIT=myFile,FMT='(A100)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 CALL StrStuct()
!   !IF (myid.eq.showid) write(*,*) myid,iAt,iEq,bOK
  IF (bOK) THEN

  READ(string(1:iAt-1),*) cVar

!   !IF (myid.eq.showid) write(*,*) myid,cVar,iAt,iEq,string
  IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN

   READ(string(iAt+1:iEq-1),*) cPar

!    !IF (myid.eq.showid) write(*,*) myid,cPar
   SELECT CASE (TRIM(ADJUSTL(cPar)))

    CASE ("MGMinLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinLev
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MinLev)

    CASE ("MGMedLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MedLev
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MedLev)

    CASE ("MGMinIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MinIterCycle)

    CASE ("MGMaxIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MaxIterCycle
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%MaxIterCycle)

    CASE ("MGSmoothSteps")
    READ(string(iEq+1:),*) myParam%MGprmIn%nSmootherSteps
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%nSmootherSteps)

    CASE ("MGIterCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%nIterCoarse
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%nIterCoarse)

    CASE ("MGDefImprCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%DefImprCoarse
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%DefImprCoarse)

    CASE ("MGCrsRelaxPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%CrsRelaxPrm
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%CrsRelaxPrm)

    CASE ("MGCriterion1")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion1
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%Criterion1)

    CASE ("MGCriterion2")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion2
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%Criterion2)

    CASE ("MGCycType")
    READ(string(iEq+1:),*) param
    myParam%MGprmIn%CycleType = TRIM(ADJUSTL(param))
    IF (myid.eq.showid) write(mterm,'(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType
    IF (myid.eq.showid) write(mfile,'(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType

    CASE ("MGCrsSolverType")
    READ(string(iEq+1:),*) myParam%MGprmIn%CrsSolverType
    call write_param_int(mfile,cVar,cPar,out_string,myParam%MGprmIn%CrsSolverType)
    
    CASE ("MGRelaxPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%RLX
    call write_param_real(mfile,cVar,cPar,out_string,myParam%MGprmIn%RLX)

   END SELECT

  END IF
 END IF
END DO

!IF (myid.eq.showid) write(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
!IF (myid.eq.showid) write(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

CLOSE (myFile)

CONTAINS

SUBROUTINE StrStuct()
IMPLICIT NONE
INTEGER i,n

n = len(string)
iAt = 0
iEq = 0
DO i=1,n
 IF (string(i:i).EQ. '@') iAt = i
 IF (string(i:i).EQ. '=') iEq = i
END DO

bOk=.FALSE.
IF (iAt.ne.0.AND.iEq.ne.0) bOk=.TRUE.

END SUBROUTINE StrStuct

END SUBROUTINE GetPresParameters
!
! ----------------------------------------------
!
SUBROUTINE GetPhysiclaParameters(Props,cName,mfile)
IMPLICIT NONE
TYPE(tProperties) Props
INTEGER mfile
CHARACTER*7 cName

INTEGER iEnd,iAt,iEq,istat
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile=90909090
  
OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action='read',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open data file: ",myDataFile  
  stop          
end if

IF (myid.eq.showid) WRITE(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
IF (myid.eq.showid) WRITE(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

DO
 READ (UNIT=myFile,FMT='(A100)',IOSTAT=iEnd) string
 IF (iEnd.EQ.-1) EXIT
 CALL StrStuct()
  IF (bOK) THEN

  READ(string(1:iAt-1),*) cVar

  IF (TRIM(ADJUSTL(cVar)).EQ.TRIM(ADJUSTL(cName))) THEN

   READ(string(iAt+1:iEq-1),*) cPar
   SELECT CASE (TRIM(ADJUSTL(cPar)))

    CASE ("Material")
    READ(string(iEq+1:),*) Props%Material
    IF (myid.eq.showid) write(mterm,'(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Material
    IF (myid.eq.showid) write(mfile,'(A,A)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Material
    CASE ("Gravity")
    READ(string(iEq+1:),*) Props%Gravity
    IF (myid.eq.showid) write(mterm,'(A,3E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Gravity
    IF (myid.eq.showid) write(mfile,'(A,3E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Gravity
    CASE ("Density")
    READ(string(iEq+1:),*) Props%Density
    IF (myid.eq.showid) write(mterm,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Density
    IF (myid.eq.showid) write(mfile,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Density
    CASE ("Viscosity")
    READ(string(iEq+1:),*) Props%Viscosity
    IF (myid.eq.showid) write(mterm,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Viscosity
    IF (myid.eq.showid) write(mfile,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Viscosity
    CASE ("DiffCoeff")
    READ(string(iEq+1:),*) Props%DiffCoeff
    IF (myid.eq.showid) write(mterm,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%DiffCoeff
    IF (myid.eq.showid) write(mfile,'(A,2E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%DiffCoeff
    CASE ("Sigma")
    READ(string(iEq+1:),*) Props%Sigma
    IF (myid.eq.showid) write(mterm,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Sigma
    IF (myid.eq.showid) write(mfile,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%Sigma
    CASE ("DiracEps")
    READ(string(iEq+1:),*) Props%DiracEps
    IF (myid.eq.showid) write(mterm,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%DiracEps
    IF (myid.eq.showid) write(mfile,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%DiracEps
    CASE ("PowerLawExp")
    READ(string(iEq+1:),*) Props%PowerLawExp
    IF (myid.eq.showid) write(mterm,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%PowerLawExp
    IF (myid.eq.showid) write(mfile,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%PowerLawExp
    CASE ("ViscoLambda")
    READ(string(iEq+1:),*) Props%ViscoLambda
    IF (myid.eq.showid) write(mterm,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ViscoLambda
    IF (myid.eq.showid) write(mfile,'(A,E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ViscoLambda
    CASE ("ViscoModel")
    READ(string(iEq+1:),*) Props%ViscoModel
    IF (myid.eq.showid) write(mterm,'(A,I2)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ViscoModel
    IF (myid.eq.showid) write(mfile,'(A,I2)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ViscoModel
    CASE ("ForceScale")
    READ(string(iEq+1:),*) Props%ForceScale
    IF (myid.eq.showid) write(mterm,'(A,6E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ForceScale
    IF (myid.eq.showid) write(mfile,'(A,6E16.8)') cVar//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",Props%ForceScale
  END SELECT

  END IF
 END IF
END DO

IF (myid.eq.showid) WRITE(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
IF (myid.eq.showid) WRITE(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

CLOSE (myFile)

CONTAINS

SUBROUTINE StrStuct()
IMPLICIT NONE
INTEGER i,n

n = len(string)
iAt = 0
iEq = 0
DO i=1,n
 IF (string(i:i).EQ. '@') iAt = i
 IF (string(i:i).EQ. '=') iEq = i
END DO

bOk=.FALSE.
IF (iAt.ne.0.AND.iEq.ne.0) bOk=.TRUE.

END SUBROUTINE StrStuct

END SUBROUTINE GetPhysiclaParameters
!
! ----------------------------------------------
!
!SUBROUTINE Correct_myQ2Coor(dcorvg,kvert,karea,kedge)
!REAL*8  dcorvg(3,*)
!INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!REAL*8 PX,PY,PZ
!INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,iedge1,iedge2
!INTEGER iat1,iat2,iat3,iat4,iat5,iat6
!INTEGER iBndr
!INTEGER NeighE(2,12),NeighA(4,6),NeighU(4,6)
!DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!DATA NeighU/1,2,3,4,  1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/
!
!! return
!k=1
!DO i=1,nel
! DO j=1,6
!  IF (k.eq.karea(j,i)) THEN
!   ivt1 = kvert(NeighA(1,j),i)
!   ivt2 = kvert(NeighA(2,j),i)
!   ivt3 = kvert(NeighA(3,j),i)
!   ivt4 = kvert(NeighA(4,j),i)
!
!   IF ((     myBoundary%LS_zero(ivt1).and.myBoundary%LS_zero(ivt2)).and.&
!       .not.(myBoundary%LS_zero(ivt3).and.myBoundary%LS_zero(ivt4))) then
!    iedge1 = neighU(1,j)
!    iedge2 = neighU(3,j)
!    PX = 0.5d0*(myQ2Coor(1,nvt+kedge(iedge1,i))+myQ2Coor(1,nvt+kedge(iedge2,i)))
!    PY = 0.5d0*(myQ2Coor(2,nvt+kedge(iedge1,i))+myQ2Coor(2,nvt+kedge(iedge2,i)))
!    PZ = 0.5d0*(myQ2Coor(3,nvt+kedge(iedge1,i))+myQ2Coor(3,nvt+kedge(iedge2,i)))
!!     WRITE(*,*) 'case: ',1
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!    myQ2Coor(:,nvt+net+k)=[PX,PY,PZ]
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!   END IF
!
!   IF ((     myBoundary%LS_zero(ivt3).and.myBoundary%LS_zero(ivt4)).and.&
!       .not.(myBoundary%LS_zero(ivt1).and.myBoundary%LS_zero(ivt2))) then
!    iedge1 = neighU(1,j)
!    iedge2 = neighU(3,j)
!    PX = 0.5d0*(myQ2Coor(1,nvt+kedge(iedge1,i))+myQ2Coor(1,nvt+kedge(iedge2,i)))
!    PY = 0.5d0*(myQ2Coor(2,nvt+kedge(iedge1,i))+myQ2Coor(2,nvt+kedge(iedge2,i)))
!    PZ = 0.5d0*(myQ2Coor(3,nvt+kedge(iedge1,i))+myQ2Coor(3,nvt+kedge(iedge2,i)))
!!     WRITE(*,*) 'case: ',2
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!    myQ2Coor(:,nvt+net+k)=[PX,PY,PZ]
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!   END IF
!
!   
!   IF ((     myBoundary%LS_zero(ivt2).and.myBoundary%LS_zero(ivt3)).and.&
!       .not.(myBoundary%LS_zero(ivt1).and.myBoundary%LS_zero(ivt4))) then
!    iedge1 = neighU(2,j)
!    iedge2 = neighU(4,j)
!    PX = 0.5d0*(myQ2Coor(1,nvt+kedge(iedge1,i))+myQ2Coor(1,nvt+kedge(iedge2,i)))
!    PY = 0.5d0*(myQ2Coor(2,nvt+kedge(iedge1,i))+myQ2Coor(2,nvt+kedge(iedge2,i)))
!    PZ = 0.5d0*(myQ2Coor(3,nvt+kedge(iedge1,i))+myQ2Coor(3,nvt+kedge(iedge2,i)))
!!     WRITE(*,*) 'case: ',3
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!    myQ2Coor(:,nvt+net+k)=[PX,PY,PZ]
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!   END IF
!
!   IF ((     myBoundary%LS_zero(ivt1).and.myBoundary%LS_zero(ivt4)).and.&
!       .not.(myBoundary%LS_zero(ivt2).and.myBoundary%LS_zero(ivt3))) then
!    iedge1 = neighU(2,j)
!    iedge2 = neighU(4,j)
!    PX = 0.5d0*(myQ2Coor(1,nvt+kedge(iedge1,i))+myQ2Coor(1,nvt+kedge(iedge2,i)))
!    PY = 0.5d0*(myQ2Coor(2,nvt+kedge(iedge1,i))+myQ2Coor(2,nvt+kedge(iedge2,i)))
!    PZ = 0.5d0*(myQ2Coor(3,nvt+kedge(iedge1,i))+myQ2Coor(3,nvt+kedge(iedge2,i)))
!!     WRITE(*,*) 'case: ',4
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!    myQ2Coor(:,nvt+net+k)=[PX,PY,PZ]
!!     WRITE(*,*) myQ2Coor(:,nvt+net+k)
!   END IF
!
!   k = k + 1
!  END IF
! END DO
!END DO
!
!! pause
!DO i=1,nel
!
! iat1 = karea(1,i)
! iat2 = karea(2,i)
! iat3 = karea(3,i)
! iat4 = karea(4,i)
! iat5 = karea(5,i)
! iat1 = karea(6,i)
!
! IF (myBoundary%LS_zero(nvt+net+iat1).and.(.not.myBoundary%LS_zero(nvt+net+iat6))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat1)+myQ2Coor(1,nvt+net+iat6))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat1)+myQ2Coor(2,nvt+net+iat6))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat1)+myQ2Coor(3,nvt+net+iat6))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
! IF (myBoundary%LS_zero(nvt+net+iat2).and.(.not.myBoundary%LS_zero(nvt+net+iat4))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat2)+myQ2Coor(1,nvt+net+iat4))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat2)+myQ2Coor(2,nvt+net+iat4))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat2)+myQ2Coor(3,nvt+net+iat4))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
! IF (myBoundary%LS_zero(nvt+net+iat3).and.(.not.myBoundary%LS_zero(nvt+net+iat5))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat3)+myQ2Coor(1,nvt+net+iat5))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat3)+myQ2Coor(2,nvt+net+iat5))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat3)+myQ2Coor(3,nvt+net+iat5))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
! IF (myBoundary%LS_zero(nvt+net+iat6).and.(.not.myBoundary%LS_zero(nvt+net+iat1))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat1)+myQ2Coor(1,nvt+net+iat6))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat1)+myQ2Coor(2,nvt+net+iat6))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat1)+myQ2Coor(3,nvt+net+iat6))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
! IF (myBoundary%LS_zero(nvt+net+iat4).and.(.not.myBoundary%LS_zero(nvt+net+iat2))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat2)+myQ2Coor(1,nvt+net+iat4))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat2)+myQ2Coor(2,nvt+net+iat4))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat2)+myQ2Coor(3,nvt+net+iat4))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
! IF (myBoundary%LS_zero(nvt+net+iat5).and.(.not.myBoundary%LS_zero(nvt+net+iat3))) then
!  PX = 0.5d0*(myQ2Coor(1,nvt+net+iat3)+myQ2Coor(1,nvt+net+iat5))
!  PY = 0.5d0*(myQ2Coor(2,nvt+net+iat3)+myQ2Coor(2,nvt+net+iat5))
!  PZ = 0.5d0*(myQ2Coor(3,nvt+net+iat3)+myQ2Coor(3,nvt+net+iat5))
!  myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
! end if
!
!END DO
!
!END SUBROUTINE Correct_myQ2Coor
!
! ----------------------------------------------
!
SUBROUTINE SetUp_myQ2Coor(dcorvg,dcorag,kvert,karea,kedge)
REAL*8  dcorvg(3,*),dcorag(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
REAL*8 PX,PY,PZ
INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
INTEGER iBndr
INTEGER NeighE(2,12),NeighA(4,6)
DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

DO i=1,nvt
 PX = dcorvg(1,I)
 PY = dcorvg(2,I)
 PZ = dcorvg(3,I)
 myQ2Coor(:,i)=[PX,PY,PZ]
END DO

k=1
DO i=1,nel
 DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
   ivt1 = kvert(NeighE(1,j),i)
   ivt2 = kvert(NeighE(2,j),i)
   PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
   PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
   PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
   myQ2Coor(:,nvt+k)=[PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO

k=1
DO i=1,nel
 DO j=1,6
  IF (k.eq.karea(j,i)) THEN
   ivt1 = kvert(NeighA(1,j),i)
   ivt2 = kvert(NeighA(2,j),i)
   ivt3 = kvert(NeighA(3,j),i)
   ivt4 = kvert(NeighA(4,j),i)
   PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
   PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
   PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
   myQ2Coor(:,nvt+net+k)=[PX,PY,PZ]
   k = k + 1
  END IF
 END DO
END DO

DO i=1,nel
 PX = 0d0
 PY = 0d0
 PZ = 0d0
 DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
 END DO
 myQ2Coor(:,nvt+net+nat+i)=[PX,PY,PZ]
END DO

END SUBROUTINE SetUp_myQ2Coor
!
! ----------------------------------------------
!
SUBROUTINE ResetTimer()

 myStat%iNonLin=0
 myStat%iLinUVW=0
 myStat%iLinP=0
 myStat%tMGUVW=0d0
 myStat%tMGP=0d0
 myStat%tDefUVW=0d0
 myStat%tDefP=0d0
 myStat%tCorrUVWP=0d0
 myStat%tGMVOut=0d0
 myStat%tDumpOut=0d0
 myStat%tSmat=0d0
 myStat%tKmat=0d0
 myStat%tDmat=0d0
 myStat%tMmat=0d0
 myStat%tCmat=0d0
 myStat%tRestUVW=0d0
 myStat%tProlUVW=0d0
 myStat%tSmthUVW=0d0
 myStat%tSolvUVW=0d0
 myStat%tRestP=0d0
 myStat%tProlP=0d0
 myStat%tSmthP=0d0
 myStat%tSolvP=0d0
 myStat%tCommP=0d0
 myStat%tCommV=0d0
 myStat%tCommS=0d0

END SUBROUTINE ResetTimer

END MODULE def_QuadScalar

