!
! ----------------------------------------------
!
SUBROUTINE Create_AMat_GenLinSc_Q1(lSc)
implicit none
TYPE(lScalarGen) :: lSc
integer NA,iFld

IF (.NOT.ALLOCATED(mg_AMat)) ALLOCATE(mg_AMat(lSc%nOfFields))

DO iFld=1,lSc%nOfFields
 IF (.NOT.ALLOCATED(mg_AMat(iFld)%Fld)) ALLOCATE(mg_AMat(iFld)%Fld(NLMIN:NLMAX))
 DO ILEV=NLMIN,NLMAX
  NA = mg_lMat(ILEV)%na
  IF (.NOT.ALLOCATED(mg_AMat(iFld)%Fld(ILEV)%a)) ALLOCATE(mg_AMat(iFld)%Fld(ILEV)%a(NA))
 END DO
END DO

END SUBROUTINE Create_AMat_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Initialize_GenLinSc_Q1(myScalar)
USE var_QuadScalar, only:knvt
TYPE(lScalarGen) myScalar
integer iFld,nFld,ndof,MaxLev

ILEV=NLMAX
CALL SETLEV(2)

plMat   => mg_lMat(ILEV)

myScalar%ndof = KNVT(ILEV)

MaxLev = myScalar%prm%MGprmIn%MaxLev - myScalar%prm%MGprmIn%MaxDifLev + 1
if (myid.eq.0) then
 ndof = KNVT(ILEV)
else
 ndof = KNVT(MaxLev)
end if

myScalar%na = plMat%na
nFld = myScalar%nOfFields

ALLOCATE(myScalar%Fld(nFld))

DO iFld = 1,nFld
 myScalar%Fld(iFld)%cName = myScalar%prm%cField(iFld)
 ALLOCATE(myScalar%Fld(iFld)%knpr(ndof))
 myScalar%Fld(iFld)%knpr = 0
 ALLOCATE(myScalar%Fld(iFld)%def(ndof))
 ALLOCATE(myScalar%Fld(iFld)%aux(ndof))
 ALLOCATE(myScalar%Fld(iFld)%rhs(ndof))
 ALLOCATE(myScalar%Fld(iFld)%val(ndof))
 ALLOCATE(myScalar%Fld(iFld)%val_old(ndof))
END DO

ALLOCATE(myScalar%knpr(nFld*ndof))

ALLOCATE(myScalar%def(NLMIN:MaxLev))
ALLOCATE(myScalar%sol(NLMIN:MaxLev))
ALLOCATE(myScalar%aux(NLMIN:MaxLev))
ALLOCATE(myScalar%rhs(NLMIN:MaxLev))
ALLOCATE(MGE011(NLMIN:MaxLev))

DO ILEV=NLMIN,MaxLev
 ALLOCATE(myScalar%def(ILEV)%x(nFld*KNVT(ILEV)))
 ALLOCATE(myScalar%sol(ILEV)%x(nFld*KNVT(ILEV)))
 ALLOCATE(myScalar%aux(ILEV)%x(nFld*KNVT(ILEV)))
 ALLOCATE(myScalar%rhs(ILEV)%x(nFld*KNVT(ILEV)))
 ALLOCATE(MGE011(ILEV)%UE(nFld))
 DO iFld = 1,nFld
  ALLOCATE(MGE011(ILEV)%UE(iFld)%x(KNVT(ILEV)))
 END DO
END DO

END SUBROUTINE Initialize_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Matdef_INSTATIONARY_GenLinSc_Q1(lSc,idef,imat)
INTEGER :: idef,imat
TYPE(lScalarGen) lSc
INTEGER i,j,iFld
REAL*8 daux
EXTERNAL E011

 if (imat.eq.1) then
  DO ILEV=NLMIN,NLMAX

   CALL SETLEV(2)

   LaplaceMat    => mg_LaplaceMat(ILEV)%a
   LMassMat      => mg_LMassMat(ILEV)%a
   ConvectionMat =>  mg_ConvMat(ILEV)%a
   plMat  => mg_lMat(ILEV)

  ! Build up the matrix
   DO iFld=1,lSc%nOfFields
    DO I=1,plMat%nu
     J = plMat%LdA(I)
     daux = LMassMat(I) - thstep*(LaplaceMat(J) + ConvectionMat(J))
     mg_AMat(iFld)%fld(ILEV)%a(J) = daux
     DO J=plMat%LdA(I)+1,plMat%LdA(I+1)-1
      daux = -thstep*(LaplaceMat(J) + ConvectionMat(J))
      mg_AMat(iFld)%fld(ILEV)%a(J) =  daux
     END DO
    END DO
   END DO
   
  END DO
 end if

 ILEV=NLMAX
 CALL SETLEV(2)

 LaplaceMat    => mg_LaplaceMat(ILEV)%a
 LMassMat      => mg_LMassMat(ILEV)%a
 ConvectionMat => mg_ConvMat(ILEV)%a
 plMat  => mg_lMat(ILEV)

 ! Build up the defect
 IF (idef.eq. 1) THEN
  DO iFld=1,lSc%nOfFields
  
   lSc%fld(iFld)%def = 0d0
   DO j=1,plMat%nu
    lSc%fld(iFld)%def(j) = lSc%fld(iFld)%def(j) + LMassMat(j)*lSc%fld(iFld)%val(j)
   END DO
   
   CALL LAX17(ConvectionMat,plMat%ColA,plMat%LdA,plMat%nu,&
   lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep,1d0)
   
   CALL LAX17(LaplaceMat,plMat%ColA,plMat%LdA,plMat%nu,&
   lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep,1d0)
   
  END DO
  
 ELSE

  DO iFld=1,lSc%nOfFields
   CALL LAX17(mg_AMat(iFld)%fld(ILEV)%a,plMat%ColA,plMat%LdA,plMat%nu,&
   lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,-1d0,1d0)
  END DO
  
 END IF
 
 IF (idef.le.0) then
  DO iFld=1,lSc%nOfFields
   CALL DefTVD_LinScalar(lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep)
  END DO
 end if


END SUBROUTINE Matdef_INSTATIONARY_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Matdef_HEATALPHA_GenLinSc_Q1(lSc,idef,imat)
INTEGER :: idef,imat
TYPE(lScalarGen) lSc
INTEGER i,j,iFld
REAL*8 daux
EXTERNAL E011

 if (imat.eq.1) then
  DO ILEV=NLMIN,NLMAX

   CALL SETLEV(2)

   HeatDiffMat   => mg_HeatDiffMat(ILEV)%a
   IF (ASSOCIATED(AlphaDiffMat)) AlphaDiffMat  => mg_AlphaDiffMat(ILEV)%a
   LMassMat      => mg_LMassMat(ILEV)%a
   ConvectionMat =>  mg_ConvMat(ILEV)%a
   plMat  => mg_lMat(ILEV)

  ! Build up the matrix
   DO iFld=1,lSc%nOfFields
    IF (TRIM(lSc%Fld(iFld)%cName).eq.'temp') then
     DO I=1,plMat%nu
      J = plMat%LdA(I)
      daux = LMassMat(I) - thstep*(HeatDiffMat(J) + ConvectionMat(J))
      mg_AMat(iFld)%fld(ILEV)%a(J) = daux
      DO J=plMat%LdA(I)+1,plMat%LdA(I+1)-1
       daux = -thstep*(HeatDiffMat(J) + ConvectionMat(J))
       mg_AMat(iFld)%fld(ILEV)%a(J) =  daux
      END DO
     END DO
    ELSE
     DO I=1,plMat%nu
      J = plMat%LdA(I)
      daux = LMassMat(I) - thstep*(AlphaDiffMat(J) + ConvectionMat(J))
      mg_AMat(iFld)%fld(ILEV)%a(J) = daux
      DO J=plMat%LdA(I)+1,plMat%LdA(I+1)-1
       daux = -thstep*(AlphaDiffMat(J) + ConvectionMat(J))
       mg_AMat(iFld)%fld(ILEV)%a(J) =  daux
      END DO
     END DO
    END IF
   END DO
  END DO
 end if

 ILEV=NLMAX
 CALL SETLEV(2)

 HeatDiffMat   => mg_HeatDiffMat(ILEV)%a
 IF (ASSOCIATED(AlphaDiffMat)) AlphaDiffMat  => mg_AlphaDiffMat(ILEV)%a
 LMassMat      => mg_LMassMat(ILEV)%a
 ConvectionMat => mg_ConvMat(ILEV)%a
 plMat  => mg_lMat(ILEV)

 ! Build up the defect
 IF (idef.eq. 1) THEN
  DO iFld=1,lSc%nOfFields
  
   lSc%fld(iFld)%def = 0d0
   DO j=1,plMat%nu
    lSc%fld(iFld)%def(j) = lSc%fld(iFld)%def(j) + LMassMat(j)*lSc%fld(iFld)%val(j)
   END DO
   
   CALL LAX17(ConvectionMat,plMat%ColA,plMat%LdA,plMat%nu,&
   lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep,1d0)
   
   IF (TRIM(lSc%Fld(iFld)%cName).eq.'temp') then
    CALL LAX17(HeatDiffMat,plMat%ColA,plMat%LdA,plMat%nu,&
    lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep,1d0)
   END IF
   IF (TRIM(lSc%Fld(iFld)%cName).eq.'alpha1') then
    CALL LAX17(AlphaDiffMat,plMat%ColA,plMat%LdA,plMat%nu,&
    lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep,1d0)
   END IF

  END DO
  
 ELSE

  DO iFld=1,lSc%nOfFields
   CALL LAX17(mg_AMat(iFld)%fld(ILEV)%a,plMat%ColA,plMat%LdA,plMat%nu,&
   lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,-1d0,1d0)
  END DO
  
 END IF
 
 IF (idef.le.0) then
  DO iFld=1,lSc%nOfFields
   CALL DefTVD_LinScalar(lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,thstep)
  END DO
 end if


END SUBROUTINE Matdef_HEATALPHA_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Matdef_STATIONARY_GenLinSc_Q1(lSc,idef,imat)
INTEGER :: idef,imat
TYPE(lScalarGen) lSc
INTEGER i,j,iFld
REAL*8 daux
EXTERNAL E011

DO ILEV=NLMIN,NLMAX

 CALL SETLEV(2)

 LaplaceMat => mg_LaplaceMat(ILEV)%a
 ConvectionMat =>  mg_ConvMat(ILEV)%a
 plMat  => mg_lMat(ILEV)

! Build up the matrix
 DO iFld=1,lSc%nOfFields
  DO I=1,plMat%nu
   J = plMat%LdA(I)
   daux = LaplaceMat(J) + ConvectionMat(J)
   mg_AMat(iFld)%fld(ILEV)%a(J) = daux
   DO J=plMat%LdA(I)+1,plMat%LdA(I+1)-1
    daux = LaplaceMat(J) + ConvectionMat(J)
    mg_AMat(iFld)%fld(ILEV)%a(J) =  daux
   END DO
  END DO
 END DO
 
END DO

ILEV=NLMAX
CALL SETLEV(2)

 LaplaceMat => mg_LaplaceMat(ILEV)%a
 ConvectionMat =>  mg_ConvMat(ILEV)%a
 plMat  => mg_lMat(ILEV)

! Build up the defect
IF (idef.eq. 1) THEN
 DO iFld=1,lSc%nOfFields
  lSc%fld(iFld)%def = 0d0 
 END DO
ELSE
 DO iFld=1,lSc%nOfFields
  lSc%fld(iFld)%def = 0d0 
  CALL LAX17(mg_AMat(iFld)%fld(ILEV)%a,plMat%ColA,plMat%LdA,plMat%nu,&
  lSc%Fld(iFld)%val,lSc%Fld(iFld)%def,-1d0,1d0)
 END DO
END IF

END SUBROUTINE Matdef_STATIONARY_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_GenLinSc_Q1(lSc,&
    resScalar,defScalar,rhsScalar)
TYPE(lScalarGen), INTENT(INOUT) :: lSc
REAL*8  resScalar(*),defScalar(*),rhsScalar(*)
integer iFld

! CALL LL21 (myScalar%rhs,myScalar%ndof,RESF)
! RESF=MAX(1D-15,RESF)
! CALL LL21 (myScalar%def,myScalar%ndof,RESU)
! resScalar = RESU/RESF
! defScalar = RESU
! rhsScalar = RESF

do iFld = 1,lSc%nOfFields
 resScalar(iFld) = 0d0 
 CALL LL21 (lSc%fld(iFld)%def,lSc%ndof,defScalar(iFld))
!  write(*,*) myid,defScalar(iFld)
 rhsScalar(iFld) = 0d0 
end do

! pause

END SUBROUTINE Resdfk_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Solve_GenLinSc_Q1_MGLinScalar(lSc,Bndry_Mat,mfile)
use var_QuadScalar, only : myStat
USE mg_LinScalar, only:mg_solver
INTEGER mfile
TYPE(lScalarGen), INTENT(INOUT), TARGET :: lSc
REAL*8 daux,nrm_U,nrm_V,nrm_W
REAL*8, ALLOCATABLE :: dnrm(:)
INTEGER ndof,iFld
EXTERNAL Bndry_Mat

 if (.not.allocated(dnrm)) allocate(dnrm(lSc%nOfFields))

 CALL ZTIME(myStat%t0)

 IF (myid.ne.0) THEN
  DO ILEV = NLMIN,NLMAX
   DO iFld=1,lSc%nOfFields
    CALL Bndry_Mat(mg_AMat(iFld)%fld(ILEV)%a,mg_lMat(ILEV)%LdA,lSc%fld(iFld)%knpr,mg_lMat(ILEV)%nu)
    CALL E011_GenLinSc_Q1_UMAT(mg_AMat(iFld)%fld(ILEV)%a,mg_lMat(ILEV)%LdA,mg_lMat(ILEV)%nu,iFld)
   END DO
  END DO
  
  daux = 0d0
  DO iFld=1,lSc%nOfFields
   CALL LL21 (lSc%fld(iFld)%def,lSc%ndof,dnrm(iFld))
!    CALL LCL1 (lSc%fld(iFld)%val,lSc%ndof)
   daux = MAX(daux,dnrm(iFld))
  END DO
 END IF

CALL COMM_Maximum(daux)


! --------------- Set up the MG driver -----------------!
IF (myid.ne.0) THEN

!  IF (.not.allocated(MyMG%AXX)) ALLOCATE(MyMG%AXX(lSc%nOfFields))
 
 MyMG%AXX => mg_AMat
 DO iFld=1,lSc%nOfFields
  ndof = lSc%ndof
  lSc%knpr((iFld-1)*ndof+1:iFld*ndof) = lSc%fld(iFld)%knpr
 END DO
 
 MyMG%L    => mg_lMat
 MyMG%D    => lSc%def
 MyMG%AUX  => lSc%aux
 MyMG%KNPR => lSc%knpr
 
 MyMG%bProlRest => lSc%bProlRest
END IF

MyMG%cVariable          = ADJUSTL(TRIM(lSc%cName))
MyMG%MinIterCycle       = lSc%prm%MGprmIn%MinIterCycle
MyMG%MaxIterCycle       = lSc%prm%MGprmIn%MaxIterCycle
MyMG%nIterCoarse        = lSc%prm%MGprmIn%nIterCoarse
MyMG%DefImprCoarse      = lSc%prm%MGprmIn%DefImprCoarse
MyMG%nSmootherSteps     = lSc%prm%MGprmIn%nSmootherSteps
MyMG%CycleType          = lSc%prm%MGprmIn%CycleType
MyMG%Criterion1         = lSc%prm%MGprmIn%Criterion1
MyMG%Criterion2         = lSc%prm%MGprmIn%Criterion2*daux
MyMG%RLX                = lSc%prm%MGprmIn%RLX
MyMG%MinLev             = lSc%prm%MGprmIn%MinLev
MyMG%MedLev             = lSc%prm%MGprmIn%MedLev
MyMG%MaxDifLev          = lSc%prm%MGprmIn%MaxDifLev
MyMG%nOfFields          = lSc%nOfFields
MyMG%nOfSubsystemEqs    = lSc%ndof
daux = DBLE(NLMAX)
CALL COMM_Maximum(daux)
MyMG%MaxLev             = NINT(daux)

!-------------------  U - Component -------------------!
IF (myid.ne.0) THEN
 DO iFld=1,lSc%nOfFields
  ndof = lSc%ndof
  lSc%sol(NLMAX)%x((iFld-1)*ndof+1:(iFld)*ndof) = lSc%fld(iFld)%val(1:ndof)
  lSc%rhs(NLMAX)%x((iFld-1)*ndof+1:(iFld)*ndof) = lSc%fld(iFld)%aux(1:ndof)
  lSc%def(NLMAX)%x((iFld-1)*ndof+1:(iFld)*ndof) = 0d0
 END DO
 MyMG%X    => lSc%sol
 MyMG%B    => lSc%rhs
END IF

CALL MG_Solver(mfile,mterm)

IF (myid.ne.0) THEN
 ndof = lSc%ndof
 DO iFld=1,lSc%nOfFields
  lSc%fld(iFld)%val(1:ndof) = lSc%sol(NLMAX)%x((iFld-1)*ndof+1:(iFld)*ndof)
 END DO

 lSc%prm%MGprmOut(1)%UsedIterCycle = myMG%UsedIterCycle
 lSc%prm%MGprmOut(1)%nIterCoarse   = myMG%nIterCoarse
 lSc%prm%MGprmOut(1)%DefInitial    = myMG%DefInitial
 lSc%prm%MGprmOut(1)%DefFinal      = myMG%DefFinal
 lSc%prm%MGprmOut(1)%RhoMG1        = myMG%RhoMG1
 lSc%prm%MGprmOut(1)%RhoMG2        = myMG%RhoMG2

! Update the solution
 DO iFld=1,lSc%nOfFields
  CALL LLC1(lSc%fld(iFld)%val_old,lSc%fld(iFld)%val,&
       lSc%ndof,1D0,1D0)
 END DO

END IF

CALL ZTIME(myStat%t1)
myStat%tMGUVW = myStat%tMGUVW + (myStat%t1-myStat%t0)
myStat%iLinUVW = myStat%iLinUVW + lSc%prm%MGprmOut(1)%UsedIterCycle &
                                + lSc%prm%MGprmOut(2)%UsedIterCycle &
                                + lSc%prm%MGprmOut(3)%UsedIterCycle

 ! write(*,*) 'reaching MG ...', daux
 ! pause

END SUBROUTINE Solve_GenLinSc_Q1_MGLinScalar
!
! ----------------------------------------------
!
SUBROUTINE Protocol_GenLinSc_Q1(mfile,lSc,nINL,&
           ResScalar,DefScalar,RhsScalar,cTitle)
TYPE(lScalarGen), INTENT(INOUT) :: lSc
INTEGER nINL,mfile
INTEGER i,length,iFld
REAL*8 ResScalar(*),DefScalar(*),RhsScalar(*)
CHARACTER, OPTIONAL:: cTitle*(*)

IF (myid.eq.showID) THEN
length =  LEN(lSc%cName)

IF (PRESENT(cTitle)) THEN
 length = LEN(cTitle)
 IF (MOD(length,2).eq.1) length = length + 1
 length = (104-length)/2
END IF

IF (nINL.EQ.0) THEN
 IF (PRESENT(cTitle)) THEN
  WRITE(*,4)
  WRITE(*,*) cTitle
!  WRITE(MFILE,4) cTitle
 ELSE
  WRITE(*,5)
!  WRITE(MFILE,5)
 END IF

 WRITE(MTERM,'(A8$)') "INL"
 WRITE(MFILE,'(A8$)') "INL"
 DO iFld=1,lSc%nOfFields
  WRITE(MTERM,'(2X,A14$)') TRIM(lSc%fld(iFld)%cName)
  WRITE(MFILE,'(2X,A14$)') TRIM(lSc%fld(iFld)%cName)
 END DO
 WRITE(MTERM,*)
 WRITE(MFILE,*)
 
 WRITE(MTERM,5)
 WRITE(MFILE,5)

 WRITE(MTERM,'(A8$)') "Criteria"
 WRITE(MFILE,'(A8$)') "Criteria"
 DO iFld=1,lSc%nOfFields
  WRITE(MTERM,'(6X,ES10.3$)') RhsScalar(iFld)
  WRITE(MFILE,'(6X,ES10.3$)') RhsScalar(iFld)
 END DO
 WRITE(MTERM,*)
 WRITE(MFILE,*)

 WRITE(MTERM,5)
 WRITE(MFILE,5)
 
 WRITE(MTERM,'(I8$)') 0
 WRITE(MFILE,'(I8$)') 0
 DO iFld=1,lSc%nOfFields
  WRITE(MTERM,'(6X,ES10.3$)') DefScalar(iFld)
  WRITE(MFILE,'(6X,ES10.3$)') DefScalar(iFld)
 END DO
 WRITE(MTERM,*)
 WRITE(MFILE,*)
 
ELSE

 WRITE(MTERM,'(I8$)') nINL
 WRITE(MFILE,'(I8$)') nINL
 DO iFld=1,lSc%nOfFields
  WRITE(MTERM,'(6X,ES10.3$)') DefScalar(iFld)
  WRITE(MFILE,'(6X,ES10.3$)') DefScalar(iFld)
 END DO
 WRITE(MTERM,'(2I5,2XES10.3)') &
 lSc%prm%MGprmOut(1)%UsedIterCycle,lSc%prm%MGprmOut(1)%nIterCoarse,lSc%prm%MGprmOut(1)%RhoMG1
 WRITE(MFILE,'(2I5,2XES10.3)') &
 lSc%prm%MGprmOut(1)%UsedIterCycle,lSc%prm%MGprmOut(1)%nIterCoarse,lSc%prm%MGprmOut(1)%RhoMG1
END IF

END IF

5  FORMAT(104('-'))
4  FORMAT(104('-'))


END SUBROUTINE Protocol_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE Get_GenLinSc_Q1_Parameters(myParam,cName,mfile)
IMPLICIT NONE
TYPE(tParam) myParam
INTEGER mfile
CHARACTER*(*) cName
integer :: iEnd,iAt,iEq
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile = 990076
integer :: istat,iFld

OPEN (UNIT=myFile,FILE=TRIM(ADJUSTL(myDataFile)),action='read',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open data file: ",myDataFile  
  stop          
end if

IF (myid.eq.showid) WRITE(mfile,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))
IF (myid.eq.showid) WRITE(mterm,'(47("-"),A10,47("-"))') TRIM(ADJUSTL(cName))

myParam%MGprmIn%MaxLev = NLMAX + 1 
myParam%MGprmIn%MaxDifLev = myParam%MGprmIn%MaxLev - NLMAX

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

    CASE ("Equation")
    READ(string(iEq+1:),*) param
    myParam%cEquation = TRIM(ADJUSTL(param))
    call inip_toupper_replace(myParam%cEquation)
    IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%cEquation
    IF (myid.eq.showid) write(mfile,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%cEquation

    CASE ("nOfFields")
    READ(string(iEq+1:),*) myParam%nOfFields
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%nOfFields
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%nOfFields
    ALLOCATE(myParam%cField(myParam%nOfFields))

    CASE ("Fields")
    READ(string(iEq+1:),*) myParam%cField
!    myParam%cField = TRIM(ADJUSTL(param))
    IF (myid.eq.showid) write(mterm,'(A$)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= "
    IF (myid.eq.showid) write(mfile,'(A$)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= "
    do ifld=1,myParam%nOfFields-1
     IF (myid.eq.showid) write(mterm,'(A$)') TRIM(myParam%cField(iFld))//","
     IF (myid.eq.showid) write(mfile,'(A$)') TRIM(myParam%cField(iFld))//","
    END DO
    IF (myid.eq.showid) write(mterm,'(A)') TRIM(myParam%cField(myParam%nOfFields))
    IF (myid.eq.showid) write(mfile,'(A)') TRIM(myParam%cField(myParam%nOfFields))
    
    CASE ("epsCrit")
    READ(string(iEq+1:),*) myParam%epsCrit
!    IF (myid.eq.showid) write(mterm,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%epsCrit
!    IF (myid.eq.showid) write(mfile,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%epsCrit
    CASE ("defCrit")
    READ(string(iEq+1:),*) myParam%defCrit
!    IF (myid.eq.showid) write(mterm,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%defCrit
!    IF (myid.eq.showid) write(mfile,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%defCrit
    CASE ("MinDef")
    READ(string(iEq+1:),*) myParam%MinDef
!    IF (myid.eq.showid) write(mterm,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MinDef
!    IF (myid.eq.showid) write(mfile,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MinDef
!     CASE ("SolvType")
!     READ(string(iEq+1:),*) myParam%SolvType
!     IF (myid.eq.showid) write(mterm,'(A,I)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%SolvType
!     IF (myid.eq.showid) write(mfile,'(A,I)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%SolvType
    CASE ("NLmin")
    READ(string(iEq+1:),*) myParam%NLmin
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%NLmin
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%NLmin
    CASE ("NLmax")
    READ(string(iEq+1:),*) myParam%NLmax
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%NLmax
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%NLmax
    CASE ("MGMinLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinLev
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MinLev
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MinLev
    CASE ("MGMedLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MedLev
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MedLev
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MedLev
    CASE ("MGMaxLev")
    READ(string(iEq+1:),*) myParam%MGprmIn%MaxLev
    IF (myid.eq.showid) write(mterm,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MaxLev
    IF (myid.eq.showid) write(mfile,'(A,I0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MaxLev
    myParam%MGprmIn%MaxDifLev = myParam%MGprmIn%MaxLev - NLMAX
    CASE ("MGMinIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MinIterCycle
    IF (myid.eq.showid) write(mterm,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MinIterCycle
    IF (myid.eq.showid) write(mfile,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MinIterCycle
    CASE ("MGMaxIterCyc")
    READ(string(iEq+1:),*) myParam%MGprmIn%MaxIterCycle
    IF (myid.eq.showid) write(mterm,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MaxIterCycle
    IF (myid.eq.showid) write(mfile,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%MaxIterCycle
    CASE ("MGSmoothSteps")
    READ(string(iEq+1:),*) myParam%MGprmIn%nSmootherSteps
    IF (myid.eq.showid) write(mterm,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%nSmootherSteps
    IF (myid.eq.showid) write(mfile,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%nSmootherSteps
    CASE ("MGIterCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%nIterCoarse
    IF (myid.eq.showid) write(mterm,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%nIterCoarse
    IF (myid.eq.showid) write(mfile,'(A,I0)') &
      TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%nIterCoarse
    CASE ("MGDefImprCoarse")
    READ(string(iEq+1:),*) myParam%MGprmIn%DefImprCoarse
!    IF (myid.eq.showid) write(mterm,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%DefImprCoarse
!    IF (myid.eq.showid) write(mfile,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%DefImprCoarse
    CASE ("MGCriterion1")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion1
!    IF (myid.eq.showid) write(mterm,'(A,E0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%Criterion1
!    IF (myid.eq.showid) write(mfile,'(A,E0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%Criterion1
    CASE ("MGCriterion2")
    READ(string(iEq+1:),*) myParam%MGprmIn%Criterion2
!    IF (myid.eq.showid) write(mterm,'(A,E0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%Criterion2
!    IF (myid.eq.showid) write(mfile,'(A,E0)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%Criterion2
    CASE ("MGCycType")
    READ(string(iEq+1:),*) param
    myParam%MGprmIn%CycleType = TRIM(ADJUSTL(param))
    IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType
    IF (myid.eq.showid) write(mfile,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CycleType
!     CASE ("MGCrsSolverType")
!     READ(string(iEq+1:),*) myParam%MGprmIn%CrsSolverType
!     IF (myid.eq.showid) write(mterm,'(A,I)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CrsSolverType
!     IF (myid.eq.showid) write(mfile,'(A,I)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%CrsSolverType
    CASE ("MGRelaxPrm")
    READ(string(iEq+1:),*) myParam%MGprmIn%RLX
!    IF (myid.eq.showid) write(mterm,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%RLX
!    IF (myid.eq.showid) write(mfile,'(A,E)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%MGprmIn%RLX
    
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

END SUBROUTINE Get_GenLinSc_Q1_Parameters
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_AFCStruct()
INTEGER I,ILOC

ILEV = NLMAX

plMat  => mg_lMat(ILEV)

AFC%nedge = (plMat%na-plMat%nu)/2
AFC%nu = plMat%nu

ALLOCATE(AFC%iaux(AFC%nu))
ALLOCATE(AFC%isep(AFC%nu))

ALLOCATE(AFC%inod (AFC%nedge))
ALLOCATE(AFC%jnod (AFC%nedge))
ALLOCATE(AFC%aedge(AFC%nedge))

ALLOCATE(AFC%pp(AFC%nu))
ALLOCATE(AFC%pm(AFC%nu))
ALLOCATE(AFC%qp(AFC%nu))
ALLOCATE(AFC%qm(AFC%nu))

DO I=1,AFC%nu
  ILOC=plMat%LdA(I)
1 ILOC=ILOC+1 
  if(iloc .le. plMat%na)then
    IF (plMat%ColA(ILOC).LT.I.AND.ILOC.LT.plMat%LdA(I+1)) GOTO 1
    AFC%isep(I)=ILOC-1
  end if
END DO

END SUBROUTINE Create_GenLinSc_Q1_AFCStruct
!
! ----------------------------------------------
!
SUBROUTINE InitAFC_GenLinSc_Q1()

 ILEV = NLMAX
 
 ConvectionMat =>  mg_ConvMat(ILEV)%a
 plMat  => mg_lMat(ILEV)

 CALL AFC_LinScalar(ConvectionMat,plMat%ColA,plMat%LdA,plMat%nu,&
      AFC%isep,AFC%iaux,AFC%inod,AFC%jnod,AFC%aedge)
      
END SUBROUTINE InitAFC_GenLinSc_Q1
!
! ----------------------------------------------
!
SUBROUTINE ProlongateSolution_GenLinSc_Q1_sub(lSc)

TYPE(lScalarGen), INTENT(INOUT), TARGET :: lSc
INTEGER J,K,ndof,iFld

IF (myid.ne.0) THEN

 NLMAX = NLMAX + 1 
 
 ILEV=NLMAX
 CALL SETLEV(2)

 J = NLMAX
 K = NLMAX-1
 mgLev = J

 IF(myid.eq.showid) WRITE(*,'(A,I0,A,I0)') "Prolongation of GenLinSc-components from L",K," to L",J

 MyMG%cVariable = "Temper"
 myMG%nOfFields = lSc%nOfFields

 IF (myid.ne.0) THEN
  DO iFld=1,lSc%nOfFields
   ndof = KNVT(K)
   lSc%sol(K)%x((iFld-1)*ndof+1:(iFld)*ndof) = lSc%fld(iFld)%val(1:ndof)
   lSc%knpr((iFld-1)*ndof+1:iFld*ndof) = 0
  END DO
 END IF

 MyMG%X    => lSc%sol
 MyMG%AUX  => lSc%aux
 MyMG%KNPR => lSc%knpr

 CALL mgProlongation()

 do iFld=1,lSc%nOfFields
  ndof = KNVT(J)
  lSc%Fld(iFld)%Val(1:ndof) = lSc%aux(J)%x((iFld-1)*ndof+1:(iFld)*ndof)
 end do
 
 NLMAX = NLMAX - 1 
 
END IF

END SUBROUTINE ProlongateSolution_GenLinSc_Q1_sub
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_DiffCoeff(lSc)
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  use Sigma_User, only: mySegmentIndicator,myMultiMat
  EXTERNAL E011
  TYPE(lScalarGen) :: lSc
  integer :: i,iel,ivt,iSeg,ndof,iMat,iFld,iMinMat
  REAL*8 :: daux,Cp,Rho,Lambda,Temp,dMinField

  if (.not.allocated(mg_DiffCoeff)) allocate(mg_DiffCoeff(NLMIN:NLMAX))

  do ilev=NLMAX,NLMIN,-1
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_DiffCoeff(ILEV)%x))then
      ALLOCATE(mg_DiffCoeff(ILEV)%x(mg_mesh%level(ilev)%nel))
    end if

    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMAX) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [Dc]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF

    ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net + mg_mesh%level(ilev)%nat
    mg_DiffCoeff(ILEV)%x = 0d0
    
    if (ilev.eq.nlmax) then
     DO iel=1,mg_mesh%level(ilev)%nel
      do i=1,8
       ivt = mg_mesh%level(ilev)%kvert(i,iel)
       iSeg = INT(mySegmentIndicator(2,ivt))
       IF (iSeg.eq.0) THEN
        daux = 0.188d0 ! cm2/s !Steel
       ELSE
!        iMat = lSc%nOfFields
        dMinField = -1d0
        iMat = 0
        do iFld=2,lSc%nOfFields
         if (lSc%Fld(iFld)%Val(ivt).gt.0.5d0) then
          iMat = iFld - 1
         else
          if (dMinField.lt.lSc%Fld(iFld)%Val(ivt)) then
           dMinField = lSc%Fld(iFld)%Val(ivt)
           iMinMat = iFld - 1
          end if
         end if
        end do
        if (iMat.eq.0) iMat = iMinMat
        Temp   = lSc%Fld(1)%Val(ivt)
        Lambda = 1e0*(myMultiMat%Mat(iMat)%Thermodyn%lambda  + Temp*myMultiMat%Mat(iMat)%Thermodyn%LambdaSteig  ) ! W/m/K
        Cp     = 1e3*(myMultiMat%Mat(iMat)%Thermodyn%Cp      + Temp*myMultiMat%Mat(iMat)%Thermodyn%CpSteig      ) ! J/kg/K
        Rho    = 1e3*(myMultiMat%Mat(iMat)%Thermodyn%density - Temp*myMultiMat%Mat(iMat)%Thermodyn%densitySteig ) ! kg/m3
        daux   = 1e4*Lambda/(Rho*Cp) ! cm2/s
       END IF
       mg_DiffCoeff(ILEV)%x(iel) = mg_DiffCoeff(ILEV)%x(iel) + 0.125d0*daux
      end do
     end do
     
    else
    
     DO iel=1,mg_mesh%level(ilev)%nel
      mg_DiffCoeff(ILEV)%x(iel) = 0.125d0*mg_DiffCoeff(ILEV+1)%x(iel)
      do i=1,7
       ivt = mg_mesh%level(ilev)%nel + (iel-1)*7 + 1
       daux = mg_DiffCoeff(ILEV+1)%x(ivt)
       mg_DiffCoeff(ILEV)%x(iel) = mg_DiffCoeff(ILEV)%x(iel) + 0.125d0*daux
      end do
     end do
    end if
    
  end do ! ilev

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

END SUBROUTINE Create_GenLinSc_Q1_DiffCoeff
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_Heat_Diffusion()
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  EXTERNAL E011
  integer :: i

  IF (.not.allocated(mg_HeatDiffMat)) ALLOCATE(mg_HeatDiffMat(NLMIN:NLMAX))

  do ilev=NLMIN,NLMAX
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_HeatDiffMat(ILEV)%a))then
      ALLOCATE(mg_HeatDiffMat(ILEV)%a(mg_lMat(ILEV)%na))
    end if

    HeatDiffMat=>mg_HeatDiffMat(ILEV)%a
    plMat=>mg_lMat(ILEV)

    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMIN) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [D]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF

    HeatDiffMat=0d0

    CALL DIE_DiffMatQ1(HeatDiffMat,mg_DiffCoeff(ILEV)%x,&
    plMat%nu,plMat%ColA,plMat%LdA,&
    mg_mesh%level(ilev)%kvert,&
    mg_mesh%level(ilev)%karea,&
    mg_mesh%level(ilev)%kedge,&
    mg_mesh%level(ilev)%dcorvg,&
    E011)
    
  end do

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

END SUBROUTINE Create_GenLinSc_Q1_Heat_Diffusion
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_Alpha_Diffusion(dAlpha)
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  EXTERNAL E011
  REAL*8 dAlpha
  integer :: i

  IF (.not.allocated(mg_AlphaDiffMat)) ALLOCATE(mg_AlphaDiffMat(NLMIN:NLMAX))

  do ilev=NLMIN,NLMAX
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_AlphaDiffMat(ILEV)%a))then
      ALLOCATE(mg_AlphaDiffMat(ILEV)%a(mg_lMat(ILEV)%na))
    end if

    AlphaDiffMat=>mg_AlphaDiffMat(ILEV)%a
    plMat=>mg_lMat(ILEV)

    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMIN) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [D]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF

    AlphaDiffMat=0d0

    CALL CnstDiffMatQ1(AlphaDiffMat,&
      plMat%nu,plMat%ColA,plMat%LdA,& 
      mg_mesh%level(ilev)%kvert,&
      mg_mesh%level(ilev)%karea,&
      mg_mesh%level(ilev)%kedge,&
      mg_mesh%level(ilev)%dcorvg,&
      dAlpha,E011)

  end do

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

END SUBROUTINE Create_GenLinSc_Q1_Alpha_Diffusion
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_Diffusion(dAlpha)
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  EXTERNAL E011
  REAL*8 dAlpha
  integer :: i

  IF (.not.allocated(mg_LaplaceMat)) ALLOCATE(mg_LaplaceMat(NLMIN:NLMAX))

  do ilev=NLMIN,NLMAX
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_LaplaceMat(ILEV)%a))then
      ALLOCATE(mg_LaplaceMat(ILEV)%a(mg_lMat(ILEV)%na))
    end if

    LaplaceMat=>mg_LaplaceMat(ILEV)%a
    plMat=>mg_lMat(ILEV)

    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMIN) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [D]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF

    LaplaceMat=0d0

    CALL CnstDiffMatQ1(LaplaceMat,&
      plMat%nu,plMat%ColA,plMat%LdA,& 
      mg_mesh%level(ilev)%kvert,&
      mg_mesh%level(ilev)%karea,&
      mg_mesh%level(ilev)%kedge,&
      mg_mesh%level(ilev)%dcorvg,&
      dAlpha,E011)

  end do

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

END SUBROUTINE Create_GenLinSc_Q1_Diffusion
!
! ----------------------------------------------
!
SUBROUTINE Create_GenLinSc_Q1_Mass(dAlpha)
  USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
  EXTERNAL E011
  REAL*8 dAlpha, DML
  integer :: i,J

  IF (.not.allocated(mg_MassMat)) ALLOCATE(mg_MassMat(NLMIN:NLMAX))

  do ilev=NLMIN,NLMAX
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_MassMat(ILEV)%a))then
      ALLOCATE(mg_MassMat(ILEV)%a(mg_lMat(ILEV)%na))
    end if

    MassMat=>mg_MassMat(ILEV)%a
    plMat=>mg_lMat(ILEV)

    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMIN) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [M]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF

    MassMat=0d0

    CALL CnstMassMatQ1(MassMat,&
      plMat%nu,plMat%ColA,plMat%LdA,& 
      mg_mesh%level(ilev)%kvert,&
      mg_mesh%level(ilev)%karea,&
      mg_mesh%level(ilev)%kedge,&
      mg_mesh%level(ilev)%dcorvg,&
      dAlpha,E011)

  end do

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

  IF (.not.allocated(mg_LMassMat)) ALLOCATE(mg_LMassMat(NLMIN:NLMAX))

  do ilev=NLMIN,NLMAX
  
    CALL SETLEV(2)

    IF (.NOT.ALLOCATED(mg_LMassMat(ILEV)%a))then
      ALLOCATE(mg_LMassMat(ILEV)%a(mg_lMat(ILEV)%nu))
    end if

    MassMat    => mg_MassMat(ILEV)%a
    LMassMat   => mg_LMassMat(ILEV)%a
    plMat      => mg_lMat(ILEV)
    
    IF (myid.eq.showID) THEN
     IF (ILEV.EQ.NLMIN) THEN
      WRITE(MTERM,'(A,I1,A)', advance='no') " [ML]: [", ILEV,"]"
     ELSE
      WRITE(MTERM,'(A,I1,A)', advance='no') ", [",ILEV,"]"
     END IF
    END IF
    
    DO I=1,plMat%nu
     DML = 0d0
     DO J=plMat%LdA(I),plMat%LdA(I+1)-1
      DML = DML + MassMat(J)
     END DO
     LMassMat(I) = DML
    END DO
    
  END DO

  IF (myid.eq.showID) WRITE(MTERM,'(A)', advance='yes') " |"

END SUBROUTINE Create_GenLinSc_Q1_Mass
