MODULE def_LinScalar

USE def_FEAT
USE PP3D_MPI, ONLY:E011Sum,E011Mat,myid,showID,MGE011,comm_maximum
USE var_QuadScalar, ONLY: myMG,myDataFile,&
                          mg_mesh, mg_Matrix, TMatrix, mg_dVector,&
                          tMGParamIn, tMGParamOut
USE mg_LinScalar, only:mg_solver


IMPLICIT NONE

TYPE tParam
 Character*20 :: cEquation

 REAL*8  defCrit,epsCrit,MinDef

 INTEGER NLmin,NLmax

 INTEGER SolvIter,SolvType,iMass

 LOGICAL AFC
 TYPE(tMGParamIn) :: MGprmIn
 TYPE(tMGParamOut):: MGprmOut(3)
END TYPE tParam

TYPE mg_vector
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: x
END TYPE mg_vector

TYPE lScalar3
 CHARACTER cName*7
 INTEGER :: ndof,na
 LOGICAL :: bProlRest=.FALSE. 

 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knprX,knprY,knprZ
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valX_old,valY_old,valZ_old
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: valX,valY,valZ
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: defX,defY,defZ
 TYPE(mg_dVector), DIMENSION(:),ALLOCATABLE :: def,sol,aux,rhs
 TYPE(tParam) :: prm
 
END TYPE

TYPE lScalar
 CHARACTER cName*7
 INTEGER :: ndof,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: knpr
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aux,rhs,def,val_old
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: src,snk
 TYPE(mg_vector), DIMENSION(:),ALLOCATABLE :: val
 TYPE(tParam) :: prm
END TYPE

TYPE tAFC
 INTEGER :: iedge,nedge,nu
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: inod,jnod,iaux,isep
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: aedge
 REAL*8  , DIMENSION(:)  , ALLOCATABLE :: pp,pm,qm,qp
END TYPE
TYPE(tAFC) AFC

TYPE TlMatrix
 INTEGER :: nu,na
 INTEGER , DIMENSION(:)  , ALLOCATABLE :: ColA,LdA
END TYPE

TYPE(TlMatrix) :: lMat

!---------------------------------------------------------------------

REAL*8  , DIMENSION(:)  , ALLOCATABLE :: Mmat,MLmat,Kmat,Dmat
REAL*4  , DIMENSION(:)  , ALLOCATABLE :: Amat
REAL*4  , DIMENSION(:)  , ALLOCATABLE :: AmatX,AmatY,AmatZ

!---------------------------------------------------------------------

! A general Matrix structure for E011
TYPE(TMatrix), DIMENSION(:), ALLOCATABLE, TARGET :: mg_lMat

! A pointer for E011 matrix level to be used as a shorthand 
TYPE(TMatrix), POINTER :: plMat

!---------------------------------------------------------------------

! Pointer to the entries of the multi-level matrices 
REAL*8  , DIMENSION(:)  , POINTER :: LaplaceMat, A11mat, A22mat, A33mat

TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_A11mat,mg_A22mat,mg_A33mat

! A multi-level Matrix entries struct 
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_LaplaceMat

TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_E011Prol,mg_E011Rest
TYPE(TMatrix)   , DIMENSION(:)  , ALLOCATABLE,  TARGET :: mg_E011ProlM,mg_E011RestM


CONTAINS
!
! ----------------------------------------------
!
SUBROUTINE Create_AFCStruct()
INTEGER I,ILOC

AFC%nedge = (lMat%na-lMat%nu)/2
AFC%nu = lMat%nu

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
  ILOC=lMat%LdA(I)
1 ILOC=ILOC+1 
  IF (lMat%ColA(ILOC).LT.I.AND.ILOC.LT.lMat%LdA(I+1)) GOTO 1
  AFC%isep(I)=ILOC-1
END DO

END SUBROUTINE Create_AFCStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_MatStruct()
INTEGER iSymm,nERow
INTEGER , DIMENSION(:)  , ALLOCATABLE :: TempColA
INTEGER I,J

EXTERNAL E011,coefst

ILEV=NLMAX
CALL SETLEV(2)

ALLOCATE(TempColA(100*mg_mesh%level(ilev)%nvt))
ALLOCATE(lMat%LdA(mg_mesh%level(ilev)%nvt+1))
lMat%nu = mg_mesh%level(ilev)%nvt
lMat%na = 100*mg_mesh%level(ilev)%nvt
iSymm =   0
nERow =   100

CALL AP7(TempColA,lMat%LdA,lMat%na,lMat%nu,E011,iSymm,nERow,&
         mg_mesh%level(ilev)%kvert,&
         mg_mesh%level(ilev)%kedge,&
         mg_mesh%level(ilev)%karea)

ALLOCATE(lMat%ColA(lMat%na))
lMat%ColA(:) = TempColA(1:lMat%na)
DEALLOCATE(TempColA)

! if (myid.eq.1) OPEN(987,FILE='#data/MAT1.txt')
! if (myid.eq.2) OPEN(987,FILE='#data/MAT2.txt')
! if (myid.eq.3) OPEN(987,FILE='#data/MAT3.txt')
! 
! DO I=1,NVT
! WRITE(987,*) "---",lMat%LdA(I)
!  DO J=lMat%LdA(I),lMat%LdA(I+1)-1
!   WRITE(987,*) I,J,lMat%ColA(J)
!  END DO
! END DO
! CLOSE(987)

END SUBROUTINE Create_MatStruct
!
! ----------------------------------------------
!
SUBROUTINE Create_MassMat()
INTEGER nbloc
PARAMETER (NBLOC = 1)
INTEGER, DIMENSION(2,nnab,NBLOC) :: kabst
INTEGER, DIMENSION(NNDER,NNDER,NBLOC) :: coecon
INTEGER, DIMENSION(NBLOC) :: koff
LOGICAL :: bcon
INTEGER icub,lint,kabstn
INTEGER iSymm,nERow
INTEGER I,J

EXTERNAL E011,coefst

ILEV=NLMAX
CALL SETLEV(2)

 bcon=.TRUE.
 kabst (1,1,1)=1
 kabst (2,1,1)=1
 kabstn       =1
 icub  = 7
 lint = 0
 coecon = 0
 koff = 0

ALLOCATE(Mmat(lMat%na))

CALL AB07(Mmat,lMat%ColA,lMat%LdA,lMat%na,lMat%nu,nbloc,&
          koff,&
          mg_mesh%level(ilev)%kvert,&
          mg_mesh%level(ilev)%kedge,&
          mg_mesh%level(ilev)%karea,&
          mg_mesh%level(ilev)%dcorvg,&
          E011,coefst,&
          bcon,coecon,kabst,kabstn,icub,ISymm,lint)


! if (myid.eq.1) OPEN(987,FILE='#data/MMAT1.txt')
! if (myid.eq.2) OPEN(987,FILE='#data/MMAT2.txt')
! if (myid.eq.3) OPEN(987,FILE='#data/MMAT3.txt')
! 
! DO I=1,NVT
! WRITE(987,*) "---",lMat%LdA(I)
!  DO J=lMat%LdA(I),lMat%LdA(I+1)-1
!   WRITE(987,*) I,J,mMat(J)
!  END DO
! END DO
! CLOSE(987)

END SUBROUTINE Create_MassMat
!
! ----------------------------------------------
!
SUBROUTINE Create_LMassMat()
INTEGER I,J
REAL*8 DML

ILEV=NLMAX
CALL SETLEV(2)

ALLOCATE(MLmat(lMat%nu))

DO I=1,lMat%nu
 DML = 0d0
 DO J=lMat%LdA(I),lMat%LdA(I+1)-1
  DML = DML + Mmat(J)
 END DO
 MLmat(I) = DML
END DO

END SUBROUTINE Create_LMassMat
!
! ----------------------------------------------
!
SUBROUTINE Create_LKonvMat()

ALLOCATE(Kmat(lMat%na))

END SUBROUTINE Create_LKonvMat
!
! ----------------------------------------------
!
SUBROUTINE Create_AMat()

ALLOCATE(Amat(lMat%na))

END SUBROUTINE Create_AMat
!
! ----------------------------------------------
!
SUBROUTINE Create_DiffMat(Alpha)
EXTERNAL E011
REAL*8 Alpha(*)

ALLOCATE(DMat(lMat%na))
DMat=0d0

ILEV=NLMAX
CALL SETLEV(2)

CALL DiffMatQ1(Alpha,DMat,lMat%nu,lMat%ColA,lMat%LdA,&
               mg_mesh%level(ilev)%kvert,&
               mg_mesh%level(ilev)%karea,&
               mg_mesh%level(ilev)%kedge,&
               mg_mesh%level(ilev)%dcorvg,&
               E011)

END SUBROUTINE Create_DiffMat
!
! ----------------------------------------------
!
SUBROUTINE Initialize(myScalar)
TYPE(lScalar) myScalar

ILEV=NLMAX
CALL SETLEV(2)

myScalar%ndof = mg_mesh%level(ilev)%nvt
myScalar%na = lMat%na

ALLOCATE(myScalar%knpr(myScalar%ndof))
ALLOCATE(myScalar%aux(myScalar%ndof))
ALLOCATE(myScalar%rhs(myScalar%ndof))
ALLOCATE(myScalar%def(myScalar%ndof))
ALLOCATE(myScalar%val_old(myScalar%ndof))

ALLOCATE(myScalar%val(NLMIN:NLMAX))
DO ILEV=NLMIN,NLMAX
 ALLOCATE(myScalar%val(ILEV)%x(mg_mesh%level(ilev)%nvt))
END DO
! ALLOCATE(myScalar%src(myScalar%ndof))
! ALLOCATE(myScalar%snk(myScalar%ndof))

END SUBROUTINE Initialize
!
! ----------------------------------------------
!
SUBROUTINE Create_Knpr(myScalar_Knpr)
EXTERNAL myScalar_Knpr

ILEV=NLMAX
CALL SETLEV(2)

CALL myScalar_Knpr(mg_mesh%level(ilev)%dcorvg)

END SUBROUTINE Create_Knpr
!
! ----------------------------------------------
!
SUBROUTINE InitCond(myScalar_InitCond)
EXTERNAL myScalar_InitCond

ILEV=NLMAX
CALL SETLEV(2)

CALL myScalar_InitCond(mg_mesh%level(ilev)%dcorvg)

END SUBROUTINE InitCond
!
! ----------------------------------------------
!
SUBROUTINE Matdef_general_LinScalar(myScalar,idef,imat)
INTEGER :: idef,imat
TYPE(lScalar) myScalar

 ! Build up the matrix
 IF (imat.eq.1) THEN
   Amat=REAL(-thstep*(Kmat+DMat))
   Amat(lMat%LdA(1:lMat%nu))=Amat(lMat%LdA(1:lMat%nu))+REAL(MLmat)
 END IF

 ! Build up the defect
 IF (idef.eq. 1) THEN
  IF (myScalar%prm%iMAss.EQ.1) THEN
     myScalar%def = MLmat*myScalar%val(NLMAX)%x
     CALL LAX17(Kmat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%val(NLMAX)%x,myScalar%def,thstep,1d0)
     CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%val(NLMAX)%x,myScalar%def,thstep,1d0)
  END IF
  IF (myScalar%prm%iMAss.EQ.2) THEN
     CALL LAX17(Mmat+thstep*(Kmat+Dmat),lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%val(NLMAX)%x,myScalar%def,1d0,0d0)
  END IF
 ELSE
  IF (myScalar%prm%iMAss.EQ.1) THEN
     CALL LAX37(Amat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%val(NLMAX)%x,myScalar%def,-1d0,1d0)
  END IF
  IF (myScalar%prm%iMAss.EQ.2) THEN
     CALL LAX17(-Mmat+thstep*(Kmat+Dmat),lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%val(NLMAX)%x,myScalar%def,1d0,1d0)
  END IF
 END IF

 ! Perform Algebraic Flux Correction (if needed)
 IF (myScalar%prm%AFC) THEN
   CALL DefTVD_LinScalar(myScalar%val(NLMAX)%x,myScalar%def,THSTEP)
 END IF

END SUBROUTINE Matdef_general_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_General_LinScalar(myScalar,&
           resScalar,defScalar,rhsScalar)
TYPE(lScalar), INTENT(INOUT) :: myScalar
REAL*8  resScalar,defScalar,rhsScalar,RESF,RESU

CALL LL21 (myScalar%rhs,myScalar%ndof,RESF)
RESF=MAX(1D-15,RESF)

CALL LL21 (myScalar%def,myScalar%ndof,RESU)

resScalar = RESU/RESF
defScalar = RESU
rhsScalar = RESF
!write(*,*) RESU,RESF

END SUBROUTINE Resdfk_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Solve_General_LinScalar(myScalar,knpr,Bndry_Val,Bndry_Mat)
INTEGER KNPR(*)
INTEGER iLinIter
REAL*8 DefInit,DefCurrent
TYPE(lScalar), INTENT(INOUT) :: myScalar
EXTERNAL Bndry_Val,Bndry_Mat

IF (myid.ne.0) THEN
 CALL Bndry_Mat(Amat,lMat%LdA,myScalar%knpr)

 CALL E011Mat(Amat,lMat%LdA,lMat%nu)

 CALL LCL1 (myScalar%val(NLMAX)%x,myScalar%ndof)
END IF

CALL GetDefNorm(Amat,lMat%ColA,lMat%LdA,myScalar%val(NLMAX)%x,&
     myScalar%def,myScalar%aux,myScalar%ndof,DefInit)

DO iLinIter=1,5
 IF (myid.ne.0) THEN
  IF (myScalar%prm%SolvType.EQ.1) THEN
   CALL SSORSolver(Amat,lMat%ColA,lMat%LdA,&
        myScalar%val(NLMAX)%x,myScalar%def,myScalar%aux,&
        KNPR,myScalar%ndof,1*myScalar%prm%SolvIter,0.7d0)
  ENDIF
  IF (myScalar%prm%SolvType.EQ.2) THEN
   CALL JacobiSolver(Amat,lMat%ColA,lMat%LdA,&
        myScalar%val(NLMAX)%x,myScalar%def,myScalar%aux,&
        myScalar%ndof,4*myScalar%prm%SolvIter,0.7d0)
  ENDIF

 END IF
 CALL GetDefNorm(Amat,lMat%ColA,lMat%LdA,myScalar%val(NLMAX)%x,&
      myScalar%def,myScalar%aux,myScalar%ndof,DefCurrent)
 IF (DefCurrent/DefInit.LT.0.09d0) GOTO 1
END DO

1 CONTINUE

! IF (myid.eq.1) WRITE(*,'(I4,(3D12.4))') iLinIter,DefInit,DefCurrent,DefCurrent/DefInit

IF (myid.ne.0) THEN
 ! Update the solution
 CALL LLC1(myScalar%val_old,myScalar%val(NLMAX)%x,&
      myScalar%ndof,1D0,1D0)

 ! Set dirichlet boundary conditions on the solution
 CALL Bndry_Val(mg_mesh%level(NLMAX)%dcorvg)

END IF

END SUBROUTINE Solve_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE GMVoutput(value)
INTEGER I
REAL*8 , INTENT(IN) :: value(3,*)

 if (myid.eq.1) OPEN(987,FILE='#gmv/VAL1.gmv')
 if (myid.eq.2) OPEN(987,FILE='#gmv/VAL2.gmv')
 if (myid.eq.3) OPEN(987,FILE='#gmv/VAL3.gmv')
 WRITE(987,'(A)') 'gmvinput ascii'
 if (myid.eq.1) WRITE(987,'(A)') 'nodes fromfile "msh_01.gmv"'
 if (myid.eq.2) WRITE(987,'(A)') 'nodes fromfile "msh_02.gmv"'
 if (myid.eq.3) WRITE(987,'(A)') 'nodes fromfile "msh_03.gmv"'
 if (myid.eq.1) WRITE(987,'(A)') 'cells fromfile "msh_01.gmv"'
 if (myid.eq.2) WRITE(987,'(A)') 'cells fromfile "msh_02.gmv"'
 if (myid.eq.3) WRITE(987,'(A)') 'cells fromfile "msh_03.gmv"'
 WRITE(987,'(A)') 'variable'
 WRITE(987,'(A)') 'norm_x 1'
 DO I=1,lMat%nu
   WRITE(987,*) value(1,i)
 END DO
 WRITE(987,'(A)') 'norm_y 1'
 DO I=1,lMat%nu
   WRITE(987,*) value(2,i)
 END DO
 WRITE(987,'(A)') 'norm_z 1'
 DO I=1,lMat%nu
   WRITE(987,*) value(3,i)
 END DO
 WRITE(987,'(A)') 'norm_D 1'
 DO I=1,lMat%nu
   WRITE(987,*) SQRT(value(1,i)**2+value(2,i)**2+value(3,i)**2)
 END DO
 WRITE(987,'(A)') 'endvars'
 WRITE(987,'(A)') 'probtime    1.0'
 WRITE(987,'(A)') 'endgmv'
 CLOSE(987)

END SUBROUTINE GMVoutput
!
! ----------------------------------------------
!
SUBROUTINE InitAFC_General_LinScalar()

CALL AFC_LinScalar(Kmat,lMat%ColA,lMat%LdA,lMat%nu,&
     AFC%isep,AFC%iaux,AFC%inod,AFC%jnod,AFC%aedge)

END SUBROUTINE InitAFC_General_LinScalar
!
! ----------------------------------------------
!
SUBROUTINE Protocol_linScalar(mfile,myScalar,nINL,&
           ResScalar,DefScalar,RhsScalar,cTitle)
TYPE(lscalar), INTENT(INOUT) :: myScalar
INTEGER nINL,mfile
INTEGER i,length
REAL*8 ResScalar,DefScalar,RhsScalar
CHARACTER C1*14,C2*14,C3*14
CHARACTER, OPTIONAL:: cTitle*(*)

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
  WRITE(MTERM,4) cTitle
  WRITE(MFILE,4) cTitle
 ELSE
  WRITE(MTERM,5)
  WRITE(MFILE,5)
 END IF
 WRITE(MTERM,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(MFILE,'(A8,5(2X,A14))') "INL",TRIM(C1),TRIM(C2),TRIM(C3)
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(A8,6XA10,5(6X,D11.4))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(MFILE,'(A8,6XA10,5(6X,D11.4))') "Criteria"," ",DefScalar*myScalar%prm%defCrit,RhsScalar
 WRITE(MTERM,5)
 WRITE(MFILE,5)
 WRITE(MTERM,'(I8,5(6X,D11.4))') 0,ResScalar,DefScalar
 WRITE(MFILE,'(I8,5(6X,D11.4))') 0,ResScalar,DefScalar
ELSE
 WRITE(MTERM,'(I8,5(6X,D11.4))') nINL,ResScalar,DefScalar,RhsScalar
 WRITE(MFILE,'(I8,5(6X,D11.4))') nINL,ResScalar,DefScalar,RhsScalar
END IF

END IF

5  FORMAT(80('-'))
4  FORMAT(80('-'))


END SUBROUTINE Protocol_linScalar
!
! ----------------------------------------------
!
SUBROUTINE Solve_Displacement_Q1(myScalar,knpr,Bndry_Val,Bndry_Mat)
INTEGER KNPR(*)
INTEGER iLinIter,i,ndof
REAL*8 DefInit,DefCurrent,DefInitX,DefInitY,DefInitZ
TYPE(lScalar3), INTENT(INOUT), TARGET :: myScalar
EXTERNAL Bndry_Val,Bndry_Mat

ILEV=NLMAX
CALL SETLEV(2)

! if (myid.eq.0) GOTO 191

IF (myid.ne.0) THEN

 ndof = myScalar%ndof

 LaplaceMat => mg_LaplaceMat(ILEV)%a
! A11Mat => mg_A11Mat(ILEV)%a
! A22Mat => mg_A22Mat(ILEV)%a
! A33Mat => mg_A33Mat(ILEV)%a
! lMat   => mg_lMat(ILEV)

! CALL Bndry_Mat(A11Mat,lMat%LdA,myScalar%knprX)
! CALL Bndry_Mat(A22Mat,lMat%LdA,myScalar%knprY)
! CALL Bndry_Mat(A33Mat,lMat%LdA,myScalar%knprZ)
!        
! CALL E011_UMAT(A11mat,lMat%LdA,lMat%nu,1)
! CALL E011_UMAT(A22mat,lMat%LdA,lMat%nu,2)
! CALL E011_UMAT(A33mat,lMat%LdA,lMat%nu,3)
!
! DO i=1,ndof
!  if (myScalar%knprX(i).eq.1) myScalar%defX(i) = 0d0
!  if (myScalar%knprY(i).eq.1) myScalar%defY(i) = 0d0
!  if (myScalar%knprZ(i).eq.1) myScalar%defZ(i) = 0d0
! END DO
!
! CALL LCL1 (myScalar%valX,ndof)
! CALL LCL1 (myScalar%valY,ndof)
! CALL LCL1 (myScalar%valZ,ndof)
!
! myScalar%sol(NLMAX)%x(0*ndof+1:1*ndof) = myScalar%ValX
! myScalar%sol(NLMAX)%x(1*ndof+1:2*ndof) = myScalar%ValY
! myScalar%sol(NLMAX)%x(2*ndof+1:3*ndof) = myScalar%ValZ
! MyMG%X    => myScalar%sol
!
! myScalar%rhs(NLMAX)%x(0*ndof+1:1*ndof) = myScalar%defX
! myScalar%rhs(NLMAX)%x(1*ndof+1:2*ndof) = myScalar%defY
! myScalar%rhs(NLMAX)%x(2*ndof+1:3*ndof) = myScalar%defZ
! MyMG%B    => myScalar%rhs
!
! CALL GetDefNorm(A11Mat,lMat%ColA,lMat%LdA,myScalar%valX,&
!      myScalar%defX,myScalar%aux(NLMAX)%x,ndof,DefInitX)
! CALL GetDefNorm(A22Mat,lMat%ColA,lMat%LdA,myScalar%valY,&
!      myScalar%defY,myScalar%aux(NLMAX)%x,ndof,DefInitY)
! CALL GetDefNorm(A33Mat,lMat%ColA,lMat%LdA,myScalar%valZ,&
!      myScalar%defZ,myScalar%aux(NLMAX)%x,ndof,DefInitZ)
     
END IF

!CALL COMM_Maximum(DefInitX)
!CALL COMM_Maximum(DefInitY)
!CALL COMM_Maximum(DefInitZ)
!
!DefInit=MAX(DefInitX,DefInitY,DefInitZ)
!
!DO iLinIter=1,5
! IF (myid.ne.0) THEN
!
!   CALL SSORSolverXXX(A11Mat,A22Mat,A33Mat,lMat%ColA,lMat%LdA,&
!        myScalar%sol(ILEV)%x,myScalar%rhs(ILEV)%x,myScalar%aux(ILEV)%x,&
!        KNPR,ndof,myScalar%prm%MGprmIn%nIterCoarse, 0.7d0)
!   CALL GetDefNorm(A11Mat,lMat%ColA,lMat%LdA,myScalar%valX,&
!        myScalar%defX,myScalar%aux(NLMAX)%x,ndof,DefInitX)
!   CALL GetDefNorm(A22Mat,lMat%ColA,lMat%LdA,myScalar%valY,&
!        myScalar%defY,myScalar%aux(NLMAX)%x,ndof,DefInitY)
!   CALL GetDefNorm(A33Mat,lMat%ColA,lMat%LdA,myScalar%valZ,&
!        myScalar%defZ,myScalar%aux(NLMAX)%x,ndof,DefInitZ)
!  END IF
!  
! CALL COMM_Maximum(DefInitX)
! CALL COMM_Maximum(DefInitY)
! CALL COMM_Maximum(DefInitZ)
! DefCurrent=MAX(DefInitX,DefInitY,DefInitZ)
! IF (DefCurrent/DefInit.LT.0.0001d0) GOTO 1
!
!END DO
!
!1 CONTINUE
!
!IF (myid.ne.0) THEN
! ! Update the solution
! myScalar%ValX = myScalar%sol(NLMAX)%x(0*ndof+1:1*ndof)
! myScalar%ValY = myScalar%sol(NLMAX)%x(1*ndof+1:2*ndof)
! myScalar%ValZ = myScalar%sol(NLMAX)%x(2*ndof+1:3*ndof)
!
! CALL LLC1(myScalar%valX_old,myScalar%valX,&
!      myScalar%ndof,1D0,1D0)
! CALL LLC1(myScalar%valY_old,myScalar%valY,&
!      myScalar%ndof,1D0,1D0)
! CALL LLC1(myScalar%valZ_old,myScalar%valZ,&
!      myScalar%ndof,1D0,1D0)
!END IF

END SUBROUTINE Solve_Displacement_Q1
!
! ----------------------------------------------
!
SUBROUTINE GetDispParameters(myParam,cName,mfile)
IMPLICIT NONE
TYPE(tParam) myParam
INTEGER mfile
CHARACTER*7 cName
integer :: iEnd,iAt,iEq
CHARACTER string*100,param*50,cVar*7,cPar*50
LOGICAL bOK
INTEGER :: myFile = 990076
integer :: istat

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
    IF (myid.eq.showid) write(mterm,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%cEquation
    IF (myid.eq.showid) write(mfile,'(A,A)') TRIM(ADJUSTL(cVar))//" - "//TRIM(ADJUSTL(cPar))//" "//"= ",myParam%cEquation
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

END SUBROUTINE GetDispParameters

!
! ----------------------------------------------
!

include 'LinSc_def_extension.f90'

END MODULE def_LinScalar

