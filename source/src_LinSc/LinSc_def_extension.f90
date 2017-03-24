!
! ----------------------------------------------
!
SUBROUTINE Create_AFCStruct_Q1()
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

END SUBROUTINE Create_AFCStruct_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_MatStruct_Q1()
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

END SUBROUTINE Create_MatStruct_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_MassMat_Q1()
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

END SUBROUTINE Create_MassMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_LMassMat_Q1()
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

END SUBROUTINE Create_LMassMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_LKonvMat_Q1()

ALLOCATE(Kmat(lMat%na))

END SUBROUTINE Create_LKonvMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_AMat_Q1()

ALLOCATE(AmatX(lMat%na))
ALLOCATE(AmatY(lMat%na))
ALLOCATE(AmatZ(lMat%na))

END SUBROUTINE Create_AMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_DiffMat_Q1(Alpha)
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

END SUBROUTINE Create_DiffMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_NewDiffMat_Q1(DCOOR,dAlpha)
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier
EXTERNAL E011
REAL*8 dCOOR(3,*)
REAL*8 dAlpha
integer :: i

IF (.not.allocated(DMat)) ALLOCATE(DMat(lMat%na))
DMat=0d0

ILEV=NLMAX
CALL SETLEV(2)

CALL CnstDiffMatQ1(dAlpha,DMat,lMat%nu,lMat%ColA,lMat%LdA,&
mg_mesh%level(ilev)%kvert,&
mg_mesh%level(ilev)%karea,&
mg_mesh%level(ilev)%kedge,&
dCOOR,&
E011)

END SUBROUTINE Create_NewDiffMat_Q1
!
! ----------------------------------------------
!
SUBROUTINE Initialize_Q1(myScalar)
USE var_QuadScalar, only:knvt
TYPE(lScalar3) myScalar

ILEV=NLMAX
CALL SETLEV(2)

myScalar%ndof = KNVT(ilev)
myScalar%na = lMat%na

ALLOCATE(myScalar%knprX(myScalar%ndof))
ALLOCATE(myScalar%knprY(myScalar%ndof))
ALLOCATE(myScalar%knprZ(myScalar%ndof))

ALLOCATE(myScalar%aux(myScalar%ndof))
ALLOCATE(myScalar%rhs(myScalar%ndof))
ALLOCATE(myScalar%defX(myScalar%ndof))
ALLOCATE(myScalar%defY(myScalar%ndof))
ALLOCATE(myScalar%defZ(myScalar%ndof))
ALLOCATE(myScalar%valX_old(myScalar%ndof))
ALLOCATE(myScalar%valY_old(myScalar%ndof))
ALLOCATE(myScalar%valZ_old(myScalar%ndof))

ALLOCATE(myScalar%valX(NLMIN:NLMAX))
ALLOCATE(myScalar%valY(NLMIN:NLMAX))
ALLOCATE(myScalar%valZ(NLMIN:NLMAX))
DO ILEV=NLMIN,NLMAX
 ALLOCATE(myScalar%valX(ILEV)%x(KNVT(ilev)))
 ALLOCATE(myScalar%valY(ILEV)%x(KNVT(ilev)))
 ALLOCATE(myScalar%valZ(ILEV)%x(KNVT(ilev)))
END DO

END SUBROUTINE Initialize_Q1
!
! ----------------------------------------------
!
SUBROUTINE Create_Knpr_Q1(myScalar_Knpr)
EXTERNAL myScalar_Knpr

ILEV=NLMAX
CALL SETLEV(2)

CALL myScalar_Knpr(mg_mesh%level(ilev)%dcorvg)

END SUBROUTINE Create_Knpr_Q1
!
! ----------------------------------------------
!
SUBROUTINE InitCond_Q1(myScalar_InitCond)
EXTERNAL myScalar_InitCond

ILEV=NLMAX
CALL SETLEV(2)

CALL myScalar_InitCond(mg_mesh%level(ilev)%dcorvg)

END SUBROUTINE InitCond_Q1
!
! ----------------------------------------------
!
SUBROUTINE Matdef_general_LinScalar_Q1(myScalar,idef,imat)
INTEGER :: idef,imat
TYPE(lScalar3) myScalar

 ! Build up the matrix
 IF (imat.eq.1) THEN
   AmatX=REAL(-thstep*(DMat))
   AmatY=REAL(-thstep*(DMat))
   AmatZ=REAL(-thstep*(DMat))
!    Amat=REAL(-thstep*(Kmat+DMat))
!    Amat(lMat%LdA(1:lMat%nu))=Amat(lMat%LdA(1:lMat%nu))+REAL(MLmat)
 END IF

 ! Build up the defect
 IF (idef.eq. 1) THEN
  IF (myScalar%prm%iMAss.EQ.1) THEN
     myScalar%defX = 0d0 !MLmat*myScalar%val(NLMAX)%x
     myScalar%defY = 0d0 !MLmat*myScalar%val(NLMAX)%x
     myScalar%defZ = 0d0 !MLmat*myScalar%val(NLMAX)%x
!      CALL LAX17(Kmat,lMat%ColA,lMat%LdA,lMat%nu,&
!      myScalar%val(NLMAX)%x,myScalar%def,thstep,1d0)
     CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valX(NLMAX)%x,myScalar%defX,thstep,1d0)
     CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valY(NLMAX)%x,myScalar%defY,thstep,1d0)
     CALL LAX17(Dmat,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valZ(NLMAX)%x,myScalar%defZ,thstep,1d0)
  END IF
!   IF (myScalar%prm%iMAss.EQ.2) THEN
!      CALL LAX17(Mmat+thstep*(Kmat+Dmat),lMat%ColA,lMat%LdA,lMat%nu,&
!      myScalar%valX(NLMAX)%x,myScalar%def,1d0,0d0)
!   END IF
 ELSE
  IF (myScalar%prm%iMAss.EQ.1) THEN
     CALL LAX37(AmatX,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valX(NLMAX)%x,myScalar%defX,-1d0,1d0)
     CALL LAX37(AmatY,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valY(NLMAX)%x,myScalar%defY,-1d0,1d0)
     CALL LAX37(AmatZ,lMat%ColA,lMat%LdA,lMat%nu,&
     myScalar%valZ(NLMAX)%x,myScalar%defZ,-1d0,1d0)
  END IF
!   IF (myScalar%prm%iMAss.EQ.2) THEN
!      CALL LAX17(-Mmat+thstep*(Kmat+Dmat),lMat%ColA,lMat%LdA,lMat%nu,&
!      myScalar%val(NLMAX)%x,myScalar%def,1d0,1d0)
!   END IF
 END IF

 ! Perform Algebraic Flux Correction (if needed)
 IF (myScalar%prm%AFC) THEN
!    CALL DefTVD_LinScalar(myScalar%val(NLMAX)%x,myScalar%def,THSTEP)
 END IF

END SUBROUTINE Matdef_general_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE Resdfk_General_LinScalar_Q1(myScalar,&
           resScalar,defScalar,rhsScalar)
TYPE(lScalar3), INTENT(INOUT) :: myScalar
REAL*8  resScalar,defScalar,rhsScalar,RESF,RESU

CALL LL21 (myScalar%rhs,myScalar%ndof,RESF)
RESF=MAX(1D-15,RESF)

CALL LL21 (myScalar%defX,myScalar%ndof,RESU)

resScalar = RESU/RESF
defScalar = RESU
rhsScalar = RESF
!write(*,*) RESU,RESF

END SUBROUTINE Resdfk_General_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE Solve_General_LinScalar_Q1(myScalar,knpr,Bndry_Val,Bndry_Mat)
INTEGER KNPR(*)
INTEGER iLinIter,i
REAL*8 DefInit,DefCurrent
TYPE(lScalar3), INTENT(INOUT) :: myScalar
EXTERNAL Bndry_Val,Bndry_Mat

IF (myid.ne.0) THEN
 CALL Bndry_Mat(AmatX,lMat%LdA,myScalar%knprX)

 CALL E011Mat(AmatX,lMat%LdA,lMat%nu)

 DO i=1,myScalar%ndof
  if (myScalar%knprX(i).eq.1) myScalar%defX(i) = 0d0
 END DO

 CALL LCL1 (myScalar%valX(NLMAX)%x,myScalar%ndof)
END IF

IF (myid.eq.1) WRITE(*,'(A,2ES12.4)')  'solving X: ',DefInit,DefCurrent

CALL GetDefNorm(AmatX,lMat%ColA,lMat%LdA,myScalar%valX(NLMAX)%x,&
     myScalar%defX,myScalar%aux,myScalar%ndof,DefInit)

DO iLinIter=1,5
 IF (myid.ne.0) THEN
  IF (myScalar%prm%SolvType.EQ.1) THEN
   CALL SSORSolver(AmatX,lMat%ColA,lMat%LdA,&
        myScalar%valX(NLMAX)%x,myScalar%defX,myScalar%aux,&
        KNPR,myScalar%ndof,1*myScalar%prm%SolvIter,0.7d0)
  ENDIF
  IF (myScalar%prm%SolvType.EQ.2) THEN
   CALL JacobiSolver(AmatX,lMat%ColA,lMat%LdA,&
        myScalar%valX(NLMAX)%x,myScalar%defX,myScalar%aux,&
        myScalar%ndof,4*myScalar%prm%SolvIter,0.7d0)
  ENDIF

 END IF
 CALL GetDefNorm(AmatX,lMat%ColA,lMat%LdA,myScalar%valX(NLMAX)%x,&
      myScalar%defX,myScalar%aux,myScalar%ndof,DefCurrent)
 IF (DefCurrent/DefInit.LT.0.09d0) GOTO 1
END DO

1 CONTINUE


IF (myid.eq.1) WRITE(*,'(A,2ES12.4)')  'solving Y: ',DefInit,DefCurrent

IF (myid.ne.0) THEN
 CALL Bndry_Mat(AmatY,lMat%LdA,myScalar%knprY)

 CALL E011Mat(AmatY,lMat%LdA,lMat%nu)

 DO i=1,myScalar%ndof
  if (myScalar%knprY(i).eq.1) myScalar%defY(i) = 0d0
 END DO

 CALL LCL1 (myScalar%valY(NLMAX)%x,myScalar%ndof)
END IF

CALL GetDefNorm(AmatY,lMat%ColA,lMat%LdA,myScalar%valY(NLMAX)%x,&
     myScalar%defY,myScalar%aux,myScalar%ndof,DefInit)

DO iLinIter=1,5
 IF (myid.ne.0) THEN
  IF (myScalar%prm%SolvType.EQ.1) THEN
   CALL SSORSolver(AmatY,lMat%ColA,lMat%LdA,&
        myScalar%valY(NLMAX)%x,myScalar%defY,myScalar%aux,&
        KNPR,myScalar%ndof,1*myScalar%prm%SolvIter,0.7d0)
  ENDIF
  IF (myScalar%prm%SolvType.EQ.2) THEN
   CALL JacobiSolver(AmatY,lMat%ColA,lMat%LdA,&
        myScalar%valY(NLMAX)%x,myScalar%defY,myScalar%aux,&
        myScalar%ndof,4*myScalar%prm%SolvIter,0.7d0)
  ENDIF

 END IF
 CALL GetDefNorm(AmatY,lMat%ColA,lMat%LdA,myScalar%valY(NLMAX)%x,&
      myScalar%defY,myScalar%aux,myScalar%ndof,DefCurrent)
 IF (DefCurrent/DefInit.LT.0.09d0) GOTO 2
END DO

2 CONTINUE

IF (myid.eq.1) WRITE(*,'(A,2ES12.4)')  'solving Z: ',DefInit,DefCurrent

IF (myid.ne.0) THEN
 CALL Bndry_Mat(AmatZ,lMat%LdA,myScalar%knprZ)

 CALL E011Mat(AmatZ,lMat%LdA,lMat%nu)

 DO i=1,myScalar%ndof
  if (myScalar%knprZ(i).eq.1) myScalar%defZ(i) = 0d0
 END DO

 CALL LCL1 (myScalar%valZ(NLMAX)%x,myScalar%ndof)
END IF

CALL GetDefNorm(AmatZ,lMat%ColA,lMat%LdA,myScalar%valZ(NLMAX)%x,&
     myScalar%defZ,myScalar%aux,myScalar%ndof,DefInit)

DO iLinIter=1,5
 IF (myid.ne.0) THEN
  IF (myScalar%prm%SolvType.EQ.1) THEN
   CALL SSORSolver(AmatZ,lMat%ColA,lMat%LdA,&
        myScalar%valZ(NLMAX)%x,myScalar%defZ,myScalar%aux,&
        KNPR,myScalar%ndof,1*myScalar%prm%SolvIter,0.7d0)
  ENDIF
  IF (myScalar%prm%SolvType.EQ.2) THEN
   CALL JacobiSolver(AmatZ,lMat%ColA,lMat%LdA,&
        myScalar%valZ(NLMAX)%x,myScalar%defZ,myScalar%aux,&
        myScalar%ndof,4*myScalar%prm%SolvIter,0.7d0)
  ENDIF

 END IF
 CALL GetDefNorm(AmatZ,lMat%ColA,lMat%LdA,myScalar%valZ(NLMAX)%x,&
      myScalar%defZ,myScalar%aux,myScalar%ndof,DefCurrent)
 IF (DefCurrent/DefInit.LT.0.09d0) GOTO 3
END DO

3 CONTINUE

! IF (myid.eq.1) WRITE(*,'(I4,(3D12.4))') iLinIter,DefInit,DefCurrent,DefCurrent/DefInit

IF (myid.ne.0) THEN
 ! Update the solution
 CALL LLC1(myScalar%valX_old,myScalar%valX(NLMAX)%x,&
      myScalar%ndof,1D0,1D0)
 CALL LLC1(myScalar%valY_old,myScalar%valY(NLMAX)%x,&
      myScalar%ndof,1D0,1D0)
 CALL LLC1(myScalar%valZ_old,myScalar%valZ(NLMAX)%x,&
      myScalar%ndof,1D0,1D0)
 ! Set dirichlet boundary conditions on the solution
 ! CALL Bndry_Val(mg_mesh%level(NLMAX)%dcorvg)
END IF

END SUBROUTINE Solve_General_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE GMVoutput_Q1(value)
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

END SUBROUTINE GMVoutput_Q1
!
! ----------------------------------------------
!
SUBROUTINE InitAFC_General_LinScalar_Q1()

CALL AFC_LinScalar(Kmat,lMat%ColA,lMat%LdA,lMat%nu,&
     AFC%isep,AFC%iaux,AFC%inod,AFC%jnod,AFC%aedge)

END SUBROUTINE InitAFC_General_LinScalar_Q1
!
! ----------------------------------------------
!
SUBROUTINE Protocol_linScalar_Q1(mfile,myScalar,nINL,&
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


END SUBROUTINE Protocol_linScalar_Q1

