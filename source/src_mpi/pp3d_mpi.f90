MODULE PP3D_MPI
  USE OctTreeSearch
  IMPLICIT NONE

  include 'mpif.h'
  INTEGER IERR
  INTEGER Variable,iunit
  INTEGER:: ShowID = 1
  INTEGER myid,numnodes,subnodes,MPI_COMM_SUBS,MASTER,MPI_COMM_SUBGROUP
  INTEGER NParNodes(9)
  CHARACTER nodefile*60

  TYPE TE013ST
    INTEGER :: Num = 0
    INTEGER, DIMENSION(2) :: nElems = (/0, 0/)
    INTEGER, DIMENSION(2) :: nEntries = (/0, 0/)
    INTEGER,DIMENSION (:,:), ALLOCATABLE :: VertLink
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: SDVect,RDVect
    REAL*4 , DIMENSION(:)  , ALLOCATABLE :: SVVect,RVVect
    LOGICAL, DIMENSION(:)  , ALLOCATABLE :: SBVect,RBVect
  END TYPE TE013ST
  TYPE(TE013ST), DIMENSION(:), ALLOCATABLE :: E013ST
  TYPE(TE013ST), DIMENSION(:), ALLOCATABLE :: E013SP
  REAL*8 , DIMENSION(:)  , ALLOCATABLE     :: E013_UE

  TYPE TMGE013
    TYPE(TE013ST), DIMENSION(:), ALLOCATABLE :: ST
    TYPE(TE013ST), DIMENSION(:), ALLOCATABLE :: SP
    REAL*8 , DIMENSION(:)  , ALLOCATABLE     :: UE,UE11,UE22,UE33
    REAL*8 , DIMENSION(:)  , ALLOCATABLE     :: CRSVect
  END TYPE

  TYPE(TMGE013), DIMENSION(:), ALLOCATABLE :: MGE013

  TYPE TDBLVect
   REAL*8 , DIMENSION(:)  , ALLOCATABLE     :: x
  END TYPE

  TYPE TMGE011
   REAL*8 , DIMENSION(:)  , ALLOCATABLE     :: UE11,UE22,UE33
   TYPE(TDBLVect), DIMENSION(:), ALLOCATABLE :: UE
  END TYPE
  TYPE(TMGE011), DIMENSION(:), ALLOCATABLE :: MGE011

  TYPE TE012ST
    INTEGER :: Num = 0
    INTEGER :: Neigh = 0
    INTEGER,DIMENSION (:,:), ALLOCATABLE :: FaceLink
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: SDVect,RDVect
    REAL*4 , DIMENSION(:)  , ALLOCATABLE :: SVVect,RVVect
    INTEGER, DIMENSION(:)  , ALLOCATABLE :: SKVect,RKVect
    LOGICAL, DIMENSION(:)  , ALLOCATABLE :: SBVect,RBVect
  END TYPE TE012ST
  TYPE(TE012ST), DIMENSION(:), ALLOCATABLE :: E012ST

  TYPE TE011ST
    INTEGER :: Num = 0
    INTEGER,DIMENSION (:,:), ALLOCATABLE :: VertLink
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: SDVect,RDVect
    REAL*4 , DIMENSION(:)  , ALLOCATABLE :: SVVect,RVVect
    LOGICAL, DIMENSION(:)  , ALLOCATABLE :: SBVect,RBVect
  END TYPE TE011ST
  TYPE(TE011ST), DIMENSION(:), ALLOCATABLE :: E011ST
  REAL*4 , DIMENSION(:)  , ALLOCATABLE     :: E011_UE

  TYPE TVert
    INTEGER :: Num = 0
    INTEGER,DIMENSION (:)  , ALLOCATABLE :: iCoor
    REAL*8, DIMENSION (:,:), ALLOCATABLE :: dCoor
  END TYPE TVert
  TYPE(TVert), DIMENSION(:), ALLOCATABLE :: CoorST,CoorSP

  TYPE Tcoarse
    INTEGER, DIMENSION (:), ALLOCATABLE :: myELEMLINK,myVERTLINK
    INTEGER, DIMENSION (:,:), ALLOCATABLE :: pELEMLINK,pFACELINK,pVERTLINK
    REAL*8, DIMENSION (:),   ALLOCATABLE :: pDX,DX
    INTEGER, DIMENSION (:), ALLOCATABLE :: pNEL,pNVT
    INTEGER :: pElem, pFace,pVert
  END TYPE Tcoarse
  TYPE(Tcoarse) :: coarse

  TYPE Tside
    INTEGER :: Neigh = 0
    INTEGER :: Num = 0
    INTEGER :: i = 0
    INTEGER, DIMENSION(:)  , ALLOCATABLE :: SideLink
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemLink,FaceLink
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ElemLin1,FaceLin1,ElemLin2,FaceLin2
    REAL*8 , DIMENSION(:,:), ALLOCATABLE :: CoragLinkX,CoragLinkY,CoragLinkZ
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: PE
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: SDVect,RDVect
    REAL*4 , DIMENSION(:)  , ALLOCATABLE :: SVVect,RVVect
  END TYPE Tside
  TYPE TparST
    INTEGER :: NeighNum
    REAL*8 , DIMENSION(:)  , ALLOCATABLE :: UE
    TYPE(Tside), DIMENSION(:), ALLOCATABLE :: parST
  END TYPE TparST

  TYPE(TparST), DIMENSION(:), ALLOCATABLE :: mg_mpi
  INTEGER NLMINp

  integer, allocatable :: VerticeCommunicationScheme(:)
  
  TYPE tmQ2
   INTEGER, ALLOCATABLE :: x(:)
  END TYPE tmQ2
  TYPE (tmQ2), ALLOCATABLE :: mQ2(:)

  REAL*8 :: DEpsPrec = 1d-5
  REAL*8 :: dZPeriodicLength,dPeriodicity(3)=[1d9,1d9,1d9]
  ! -------------- Subroutines -------------------
CONTAINS
  ! ----------------------------------------------
  SUBROUTINE INIT_MPI !ok
    implicit none

    INTEGER I
    INTEGER NULLNODE(1)
    INTEGER MPI_COMM_SUBS0,MPI_COMM_SUBS1

    MASTER=0
    NULLNODE(1)=MASTER

    CALL MPI_INIT(IERR)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,IERR)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,numnodes,IERR)
    subnodes=numnodes-1

    CALL MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_COMM_SUBS0,IERR)
    CALL MPI_GROUP_EXCL(MPI_COMM_SUBS0,1,NULLNODE,MPI_COMM_SUBS1,IERR)
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_COMM_SUBS1,MPI_COMM_SUBS,IERR)

    if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

    ! DO I=1,numnodes
    !    iunit = 110+I
    !    nodefile = '#data/proc     '
    !    IF (numnodes.lt.10) THEN
    !     WRITE(nodefile(11:16),'(I1,I1,A4)') 0,I,'.txt'
    !    ELSE
    !     WRITE(nodefile(11:16),'(I2,A4)') I,'.txt'
    !    END IF
    ! 
    !    IF (myid.eq.i-1) THEN
    !     OPEN (iunit,FILE=nodefile)
    !     WRITE(iunit,*) "Report file ..."
    !     WRITE(iunit,*) myid+1,"out of",numnodes
    !     CLOSE(iunit)
    !    END IF
    ! 
    ! END DO

    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  END SUBROUTINE INIT_MPI
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CREATECOMM(ILEV,NAT,NEL,NVT,DCORAG,DCORVG,&
    KADJ,KAREA,KVERT,BLIN) !ok
  implicit none
  REAL*8  DCORAG(3,*),DCORVG(3,*)
  INTEGER KAREA(6,*),KADJ(6,*),KVERT(8,*)
  INTEGER NAT,NEL,NVT,ILEV
  INTEGER IAT,IMT,IEL
  INTEGER IA1,IA2,IA3,IA4,IM1,IM2,IM3,IM4,IE1,IE2,IE3,IE4
  INTEGER IEL1,IEL2,IEL3,IEL4,IEL5,IEL6,IEL7,IEL8
  INTEGER I,J,JI,K
  INTEGER pID,pJD,nSIZE,nEIGH

  REAL*8,  ALLOCATABLE :: dCoor(:,:)
  INTEGER, ALLOCATABLE :: iAux(:),iCoor(:)
  INTEGER IV1,IV2,IV3,IV4,IVT,IVT1,IVT2,IVT3,IVT4
  INTEGER jAux,pNVT
  REAL*8  P1X,P2X,P1Y,P2Y,P1Z,P2Z,pPoint(3),dist
  LOGICAL BLIN
  REAL :: txt1, txt0

  NLMINp = ILEV

  if(myid.eq.0)then
    goto 555
  endif

  mg_mpi(ILEV)%NeighNum = mg_mpi(ILEV-1)%NeighNum
  ALLOCATE (mg_mpi(ILEV)%parST(mg_mpi(ILEV)%NeighNum))
  ALLOCATE (mg_mpi(ILEV)%UE(NAT))

  !  DO pID=1,mg_mpi(ILEV-1)%NeighNum
  !  WRITE(*,*) (mg_mpi(ILEV-1)%parST(pID)%ElemLink(1,I),I=1,mg_mpi(ILEV-1)%parST(pID)%Num)
  !  WRITE(*,*) (mg_mpi(ILEV-1)%parST(pID)%FaceLink(1,I),I=1,mg_mpi(ILEV-1)%parST(pID)%Num)
  !  WRITE(*,*) (mg_mpi(ILEV-1)%parST(pID)%SideLink(  I),I=1,mg_mpi(ILEV-1)%parST(pID)%Num)
  !  END DO

  !  WRITE(*,*)myid," has arrived ...."
  !  STOP
  DO pID=1,mg_mpi(ILEV)%NeighNum
  nSIZE = mg_mpi(ILEV-1)%parST(pID)%Num*4
  mg_mpi(ILEV)%parST(pID)%Num = nSIZE
  mg_mpi(ILEV)%parST(pID)%Neigh=mg_mpi(ILEV-1)%parST(pID)%Neigh
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%ElemLink(2,nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%FaceLink(2,nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%SideLink(  nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%CoragLinkX(2,nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%CoragLinkY(2,nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%CoragLinkZ(2,nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%PE(nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%SDVect(nSIZE),mg_mpi(ILEV)%parST(pID)%RDVect(nSIZE))
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%SVVect(nSIZE),mg_mpi(ILEV)%parST(pID)%RVVect(nSIZE))

  ALLOCATE(mg_mpi(ILEV)%parST(pID)%ElemLin1(2,nSIZE)) ! sorted elems
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%FaceLin1(2,nSIZE)) ! sorted faces
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%ElemLin2(2,nSIZE)) ! sorted elems
  ALLOCATE(mg_mpi(ILEV)%parST(pID)%FaceLin2(2,nSIZE)) ! sorted faces
  END DO
  
  DO pID=1,mg_mpi(ILEV-1)%NeighNum
  !  if (myid.eq.1) write(*,*) mg_mpi(ILEV-1)%parST(pID)%Neigh
  DO I=1,mg_mpi(ILEV-1)%parST(pID)%Num
  IEL = mg_mpi(ILEV-1)%parST(pID)%ElemLink(1,I)
  IMT = mg_mpi(ILEV-1)%parST(pID)%FaceLink(1,I)
  IAT = mg_mpi(ILEV-1)%parST(pID)%SideLink(I)

  !   write(*,*) IAT,IMT,IEL
  IEL1=IEL
  IEL2=KADJ(3,IEL1)
  IEL3=KADJ(3,IEL2)
  IEL4=KADJ(3,IEL3)
  IEL5=KADJ(6,IEL1)
  IEL6=KADJ(6,IEL2)
  IEL7=KADJ(6,IEL3)
  IEL8=KADJ(6,IEL4)

  IF (IAT.EQ.1) THEN
    IA1 = 1;  IE1 = IEL1; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 1;  IE2 = IEL2; IM2 = KAREA(IA2,IE2) !ok
    IA3 = 1;  IE3 = IEL3; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 1;  IE4 = IEL4; IM4 = KAREA(IA4,IE4) !ok
  END IF

  IF (IAT.EQ.2) THEN
    IA1 = 2;  IE1 = IEL1; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 5;  IE2 = IEL2; IM2 = KAREA(IA2,IE2) !ok 
    IA3 = 2;  IE3 = IEL5; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 5;  IE4 = IEL6; IM4 = KAREA(IA4,IE4) !ok
  END IF

  IF (IAT.EQ.3) THEN
    IA1 = 2;  IE1 = IEL2; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 5;  IE2 = IEL3; IM2 = KAREA(IA2,IE2) !ok
    IA3 = 2;  IE3 = IEL6; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 5;  IE4 = IEL7; IM4 = KAREA(IA4,IE4) !ok
  END IF

  IF (IAT.EQ.4) THEN
    IA1 = 2;  IE1 = IEL3; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 5;  IE2 = IEL4; IM2 = KAREA(IA2,IE2) !ok
    IA3 = 2;  IE3 = IEL7; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 5;  IE4 = IEL8; IM4 = KAREA(IA4,IE4) !ok
  END IF

  IF (IAT.EQ.5) THEN
    IA1 = 2;  IE1 = IEL4; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 5;  IE2 = IEL1; IM2 = KAREA(IA2,IE2) !ok
    IA3 = 2;  IE3 = IEL8; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 5;  IE4 = IEL5; IM4 = KAREA(IA4,IE4) !ok
  END IF

  IF (IAT.EQ.6) THEN
    IA1 = 1;  IE1 = IEL5; IM1 = KAREA(IA1,IE1) !ok
    IA2 = 1;  IE2 = IEL6; IM2 = KAREA(IA2,IE2) !ok
    IA3 = 1;  IE3 = IEL7; IM3 = KAREA(IA3,IE3) !ok
    IA4 = 1;  IE4 = IEL8; IM4 = KAREA(IA4,IE4) !ok
  END IF

  !   write(*,*) IA1,IA2,IA3,IA4,IM1,IM2,IM3,IM4,IE1,IE2,IE3,IE4
  mg_mpi(ILEV)%parST(pID)%ElemLink(1,4*(I-1)+1) = IE1
  mg_mpi(ILEV)%parST(pID)%ElemLink(1,4*(I-1)+2) = IE2
  mg_mpi(ILEV)%parST(pID)%ElemLink(1,4*(I-1)+3) = IE3
  mg_mpi(ILEV)%parST(pID)%ElemLink(1,4*(I-1)+4) = IE4

  mg_mpi(ILEV)%parST(pID)%FaceLink(1,4*(I-1)+1) = IM1
  mg_mpi(ILEV)%parST(pID)%FaceLink(1,4*(I-1)+2) = IM2
  mg_mpi(ILEV)%parST(pID)%FaceLink(1,4*(I-1)+3) = IM3
  mg_mpi(ILEV)%parST(pID)%FaceLink(1,4*(I-1)+4) = IM4

  mg_mpi(ILEV)%parST(pID)%SideLink(  4*(I-1)+1) = IA1
  mg_mpi(ILEV)%parST(pID)%SideLink(  4*(I-1)+2) = IA2
  mg_mpi(ILEV)%parST(pID)%SideLink(  4*(I-1)+3) = IA3
  mg_mpi(ILEV)%parST(pID)%SideLink(  4*(I-1)+4) = IA4

  mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,4*(I-1)+1) = DCORAG(1,IM1)
  mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,4*(I-1)+2) = DCORAG(1,IM2)
  mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,4*(I-1)+3) = DCORAG(1,IM3)
  mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,4*(I-1)+4) = DCORAG(1,IM4)

  mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,4*(I-1)+1) = DCORAG(2,IM1)
  mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,4*(I-1)+2) = DCORAG(2,IM2)
  mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,4*(I-1)+3) = DCORAG(2,IM3)
  mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,4*(I-1)+4) = DCORAG(2,IM4)

  mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,4*(I-1)+1) = DCORAG(3,IM1)
  mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,4*(I-1)+2) = DCORAG(3,IM2)
  mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,4*(I-1)+3) = DCORAG(3,IM3)
  mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,4*(I-1)+4) = DCORAG(3,IM4)

  !    IF (IAT.EQ.6) THEN
  !      write(*,'(4(1xG12.5))') DCORAG(1,IM1),DCORAG(1,IM2),DCORAG(1,IM3),DCORAG(1,IM4)
  !      write(*,'(4(1xG12.5))') DCORAG(2,IM1),DCORAG(2,IM2),DCORAG(2,IM3),DCORAG(2,IM4)
  !      write(*,'(4(1xG12.5))') DCORAG(3,IM1),DCORAG(3,IM2),DCORAG(3,IM3),DCORAG(3,IM4)
  !    END IF
  !    if (myid.eq.1) then
  !     write(*,*) DCORAG(1,IM1),DCORAG(2,IM1),DCORAG(3,IM1),IAT
  !     write(*,*) DCORAG(1,IM2),DCORAG(2,IM2),DCORAG(3,IM2),IAT
  !     write(*,*) DCORAG(1,IM3),DCORAG(2,IM3),DCORAG(3,IM3),IAT
  !     write(*,*) DCORAG(1,IM4),DCORAG(2,IM4),DCORAG(3,IM4),IAT
  !    end if

  END DO

  !  pause
  CALL SORTALL(mg_mpi(ILEV)%parST(pID)%FaceLink(1,:)  ,&
    mg_mpi(ILEV)%parST(pID)%ElemLink(1,:)  ,&
    mg_mpi(ILEV)%parST(pID)%SideLink(:)    ,&
    mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,:),&
    mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,:),&
    mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,:),&
    mg_mpi(ILEV)%parST(pID)%Num)

  END DO

  DO pID=1,subnodes
  IF (myid.NE.pID) THEN
    !   WRITE(*,*)myid,"'s neighnum",mg_mpi(ILEV)%NeighNum
    DO pJD=1,mg_mpi(ILEV)%NeighNum
    IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkX(1,:),nSIZE,pID)
      CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkY(1,:),nSIZE,pID)
      CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkZ(1,:),nSIZE,pID)
    END IF
    END DO
  ELSE
    DO pJD=1,mg_mpi(ILEV)%NeighNum
    nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
    nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh
    CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkX(2,:),nSize,nEIGH)
    CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkY(2,:),nSize,nEIGH)
    CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%CoragLinkZ(2,:),nSize,nEIGH)
    END DO
  END IF
  END DO

  DO pID=1,mg_mpi(ILEV)%NeighNum
  nSIZE = mg_mpi(ILEV)%parST(pID)%Num
  !   write(*,*)myid,"|",(mg_mpi(ILEV)%parST(pID)%ElemLink(1,I),I=1,nSIZE)
  DO I=1,nSIZE
  !      if (myid.eq.1) then
  !      WRITE(*,'(6(G12.4))') mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,I),&
  !                            mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,I),&
  !                            mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,I),&
  !                            mg_mpi(ILEV)%parST(pID)%CoragLinkX(2,I),&
  !                            mg_mpi(ILEV)%parST(pID)%CoragLinkY(2,I),&
  !                            mg_mpi(ILEV)%parST(pID)%CoragLinkZ(2,I)
  !      end if
  DO J=1,nSIZE
     IF (((ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkX(2,J)).LT.DEpsPrec).OR.&
           (ABS(ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkX(2,J))-dPeriodicity(1)).LT.DEpsPrec)).AND.&
         ((ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkY(2,J)).LT.DEpsPrec).OR.&
           (ABS(ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkY(2,J))-dPeriodicity(2)).LT.DEpsPrec)).AND.&
         ((ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkZ(2,J)).LT.DEpsPrec).OR.&
           (ABS(ABS(mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,I)-mg_mpi(ILEV)%parST(pID)%CoragLinkZ(2,J))-dPeriodicity(3)).LT.DEpsPrec))) THEN

  !        if (myid.eq.1.and.ilev.eq.2) write(*,*) I,J
  mg_mpi(ILEV)%parST(pID)%FaceLink(2,J) = mg_mpi(ILEV)%parST(pID)%FaceLink(1,I)
  mg_mpi(ILEV)%parST(pID)%ElemLink(2,J) = mg_mpi(ILEV)%parST(pID)%ElemLink(1,I)
  GOTO 2
END IF
END DO
WRITE(*,*) "PROBLEM with parallel assignement ...",myid,ilev
WRITE(*,*) mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,I),&
  mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,I),&
  mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,I)

2   CONTINUE 
END DO
END DO

!   DO pID=1,mg_mpi(ILEV)%NeighNum
!   nSIZE = mg_mpi(ILEV)%parST(pID)%Num
!   !   write(*,*)myid,"|",(mg_mpi(ILEV)%parST(pID)%ElemLink(1,I),I=1,nSIZE)
!   CALL InitOctTreeC(mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,:),&
!                     mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,:),&
!                     mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,:),&
!                     nsize)
!   DO I=1,nSIZE
!   pPoint = [mg_mpi(ILEV)%parST(pID)%CoragLinkX(2,i),&
!             mg_mpi(ILEV)%parST(pID)%CoragLinkY(2,i),&
!             mg_mpi(ILEV)%parST(pID)%CoragLinkZ(2,i)]
!   CALL FindInOctTreeC(mg_mpi(ILEV)%parST(pID)%CoragLinkX(1,:),&
!                       mg_mpi(ILEV)%parST(pID)%CoragLinkY(1,:),&
!                       mg_mpi(ILEV)%parST(pID)%CoragLinkZ(1,:),&
!                       nsize,pPoint,J,DIST)
!       IF (J.lt.0) then
!        WRITE(*,*) I,"PROBLEM of parallel assignement ..."
!       end if
!       IF (DIST.LT.DEpsPrec) THEN
!        mg_mpi(ILEV)%parST(pID)%FaceLink(2,J) = mg_mpi(ILEV)%parST(pID)%FaceLink(1,I)
!        mg_mpi(ILEV)%parST(pID)%ElemLink(2,J) = mg_mpi(ILEV)%parST(pID)%ElemLink(1,I)
!       END IF
!   END DO
!     CALL FreeOctTree()
!   END DO

DEALLOCATE(mg_mpi(ILEV)%parST(1)%CoragLinkX)
DEALLOCATE(mg_mpi(ILEV)%parST(1)%CoragLinkY)
DEALLOCATE(mg_mpi(ILEV)%parST(1)%CoragLinkZ)

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

555 IF (.NOT.BLIN) RETURN

IF (myid.NE.MASTER) THEN

  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
  CALL ZTIME(txt0)

  ALLOCATE (iAux(1:nvt))
  ALLOCATE (dCoor(3,1:nvt))
  ALLOCATE (iCoor(1:nvt))

  ALLOCATE(CoorSP(subnodes))
  CoorSP(:)%Num = 0

  DO pID=1,mg_mpi(ILEV)%NeighNum
  nSIZE = mg_mpi(ILEV)%parST(pID)%Num
  pJD   = mg_mpi(ILEV)%parST(pID)%Neigh

  IF (nSIZE.EQ.0) THEN
    CYCLE
  END IF

  iAux = 0
  jAux = 0
  DO I=1,nSIZE
  IEL = mg_mpi(ILEV)%parST(pID)%ElemLink(1,I)
  IAT = mg_mpi(ILEV)%parST(pID)%SideLink(  I)

  IF (IAT.eq.1) THEN
    IV1 = 1; IV2 = 2; IV3 = 3; IV4 = 4;
  END IF

  IF (IAT.eq.2) THEN
    IV1 = 1; IV2 = 2; IV3 = 6; IV4 = 5;
  END IF

  IF (IAT.eq.3) THEN
    IV1 = 2; IV2 = 3; IV3 = 7; IV4 = 6;
  END IF

  IF (IAT.eq.4) THEN
    IV1 = 3; IV2 = 4; IV3 = 8; IV4 = 7;
  END IF

  IF (IAT.eq.5) THEN
    IV1 = 4; IV2 = 1; IV3 = 5; IV4 = 8;
  END IF

  IF (IAT.eq.6) THEN
    IV1 = 5; IV2 = 6; IV3 = 7; IV4 = 8;
  END IF

  IVT1 = KVERT(IV1,IEL); IVT2 = KVERT(IV2,IEL)
  IVT3 = KVERT(IV3,IEL); IVT4 = KVERT(IV4,IEL)

  IF (iAUX(IVT1).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT2).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT3).EQ.0) jAux = jAux + 1
  IF (iAUX(IVT4).EQ.0) jAux = jAux + 1
  iAUX(IVT1) = 1
  iAUX(IVT2) = 1
  iAUX(IVT3) = 1
  iAUX(IVT4) = 1

  END DO

  ALLOCATE(CoorSP(pJD)%dCoor(3,jAux))
  ALLOCATE(CoorSP(pJD)%iCoor(  jAux))
  CoorSP(pJD)%Num = jAux

  jAux = 0
  DO IVT=1,NVT
  IF (iAux(IVT).EQ.1) THEN
    jAux = jAux + 1
    CoorSP(pJD)%dCoor(1,jAux) = DCORVG(1,IVT)
    CoorSP(pJD)%dCoor(2,jAux) = DCORVG(2,IVT)
    CoorSP(pJD)%dCoor(3,jAux) = DCORVG(3,IVT)
    CoorSP(pJD)%iCoor(  jAux) = IVT
  END IF
  END DO

  END DO

  ALLOCATE(CoorST(subnodes))
  DO pID=1,subnodes
  IF (myid.eq.pID) THEN
    DO pJD=1,subnodes
    IF (myid.NE.pJD) THEN
      CALL RECVI_myMPI(pNVT,pJD)
      CoorST(pJD)%Num = pNVT
      IF (pNVT.GT.0) THEN
        ALLOCATE(CoorST(pJD)%dCoor(3,pNVT))
        ALLOCATE(CoorST(pJD)%iCoor(  pNVT))
        CALL RECVD_myMPI(CoorST(pJD)%dCoor,3*pNVT,pJD)
        CALL RECVK_myMPI(CoorST(pJD)%iCoor,  pNVT,pJD)
      END IF
      !    ELSE
      !     CoorST(pJD)%Num = jAux
      !     ALLOCATE(CoorST(pJD)%dCoor(3,jAux))
      !     ALLOCATE(CoorST(pJD)%iCoor(  jAux))
      !     CoorST(pJD)%dCoor(:,:) = dCoor(:,1:jAux)
      !     CoorST(pJD)%iCoor(  :) = iCoor(  1:jAux)
    END IF
    END DO
  ELSE
    jAux = CoorSP(pID)%Num
    CALL SENDI_myMPI(jAux ,pID)
    IF (jAux.GT.0) THEN
      CALL SENDD_myMPI(CoorSP(pID)%dCoor,3*jAux,pID)
      CALL SENDK_myMPI(CoorSP(pID)%iCoor,  jAux,pID)
    END IF
  END IF

  END DO

  ALLOCATE(E011ST(subnodes))

  DO pID=1,subnodes
  E011ST(pID)%Num = CoorSP(pID)%Num
  END DO

  ALLOCATE(E011_UE(NVT))

  DO pID=1,subnodes
  jAux = 0
  IF (E011ST(pID)%Num.GT.0.AND.pID.NE.myid) THEN
    ALLOCATE(E011ST(pID)%VertLink(2,E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%SVVect  (  E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%RVVect  (  E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%SDVect  (  E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%RDVect  (  E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%SBVect  (  E011ST(pID)%Num))
    ALLOCATE(E011ST(pID)%RBVect  (  E011ST(pID)%Num))

!     CALL InitOctTree(CoorSP(pID)%dCoor,CoorSP(pID)%Num)
    
    DO I=1,CoorST(pID)%Num
!      CALL FindInOctTree(CoorSP(pID)%dCoor,CoorSP(pID)%Num,CoorST(pID)%dCoor(:,I),J,DIST)
!       IF (J.lt.0) then
!        WRITE(*,*) I,"PROBLEM of parallel assignement SP/ST..."
!       end if
!       IF (DIST.LT.DEpsPrec) THEN
!        jAux = jAux + 1
!        E011ST(pID)%VertLink(1,jAux) = CoorSP(pID)%iCoor(J)
!        E011ST(pID)%VertLink(2,jAux) = CoorSP(pID)%iCoor(J)
!       END IF
     
    P1X = CoorST(pID)%dCoor(1,I)
    P1Y = CoorST(pID)%dCoor(2,I)
    P1Z = CoorST(pID)%dCoor(3,I)
    DO J=1,CoorSP(pID)%Num
    P2X = CoorSP(pID)%dCoor(1,J)
    P2Y = CoorSP(pID)%dCoor(2,J)
    P2Z = CoorSP(pID)%dCoor(3,J)
    IF ((ABS(P1X-P2X).LT.DEpsPrec).AND.&
      (ABS(P1Y-P2Y).LT.DEpsPrec).AND.&
      (ABS(P1Z-P2Z).LT.DEpsPrec)) THEN
    jAux = jAux + 1
    E011ST(pID)%VertLink(1,jAux) = CoorSP(pID)%iCoor(J)
    E011ST(pID)%VertLink(2,jAux) = CoorSP(pID)%iCoor(J)
    EXIT
  END IF
  END DO
  END DO
  CALL SORT1D(E011ST(pID)%VertLink(1,:),E011ST(pID)%Num)
!   CALL FreeOctTree()
END IF
END DO
! 
CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
CALL ZTIME(txt1)
IF (myid.eq.1) write(*,*) "measured full time:", txt1-txt0

DEALLOCATE (CoorSP,CoorST)
DEALLOCATE (iAux,dCoor,iCoor)

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CREATECOMM
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011Sum(FX) !ok
USE def_feat, ONLY: ILEV,NLMIN,NLMAX
USE var_QuadScalar, only: knvt
IMPLICIT NONE

REAL*8  FX(*)
INTEGER I,pID,pJD,nSIZE,nEIGH,NU

IF (myid.ne.MASTER) THEN

 DO pID=1,subnodes
  IF (myid.NE.pID) THEN
   DO pJD=1,subnodes
    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
     nSIZE = E011ST(pJD)%Num

     NU = KNVT(ILEV)
     
     DO I=1,nSIZE
      IF (E011ST(pJD)%VertLink(1,I).le.NU) THEN
       E011ST(pJD)%SDVect(I) = FX(E011ST(pJD)%VertLink(1,I))
      end if
     END DO

     CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)

    END IF
   END DO
  ELSE
   DO pJD=1,subnodes
     IF (E011ST(pJD)%Num.GT.0) THEN
      nSIZE = E011ST(pJD)%Num

      CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)

     END IF
   END DO
  END IF
 END DO

 DO pJD=1,subnodes

   NU = KNVT(ILEV)
   IF (ILEV.lt.NLMIN.OR.ILEV.GT.NLMAX+1) THEN
    WRITE(*,*) myid,ILEV,'problem'
    pause
   END IF
   IF (E011ST(pJD)%Num.GT.0) THEN
     nSIZE = E011ST(pJD)%Num

     DO I=1,nSIZE
      IF (E011ST(pJD)%VertLink(2,I).le.NU) THEN
       FX(E011ST(pJD)%VertLink(2,I)) = &
       FX(E011ST(pJD)%VertLink(2,I)) +  E011ST(pJD)%RDVect(I)
      END IF
     END DO

   END IF
 END DO

END IF

if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011Sum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011_ABSMIN(FX) !ok

  REAL*8  FX(*)
  INTEGER I,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,subnodes
      IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        DO I=1,nSIZE
        E011ST(pJD)%SDVect(I) = FX(E011ST(pJD)%VertLink(1,I))
        END DO

        CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,subnodes
      IF (E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)

      END IF
      END DO
    END IF
    END DO

    DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN
      nSIZE = E011ST(pJD)%Num

      DO I=1,nSIZE
      IF (E011ST(pJD)%RDVect(I).LT.0D0) THEN
        FX(E011ST(pJD)%VertLink(2,I)) = &
          MAX(FX(E011ST(pJD)%VertLink(2,I)),E011ST(pJD)%RDVect(I))
      END IF
      IF (E011ST(pJD)%RDVect(I).GT.0D0) THEN
        FX(E011ST(pJD)%VertLink(2,I)) = &
          MIN(FX(E011ST(pJD)%VertLink(2,I)),E011ST(pJD)%RDVect(I))
      END IF
      END DO

    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011_ABSMIN
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011Knpr(FX) !ok

  LOGICAL  FX(*)
  INTEGER I,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,subnodes
      IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        DO I=1,nSIZE
        E011ST(pJD)%SBVect(I) = FX(E011ST(pJD)%VertLink(1,I))
        END DO

        CALL SENDB_myMPI(E011ST(pJD)%SBVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,subnodes
      IF (E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        CALL RECVB_myMPI(E011ST(pJD)%RBVect,nSIZE,pJD)

      END IF
      END DO
    END IF
    END DO

    DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN
      nSIZE = E011ST(pJD)%Num

      DO I=1,nSIZE
      FX(E011ST(pJD)%VertLink(2,I)) = &
        (FX(E011ST(pJD)%VertLink(2,I)).OR.E011ST(pJD)%RBVect(I))
      END DO

    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011Knpr
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011Mean(FX,DM) !ok

  REAL*8  FX(*),DM(*)
  INTEGER I,pID,pJD,nSIZE,nEIGH
  REAL*8  DWEIGHT

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,subnodes
      IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        DO I=1,nSIZE
        E011ST(pJD)%SDVect(I) = &
          FX(E011ST(pJD)%VertLink(1,I))*DM(E011ST(pJD)%VertLink(1,I))
        E011ST(pJD)%SVVect(I) = REAL(DM(E011ST(pJD)%VertLink(1,I)))
        END DO

        CALL SENDD_myMPI(E011ST(pJD)%SDVect,nSIZE,pID)
        CALL SENDV_myMPI(E011ST(pJD)%SVVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,subnodes
      IF (E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        CALL RECVD_myMPI(E011ST(pJD)%RDVect,nSIZE,pJD)
        CALL RECVV_myMPI(E011ST(pJD)%RVVect,nSIZE,pJD)

      END IF
      END DO
    END IF
    END DO

    DO pJD=1,subnodes
    IF (E011ST(pJD)%Num.GT.0) THEN
      nSIZE = E011ST(pJD)%Num

      DO I=1,nSIZE
      DWEIGHT = &
        DBLE(DM(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RVVect(I))

      FX(E011ST(pJD)%VertLink(2,I)) = &
        (FX(E011ST(pJD)%VertLink(2,I))*DM(E011ST(pJD)%VertLink(2,I)) +&
        E011ST(pJD)%RDVect(I))/DWEIGHT

      END DO

    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011Mean
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CommSum(FX,ILEV) !ok

  REAL*8  FX(*)
  INTEGER ILEV
  INTEGER I,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SDVect(I) = FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(1,I))
        END DO

        CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%RDVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

    DO pJD=1,mg_mpi(ILEV)%NeighNum
    nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
    nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

    DO I=1,nSIZE
    !     if (myid.eq.1) WRITE(*,*) mg_mpi(ILEV)%parST(pJD)%RDVect(I),&
    !     mg_mpi(ILEV)%parST(pJD)%FaceLink(1,I),nSIZE
    FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) = &
      FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) + mg_mpi(ILEV)%parST(pJD)%RDVect(I)
    END DO

    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CommSum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E012_CommVal(FX,KX) !ok

  REAL*8  FX(4,*)
  INTEGER KX(5,*)
  INTEGER III,II,IJ,I,J,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    III=0
    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (pID.EQ.E012ST(pJD)%Neigh) THEN
        nSIZE = E012ST(pJD)%Num

        DO I=1,nSIZE
        III=III+1
        IJ = KX(1,III)
        DO J=1,4
        II = 4*(I-1) + J
        E012ST(pJD)%SDVect(II) = FX(J,IJ)
        END DO
        END DO

        !      write(*,*) "S", myid, pID
        CALL SENDD_myMPI(E012ST(pJD)%SDVect,4*nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (E012ST(pJD)%Num.GT.0) THEN
        nSIZE = E012ST(pJD)%Num
        nEIGH = E012ST(pJD)%Neigh

        !       write(*,*) "R", myid, nEIGH
        CALL RECVD_myMPI(E012ST(pJD)%RDVect,4*nSIZE,nEIGH)

      END IF
      END DO
    END IF
    END DO

    III=0
    DO pJD=1,subnodes
    IF (E012ST(pJD)%Num.GT.0) THEN
      nSIZE = E012ST(pJD)%Num

      DO I=1,nSIZE
      III=III+1
      IJ = KX(3,III)
      DO J=1,4
      II = 4*(I-1) + J
      FX(J,IJ) = E012ST(pJD)%RDVect(II)
      END DO
      END DO

    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E012_CommVal
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E012_CommSlice(FX,NX) !ok

  INTEGER NX
  INTEGER  FX(*)
  INTEGER III,II,IJ,I,J,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (pID.EQ.E012ST(pJD)%Neigh) THEN

        IF (myid.LT.pID) THEN
          E012ST(pJD)%SKVect(1) = FX(3   )
          E012ST(pJD)%SKVect(2) = FX(4   )
        ELSE
          E012ST(pJD)%SKVect(1) = FX(NX+1)
          E012ST(pJD)%SKVect(2) = FX(NX+2)
        END IF

        !      write(*,*) myid,pID,pJD
        !      IF (myid.eq.1) write(*,*) myid,pID,":S:",E012ST(pJD)%SKVect(1),E012ST(pJD)%SKVect(2)
        CALL SENDK_myMPI(E012ST(pJD)%SKVect,2,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (E012ST(pJD)%Num.GT.0) THEN
        nEIGH = E012ST(pJD)%Neigh

        CALL RECVK_myMPI(E012ST(pJD)%RKVect,2,nEIGH)

        !       IF (myid.eq.2) write(*,*) myid,nEIGH,":R:",E012ST(pJD)%RKVect(1),E012ST(pJD)%RKVect(2)
        IF (myid.LT.nEIGH) THEN
          FX(   1) = E012ST(pJD)%RKVect(1)
          FX(   2) = E012ST(pJD)%RKVect(2)
        ELSE
          FX(NX+3) = E012ST(pJD)%RKVect(1)
          FX(NX+4) = E012ST(pJD)%RKVect(2)
        END IF

      END IF
      END DO
    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E012_CommSlice
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E012_CommInt(FX,KX) !ok

  INTEGER FX(*)
  INTEGER KX(5,*)
  INTEGER III,II,IJ,I,J,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    III=0
    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (pID.EQ.E012ST(pJD)%Neigh) THEN
        nSIZE = E012ST(pJD)%Num

        DO I=1,nSIZE
        III=III+1
        IJ = KX(1,III)
        E012ST(pJD)%SKVect(I) = FX(IJ)
        END DO

        !      write(*,*) "S", myid, pID
        CALL SENDK_myMPI(E012ST(pJD)%SKVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (E012ST(pJD)%Num.GT.0) THEN
        nSIZE = E012ST(pJD)%Num
        nEIGH = E012ST(pJD)%Neigh

        !       write(*,*) "R", myid, nEIGH
        CALL RECVK_myMPI(E012ST(pJD)%RKVect,nSIZE,nEIGH)

      END IF
      END DO
    END IF
    END DO

    III=0
    DO pJD=1,subnodes
    IF (E012ST(pJD)%Num.GT.0) THEN
      nSIZE = E012ST(pJD)%Num

      DO I=1,nSIZE
      III=III+1
      IJ = KX(3,III)
      FX(IJ) = E012ST(pJD)%RKVect(I)
      END DO

    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E012_CommInt
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E012_CommSetUp(iSwitch,KX,II) !ok

  INTEGER iSwitch
  INTEGER KX(5,*)
  INTEGER I,II,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    IF (iSwitch.EQ.1) THEN
      II=0
      DO pID=1,mg_mpi(NLMINp)%NeighNum
      II = II + E012ST(pID)%Num
      END DO
      RETURN
    END IF

    II=0
    DO pID=1,subnodes
    !    write(*,*) "#",myid,E012ST(pID)%Neigh,E012ST(pID)%Num
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(NLMINp)%NeighNum
      IF (pID.EQ.E012ST(pJD)%Neigh) THEN
        nSIZE = E012ST(pJD)%Num

        DO I=1,nSIZE
        II = II + 1
        KX(1,II) = E012ST(pJD)%FaceLink(1,I)
        KX(2,II) = II
        KX(3,II) = E012ST(pJD)%FaceLink(2,I)
        KX(4,II) = II
        END DO

      END IF
      END DO
    END IF
    END DO

    CALL SORT2D(KX(1,1:II),KX(2,1:II),II)
    CALL SORT2D(KX(3,1:II),KX(4,1:II),II)
    DO I=1,II
    KX(5,I) = KX(1,I)
    KX(1,I) = I
    KX(3,I) = I
    END DO
    CALL SORT2D(KX(2,1:II),KX(1,1:II),II)
    CALL SORT2D(KX(4,1:II),KX(3,1:II),II)

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E012_CommSetUp
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CommSumHalf(FX,ILEV) !ok

  REAL*8  FX(*)
  INTEGER ILEV
  INTEGER I,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SDVect(I) = 0.5d0*FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(1,I))
        END DO

        CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%RDVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

    DO pJD=1,mg_mpi(ILEV)%NeighNum
    nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
    nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

    DO I=1,nSIZE
    FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) = &
      0.5d0*FX(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) + mg_mpi(ILEV)%parST(pJD)%RDVect(I)
    END DO
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CommSumHalf
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CommPres(FX,PX,ILEV) !ok

  REAL*8  FX(*),PX(*)
  INTEGER ILEV
  INTEGER I,pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SDVect(I)=&
          PX(mg_mpi(ILEV)%parST(pJD)%ElemLink(1,I))*mg_mpi(ILEV)%parST(pJD)%PE(I)
        END DO

        CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%RDVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

    DO pJD=1,mg_mpi(ILEV)%NeighNum
    nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
    nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

    DO I=1,nSIZE
    FX(mg_mpi(ILEV)%parST(pJD)%ElemLink(2,I)) = &
      FX(mg_mpi(ILEV)%parST(pJD)%ElemLink(2,I)) - mg_mpi(ILEV)%parST(pJD)%RDVect(I)
    END DO

    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CommPres
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011MAT(A,KLDA,NU) !ok

  REAL*4 A(*)
  INTEGER KLDA(*),NU
  INTEGER I,J
  INTEGER pID,pJD,nSIZE

  IF (myid.ne.MASTER) THEN

    DO I=1,NU
    E011_UE(I)=A(KLDA(I))
    ENDDO

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,subnodes
      IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        DO I=1,nSIZE
        E011ST(pJD)%SVVect(I) = E011_UE(E011ST(pJD)%VertLink(1,I))
        END DO

        CALL SENDV_myMPI(E011ST(pJD)%SVVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,subnodes
      IF (E011ST(pJD)%Num.GT.0) THEN

        nSIZE = E011ST(pJD)%Num
        CALL RECVV_myMPI(E011ST(pJD)%RVVect,nSIZE,pJD)

      END IF

      END DO
    END IF
    END DO

    DO pJD=1,subnodes

    IF (E011ST(pJD)%Num.GT.0) THEN

      nSIZE = E011ST(pJD)%Num

      DO I=1,nSIZE
      E011_UE(E011ST(pJD)%VertLink(2,I)) = &
        E011_UE(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RVVect(I)
      END DO

    END IF

    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011MAT
!
!
!
!SUBROUTINE E011MAT(A,KLDA,NU) !ok
!IMPLICIT NONE
!
!REAL*8 A(*)
!INTEGER KLDA(*),NU
!INTEGER I,J
!INTEGER pID,pJD,nSIZE
!
!IF (myid.ne.MASTER) THEN
!
! DO I=1,NU
!  E011_UE(I)=A(KLDA(I))
! ENDDO
!
! DO pID=1,subnodes
!  IF (myid.NE.pID) THEN
!   DO pJD=1,subnodes
!    IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!      E011ST(pJD)%SVVect(I) = E011_UE(E011ST(pJD)%VertLink(1,I))
!     END DO
!
!     CALL SENDV_myMPI(E011ST(pJD)%SVVect,nSIZE,pID)
!
!    END IF
!   END DO
!  ELSE
!   DO pJD=1,subnodes
!    IF (E011ST(pJD)%Num.GT.0) THEN
!
!     nSIZE = E011ST(pJD)%Num
!     CALL RECVV_myMPI(E011ST(pJD)%RVVect,nSIZE,pJD)
!
!    END IF
!
!   END DO
!  END IF
! END DO
!
! DO pJD=1,subnodes
!
!   IF (E011ST(pJD)%Num.GT.0) THEN
!
!     nSIZE = E011ST(pJD)%Num
!
!     DO I=1,nSIZE
!       E011_UE(E011ST(pJD)%VertLink(2,I)) = &
!       E011_UE(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RVVect(I)
!     END DO
!
!   END IF
!
! END DO
!
!END IF
!
!if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
!
!END SUBROUTINE E011MAT

! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE E011DMAT(A,KLDA,NU) !ok

  REAL*8 A(*)
  INTEGER KLDA(*),NU
  INTEGER I,J
  INTEGER pID,pJD,nSIZE

  IF (myid.ne.MASTER) THEN

    DO I=1,NU
    E011_UE(I)=REAL(A(KLDA(I)))
    ENDDO

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,subnodes
      IF (pID.EQ.pJD.AND.E011ST(pJD)%Num.GT.0) THEN
        nSIZE = E011ST(pJD)%Num

        DO I=1,nSIZE
        E011ST(pJD)%SVVect(I) = E011_UE(E011ST(pJD)%VertLink(1,I))
        END DO

        CALL SENDV_myMPI(E011ST(pJD)%SVVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,subnodes
      IF (E011ST(pJD)%Num.GT.0) THEN

        nSIZE = E011ST(pJD)%Num
        CALL RECVV_myMPI(E011ST(pJD)%RVVect,nSIZE,pJD)

      END IF

      END DO
    END IF
    END DO

    DO pJD=1,subnodes

    IF (E011ST(pJD)%Num.GT.0) THEN

      nSIZE = E011ST(pJD)%Num

      DO I=1,nSIZE
      E011_UE(E011ST(pJD)%VertLink(2,I)) = &
        E011_UE(E011ST(pJD)%VertLink(2,I)) + E011ST(pJD)%RVVect(I)
      END DO

    END IF

    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE E011DMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMMMAT(A,KLDA,NU,ILEV) !ok

  REAL*4 A(*)
  INTEGER KLDA(*),NU,ILEV
  INTEGER I,J
  INTEGER pID,pJD,nSIZE,nEIGH

  IF (myid.ne.MASTER) THEN

    DO I=1,NU
    mg_mpi(ILEV)%UE(I)=A(KLDA(I))
    ENDDO

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SVVect(I) = mg_mpi(ILEV)%UE(mg_mpi(ILEV)%parST(pJD)%FaceLink(1,I))
        END DO

        CALL SENDV_myMPI(mg_mpi(ILEV)%parST(pJD)%SVVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVV_myMPI(mg_mpi(ILEV)%parST(pJD)%RVVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

    DO pJD=1,mg_mpi(ILEV)%NeighNum
    nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
    nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

    DO I=1,nSIZE
    mg_mpi(ILEV)%UE(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) = &
      mg_mpi(ILEV)%UE(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I)) + mg_mpi(ILEV)%parST(pJD)%RVVect(I)
    END DO

    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE COMMMAT
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMMDISTRU(U,NU) !ok
  REAL*8 U(*)
  INTEGER NU
  INTEGER pNAT,NAT
  INTEGER I,pID

  NAT=NU
  IF (myid.eq.0) THEN

    DO pID=1,subnodes
    CALL RECVI_myMPI(pNAT,pID)
    CALL RECVD_myMPI(coarse%pDX,pNAT,pID)
    DO I=1,pNAT
    U(coarse%pFACELINK(pID,I))=coarse%pDX(I)
    END DO
    END DO

  ELSE

    DO I=1,NAT
    coarse%DX(I)=U(I)
    END DO
    CALL SENDI_myMPI(NAT,0)
    CALL SENDD_myMPI(coarse%DX,NAT,0)

  END IF

END SUBROUTINE COMMDISTRU
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------

SUBROUTINE COMMDEX(FX,FB,FD,NEQ,RHO) !ok

  REAL*8  FX(*),FB(*),FD(*),RHO,pRHO
  INTEGER NEQ,pNAT,NAT,NEL,pNEL
  INTEGER I,pID

  IF (VARIABLE.EQ.1) THEN
    NAT=NEQ
    IF (myid.eq.0) THEN

      !       Collect the RHS and X from the subdomains
      DO pID=1,subnodes
      CALL RECVI_myMPI(pNAT,pID)
      CALL RECVD_myMPI(coarse%pDX,pNAT,pID)
      DO I=1,pNAT
      FX(coarse%pFACELINK(pID,I))=coarse%pDX(I)
      FB(coarse%pFACELINK(pID,I))=coarse%pDX(I)
      END DO
      END DO

      !       Solve the system on coarse grid
!       CALL YEXU(FX,FB,FD,NAT,RHO)
      !  WRITE(*,*) "solve ...",RHO

      DO pID=1,subnodes
      CALL RECVI_myMPI(pNAT,pID)
      DO I=1,pNAT
      coarse%pDX(I)=FX(coarse%pFACELINK(pID,I))
      !          WRITE(*,*) "a",pDX(I)
      END DO
      !         PAUSE
      CALL SENDD_myMPI(coarse%pDX,pNAT,pID)
      CALL SENDDD_myMPI(RHO,pID)
      END DO
    ELSE
      DO I=1,NAT
      coarse%DX(I)=FX(I)
      END DO
      CALL SENDI_myMPI(NAT,0)
      CALL SENDD_myMPI(coarse%DX,NAT,0)

      CALL SENDI_myMPI(NAT,0)
      CALL RECVD_myMPI(coarse%DX,NAT,0)
      CALL RECVDD_myMPI(pRHO,0)
      DO I=1,NAT
      FX(I)=coarse%DX(I)
      END DO
      RHO = pRHO
    END IF

  ELSEIF (VARIABLE.EQ.2) THEN

    NEL=NEQ
    IF (myid.eq.0) THEN

      !       Collect the RHS and X from the subdomains
      DO pID=1,subnodes
      CALL RECVI_myMPI(pNEL,pID)
      CALL RECVD_myMPI(coarse%pDX,pNEL,pID)
      DO I=1,pNEL
      FX(coarse%pELEMLINK(pID,I))=coarse%pDX(I)
      FB(coarse%pELEMLINK(pID,I))=coarse%pDX(I)
      END DO
      END DO

      !       Solve the system on coarse grid
!       CALL YEXP(FX,FB,FD,NEL,RHO)
      !         WRITE(*,*) "solve ...",RHO

      DO pID=1,subnodes
      CALL RECVI_myMPI(pNEL,pID)
      DO I=1,pNEL
      coarse%pDX(I)=FX(coarse%pELEMLINK(pID,I))
      !          WRITE(*,*) "a",pDX(I)
      END DO
      !         PAUSE
      CALL SENDD_myMPI(coarse%pDX,pNEL,pID)
      CALL SENDDD_myMPI(RHO,pID)
      END DO
    ELSE
      DO I=1,NEL
      coarse%DX(I)=FX(I)
      END DO
      CALL SENDI_myMPI(NEL,0)
      CALL SENDD_myMPI(coarse%DX,NEL,0)

      CALL SENDI_myMPI(NEL,0)
      CALL RECVD_myMPI(coarse%DX,NEL,0)
      CALL RECVDD_myMPI(pRHO,0)
      DO I=1,NEL
      FX(I)=coarse%DX(I)
      END DO
      RHO = pRHO
    END IF

  END IF

  ! WRITE(*,*) "ok-",VARIABLE,myid
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE COMMDEX
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_MGComplete(IMGComplete) !ok

  INTEGER IMGComplete
  INTEGER iMG,piMG,pID

  ! if (VARIABLE.EQ.2) write(*,*) "complete?",myid,IMGComplete
  IF (myid.eq.MASTER) THEN

    iMG=1
    DO pID=1,subnodes
    CALL RECVI_myMPI(piMG,pID)
    !  write(*,*) pID,piMG
    iMG=iMG*pIMG
    END DO
    !write(*,*) "-----------"

    pIMG=iMG
    DO pID=1,subnodes
    CALL SENDI_myMPI(piMG,pID)
    END DO
    IMGComplete=iMG

  ELSE
    piMG=IMGComplete
    CALL SENDI_myMPI(piMG,0)
    CALL RECVI_myMPI(piMG,0)
    IMGComplete=piMG
  END IF

  ! if (VARIABLE.EQ.2) write(*,*) "complete!",myid,IMGComplete
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE COMM_MGComplete
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_Maximumn(DVAL,NN)
USE var_QuadScalar, ONLY :  myStat,iCommSwitch,myTimer
INTEGER NN
REAL*8 DVAL(NN)
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

if (iCommSwitch.le.2) CALL COMM_MaximumN_OLD(DVAL,NN)
if (iCommSwitch.ge.3) CALL COMM_MaximumN_NEW(DVAL,NN)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(2)%n = myTimer(2)%n + 1
myTimer(2)%t = myTimer(2)%t + dble(tt1-tt0)

END SUBROUTINE COMM_Maximumn
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_MaximumX(value) !ok
USE var_QuadScalar, only :myStat,iCommSwitch,myTimer
REAL*8 value
INTEGER NN
REAL*8 DVAL(1)
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

NN = 1
dval(1) = value 
if (iCommSwitch.le.2) CALL COMM_MaximumN_OLD(DVAL,NN)
if (iCommSwitch.ge.3) CALL COMM_MaximumN_NEW(DVAL,NN)

value = DVAL(1)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(6)%n = myTimer(6)%n + 1
myTimer(6)%t = myTimer(6)%t + dble(tt1-tt0)

END SUBROUTINE COMM_MaximumX
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_Maximum(value) !ok
USE var_QuadScalar, only :myStat,iCommSwitch,myTimer
REAL*8 value
INTEGER NN
REAL*8 DVAL(1)
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

NN = 1
dval(1) = value 
if (iCommSwitch.le.2) CALL COMM_MaximumN_OLD(DVAL,NN)
if (iCommSwitch.ge.3) CALL COMM_MaximumN_NEW(DVAL,NN)

value = DVAL(1)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(1)%n = myTimer(1)%n + 1
myTimer(1)%t = myTimer(1)%t + dble(tt1-tt0)

END SUBROUTINE COMM_Maximum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_Maximum_old(value) !ok

  REAL*8  value
  REAL*8  Val,pVal
  INTEGER pID

  ! if (VARIABLE.EQ.2) write(*,*) "complete?",myid,IMGComplete
  IF (myid.eq.MASTER) THEN

    Val=-1D99
    DO pID=1,subnodes
    CALL RECVDD_myMPI(pVal,pID)
    !  write(*,*) pID,piMG
    Val=MAX(Val,pVal)
    END DO
    !write(*,*) "-----------"

    pVal=Val
    DO pID=1,subnodes
    CALL SENDDD_myMPI(pVal,pID)
    END DO
    value=Val

  ELSE
    pVal=value
    CALL SENDDD_myMPI(pVal,0)
    CALL RECVDD_myMPI(pVal,0)
    value=pVal
  END IF

  ! if (VARIABLE.EQ.2) write(*,*) "complete!",myid,IMGComplete
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE COMM_Maximum_old
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_Minimum(value) !ok

REAL*8  value
REAL*8  Val,pVal
INTEGER pID

! if (VARIABLE.EQ.2) write(*,*) "complete?",myid,IMGComplete
IF (myid.eq.MASTER) THEN

 Val=+1D99
 DO pID=1,subnodes
  CALL RECVDD_myMPI(pVal,pID)
!  write(*,*) pID,piMG
  Val=MIN(Val,pVal)
 END DO
 !write(*,*) "-----------"

 pVal=Val
 DO pID=1,subnodes
  CALL SENDDD_myMPI(pVal,pID)
 END DO
 value=Val

ELSE
  pVal=value
  CALL SENDDD_myMPI(pVal,0)
  CALL RECVDD_myMPI(pVal,0)
  value=pVal
END IF

! if (VARIABLE.EQ.2) write(*,*) "complete!",myid,IMGComplete
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE COMM_Minimum
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_Minimumn(DVAL,NN)
  INTEGER NN
  REAL*8 DVAL(NN)
  REAL*8, ALLOCATABLE :: pVal(:)
  INTEGER pID,i

  ALLOCATE (pVal(NN))

  IF (myid.eq.MASTER) THEN
    DVAL = +1d30
    DO pID=1,subnodes
    CALL RECVD_myMPI(pVAL,NN,pID)
    DO i=1,NN
    DVal(i) = MIN(DVAL(i),pVal(i))
    END DO
    END DO

    pVAL = DVAL
    DO pID=1,subnodes
    CALL SENDD_myMPI(pVAL,NN,pID)
    END DO

  ELSE
    pVAL=DVAL
    CALL SENDD_myMPI(pVAL,NN,0)
    CALL RECVD_myMPI(pVAL,NN,0)
    DVAL = pVal
  END IF

  DEALLOCATE (pVal)

END SUBROUTINE COMM_Minimumn
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
! SUBROUTINE COMM_NLComplete(INLComplete)
! 
!   INTEGER INLComplete
!   INTEGER iNL,piNL,pID
! 
!   ! write(*,*) "complete?",myid,INLComplete
!   IF (myid.eq.MASTER) THEN
! 
!     iNL=1
!     DO pID=1,subnodes
!     CALL RECVI_myMPI(piNL,pID)
!     iNL=iNL*pINL
!     END DO
! 
!     pINL=iNL
!     DO pID=1,subnodes
!     CALL SENDI_myMPI(piNL,pID)
!     END DO
!     INLComplete=iNL
! 
!   ELSE
!     piNL=INLComplete
!     CALL SENDI_myMPI(piNL,0)
!     CALL RECVI_myMPI(piNL,0)
!     INLComplete=piNL
!   END IF
! 
!   ! write(*,*) "complete!",myid,INLComplete
!   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
! 
! END SUBROUTINE COMM_NLComplete
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_NLComplete(INLComplete)
USE var_QuadScalar, ONLY :  myStat,iCommSwitch,myTimer
INTEGER INLComplete
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

if (iCommSwitch.le.2) CALL COMM_INLN_OLD(INLComplete)
if (iCommSwitch.ge.3) CALL COMM_INLN_NEW(INLComplete)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(5)%n = myTimer(5)%n + 1
myTimer(5)%t = myTimer(5)%t + dble(tt1-tt0)

END SUBROUTINE COMM_NLComplete
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_SUMMN(DVAL,NN)
USE var_QuadScalar, ONLY :  myStat,iCommSwitch,myTimer
INTEGER NN
REAL*8 DVAL(NN)
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

if (iCommSwitch.le.2) CALL COMM_SUMMN_OLD(DVAL,NN)
if (iCommSwitch.ge.3) CALL COMM_SUMMN_NEW(DVAL,NN)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(4)%n = myTimer(4)%n + 1
myTimer(4)%t = myTimer(4)%t + dble(tt1-tt0)


END SUBROUTINE COMM_SUMMN
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE COMM_SUMM(value)
USE var_QuadScalar, ONLY :  myStat,iCommSwitch,myTimer
REAL*8 value
INTEGER NN
REAL*8 DVAL(1)
REAL*4  tt1,tt0

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt0)

NN = 1
dval(1) = value 
if (iCommSwitch.le.2) CALL COMM_SUMMN_OLD(DVAL,NN)
if (iCommSwitch.ge.3) CALL COMM_SUMMN_NEW(DVAL,NN)

value = DVAL(1)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL ztime(tt1)
myStat%tCommS = myStat%tCommS + dble(tt1-tt0)
myTimer(3)%n = myTimer(3)%n + 1
myTimer(3)%t = myTimer(3)%t + dble(tt1-tt0)

END SUBROUTINE COMM_SUMM
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CommBMul(PX,ILEV) !ok

  REAL*8  PX(*)
  INTEGER ILEV
  INTEGER I
  INTEGER pID,pJD,nSIZE,nEIGH

  IF (myid.ne.0) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SDVect(I) = PX(mg_mpi(ILEV)%parST(pJD)%ElemLink(2,I))
        END DO

        CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%RDVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

  END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CommBMul
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE CommBMMul(PX,ILEV)

  REAL*8  PX(*)
  INTEGER ILEV
  INTEGER I
  INTEGER pID,pJD,nSIZE,nEIGH

  IF (myid.ne.0) THEN

    DO pID=1,subnodes
    IF (myid.NE.pID) THEN
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      IF (pID.EQ.mg_mpi(ILEV)%parST(pJD)%Neigh) THEN
        nSIZE = mg_mpi(ILEV)%parST(pJD)%Num

        DO I=1,nSIZE
        mg_mpi(ILEV)%parST(pJD)%SDVect(I) = PX(mg_mpi(ILEV)%parST(pJD)%FaceLink(2,I))
        END DO

        CALL SENDD_myMPI(mg_mpi(ILEV)%parST(pJD)%SDVect,nSIZE,pID)

      END IF
      END DO
    ELSE
      DO pJD=1,mg_mpi(ILEV)%NeighNum
      nSIZE = mg_mpi(ILEV)%parST(pJD)%Num
      nEIGH = mg_mpi(ILEV)%parST(pJD)%Neigh

      CALL RECVD_myMPI(mg_mpi(ILEV)%parST(pJD)%RDVect,nSIZE,nEIGH)

      END DO
    END IF
    END DO

  END IF

  !  DO I=1,NParNodes(ILEV)
  !   mg_mpi(ILEV)%parST(1)%SDVect(I)=PX(mg_mpi(ILEV)%parST(1)%FaceLink(2,I))
  !  END DO
  ! 
  !  IF (myid.eq.1) THEN
  !   CALL SENDD_myMPI(mg_mpi(ILEV)%parST(1)%SDVect,NParNodes(ILEV),2)
  !   CALL RECVD_myMPI(mg_mpi(ILEV)%parST(1)%RDVect,NParNodes(ILEV),2)
  !  END IF
  ! 
  !  IF (myid.eq.2) THEN
  !   CALL RECVD_myMPI(mg_mpi(ILEV)%parST(1)%RDVect,NParNodes(ILEV),1)
  !   CALL SENDD_myMPI(mg_mpi(ILEV)%parST(1)%SDVect,NParNodes(ILEV),1)
  !  END IF
  ! 
  ! END IF

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE CommBMMul
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
! SUBROUTINE CommSortPre(KTR1,KTR2,NEQ,IPAR,ILEV)
!   INTEGER KTR1(*),KTR2(*)
!   INTEGER IPAR,IEQ,NEQ,ILEV
!   INTEGER pID,I,IAUX
! 
!   DO pID=1,mg_mpi(ILEV)%NeighNum
!   DO I=1,mg_mpi(ILEV)%parST(pID)%Num
!   IAUX = KTR1(mg_mpi(ILEV)%parST(pID)%FaceLink(1,I))
!   mg_mpi(ILEV)%parST(pID)%FaceLin1(1,I) = IAUX
!   IAUX = KTR1(mg_mpi(ILEV)%parST(pID)%FaceLink(2,I))
!   mg_mpi(ILEV)%parST(pID)%FaceLin1(2,I) = IAUX
! 
!   IAUX = KTR2(mg_mpi(ILEV)%parST(pID)%FaceLink(1,I))
!   mg_mpi(ILEV)%parST(pID)%FaceLin2(1,I) = IAUX
!   IAUX = KTR2(mg_mpi(ILEV)%parST(pID)%FaceLink(2,I))
!   mg_mpi(ILEV)%parST(pID)%FaceLin2(2,I) = IAUX
!   END DO
!   END DO
! 
!   ! IF (IPAR.EQ.1) THEN
!   !  DO IEQ=1,NEQ
!   !   DD(IEQ)=DX(KTR1(IEQ))
!   !  END DO
!   ! ELSE
!   !  DO IEQ=1,NEQ
!   !   DD(IEQ)=DX(KTR2(IEQ))
!   !  END DO
!   ! ENDIF
! 
! END SUBROUTINE CommSortPre
! ----------------------------------------------
! ----------------------------------------------
! ----------------------------------------------
SUBROUTINE myMPI_Barrier()

  if (myid.ne.MASTER) CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)

END SUBROUTINE myMPI_Barrier












! --------------- BASIC ROUTINES FOR DATA EXCHANGE ------------
! -------------------------------------------------------------
! ---------------------- BARRIER ------------------------------
SUBROUTINE Barrier_myMPI()

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE Barrier_myMPI

! -------------------------------------------------------------
! -------------------- BCAST of INTEGERs- ---------------------
SUBROUTINE ShareValueK_myMPI(KDATA,nData,iProc)

  INTEGER KDATA(*)
  INTEGER iProc,nData

  CALL MPI_BCAST(KDATA,nData,MPI_INTEGER,iProc,MPI_COMM_WORLD,IERR)

END SUBROUTINE ShareValueK_myMPI
! -------------------------------------------------------------
! ---------------------- BCAST of DOUBLEs ---------------------
SUBROUTINE ShareValueD_myMPI(DDATA,nData,iProc)

  REAL*8 DDATA(*)
  INTEGER iProc,nData

  CALL MPI_BCAST(DDATA,nData,MPI_DOUBLE_PRECISION,iProc,MPI_COMM_WORLD,IERR)

END SUBROUTINE ShareValueD_myMPI
! -------------------------------------------------------------
! ---------------------- BCAST of SINGLEs ---------------------
SUBROUTINE ShareValueV_myMPI(VDATA,nData,iProc)

  REAL*4  VDATA(*)
  INTEGER iProc,nData

  CALL MPI_BCAST(VDATA,nData,MPI_REAL,iProc,MPI_COMM_WORLD,IERR)

END SUBROUTINE ShareValueV_myMPI
! -------------------------------------------------------------
! ---------------------- BCAST of CHARs -----------------------
SUBROUTINE ShareValueC_myMPI(CDATA,nData,iProc)

  CHARACTER CDATA(*)
  INTEGER iProc,nData

  CALL MPI_BCAST(CDATA,nData,MPI_CHARACTER,iProc,MPI_COMM_WORLD,IERR)

END SUBROUTINE ShareValueC_myMPI



! -------------------------------------------------------------
! ------------------ SEND/RECV of LOGICALs --------------------
SUBROUTINE SENDB_myMPI(DDATA,N,DEST)
implicit none
LOGICAL DDATA(*)
INTEGER N
integer DEST
INTEGER ITAG

ITAG=7

CALL MPI_SSEND(DDATA,N,MPI_LOGICAL,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDB_myMPI
! -------------------------------------------------------------
SUBROUTINE RECVB_myMPI(DDATA,N,SRC)
implicit none
LOGICAL DDATA(*)
INTEGER N
integer SRC
INTEGER ITAG
integer st(MPI_STATUS_SIZE)

ITAG=7

CALL MPI_RECV(DDATA,N,MPI_LOGICAL,SRC,ITAG,MPI_COMM_WORLD,ST,IERR)

END SUBROUTINE RECVB_myMPI
! ------------------ SEND/RECV of DOUBLEs ---------------------
SUBROUTINE SENDD_myMPI(DDATA,N,DEST)
implicit none
REAL*8  DDATA(*)
INTEGER N
integer DEST
INTEGER ITAG

ITAG=2

CALL MPI_SSEND(DDATA,N,MPI_DOUBLE_PRECISION,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDD_myMPI



SUBROUTINE RECVD_myMPI(DDATA,N,SRC)
implicit none
REAL*8 DDATA(*)
INTEGER N
integer SRC
INTEGER ITAG
integer st(MPI_STATUS_SIZE)

ITAG=2
CALL MPI_RECV(DDATA,N,MPI_DOUBLE_PRECISION,SRC,ITAG,MPI_COMM_WORLD,st,IERR)

END SUBROUTINE RECVD_myMPI



!     ------------ SEND/RECV of REALs ---------------
SUBROUTINE SENDV_myMPI(DDATA,N,DEST)
implicit none
REAL*4  DDATA(*)
INTEGER N
integer DEST
INTEGER ITAG

ITAG=3

CALL MPI_SSEND(DDATA,N,MPI_REAL,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDV_myMPI



SUBROUTINE RECVV_myMPI(DDATA,N,SRC)
implicit none
REAL*4 DDATA(*)
INTEGER N
integer SRC
INTEGER ITAG
integer St(MPI_STATUS_SIZE)

ITAG=3

CALL MPI_RECV(DDATA,N,MPI_REAL,SRC,ITAG,MPI_COMM_WORLD,St,IERR)

END SUBROUTINE RECVV_myMPI



!     ------------ SEND/RECV of INTEGERs ---------------
SUBROUTINE SENDK_myMPI(DDATA,N,DEST)
implicit none
INTEGER DDATA(*)
INTEGER N
integer DEST
INTEGER ITAG

ITAG=4

CALL MPI_SSEND(DDATA,N,MPI_INTEGER,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDK_myMPI



SUBROUTINE RECVK_myMPI(DDATA,N,SRC)
implicit none
INTEGER DDATA(*)
INTEGER N
integer SRC
INTEGER ITAG
integer St(MPI_STATUS_SIZE)

ITAG=4

CALL MPI_RECV(DDATA,N,MPI_INTEGER,SRC,ITAG,MPI_COMM_WORLD,ST,IERR)

END SUBROUTINE RECVK_myMPI



!     ------------ SEND/RECV of INTEGERs ---------------
SUBROUTINE SENDI_myMPI(DDATA,DEST)
implicit none
INTEGER DDATA
INTEGER DEST
INTEGER ITAG
INTEGER SBUFF

ITAG=0
CALL MPI_SSEND(DDATA,1,MPI_INTEGER,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDI_myMPI



SUBROUTINE RECVI_myMPI(DDATA,SRC)
implicit none
INTEGER DDATA
INTEGER SRC
!locals
INTEGER ITAG
integer St(MPI_STATUS_SIZE)
integer ncount

ncount=1
itag=0
CALL MPI_RECV(DDATA,1,MPI_INTEGER,SRC,ITAG,MPI_COMM_WORLD,St,IERR)

END SUBROUTINE RECVI_myMPI


! ------------------ SEND/RECV of DOUBLEs ---------------------
SUBROUTINE SENDDD_myMPI(DDATA,DEST)
implicit none
REAL*8  DDATA
INTEGER DEST
INTEGER ITAG
REAL*8  SBUFF

ITAG=6

SBUFF=DDATA
CALL MPI_SSEND(SBUFF,1,MPI_DOUBLE_PRECISION,DEST,ITAG,MPI_COMM_WORLD,IERR)

END SUBROUTINE SENDDD_myMPI



SUBROUTINE RECVDD_myMPI(DDATA,SRC)
implicit none
REAL*8  DDATA
INTEGER SRC
INTEGER ITAG
integer ST(MPI_STATUS_SIZE)
REAL*8  RBUFF

ITAG=6

CALL MPI_RECV(RBUFF,1,MPI_DOUBLE_PRECISION,SRC,ITAG,MPI_COMM_WORLD,ST,IERR)
DDATA=RBUFF
END SUBROUTINE RECVDD_myMPI


SUBROUTINE SORT1D(LW,N)
  INTEGER LW(N),LWA
  INTEGER I,J,N

  DO I=2,N
  DO J=N,I,-1
  IF (LW(J).LT.LW(J-1)) THEN
    LWA     = LW(J)
    LW(J)   = LW(J-1)
    LW(J-1) = LWA
  END IF
  END DO
  END DO

END SUBROUTINE SORT1D


SUBROUTINE SORT2D(LW,KW,N)
  INTEGER LW(N),KW(N),LWA,KWA
  INTEGER I,J,N

  DO I=2,N
  DO J=N,I,-1
  IF (LW(J).LT.LW(J-1)) THEN
    LWA     = LW(J)
    KWA     = KW(J)
    LW(J)   = LW(J-1)
    KW(J)   = KW(J-1)
    LW(J-1) = LWA
    KW(J-1) = KWA
  END IF
  END DO
  END DO

  ! DO I=1,N
  !  write(*,*) I,LW(I),KW(I)
  ! END DO

END SUBROUTINE SORT2D


SUBROUTINE SORTALL(LW,KW,MW,x,y,z,N)
  INTEGER LW(N),KW(N),MW(N),LWA,KWA,MWA
  REAL*8  x(N),y(N),z(N),xd,yd,zd
  INTEGER I,J,N

  DO I=2,N
  DO J=N,I,-1
  IF (LW(J).LT.LW(J-1)) THEN
    LWA     = LW(J)
    KWA     = KW(J)
    MWA     = MW(J)
    xd      = x(J)
    yd      = y(J)
    zd      = z(J)
    LW(J)   = LW(J-1)
    KW(J)   = KW(J-1)
    MW(J)   = MW(J-1)
    x(J)    = x(J-1)
    y(J)    = y(J-1)
    z(J)    = z(J-1)
    LW(J-1) = LWA
    KW(J-1) = KWA
    MW(J-1) = MWA
    x(J-1)  = xd
    y(J-1)  = yd
    z(J-1)  = zd
  END IF
  END DO
  END DO

END SUBROUTINE SORTALL
!-------------------------------------------------------------
!
!
!-------------------------------------------------------------
SUBROUTINE Reduce_myMPI(localMax, totalMax)

  REAL*8 :: localMax 
  REAL*8 :: totalMax 
  integer :: error_indicator

  totalMax = 0.0
  call MPI_Reduce(localMax, totalMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, error_indicator)
  CALL MPI_BARRIER(MPI_COMM_WORLD, error_indicator)

  if (myid == 0) then
    write(*,'(A,E12.6)')'Total Max value: ', totalMax
  end if

END SUBROUTINE Reduce_myMPI
!-------------------------------------------------------------
!
!
!-------------------------------------------------------------
SUBROUTINE Sum_myMPI(localMax, totalMax)

  REAL*8 :: localMax 
  REAL*8 :: totalMax 
  integer :: error_indicator

  totalMax = 0.0
  call MPI_Reduce(localMax, totalMax, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, error_indicator)
  CALL MPI_BARRIER(MPI_COMM_WORLD, error_indicator)

  if (myid == 0) then
    print *,'Global lubrication', totalMax
  end if

END SUBROUTINE Sum_myMPI
!-------------------------------------------------------------
!
!
!-------------------------------------------------------------
SUBROUTINE SynchronizeValue_myMPI(localMax, totalMax)
  implicit none

  ! Parameters
  REAL*8 :: localMax 
  REAL*8 :: totalMax 

  ! Local variables
  integer :: error_indicator
  integer :: root = 0

  totalMax = 0.0

  ! Reduce to find the maximum
  call MPI_Reduce(localMax, totalMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, MPI_COMM_WORLD, error_indicator)

  ! Now bcast the maximum to synchronize the value among the processes
  call MPI_Bcast(totalMax, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, error_indicator) 

  CALL MPI_BARRIER(MPI_COMM_WORLD, error_indicator)

  if (myid == 0) then
    write(*,'(A,E12.6)')'Total Max value: ', totalMax
  end if

END SUBROUTINE SynchronizeValue_myMPI

END MODULE
