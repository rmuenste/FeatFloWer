      SUBROUTINE AddDiffPrec(DMP,RHS,VAL,NEL,dF)
      USE PP3D_MPI, ONLY:myid
      implicit none
C
      INTEGER NEL
      REAL*8 DMP(4,4,*),RHS(*),VAL(*),dF
      !
      integer i4,j4,i,j,iel

!       if (myid.eq.1) WRITE(*,*) 'RHS: ',RHS(1:4*NEL)

      DO iel=1,nel
       DO i4=1,4
        i=4*(iel-1)+i4
        
        DO j4=1,4
         j=4*(iel-1)+j4
         VAL(i) = VAL(i) + dF*DMP(i4,j4,iel)*RHS(j)
        END DO
       END DO
      END DO
      
!       if (myid.eq.1) WRITE(*,*) 'VAL: ',VAL(1:4*NEL)
      END
************************************************************************
      SUBROUTINE Get_BxiMP(DMP,BX,BY,BZ,iBX,iBY,iBZ,
     *           KLD,KCOL,NEL,NDOF,dF)
************************************************************************
      USE PP3D_MPI, ONLY:myid
      implicit none
C
      REAL*8 DMP(4,4,*)
      REAL*8 BX(*),BY(*),BZ(*),iBX(*),iBY(*),iBZ(*)
      real*8 dF
      integer KCOL(*),KLD(*)
      integer NEL,NDOF
C      
      integer idof,iel,jel,iEntry,jEntry,iaux,jaux,iComp,jComp
      real*8 dentry(3),MP_Entry
      
      
      DO idof=1,ndof
       do iEntry=kld(idof),kld(idof+1)-1
       
        iaux = kcol(iEntry)
        iComp = mod(iaux-1,4) + 1
        iel = (iaux - iComp)/4 + 1
        
        dentry = 0d0
!         dentry = [BX(iEntry),BY(iEntry),BZ(iEntry)]

        do jEntry=kld(idof),kld(idof+1)-1
        
         jaux = kcol(jEntry)
         jComp = mod(jaux-1,4) + 1
         jel = (jaux - jComp)/4 + 1
         
         if (jel.eq.iel) THEN
!           if (myid.eq.1) write(*,'(a,3I12)') 'entry: ',iComp,jComp,iel
         
          MP_Entry = DMP(jComp,iComp,iel)
          dentry(1) = dentry(1) + MP_entry*BX(jEntry)
          dentry(2) = dentry(2) + MP_entry*BY(jEntry)
          dentry(3) = dentry(3) + MP_entry*BZ(jEntry)
         end if
         
        end do
        
        iBX(iEntry) = dF*dentry(1)
        iBY(iEntry) = dF*dentry(2)
        iBZ(iEntry) = dF*dentry(3)
        
!         if (myid.eq.1) write(*,'(a,3es12.4)') 'dentry: ',dentry 
       
       end do
      ENDdo
      
      END 
C      
************************************************************************
      SUBROUTINE InvertP1Mat(DA,iDA,nel)
************************************************************************
      USE PP3D_MPI, ONLY:myid
      implicit none
      
      REAL*8 DA(4,4,*),iDA(4,4,*)
      integer nel
      
      integer i,errorflag
   
      DO I=1,NEL
       CALL FINDInv(DA(:,:,I), iDA(:,:,I), 4, errorflag)
       if (errorflag.ne.0) write(*,*) 'error:', myid, errorflag
!        if (myid.eq.1) write(*,'(4ES12.4)') DA(1,:,I)
!        if (myid.eq.1) write(*,'(4ES12.4)') DA(2,:,I)
!        if (myid.eq.1) write(*,'(4ES12.4)') DA(3,:,I)
!        if (myid.eq.1) write(*,'(4ES12.4)') DA(4,:,I)
!        if (myid.eq.1) write(*,'(4ES12.4)') 
      END DO
!        pause
      END 
C
************************************************************************
      SUBROUTINE Build_MMatP1(DA,KVERT,KAREA,KEDGE,
     *           DCORVG,ICUB,ELE)
************************************************************************
*     Discrete mass-matrix operator: P1 elements
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      REAL*8 DAUX
      DIMENSION DA(4,4,*)
      DIMENSION DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      SAVE
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      BDER(4)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      DAUX = 1d0
C
      ILINT=0
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
!       write(*,*) KDFL(1),KDFL(2),KDFL(3),KDFL(4)
!       write(*,*) KDFG(1),KDFG(2),KDFG(3),KDFG(4)
!       stop
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DX0I = 0d0
      DY0I = 0d0
      DZ0I = 0d0
      DO 120 IVE=1,NVE
      IP=KVERT(IVE,IEL)
      KVE(IVE)=IP
      DX(IVE)=DCORVG(1,IP)
      DY(IVE)=DCORVG(2,IP)
      DZ(IVE)=DCORVG(3,IP)
      DX0I = DX0I + 0.125d0*DCORVG(1,IP)
      DY0I = DY0I + 0.125d0*DCORVG(2,IP)
      DZ0I = DZ0I + 0.125d0*DCORVG(3,IP)
120   CONTINUE
C
      DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
      DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
      DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
      DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
      DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
      DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
      DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
      DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
      DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
      DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
      DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
      DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
      DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
      DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
      DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
      DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
      DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
      DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
      DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
      DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
      DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
      DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
      DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
      DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
      DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
      DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
      DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
      DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
      DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
      DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
      DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
!     Computation of the cartesian coordiante of the cubature point
      XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
      YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
      ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
C
      CALL ELE(XX-DX0I,YY-DY0I,ZZ-DZ0I,0)
C     
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HBASJ*HBASJ
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH=HBASJ*HBASI
        ENDIF
        DA(JDOFE,IDOFE,IEL)=DA(JDOFE,IDOFE,IEL) + DAUX*OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE BuildMRhoMat(DENS,DA,NA,KCOLA,KLDA,KVERT,KAREA,
     *                  KEDGE,DCORVG,ELE)
************************************************************************
*     Discrete diffusion operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DENS(*),DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      SAVE
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
      BDER(2:4)=.FALSE.
C     
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
!       ILINT=0!KINT(IEL)
      DENSITY = DENS(IEL)
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ILINT.eq.2) THEN ! Q2
      DO JDOFE=1,27
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q2(JDFL,4,ICUBP)
      END DO
      END IF
      IF (Transform%ILINT.eq.1) THEN ! Q1
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
      END IF
C      
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=DENSITY*(HBASJ*HBASJ)
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH=DENSITY*(HBASJ*HBASI)
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)+DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE BuildMMat(DA,NA,KCOLA,KLDA,KVERT,KAREA,
     *                  KEDGE,DCORVG,ELE)
************************************************************************
*     Discrete diffusion operator: Q1 elements
*-----------------------------------------------------------------------
      USE PP3D_MPI, ONLY:myid
      USE var_QuadScalar, ONLY : transform
C     
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION DA(*)
      DIMENSION KCOLA(*),KLDA(*),DCORVG(NNDIM,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q2(27,4,NNCUBP),DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
C     --------------------------------------------------------------------------
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      SAVE
C
      DO 1 I= 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      ICUB=9
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
       CALL E013A(XI1,XI2,XI3,DHELP_Q2,ICUBP)
      END DO
C
      BDER(2)=.FALSE.
      BDER(3)=.FALSE.
      BDER(4)=.FALSE.
C     
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      IF (Transform%ILINT.eq.2) THEN ! Q2
      DO JDOFE=1,27
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q2(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q2(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q2(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q2(JDFL,4,ICUBP)
      END DO
      END IF
      IF (Transform%ILINT.eq.1) THEN ! Q1
      DO JDOFE=1,8
       JDFL=KDFL(JDOFE)
       JDFG=KDFG(JDOFE)
       DPP(:) = DCORVG(:,JDFG)
       DJAC(1,1)= DJAC(1,1) +  DPP(1)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(2,1)= DJAC(2,1) +  DPP(2)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(3,1)= DJAC(3,1) +  DPP(3)*DHELP_Q1(JDFL,2,ICUBP)
       DJAC(1,2)= DJAC(1,2) +  DPP(1)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(2,2)= DJAC(2,2) +  DPP(2)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(3,2)= DJAC(3,2) +  DPP(3)*DHELP_Q1(JDFL,3,ICUBP)
       DJAC(1,3)= DJAC(1,3) +  DPP(1)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(2,3)= DJAC(2,3) +  DPP(2)*DHELP_Q1(JDFL,4,ICUBP)
       DJAC(3,3)= DJAC(3,3) +  DPP(3)*DHELP_Q1(JDFL,4,ICUBP)
      END DO
      END IF
C      
      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *     -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *     +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ=DBAS(1,JDOFEH,1)
C
       DO 240 IDOFE=1,IDFL
        IF (IDOFE.EQ.JDOFE) THEN
         AH=HBASJ*HBASJ
        ELSE
         IDOFEH=KDFL(IDOFE)
         HBASI=DBAS(1,IDOFEH,1)
         AH=HBASJ*HBASI
        ENDIF
        DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230   CONTINUE
C
200   CONTINUE
C
      DO 300 JDOFE=1,IDFL
      DO 300 IDOFE=1,IDFL
        IA    =KENTRY(JDOFE,IDOFE)
        DA(IA)=DA(IA)+DENTRY(JDOFE,IDOFE)
300   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
C
      SUBROUTINE Q1ElemShapeFunction(X1,X2,X3,DPHI)
      REAL*8 X1,X2,X3,DPHI(8)
      REAL*8 :: Q8=0.125D0

      DPHI(1)=Q8*(1D0-X1)*(1D0-X2)*(1D0-X3)
      DPHI(2)=Q8*(1D0+X1)*(1D0-X2)*(1D0-X3)
      DPHI(3)=Q8*(1D0+X1)*(1D0+X2)*(1D0-X3)
      DPHI(4)=Q8*(1D0-X1)*(1D0+X2)*(1D0-X3)
      DPHI(5)=Q8*(1D0-X1)*(1D0-X2)*(1D0+X3)
      DPHI(6)=Q8*(1D0+X1)*(1D0-X2)*(1D0+X3)
      DPHI(7)=Q8*(1D0+X1)*(1D0+X2)*(1D0+X3)
      DPHI(8)=Q8*(1D0-X1)*(1D0+X2)*(1D0+X3)

      END
