      SUBROUTINE DIFFQ2_Elem(DCORVG,ELE)
      use UMFPackSolver, only : myUmfPack_Factorize,myUmfPack_Solve,myUmfPack_Free
      use types, only : tmatrix
      
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

      CHARACTER SUB*6,FMT*15,CPARAM*120

      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)

      INTEGER   KVERT(nnve),KAREA(nnae),KEDGE(nnee)
      DATA      KVERT/1,2,3,4,5,6,7,8/
      DATA      KEDGE/1,2,3,4,5,6,7,8,9,10,11,12/
      DATA      KAREA/1,2,3,4,5,6/
      
      REAL*8    DCORVG(nndim,nnve)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      REAL*8    DQ1BAS(8)
      
      TYPE(tMatrix) dMat
      REAL*8, ALLOCATABLE :: SOL(:),RHS(:)

!     --------------------------- Transformation -------------------------------
      REAL*8    DHELP_Q1(8,4,NNCUBP)
      REAL*8    DPP(3)
!     --------------------------------------------------------------------------
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                      DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                      IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                      NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL

      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
                     IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
                     ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
                     INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
                     ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
                     IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

      SAVE

      DO  I= 1,NNDER
       BDER(I)=.FALSE.
      end do

      DO I=1,4
        BDER(I)=.TRUE.
      end do

      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=27
      write(*,*) 'IELTYP',IELTYP,idfl

      ICUB=9
      CALL CB3H(ICUB)

      DO ICUBP=1,NCUBP
       XI1=DXI(ICUBP,1)
       XI2=DXI(ICUBP,2)
       XI3=DXI(ICUBP,3)
       CALL E011A(XI1,XI2,XI3,DHELP_Q1,ICUBP)
      END DO
      
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
      
      iel = 1
      
      dVisc = 1d0
     
! Evaluation of coordinates of the vertices
      DO IVE=1,NVE
       DX(IVE)=DCORVG(1,ive)
       DY(IVE)=DCORVG(2,ive)
       DZ(IVE)=DCORVG(3,ive)
      END DO

!  Loop over all cubature points
      DO ICUBP=1,NCUBP
      
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
      
! Jacobian of the (trilinear,triquadratic,or simple) mapping onto the reference element
      DJAC=0d0
      DO JDFL=1,8
       DPP(:) = DCORVG(:,JDFL)
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

      DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3)) &
           -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3)) &
           +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
           
      OM=DOMEGA(ICUBP)*ABS(DETJ)

      CALL ELE(XI1,XI2,XI3,-3)

!         write(*,*) HBASJ1,HBASJ2,HBASJ3,HBASJ4
!     Summing up over all pairs of multiindices
      DO JDOFEH=1,IDFL
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
      
!          write(*,*) idfl
!          write(*,*) HBASJ1,HBASJ2,HBASJ3,HBASJ4
       DO IDOFEH=1,IDFL
        IF (IDOFEH.EQ.JDOFEH) THEN
         AH=dVisc*(HBASJ2*HBASJ2+HBASJ3*HBASJ3+HBASJ4*HBASJ4)
        ELSE
         HBASI2=DBAS(1,IDOFEH,2)
         HBASI3=DBAS(1,IDOFEH,3)
         HBASI4=DBAS(1,IDOFEH,4)
         AH=dVisc*(HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4)
        ENDIF
        DENTRY(JDOFEH,IDOFEH)=DENTRY(JDOFEH,IDOFEH) + OM*AH
       END DO
      END DO

      END DO
      
      dMat%nu = 27
      dMat%na = 27*27
      
      allocate(dMat%LdA(dMat%nu+1))
      allocate(dMat%ColA(dMat%na))
      allocate(RHS(dMat%nu))
      allocate(SOL(dMat%nu))
      RHS= 0d0 
      SOL(9) = 0.15d0
      SOL(10) = 0.15d0
      SOL(11) = 0.15d0
      SOL(12) = 0.1d0
      SOL(21) = 0.175d0
      
      dMat%LdA(1) = 1
      DO i=1,27
       dMat%LdA(i+1) = dMat%LdA(i) + 27
       jj = 0
       DO j=dMat%LdA(i),dMat%LdA(i+1)-1
        jj = jj + 1
        dMat%ColA(j) = jj
       end do
      end do
      
      !       write(*,*) DENTRY
!       pause
!       write(*,*) dMat%ColA
!       pause
!       write(*,*) dMat%LdA
      DO i=1,8
       DENTRY(i,1:27) = 0d0
       DENTRY(i,i) = 1d0
      end do
      DO i=9,12
       DENTRY(i,1:27) = 0d0
       DENTRY(i,i) = 1d0
      end do
      i=21
      DENTRY(i,1:27) = 0d0
      DENTRY(i,i) = 1d0
       
      DO i=1,27
       DO j=dMat%LdA(i),dMat%LdA(i+1)-1
        jj = dMat%ColA(j)
        RHS(i) = RHS(i) - DENTRY(i,jj)*SOL(jj)
       end do
       write(*,*) i,rhs(i)
      end do
      
      DO i=1,8 
       RHS(i) = 0d0
      END DO
      RHS(9) = 0d0
      RHS(10) = 0d0
      RHS(11) = 0d0
      RHS(12) = 0d0
      RHS(21) = 0d0
      SOL = 0d0

      CALL myUmfPack_Factorize(DENTRY,dMat)
      
      CALL myUmfPack_Solve(SOL,RHS,DENTRY,dMat,1)
      
      WRITE(*,*) 'SOL'
      DO i=1,27
       write(*,*) i,sol(i)
      end do
      
      CALL myUmfPack_Free()

END subroutine
