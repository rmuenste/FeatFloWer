************************************************************************
      SUBROUTINE BDRNEU(KABD,KNPR,DCORVG,INEUM,KAREA,KADJ,KVERT,KEBD,
     *                  KVBD,KFBD,AFBD,NFBD,ILEV)
C
************************************************************************
*    Purpose:  sets the DIRICHLET- and NEUMANN components
*-----------------------------------------------------------------------
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KABD(*),KEBD(*),KAREA(6,*),KADJ(6,*)
      DIMENSION KNPR(*),DCORVG(3,*),KVERT(8,*),KVBD(*)
      DIMENSION KFBD(*),AFBD(*)
	REAL*8 NFBD(3,*)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      DIMENSION KAUX(1)
      SAVE
C
C
      IF (ilev.gt.1.and.myid.ne.MASTER) THEN 
       DO IUU=1,mg_mpi(ILEV)%NeighNum                            ! PARALLEL
        mg_mpi(ILEV)%parST(IUU)%i=1                              ! PARALLEL
       END DO                                                    ! PARALLEL
      END IF

      NABD =0
      NVBD =0
      INEUM=0
C
      DO 1 IEL=1,NEL
C
      DO 2 IAT=1,6
      IAR =KAREA(IAT,IEL)
      INPR=KNPR(NVT+IAR)
      IF (INPR.EQ.0) GOTO 2
C
      NABD=NABD+1
      KABD(NABD)=IAR
      KEBD(NABD)=IEL
      KAUX(NABD)=IAT
C
2     CONTINUE
C
1     CONTINUE
C
      DO 10 IABD=1,NABD
C
      IEL=KEBD(IABD) 
      IAT=KAUX(IABD) 
      KFBD(IABD)=IAT
C
      PXC=0D0
      PYC=0D0
      PZC=0D0
C
      DO 20 IVE=1,8
      IVT=KVERT(IVE,IEL)
      PXC=PXC+DCORVG(1,IVT)
      PYC=PYC+DCORVG(2,IVT)
      PZC=PZC+DCORVG(3,IVT)
 20   CONTINUE
C
C *** PXC,PYC,PZC are coordinates of the center of the element
      PXC=0.125D0*PXC
      PYC=0.125D0*PYC
      PZC=0.125D0*PZC
C
      IF (IAT.EQ.1) THEN
       IVT1=KVERT(1,IEL)
       IVT2=KVERT(2,IEL)
       IVT3=KVERT(3,IEL)
       IVT4=KVERT(4,IEL)
       GOTO 30
      ENDIF
C
      IF (IAT.EQ.2) THEN
       IVT1=KVERT(1,IEL)
       IVT2=KVERT(2,IEL)
       IVT3=KVERT(6,IEL)
       IVT4=KVERT(5,IEL)
       GOTO 30
      ENDIF
C
      IF (IAT.EQ.3) THEN
       IVT1=KVERT(2,IEL)
       IVT2=KVERT(3,IEL)
       IVT3=KVERT(7,IEL)
       IVT4=KVERT(6,IEL)
       GOTO 30
      ENDIF
C
      IF (IAT.EQ.4) THEN
       IVT1=KVERT(3,IEL)
       IVT2=KVERT(4,IEL)
       IVT3=KVERT(8,IEL)
       IVT4=KVERT(7,IEL)
       GOTO 30
      ENDIF
C
      IF (IAT.EQ.5) THEN
       IVT1=KVERT(4,IEL)
       IVT2=KVERT(1,IEL)
       IVT3=KVERT(5,IEL)
       IVT4=KVERT(8,IEL)
       GOTO 30
      ENDIF
C
      IF (IAT.EQ.6) THEN
       IVT1=KVERT(5,IEL)
       IVT2=KVERT(6,IEL)
       IVT3=KVERT(7,IEL)
       IVT4=KVERT(8,IEL)
       GOTO 30
      ENDIF
C
C
30    P1X=DCORVG(1,IVT1)
      P1Y=DCORVG(2,IVT1)
      P1Z=DCORVG(3,IVT1)
      P2X=DCORVG(1,IVT2)
      P2Y=DCORVG(2,IVT2)
      P2Z=DCORVG(3,IVT2)
      P3X=DCORVG(1,IVT3)
      P3Y=DCORVG(2,IVT3)
      P3Z=DCORVG(3,IVT3)
      P4X=DCORVG(1,IVT4)
      P4Y=DCORVG(2,IVT4)
      P4Z=DCORVG(3,IVT4)
C
      PXA=(P1X+P2X+P3X+P4X)*0.25D0
      PYA=(P1Y+P2Y+P3Y+P4Y)*0.25D0
      PZA=(P1Z+P2Z+P3Z+P4Z)*0.25D0
C
      AX2=P2X-P1X
      AY2=P2Y-P1Y
      AZ2=P2Z-P1Z
      AY3=P3Y-P1Y
      AX3=P3X-P1X
      AZ3=P3Z-P1Z
      AY4=P4Y-P1Y
      AX4=P4X-P1X
      AZ4=P4Z-P1Z
C
      AX =PXC-PXA
      AY =PYC-PYA
      AZ =PZC-PZA
C
      DNX=(AY3*AZ2)-(AZ3*AY2)
      DNY=(AZ3*AX2)-(AX3*AZ2)
      DNZ=(AX3*AY2)-(AY3*AX2)
      DNAR1=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
      DNX=(AY4*AZ3)-(AZ4*AY3)
      DNY=(AZ4*AX3)-(AX4*AZ3)
      DNZ=(AX4*AY3)-(AY4*AX3)
      DNAR2=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
      DNA=0.5D0*(DNAR1+DNAR2)
      DNX=DNX/DNAR2
      DNY=DNY/DNAR2
      DNZ=DNZ/DNAR2
C
      DHN=DNX*AX+DNY*AY+DNZ*AZ
      IF (DHN.GT.0D0) THEN
       DNX=-DNX
       DNY=-DNY
       DNZ=-DNZ
      ENDIF
C
      AFBD(  IABD)=DNA
      NFBD(1,IABD)=DNX
      NFBD(2,IABD)=DNY
      NFBD(3,IABD)=DNZ
C
      CALL NEUDAT(PXA,PYA,PZA,TIMENS,INPR,1)

      IF (ilev.gt.1.and.myid.ne.MASTER) THEN 
       DO pID=1,mg_mpi(ILEV)%NeighNum                                  ! PARALLEL
        IUU = mg_mpi(ILEV)%parST(pID)%i                                ! PARALLEL
        IF (mg_mpi(ILEV)%parST(pID)%FaceLink(1,IUU).EQ.KABD(IABD))THEN ! PARALLEL
         mg_mpi(ILEV)%parST(pID)%i=IUU+1                               ! PARALLEL
         INPR=0                                                        ! PARALLEL
         GOTO 9                                                        ! PARALLEL
        END IF                                                         ! PARALLEL
       END DO                                                          ! PARALLEL
      END IF
C
9     KNPR(NVT+KABD(IABD))=INPR
      IF (INPR.EQ.0) INEUM=1
C
 10   CONTINUE
C
      DO 40 IVT=1,NVT
        IF (KNPR(IVT).EQ.0) GOTO 40
        NVBD=NVBD+1 
        KVBD(NVBD)=IVT
        PX=DCORVG(1,IVT)
        PY=DCORVG(2,IVT)
        PZ=DCORVG(3,IVT)
        CALL NEUDAT(PX,PY,PZ,TIMENS,INPR,0)
        KNPR(IVT)=INPR
 40   CONTINUE  
C
C
      IF (ilev.gt.1.and.myid.ne.MASTER) THEN 
      DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
       IF (mg_mpi(ILEV)%parST(pID)%i-1.NE.                          ! PARALLEL
     *     mg_mpi(ILEV)%parST(pID)%Num) WRITE(*,*)                  ! PARALLEL
     *     "Problem with",myid,"--",mg_mpi(ILEV)%parST(pID)%Neigh   ! PARALLEL
      END DO
      END IF

      END
C
************************************************************************
      SUBROUTINE BDRSET(DU1,DU2,DU3,KABD,NABD,DCORVG,
     *                  KVERT,KAREA,KELBD,UE,KNPR,NFBD,NVT)
************************************************************************
*    Purpose:  updates the solution vector (DU1,DU2,DU3) and the right 
*              hand side (DF1,DF2,DF3) for all DIRICHLET boundary nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DU1(*),DU2(*),DU3(*),DN(3)
      DIMENSION KABD(*),KELBD(*),DCORVG(3,*)
	REAL*8 NFBD(3,*)
      DIMENSION KVERT(8,*),KAREA(6,*),KNPR(*)
      EXTERNAL UE
C
C
      DO 1 IAT=1,NABD
C
      IABD=KABD(IAT)
      INPR=KNPR(NVT+IABD)
C
      IF (INPR.LT.0) THEN
       DN(1)=NFBD(1,IAT)
       DN(2)=NFBD(2,IAT)
       DN(3)=NFBD(3,IAT)
       UN=DU1(IABD)*DN(1)+DU2(IABD)*DN(2)+DU3(IABD)*DN(3)
       DU1(IABD)=DU1(IABD)-UN*DN(1)
       DU2(IABD)=DU2(IABD)-UN*DN(2)
       DU3(IABD)=DU3(IABD)-UN*DN(3)
      ENDIF
      IF (INPR.LE.0) GOTO 1
C
      IELBD=KELBD(IAT)
      IVT1=KVERT(1,IELBD)
      IVT2=KVERT(2,IELBD)
      IVT3=KVERT(3,IELBD)
      IVT4=KVERT(4,IELBD)
      IVT5=KVERT(5,IELBD)
      IVT6=KVERT(6,IELBD)
      IVT7=KVERT(7,IELBD)
      IVT8=KVERT(8,IELBD)
C
      XV1=DCORVG(1,IVT1)
      YV1=DCORVG(2,IVT1)
      ZV1=DCORVG(3,IVT1)
      XV2=DCORVG(1,IVT2)
      YV2=DCORVG(2,IVT2)
      ZV2=DCORVG(3,IVT2)
      XV3=DCORVG(1,IVT3)
      YV3=DCORVG(2,IVT3)
      ZV3=DCORVG(3,IVT3)
      XV4=DCORVG(1,IVT4)
      YV4=DCORVG(2,IVT4)
      ZV4=DCORVG(3,IVT4)
      XV5=DCORVG(1,IVT5)
      YV5=DCORVG(2,IVT5)
      ZV5=DCORVG(3,IVT5)
      XV6=DCORVG(1,IVT6)
      YV6=DCORVG(2,IVT6)
      ZV6=DCORVG(3,IVT6)
      XV7=DCORVG(1,IVT7)
      YV7=DCORVG(2,IVT7)
      ZV7=DCORVG(3,IVT7)
      XV8=DCORVG(1,IVT8)
      YV8=DCORVG(2,IVT8)
      ZV8=DCORVG(3,IVT8)
C
      IF (IABD.EQ.KAREA(1,IELBD)) THEN
       PX=(XV1+XV2+XV3+XV4)*0.25D0
       PY=(YV1+YV2+YV3+YV4)*0.25D0
       PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(2,IELBD)) THEN
       PX=(XV1+XV2+XV6+XV5)*0.25D0
       PY=(YV1+YV2+YV6+YV5)*0.25D0
       PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(3,IELBD)) THEN
       PX=(XV2+XV3+XV6+XV7)*0.25D0
       PY=(YV2+YV3+YV6+YV7)*0.25D0
       PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(4,IELBD)) THEN
       PX=(XV3+XV4+XV8+XV7)*0.25D0
       PY=(YV3+YV4+YV8+YV7)*0.25D0
       PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(5,IELBD)) THEN
       PX=(XV4+XV1+XV5+XV8)*0.25D0
       PY=(YV4+YV1+YV5+YV8)*0.25D0
       PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(6,IELBD)) THEN
       PX=(XV5+XV6+XV7+XV8)*0.25D0
       PY=(YV5+YV6+YV7+YV8)*0.25D0
       PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
      ENDIF
C
      DU1(IABD)=UE(PX,PY,PZ,1)
      DU2(IABD)=UE(PX,PY,PZ,2)
      DU3(IABD)=UE(PX,PY,PZ,3)
C
1     CONTINUE
C
      END

C
C
************************************************************************
      SUBROUTINE BDRYA(VA,KCOL,KLD,KABD,NABD,KNPR,NVT)
************************************************************************
*    Purpose:  updates the matrix entries for all DIRICHLET boundary
*              nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 VA
      DIMENSION VA(*),KCOL(*),KLD(*),KABD(*),KNPR(*)
C
      DO 1 IAT=1,NABD
      IABD=KABD(IAT)
      INPR=KNPR(NVT+IABD)
      IF (INPR.EQ.0) GOTO 1
C
      DO 2 ICOL=KLD(IABD)+1,KLD(IABD+1)-1
2     VA(ICOL)=0E0
C
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE BDRYB(B1,B2,B3,KCOLB,KLDB,DNORM,KABD,NABD,KNPR,NVT)
************************************************************************
*     Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*8 B1,B2,B3
C
      DIMENSION B1(*),B2(*),B3(*),KCOLB(*),KLDB(*)
      DIMENSION KABD(*),KNPR(*),DNORM(3,*)


      DO 1 IAT=1,NABD
C
       IAREA=KABD(IAT)
       INPR=KNPR(NVT+IAREA)
       IF (INPR.GE.0) GOTO 1
C
       DNX=DNORM(1,IAT)
       DNY=DNORM(2,IAT)
       DNZ=DNORM(3,IAT)
C
       ILDB=KLDB(IAREA)
       DAUX=B1(ILDB)*DNX+B2(ILDB)*DNY+B3(ILDB)*DNZ
       B1(ILDB)=B1(ILDB)-DAUX*DNX
       B2(ILDB)=B2(ILDB)-DAUX*DNY
       B3(ILDB)=B3(ILDB)-DAUX*DNZ
C
 1     CONTINUE
C
      END
C
C
************************************************************************
      SUBROUTINE BDRY0(DD,KABD,NABD,KNPR,NVT)
************************************************************************
*    Purpose:  sets the DIRICHLET-components of the vector (D1,D2,D3) to
*              zero
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DD(*),KABD(*),KNPR(*)
C
      DO 1 IAT=1,NABD
      IABD=KABD(IAT)
      INPR=KNPR(NVT+IABD)
      IF (INPR.NE.0) DD(IABD)=0D0
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE BDRYS(U1,U2,U3,D1,D2,D3,VA,KLD,KABD,NABD,KNPR,NFBD,NVT,
     *                 THSTEP,IDEF)
************************************************************************
*     Purpose: Assembly of the surface stress integral
*-----------------------------------------------------------------------
C  
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*4 VA
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*),DN(3)
      DIMENSION VA(*),KLD(*),KABD(*),KNPR(*)
	REAL*8 NFBD(3,*)
C
      IF (ABS(THSTEP).LE.1D-8) RETURN
C
      DO 1 IAT=1,NABD
C
       IABD=KABD(IAT)
       INPR=KNPR(NVT+IABD)
       IF (INPR.EQ.0) GOTO 1
C
C      Impose Dirichlet boundary conditions
       IF (INPR.EQ.1.AND.IDEF.GE.1) THEN
        DO 2 ICOL=KLD(IABD)+1,KLD(IABD+1)-1
 2         VA(ICOL)=0E0
        D1(IABD)=0D0
        D2(IABD)=0D0
        D3(IABD)=0D0
        GOTO 1
       ENDIF
C
C      Add tangential stress given by symmetry BC
C      Project the defect onto the tangent plane
       IF (IDEF.GE.1) THEN
        DO 3 ICOL=KLD(IABD)+1,KLD(IABD+1)-1
 3           VA(ICOL)=0E0
        DN(1)=NFBD(1,IAT)
        DN(2)=NFBD(2,IAT)
        DN(3)=NFBD(3,IAT)
        DAUX=D1(IABD)*DN(1)+D2(IABD)*DN(2)+D3(IABD)*DN(3)
        D1(IABD)=D1(IABD)-DAUX*DN(1)
        D2(IABD)=D2(IABD)-DAUX*DN(2)
        D3(IABD)=D3(IABD)-DAUX*DN(3)
       ENDIF
C
 1    CONTINUE
C
      END
