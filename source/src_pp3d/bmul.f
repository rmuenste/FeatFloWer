************************************************************************
      SUBROUTINE BBUILD(KVERT,KAREA,KADJ,DCORVG,B1,B2,B3,KCOLB,KLDB,
     *                  KABD,KERIOA,KERIOB,KCROSS,NB)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL*8 B1,B2,B3
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KADJ(NNAE,*),DCORVG(3,*)
      DIMENSION B1(*),B2(*),B3(*),KCOLB(*),KLDB(*),KCROSS(*)
      DIMENSION KABD(*),KERIOA(*),KERIOB(NEL,2)
      SAVE
C
C
C
      ILD=0
C
      DO 10 IEL=1,NEL
C
      IELPERFOUND = 0
C
      PXC=0D0
      PYC=0D0
      PZC=0D0
C
      DO 12 IVE=1,8
      IVT=KVERT(IVE,IEL)
      PXC=PXC+DCORVG(1,IVT)
      PYC=PYC+DCORVG(2,IVT)
      PZC=PZC+DCORVG(3,IVT)
12    CONTINUE
C
C *** PXC,PYC,PZC are coordinates of the center of the element
      PXC=0.125D0*PXC
      PYC=0.125D0*PYC
      PZC=0.125D0*PZC
C
C
      DO 20 IAT=1,6
      IAREA = KAREA(IAT,IEL)
      IF (KCROSS(IAREA).EQ.0) THEN
       IADJ=KADJ(IAT,IEL)
       IF ((IADJ.GT.0).AND.(IADJ.LT.IEL)) GOTO 20
      ELSE
       IELPERFOUND = IELPERFOUND + 1
       IF (IELPERFOUND.EQ.1) THEN
        IADJ = ABS(KERIOB(IEL,1))
       ELSE
        IADJ = ABS(KERIOB(IEL,2))
       END IF
      END IF
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
      DFAC=0.5D0*(DNAR1+DNAR2)/DNAR2
      DNX =DFAC*DNX
      DNY =DFAC*DNY
      DNZ =DFAC*DNZ
C
      DHN=DNX*AX+DNY*AY+DNZ*AZ
      IF (DHN.LT.0D0) THEN
       DNX=-DNX
       DNY=-DNY
       DNZ=-DNZ
      ENDIF
C
C
      IAREA=KAREA(IAT,IEL)
      ILD  =ILD+1
      KLDB (IAREA)=ILD
      KCOLB(ILD ) =IEL
      B1   (ILD ) =DNX
      B2   (ILD ) =DNY
      B3   (ILD ) =DNZ
C
      IF (IADJ.GT.0) THEN
       ILD=ILD+1
       KCOLB(ILD)= IADJ
       B1   (ILD)=-DNX
       B2   (ILD)=-DNY
       B3   (ILD)=-DNZ
      ENDIF
C
C
20    CONTINUE
C
10    CONTINUE
C
      NB         =ILD
      KLDB(NAT+1)=ILD+1
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE BMUL1 (KVERT,KAREA,KADJ,DCORVG,DP,DF1,DF2,DF3,KERIOA,
     *                  KERIOB,KCROSS,KABD,A1,A2,ILEV)
************************************************************************
      USE PP3D_MPI
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KADJ(NNAE,*),DCORVG(3,*)
      DIMENSION DP(*),DF1(*),DF2(*),DF3(*),KCROSS(*),KABD(*)
      DIMENSION KERIOA(*),KERIOB(NEL,2)
      SAVE
C
C
      IF (ABS(A1).LE.1D-8) RETURN
C
      DO IUU=1,mg_mpi(ILEV)%NeighNum                            ! PARALLEL
       mg_mpi(ILEV)%parST(IUU)%i=1                              ! PARALLEL
      END DO                                                    ! PARALLEL
C
      CALL CommBMul(DP,ILEV)                                    ! PARALLEL
C
      DO 10 IEL=1,NEL
C
      IELPERFOUND = 0
C
      PXC=0D0
      PYC=0D0
      PZC=0D0
C
      DO 12 IVE=1,8
      IVT=KVERT(IVE,IEL)
      PXC=PXC+DCORVG(1,IVT)
      PYC=PYC+DCORVG(2,IVT)
      PZC=PZC+DCORVG(3,IVT)
12    CONTINUE
C
C *** PXC,PYC,PZC are coordinates of the center of the element
      PXC=0.125D0*PXC
      PYC=0.125D0*PYC
      PZC=0.125D0*PZC
C
C
      DO 20 IAT=1,6
      IADJ=KADJ(IAT,IEL)
      IF (KCROSS(KAREA(IAT,IEL)).EQ.0) THEN
       IF ((IADJ.GT.0).AND.(IADJ.LT.IEL)) GOTO 20
      END IF
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
      DFAC=0.5D0*(DNAR1+DNAR2)/DNAR2
      DNX =DFAC*DNX
      DNY =DFAC*DNY
      DNZ =DFAC*DNZ
C
      DHN=DNX*AX+DNY*AY+DNZ*AZ
      IF (DHN.LT.0D0) THEN
       DNX=-DNX
       DNY=-DNY
       DNZ=-DNZ
      ENDIF
C
      IAREA=KAREA(IAT,IEL)
      IF (KCROSS(IAREA).NE.0) THEN
       IELPERFOUND = IELPERFOUND + 1
       IF (IELPERFOUND.EQ.1) THEN
        IADJ = ABS(KERIOB(IEL,1))
       ELSE
        IADJ = ABS(KERIOB(IEL,2))
       END IF
      END IF
      IF (IADJ.GT.0) THEN
       IF (KCROSS(IAREA).EQ.0) THEN
        IF (A2.EQ.0D0) THEN
         DPDIFF=DP(IEL)-DP(IADJ)     !!!! CORRECTION
        ELSE
         DPDIFF=DP(IEL)-DP(IADJ)     !!!! GRAD
        END IF
       ELSE
        IF (A2.EQ.0D0) THEN
         DPDIFF=DP(IEL)-DP(IADJ)    !!!! CORRECTION                  ! PERIODIC
        ELSE                                                         ! PERIODIC
         IF (KCROSS(IAREA).LT.0) THEN                                ! PERIODIC
          DPDIFF=0.0D0               !!!! GRAD                       ! PERIODIC
         ELSE                                                        ! PERIODIC
          DPDIFF=(DP(IEL)-DP(IADJ))  !!!! GRAD                       ! PERIODIC
         END IF
        END IF
       END IF
      ELSE
       DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
        IUU = mg_mpi(ILEV)%parST(pID)%i                              ! PARALLEL
        IF (mg_mpi(ILEV)%parST(pID)%FaceLink(1,IUU).EQ.IAREA) THEN   ! PARALLEL
         DP_ADJ = mg_mpi(ILEV)%parST(pID)%RDVect(IUU)                ! PARALLEL
         DPDIFF=DP(IEL)-(1D0-A2)*DP_ADJ                              ! PARALLEL
         mg_mpi(ILEV)%parST(pID)%i=IUU+1                             ! PARALLEL
!         WRITE(*,*)myid," parallel", mg_mpi(ILEV)%parST(pID)%Neigh
         GOTO 66                                                     ! PARALLEL
        END IF                                                       ! PARALLEL
       END DO                                                        ! PARALLEL
       CALL NEUDAT(PXA,PYA,PZA,TIMENS,INPR,0)
       IF (INPR.LT.0.OR.KCROSS(IAREA).NE.0) THEN
         DPDIFF = 0.0D0
       ELSE
         DPDIFF=DP(IEL)
       END IF
66     CONTINUE                                                      ! PARALLEL
      END IF
      DF1(IAREA)=A2*DF1(IAREA)+A1*DNX*DPDIFF
      DF2(IAREA)=A2*DF2(IAREA)+A1*DNY*DPDIFF
      DF3(IAREA)=A2*DF3(IAREA)+A1*DNZ*DPDIFF
C
20    CONTINUE
C
10    CONTINUE
C
C
C
      DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
       IF (mg_mpi(ILEV)%parST(pID)%i-1.NE.                          ! PARALLEL
     *     mg_mpi(ILEV)%parST(pID)%Num) WRITE(*,*)                  ! PARALLEL
     *     "Problem with",myid,"--",mg_mpi(ILEV)%parST(pID)%Neigh   ! PARALLEL
      END DO

      END
C
C
C
************************************************************************
      SUBROUTINE BTMUL1(KVERT,KAREA,KADJ,DCORVG,DU1,DU2,DU3,DFP,
     *                  NEL,NVT,NAT,A1)
************************************************************************
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KADJ(NNAE,*),DCORVG(3,*)
      DIMENSION DFP(*),DU1(*),DU2(*),DU3(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      DO 10 IEL=1,NEL
      DFP(IEL)=0D0
C
      PXC=0D0
      PYC=0D0
      PZC=0D0
C
      DO 12 IVE=1,8
      IVT=KVERT(IVE,IEL)
      PXC=PXC+DCORVG(1,IVT)
      PYC=PYC+DCORVG(2,IVT)
      PZC=PZC+DCORVG(3,IVT)
12    CONTINUE
C
C *** PXC,PYC,PZC are coordinates of the center of the element
      PXC=0.125D0*PXC
      PYC=0.125D0*PYC
      PZC=0.125D0*PZC
C
C
      DO 20 IAT=1,6
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
      DFAC=0.5D0*(DNAR1+DNAR2)/DNAR2
      DNX =DFAC*DNX
      DNY =DFAC*DNY
      DNZ =DFAC*DNZ
C
      DHN=DNX*AX+DNY*AY+DNZ*AZ
      IF (DHN.LT.0D0) THEN
       DNX=-DNX
       DNY=-DNY
       DNZ=-DNZ
      ENDIF
C
C
      IAREA=KAREA(IAT,IEL)
      DFP(IEL)=DFP(IEL)+A1*(DNX*DU1(IAREA)+DNY*DU2(IAREA)+
     *                      DNZ*DU3(IAREA))
C
20    CONTINUE
C
!      if (myid.eq.1) WRITE(*,*) DFP(IEL)
10    CONTINUE

!       if (myid.eq.2) then
!       ILEV=2
!       DO pID=1,mg_mpi(ILEV)%NeighNum
!       DO I=1,mg_mpi(ILEV)%parST(pID)%Num
!        WRITE(*,*) DFP(mg_mpi(ILEV)%parST(pID)%ElemLink(1,I))
!       END DO
!       END DO
!       end if
! 
!       STOP
C
C
C
      END
C
C
