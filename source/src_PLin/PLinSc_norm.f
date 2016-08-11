************************************************************************
      SUBROUTINE Build_NormQP1(DNORM,KVERT,DCORVG,NEL)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION B1,B2,B3
C
      DIMENSION KVERT(8,*),DCORVG(3,*)
      DIMENSION DNORM(6,3,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      DO 10 IEL=1,NEL
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
       IVT3=KVERT(6,IEL)
       IVT4=KVERT(7,IEL)
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
      DNORM(IAT,1,IEL ) =DNX
      DNORM(IAT,2,IEL ) =DNY
      DNORM(IAT,3,IEL ) =DNZ
C
20    CONTINUE
C
10    CONTINUE
C
C
C
      END
