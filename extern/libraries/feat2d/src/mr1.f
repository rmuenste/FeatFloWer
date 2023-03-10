************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* MR1nn                                                                *
*                                                                      *
* Purpose  Multigrid restriction (version 1)                           *
*          nn refers to elementtype                                    *
*          Double precision version                                    *
*          The meshes are assumed to be generated by SA1 or SB1,       *
*          same numbering conventions for vertices and elements        *
*          Uniform refinement assumed                                  *
*                                                                      *
* Subroutines/functions called   LCP1                                  *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DU2      R*8    Vector on fine grid                                  *
* KVERT1   I*4    1 fine grid information                              *
* KVERT2   I*4    2 fine grid information                              *
* NEL1     I*4                                                         *
* KMAVT    I*4                                                         *
* KMAADJ   I*4    Macro-element information                            *
* KMAVE    I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DU1      R*8    Vector on coarse grid                                *
*                                                                      *
************************************************************************
C
      SUBROUTINE MR111(DU2,DU1,NEL1,KVERT2,KVERT1,KMAVT,KMAADJ,KMAVE)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B) 
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      PARAMETER (Q2=.5D0,Q4=.25D0)
      DIMENSION DU1(*),DU2(*)
      DIMENSION KVERT1(NNVE,*),KVERT2(NNVE,*)
      DIMENSION KMAVT(NNVE,*),KMAADJ(NNVE,*),KMAVE(2,*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MACROD/ NMAEL,NMAVT,NMAEDG,NMAVE,NMAVEL,NMABCT,NMAVBD
      SAVE /ERRCTL/,/CHAR/,/MACROD/
C
      SUB='MR111'
      IF (ICHECK.GE.998) CALL OTRC('MR111 ','04/12/91')
C

C *** NFINE1 is the degree of refinement on the coarse grid
      NFINE1=NINT(SQRT(REAL(NEL1/NMAEL)))
      NFINE2=2*NFINE1
      NF11=NFINE1-1
      NF21=NFINE2-1
C
C *** Step 1 - copy values for macro-vertices
      CALL LCP1(DU2,DU1,NMAVT)
C
      DO 11 IMAEL=1,NMAEL
      IEL1=(IMAEL-1)*NFINE1**2
      IEL2=(IMAEL-1)*NFINE2**2
C
      DO 12 IVE=1,4
      IF (IVE.EQ.1) THEN
       IEL11=IEL1+1
       IEL21=IEL2+1
      ELSE IF (IVE.EQ.2) THEN
       IEL11=IEL1+NFINE1
       IEL21=IEL2+NFINE2
      ELSE IF (IVE.EQ.3) THEN
       IEL11=IEL1+NFINE1**2
       IEL21=IEL2+NFINE2**2
      ELSE
       IEL11=IEL1+NFINE1**2-NFINE1+1
       IEL21=IEL2+NFINE2**2-NFINE2+1
      ENDIF
      IVT1=KVERT1(IVE,IEL11)
      IVT21=KVERT2(MOD(IVE,4)+1,IEL21)
      IVT22=KVERT2(MOD(IVE+1,4)+1,IEL21)
      IVT23=KVERT2(MOD(IVE+2,4)+1,IEL21)
      DU1(IVT1)=DU1(IVT1)+
     *          Q4*(DU2(IVT21)+DU2(IVT22)+DU2(IVT23))
      IF (KMAADJ(IVE,IMAEL).EQ.0) DU1(IVT1)=DU1(IVT1)+Q4*DU2(IVT21)
      IF (KMAADJ(MOD(IVE+2,4)+1,IMAEL).EQ.0)
     *          DU1(IVT1)=DU1(IVT1)+Q4*DU2(IVT23)
C
12    CONTINUE
11    CONTINUE
C
C
      DO 2 IMAEDG=1,NMAEDG
      IVT1=NMAVT+(IMAEDG-1)*NF11
      IVT2=NMAVT+(IMAEDG-1)*NF21
      DO 2 I=1,NF11
      IVT11=IVT1+I
      IVT21=IVT2+2*I
      DU1(IVT11)=DU2(IVT21)+Q2*(DU2(IVT21-1)+DU2(IVT21+1))
2     CONTINUE 
C
      DO 21 IMAEL=1,NMAEL
      IEL1=(IMAEL-1)*NFINE1**2
      IEL2=(IMAEL-1)*NFINE2**2
      DO 21 I=1,NF11
C
      DO 22 IVE=1,4
      IF (IVE.EQ.1) THEN
       IEL11=IEL1+I
       IEL21=IEL2+2*I
       JEL=1
      ELSE IF (IVE.EQ.2) THEN
       IEL11=IEL1+I*NFINE1
       IEL21=IEL2+2*I*NFINE2
       JEL=NFINE2
      ELSE IF (IVE.EQ.3) THEN
       IEL11=IEL1+NFINE1**2-I+1
       IEL21=IEL2+NFINE2**2-2*I+1
       JEL=-1
      ELSE 
       IEL11=IEL1+NFINE1**2-I*NFINE1+1
       IEL21=IEL2+NFINE2**2-2*I*NFINE2+1
       JEL=-NFINE2
      ENDIF
      IVT1=KVERT1(MOD(IVE,4)+1,IEL11)
      IVT21=KVERT2(MOD(IVE+2,4)+1,IEL21)
      IVT22=KVERT2(MOD(IVE+1,4)+1,IEL21)
      IVT23=KVERT2(MOD(IVE+1,4)+1,IEL21+JEL)
      DU1(IVT1)=DU1(IVT1)+
     *          Q2*DU2(IVT22)+Q4*(DU2(IVT21)+DU2(IVT23))
C
22    CONTINUE
21    CONTINUE
C
C
      IVT1=NMAVT+NMAEDG*NF11
      IVT2=NMAVT+NMAEDG*NF21
C
      DO 3 IMAEL=1,NMAEL
C *** Step 3 - define values interior to macro-elements
      IVT11=IVT1+(IMAEL-1)*NF11**2
      IVT21=IVT2+(IMAEL-1)*NF21**2
      DO 31 I=1,NF11
      IVT11A=IVT11+(I-1)*NF11
      IVT21A=IVT21+(2*I-1)*NF21
      DO 31 J=1,NF11
      IVT11B=IVT11A+J
      IVT21B=IVT21A+2*J
      DU1(IVT11B)=DU2(IVT21B)+
     *            Q2*(DU2(IVT21B+1)+DU2(IVT21B-1)+
     *                DU2(IVT21B+NF21)+DU2(IVT21B-NF21))+
     *            Q4*(DU2(IVT21B+NF21+1)+DU2(IVT21B+NF21-1)+
     *                DU2(IVT21B-NF21+1)+DU2(IVT21B-NF21-1))
31     CONTINUE
C
3     CONTINUE
C
      END
