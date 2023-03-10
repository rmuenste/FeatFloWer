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
* E032                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  32        *
*          nonconforming element with 5 local degrees of freedom       *
*          (Han element)                                               *
*                                                                      *
*          A*(5X**4 - 3X**2) + B*(5Y**4-3Y**2) +                       *
*           + C * X + D * Y + E                                        *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8    Evaluation point in cartesian coordinates            *
* XI2      R*8    with respect to the unit square                      *
* IPAR     I*4    Switch for desired action                            *
*                  0  Calculation of the desired values of DBAS        *
*                 -1  Set number of element                            *
*                 -2  Calculate values on the reference element        *
*                     for all cubature points and save them            *
*                 -3  same as 0, but use the values saved before       *
*                                                                      *
* BDER     LOG    Derivative J is calculated if BDER(J).EQ..TRUE.      *
*                 Multiindices are enumbered from 1 to 6               *
* DJAC     R*8    Jacobian                                             *
* DETJ     R*8    Determinant of the jacobian                          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DBAS     R*8    Array of function values and derivatives for         *
*                 all basis functions (block /ELEM/)                   *
* IPAR     I*4    Set to  32  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E032(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=5)
      DIMENSION DHELP(NDOFL,2)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      SUB='E032'
      IF (ICHECK.GE.998) CALL OTRC('E032  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
1     IPAR=32
2     GOTO 99999
C
3     CONTINUE
10    IF (ICHECK.EQ.0) GOTO 3200
      IF (DETJ.LT.0D0) CALL WERR(-130,'E032  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E032  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E032  ')
      IF (IER.NE.0) GOTO 99999
C
3200  IF (.NOT.BDER(1)) GOTO 101
C *** Function values
      DBAS(1,1)=-0.5D0*XI2+0.25D0*(5D0*(XI2**4)-3D0*(XI2**2))
      DBAS(2,1)= 0.5D0*XI1+0.25D0*(5D0*(XI1**4)-3D0*(XI1**2))
      DBAS(3,1)= 0.5D0*XI2+0.25D0*(5D0*(XI2**4)-3D0*(XI2**2))
      DBAS(4,1)=-0.5D0*XI1+0.25D0*(5D0*(XI1**4)-3D0*(XI1**2))
      DBAS(5,1)= 1D0-0.5D0*(5D0*(XI1**4)-3D0*(XI1**2))-
     *           0.5D0*(5D0*(XI2**4)-3D0*(XI2**2))

101   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C *** First order derivatives
      XJ=1D0/DETJ
      DHELP(1,1)= 0D0
      DHELP(2,1)= 0.5D0+5D0*(XI1**3)-1.5D0*XI1
      DHELP(3,1)= 0D0
      DHELP(4,1)=-0.5D0+5D0*(XI1**3)-1.5D0*XI1
      DHELP(5,1)=(-10D0*(XI1**3)+3D0*XI1)*1D0
      DHELP(1,2)=-0.5D0+5D0*(XI2**3)-1.5D0*XI2
      DHELP(2,2)= 0D0
      DHELP(3,2)= 0.5D0+5D0*(XI2**3)-1.5D0*XI2
      DHELP(4,2)= 0D0
      DHELP(5,2)=(-10D0*(XI2**3)+3D0*XI2)*1D0
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,2)= XJ*(DJAC(2,2)*DHELP(1,1)-DJAC(2,1)*DHELP(1,2))
      DBAS(2,2)= XJ*(DJAC(2,2)*DHELP(2,1)-DJAC(2,1)*DHELP(2,2))
      DBAS(3,2)= XJ*(DJAC(2,2)*DHELP(3,1)-DJAC(2,1)*DHELP(3,2))
      DBAS(4,2)= XJ*(DJAC(2,2)*DHELP(4,1)-DJAC(2,1)*DHELP(4,2))
      DBAS(5,2)= XJ*(DJAC(2,2)*DHELP(5,1)-DJAC(2,1)*DHELP(5,2))
102   IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)=-XJ*(DJAC(1,2)*DHELP(1,1)-DJAC(1,1)*DHELP(1,2))
      DBAS(2,3)=-XJ*(DJAC(1,2)*DHELP(2,1)-DJAC(1,1)*DHELP(2,2))
      DBAS(3,3)=-XJ*(DJAC(1,2)*DHELP(3,1)-DJAC(1,1)*DHELP(3,2))
      DBAS(4,3)=-XJ*(DJAC(1,2)*DHELP(4,1)-DJAC(1,1)*DHELP(4,2))
      DBAS(5,3)=-XJ*(DJAC(1,2)*DHELP(5,1)-DJAC(1,1)*DHELP(5,2))
C
99999 END
