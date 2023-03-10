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
* XORA                                                                 *
*                                                                      *
* Purpose  Get arrays, previously stored by XOWA, back on DWORK        *
*                                                                      *
* Subroutines/functions called  ORA0, ZNEW, ZLEN, ZTYPE, ZDISP         *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* LNR     I*4    Number of array to be stored                          *
* ARR     C*6    Name of array                                         *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
*                CFILE.EQ.'SCRATCH' means temporary file               *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -106  Wrong value of  LNR                             *
*                -110  MFILE  exceeds maximum range 16,...,80          *
*                -110  Error while reading from unit  MFILE            *
*                -112  File  CFILE  could not be opened                *
*                -113  Wrong value of  ITYPE  read from  MFILE         *
*                                                                      *
*************************************************************************
C
      SUBROUTINE XORA(LNR,ARR,MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFILE*(*),CFORM*15
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XORA'
      IF (ICHECK.GE.997) CALL OTRC('XORA  ','12/11/89')
C
      IER=0
      BWARN=.FALSE.
      BFMT=IFMT.EQ.1
C
C *** Valid array number - 0 means that array does not yet exist
C
      IF ((LNR.LT.0).OR.(LNR.GT.NNARR)) THEN
       WRITE (CPARAM,'(A6)') ARR
       CALL WERR(-106,'XORA  ')
       GOTO 99999
      ENDIF
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,IFMT)
      IF (IER.NE.0) GOTO 99999
C
      IF (BFMT) THEN
       READ (MFILE,'(2A15,2I15)',ERR=99997,END=99997)
     *       ARR,CFORM,ITYPE,ILEN
      ELSE
       READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
      ENDIF
C
      IF ((ITYPE.LT.1).OR.(ITYPE.GT.3).OR.(ILEN.LT.1)) GOTO 99997
C
      IF (LNR.EQ.0) THEN
       CALL ZNEW(ILEN,-ITYPE,LNR,ARR)
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZTYPE(LNR,JTYPE)
       IF (ITYPE.NE.JTYPE) THEN
        WRITE (CPARAM,'(A6,I15)') ARR,MFILE
        CALL WERR(-113,'XORA  ')
        GOTO 99999
       ELSE
        CALL ZLEN(LNR,JLEN)
        IF (JLEN.LT.ILEN) THEN
         CALL ZDISP(0,LNR,ARR)
         CALL ZNEW(ILEN,-ITYPE,LNR,ARR)
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
      ENDIF
C
      L1=L(LNR)
      IF (BFMT) THEN
       GOTO (111,222,333),ITYPE
111    READ(MFILE,CFORM,ERR=99997,END=99997) (DWORK(L1+I),I=0,ILEN-1)
       GOTO 555
222    READ(MFILE,CFORM,ERR=99997,END=99997) (VWORK(L1+I),I=0,ILEN-1)
       GOTO 555
333    READ(MFILE,CFORM,ERR=99997,END=99997) (KWORK(L1+I),I=0,ILEN-1)
555    CONTINUE
      ELSE
       IF (ITYPE.GT.1) L1=(L1+1)/2
       CALL ORA0(DWORK(L1),ILEN8,MFILE,JRECL)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'XORA  ')
      GOTO 99999
C
99997 WRITE (CPARAM,'(A6,2I15)') MFILE
      CALL WERR(-110,'XORA  ')
C
99999 END
C
C
C
      SUBROUTINE ORA0(DX,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('ORA0  ','12/11/89')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
1     READ (MFILE,ERR=99998,END=99998) (DX(J1+I),I=1,IRECL8)
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      READ (MFILE,ERR=99998,END=99998) (DX(J1+I),I=1,MOD(N,IRECL8))
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORA0  ')
C
99999 END
