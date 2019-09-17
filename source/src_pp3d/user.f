C     Processing of Ctrl C signals
      SUBROUTINE USER_OPT(tout)
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2

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

      DOUBLE PRECISION VRPARM,NY,NY0
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
      COMMON /LESPAR/ D_SMAG,I_SMAG

      CHARACTER prompt*3,cfile*60

      WRITE(*,10)
      PRINT '(" PP3D_LES suspended due to user interuption... "/
     *" Interactive session - possible options: "/
     *" * Write out the solution vector in raw format           - w"/
     *" * Read in the solution vector in raw format             - r"/
     *" * Write out the solution vector in GMV format           - wg"/
     *" * Change the tolerance for the defect (EPSUD)           - ed"/
     *" * Change the tolerance for the relative changes (EPSUR) - er"/
     *" * Change the viscosity coefficient                      - d"/
     *" * Change the value of the actual time                   - t"/
     *" * Change the time step value                            - td"/
     *" * Change maximum time value                             - tm"/
     *" * Change the GMV output frequency                       - tg"/
     *" * Quit/Stop the program                                 - q"/
     *" * Continue the program/quit this interactive section    - c")'

1     CONTINUE 

      PRINT '("Waiting for input ...")'
      READ(*,'(A)') prompt

      IF     ((prompt.EQ."TD ").or.(prompt.EQ."td ")) THEN
          WRITE(*,'(A10,G12.5)') " tstep = ", tstep
          READ (*,*,ERR=2) tstep
      ELSEIF ((prompt.EQ."T  ").or.(prompt.EQ."t  ")) THEN
          WRITE(*,'(A10,G12.5)') " timens = ", timens
          READ (*,*,ERR=2) timens
      ELSEIF ((prompt.EQ."TM ").or.(prompt.EQ."tm ")) THEN
          WRITE(*,'(A10,G12.5)') " timemx = ", timemx
          READ (*,*,ERR=2) timemx
      ELSEIF ((prompt.EQ."ED ").or.(prompt.EQ."ed ")) THEN
          WRITE(*,'(A10,G12.5)') " epsud = ", epsud
          READ (*,*,ERR=2) epsud
      ELSEIF ((prompt.EQ."ER ").or.(prompt.EQ."er ")) THEN
          WRITE(*,'(A10,G12.5)') " epsur = ", epsur
          READ (*,*,ERR=2) epsur
      ELSEIF ((prompt.EQ."DS ").or.(prompt.EQ."ds ")) THEN
          WRITE(*,'(A10,G12.5)') " dsmag = ", d_smag
          READ (*,*,ERR=2) d_smag
      ELSEIF ((prompt.EQ."D  ").or.(prompt.EQ."d  ")) THEN
          WRITE(*,'(A10,G12.5)') " re = ", re
          READ (*,*,ERR=2) re
          ny0=ny
          ny = 1D0/re
          CALL ReCompute_RE(ny0)
      ELSEIF ((prompt.EQ."TG ").or.(prompt.EQ."tg ")) THEN
          WRITE(*,'(A10,G12.5)') " dtgmv = ", dtgmv
          READ (*,*,ERR=2) dtgmv
      ELSEIF ((prompt.EQ."W  ").or.(prompt.EQ."w  ")) THEN
          WRITE(*,*)
     *    "Enter the file name for writing out the scratch file:"
          READ (*,'(A60)') cfile
          CALL dump_out(cfile)
      ELSEIF ((prompt.EQ."R  ").or.(prompt.EQ."r  ")) THEN
          WRITE(*,*)
     *    " Enter the file name for reading in the scratch file:"
          READ (*,'(A60)') cfile
	  CALL dump_in ('#data/'//cfile)
      ELSEIF ((prompt.EQ."WG ").or.(prompt.EQ."wg ")) THEN
          aux=tout
          tout=timens
          CALL POST_OUT
          tout=aux
      ELSEIF ((prompt.EQ."C  ").or.(prompt.EQ."c  ")) THEN
          GOTO 4
      ELSEIF ((prompt.EQ."Q  ").or.(prompt.EQ."q  ")) THEN
5         WRITE(*,*) " Are You sure, You want to exit? (y,n)"
          READ(*,'(A)') prompt
          IF     ((prompt.EQ."Y  ").or.(prompt.EQ."y  ")) THEN
            STOP
          ELSEIF ((prompt.EQ."N  ").or.(prompt.EQ."n  ")) THEN
          ELSE
           GOTO 5
          END IF
      ELSE
        WRITE(*,'(" Invalid option! " )')
      END IF

      GOTO 3

2     WRITE(*,'(" Repeat Your choice, 
     *the variable was NOT saved ... ")')

3     GOTO 1

4     WRITE(*,10)

10    FORMAT(80('='))
      END

      SUBROUTINE CHECKER(i)
      INTEGER i
      LOGICAL bexists
	CHARACTER CFILE*(25)
	DATA CFILE/'#data/checker.dat'/

      i = 0
      INQUIRE (FILE=CFILE, EXIST = bexists) 
      IF (bexists)  THEN
      OPEN (UNIT=1,FILE=CFILE)
      READ (1,*) i
      BACKSPACE 1
      WRITE (1,'(I1)') 0
      CLOSE(1)
      ELSE
      RETURN
      END IF

      END


      SUBROUTINE ReCompute_RE(NY0)
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)

      PARAMETER (NNARR=299,NNAB=21,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

      DOUBLE PRECISION VRPARM,NY,NY0
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      IF (IPRECA.NE.4) THEN
       DO ILEV=NLMIN,NLMAX
        ISETLV=2
        CALL SETLEV(ISETLV)
        DO I=1,NA
         DWORK(KST1+I-1) = DWORK(KST1+I-1)*NY/NY0
        END DO
       END DO
      ELSE
        WRITE(*,*) "Diffusion operatod can not be rebuilded..."
        PAUSE
      END IF
C
      END