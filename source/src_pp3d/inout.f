      SUBROUTINE POST_OUT

!       USE module_levelset, ONLY: levelset,normal
!       USE LinScalar, ONLY: LevSet
      USE Transport_UxyzP_Q2P1, ONLY: QuadSc

      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
C *** Standard COMMON blocks
      PARAMETER       (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD,LVOL
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /EXTRA/ KLKX(NNLEV),KLKY(NNLEV),KLKZ(NNLEV),KLDK(NNLEV),
     *               LRHS,LRHS1,LRHS2,LRHS3,LNUT,KLNUT(NNLEV),
     *               LISEP,LIAUX,LINOD,LJNOD,LAEDGE,
     *               LF0U1,LF0U2,LF0U3,LF1U1,LF1U2,LF1U3
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGIEL/  KLINT(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGPER/  KKERI(NNLEV),KKERP(NNLEV),KKERV(NNLEV),
     *                NVERPER(NNLEV),KKROS(NNLEV),LKERI,LKERP,LKROS
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IGID,IGMV,IFINIT
C *** user COMMON block
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
      COMMON /DELTA/ KLDELT,iLDELT,LDSTRS
      COMMON /LESPAR/D_SMAG,I_SMAG
C
      INTEGER OUTPUTS(2),PROFILES(4)
      DATA PROFILES/1,1,1,1/
      SAVE 
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL UE,EM30,EM31,E030,E031
C	
      ILEV=NLEV
      CALL SETLEV(2)

      IF (PROFILES(2).EQ.1) THEN
C     Extracting the turbulent viscosities by Smagorinsky LES
C-----------------------------------------------------------------------
      CALL ZTIME(TTT0)
C
!       IF (IELT.EQ.0 .OR. IELT.EQ.2) 
!      *   CALL GETGRAD(DWORK(KU1),DWORK(KU2),DWORK(KU3),KWORK(L(LVERT)),
!      *        KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),
!      *        DWORK(L(LDSTRS)),EM31)
!       IF (IELT.EQ.1 .OR. IELT.EQ.3) 
!      *   CALL GETGRAD(DWORK(KU1),DWORK(KU2),DWORK(KU3),KWORK(L(LVERT)),
!      *        KWORK(L(LAREA)),KWORK(L(LEDGE)),DWORK(L(LCORVG)),
!      *        DWORK(L(LDSTRS)),EM30)
! C
!       DO IEL=1,NEL
!        DELTA=VWORK(L(LVOL)+IEL-1)**0.333333333D0
!        DWORK(L(LNUT)+IEL-1)=((D_SMAG*DELTA)**2D0)*DWORK(L(LDSTRS)+IEL-1)
!       END DO
C
      CALL ZTIME(TTT1)
      TTLES=TTLES+TTT1-TTT0
C-----------------------------------------------------------------------
      END IF
C
      CALL ZNEW(NEL,1,LW1 ,'WU  ')
      CALL ZNEW(NEL,1,LW2 ,'WV  ')
      CALL ZNEW(NEL,1,LW3 ,'WW  ')

      CALL ZNEW(NVT,1,LAUX,'AUX ')
      CALL ZNEW(NVT,2,LVU ,'VU  ')
      CALL ZNEW(NVT,2,LVV ,'VV  ')
      CALL ZNEW(NVT,2,LVW ,'VW  ')
      CALL ZNEW(NVT,2,LVP ,'VP  ')
      CALL ZNEW(NVT,2,LVNU,'VNU ')
      CALL ZNEW(NVT,2,LVRT,'VORT')
C
      CALL ZTIME(TTT0)
C------------------------------------------------------------------------

      IF (PROFILES(4).EQ.1) THEN
      ! Interpolate the vorticities to the conforming space
      IF (IELT.EQ.0 .OR. IELT.EQ.2)
     * CALL VORT(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *           DWORK(L(LW1)),DWORK(L(LW2)),DWORK(L(LW3)),
     *           KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *           DWORK(L(LCORVG)),EM31)
C
      IF (IELT.EQ.1 .OR. IELT.EQ.3)
     * CALL VORT(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *          DWORK(L(LW1)),DWORK(L(LW2)),DWORK(L(LW3)),
     *          KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *          DWORK(L(LCORVG)),EM30)
C
      CALL INTPV(DWORK(L(LW1)),VWORK(L(LVU)),KWORK(L(LAUX)),
     *           VWORK(L(LVOL)),KWORK(L(LVERT)))
      CALL INTPV(DWORK(L(LW2)),VWORK(L(LVV)),KWORK(L(LAUX)),
     *           VWORK(L(LVOL)),KWORK(L(LVERT)))
      CALL INTPV(DWORK(L(LW3)),VWORK(L(LVW)),KWORK(L(LAUX)),
     *           VWORK(L(LVOL)),KWORK(L(LVERT)))
C
      DO IVT=1,NVT
       VWORK(L(LVRT)+IVT-1)=SQRT(VWORK(L(LVU)+IVT-1)**2+
     * VWORK(L(LVV)+IVT-1)**2+VWORK(L(LVW)+IVT-1)**2) 
      END DO
      END IF
C
C------------------------------------------------------------------------
      IF (PROFILES(1).EQ.1) THEN
      ! Interpolate velocities to the conforming space
!       CALL INTUV(normal%l(ILEV)%val(1),normal%l(ILEV)%val(NAT+1),
!      *     normal%l(ILEV)%val(2*NAT+1),
      CALL INTUV(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *           VWORK(L(LVU)),VWORK(L(LVV)),VWORK(L(LVW)),
     *           KWORK(L(LAUX)),NVT,NEL,NVBD,KWORK(L(LAREA)),
     *           KWORK(L(LVERT)),KWORK(L(LNPR)),
     *           KWORK(L(LVBD)),DWORK(L(LCORVG)),UE)
      CALL PERGMV(VWORK(L(LVU)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
      CALL PERGMV(VWORK(L(LVV)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
      CALL PERGMV(VWORK(L(LVW)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
      END IF
C------------------------------------------------------------------------
      IF (PROFILES(2).EQ.1) THEN
      ! Interpolate the numerical diffussion fluxes to the conforming space
!         DO IVT=1,NVT
!          VWORK(L(LVNU)+IVT-1)=REAL(DWORK(L(LNUT)+IVT-1))
!         END DO

!         DO IVT=1,NVT
!          VWORK(L(LVNU)+IVT-1)=REAL(QuadSc%valU(ILEV)%x(ivt))
!         END DO

!       CALL INTPV(DWORK(L(LNUT)),VWORK(L(LVNU)),VWORK(L(LAUX)),
!      *           VWORK(L(LVOL)),KWORK(L(LVERT)))
!       CALL PERGMV(VWORK(L(LVNU)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
!
!       CALL INTNU(levelset%l(ILEV)%val,VWORK(L(LVNU)),
!      *           KWORK(L(LAUX)),NVT,NEL,NVBD,KWORK(L(LAREA)),
!      *           KWORK(L(LVERT)),KWORK(L(LNPR)),
!      *           KWORK(L(LVBD)),DWORK(L(LCORVG)),UE)
!       CALL PERGMV(VWORK(L(LVNU)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
      END IF
C------------------------------------------------------------------------
      IF (PROFILES(3).EQ.1) THEN
      ! Interpolate the pressure to the conforming space
      CALL INTPV(DWORK(KP),VWORK(L(LVP)),VWORK(L(LAUX)),
     *           VWORK(L(LVOL)),KWORK(L(LVERT)))
      CALL PERGMV(VWORK(L(LVP)),KWORK(L(KKERV(ILEV))),NVERPER(ILEV))
      END IF
C------------------------------------------------------------------------
      CALL ZTIME(TTT1)
      TTPROJ=TTPROJ+TTT1-TTT0
C
C========================================================================
      ILEV=MAX(IGMV,iGiD)
      CALL SETLEV(2)
      OUTPUTS(1)=IGMV
      OUTPUTS(2)=iGiD
C
      CALL POSTPROC(VWORK(L(LVU)),VWORK(L(LVV)),VWORK(L(LVW)),
     *           VWORK(L(LVP)),VWORK(L(LVNU)),VWORK(L(LVRT)),
     *           DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *           NEL,NVT,TIMENS,OUTPUTS,PROFILES)
C
C========================================================================
C
      CALL ZDISP(0,LW1 ,'WU  ')
      CALL ZDISP(0,LW2 ,'WV  ')
      CALL ZDISP(0,LW3 ,'WW  ')
C
      CALL ZDISP(0,LAUX,'AUX ')
      CALL ZDISP(0,LVU ,'VU  ')
      CALL ZDISP(0,LVV ,'VV  ')
      CALL ZDISP(0,LVW ,'VW  ')
      CALL ZDISP(0,LVP ,'VP  ')
      CALL ZDISP(0,LVNU,'VNU ')
      CALL ZDISP(0,LVRT,'VORT')
C
      ILEV=NLEV
      CALL SETLEV(2)
C
      END 
C
C
C
      SUBROUTINE DUMP_OUT(CFILE)
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILE*15
C *** Standard COMMON blocks
      PARAMETER       (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      SAVE 
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
	
      ILEV=NLEV
      CALL SETLEV(2)
	CALL WRITE_DUMP(TIMENS,DWORK(KU1),DWORK(KU2),
     *                DWORK(KU3),DWORK(KP),CFILE)

      END 


      SUBROUTINE DUMP_IN(CFILE)
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILE*(*)
C *** Standard COMMON blocks
      PARAMETER       (NNARR=299,NNLEV=9,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      SAVE 
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))

      ILEV=NLEV
      CALL SETLEV(2)
      CALL READ_DUMP(TIMENS,DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KP),CFILE)
      END 

      ! Unformatted output of data

      SUBROUTINE READ_DUMP(time,du_x,du_y,du_z,dp,cdump)
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      CHARACTER cdump*(*),CIFile*15
      REAL*8 du_x(*),du_y(*),du_z(*),dp(*)
      REAL*8 time
      INTEGER i,IFl

      IF (myid.EQ.MASTER) RETURN

      IFl=LEN(TRIM(ADJUSTL(cdump)))
      CIFile=TRIM(ADJUSTL(cdump))
      IF (myid.lt.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I1,I1)') '_',0,myid
      IF (myid.ge.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I2)') '_',myid
      OPEN (UNIT=1,FILE=CIFile,STATUS="OLD",FORM="FORMATTED")

      ! Hydrodynamic variables
      READ(1,'(G12.6)') timens
      DO I=1,NAT
       READ(1,'(3(1XG12.6))') du_x(i),du_y(i),du_z(i)
      END DO
      DO I=1,NEL
       READ(1,'(G12.6)') dp  (i)
      END DO
!       READ(UNIT=1) time,(du_x(i),i=1,NAT),
!      *                   (du_y(i),i=1,NAT),
!      *                   (du_z(i),i=1,NAT),
!      *                   (dp  (i),i=1,NEL)
      IF (timens.LE.1d-10) timens=time

      CLOSE(1)

      END


      ! Unformatted output of data

      SUBROUTINE WRITE_DUMP(time,du_x,du_y,du_z,dp,cdump)
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      REAL*8 du_x(*),du_y(*),du_z(*),dp(*)
      REAL*8 TIME
      CHARACTER cdump*15
      CHARACTER COFile*15
      INTEGER i,itwx,ifilen
      COMMON /NSSAV/  INSAV,INSAVN
C
      DATA ifilen/0/

      IF (myid.EQ.MASTER) RETURN

      IF (cdump.NE."DX             ") THEN
       OPEN (UNIT=2,FILE='#data/'//cdump,FORM="FORMATTED")
      ELSE
       ifilen=ifilen+1
       itwx=MOD(ifilen+insavn-1,insavn)+1
       COFile='#ns/DX        '
       IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1,A1)') 0,itwx,'_'
       IF (itwx.ge.10) WRITE(COFile(7:9),'(I2,A1)') itwx,'_'
       IF (myid.lt.10) WRITE(COFile(10:11),'(I1,I1)') 0,myid
       IF (myid.ge.10) WRITE(COFile(10:11),'(I2)') myid
       OPEN (UNIT=2,FILE=COFile,FORM="FORMATTED")
      END IF

      ! Hydrodynamic variables
!       WRITE(UNIT=2) timens,(du_x(i),i=1,NAT),
!      *                       (du_y(i),i=1,NAT),
!      *                       (du_z(i),i=1,NAT),
!      *                       (dp  (i),i=1,NEL)
      WRITE(2,'(G12.6)') timens
      DO I=1,NAT
       WRITE(2,'(3(1XG12.6))') du_x(i),du_y(i),du_z(i)
      END DO
      DO I=1,NEL
       WRITE(2,'(G12.6)') dp  (i)
      END DO

      CLOSE(2)

      END
C
C
C
      SUBROUTINE read_prol(time,ilevH,iintu,cdump,
     *                     DUL,nuL,npL,KAREAL,KADJL,
     *                     DUH,nuH,npH,KAREAH,KADJH)
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
      CHARACTER cdump*15,CIFile*15
      REAL*8 DUL(*),DUH(*)
      INTEGER KAREAL(6,*),KAREAH(6,*),KADJL(6,*),KADJH(6,*)
      INTEGER IFl
C
      IF (myid.EQ.MASTER) RETURN

      IFl=LEN(TRIM(ADJUSTL(cdump)))
      CIFile=TRIM(ADJUSTL(cdump))
      IF (myid.lt.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I1,I1)') '_',0,myid
      IF (myid.ge.10) WRITE(CIFile(iFl+1:iFl+3),'(A1,I2)') '_',myid

      ku1L = 1
      ku2L = ku1L + nuL
      ku3L = ku2L + nuL
      kpL  = ku3L + nuL
      ku1H = 1
      ku2H = ku1H + nuH
      ku3H = ku2H + nuH
      kpH  = ku3H + nuH
C
      OPEN (UNIT=1,FILE=CIFile,FORM="FORMATTED")
      ! Hydrodynamic variables
!       READ(UNIT=1) time,(DUL(ku1L+i-1),i=1,nuL),
!      *                   (DUL(ku2L+i-1),i=1,nuL),
!      *                   (DUL(ku3L+i-1),i=1,nuL),
!      *                   (DUL(kpL +i-1),i=1,npL)
      READ(1,'(G12.6)') timens
      DO I=1,nuL
       READ(1,'(3(1XG12.6))') DUL(ku1L+i-1),DUL(ku2L+i-1),DUL(ku3L+i-1)
      END DO
      DO I=1,npL
       READ(1,'(G12.6)') DUL(kpL +i-1)
      END DO

      IF (timens.LE.1d-10) timens=time
      CLOSE(1)
C
      CALL PROLU2(DUL(ku1L),DUL(ku2L),DUL(ku3L),DUL(kpL),
     *            DUH(ku1H),DUH(ku2H),DUH(ku3H),DUH(kpH),
     *            KAREAL,KADJL,nuL,npL,
     *            KAREAH,KADJH,nuH,npH,iintu)
C
      END
C
C
C
      SUBROUTINE POSTPROC(vu,vv,vw,vp,vnu,vrt,xyz,
     *           kvert,nel,nvt,t,o,p)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,D-H,O-U,W-Z),LOGICAL(B)
        REAL*8 xyz(3,*),t
        REAL*4 vu(*),vv(*),vw(*),vp(*),vnu(*),vrt(*)
        INTEGER kvert(8,*)
        CHARACTER cf*(20),cmm*(20),cm*(20)
        INTEGER mf,nel,nvt,o(2),p(4)
        DATA i_Out /0/,mf /1/
C
C
       i_Out=i_Out+1
       IF (o(1).GT.0) THEN     ! gmv

       cf="#gmv/               "
       IF(myid.lt.10) WRITE(cf(6:12),'(A4,I1,I1,A1)') 'res_',0,myid,'_'
       IF(myid.ge.10) WRITE(cf(6:12),'(A4,I2,A1)') 'res_',myid,'_'
C
       IF     ((i_Out.GE.0   ).AND.(i_Out.LT.10   )) THEN
        WRITE(cf(13:15),'(A3)') "000"
        WRITE(cf(16:20),'(I1,A4)') i_Out,".gmv"
       ELSEIF ((i_Out.GE.10  ).AND.(i_Out.LT.100  )) THEN
        WRITE(cf(13:14),'(A2)') "00"
        WRITE(cf(15:20),'(I2,A4)') i_Out,".gmv"
       ELSEIF ((i_Out.GE.100 ).AND.(i_Out.LT.1000 )) THEN
        WRITE(cf(13:13),'(A1)') "0"
        WRITE(cf(14:20),'(I3,A4)') i_Out,".gmv"
       ELSEIF ((i_Out.GE.1000).AND.(i_Out.LT.10000)) THEN
        WRITE(cf(13:20),'(I4,A4)') i_Out,".gmv"
       ELSEIF (i_Out.GE.10000)                       THEN
        STOP
       END IF

       cmm="msh                 "
       IF(myid.lt.10) WRITE(cmm(4:10),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
       IF(myid.ge.10) WRITE(cmm(4:10),'(A1,I2,A4)') '_',myid,".gmv"

       cm="#gmv/msh             "
       IF(myid.lt.10) WRITE(cm(9:15),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
       IF(myid.ge.10) WRITE(cm(9:15),'(A1,I2,A4)') '_',myid,".gmv"

!       WRITE(*,*) myid,cf,cm,cmm

       if (myid.ne.MASTER) THEN
       IF (i_Out.EQ.1) CALL XGMVMS(mf,cm,nel,nvt,kvert,xyz)
        CALL XGMV3D(mf,cf,cmm,nel,nvt,vu,vv,vw,vp,vnu,vrt,t,p)
       END IF
       END IF
       IF (o(2).GT.0) THEN     ! GiD
       cf="#GiD/flavia.res"
       cm="#GiD/flavia.msh"

       IF (i_Out.EQ.1) CALL XGIDMS(mf,cm,cf,nel,nvt,kvert,xyz)
        CALL XGID3D(mf,cf,nel,nvt,vu,vv,vw,vp,vnu,vrt,i_Out,p)
       END IF

       END
