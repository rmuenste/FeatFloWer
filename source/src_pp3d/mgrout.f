************************************************************************
*    M011   (edited from M010)                                         *
************************************************************************
* FINIITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* M011                                                                 *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          multigrid iteration                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LL21 , LLC1                            *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting address of vectors containing the           *
* DB       R*8    solution and the right hand side, DD is used as      *
* DD       R*8    auxiliary vector only                                *
* KOFFX           The actual starting address of DX on level ILEV      *
* KOFFB           is DX(1+KOFFX(ILEV)) (analogously for DB and DD)     *
* KOFFD           Total space required for all vectors is              *
*                 KNEQ(NLMIN)+...+KNEQ(NLMAX)                          *
*                 DX(1+KOFFX(NLMAX)) contains initial solution         *
*                 DB(1+KOFFB(NLMAX)) contains right hand side          *
* KNEQ     I*4    Number of equations for all levels                   *
* NLMAX    I*4    Iteration uses levels NLMIN to NLMAX,                *
* NLMIN    I*4    NLMAX  is the finest level                           *
* NIT      I*4    Maximum number of iterations                         *
*                 Iteration completed after reaching the finest level  *
* EPS1      R*8   Desired precision                                    *
*                 Stop if !!DEFn!! < EPS1 !!DEF0!!                     *
* EPS2      R*8   Desired precision                                    *
*                 Stop if !!DEF!! < EPS2                               *
* KPRSM    I*4    Number of pre -smoothing steps for all levels        *
* KPOSM    I*4    Number of post-smoothing steps for all levels        *
* ICYCLE   I*4    <0: special cycle types (not yet implemented)        *
*                 =0: F-Cycle                                          *
*                 =1: V-Cycle                                          *
*                 =2: W-Cycle                                          *
*                 >2: Cycle of higher order                            *
* DAX      SUBR   CALL DAX(DX,DAX,NEQ,A1,A2)                           *
*                 Returns  DAX := A1*A*DX+A2*DAX                       *
* DPROL    SUBR   DPROL(DX1,DX2)                                       *
*                 Returns  DX2 := Prolongation(DX1)  to higher level   *
* DSTEP    SUBR   DPROL(DX,DD,DB,DSTEPP)                               *
*                 Returns  DSTEPP := optimal step size for correction  *
* DREST    SUBR   DREST(DD1,DD2)                                       *
*                 Returns  DD2 := Restriction(DD1)  to lower level     *
* DPRSM    SUBR   DPRSM(DX,DB,DD,NEQ,NPRSM)                            *
*                 Returns  DX  after  NPRSM:=KPRSM(ILEV)               *
*                 pre-smoothing steps, DD  is used as auxiliary vector *
* DPOSM    SUBR   Same as above, used for post-smoothing               *
* DEX      SUBR   DEX(DX,DB,DD,NEQ)                                    *
*                 Returns exact solution                               *
* DBC      SUBR   DBC(DX,NEQ)                                          *
*                 Copies boundary data onto components of  DX          *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector on DX(KOFFX(NLMAX))                  *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
* RHOLMG   R*8    Multigrid convergence rate                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE  M011 (DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DEF,DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXA,DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,
     *                  BMGEND)
C
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      SAVE 
C
C
      SUB='M011  '
      IF (ICHECK.GE.997) CALL OTRC('M011  ','05/13/91')
      IER=0
C
C
      BMGEND=.FALSE.
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0
C
      BTIME=IMTIME.GT.0
      IF (BTIME) THEN
       IF (IMTIME.EQ.1) THEN
        TTMG=0D0
        TTS=0D0
        TTE=0D0
        TTD=0D0
        TTP=0D0
        TTR=0D0
       ENDIF
       CALL ZTIME(TTMG0)
      ENDIF
C
C
      NIT0=MAX(ITE,0)
      IF (ICYCLE.LT.0) THEN
       CALL WERR(-181,'M011  ')
       GOTO 99999
      ENDIF
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
       CALL WERR(-182,'M011  ')
       GOTO 99999
      ENDIF
      ILEV=NLMAX
C
C *** special case - zero solution
      IMGComplete=0                                                  ! PARALLEL
      ISTARTIMGComplete=0                                                  ! PARALLEL
      IF (myid.NE.MASTER) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LLI1(DB(1+KOFFB(NLMAX)),KNEQ(NLMAX),DMAX,INDMAX)
       IF (DMAX.LE.1D-20) THEN
        CALL LCP1(DB(1+KOFFB(NLMAX)),DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
        RHOLMG=0D0
        IMGComplete = 1
        ISTARTIMGComplete = 1
!        write(*,*) myid
       ENDIF
!       write(*,*) myid,IMGComplete
       CALL COMM_MGComplete(IMGComplete)                             ! PARALLEL
       IF (IMGComplete.EQ.1) GOTO 99999                              ! PARALLEL
      ENDIF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
C *** level counts
!       WRITE(*,*) "subs",myid
      KIT0(NLMAX)=1
      DO 2 ILEV=NLMIN+1,NLMAX-1
      IF (ICYCLE.EQ.0) THEN
       KIT0(ILEV)=2
      ELSE
       KIT0(ILEV)=ICYCLE
      ENDIF
2     CONTINUE
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
C    ----- calculate the initial defect in 3 steps (due to communication) ----
        CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),  ! PARALLEL
     *           -1D0,0D0)                                           ! PARALLEL
        IF (VARIABLE.EQ.1) THEN                                      ! PARALLEL
         CALL CommSum(DD(1+KOFFD(NLMAX)),NLMAX)                      ! PARALLEL
        ELSE                                                         ! PARALLEL
          CALL CommPres(DD(1+KOFFD(NLMAX)),DX(1+KOFFX(NLMAX)),NLMAX) ! PARALLEL
        END IF                                                       ! PARALLEL
        CALL LLC1 (DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),            ! PARALLEL
     *             KNEQ(NLMAX),1d0,1D0)                              ! PARALLEL
C    ------------------------------------------------------------------
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
!       if (myid.eq.1) write(*,*) myid,DEF
       DEFOLD=DEF
C *** FD  is considered as initial defect
       FD=DEF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       IF (DEF.LE.EPS2.AND..NOT.BREL) THEN
        ITE=0
        IMGComplete = 1                                              ! PARALLEL
       ENDIF

!      WRITE(*,*) "myid",myid,IMGComplete
      CALL COMM_MGComplete(IMGComplete)                              ! PARALLEL
      IF (IMGComplete.EQ.1) GOTO 1000                                ! PARALLEL
      ENDIF
C
C
C *** Start multigrid iteration
      DO 100 ITE=1,NIT

      IF (ISTARTIMGComplete.EQ.1) IMGComplete = 1                    ! PARALLEL
C
C *** initialize level counts for all levels
      DO 101 ILEV=NLMIN,NLMAX
101   KIT(ILEV)=KIT0(ILEV)
C
      ILEV=NLMAX
110   IF (ILEV.NE.NLMIN) THEN
C
        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0) 
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))
C
        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF
C
        IF (BTIME) CALL ZTIME(TTD0)
C    ----- calculate defect in three steps (due to communication) ----
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),     ! PARALLEL
     *           -1D0,0D0)                                           ! PARALLEL
        IF (VARIABLE.EQ.1) THEN                                      ! PARALLEL
         CALL CommSum(DD(1+KOFFD(ILEV)),ILEV)                        ! PARALLEL
        ELSE                                                         ! PARALLEL
         CALL CommPres(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),ILEV)     ! PARALLEL
        END IF                                                       ! PARALLEL
        CALL LLC1 (DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),              ! PARALLEL
     *             KNEQ(ILEV),1d0,1D0)                               ! PARALLEL
       CALL LL21(DD(1+KOFFD(ILEV)),KNEQ(ILEV),DEF)
!       if (myid.eq.2.and.ilev.eq.nlmax) write(*,*) "1",DEF
!        STOP
C    ------------------------------------------------------------------
C
        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF
C
       ILEV=ILEV-1
C
C ***  restriction of defect
       IF (BTIME) CALL ZTIME(TTR0)
       CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
       IF (VARIABLE.EQ.1) CALL CommSumHalf(DB(1+KOFFB(ILEV)),ILEV)   ! PARALLEL
C
       IF (BTIME) THEN
        CALL ZTIME(TTR1)
        TTR=TTR+TTR1-TTR0
       ENDIF
C
C ***  choose zero as initial vector on lower level
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
       CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       GOTO 110
C
      ENDIF
C
C *** exact solution on lowest level
C
      IF (BTIME) CALL ZTIME(TTE0)
      CALL LCP1(DB(1+KOFFB(NLMIN)),DX(1+KOFFX(NLMIN)),KNEQ(NLMIN))
      CALL CommDex(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),            ! PARALLEL
     *             DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),RHOLMG)            ! PARALLEL
C
      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF
C
130   IF (ILEV.NE.NLMAX) THEN
       ILEV=ILEV+1
C
C *** DPROL  returns  DD:=PROL(DX)
      IF (BTIME) CALL ZTIME(TTP0)
      CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
      IF (VARIABLE.EQ.1) CALL CommSumHalf(DD(1+KOFFD(ILEV)),ILEV)      ! PARALLEL
      CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))
C
      CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *           DB(1+KOFFB(ILEV)),KNEQ(ILEV),DSTEPP)
      CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),
     *          DSTEPP,1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTP1)
       TTP=TTP+TTP1-TTP0
      ENDIF
C
C
C    ----- calculate defect in three steps (due to communication) ----
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),  ! PARALLEL
     *           -1D0,0D0)                                           ! PARALLEL
        IF (VARIABLE.EQ.1) THEN                                      ! PARALLEL
         CALL CommSum(DD(1+KOFFD(ILEV)),ILEV)                      ! PARALLEL
        ELSE                                                         ! PARALLEL
         CALL CommPres(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),ILEV)  ! PARALLEL
        END IF                                                       ! PARALLEL
        CALL LLC1 (DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),            ! PARALLEL
     *             KNEQ(ILEV),1d0,1D0)                              ! PARALLEL
       CALL LL21(DD(1+KOFFD(ILEV)),KNEQ(ILEV),DEF)
!       if (myid.eq.1) write(*,*) "2",DEF
C    ------------------------------------------------------------------

C *** Post-smoothing
      IF (BTIME) CALL ZTIME(TTS0)
      IF (KPOSM(ILEV).GT.0) 
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *            DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPOSM(ILEV))
C
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTS0
      ENDIF
C
       KIT(ILEV)=KIT(ILEV)-1
       IF (KIT(ILEV).EQ.0) THEN
        IF (ICYCLE.EQ.0) THEN
         KIT(ILEV)=1
        ELSE
         KIT(ILEV)=KIT0(ILEV)
        ENDIF
        GOTO 130
       ELSE
        CALL COMM_MGComplete(IMGComplete)                            ! PARALLEL
        IF (ISTARTIMGComplete.EQ.1) IMGComplete = 1                  ! PARALLEL
        GOTO 110
       ENDIF
C
      ENDIF
C
140   IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
C    ----- calculate defect in three steps (due to communication) ----
        CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),  ! PARALLEL
     *           -1D0,0D0)                                           ! PARALLEL
        IF (VARIABLE.EQ.1) THEN                                      ! PARALLEL
         CALL CommSum(DD(1+KOFFD(NLMAX)),NLMAX)                      ! PARALLEL
        ELSE                                                         ! PARALLEL
         CALL CommPres(DD(1+KOFFD(NLMAX)),DX(1+KOFFX(NLMAX)),NLMAX)  ! PARALLEL
        END IF                                                       ! PARALLEL
        CALL LLC1 (DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),            ! PARALLEL
     *             KNEQ(NLMAX),1d0,1D0)                              ! PARALLEL
C    ------------------------------------------------------------------
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
!       if (myid.eq.2) write(*,*) "3",DEF
       IF (MT.GE.9) WRITE (6,10001) ITE,DEF
       IF (BMSG2) THEN
        MT=MT0
        WRITE (CPARAM,'(I15,D25.16)') ITE,DEF
        CALL OMSG(73,'M011  ')
        WRITE (MPROT,10001) ITE,DEF
        WRITE (*,10001) ITE,DEF
        MT=0
       ENDIF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
C=======================================================================
C ***  Unexpected STOP !!!
C=======================================================================
       IF ((DEF/DEFOLD.GT.1D4).OR.(DEF/FD.GT.5D3)) THEN
!        BMGEND=.TRUE.
         IMGComplete = 1                                             ! PARALLEL
!      WRITE(*,*) "critical:",myid,IMGComplete
 !       GOTO 1000
       ENDIF
C
       DEFOLD=DEF
       IF (BREL) THEN
        IF (DEF.LE.FD*EPS1.AND.DEF.LE.EPS2.AND.ITE.GE.NIT0) THEN
!          if (variable.eq.2) write(*,*) "ok",myid
         IMGComplete = 1                                             ! PARALLEL
        END IF
        IF (DEF.LE.1d-14) THEN ! This will be the precision of MG (EPS2)
         IMGComplete = 1                                             ! PARALLEL
        END IF
       ELSE
        IF (DEF.LE.EPS2) THEN
!          if (variable.eq.2) write(*,*) "ok",myid
         IMGComplete = 1                                             ! PARALLEL
        END IF
       ENDIF
      ENDIF
C
      IF (ITE.GE.NIT) IMGComplete = 1                                ! PARALLEL
!       if (variable.eq.2) write(*,*) "not ok",IMGComplete,NIT0
!      if (myid.eq.1) write(*,*) "converged?",IMGComplete,DEF,brel
      CALL COMM_MGComplete(IMGComplete)                              ! PARALLEL
      IF (IMGComplete.EQ.1) GOTO 1000                                ! PARALLEL
      IF (ISTARTIMGComplete.EQ.1) IMGComplete = 1                    ! PARALLEL
C
100   CONTINUE
C
      ITE=NIT
      IF (IDEFMG.EQ.1) THEN
       MT=MT0
       RHOLMG=DEF/FD
       RHOLMG=RHOLMG**(1D0/DBLE(NIT))
       IF ((RHOLMG).GE.1D0) BMGEND=.TRUE.
       IF (MT.GE.9) WRITE (6,10000) NIT,DEF,FD,DEF/FD,RHOLMG
       IF (BMSG2) THEN
        WRITE (CPARAM,'(I15,4D25.16)') NIT,DEF,FD,DEF/FD,RHOLMG
        CALL OMSG(72,'M011  ')
        WRITE (MPROT,10000) NIT,DEF,FD,DEF/FD,RHOLMG
        WRITE (MPROT,*) 'IER=1 in M011'
       ENDIF
      ELSE
       RHOLMG=DBLE(NIT)
      ENDIF
C
      GOTO 99999
C
1000  IER=0
      MT=MT0
      RHOLMG=0.D0
      IF (FD.GE.1D-70) RHOLMG=(DEF/FD)**(1D0/DBLE(ITE))
!      if (myid.eq.1) write(*,*) RHOLMG,ITE
      IF (MT.GE.9) WRITE (6,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
       CALL OMSG(72,'M011  ')
       WRITE (CPARAM,'(D25.16)')  RHOLMG
       CALL OMSG(76,'M011  ')
       WRITE (MPROT,10002) ITE,DEF,FD,DEF/FD,RHOLMG
      ENDIF
C
99999 MT=MT0 
      IF (BTIME) THEN
       CALL ZTIME(TTMG1)
       TTMG=TTMG+TTMG1-TTMG0
      ENDIF
C
10000 FORMAT (//' CONVERGENCE OF MG-METHOD FAILED'/
     *        ' NORM OF DEFECT AFTER',I6,' ITERATIONS',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
10001 FORMAT (' ITERATION',I6,5X,'  RESIDUUM =',D12.3)
10002 FORMAT (' NUMBER OF ITERATIONS    ',I10/
     *        ' NORM OF DEFECT          ',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
C
C
C
!      if (myid.eq.2) write(*,*) "ok...."
      if (ITE.EQ.NIT.AND.variable.eq.1)
     *    write(*,*) myid, "maximum iterations"
!      if (variable.eq.2)write(*,*) myid, DEF
      END
C
C
C
************************************************************************
      SUBROUTINE   PROLU2  (DU1,DV1,DW1,DP1,DU2,DV2,DW2,DP2,
     *                     KAREA1,KADJ1,NAT1,NEL1,
     *                     KAREA2,KADJ2,NAT2,NEL2,IINTU)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DW1,DP1) 
*                to the fine grid vector (DU2,DV2,DW2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DW1,DP1           - coarse grid vector
*      KAREA1,KADJ1,NAT1,NEL1,  - data of the coarse grid
*      KAREA2,KADJ2,NAT2,NEL2,  - data of the fine grid
*
*    Output:
*      DU2,DV2,DW2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NNAE=6)
      DIMENSION DU1(1),DV1(1),DW1(1),DP1(1),DU2(1),DV2(1),DW2(1),DP2(1),
     *          KAREA1(NNAE,*),KAREA2(NNAE,*),
     *          KADJ1(NNAE,*),KADJ2(NNAE,*)
      SAVE
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.1) THEN
       A1= 11D0/24D0
       A2=  7D0/48D0
       A3=- 5D0/48D0
       A4=- 1D0/24D0
       A5= 11D0/24D0
       A6=  1D0/12D0
       A7=- 1D0/24D0
      ENDIF
C
      IF (IINTU.EQ.2) THEN
       A1= 0.500D0
       A2= 0.125D0
       A3=-0.125D0
       A4= 0.000D0
       A5= 0.500D0
       A6= 0.000D0
       A7= 0.000D0
      ENDIF
C
C
C
C *** Zero initialization of (DU2,DV2,DW2)
      CALL  LCL1 (DU2,NAT2)
      CALL  LCL1 (DV2,NAT2)
      CALL  LCL1 (DW2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DUH1=DU1(KAREA1(1,IEL1))
      DUH2=DU1(KAREA1(2,IEL1))
      DUH3=DU1(KAREA1(3,IEL1))
      DUH4=DU1(KAREA1(4,IEL1))
      DUH5=DU1(KAREA1(5,IEL1))
      DUH6=DU1(KAREA1(6,IEL1))

C
      DVH1=DV1(KAREA1(1,IEL1))
      DVH2=DV1(KAREA1(2,IEL1))
      DVH3=DV1(KAREA1(3,IEL1))
      DVH4=DV1(KAREA1(4,IEL1))
      DVH5=DV1(KAREA1(5,IEL1))
      DVH6=DV1(KAREA1(6,IEL1))
C
      DWH1=DW1(KAREA1(1,IEL1))
      DWH2=DW1(KAREA1(2,IEL1))
      DWH3=DW1(KAREA1(3,IEL1))
      DWH4=DW1(KAREA1(4,IEL1))
      DWH5=DW1(KAREA1(5,IEL1))
      DWH6=DW1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
      DP2(IELH5)=DPH
      DP2(IELH6)=DPH
      DP2(IELH7)=DPH
      DP2(IELH8)=DPH
C
      B11=A1*DUH1
      B12=A2*DUH1
      B13=A3*DUH1
      B14=A4*DUH1
      B15=A5*DUH1
      B16=A6*DUH1
      B17=A7*DUH1
C
      C11=A1*DUH2
      C12=A2*DUH2
      C13=A3*DUH2
      C14=A4*DUH2
      C15=A5*DUH2
      C16=A6*DUH2
      C17=A7*DUH2
C
      D11=A1*DUH3
      D12=A2*DUH3
      D13=A3*DUH3
      D14=A4*DUH3
      D15=A5*DUH3
      D16=A6*DUH3
      D17=A7*DUH3
C
      E11=A1*DUH4
      E12=A2*DUH4
      E13=A3*DUH4
      E14=A4*DUH4
      E15=A5*DUH4
      E16=A6*DUH4
      E17=A7*DUH4
C
      F11=A1*DUH5
      F12=A2*DUH5
      F13=A3*DUH5
      F14=A4*DUH5
      F15=A5*DUH5
      F16=A6*DUH5
      F17=A7*DUH5
C
      G11=A1*DUH6
      G12=A2*DUH6
      G13=A3*DUH6
      G14=A4*DUH6
      G15=A5*DUH6
      G16=A6*DUH6
      G17=A7*DUH6
      R1=B11+G14
      R2=C11+E14
      R3=D11+F14
      R4=C14+E11
      R5=F11+D14
      R6=B14+G11
      R7=B15+G17
      R8=B16+G16
      R9=B17+G15
C
C
C
C
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH1))=DU2(KAREA2(1,IELH1))+
     *                      (R1+C12+D13+
     *                      E13+F12)
       DU2(KAREA2(1,IELH2))=DU2(KAREA2(1,IELH2))+
     *                      (R1+C12+D12+
     *                      E13+F13)
       DU2(KAREA2(1,IELH3))=DU2(KAREA2(1,IELH3))+
     *                      (R1+C13+D12+
     *                      E12+F13)
       DU2(KAREA2(1,IELH4))=DU2(KAREA2(1,IELH4))+
     *                      (R1+C13+D13+
     *                      E12+F12)
      ELSE
       DU2(KAREA2(1,IELH1))=DUH1
       DU2(KAREA2(1,IELH2))=DUH1
       DU2(KAREA2(1,IELH3))=DUH1
       DU2(KAREA2(1,IELH4))=DUH1
      ENDIF
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH1))=DU2(KAREA2(2,IELH1))+
     *                      (B12+R2+D13+
     *                      F12+G13)
       DU2(KAREA2(5,IELH2))=DU2(KAREA2(5,IELH2))+
     *                      (B12+R2+D12+
     *                      F13+G13)
       DU2(KAREA2(5,IELH6))=DU2(KAREA2(5,IELH6))+
     *                      (B13+R2+D12+
     *                      F13+G12)
       DU2(KAREA2(2,IELH5))=DU2(KAREA2(2,IELH5))+
     *                      (B13+R2+D13+
     *                      F12+G12)
      ELSE
       DU2(KAREA2(2,IELH1))=DUH2
       DU2(KAREA2(5,IELH2))=DUH2
       DU2(KAREA2(5,IELH6))=DUH2
       DU2(KAREA2(2,IELH5))=DUH2
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH2))=DU2(KAREA2(2,IELH2))+
     *                      (B12+C12+R3+
     *                      E13+G13)
       DU2(KAREA2(5,IELH3))=DU2(KAREA2(5,IELH3))+
     *                      (B12+C13+R3+
     *                      E12+G13)
       DU2(KAREA2(5,IELH7))=DU2(KAREA2(5,IELH7))+
     *                      (B13+C13+R3+
     *                      E12+G12)
       DU2(KAREA2(2,IELH6))=DU2(KAREA2(2,IELH6))+
     *                      (B13+C12+R3+
     *                      E13+G12)
      ELSE
       DU2(KAREA2(2,IELH2))=DUH3
       DU2(KAREA2(5,IELH3))=DUH3
       DU2(KAREA2(5,IELH7))=DUH3
       DU2(KAREA2(2,IELH6))=DUH3
      ENDIF
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH3))=DU2(KAREA2(2,IELH3))+
     *                      (B12+R4+D12+
     *                      F13+G13)
       DU2(KAREA2(5,IELH4))=DU2(KAREA2(5,IELH4))+
     *                      (B12+R4+D13+
     *                      F12+G13)
       DU2(KAREA2(5,IELH8))=DU2(KAREA2(5,IELH8))+
     *                      (B13+R4+D13+
     *                      F12+G12)
       DU2(KAREA2(2,IELH7))=DU2(KAREA2(2,IELH7))+
     *                      (B13+R4+D12+
     *                      F13+G12)
      ELSE
       DU2(KAREA2(2,IELH3))=DUH4
       DU2(KAREA2(5,IELH4))=DUH4
       DU2(KAREA2(5,IELH8))=DUH4
       DU2(KAREA2(2,IELH7))=DUH4
      ENDIF
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DU2(KAREA2(2,IELH4))=DU2(KAREA2(2,IELH4))+
     *                      (B12+C13+R5+
     *                      E12+G13)
       DU2(KAREA2(5,IELH1))=DU2(KAREA2(5,IELH1))+
     *                      (B12+C12+R5+
     *                      E13+G13)
       DU2(KAREA2(5,IELH5))=DU2(KAREA2(5,IELH5))+
     *                      (B13+C12+R5+
     *                      E13+G12)
       DU2(KAREA2(2,IELH8))=DU2(KAREA2(2,IELH8))+
     *                      (B13+C13+R5+
     *                      E12+G12)
      ELSE
       DU2(KAREA2(2,IELH4))=DUH5
       DU2(KAREA2(5,IELH1))=DUH5
       DU2(KAREA2(5,IELH5))=DUH5
       DU2(KAREA2(2,IELH8))=DUH5
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DU2(KAREA2(1,IELH5))=DU2(KAREA2(1,IELH5))+
     *                      (R6+C12+D13+
     *                      E13+F12)
       DU2(KAREA2(1,IELH6))=DU2(KAREA2(1,IELH6))+
     *                      (R6+C12+D12+
     *                      E13+F13)
       DU2(KAREA2(1,IELH7))=DU2(KAREA2(1,IELH7))+
     *                      (R6+C13+D12+
     *                      E12+F13)
       DU2(KAREA2(1,IELH8))=DU2(KAREA2(1,IELH8))+
     *                      (R6+C13+D13+
     *                      E12+F12)
      ELSE
       DU2(KAREA2(1,IELH5))=DUH6
       DU2(KAREA2(1,IELH6))=DUH6
       DU2(KAREA2(1,IELH7))=DUH6
       DU2(KAREA2(1,IELH8))=DUH6
      ENDIF
C
C
      DU2(KAREA2(3,IELH1))=R7+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH2))=R7+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH3))=R7+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH4))=R7+C16+D17+
     *                     E16+F15
C
      DU2(KAREA2(6,IELH1))=R8+C15+D17+
     *                     E17+F15
      DU2(KAREA2(6,IELH2))=R8+C15+D15+
     *                     E17+F17
      DU2(KAREA2(6,IELH3))=R8+C17+D15+
     *                     E15+F17
      DU2(KAREA2(6,IELH4))=R8+C17+D17+
     *                     E15+F15
C
      DU2(KAREA2(3,IELH5))=R9+C15+D16+
     *                     E17+F16
      DU2(KAREA2(3,IELH6))=R9+C16+D15+
     *                     E16+F17
      DU2(KAREA2(3,IELH7))=R9+C17+D16+
     *                     E15+F16
      DU2(KAREA2(3,IELH8))=R9+C16+D17+
     *                     E16+F15
C
      B21=A1*DVH1
      B22=A2*DVH1
      B23=A3*DVH1
      B24=A4*DVH1
      B25=A5*DVH1
      B26=A6*DVH1
      B27=A7*DVH1
C
      C21=A1*DVH2
      C22=A2*DVH2
      C23=A3*DVH2
      C24=A4*DVH2
      C25=A5*DVH2
      C26=A6*DVH2
      C27=A7*DVH2
C
      D21=A1*DVH3
      D22=A2*DVH3
      D23=A3*DVH3
      D24=A4*DVH3
      D25=A5*DVH3
      D26=A6*DVH3
      D27=A7*DVH3
C
      E21=A1*DVH4
      E22=A2*DVH4
      E23=A3*DVH4
      E24=A4*DVH4
      E25=A5*DVH4
      E26=A6*DVH4
      E27=A7*DVH4
C
      F21=A1*DVH5
      F22=A2*DVH5
      F23=A3*DVH5
      F24=A4*DVH5
      F25=A5*DVH5
      F26=A6*DVH5
      F27=A7*DVH5
C
      G21=A1*DVH6
      G22=A2*DVH6
      G23=A3*DVH6
      G24=A4*DVH6
      G25=A5*DVH6
      G26=A6*DVH6
      G27=A7*DVH6
C
C
      S1=B21+G24
      S2=C21+E24
      S3=D21+F24
      S4=C24+E21
      S5=D24+F21
      S6=B24+G21
      S7=B25+G27
      S8=B26+G26
      S9=B27+G25
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH1))=DV2(KAREA2(1,IELH1))+
     *                      (S1+C22+D23+
     *                      E23+F22)
       DV2(KAREA2(1,IELH2))=DV2(KAREA2(1,IELH2))+
     *                      (S1+C22+D22+
     *                      E23+F23)
       DV2(KAREA2(1,IELH3))=DV2(KAREA2(1,IELH3))+
     *                      (S1+C23+D22+
     *                      E22+F23)
       DV2(KAREA2(1,IELH4))=DV2(KAREA2(1,IELH4))+
     *                      (S1+C23+D23+
     *                      E22+F22)
      ELSE
       DV2(KAREA2(1,IELH1))=DVH1
       DV2(KAREA2(1,IELH2))=DVH1
       DV2(KAREA2(1,IELH3))=DVH1
       DV2(KAREA2(1,IELH4))=DVH1
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH1))=DV2(KAREA2(2,IELH1))+
     *                      (B22+S2+D23+
     *                      F22+G23)
       DV2(KAREA2(5,IELH2))=DV2(KAREA2(5,IELH2))+
     *                      (B22+S2+D22+
     *                      F23+G23)
       DV2(KAREA2(5,IELH6))=DV2(KAREA2(5,IELH6))+
     *                      (B23+S2+D22+
     *                      F23+G22)
       DV2(KAREA2(2,IELH5))=DV2(KAREA2(2,IELH5))+
     *                      (B23+S2+D23+
     *                      F22+G22)
      ELSE
       DV2(KAREA2(2,IELH1))=DVH2
       DV2(KAREA2(5,IELH2))=DVH2
       DV2(KAREA2(5,IELH6))=DVH2
       DV2(KAREA2(2,IELH5))=DVH2
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH2))=DV2(KAREA2(2,IELH2))+
     *                      (B22+C22+S3+
     *                      E23+G23)
       DV2(KAREA2(5,IELH3))=DV2(KAREA2(5,IELH3))+
     *                      (B22+C23+S3+
     *                      E22+G23)
       DV2(KAREA2(5,IELH7))=DV2(KAREA2(5,IELH7))+
     *                      (B23+C23+S3+
     *                      E22+G22)
       DV2(KAREA2(2,IELH6))=DV2(KAREA2(2,IELH6))+
     *                      (B23+C22+S3+
     *                      E23+G22)
      ELSE
       DV2(KAREA2(2,IELH2))=DVH3
       DV2(KAREA2(5,IELH3))=DVH3
       DV2(KAREA2(5,IELH7))=DVH3
       DV2(KAREA2(2,IELH6))=DVH3
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH3))=DV2(KAREA2(2,IELH3))+
     *                      (B22+S4+D22+
     *                      F23+G23)
       DV2(KAREA2(5,IELH4))=DV2(KAREA2(5,IELH4))+
     *                      (B22+S4+D23+
     *                      F22+G23)
       DV2(KAREA2(5,IELH8))=DV2(KAREA2(5,IELH8))+
     *                      (B23+S4+D23+
     *                      F22+G22)
       DV2(KAREA2(2,IELH7))=DV2(KAREA2(2,IELH7))+
     *                      (B23+S4+D22+
     *                      F23+G22)
      ELSE
       DV2(KAREA2(2,IELH3))=DVH4
       DV2(KAREA2(5,IELH4))=DVH4
       DV2(KAREA2(5,IELH8))=DVH4
       DV2(KAREA2(2,IELH7))=DVH4
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DV2(KAREA2(2,IELH4))=DV2(KAREA2(2,IELH4))+
     *                      (B22+C23+S5+
     *                      E22+G23)
       DV2(KAREA2(5,IELH1))=DV2(KAREA2(5,IELH1))+
     *                      (B22+C22+S5+
     *                      E23+G23)
       DV2(KAREA2(5,IELH5))=DV2(KAREA2(5,IELH5))+
     *                      (B23+C22+S5+
     *                      E23+G22)
       DV2(KAREA2(2,IELH8))=DV2(KAREA2(2,IELH8))+
     *                      (B23+C23+S5+
     *                      E22+G22)
      ELSE
       DV2(KAREA2(2,IELH4))=DVH5
       DV2(KAREA2(5,IELH1))=DVH5
       DV2(KAREA2(5,IELH5))=DVH5
       DV2(KAREA2(2,IELH8))=DVH5
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DV2(KAREA2(1,IELH5))=DV2(KAREA2(1,IELH5))+
     *                      (S6+C22+D23+
     *                      E23+F22)
       DV2(KAREA2(1,IELH6))=DV2(KAREA2(1,IELH6))+
     *                      (S6+C22+D22+
     *                      E23+F23)
       DV2(KAREA2(1,IELH7))=DV2(KAREA2(1,IELH7))+
     *                      (S6+C23+D22+
     *                      E22+F23)
       DV2(KAREA2(1,IELH8))=DV2(KAREA2(1,IELH8))+
     *                      (S6+C23+D23+
     *                      E22+F22)
      ELSE
       DV2(KAREA2(1,IELH5))=DVH6
       DV2(KAREA2(1,IELH6))=DVH6
       DV2(KAREA2(1,IELH7))=DVH6
       DV2(KAREA2(1,IELH8))=DVH6
      ENDIF
C
C
      DV2(KAREA2(3,IELH1))=S7+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH2))=S7+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH3))=S7+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH4))=S7+C26+D27+
     *                     E26+F25
C
      DV2(KAREA2(6,IELH1))=S8+C25+D27+
     *                     E27+F25
      DV2(KAREA2(6,IELH2))=S8+C25+D25+
     *                     E27+F27
      DV2(KAREA2(6,IELH3))=S8+C27+D25+
     *                     E25+F27
      DV2(KAREA2(6,IELH4))=S8+C27+D27+
     *                     E25+F25
C
      DV2(KAREA2(3,IELH5))=S9+C25+D26+
     *                     E27+F26
      DV2(KAREA2(3,IELH6))=S9+C26+D25+
     *                     E26+F27
      DV2(KAREA2(3,IELH7))=S9+C27+D26+
     *                     E25+F26
      DV2(KAREA2(3,IELH8))=S9+C26+D27+
     *                     E26+F25
C
      B31=A1*DWH1
      B32=A2*DWH1
      B33=A3*DWH1
      B34=A4*DWH1
      B35=A5*DWH1
      B36=A6*DWH1
      B37=A7*DWH1
C
      C31=A1*DWH2
      C32=A2*DWH2
      C33=A3*DWH2
      C34=A4*DWH2
      C35=A5*DWH2
      C36=A6*DWH2
      C37=A7*DWH2
C
      D31=A1*DWH3
      D32=A2*DWH3
      D33=A3*DWH3
      D34=A4*DWH3
      D35=A5*DWH3
      D36=A6*DWH3
      D37=A7*DWH3
C
      E31=A1*DWH4
      E32=A2*DWH4
      E33=A3*DWH4
      E34=A4*DWH4
      E35=A5*DWH4
      E36=A6*DWH4
      E37=A7*DWH4
C
      F31=A1*DWH5
      F32=A2*DWH5
      F33=A3*DWH5
      F34=A4*DWH5
      F35=A5*DWH5
      F36=A6*DWH5
      F37=A7*DWH5
C
      G31=A1*DWH6
      G32=A2*DWH6
      G33=A3*DWH6
      G34=A4*DWH6
      G35=A5*DWH6
      G36=A6*DWH6
      G37=A7*DWH6
C
C
      T1=B31+G34
      T2=C31+E34
      T3=D31+F34
      T4=C34+E31
      T5=D34+F31
      T6=B34+G31
      T7=B35+G37
      T8=B36+G36
      T9=B37+G35
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH1))=DW2(KAREA2(1,IELH1))+
     *                      (T1+C32+D33+
     *                      E33+F32)
       DW2(KAREA2(1,IELH2))=DW2(KAREA2(1,IELH2))+
     *                      (T1+C32+D32+
     *                      E33+F33)
       DW2(KAREA2(1,IELH3))=DW2(KAREA2(1,IELH3))+
     *                      (T1+C33+D32+
     *                      E32+F33)
       DW2(KAREA2(1,IELH4))=DW2(KAREA2(1,IELH4))+
     *                      (T1+C33+D33+
     *                      E32+F32)
      ELSE
       DW2(KAREA2(1,IELH1))=DWH1
       DW2(KAREA2(1,IELH2))=DWH1
       DW2(KAREA2(1,IELH3))=DWH1
       DW2(KAREA2(1,IELH4))=DWH1
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH1))=DW2(KAREA2(2,IELH1))+
     *                      (B32+T2+D33+
     *                      F32+G33)
       DW2(KAREA2(5,IELH2))=DW2(KAREA2(5,IELH2))+
     *                      (B32+T2+D32+
     *                      F33+G33)
       DW2(KAREA2(5,IELH6))=DW2(KAREA2(5,IELH6))+
     *                      (B33+T2+D32+
     *                      F33+G32)
       DW2(KAREA2(2,IELH5))=DW2(KAREA2(2,IELH5))+
     *                      (B33+T2+D33+
     *                      F32+G32)
      ELSE
       DW2(KAREA2(2,IELH1))=DWH2
       DW2(KAREA2(5,IELH2))=DWH2
       DW2(KAREA2(5,IELH6))=DWH2
       DW2(KAREA2(2,IELH5))=DWH2
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH2))=DW2(KAREA2(2,IELH2))+
     *                      (B32+C32+T3+
     *                      E33+G33)
       DW2(KAREA2(5,IELH3))=DW2(KAREA2(5,IELH3))+
     *                      (B32+C33+T3+
     *                      E32+G33)
       DW2(KAREA2(5,IELH7))=DW2(KAREA2(5,IELH7))+
     *                      (B33+C33+T3+
     *                      E32+G32)
       DW2(KAREA2(2,IELH6))=DW2(KAREA2(2,IELH6))+
     *                      (B33+C32+T3+
     *                      E33+G32)
      ELSE
       DW2(KAREA2(2,IELH2))=DWH3
       DW2(KAREA2(5,IELH3))=DWH3
       DW2(KAREA2(5,IELH7))=DWH3
       DW2(KAREA2(2,IELH6))=DWH3
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH3))=DW2(KAREA2(2,IELH3))+
     *                      (B32+T4+D32+
     *                      F33+G33)
       DW2(KAREA2(5,IELH4))=DW2(KAREA2(5,IELH4))+
     *                      (B32+T4+D33+
     *                      F32+G33)
       DW2(KAREA2(5,IELH8))=DW2(KAREA2(5,IELH8))+
     *                      (B33+T4+D33+
     *                      F32+G32)
       DW2(KAREA2(2,IELH7))=DW2(KAREA2(2,IELH7))+
     *                      (B33+T4+D32+
     *                      F33+G32)
      ELSE
       DW2(KAREA2(2,IELH3))=DWH4
       DW2(KAREA2(5,IELH4))=DWH4
       DW2(KAREA2(5,IELH8))=DWH4
       DW2(KAREA2(2,IELH7))=DWH4
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DW2(KAREA2(2,IELH4))=DW2(KAREA2(2,IELH4))+
     *                      (B32+C33+T5+
     *                      E32+G33)
       DW2(KAREA2(5,IELH1))=DW2(KAREA2(5,IELH1))+
     *                      (B32+C32+T5+
     *                      E33+G33)
       DW2(KAREA2(5,IELH5))=DW2(KAREA2(5,IELH5))+
     *                      (B33+C32+T5+
     *                      E33+G32)
       DW2(KAREA2(2,IELH8))=DW2(KAREA2(2,IELH8))+
     *                      (B33+C33+T5+
     *                      E32+G32)
      ELSE
       DW2(KAREA2(2,IELH4))=DWH5
       DW2(KAREA2(5,IELH1))=DWH5
       DW2(KAREA2(5,IELH5))=DWH5
       DW2(KAREA2(2,IELH8))=DWH5
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DW2(KAREA2(1,IELH5))=DW2(KAREA2(1,IELH5))+
     *                      (T6+C32+D33+
     *                      E33+F32)
       DW2(KAREA2(1,IELH6))=DW2(KAREA2(1,IELH6))+
     *                      (T6+C32+D32+
     *                      E33+F33)
       DW2(KAREA2(1,IELH7))=DW2(KAREA2(1,IELH7))+
     *                      (T6+C33+D32+
     *                      E32+F33)
       DW2(KAREA2(1,IELH8))=DW2(KAREA2(1,IELH8))+
     *                      (T6+C33+D33+
     *                      E32+F32)
      ELSE
       DW2(KAREA2(1,IELH5))=DWH6
       DW2(KAREA2(1,IELH6))=DWH6
       DW2(KAREA2(1,IELH7))=DWH6
       DW2(KAREA2(1,IELH8))=DWH6
      ENDIF
C
C
      DW2(KAREA2(3,IELH1))=T7+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH2))=T7+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH3))=T7+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH4))=T7+C36+D37+
     *                     E36+F35
C
      DW2(KAREA2(6,IELH1))=T8+C35+D37+
     *                     E37+F35
      DW2(KAREA2(6,IELH2))=T8+C35+D35+
     *                     E37+F37
      DW2(KAREA2(6,IELH3))=T8+C37+D35+
     *                     E35+F37
      DW2(KAREA2(6,IELH4))=T8+C37+D37+
     *                     E35+F35
C
      DW2(KAREA2(3,IELH5))=T9+C35+D36+
     *                     E37+F36
      DW2(KAREA2(3,IELH6))=T9+C36+D35+
     *                     E36+F37
      DW2(KAREA2(3,IELH7))=T9+C37+D36+
     *                     E35+F36
      DW2(KAREA2(3,IELH8))=T9+C36+D37+
     *                     E36+F35
C
C
C
10    CONTINUE
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRU  (DU1,DV1,DW1,DU2,DV2,DW2,
     *                     KVERT1,KAREA1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KAREA2,KADJ2,NEQ2,NEL2,NVT2,
     *                     AVOL,IVEL)
************************************************************************
*    Purpose:  Restricts the  fine grid solution vector (DU2,DV2,DW2) to
*              the coarse grid solution vector (DU1,DV1,DW1)
*
*    Remark:   i.) RESTRU works for all cases of boundary conditions
*              ii.) weighted interpolation included
*                
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DW2               - fine grid solution vector
*      KVERT1,KAREA1,..,NVT1  - data of the coarse grid
*      KVERT2,KAREA2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1,DV1               - coarse grid solution vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*4 AVOL
C
      PARAMETER (NNVE=8,NNAE=6)
      PARAMETER (A1=0.125D0,A2=0.25D0,A3=0.08333333D0)
      PARAMETER (R1=0.25D0,R2=0.5D0,R3=0.16666666D0)
C
      DIMENSION DU1(*),DV1(*),DW1(*),  DU2(*),DV2(*),DW2(*),
     *          KVERT1(NNVE,*),KAREA1(NNAE,*),KADJ1(NNAE,*),
     *          KVERT2(NNVE,*),KAREA2(NNAE,*),KADJ2(NNAE,*),
     *          AVOL(*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C
C
      I1=KAREA2(1,IELH1)
      I2=KAREA2(2,IELH1)
      I3=KAREA2(3,IELH1)
      I4=KAREA2(4,IELH1)
      I5=KAREA2(5,IELH1)
      I6=KAREA2(6,IELH1)
      I7=KAREA2(1,IELH2)
      I8=KAREA2(2,IELH2)
      I9=KAREA2(3,IELH2)
      I10=KAREA2(5,IELH2)
      I11=KAREA2(6,IELH2)
      I12=KAREA2(1,IELH3)
      I13=KAREA2(2,IELH3)
      I14=KAREA2(3,IELH3)
      I15=KAREA2(5,IELH3)
      I16=KAREA2(6,IELH3)
      I17=KAREA2(1,IELH4)
      I18=KAREA2(2,IELH4)
      I19=KAREA2(5,IELH4)
      I20=KAREA2(6,IELH4)
      I21=KAREA2(1,IELH5)
      I22=KAREA2(2,IELH5)
      I23=KAREA2(3,IELH5)
      I24=KAREA2(4,IELH5)
      I25=KAREA2(5,IELH5)
      I26=KAREA2(1,IELH6)
      I27=KAREA2(2,IELH6)
      I28=KAREA2(3,IELH6)
      I29=KAREA2(5,IELH6)
      I30=KAREA2(1,IELH7)
      I31=KAREA2(2,IELH7)
      I32=KAREA2(3,IELH7)
      I33=KAREA2(5,IELH7)
      I34=KAREA2(1,IELH8)
      I35=KAREA2(2,IELH8)
      I36=KAREA2(5,IELH8)
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
      DUH13=DU2(I13)
      DUH14=DU2(I14)
      DUH15=DU2(I15)
      DUH16=DU2(I16)
      DUH17=DU2(I17)
      DUH18=DU2(I18)
      DUH19=DU2(I19)
      DUH20=DU2(I20)
      DUH21=DU2(I21)
      DUH22=DU2(I22)
      DUH23=DU2(I23)
      DUH24=DU2(I24)
      DUH25=DU2(I25)
      DUH26=DU2(I26)
      DUH27=DU2(I27)
      DUH28=DU2(I28)
      DUH29=DU2(I29)
      DUH30=DU2(I30)
      DUH31=DU2(I31)
      DUH32=DU2(I32)
      DUH33=DU2(I33)
      DUH34=DU2(I34)
      DUH35=DU2(I35)
      DUH36=DU2(I36)
C
C
      DVH1= DV2(I1)
      DVH2= DV2(I2)
      DVH3= DV2(I3)
      DVH4= DV2(I4)
      DVH5= DV2(I5)
      DVH6= DV2(I6)
      DVH7= DV2(I7)
      DVH8= DV2(I8)
      DVH9= DV2(I9)
      DVH10=DV2(I10)
      DVH11=DV2(I11)
      DVH12=DV2(I12)
      DVH13=DV2(I13)
      DVH14=DV2(I14)
      DVH15=DV2(I15)
      DVH16=DV2(I16)
      DVH17=DV2(I17)
      DVH18=DV2(I18)
      DVH19=DV2(I19)
      DVH20=DV2(I20)
      DVH21=DV2(I21)
      DVH22=DV2(I22)
      DVH23=DV2(I23)
      DVH24=DV2(I24)
      DVH25=DV2(I25)
      DVH26=DV2(I26)
      DVH27=DV2(I27)
      DVH28=DV2(I28)
      DVH29=DV2(I29)
      DVH30=DV2(I30)
      DVH31=DV2(I31)
      DVH32=DV2(I32)
      DVH33=DV2(I33)
      DVH34=DV2(I34)
      DVH35=DV2(I35)
      DVH36=DV2(I36)
C
      DWH1= DW2(I1)
      DWH2= DW2(I2)
      DWH3= DW2(I3)
      DWH4= DW2(I4)
      DWH5= DW2(I5)
      DWH6= DW2(I6)
      DWH7= DW2(I7)
      DWH8= DW2(I8)
      DWH9= DW2(I9)
      DWH10=DW2(I10)
      DWH11=DW2(I11)
      DWH12=DW2(I12)
      DWH13=DW2(I13)
      DWH14=DW2(I14)
      DWH15=DW2(I15)
      DWH16=DW2(I16)
      DWH17=DW2(I17)
      DWH18=DW2(I18)
      DWH19=DW2(I19)
      DWH20=DW2(I20)
      DWH21=DW2(I21)
      DWH22=DW2(I22)
      DWH23=DW2(I23)
      DWH24=DW2(I24)
      DWH25=DW2(I25)
      DWH26=DW2(I26)
      DWH27=DW2(I27)
      DWH28=DW2(I28)
      DWH29=DW2(I29)
      DWH30=DW2(I30)
      DWH31=DW2(I31)
      DWH32=DW2(I32)
      DWH33=DW2(I33)
      DWH34=DW2(I34)
      DWH35=DW2(I35)
      DWH36=DW2(I36)
C
      IF (IVEL.EQ.1) THEN
      AVOL11=AVOL(IEL1)
      ENDIF
C
      IF (IVEL.EQ.0) THEN
      AVOL11=1D0
      ENDIF
C
C *** The area IA1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IA1)=(A1*(DUH1+DUH7+DUH12+DUH17)+
     *    A2*(DUH3+DUH14+DUH4+DUH9)-
     *    A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11
        DV1(IA1)=(A1*(DVH1+DVH7+DVH12+DVH17)+
     *    A2*(DVH3+DVH14+DVH4+DVH9)-
     *    A3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11
        DW1(IA1)=(A1*(DWH1+DWH7+DWH12+DWH17)+
     *    A2*(DWH3+DWH14+DWH4+DWH9)-
     *    A3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11
      ELSE
        DU1(IA1)=DU1(IA1)+
     *    (A1*(DUH1+DUH7+DUH12+DUH17)+A2*(DUH3+DUH14+DUH4+DUH9)-
     *    A3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11
        DV1(IA1)=DV1(IA1)+
     *    (A1*(DVH1+DVH7+DVH12+DVH17)+A2*(DVH3+DVH14+DVH4+DVH9)-
     *    A3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11
        DW1(IA1)=DW1(IA1)+
     *    (A1*(DWH1+DWH7+DWH12+DWH17)+A2*(DWH3+DWH14+DWH4+DWH9)-
     *    A3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11
      ENDIF
      ELSE
C     case of a boundary edge
        DU1(IA1)=R1*(DUH1+DUH7+DUH12+DUH17)+
     *     R2*(DUH3+DUH14+DUH4+DUH9)-
     *    R3*(DUH2+DUH10+DUH8+DUH15+DUH13+DUH19+DUH18+
     *        DUH5+DUH6+DUH11+
     *        DUH16+DUH20)
        DV1(IA1)=R1*(DVH1+DVH7+DVH12+DVH17)+
     *     R2*(DVH3+DVH14+DVH4+DVH9)-
     *    R3*(DVH2+DVH10+DVH8+DVH15+DVH13+DVH19+DVH18+
     *        DVH5+DVH6+DVH11+
     *        DVH16+DVH20)
        DW1(IA1)=R1*(DWH1+DWH7+DWH12+DWH17)+
     *    R2*(DWH3+DWH14+DWH4+DWH9)-
     *    R3*(DWH2+DWH10+DWH8+DWH15+DWH13+DWH19+DWH18+
     *        DWH5+DWH6+DWH11+
     *        DWH16+DWH20)
      ENDIF
C
C *** The area IA2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IA2)=(A1*(DUH2+DUH10+DUH29+DUH22)+
     *     A2*(DUH3+DUH23+DUH6+DUH11)-
     *     A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28))*AVOL11  
        DV1(IA2)=(A1*(DVH2+DVH10+DVH29+DVH22)+
     *     A2*(DVH3+DVH23+DVH6+DVH11)-
     *     A3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28))*AVOL11  
        DW1(IA2)=(A1*(DWH2+DWH10+DWH29+DWH22)+
     *     A2*(DWH3+DWH23+DWH6+DWH11)-
     *     A3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28))*AVOL11  
      ELSE
        DU1(IA2)=DU1(IA2)+
     *     (A1*(DUH2+DUH10+DUH29+DUH22)+A2*(DUH3+DUH23+DUH6+DUH11)-
     *     A3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28))*AVOL11  
        DV1(IA2)=DV1(IA2)+
     *     (A1*(DVH2+DVH10+DVH29+DVH22)+A2*(DVH3+DVH23+DVH6+DVH11)-
     *     A3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28))*AVOL11  
        DW1(IA2)=DW1(IA2)+
     *     (A1*(DWH2+DWH10+DWH29+DWH22)+A2*(DWH3+DWH23+DWH6+DWH11)-
     *     A3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28))*AVOL11
       ENDIF
       ELSE  
C     case of a boundary area
        DU1(IA2)=R1*(DUH2+DUH10+DUH29+DUH22)+
     *     R2*(DUH3+DUH23+DUH6+DUH11)-
     *     R3*(DUH1+DUH7+DUH8+DUH27+DUH26+DUH21+DUH25+
     *         DUH5+DUH4+DUH9+
     *         DUH24+DUH28)  
        DV1(IA2)=R1*(DVH2+DVH10+DVH29+DVH22)+
     *     R2*(DVH3+DVH23+DVH6+DVH11)-
     *     R3*(DVH1+DVH7+DVH8+DVH27+DVH26+DVH21+DVH25+
     *         DVH5+DVH4+DVH9+
     *         DVH24+DVH28)  
        DW1(IA2)=R1*(DWH2+DWH10+DWH29+DWH22)+
     *     R2*(DWH3+DWH23+DWH6+DWH11)-
     *     R3*(DWH1+DWH7+DWH8+DWH27+DWH26+DWH21+DWH25+
     *         DWH5+DWH4+DWH9+
     *         DWH24+DWH28)  
      ENDIF
C
C *** The edge IA3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IA3)=(A1*(DUH8+DUH15+DUH33+DUH27)+
     *     A2*(DUH11+DUH16+DUH9+DUH28)-
     *     A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10))*AVOL11  
        DV1(IA3)=(A1*(DVH8+DVH15+DVH33+DVH27)+
     *     A2*(DVH11+DVH16+DVH9+DVH28)-
     *     A3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10))*AVOL11  
        DW1(IA3)=(A1*(DWH8+DWH15+DWH33+DWH27)+
     *     A2*(DWH11+DWH16+DWH9+DWH28)-
     *     A3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10))*AVOL11  
      ELSE
        DU1(IA3)= DU1(IA3)+
     *     (A1*(DUH8+DUH15+DUH33+DUH27)+A2*(DUH11+DUH16+DUH9+DUH28)-
     *     A3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10))*AVOL11  
        DV1(IA3)=DV1(IA3)+
     *     (A1*(DVH8+DVH15+DVH33+DVH27)+A2*(DVH11+DVH16+DVH9+DVH28)-
     *     A3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10))*AVOL11  
        DW1(IA3)=DW1(IA3)+
     *     (A1*(DWH8+DWH15+DWH33+DWH27)+A2*(DWH11+DWH16+DWH9+DWH28)-
     *     A3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10))*AVOL11
        ENDIF
        ELSE  
C     case of a boundary area 
        DU1(IA3)=R1*(DUH8+DUH15+DUH33+DUH27)+
     *     R2*(DUH11+DUH16+DUH9+DUH28)-
     *     R3*(DUH3+DUH14+DUH32+DUH23+DUH7+DUH12+DUH13+
     *         DUH31+DUH30+DUH26+
     *          DUH29+DUH10)  
        DV1(IA3)=R1*(DVH8+DVH15+DVH33+DVH27)+
     *     R2*(DVH11+DVH16+DVH9+DVH28)-
     *     R3*(DVH3+DVH14+DVH32+DVH23+DVH7+DVH12+DVH13+
     *         DVH31+DVH30+DVH26+
     *          DVH29+DVH10)  
        DW1(IA3)=R1*(DWH8+DWH15+DWH33+DWH27)+
     *     R2*(DWH11+DWH16+DWH9+DWH28)-
     *     R3*(DWH3+DWH14+DWH32+DWH23+DWH7+DWH12+DWH13+
     *         DWH31+DWH30+DWH26+
     *          DWH29+DWH10)  
      ENDIF
C
C *** The area IA4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IA4)=(A1*(DUH13+DUH31+DUH36+DUH19)+
     *   A2*(DUH32+DUH14+DUH16+DUH20)-
     *   A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24))*AVOL11  
        DV1(IA4)=(A1*(DVH13+DVH31+DVH36+DVH19)+
     *   A2*(DVH32+DVH14+DVH16+DVH20)-
     *   A3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24))*AVOL11  
        DW1(IA4)=(A1*(DWH13+DWH31+DWH36+DWH19)+
     *   A2*(DWH32+DWH14+DWH16+DWH20)-
     *   A3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24))*AVOL11  
      ELSE
        DU1(IA4)=DU1(IA4)+
     *   (A1*(DUH13+DUH31+DUH36+DUH19)+A2*(DUH32+DUH14+DUH16+DUH20)-
     *   A3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24))*AVOL11  
        DV1(IA4)= DV1(IA4)+
     *   (A1*(DVH13+DVH31+DVH36+DVH19)+A2*(DVH32+DVH14+DVH16+DVH20)-
     *   A3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24))*AVOL11  
        DW1(IA4)=DW1(IA4)+
     *   (A1*(DWH13+DWH31+DWH36+DWH19)+A2*(DWH32+DWH14+DWH16+DWH20)-
     *   A3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24))*AVOL11
        ENDIF
        ELSE  
C     case of a boundary area
        DU1(IA4)=R1*(DUH13+DUH31+DUH36+DUH19)+
     *   R2*(DUH32+DUH14+DUH16+DUH20)-
     *   R3*(DUH17+DUH12+DUH15+DUH33+DUH30+DUH34+DUH35+
     *       DUH18+DUH4+DUH9+
     *       DUH28+DUH24)  
        DV1(IA4)=R1*(DVH13+DVH31+DVH36+DVH19)+
     *   R2*(DVH32+DVH14+DVH16+DVH20)-
     *   R3*(DVH17+DVH12+DVH15+DVH33+DVH30+DVH34+DVH35+
     *       DVH18+DVH4+DVH9+
     *       DVH28+DVH24)  
        DW1(IA4)=R1*(DWH13+DWH31+DWH36+DWH19)+
     *   R2*(DWH32+DWH14+DWH16+DWH20)-
     *   R3*(DWH17+DWH12+DWH15+DWH33+DWH30+DWH34+DWH35+
     *       DWH18+DWH4+DWH9+
     *       DWH28+DWH24)  
      ENDIF
C
C *** The area IA5
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(5,IEL1).GT.IEL1) THEN
        DU1(IA5)=(A1*(DUH5+DUH18+DUH35+DUH25)+
     *    A2*(DUH6+DUH20+DUH4+DUH24)-
     *    A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23))*AVOL11  
        DV1(IA5)=(A1*(DVH5+DVH18+DVH35+DVH25)+
     *    A2*(DVH6+DVH20+DVH4+DVH24)-
     *    A3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23))*AVOL11  
        DW1(IA5)=(A1*(DWH5+DWH18+DWH35+DWH25)+
     *    A2*(DWH6+DWH20+DWH4+DWH24)-
     *    A3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23))*AVOL11  
      ELSE
        DU1(IA5)=DU1(IA5)+
     *    (A1*(DUH5+DUH18+DUH35+DUH25)+A2*(DUH6+DUH20+DUH4+DUH24)-
     *    A3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23))*AVOL11  
        DV1(IA5)=DV1(IA5)+
     *    (A1*(DVH5+DVH18+DVH35+DVH25)+A2*(DVH6+DVH20+DVH4+DVH24)-
     *    A3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23))*AVOL11  
        DW1(IA5)= DW1(IA5)+
     *    (A1*(DWH5+DWH18+DWH35+DWH25)+A2*(DWH6+DWH20+DWH4+DWH24)-
     *    A3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23))*AVOL11
       ENDIF
       ELSE  
C     case of a boundary area
        DU1(IA5)=R1*(DUH5+DUH18+DUH35+DUH25)+
     *    R2*(DUH6+DUH20+DUH4+DUH24)-
     *    R3*(DUH1+DUH17+DUH19+DUH36+DUH34+DUH21+DUH22+
     *        DUH2+DUH3+DUH14+
     *        DUH32+DUH23)  
        DV1(IA5)=R1*(DVH5+DVH18+DVH35+DVH25)+
     *    R2*(DVH6+DVH20+DVH4+DVH24)-
     *    R3*(DVH1+DVH17+DVH19+DVH36+DVH34+DVH21+DVH22+
     *        DVH2+DVH3+DVH14+
     *        DVH32+DVH23)  
        DW1(IA5)=R1*(DWH5+DWH18+DWH35+DWH25)+
     *    R2*(DWH6+DWH20+DWH4+DWH24)-
     *    R3*(DWH1+DWH17+DWH19+DWH36+DWH34+DWH21+DWH22+
     *        DWH2+DWH3+DWH14+
     *        DWH32+DWH23)  
      ENDIF
C
C *** The area IA6
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(6,IEL1).GT.IEL1) THEN
        DU1(IA6)=(A1*(DUH26+DUH30+DUH34+DUH21)+
     *    A2*(DUH23+DUH32+DUH24+DUH28)-
     *    A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11  
        DV1(IA6)=(A1*(DVH26+DVH30+DVH34+DVH21)+
     *    A2*(DVH23+DVH32+DVH24+DVH28)-
     *    A3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11  
        DW1(IA6)=(A1*(DWH26+DWH30+DWH34+DWH21)+
     *    A2*(DWH23+DWH32+DWH24+DWH28)-
     *    A3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11  
      ELSE
        DU1(IA6)= DU1(IA6)+
     *    (A1*(DUH26+DUH30+DUH34+DUH21)+A2*(DUH23+DUH32+DUH24+DUH28)-
     *    A3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20))*AVOL11  
        DV1(IA6)=DV1(IA6)+
     *    (A1*(DVH26+DVH30+DVH34+DVH21)+A2*(DVH23+DVH32+DVH24+DVH28)-
     *    A3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20))*AVOL11  
        DW1(IA6)=DW1(IA6)+
     *    (A1*(DWH26+DWH30+DWH34+DWH21)+A2*(DWH23+DWH32+DWH24+DWH28)-
     *    A3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20))*AVOL11 
        ENDIF
        ELSE 
C     case of a boundary area
        DU1(IA6)=R1*(DUH26+DUH30+DUH34+DUH21)+
     *    R2*(DUH23+DUH32+DUH24+DUH28)-
     *    R3*(DUH22+DUH29+DUH27+DUH33+DUH31+DUH36+DUH21+
     *        DUH25+DUH6+DUH11+
     *        DUH16+DUH20)  
        DV1(IA6)=R1*(DVH26+DVH30+DVH34+DVH21)+
     *    R2*(DVH23+DVH32+DVH24+DVH28)-
     *    R3*(DVH22+DVH29+DVH27+DVH33+DVH31+DVH36+DVH21+
     *        DVH25+DVH6+DVH11+
     *        DVH16+DVH20)  
        DW1(IA6)=R1*(DWH26+DWH30+DWH34+DWH21)+
     *    R2*(DWH23+DWH32+DWH24+DWH28)-
     *    R3*(DWH22+DWH29+DWH27+DWH33+DWH31+DWH36+DWH21+
     *        DWH25+DWH6+DWH11+
     *        DWH16+DWH20)  
      ENDIF
C
10    CONTINUE
C
      IF (IVEL.EQ.1) THEN
      DO 555 IEL1=1,NEL1
      DO 666 I=1,6
      IADJ=KADJ1(I,IEL1)
      IA=KAREA1(I,IEL1)
      IF (IADJ.EQ.0) THEN
      GOTO 666
      ENDIF
      IF (IADJ.GT.IEL1) THEN
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DU1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      DV1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DV1(IA)
      DW1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DW1(IA)
      ELSE
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DU1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DU1(IA)
      DV1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DV1(IA)
      DW1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DW1(IA)
      ENDIF
666   CONTINUE
555   CONTINUE
      ENDIF
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRD  (DNU1,DNU2,
     *                     KVERT1,KAREA1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KAREA2,KADJ2,NEQ2,NEL2,NVT2,
     *                     AVOL,IVEL)
************************************************************************
*    Purpose:  Restricts the  fine grid solution vector DNU2 to
*              the coarse grid solution vector DNU1
*
*    Remark:   i.) RESTRU works for all cases of boundary conditions
*              ii.) weighted interpolation included
*                
*-----------------------------------------------------------------------
*    Input:
*      DNU2               - fine grid solution vector
*      KVERT1,KAREA1,..,NVT1  - data of the coarse grid
*      KVERT2,KAREA2,..,NVT2  - data of the fine grid
*
*    Output:
*      DNU1               - coarse grid solution vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      REAL*4 AVOL
C
      PARAMETER (NNVE=8,NNAE=6)
      PARAMETER (A1=0.125D0,A2=0.25D0,A3=0.08333333D0)
      PARAMETER (R1=0.25D0,R2=0.5D0,R3=0.16666666D0)
C
      DIMENSION KVERT1(NNVE,*),KAREA1(NNAE,*),KADJ1(NNAE,*),
     *          KVERT2(NNVE,*),KAREA2(NNAE,*),KADJ2(NNAE,*),
     *          AVOL(*),DNU1(*),DNU2(*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C
C
      I1=KAREA2(1,IELH1)
      I2=KAREA2(2,IELH1)
      I3=KAREA2(3,IELH1)
      I4=KAREA2(4,IELH1)
      I5=KAREA2(5,IELH1)
      I6=KAREA2(6,IELH1)
      I7=KAREA2(1,IELH2)
      I8=KAREA2(2,IELH2)
      I9=KAREA2(3,IELH2)
      I10=KAREA2(5,IELH2)
      I11=KAREA2(6,IELH2)
      I12=KAREA2(1,IELH3)
      I13=KAREA2(2,IELH3)
      I14=KAREA2(3,IELH3)
      I15=KAREA2(5,IELH3)
      I16=KAREA2(6,IELH3)
      I17=KAREA2(1,IELH4)
      I18=KAREA2(2,IELH4)
      I19=KAREA2(5,IELH4)
      I20=KAREA2(6,IELH4)
      I21=KAREA2(1,IELH5)
      I22=KAREA2(2,IELH5)
      I23=KAREA2(3,IELH5)
      I24=KAREA2(4,IELH5)
      I25=KAREA2(5,IELH5)
      I26=KAREA2(1,IELH6)
      I27=KAREA2(2,IELH6)
      I28=KAREA2(3,IELH6)
      I29=KAREA2(5,IELH6)
      I30=KAREA2(1,IELH7)
      I31=KAREA2(2,IELH7)
      I32=KAREA2(3,IELH7)
      I33=KAREA2(5,IELH7)
      I34=KAREA2(1,IELH8)
      I35=KAREA2(2,IELH8)
      I36=KAREA2(5,IELH8)
C
      DNUH1= DNU2(I1)
      DNUH2= DNU2(I2)
      DNUH3= DNU2(I3)
      DNUH4= DNU2(I4)
      DNUH5= DNU2(I5)
      DNUH6= DNU2(I6)
      DNUH7= DNU2(I7)
      DNUH8= DNU2(I8)
      DNUH9= DNU2(I9)
      DNUH10=DNU2(I10)
      DNUH11=DNU2(I11)
      DNUH12=DNU2(I12)
      DNUH13=DNU2(I13)
      DNUH14=DNU2(I14)
      DNUH15=DNU2(I15)
      DNUH16=DNU2(I16)
      DNUH17=DNU2(I17)
      DNUH18=DNU2(I18)
      DNUH19=DNU2(I19)
      DNUH20=DNU2(I20)
      DNUH21=DNU2(I21)
      DNUH22=DNU2(I22)
      DNUH23=DNU2(I23)
      DNUH24=DNU2(I24)
      DNUH25=DNU2(I25)
      DNUH26=DNU2(I26)
      DNUH27=DNU2(I27)
      DNUH28=DNU2(I28)
      DNUH29=DNU2(I29)
      DNUH30=DNU2(I30)
      DNUH31=DNU2(I31)
      DNUH32=DNU2(I32)
      DNUH33=DNU2(I33)
      DNUH34=DNU2(I34)
      DNUH35=DNU2(I35)
      DNUH36=DNU2(I36)
C
      IF (IVEL.EQ.1) THEN
      AVOL11=AVOL(IEL1)
      ENDIF
C
      IF (IVEL.EQ.0) THEN
      AVOL11=1D0
      ENDIF
C
C *** The area IA1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DNU1(IA1)=(A1*(DNUH1+DNUH7+DNUH12+DNUH17)+
     *    A2*(DNUH3+DNUH14+DNUH4+DNUH9)-
     *    A3*(DNUH2+DNUH10+DNUH8+DNUH15+DNUH13+DNUH19+DNUH18+
     *        DNUH5+DNUH6+DNUH11+
     *        DNUH16+DNUH20))*AVOL11  
      ELSE
        DNU1(IA1)=DNU1(IA1)+
     *    (A1*(DNUH1+DNUH7+DNUH12+DNUH17)+A2*(DNUH3+DNUH14+DNUH4+DNUH9)-
     *    A3*(DNUH2+DNUH10+DNUH8+DNUH15+DNUH13+DNUH19+DNUH18+
     *        DNUH5+DNUH6+DNUH11+
     *        DNUH16+DNUH20))*AVOL11
      ENDIF
      ELSE
C     case of a boundary edge
        DNU1(IA1)=R1*(DNUH1+DNUH7+DNUH12+DNUH17)+
     *     R2*(DNUH3+DNUH14+DNUH4+DNUH9)-
     *    R3*(DNUH2+DNUH10+DNUH8+DNUH15+DNUH13+DNUH19+DNUH18+
     *        DNUH5+DNUH6+DNUH11+
     *        DNUH16+DNUH20)
      ENDIF
C
C *** The area IA2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DNU1(IA2)=(A1*(DNUH2+DNUH10+DNUH29+DNUH22)+
     *     A2*(DNUH3+DNUH23+DNUH6+DNUH11)-
     *     A3*(DNUH1+DNUH7+DNUH8+DNUH27+DNUH26+DNUH21+DNUH25+
     *         DNUH5+DNUH4+DNUH9+
     *         DNUH24+DNUH28))*AVOL11  
      ELSE
        DNU1(IA2)=DNU1(IA2)+(A1*(DNUH2+DNUH10+DNUH29+DNUH22)+
     *     A2*(DNUH3+DNUH23+DNUH6+DNUH11)-
     *     A3*(DNUH1+DNUH7+DNUH8+DNUH27+DNUH26+DNUH21+DNUH25+
     *         DNUH5+DNUH4+DNUH9+
     *         DNUH24+DNUH28))*AVOL11  
       ENDIF
       ELSE  
C     case of a boundary area
        DNU1(IA2)=R1*(DNUH2+DNUH10+DNUH29+DNUH22)+
     *     R2*(DNUH3+DNUH23+DNUH6+DNUH11)-
     *     R3*(DNUH1+DNUH7+DNUH8+DNUH27+DNUH26+DNUH21+DNUH25+
     *         DNUH5+DNUH4+DNUH9+
     *         DNUH24+DNUH28)  
      ENDIF
C
C *** The edge IA3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner area
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DNU1(IA3)=(A1*(DNUH8+DNUH15+DNUH33+DNUH27)+
     *     A2*(DNUH11+DNUH16+DNUH9+DNUH28)-
     *     A3*(DNUH3+DNUH14+DNUH32+DNUH23+DNUH7+DNUH12+DNUH13+
     *         DNUH31+DNUH30+DNUH26+
     *          DNUH29+DNUH10))*AVOL11  
      ELSE
        DNU1(IA3)=DNU1(IA3)+(A1*(DNUH8+DNUH15+DNUH33+DNUH27)+
     *     A2*(DNUH11+DNUH16+DNUH9+DNUH28)-
     *     A3*(DNUH3+DNUH14+DNUH32+DNUH23+DNUH7+DNUH12+DNUH13+
     *         DNUH31+DNUH30+DNUH26+
     *          DNUH29+DNUH10))*AVOL11  
        ENDIF
        ELSE  
C     case of a boundary area 
        DNU1(IA3)=R1*(DNUH8+DNUH15+DNUH33+DNUH27)+
     *     R2*(DNUH11+DNUH16+DNUH9+DNUH28)-
     *     R3*(DNUH3+DNUH14+DNUH32+DNUH23+DNUH7+DNUH12+DNUH13+
     *         DNUH31+DNUH30+DNUH26+
     *          DNUH29+DNUH10)  
      ENDIF
C
C *** The area IA4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DNU1(IA4)=(A1*(DNUH13+DNUH31+DNUH36+DNUH19)+
     *   A2*(DNUH32+DNUH14+DNUH16+DNUH20)-
     *   A3*(DNUH17+DNUH12+DNUH15+DNUH33+DNUH30+DNUH34+DNUH35+
     *       DNUH18+DNUH4+DNUH9+
     *       DNUH28+DNUH24))*AVOL11  
      ELSE
        DNU1(IA4)= DNU1(IA4)+(A1*(DNUH13+DNUH31+DNUH36+DNUH19)+
     *   A2*(DNUH32+DNUH14+DNUH16+DNUH20)-
     *   A3*(DNUH17+DNUH12+DNUH15+DNUH33+DNUH30+DNUH34+DNUH35+
     *       DNUH18+DNUH4+DNUH9+
     *       DNUH28+DNUH24))*AVOL11  
        ENDIF
        ELSE  
C     case of a boundary area
        DNU1(IA4)=R1*(DNUH13+DNUH31+DNUH36+DNUH19)+
     *   R2*(DNUH32+DNUH14+DNUH16+DNUH20)-
     *   R3*(DNUH17+DNUH12+DNUH15+DNUH33+DNUH30+DNUH34+DNUH35+
     *       DNUH18+DNUH4+DNUH9+
     *       DNUH28+DNUH24) 
      ENDIF
C
C *** The area IA5
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(5,IEL1).GT.IEL1) THEN
        DNU1(IA5)=(A1*(DNUH5+DNUH18+DNUH35+DNUH25)+
     *    A2*(DNUH6+DNUH20+DNUH4+DNUH24)-
     *    A3*(DNUH1+DNUH17+DNUH19+DNUH36+DNUH34+DNUH21+DNUH22+
     *        DNUH2+DNUH3+DNUH14+
     *        DNUH32+DNUH23))*AVOL11  
      ELSE
        DNU1(IA5)=DNU1(IA5)+(A1*(DNUH5+DNUH18+DNUH35+DNUH25)+
     *    A2*(DNUH6+DNUH20+DNUH4+DNUH24)-
     *    A3*(DNUH1+DNUH17+DNUH19+DNUH36+DNUH34+DNUH21+DNUH22+
     *        DNUH2+DNUH3+DNUH14+
     *        DNUH32+DNUH23))*AVOL11  
       ENDIF
       ELSE  
C     case of a boundary area
        DNU1(IA5)=R1*(DNUH5+DNUH18+DNUH35+DNUH25)+
     *    R2*(DNUH6+DNUH20+DNUH4+DNUH24)-
     *    R3*(DNUH1+DNUH17+DNUH19+DNUH36+DNUH34+DNUH21+DNUH22+
     *        DNUH2+DNUH3+DNUH14+
     *        DNUH32+DNUH23) 
      ENDIF
C
C *** The area IA6
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
C     case of an inner area
       IF (KADJ1(6,IEL1).GT.IEL1) THEN
        DNU1(IA6)=(A1*(DNUH26+DNUH30+DNUH34+DNUH21)+
     *    A2*(DNUH23+DNUH32+DNUH24+DNUH28)-
     *    A3*(DNUH22+DNUH29+DNUH27+DNUH33+DNUH31+DNUH36+DNUH21+
     *        DNUH25+DNUH6+DNUH11+
     *        DNUH16+DNUH20))*AVOL11  
      ELSE
        DNU1(IA6)=DNU1(IA6)+(A1*(DNUH26+DNUH30+DNUH34+DNUH21)+
     *    A2*(DNUH23+DNUH32+DNUH24+DNUH28)-
     *    A3*(DNUH22+DNUH29+DNUH27+DNUH33+DNUH31+DNUH36+DNUH21+
     *        DNUH25+DNUH6+DNUH11+
     *        DNUH16+DNUH20))*AVOL11  
        ENDIF
        ELSE 
C     case of a boundary area
        DNU1(IA6)=R1*(DNUH26+DNUH30+DNUH34+DNUH21)+
     *    R2*(DNUH23+DNUH32+DNUH24+DNUH28)-
     *    R3*(DNUH22+DNUH29+DNUH27+DNUH33+DNUH31+DNUH36+DNUH21+
     *        DNUH25+DNUH6+DNUH11+
     *        DNUH16+DNUH20)
      ENDIF
C
10    CONTINUE
C
      IF (IVEL.EQ.1) THEN
      DO 555 IEL1=1,NEL1
      DO 666 I=1,6
      IADJ=KADJ1(I,IEL1)
      IA=KAREA1(I,IEL1)
      IF (IADJ.EQ.0) THEN
      GOTO 666
      ENDIF
      IF (IADJ.GT.IEL1) THEN
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DNU1(IA)=(1D0/SQRT(AVOL1+AVOL2))*DNU1(IA)
      ELSE
      AVOL1=AVOL(IEL1)
      AVOL2=AVOL(IADJ)
      DNU1(IA)=(2D0/SQRT(AVOL1+AVOL2))*DNU1(IA)
      ENDIF
666   CONTINUE
555   CONTINUE
      ENDIF
C
      END
C
C
************************************************************************
      SUBROUTINE MP031(DX1,DX2,KAREA1,KAREA2,KADJ1,KADJ2,NAT1,NAT2,
     *                 NEL1,NEL2)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 11D0/24D0,A2= 7D0/48D0,A3=-5D0/48D0,A4=-1D0/24D0)
      PARAMETER (A5= 11D0/24D0,A6= 1D0/12D0,A7=-1D0/24D0)
      DIMENSION DX1(1),DX2(1),KAREA1(NNAE,1),KAREA2(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DX2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DY1=DX1(KAREA1(1,IEL1))
      DY2=DX1(KAREA1(2,IEL1))
      DY3=DX1(KAREA1(3,IEL1))
      DY4=DX1(KAREA1(4,IEL1))
      DY5=DX1(KAREA1(5,IEL1))
      DY6=DX1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+
     *                      A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A4*DY6
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+
     *                      A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A4*DY6
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+
     *                      A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A4*DY6
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+
     *                      A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A4*DY6
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+
     *                      A2*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+
     *                      A2*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+
     *                      A3*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A2*DY6
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+
     *                      A3*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+
     *                      A2*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A3*DY6
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+
     *                      A2*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A3*DY6
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+
     *                      A3*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A2*DY6
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+
     *                      A3*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+
     *                      A2*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+
     *                      A2*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+
     *                      A3*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A2*DY6
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+
     *                      A3*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+
     *                      A2*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+
     *                      A2*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+
     *                      A3*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A2*DY6
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+
     *                      A3*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+
     *                      A4*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+
     *                      A4*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+
     *                      A4*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+
     *                      A4*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6
      ENDIF
C
C
      DX2(KAREA2(3,IELH1))=A5*DY1+A5*DY2+A6*DY3+A7*DY4+A6*DY5+A7*DY6
      DX2(KAREA2(3,IELH2))=A5*DY1+A6*DY2+A5*DY3+A6*DY4+A7*DY5+A7*DY6
      DX2(KAREA2(3,IELH3))=A5*DY1+A7*DY2+A6*DY3+A5*DY4+A6*DY5+A7*DY6
      DX2(KAREA2(3,IELH4))=A5*DY1+A6*DY2+A7*DY3+A6*DY4+A5*DY5+A7*DY6
C
      DX2(KAREA2(6,IELH1))=A6*DY1+A5*DY2+A7*DY3+A7*DY4+A5*DY5+A6*DY6
      DX2(KAREA2(6,IELH2))=A6*DY1+A5*DY2+A5*DY3+A7*DY4+A7*DY5+A6*DY6
      DX2(KAREA2(6,IELH3))=A6*DY1+A7*DY2+A5*DY3+A5*DY4+A7*DY5+A6*DY6
      DX2(KAREA2(6,IELH4))=A6*DY1+A7*DY2+A7*DY3+A5*DY4+A5*DY5+A6*DY6
C
      DX2(KAREA2(3,IELH5))=A7*DY1+A5*DY2+A6*DY3+A7*DY4+A6*DY5+A5*DY6
      DX2(KAREA2(3,IELH6))=A7*DY1+A6*DY2+A5*DY3+A6*DY4+A7*DY5+A5*DY6
      DX2(KAREA2(3,IELH7))=A7*DY1+A7*DY2+A6*DY3+A5*DY4+A6*DY5+A5*DY6
      DX2(KAREA2(3,IELH8))=A7*DY1+A6*DY2+A7*DY3+A6*DY4+A5*DY5+A5*DY6
C
      IF (KADJ1(1,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+2D0*
     *                   (A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A4*DY6)
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+2D0*
     *                   (A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A4*DY6)
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+2D0*
     *                   (A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A4*DY6)
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+2D0*
     *                   (A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A4*DY6)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+2D0*
     *                   (A2*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+2D0*
     *                   (A2*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+2D0*
     *                   (A3*DY1+A1*DY2+A2*DY3+A4*DY4+A3*DY5+A2*DY6)
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+2D0*
     *                   (A3*DY1+A1*DY2+A3*DY3+A4*DY4+A2*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+2D0*
     *                   (A2*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A3*DY6)
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+2D0*
     *                   (A2*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A3*DY6)
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+2D0*
     *                   (A3*DY1+A3*DY2+A1*DY3+A2*DY4+A4*DY5+A2*DY6)
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+2D0*
     *                   (A3*DY1+A2*DY2+A1*DY3+A3*DY4+A4*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+2D0*
     *                   (A2*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+2D0*
     *                   (A2*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+2D0*
     *                   (A3*DY1+A4*DY2+A3*DY3+A1*DY4+A2*DY5+A2*DY6)
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+2D0*
     *                   (A3*DY1+A4*DY2+A2*DY3+A1*DY4+A3*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+2D0*
     *                   (A2*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+2D0*
     *                   (A2*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+2D0*
     *                   (A3*DY1+A2*DY2+A4*DY3+A3*DY4+A1*DY5+A2*DY6)
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+2D0*
     *                   (A3*DY1+A3*DY2+A4*DY3+A2*DY4+A1*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+2D0*
     *                   (A4*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6)
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+2D0*
     *                   (A4*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+2D0*
     *                   (A4*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+2D0*
     *                   (A4*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6)
      ENDIF
C
c
10    CONTINUE
C
C
      END
C
************************************************************************
      SUBROUTINE MP030(DX1,DX2,KAREA1,KAREA2,KADJ1,KADJ2,NAT1,NAT2,
     *                 NEL1,NEL2)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 0.5D0,A2= 0.125D0,A3=-0.125D0,A4= 0D0)
      PARAMETER (A5=0.5D0,A6=0D0,A7=0D0)
      DIMENSION DX1(1),DX2(1),KAREA1(NNAE,1),KAREA2(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DX2,NAT2)
C
      DO 10 IEL1=1,NEL1
C
      DY1=DX1(KAREA1(1,IEL1))
      DY2=DX1(KAREA1(2,IEL1))
      DY3=DX1(KAREA1(3,IEL1))
      DY4=DX1(KAREA1(4,IEL1))
      DY5=DX1(KAREA1(5,IEL1))
      DY6=DX1(KAREA1(6,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+
     *                      A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+
     *                      A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+
     *                      A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+
     *                      A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+
     *                      A2*DY1+A1*DY2+A3*DY3+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+
     *                      A2*DY1+A1*DY2+A2*DY3+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+
     *                      A3*DY1+A1*DY2+A2*DY3+A3*DY5+A2*DY6
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+
     *                      A3*DY1+A1*DY2+A3*DY3+A2*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+
     *                      A2*DY1+A2*DY2+A1*DY3+A3*DY4+A3*DY6
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+
     *                      A2*DY1+A3*DY2+A1*DY3+A2*DY4+A3*DY6
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+
     *                      A3*DY1+A3*DY2+A1*DY3+A2*DY4+A2*DY6
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+
     *                      A3*DY1+A2*DY2+A1*DY3+A3*DY4+A2*DY6
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+
     *                      A2*DY1+A2*DY3+A1*DY4+A3*DY5+A3*DY6
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+
     *                      A2*DY1+A3*DY3+A1*DY4+A2*DY5+A3*DY6
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+
     *                      A3*DY1+A3*DY3+A1*DY4+A2*DY5+A2*DY6
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+
     *                      A3*DY1+A2*DY3+A1*DY4+A3*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+
     *                      A2*DY1+A3*DY2+A2*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+
     *                      A2*DY1+A2*DY2+A3*DY4+A1*DY5+A3*DY6
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+
     *                      A3*DY1+A2*DY2+A3*DY4+A1*DY5+A2*DY6
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+
     *                      A3*DY1+A3*DY2+A2*DY4+A1*DY5+A2*DY6
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+
     *                      A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+
     *                      A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+
     *                      A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+
     *                      A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6
      ENDIF
C
C
      DX2(KAREA2(3,IELH1))=A5*(DY1+DY2)
      DX2(KAREA2(3,IELH2))=A5*(DY1+DY3)
      DX2(KAREA2(3,IELH3))=A5*(DY1+DY4)
      DX2(KAREA2(3,IELH4))=A5*(DY1+DY5)
C
      DX2(KAREA2(6,IELH1))=A5*(DY2+DY5)
      DX2(KAREA2(6,IELH2))=A5*(DY2+DY3)
      DX2(KAREA2(6,IELH3))=A5*(DY3+DY4)
      DX2(KAREA2(6,IELH4))=A5*(DY4+DY5)
C
      DX2(KAREA2(3,IELH5))=A5*(DY2+DY6)
      DX2(KAREA2(3,IELH6))=A5*(DY3+DY6)
      DX2(KAREA2(3,IELH7))=A5*(DY4+DY6)
      DX2(KAREA2(3,IELH8))=A5*(DY5+DY6)
C
      IF (KADJ1(1,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH1))=DX2(KAREA2(1,IELH1))+2D0*
     *                   (A1*DY1+A2*DY2+A3*DY3+A3*DY4+A2*DY5)
       DX2(KAREA2(1,IELH2))=DX2(KAREA2(1,IELH2))+2D0*
     *                   (A1*DY1+A2*DY2+A2*DY3+A3*DY4+A3*DY5)
       DX2(KAREA2(1,IELH3))=DX2(KAREA2(1,IELH3))+2D0*
     *                   (A1*DY1+A3*DY2+A2*DY3+A2*DY4+A3*DY5)
       DX2(KAREA2(1,IELH4))=DX2(KAREA2(1,IELH4))+2D0*
     *                   (A1*DY1+A3*DY2+A3*DY3+A2*DY4+A2*DY5)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH1))=DX2(KAREA2(2,IELH1))+2D0*
     *                   (A2*DY1+A1*DY2+A3*DY3+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH2))=DX2(KAREA2(5,IELH2))+2D0*
     *                   (A2*DY1+A1*DY2+A2*DY3+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH6))=DX2(KAREA2(5,IELH6))+2D0*
     *                   (A3*DY1+A1*DY2+A2*DY3+A3*DY5+A2*DY6)
       DX2(KAREA2(2,IELH5))=DX2(KAREA2(2,IELH5))+2D0*
     *                   (A3*DY1+A1*DY2+A3*DY3+A2*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH2))=DX2(KAREA2(2,IELH2))+2D0*
     *                   (A2*DY1+A2*DY2+A1*DY3+A3*DY4+A3*DY6)
       DX2(KAREA2(5,IELH3))=DX2(KAREA2(5,IELH3))+2D0*
     *                   (A2*DY1+A3*DY2+A1*DY3+A2*DY4+A3*DY6)
       DX2(KAREA2(5,IELH7))=DX2(KAREA2(5,IELH7))+2D0*
     *                   (A3*DY1+A3*DY2+A1*DY3+A2*DY4+A2*DY6)
       DX2(KAREA2(2,IELH6))=DX2(KAREA2(2,IELH6))+2D0*
     *                   (A3*DY1+A2*DY2+A1*DY3+A3*DY4+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH3))=DX2(KAREA2(2,IELH3))+2D0*
     *                   (A2*DY1+A2*DY3+A1*DY4+A3*DY5+A3*DY6)
       DX2(KAREA2(5,IELH4))=DX2(KAREA2(5,IELH4))+2D0*
     *                   (A2*DY1+A3*DY3+A1*DY4+A2*DY5+A3*DY6)
       DX2(KAREA2(5,IELH8))=DX2(KAREA2(5,IELH8))+2D0*
     *                   (A3*DY1+A3*DY3+A1*DY4+A2*DY5+A2*DY6)
       DX2(KAREA2(2,IELH7))=DX2(KAREA2(2,IELH7))+2D0*
     *                   (A3*DY1+A2*DY3+A1*DY4+A3*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN 
       DX2(KAREA2(2,IELH4))=DX2(KAREA2(2,IELH4))+2D0*
     *                   (A2*DY1+A3*DY2+A2*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH1))=DX2(KAREA2(5,IELH1))+2D0*
     *                   (A2*DY1+A2*DY2+A3*DY4+A1*DY5+A3*DY6)
       DX2(KAREA2(5,IELH5))=DX2(KAREA2(5,IELH5))+2D0*
     *                   (A3*DY1+A2*DY2+A3*DY4+A1*DY5+A2*DY6)
       DX2(KAREA2(2,IELH8))=DX2(KAREA2(2,IELH8))+2D0*
     *                   (A3*DY1+A3*DY2+A2*DY4+A1*DY5+A2*DY6)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN 
       DX2(KAREA2(1,IELH5))=DX2(KAREA2(1,IELH5))+2D0*
     *                   (A2*DY2+A3*DY3+A3*DY4+A2*DY5+A1*DY6)
       DX2(KAREA2(1,IELH6))=DX2(KAREA2(1,IELH6))+2D0*
     *                   (A2*DY2+A2*DY3+A3*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH7))=DX2(KAREA2(1,IELH7))+2D0*
     *                   (A3*DY2+A2*DY3+A2*DY4+A3*DY5+A1*DY6)
       DX2(KAREA2(1,IELH8))=DX2(KAREA2(1,IELH8))+2D0*
     *                   (A3*DY2+A3*DY3+A2*DY4+A2*DY5+A1*DY6)
      ENDIF
C
C
10    CONTINUE
C
C
      END
C
C
************************************************************************
      SUBROUTINE MR031(DF2,DF1,KAREA2,KAREA1,KADJ2,KADJ1,NAT2,NAT1,
     *                 NEL2,NEL1)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 11D0/24D0,A2= 7D0/48D0,A3=-5D0/48D0,A4=-1D0/24D0)
      PARAMETER (A5= 11D0/24D0,A6= 1D0/12D0,A7=-1D0/24D0)
      DIMENSION DF1(1),DF2(1),KAREA2(NNAE,1),KAREA1(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DF1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      DY1= DF2(KAREA2(1,IELH1))
      DY2= DF2(KAREA2(1,IELH2))
      DY3= DF2(KAREA2(1,IELH3))
      DY4= DF2(KAREA2(1,IELH4))
      DY5= DF2(KAREA2(2,IELH1))
      DY6= DF2(KAREA2(5,IELH2))
      DY7= DF2(KAREA2(5,IELH6))
      DY8= DF2(KAREA2(2,IELH5))
      DY9= DF2(KAREA2(2,IELH2))
      DY10=DF2(KAREA2(5,IELH3))
      DY11=DF2(KAREA2(5,IELH7))
      DY12=DF2(KAREA2(2,IELH6))
      DY13=DF2(KAREA2(2,IELH3))
      DY14=DF2(KAREA2(5,IELH4))
      DY15=DF2(KAREA2(5,IELH8))
      DY16=DF2(KAREA2(2,IELH7))
      DY17=DF2(KAREA2(2,IELH4))
      DY18=DF2(KAREA2(5,IELH1))
      DY19=DF2(KAREA2(5,IELH5))
      DY20=DF2(KAREA2(2,IELH8))
      DY21=DF2(KAREA2(1,IELH5))
      DY22=DF2(KAREA2(1,IELH6))
      DY23=DF2(KAREA2(1,IELH7))
      DY24=DF2(KAREA2(1,IELH8))
      DY25=DF2(KAREA2(3,IELH1))
      DY26=DF2(KAREA2(3,IELH2))
      DY27=DF2(KAREA2(3,IELH3))
      DY28=DF2(KAREA2(3,IELH4))
      DY29=DF2(KAREA2(6,IELH1))
      DY30=DF2(KAREA2(6,IELH2))
      DY31=DF2(KAREA2(6,IELH3))
      DY32=DF2(KAREA2(6,IELH4))
      DY33=DF2(KAREA2(3,IELH5))
      DY34=DF2(KAREA2(3,IELH6))
      DY35=DF2(KAREA2(3,IELH7))
      DY36=DF2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DF1(IA1)= DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A4*(DY21+DY22+DY23+DY24)
     *          +A5*(DY25+DY26+DY27+DY28)
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY33+DY34+DY35+DY36)     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DF1(IA2)= DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A4*(DY13+DY14+DY15+DY16)
     *          +A5*(DY25+DY29+DY30+DY33)
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY27+DY31+DY32+DY35)     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DF1(IA3)= DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A4*(DY17+DY18+DY19+DY20)
     *          +A5*(DY26+DY30+DY31+DY34)
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY28+DY29+DY32+DY36)     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DF1(IA4)= DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A4*(DY5+DY6+DY7+DY8)
     *          +A5*(DY27+DY31+DY32+DY35)  
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY25+DY29+DY30+DY33)
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DF1(IA5)= DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A4*(DY9+DY10+DY11+DY12)
     *          +A5*(DY28+DY29+DY32+DY36)     
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY26+DY30+DY31+DY34)
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DF1(IA6)= DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A4*(DY1+DY2+DY3+DY4)
     *          +A5*(DY33+DY34+DY35+DY36)     
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY25+DY26+DY27+DY28)
      ENDIF
C
      IF (KADJ1(1,IEL1).EQ.0) THEN
       DF1(IA1)= 2D0*(DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A4*(DY21+DY22+DY23+DY24)
     *          +A5*(DY25+DY26+DY27+DY28)
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY33+DY34+DY35+DY36))     
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN
       DF1(IA2)= 2D0*(DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A4*(DY13+DY14+DY15+DY16)
     *          +A5*(DY25+DY29+DY30+DY33)
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY27+DY31+DY32+DY35))     
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN
       DF1(IA3)= 2D0*(DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A4*(DY17+DY18+DY19+DY20)
     *          +A5*(DY26+DY30+DY31+DY34)
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY28+DY29+DY32+DY36))     
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN
       DF1(IA4)= 2D0*(DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A4*(DY5+DY6+DY7+DY8)
     *          +A5*(DY27+DY31+DY32+DY35)  
     *          +A6*(DY26+DY28+DY34+DY36)
     *          +A7*(DY25+DY29+DY30+DY33))
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN
       DF1(IA5)= 2D0*(DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A4*(DY9+DY10+DY11+DY12)
     *          +A5*(DY28+DY29+DY32+DY36)     
     *          +A6*(DY25+DY27+DY33+DY35)
     *          +A7*(DY26+DY30+DY31+DY34))
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN
       DF1(IA6)= 2D0*(DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A4*(DY1+DY2+DY3+DY4)
     *          +A5*(DY33+DY34+DY35+DY36)     
     *          +A6*(DY29+DY30+DY31+DY32)
     *          +A7*(DY25+DY26+DY27+DY28))
      ENDIF
C
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MR030(DF2,DF1,KAREA2,KAREA1,KADJ2,KADJ1,NAT2,NAT1,
     *                 NEL2,NEL1)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      PARAMETER (A1= 0.5D0,A2= 0.125D0,A3=-0.125D0,A4= 0D0)
      PARAMETER (A5=0.5D0,A6=0D0,A7=0D0)
      DIMENSION DF1(1),DF2(1),KAREA2(NNAE,1),KAREA1(NNAE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCL1(DF1,NAT1)
C
      DO 10 IEL1=1,NEL1
C
      IA1=KAREA1(1,IEL1)
      IA2=KAREA1(2,IEL1)
      IA3=KAREA1(3,IEL1)
      IA4=KAREA1(4,IEL1)
      IA5=KAREA1(5,IEL1)
      IA6=KAREA1(6,IEL1)
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C
      DY1= DF2(KAREA2(1,IELH1))
      DY2= DF2(KAREA2(1,IELH2))
      DY3= DF2(KAREA2(1,IELH3))
      DY4= DF2(KAREA2(1,IELH4))
      DY5= DF2(KAREA2(2,IELH1))
      DY6= DF2(KAREA2(5,IELH2))
      DY7= DF2(KAREA2(5,IELH6))
      DY8= DF2(KAREA2(2,IELH5))
      DY9= DF2(KAREA2(2,IELH2))
      DY10=DF2(KAREA2(5,IELH3))
      DY11=DF2(KAREA2(5,IELH7))
      DY12=DF2(KAREA2(2,IELH6))
      DY13=DF2(KAREA2(2,IELH3))
      DY14=DF2(KAREA2(5,IELH4))
      DY15=DF2(KAREA2(5,IELH8))
      DY16=DF2(KAREA2(2,IELH7))
      DY17=DF2(KAREA2(2,IELH4))
      DY18=DF2(KAREA2(5,IELH1))
      DY19=DF2(KAREA2(5,IELH5))
      DY20=DF2(KAREA2(2,IELH8))
      DY21=DF2(KAREA2(1,IELH5))
      DY22=DF2(KAREA2(1,IELH6))
      DY23=DF2(KAREA2(1,IELH7))
      DY24=DF2(KAREA2(1,IELH8))
      DY25=DF2(KAREA2(3,IELH1))
      DY26=DF2(KAREA2(3,IELH2))
      DY27=DF2(KAREA2(3,IELH3))
      DY28=DF2(KAREA2(3,IELH4))
      DY29=DF2(KAREA2(6,IELH1))
      DY30=DF2(KAREA2(6,IELH2))
      DY31=DF2(KAREA2(6,IELH3))
      DY32=DF2(KAREA2(6,IELH4))
      DY33=DF2(KAREA2(3,IELH5))
      DY34=DF2(KAREA2(3,IELH6))
      DY35=DF2(KAREA2(3,IELH7))
      DY36=DF2(KAREA2(3,IELH8))
C
      IF (KADJ1(1,IEL1).NE.0) THEN
       DF1(IA1)= DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A5*(DY25+DY26+DY27+DY28)
      ENDIF
C
C
      IF (KADJ1(2,IEL1).NE.0) THEN
       DF1(IA2)= DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A5*(DY25+DY29+DY30+DY33)
      ENDIF
C
C
      IF (KADJ1(3,IEL1).NE.0) THEN
       DF1(IA3)= DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A5*(DY26+DY30+DY31+DY34)
      ENDIF
C
C
      IF (KADJ1(4,IEL1).NE.0) THEN
       DF1(IA4)= DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A5*(DY27+DY31+DY32+DY35)  
      ENDIF
C
C
      IF (KADJ1(5,IEL1).NE.0) THEN
       DF1(IA5)= DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A5*(DY28+DY29+DY32+DY36)     
      ENDIF
C
C
      IF (KADJ1(6,IEL1).NE.0) THEN
       DF1(IA6)= DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A5*(DY33+DY34+DY35+DY36)     
      ENDIF
C
      IF (KADJ1(1,IEL1).EQ.0) THEN
       DF1(IA1)= 2D0*(DF1(IA1)
     *          +A1*(DY1+DY2+DY3+DY4)
     *          +A2*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A3*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A5*(DY25+DY26+DY27+DY28))
      ENDIF
C
C
      IF (KADJ1(2,IEL1).EQ.0) THEN
       DF1(IA2)= 2D0*(DF1(IA2)
     *          +A1*(DY5+DY6+DY7+DY8)
     *          +A2*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A3*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A5*(DY25+DY29+DY30+DY33))
      ENDIF
C
C
      IF (KADJ1(3,IEL1).EQ.0) THEN
       DF1(IA3)= 2D0*(DF1(IA3)
     *          +A1*(DY9+DY10+DY11+DY12)
     *          +A2*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A3*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A5*(DY26+DY30+DY31+DY34))
      ENDIF
C
C
      IF (KADJ1(4,IEL1).EQ.0) THEN
       DF1(IA4)= 2D0*(DF1(IA4)
     *          +A1*(DY13+DY14+DY15+DY16)
     *          +A2*(DY3+DY4+DY10+DY11+DY23+DY24+DY17+DY20)
     *          +A3*(DY1+DY2+DY9+DY12+DY21+DY22+DY18+DY19)
     *          +A5*(DY27+DY31+DY32+DY35))  
      ENDIF
C
C
      IF (KADJ1(5,IEL1).EQ.0) THEN
       DF1(IA5)= 2D0*(DF1(IA5)
     *          +A1*(DY17+DY18+DY19+DY20)
     *          +A2*(DY1+DY4+DY5+DY8+DY21+DY24+DY14+DY15)
     *          +A3*(DY2+DY3+DY6+DY7+DY22+DY23+DY13+DY16)
     *          +A5*(DY28+DY29+DY32+DY36))     
      ENDIF
C
C
      IF (KADJ1(6,IEL1).EQ.0) THEN
       DF1(IA6)= 2D0*(DF1(IA6)
     *          +A1*(DY21+DY22+DY23+DY24)
     *          +A2*(DY7+DY8+DY11+DY12+DY15+DY16+DY19+DY20)
     *          +A3*(DY5+DY6+DY9+DY10+DY13+DY14+DY17+DY18)
     *          +A5*(DY33+DY34+DY35+DY36))    
      ENDIF
C
C
C
10    CONTINUE
C
C
      END
************************************************************************
      SUBROUTINE MP010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NNAE=6)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNAE,*),KADJ2(NNAE,*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(3,IELH5)
      IELH7=KADJ2(3,IELH6)
      IELH8=KADJ2(3,IELH7)
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
      DP2(IELH5)=DPH
      DP2(IELH6)=DPH
      DP2(IELH7)=DPH
      DP2(IELH8)=DPH
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NNAE=6)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNAE,*),KADJ2(NNAE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
C
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)+
     *           DP2(IELH5)+DP2(IELH6)+DP2(IELH7)+DP2(IELH8)
C
10    CONTINUE
C
C
      END
c
c
************************************************************************
      SUBROUTINE MP511(DX1,DX2,KVERT1,KVERT2,KADJ1,KADJ2,
     *                 NVT1,NVT2,NEL1,NEL2)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DX1(1),DX2(1),KVERT1(NNVE,1),KVERT2(NNVE,1),
     *          KADJ1(NNAE,1),KADJ2(NNAE,1)
C
C
C
      CALL LCP1(DX1,DX2,NVT1)
C
      DO 10 IEL1=1,NEL1
C
      DXH1=0.5D0*DX1(KVERT1(1,IEL1))
      DXH2=0.5D0*DX1(KVERT1(2,IEL1))
      DXH3=0.5D0*DX1(KVERT1(3,IEL1))
      DXH4=0.5D0*DX1(KVERT1(4,IEL1))
      DXH5=0.5D0*DX1(KVERT1(5,IEL1))
      DXH6=0.5D0*DX1(KVERT1(6,IEL1))
      DXH7=0.5D0*DX1(KVERT1(7,IEL1))
      DXH8=0.5D0*DX1(KVERT1(8,IEL1))
c
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C
C
      IF ((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)) 
     *   DX2(KVERT2(3,IELH1))=0.5D0*(DXH1+DXH2+DXH3+DXH4)
C
      IF ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0)) 
     *   DX2(KVERT2(6,IELH1))=0.5D0*(DXH1+DXH2+DXH5+DXH6)
C
      IF ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0)) 
     *   DX2(KVERT2(6,IELH2))=0.5D0*(DXH2+DXH3+DXH6+DXH7)
C
      IF ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0)) 
     *   DX2(KVERT2(6,IELH3))=0.5D0*(DXH3+DXH4+DXH7+DXH8)
C
      IF ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0)) 
     *   DX2(KVERT2(6,IELH4))=0.5D0*(DXH1+DXH4+DXH5+DXH8)
C
      IF ((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)) 
     *   DX2(KVERT2(3,IELH5))=0.5D0*(DXH5+DXH6+DXH7+DXH8)
C
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH1))=DXH1+DXH2
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH2))=DXH2+DXH3
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH3))=DXH3+DXH4
C
      IF (((KADJ1(1,IEL1).GT.IEL1).OR.(KADJ1(1,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH4))=DXH1+DXH4
C
      IF (((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   DX2(KVERT2(5,IELH1))=DXH1+DXH5
C
      IF (((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   DX2(KVERT2(5,IELH2))=DXH2+DXH6
C
      IF (((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   DX2(KVERT2(5,IELH3))=DXH3+DXH7
C
      IF (((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   DX2(KVERT2(5,IELH4))=DXH4+DXH8
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(2,IEL1).GT.IEL1).OR.(KADJ1(2,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH5))=DXH5+DXH6
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(3,IEL1).GT.IEL1).OR.(KADJ1(3,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH6))=DXH6+DXH7
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(4,IEL1).GT.IEL1).OR.(KADJ1(4,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH7))=DXH7+DXH8
C
      IF (((KADJ1(6,IEL1).GT.IEL1).OR.(KADJ1(6,IEL1).EQ.0)).AND.
     *    ((KADJ1(5,IEL1).GT.IEL1).OR.(KADJ1(5,IEL1).EQ.0))) 
     *   DX2(KVERT2(2,IELH8))=DXH5+DXH8
C
      DX2(KVERT2(7,IELH1))=0.25D0*(DXH1+DXH2+DXH3+DXH4+
     *                             DXH5+DXH6+DXH7+DXH8)
C
10    CONTINUE
C
C
C
      END
c
************************************************************************
      SUBROUTINE MR511(DF2,DF1,KVERT2,KVERT1,KADJ2,KADJ1,
     *                 NVT2,NVT1,NEL2,NEL1)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DF1(1),DF2(1),KVERT2(NNVE,*),KVERT1(NNVE,*),
     *          KADJ2(NNAE,*),KADJ1(NNAE,*)
C
C
C
      CALL LCP1(DF2,DF1,NVT1)
      DO 10 IEL2=1,NEL2
C
      I1=KVERT2(1,IEL2)
      I3=KVERT2(3,IEL2)
      I6=KVERT2(6,IEL2)
      I7=KVERT2(7,IEL2)
      I8=KVERT2(8,IEL2)
C
      IF (KADJ2(1,IEL2).EQ.0) THEN
       A3=0.25D0
      ELSE
       A3=0.125D0
      ENDIF
C
      IF (KADJ2(2,IEL2).EQ.0) THEN
       A6=0.25D0
      ELSE
       A6=0.125D0
      ENDIF
C
      IF (KADJ2(5,IEL2).EQ.0) THEN
       A8=0.25D0
      ELSE
       A8=0.125D0
      ENDIF
C
      DF1(I1)=DF1(I1)+A3*DF2(I3)+A6*DF2(I6)+0.125D0*DF2(I7)+A8*DF2(I8)
C
10    CONTINUE
C
      IVT=NVT1
C
      DO 20 IEL1=1,NEL1
      IELH1=IEL1
      IELH2=KADJ2(3,IELH1)
      IELH3=KADJ2(3,IELH2)
      IELH4=KADJ2(3,IELH3)
      IELH5=KADJ2(6,IELH1)
      IELH6=KADJ2(6,IELH2)
      IELH7=KADJ2(6,IELH3)
      IELH8=KADJ2(6,IELH4)
C
c      IELH1=IEL1
c      IELH2=KADJ2(3,IELH1)
c      IELH3=KADJ2(3,IELH2)
c      IELH4=KADJ2(3,IELH3)
c      IELH5=KADJ2(6,IELH1)
c      IELH6=KADJ2(3,IELH5)
c      IELH7=KADJ2(3,IELH6)
c      IELH8=KADJ2(3,IELH7)
C
      I1 =KVERT2(2,IELH1)
      I2 =KVERT2(2,IELH2)
      I3 =KVERT2(2,IELH3)
      I4 =KVERT2(2,IELH4)
      I5 =KVERT2(5,IELH1)
      I6 =KVERT2(5,IELH2)
      I7 =KVERT2(5,IELH3)
      I8 =KVERT2(5,IELH4)
      I9 =KVERT2(2,IELH5)
      I10=KVERT2(2,IELH6)
      I11=KVERT2(2,IELH7)
      I12=KVERT2(2,IELH8)
C
      J1=KVERT1(1,IEL1)
      J2=KVERT1(2,IEL1)
      J3=KVERT1(3,IEL1)
      J4=KVERT1(4,IEL1)
      J5=KVERT1(5,IEL1)
      J6=KVERT1(6,IEL1)
      J7=KVERT1(7,IEL1)
      J8=KVERT1(8,IEL1)
C
      IF (I1.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I1)
       DF1(J2)=DF1(J2)+0.5D0*DF2(I1)
       IVT=IVT+1
      ENDIF
C
      IF (I2.GT.IVT) THEN
       DF1(J2)=DF1(J2)+0.5D0*DF2(I2)
       DF1(J3)=DF1(J3)+0.5D0*DF2(I2)
       IVT=IVT+1
      ENDIF
C
      IF (I3.GT.IVT) THEN
       DF1(J3)=DF1(J3)+0.5D0*DF2(I3)
       DF1(J4)=DF1(J4)+0.5D0*DF2(I3)
       IVT=IVT+1
      ENDIF
C
      IF (I4.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I4)
       DF1(J4)=DF1(J4)+0.5D0*DF2(I4)
       IVT=IVT+1
      ENDIF
C
      IF (I5.GT.IVT) THEN
       DF1(J1)=DF1(J1)+0.5D0*DF2(I5)
       DF1(J5)=DF1(J5)+0.5D0*DF2(I5)
       IVT=IVT+1
      ENDIF
C
      IF (I6.GT.IVT) THEN
       DF1(J2)=DF1(J2)+0.5D0*DF2(I6)
       DF1(J6)=DF1(J6)+0.5D0*DF2(I6)
       IVT=IVT+1
      ENDIF
C
      IF (I7.GT.IVT) THEN
       DF1(J3)=DF1(J3)+0.5D0*DF2(I7)
       DF1(J7)=DF1(J7)+0.5D0*DF2(I7)
       IVT=IVT+1
      ENDIF
C
      IF (I8.GT.IVT) THEN
       DF1(J4)=DF1(J4)+0.5D0*DF2(I8)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I8)
       IVT=IVT+1
      ENDIF
C
      IF (I9.GT.IVT) THEN
       DF1(J5)=DF1(J5)+0.5D0*DF2(I9)
       DF1(J6)=DF1(J6)+0.5D0*DF2(I9)
       IVT=IVT+1
      ENDIF
C
      IF (I10.GT.IVT) THEN
       DF1(J6)=DF1(J6)+0.5D0*DF2(I10)
       DF1(J7)=DF1(J7)+0.5D0*DF2(I10)
       IVT=IVT+1
      ENDIF
C
      IF (I11.GT.IVT) THEN
       DF1(J7)=DF1(J7)+0.5D0*DF2(I11)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I11)
       IVT=IVT+1
      ENDIF
C
      IF (I12.GT.IVT) THEN
       DF1(J5)=DF1(J5)+0.5D0*DF2(I12)
       DF1(J8)=DF1(J8)+0.5D0*DF2(I12)
       IVT=IVT+1
      ENDIF
C
20    CONTINUE
C
C
C
      END
c
c
************************************************************************
      SUBROUTINE YAXU (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DAX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DAX,NEQ)
      ENDIF
C
C
      CALL LAX37(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *           KWORK(L(KLLDA(ILEV))),NEQ,DX,DAX,A1,A2)
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DAX,NEQ)
      ENDIF
C
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE YDBCU (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX=(D1,D2,D3) to 0
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      KABD=L(KLABD(ILEV))
      NABD=KNABD(ILEV)
      KNPR=L(KLNPR(ILEV))
      NVT=KNVT(ILEV)        
C
      CALL BDRY0(DX,KWORK(KABD),NABD,KWORK(KNPR),NVT)
C
      END
c
c
c
************************************************************************
      SUBROUTINE YEXU (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX37,YIA137,YID137,YIF137
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.1) THEN
       CALL IC037(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSLU,ITE,DMPUSL,
     *            RLXSLU)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.2) THEN
C
      OMEGA=RLXSLU
C
      IREQ=5*NEQ
      BNOCON=OMEGA.GT.0D0
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE IF (OMEGA.LT.2D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,I000,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ENDIF
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.3) THEN
       CALL IF027(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSLU,ITE,DMPUSL,RLXSLU,RHO)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.4) THEN
C
      OMEGA=RLXSLU
C
      IREQ=5*NEQ
      BNOCON=.TRUE.
      IREQ=MAX(IREQ,5)
      CALL ZNEW(IREQ,1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
      L1=L(LWORK)
      L2=L1+NEQ
      L3=L2+NEQ
      L4=L3+NEQ
      L5=L4+NEQ
C
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      KXYPAR(4)=L(KLAILU(ILEV))
      DXYPAR(1)=RLXSLU
C
      CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIF137,BNOCON,
     *           DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *           DWORK(L5),1,RHO)
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YEXAU (DX,DB,DD,NEQ,RHO,EPS,ITECG)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX37,YIA137,YID137,YIF137
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.1) THEN
       CALL IC037(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NSLU,ITE,DMPUSL,
     *            RLXSLU)
       RHO  =1D0
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.2) THEN
C
      OMEGA=RLXSLU
C
      IREQ=1*NEQ
      CALL ZNEW(IREQ,1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
C
      L1=L(LWORK)
      L2=KAUX1
      L3=KAUX2
      L4=L(LD1)
      L5=L(LD2)
C
      BNOCON=OMEGA.GT.0D0
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      DXYPAR(1)=OMEGA
C
      IF (OMEGA.EQ.0D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIA137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE IF (OMEGA.LT.2D0) THEN
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YID137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ELSE
       CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,I000,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *            DWORK(L5),1,RHO)
      ENDIF
      ITECG=ITE
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
      IF (ISLU.EQ.3) THEN
       CALL IF027(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSLU,ITE,DMPUSL,RLXSLU,RHO)
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLU.EQ.4) THEN
C
      IREQ=1*NEQ
      CALL ZNEW(IREQ,1,LWORK,'WORKCG')
      IF (IER.NE.0) GOTO 99999
C
      L1=L(LWORK)
      L2=KAUX1
      L3=KAUX2
      L4=L(LD1)
      L5=L(LD2)
C
      BNOCON=.TRUE.
      KXYPAR(1)=L(KLA(ILEV))
      KXYPAR(2)=L(KLCOLA(ILEV))
      KXYPAR(3)=L(KLLDA(ILEV))
      KXYPAR(4)=L(KLAILU(ILEV))
      DXYPAR(1)=RLXSLU
C
      CALL II010(DX,DB,NEQ,NSLU,ITE,DMPUSL,YLAX37,YIF137,BNOCON,
     *           DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),
     *           DWORK(L5),1,RHO)
      ITECG=ITE
C
      IER1=IER
      CALL ZDISP(0,LWORK,'WORKCG')
      IER=IER1
      GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
c
C
C
99999 END
c
c
************************************************************************
      SUBROUTINE YPROLU (DUC,DUF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DUF(*),DUC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
      LLM1=L(KLAREA(ILEV-1))
      LLM2=L(KLAREA(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      LLA2=L(KLADJ(ILEV))
      NAT1=KNAT(ILEV-1)
      NAT2=KNAT(ILEV)
      NEL1=KNEL(ILEV-1)
      NEL2=KNEL(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTU.EQ.1) THEN
      CALL MP031(DUC,DUF,KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),
     *            KWORK(LLA2),NAT1,NAT2,NEL1,NEL2)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.2) THEN
      CALL MP030(DUC,DUF,KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),
     *            KWORK(LLA2),NAT1,NAT2,NEL1,NEL2)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YRESTU(DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      LLM1=L(KLAREA(ILEV))
      LLM2=L(KLAREA(ILEV+1))
      LLA1=L(KLADJ(ILEV))
      LLA2=L(KLADJ(ILEV+1))
      NAT1=KNAT(ILEV)
      NAT2=KNAT(ILEV+1)
      NEL2=KNEL(ILEV+1)
      NEL1=KNEL(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTU.EQ.1) THEN
      CALL MR031(DDF,DDC,KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),
     *            KWORK(LLA1),NAT2,NAT1,NEL2,NEL1)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
      IF (IINTU.EQ.2) THEN
      CALL MR030(DDF,DDC,KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),
     *            KWORK(LLA1),NAT2,NAT1,NEL2,NEL1)
       GOTO 99999
      ENDIF
C-----------------------------------------------------------------------
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSMU (DX,DB,DD,NEQ,NSM)  
************************************************************************
*
*   Purpose: - performs NSM smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C
C
      IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISMU.EQ.1) THEN
C
      CALL ParIA237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *     KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMU,ILEV)
      GOTO 999
C
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.2) THEN
C
      IF (RLXSMU.EQ.1D0) THEN
       CALL ParIB237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NSM,ILEV)
      ELSE
       CALL ParIC237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMU,ILEV)
      ENDIF
      GOTO 999
C
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.3) THEN
       CALL ParID237(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMU,ILEV)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
      IF (ISMU.EQ.4) THEN
       CALL ParIF237(VWORK(L(KLA(ILEV))),VWORK(L(KLAILU(ILEV))),
     *            KWORK(L(KLCOLA(ILEV))),KWORK(L(KLLDA(ILEV))),
     *            DX,DB,DD,NEQ,NSM,RLXSMU,ILEV)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C
999   IF (ISORTU.GT.0) THEN
       KTRA1=L(KLTRA1(ILEV))
       KTRA2=L(KLTRA2(ILEV))
       CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
       CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
       CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
      ENDIF
C
C
C 
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSTEPU(DX,DD,DB,NEQ,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          ALFMIN and ALFMAX (COMMON /RPARM/)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C *** Setting of given ALPHA
C
      IF (AMINU.EQ.AMAXU) THEN
       ALPHA=AMINU
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
C
      IF (AMINU.NE.AMAXU) THEN
C
       IF (ISORTU.GT.0) THEN
        KTRA1=L(KLTRA1(ILEV))
        KTRA2=L(KLTRA2(ILEV))
        CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
        CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
        CALL VECSRT(DD,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(L(LD1)),DD,NEQ)
       ENDIF
C
       CALL ParVYAX7(VWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *            KWORK(L(KLLDA(ILEV))),NEQ,DD,DX,DD,DYAX1,DYAX2)
       CALL LSP1(DB,DD,NEQ,DBX)
       ALPHA=(DBX-DYAX1)/DYAX2
       IF (ALPHA.LT.AMINU) ALPHA=AMINU
       IF (ALPHA.GT.AMAXU) ALPHA=AMAXU
C
       IF (ISORTU.GT.0) THEN
        KTRA1=L(KLTRA1(ILEV))
        KTRA2=L(KLTRA2(ILEV))
        CALL VECSRT(DX,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DX,NEQ)
        CALL VECSRT(DB,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DB,NEQ)
        CALL VECSRT(DD,DWORK(L(LD1)),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(L(LD1)),DD,NEQ)
       ENDIF
C
       RETURN
C
      ENDIF
C
C
      END
************************************************************************
      SUBROUTINE YAXP (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DAX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C-----------------------------------------------------------------------
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DAX,NEQ)
      ENDIF
C
C
      CALL LAX17(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *           KWORK(L(KLLDC(ILEV))),NEQ,DX,DAX,A1,A2)
      CALL TOL20A(DAX,VWORK(L(KLVOL(ILEV))),NEQ,INEUM)
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DAX,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DAX,NEQ)
      ENDIF
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YDBCP (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX to zero
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      END
c
c
c
************************************************************************
      SUBROUTINE YEXP (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX17,YIA117,YID117,YIF137
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.1) THEN
       CALL IC017(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSLP,ITE,DMPPSL,
     *            RLXSLP)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.2) THEN
C
       OMEGA=RLXSLP
C
       IREQ=4*NEQ
       BNOCON=OMEGA.LT.0D0
       IF (BNOCON) IREQ=3*NEQ
       IREQ=MAX(IREQ,4)
       CALL ZNEW(IREQ,1,LWORK,'WORKCG')
       IF (IER.NE.0) GOTO 99999
       L1=L(LWORK)
       L2=L1+NEQ
       L3=L2+NEQ
       L4=L3+NEQ
       IF (BNOCON) L4=L1
C
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       DXYPAR(1)=OMEGA
C
       IF (OMEGA.EQ.0D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YIA117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ELSE IF (OMEGA.LT.2D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YID117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ELSE
        CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,I000,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
       ENDIF
C
       IER1=IER
       CALL ZDISP(0,LWORK,'WORKCG')
       IER=IER1
       GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.3) THEN
       CALL IF037(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSLP,ITE,DMPPSL,RLXSLP,RHO)
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.4) THEN
C
       IREQ=4*NEQ
       BNOCON=.FALSE.
       IREQ=MAX(IREQ,4)
       CALL ZNEW(IREQ,1,LWORK,'WORKCG')
       IF (IER.NE.0) GOTO 99999
       L1=L(LWORK)
       L2=L1+NEQ
       L3=L2+NEQ
       L4=L3+NEQ
C
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       KXYPAR(4)=L(KLCILU(ILEV))
       DXYPAR(1)=RLXSLP
C
       CALL IE010(DX,DB,NEQ,NSLP,ITE,DMPPSL,YLAX17,YIF137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),1,RHO)
C
       IER1=IER
       CALL ZDISP(0,LWORK,'WORKCG')
       IER=IER1
       GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLVOL(ILEV))),NEQ,INEUM)
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YEXAP (DX,DB,DD,NEQ,RHO,EPS,ITECG)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
C
C *** User COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL I000,YLAX17,YIA117,YID117,YIF137
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.1) THEN
       CALL IC017(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,NEQ,NSLP,ITE,EPS,RLXSLP)
       RHO  =1D0
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
      IF (ISLP.EQ.2) THEN
C
       OMEGA=RLXSLP
C
       L1=KAUX1
       L2=KAUX2
       L3=L(LD1)
       L4=L(LD2)
C
       BNOCON=OMEGA.LT.0D0
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       DXYPAR(1)=OMEGA
C
       IF (OMEGA.EQ.0D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YIA117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ELSE IF (OMEGA.LT.2D0) THEN
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YID117,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ELSE
        CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,I000,BNOCON,
     *             DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ENDIF
       ITECG=ITE
C
       GOTO 999
C
      ENDIF 
C
C-----------------------------------------------------------------------
C
      IF (ISLP.EQ.3) THEN
       CALL IF037(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSLP,ITE,EPS,RLXSLP,RHO)
       ITECG=ITE
       GOTO 999
      ENDIF 
C
C-----------------------------------------------------------------------
C
      IF (ISLP.EQ.4) THEN
C
       L1=KAUX1
       L2=KAUX2
       L3=L(LD1)
       L4=L(LD2)
C
       BNOCON=.FALSE.
       KXYPAR(1)=L(KLC(ILEV))
       KXYPAR(2)=L(KLCOLC(ILEV))
       KXYPAR(3)=L(KLLDC(ILEV))
       KXYPAR(4)=L(KLCILU(ILEV))
       DXYPAR(1)=RLXSLP
C
       CALL IE010(DX,DB,NEQ,NSLP,ITE,EPS,YLAX17,YIF137,BNOCON,
     *            DWORK(L1),DWORK(L2),DWORK(L3),DWORK(L4),0,RHO)
       ITECG=ITE
C
       GOTO 999
C
      ENDIF 
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLVOL(ILEV))),NEQ,INEUM)
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YPROLP (DPC,DPF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DPF(*),DPC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
      LLV1=L(KLVERT(ILEV-1))
      LLV2=L(KLVERT(ILEV))
      LLM1=L(KLAREA(ILEV-1))
      LLM2=L(KLAREA(ILEV))
      LLA1=L(KLADJ(ILEV-1))
      LLA2=L(KLADJ(ILEV))
      NVT1=KNVT(ILEV-1)
      NVT2=KNVT(ILEV)
      NEL1=KNEL(ILEV-1)
      NEL2=KNEL(ILEV)
      NAT1=KNAT(ILEV-1)
      NAT2=KNAT(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.1) THEN
C
       CALL MP010(DPC,DPF,KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
       GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
c
      IF (IINTP.EQ.2) THEN
C
c      KPLC=L(LD1)
c      KPLF=L(LD2)
      CALL ZNEW(NVT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NVT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NVT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DPC,DWORK(L(LPLC)),VWORK(L(KLVOL(ILEV-1))),
     *           VWORK(L(LAUXCL)),KWORK(LLV1),NEL1,NVT1,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MP511(DWORK(L(LPLC)),DWORK(L(LPLF)),KWORK(LLV1),
     *           KWORK(LLV2),
     *           KWORK(LLA1),KWORK(LLA2),NVT1,NVT2,NEL1,NEL2)
C
      CALL ZNEW(NVT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DPF,DWORK(L(LPLF)),VWORK(L(KLVOL(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLV2),NEL2,NVT2,1)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.3) THEN
C
      CALL ZNEW(NAT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NAT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NAT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DPC,DWORK(L(LPLC)),VWORK(L(KLVOL(ILEV-1))),
     *           VWORK(L(LAUXCL)),KWORK(LLM1),NEL1,NAT1,NVT1,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MP031(DWORK(L(LPLC)),DWORK(L(LPLF)),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           NAT1,NAT2,NEL1,NEL2)
C
      CALL ZNEW(NAT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DPF,DWORK(L(LPLF)),VWORK(L(KLVOL(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLM2),NEL2,NAT2,NVT2,1)
C     
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.4) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
C
      CALL C2N2DM(DPC,DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *            NEL1,NAT1,NVT1,0)
C
      CALL MP031(DWORK(KPLC),DWORK(KPLF),
     *           KWORK(LLM1),KWORK(LLM2),KWORK(LLA1),KWORK(LLA2),
     *           NAT1,NAT2,NEL1,NEL2)
C
      CALL C2N2DM(DPF,DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *            NEL2,NAT2,NVT2,1)
C
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C
C
99999 CALL TOL20A(DPF,VWORK(L(KLVOL(ILEV))),NEL2,INEUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE YRESTP (DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON block
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
      LLV2=L(KLVERT(ILEV+1))
      LLV1=L(KLVERT(ILEV))
      LLA2=L(KLADJ(ILEV+1))
      LLA1=L(KLADJ(ILEV))
      LLM2=L(KLAREA(ILEV+1))
      LLM1=L(KLAREA(ILEV))
      NVT2=KNVT(ILEV+1)
      NVT1=KNVT(ILEV)
      NAT2=KNAT(ILEV+1)
      NAT1=KNAT(ILEV)
      NEL2=KNEL(ILEV+1)
      NEL1=KNEL(ILEV)
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.1) THEN
C
       CALL MR010(DDC,DDF,KWORK(LLA1),KWORK(LLA2),NEL1,NEL2)
       GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
      IF (IINTP.EQ.2) THEN
C
C      KPLC=L(LD1)
C      KPLF=L(LD2)
      CALL ZNEW(NVT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NVT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NVT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DDF,DWORK(L(LPLF)),VWORK(L(KLVOL(ILEV+1))),
     *           VWORK(L(LAUXCL)),KWORK(LLV2),NEL2,NVT2,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MR511(DWORK(L(LPLF)),DWORK(L(LPLC)),KWORK(LLV2),
     *           KWORK(LLV1),
     *           KWORK(LLA2),KWORK(LLA1),NVT2,NVT1,NEL2,NEL1)
C
      CALL ZNEW(NVT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2L2D(DDC,DWORK(L(LPLC)),VWORK(L(KLVOL(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLV1),NEL1,NVT1,1)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.3) THEN
C
      CALL ZNEW(NAT1,1,LPLC  ,'DPLC  ')
      CALL ZNEW(NAT2,1,LPLF  ,'DPLF  ')
      CALL ZNEW(NAT2,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DDF,DWORK(L(LPLF)),VWORK(L(KLVOL(ILEV+1))),
     *           VWORK(L(LAUXCL)),KWORK(LLM2),NEL2,NAT2,NVT2,0)
C
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL MR031(DWORK(L(LPLF)),DWORK(L(LPLC)),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           NAT2,NAT1,NEL2,NEL1)
C     
      CALL ZNEW(NAT1,2,LAUXCL,'DAUXCL')
      IF (IER.NE.0) GOTO 99999
C
      CALL C2N2D(DDC,DWORK(L(LPLC)),VWORK(L(KLVOL(ILEV))),
     *           VWORK(L(LAUXCL)),KWORK(LLM1),NEL1,NAT1,NVT1,1)
C     
      CALL ZDISP (0,LAUXCL,'DAUXCL')
      CALL ZDISP (0,LPLF  ,'DPLF  ')
      CALL ZDISP (0,LPLC  ,'DPLC  ')
      IF (IER.NE.0) GOTO 99999
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (IINTP.EQ.4) THEN
C
      KPLC=L(LD1)
      KPLF=L(LD2)
C
      CALL C2N2DM(DDF,DWORK(KPLF),KWORK(LLM2),KWORK(LLA2),
     *            NEL2,NAT2,NVT2,0)
C
      CALL MR031(DWORK(KPLF),DWORK(KPLC),
     *           KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *           NAT2,NAT1,NEL2,NEL1)
C
      CALL C2N2DM(DDC,DWORK(KPLC),KWORK(LLM1),KWORK(LLA1),
     *            NEL1,NAT1,NVT1,1)
C
      GOTO 99999
C
      ENDIF
C-----------------------------------------------------------------------
C
C
C
99999 CALL TOL20A(DDC,VWORK(L(KLVOL(ILEV))),NEL1,INEUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE YSMP (DX,DB,DD,NEQ,NSM)  
************************************************************************
*
*   Purpose: - performs NSM smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
C
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C
C
      IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.1) THEN
       CALL ParIA217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMP,ILEV)
       GOTO 999
      ENDIF
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.2) THEN
       IF (RLXSMP.EQ.1D0) THEN
        CALL ParIB217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),DX,DB,DD,NEQ,NSM,ILEV)
       ELSE
        CALL ParIC217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMP,ILEV)
       ENDIF
       GOTO 999
      ENDIF
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.3) THEN
       CALL ParID217(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),DX,DB,DD,NEQ,NSM,RLXSMP,ILEV)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      IF (ISMP.EQ.4) THEN
       CALL PArIF217(DWORK(L(KLC(ILEV))),VWORK(L(KLCILU(ILEV))),
     *            KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *            DX,DB,DD,NEQ,NSM,RLXSMP,ILEV)
       GOTO 999
      ENDIF
C-----------------------------------------------------------------------
C
C
999   IF (ISORTP.GT.0) THEN
       KTRC1=L(KLTRC1(ILEV))
       KTRC2=L(KLTRC2(ILEV))
       CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
       CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
       CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
      ENDIF
C
C
      CALL TOL20A(DX,VWORK(L(KLVOL(ILEV))),NEQ,INEUM)
C
C
C
99999 END
c
c
c
************************************************************************
      SUBROUTINE YSTEPP (DX,DD,DB,NEQ,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          ALFMIN and ALFMAX (COMMON /RPARM/)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C *** Setting of given ALPHA
C
      IF (AMINP.EQ.AMAXP) THEN
       ALPHA=AMINP
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
C
      IF (AMINP.NE.AMAXP) THEN
C
       IF (ISORTP.GT.0) THEN
        KTRC1=L(KLTRC1(ILEV))
        KTRC2=L(KLTRC2(ILEV))
        CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
        CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
        CALL VECSRT(DD ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,1)
        CALL LCP1  (DWORK(L(LDP)),DD ,NEQ)
       ENDIF
C
       CALL DYAX7(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *            KWORK(L(KLLDC(ILEV))),NEQ,DD,DX,DD,DYAX1,DYAX2)
       CALL LSP1(DB,DD,NEQ,DBX)
       ALPHA=(DBX-DYAX1)/DYAX2
       IF (ALPHA.LT.AMINP) ALPHA=AMINP
       IF (ALPHA.GT.AMAXP) ALPHA=AMAXP
C
       IF (ISORTP.GT.0) THEN
        KTRC1=L(KLTRC1(ILEV))
        KTRC2=L(KLTRC2(ILEV))
        CALL VECSRT(DX ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DX ,NEQ)
        CALL VECSRT(DB ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DB ,NEQ)
        CALL VECSRT(DD ,DWORK(L(LDP)),KWORK(KTRC1),KWORK(KTRC2),NEQ,2)
        CALL LCP1  (DWORK(L(LDP)),DD ,NEQ)
       ENDIF
C
       RETURN
C
      ENDIF
C
C
      END
