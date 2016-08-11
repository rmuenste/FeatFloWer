************************************************************************
      SUBROUTINE   PROJST  (KAREA,KCOLC,KLDC,KADJ,NEL,NAT,
     *                      NABD,NC,KABD,KERIOA,KERIOB,KCROSS)
************************************************************************
*    Purpose:  calculates the matrix positions for projection of C
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNAE=6)
      DIMENSION KCOLC(*),KLDC(*),KADJ(NNAE,*),KAREA(NNAE,*)
      DIMENSION KCROSS(*),KABD(*),KERIOA(*),KERIOB(NEL,2)
C
C
C
      NC=0
      KLDC(1)=1
C
      DO 10 IEL=1,NEL
C
      IELPERFOUND = 0
C
      NC=NC+1
      KCOLC(NC)=IEL
C
      DO 20 IAT=1,6
      IAREA=KAREA(IAT,IEL)
      IF (KCROSS(IAREA).EQ.0) THEN
       IADJ=KADJ(IAT,IEL)
      ELSE
       IELPERFOUND = IELPERFOUND + 1
       IF (IELPERFOUND.EQ.1) THEN
        IADJ = ABS(KERIOB(IEL,1))
       ELSE
        IADJ = ABS(KERIOB(IEL,2))
       END IF
      END IF
      IF (IADJ.EQ.0) GOTO 20
C
      NC=NC+1
      KCOLC(NC)=IADJ
20    CONTINUE
C
      KLDC(IEL+1)=NC+1
C
10    CONTINUE
C
C
      DO 30 IEQ=1,NEL
C
31    BSORT=.TRUE.
      DO 32 ICOL=KLDC(IEQ)+1,KLDC(IEQ+1)-2
      IF (KCOLC(ICOL).GT.KCOLC(ICOL+1)) THEN
       IHELP=KCOLC(ICOL)
       KCOLC(ICOL)=KCOLC(ICOL+1)
       KCOLC(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
32    CONTINUE
      IF (.NOT.BSORT) GOTO 31
C
30    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE   PROJMA  (DC,KCOLC,KLDC,KNPR,KAREA,KADJ,DM,DRHO,
     *             B1,B2,B3,KCOLB,KLDB,KABD,ILEV)
************************************************************************
*     Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)

      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      REAL*8 B1,B2,B3
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DC(*),KCOLC(*),KLDC(*)
      DIMENSION KNPR(*),KAREA(NNAE,*),KADJ(NNAE,*)
      DIMENSION DM(*),DRHO(*),B1(*),B2(*),B3(*),KCOLB(*),KLDB(*)
      DIMENSION KABD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: ParM
C

      IF (myid.NE.MASTER) THEN                      ! PARALLEL
       ALLOCATE (ParM(NAT))                         ! PARALLEL
       DO IAT=1,NAT                                 ! PARALLEL
        ParM(IAT) = DRHO(IAT)*DM(IAT)               ! PARALLEL
       END DO                                       ! PARALLEL
                                                    ! PARALLEL
       CALL CommSum(ParM,ILEV)                      ! PARALLEL
                                                    ! PARALLEL
       DO IU=1,mg_mpi(ILEV)%NeighNum                ! PARALLEL
        mg_mpi(ILEV)%parST(IU)%i=1                  ! PARALLEL
       END DO                                       ! PARALLEL
                                                    ! PARALLEL
      END IF                                        ! PARALLEL

      DO 10 IEL=1,NEL
C
      ILD1=KLDC(IEL)
      ILD2=KLDC(IEL+1)-1
C
      DO 20 IAT=1,6
      IAREA = KAREA(IAT,IEL)
      IADJ  = KADJ(IAT,IEL)
      DMM   = DRHO(IAREA)*DM(IAREA)
      INPR  = KNPR(IAREA+NVT)
      IF (INPR.EQ.1) GOTO 20
C
      KLOWER = KLDB(IAREA  )
      KUPPER = KLDB(IAREA+1)
C
      DO 30 ILDB1=KLOWER,KUPPER-1
      IF (KCOLB(ILDB1).EQ.IEL) GOTO 32
30    CONTINUE
      WRITE(*,*) "Indice not found in C-Matrix E",IEL,IAREA
      STOP
C
32    CONTINUE 
C
      DH1= B1(ILDB1)*B1(ILDB1)+B2(ILDB1)*B2(ILDB1)+
     *     B3(ILDB1)*B3(ILDB1)
C
      IF (INPR.EQ.2) THEN                                            ! PARALLEL
       DC(ILD1)=DC(ILD1)+DH1/ParM(IAREA)                             ! PARALLEL
       DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
        IUU = mg_mpi(ILEV)%parST(pID)%i                              ! PARALLEL
        IF (IUU.LE.mg_mpi(ILEV)%parST(pID)%Num) THEN                 ! PARALLEL
         IF (mg_mpi(ILEV)%parST(pID)%FaceLink(1,IUU).EQ.IAREA) THEN  ! PARALLEL
          mg_mpi(ILEV)%parST(pID)%PE(IUU)=-DH1/ParM(IAREA)           ! PARALLEL
          mg_mpi(ILEV)%parST(pID)%i=IUU+1                            ! PARALLEL
          GOTO 66                                                    ! PARALLEL
         END IF                                                      ! PARALLEL
        END IF                                                       ! PARALLEL
       END DO                                                        ! PARALLEL

       WRITE(*,*) myid,IAREA,"Problem in parallelization of projma"  ! PARALLEL
       DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
       WRITE(*,*)"--",mg_mpi(ILEV)%parST(pID)%i
        DO IAS=1,mg_mpi(ILEV)%parST(pID)%Num
         WRITE(*,*) mg_mpi(ILEV)%parST(pID)%FaceLink(1,IAS),
     * mg_mpi(ILEV)%parST(pID)%ElemLink(1,IAS)
        END DO
       END DO
       STOP                                                          ! PARALLEL
66     CONTINUE                                                      ! PARALLEL
      ELSE                                                           ! PARALLEL
       DC(ILD1)=DC(ILD1)+DH1/DMM                                     ! PARALLEL
      END IF                                                         ! PARALLEL
      IF (IADJ.EQ.0) GOTO 20
C
      DO 40 ILDB2=KLOWER,KUPPER-1
      IF (KCOLB(ILDB2).EQ.IADJ) GOTO 42
40    CONTINUE
      WRITE(*,*) "Indice not found in C-Matrix A",IEL,IAREA
      STOP
C
42    CONTINUE 
C
      DO 50 ILDC=ILD1+1,ILD2
      IF (KCOLC(ILDC).EQ.IADJ) GOTO 52
50    CONTINUE
      WRITE(*,*) 
     *"Indice not found in C-Matrix C",IEL,IADJ,ILD1,ILD2,IAREA
      STOP
52    CONTINUE
C
      DH2= B1(ILDB1)*B1(ILDB2)+B2(ILDB1)*B2(ILDB2)+
     *     B3(ILDB1)*B3(ILDB2)
      DC(ILDC)=DH2/DMM
C
20    CONTINUE
C
10    CONTINUE
C
      IF (myid.NE.MASTER) THEN                               ! PARALLEL
       DEALLOCATE (ParM)                                     ! PARALLEL
      END IF                                                 ! PARALLEL
C
!       IF (myid.eq.2) THEN
!       DO 200 IE=1,mg_mpi(ILEV)%parST(1)%Num
!        IEQ=mg_mpi(ILEV)%parST(1)%ElemLink(1,IE)
!        WRITE(*,*) "---------------",IEQ,mg_mpi(ILEV)%parST(1)%PE(IE)
!        DO 210 ILD=KLDC(IEQ),KLDC(IEQ+1)-1
!         WRITE(6,*) ILD,KCOLC(ILD),DC(ILD)
! 210    CONTINUE
!        WRITE(*,*) "---------------"
! 200   CONTINUE
!       END IF
      IF (myid.eq.-1) THEN
!       IF (MT.GT.9) THEN
      DO 200 IEQ=1,NEL
!       IF (KERIOB(IEQ,1).GT.0) THEN
       WRITE(*,*) "---------------",IEQ!,KERIOB(IEQ,1)
       DO 210 ILD=KLDC(IEQ),KLDC(IEQ+1)-1
        WRITE(6,*) ILD,KCOLC(ILD),DC(ILD)
210    CONTINUE
       WRITE(*,*) "---------------"
!       END IF
200   CONTINUE
!       ENDIF
      END IF
C
      END


************************************************************************
      SUBROUTINE   PROJMMA  (DC,KCOLC,KLDC,KNPR,KAREA,KADJ,DM,
     *             B1,B2,B3,KCOLB,KLDB,KABD,ILEV)
************************************************************************
*     Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      USE PP3D_MPI
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)

      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      REAL*8 B1,B2,B3
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DC(*),KCOLC(*),KLDC(*)
      DIMENSION KNPR(*),KAREA(NNAE,*),KADJ(NNAE,*)
      DIMENSION DM(*),B1(*),B2(*),B3(*),KCOLB(*),KLDB(*)
      DIMENSION KABD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      REAL*8 ,DIMENSION(:), ALLOCATABLE :: ParM
C

      IF (myid.NE.MASTER) THEN                      ! PARALLEL
       ALLOCATE (ParM(NAT))                         ! PARALLEL
       DO IAT=1,NAT                                 ! PARALLEL
        ParM(IAT) = DM(IAT)               ! PARALLEL
       END DO                                       ! PARALLEL
                                                    ! PARALLEL
       CALL CommSum(ParM,ILEV)                      ! PARALLEL
                                                    ! PARALLEL
       DO IU=1,mg_mpi(ILEV)%NeighNum                ! PARALLEL
        mg_mpi(ILEV)%parST(IU)%i=1                  ! PARALLEL
       END DO                                       ! PARALLEL
                                                    ! PARALLEL
      END IF                                        ! PARALLEL

      DO 10 IEL=1,NEL
C
      ILD1=KLDC(IEL)
      ILD2=KLDC(IEL+1)-1
C
      DO 20 IAT=1,6
      IAREA = KAREA(IAT,IEL)
      IADJ  = KADJ(IAT,IEL)
      DMM   = DM(IAREA)
      INPR  = KNPR(IAREA+NVT)
      IF (INPR.EQ.1) THEN
!        write(*,*) iel
!        GOTO 20
      END IF
C
      KLOWER = KLDB(IAREA  )
      KUPPER = KLDB(IAREA+1)
C
      DO 30 ILDB1=KLOWER,KUPPER-1
      IF (KCOLB(ILDB1).EQ.IEL) GOTO 32
30    CONTINUE
      WRITE(*,*) "Indice not found in C-Matrix E",IEL,IAREA
      STOP
C
32    CONTINUE 
C
      DH1= B1(ILDB1)*B1(ILDB1)+B2(ILDB1)*B2(ILDB1)+
     *     B3(ILDB1)*B3(ILDB1)
C
      IF (INPR.EQ.2) THEN                                            ! PARALLEL
       DC(ILD1)=DC(ILD1)+DH1/ParM(IAREA)                             ! PARALLEL
       DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
        IUU = mg_mpi(ILEV)%parST(pID)%i                              ! PARALLEL
        IF (IUU.LE.mg_mpi(ILEV)%parST(pID)%Num) THEN                 ! PARALLEL
         IF (mg_mpi(ILEV)%parST(pID)%FaceLink(1,IUU).EQ.IAREA) THEN  ! PARALLEL
          mg_mpi(ILEV)%parST(pID)%PE(IUU)=-DH1/ParM(IAREA)           ! PARALLEL
          mg_mpi(ILEV)%parST(pID)%i=IUU+1                            ! PARALLEL
          GOTO 66                                                    ! PARALLEL
         END IF                                                      ! PARALLEL
        END IF                                                       ! PARALLEL
       END DO                                                        ! PARALLEL

       WRITE(*,*) myid,IAREA,"Problem in parallelization of projma"  ! PARALLEL
       DO pID=1,mg_mpi(ILEV)%NeighNum                                ! PARALLEL
       WRITE(*,*)"--",mg_mpi(ILEV)%parST(pID)%i
        DO IAS=1,mg_mpi(ILEV)%parST(pID)%Num
         WRITE(*,*) mg_mpi(ILEV)%parST(pID)%FaceLink(1,IAS),
     * mg_mpi(ILEV)%parST(pID)%ElemLink(1,IAS)
        END DO
       END DO
       STOP                                                          ! PARALLEL
66     CONTINUE                                                      ! PARALLEL
      ELSE                                                           ! PARALLEL
       DC(ILD1)=DC(ILD1)+DH1/DMM                                     ! PARALLEL
      END IF                                                         ! PARALLEL
      IF (IADJ.EQ.0) GOTO 20
C
      DO 40 ILDB2=KLOWER,KUPPER-1
      IF (KCOLB(ILDB2).EQ.IADJ) GOTO 42
40    CONTINUE
      WRITE(*,*) "Indice not found in C-Matrix A",IEL,IAREA
      STOP
C
42    CONTINUE 
C
      DO 50 ILDC=ILD1+1,ILD2
      IF (KCOLC(ILDC).EQ.IADJ) GOTO 52
50    CONTINUE
      WRITE(*,*) 
     *"Indice not found in C-Matrix C",IEL,IADJ,ILD1,ILD2,IAREA
      STOP
52    CONTINUE
C
      DH2= B1(ILDB1)*B1(ILDB2)+B2(ILDB1)*B2(ILDB2)+
     *     B3(ILDB1)*B3(ILDB2)
      DC(ILDC)=DH2/DMM
C
20    CONTINUE
C
10    CONTINUE
C
      IF (myid.NE.MASTER) THEN                               ! PARALLEL
       DEALLOCATE (ParM)                                     ! PARALLEL
      END IF                                                 ! PARALLEL
C
!       IF (myid.eq.2) THEN
!       DO 200 IE=1,mg_mpi(ILEV)%parST(1)%Num
!        IEQ=mg_mpi(ILEV)%parST(1)%ElemLink(1,IE)
!        WRITE(*,*) "---------------",IEQ,mg_mpi(ILEV)%parST(1)%PE(IE)
!        DO 210 ILD=KLDC(IEQ),KLDC(IEQ+1)-1
!         WRITE(6,*) ILD,KCOLC(ILD),DC(ILD)
! 210    CONTINUE
!        WRITE(*,*) "---------------"
! 200   CONTINUE
!       END IF

!       IF (myid.eq.1) THEN
!       IF (MT.GT.9) THEN
      DO 200 IEQ=1,NEL
!       IF (KERIOB(IEQ,1).GT.0) THEN
      IF (IEQ.EQ.12265) THEN
       WRITE(*,*) "---------------",IEQ!,KERIOB(IEQ,1)
       DO 210 ILD=KLDC(IEQ),KLDC(IEQ+1)-1
        WRITE(6,*) ILD,KCOLC(ILD),DC(ILD)
210    CONTINUE
       WRITE(*,*) "---------------"
      END IF
200   CONTINUE
!       ENDIF
!       END IF
C
      END
