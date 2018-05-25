MODULE Sigma_User
USE PP3D_MPI, ONLY:myid,showid,subnodes,dZPeriodicLength
USE var_QuadScalar ,ONLY : bNoOutflow

IMPLICIT NONE

TYPE tSegment
  INTEGER :: nOFFfilesR=0,nOFFfilesL=0,nOFFfiles=0
  CHARACTER*200, ALLOCATABLE :: OFFfilesR(:),OFFfilesL(:),OFFfiles(:)
  INTEGER, ALLOCATABLE :: idxCgal(:)
  CHARACTER*20 ObjectType,Unit
  CHARACTER*10 name
  CHARACTER*99 ::  cF
  CHARACTER*8 ::  ART
  INTEGER ::    KNETz,N
  REAL*8  :: Ds,s,delta,Dss
  REAL*8, ALLOCATABLE :: Zknet(:)
  REAL*8 :: t,D,Alpha,StartAlpha ! t=Gangsteigung
  REAL*8 :: Min, Max,L
  REAL*8 :: ZME_DiscThick,ZME_gap_SG, ZME_gap_SS 
  INTEGER :: ZME_N
  REAL*8  :: SecProf_W, SecProf_D,SecProf_L
  INTEGER :: SecProf_N, SecProf_I
END TYPE tSegment
!------------------------------------------------------------
TYPE tSigma
!   REAL*8 :: Dz_out,Dz_in, a, L, Ds, s, delta,SegmentLength, DZz,W
  CHARACTER cType*(50),cZwickel*(50)
  REAL*8 :: Dz_out,Dz_in, a, L, SegmentLength, DZz,W
  REAL*8 :: SecStr_W,SecStr_D
  INTEGER :: NumberOfSeg, GANGZAHL,STLSeg=0
  TYPE (tSegment), ALLOCATABLE :: mySegment(:)
END TYPE tSigma

TYPE(tSigma) :: mySigma
!------------------------------------------------------------
TYPE tProcess
   REAL*8 :: Umdr, Ta, Ti, T0, Massestrom, Dreh, Angle, dPress,Phase
   REAL*8 :: MinInflowDiameter,MaxInflowDiameter
   INTEGER :: Periodicity
   REAL*8 :: dAlpha
   CHARACTER*6 :: Rotation !RECHT, LINKS
   CHARACTER*50 :: pTYPE !RECHT, LINKS
   INTEGER :: ind
END TYPE tProcess

TYPE(tProcess) :: myProcess
!------------------------------------------------------------
TYPE tRheology
   INTEGER :: Equation = 5
   INTEGER :: AtFunc
   REAL*8 :: A, B, C ! Carreau Parameter
   REAL*8 :: n, K ! Power Law
   REAL*8 :: eta_max, eta_min 
   REAL*8 :: Ts, Tb, C1, C2 ! WLF Parameter
   REAL*8 :: ViscoMin = 1e-4
   REAL*8 :: ViscoMax = 1e10
END TYPE tRheology

TYPE(tRheology) :: myRheology

!------------------------------------------------------------
TYPE tThermodyn
   REAL*8 :: density, densitySteig
   REAL*8 :: lambda, Cp, lambdaSteig,CpSteig 
   REAL*8 :: Alpha, Beta, Gamma
END TYPE tThermodyn

TYPE(tThermodyn) :: myThermodyn
!------------------------------------------------------------

TYPE tSetup
 CHARACTER*200 cMeshPath
 CHARACTER*20 cMesher
 INTEGER MeshResolution,nSolutions
 INTEGER m_nT,m_nR,m_nZ,m_nP
 LOGICAL :: bGeoTest=.FALSE.,bSendEmail=.TRUE.
END TYPE tSetup
TYPE(tSetup) :: mySetup

TYPE tOutput
 INTEGER nOf1DLayers,nRegion
END TYPE tOutput

TYPE(tOutput) :: myOutput

LOGICAL bCtrl
REAL*8 DistTolerance
LOGICAL :: bKTPRelease = .TRUE.

CONTAINS


! !-------------------------------------------------------------
! SUBROUTINE GetSigmaParameters()
! IMPLICIT NONE
! INTEGER :: myFile=377
! INTEGER iEnd,iAt,iEq,iStrEnd
! real*8 daux,myPI
! CHARACTER*8 netz
! INTEGER ik,k,iSeg,iaux,i
! CHARACTER string*100,param*25,tVar*10,cPar*50, cName*7,caux*50,cType*50
! LOGICAL bOK
! 
! ALLOCATE (mySigma%mySegment(10))
! 
! OPEN (UNIT=myFile,FILE="_data/Extrud3D.dat")
! myRheology%AtFunc = 1
! 
! mySigma%NumberOfSeg=1
! 
! READ (myFile,*) myProcess%Angle
! 
! DO
!    READ (UNIT=myFile,FMT='(A100)' ,IOSTAT=iEnd) string
!    IF (iEnd.EQ.-1) EXIT
!    CALL StrStuct()
!    IF (bOK) THEN
!    
! 
! 
!       READ(string(1:iAt-1),*) cType
!       READ(string(iAt+1:iEq-1),*) cPar
! !       IF (myid.eq.showid) write(*,'(10A)') "'",TRIM(ADJUSTL(cType)),"'","-'",TRIM(ADJUSTL(cPar)),"'"
! 
!       SELECT CASE (TRIM(ADJUSTL(cType)))
! !---------------------------------------------------------- 
!       CASE ("GeoTest")       
!         bGeoTest = .true.
!         DO i=1,5
!          IF (myid.eq.showid) write(*,*) 
!         END DO
!         IF (myid.eq.showid) write(*,*) "!!!!              warning              !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!                                   !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!                                   !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!      Only geometrical testrun     !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!                                   !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!                                   !!!! "
!         IF (myid.eq.showid) write(*,*) "!!!!              warning              !!!! "
!         DO i=1,5
!          IF (myid.eq.showid) write(*,*) 
!         END DO
! 
!       CASE ("Maschinendaten")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Maschine_Standard")
!           READ(string(iEq+1:iStrEnd-1),*) myProcess%Rotation
!           IF (myid.eq.showid) write(*,*) "Drehrichtung= ", myProcess%Rotation
!           IF (TRIM(myProcess%Rotation).eq."RECHTS") myProcess%ind=1
!           IF (TRIM(myProcess%Rotation).eq."LINKS") myProcess%ind=-1
!         END SELECT
! !---------------------------------------------------------- 
!       CASE ("Gehaeusedaten")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Zylinder_Standard")
!            READ(string(iEq+1:),*) mySigma%Dz,mySigma%a,mySigma%L
!            IF (myid.eq.showid) write(*,*) "Gehaeusedaten= ", cPar
!            IF (myid.eq.showid) write(*,*) "Dz= ",mySigma%Dz
!            IF (myid.eq.showid) write(*,*) "a= ",mySigma%a
!            IF (myid.eq.showid) write(*,*) "L= ",mySigma%L
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Elementdaten")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("0_Gang_Standard")
!            mySigma%GANGZAHL= 0
!            IF (myid.eq.showid) write(*,*) "Elementdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Gang= ",mySigma%GANGZAHL
!            READ(string(iEq+1:),*) mySigma%Ds
!            IF (myid.eq.showid) write(*,*) "InnerRadius= ",mySigma%Ds
!         CASE ("1_Gang_Standard")
!            mySigma%GANGZAHL= 1
!            IF (myid.eq.showid) write(*,*) "Elementdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Gang= ",mySigma%GANGZAHL
!            READ(string(iEq+1:),*) mySigma%Ds, mySigma%s
!            IF (myid.eq.showid) write(*,*) "Ds= ",mySigma%Ds
!            IF (myid.eq.showid) write(*,*) "Spiel SS= ", mySigma%s
!            mySigma%delta=(mySigma%Dz - mySigma%Ds)/2d0
!            IF (myid.eq.showid) write(*,*) "Spiel SG= ", mySigma%delta
!         CASE ("2_Gang_Standard")
!            mySigma%GANGZAHL= 2
!            IF (myid.eq.showid) write(*,*) "Elementdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Gang= ",mySigma%GANGZAHL
!            READ(string(iEq+1:),*) mySigma%Ds, mySigma%s
!            IF (myid.eq.showid) write(*,*) "Ds= ",mySigma%Ds
!            IF (myid.eq.showid) write(*,*) "Spiel SS= ", mySigma%s
!            mySigma%delta=(mySigma%Dz - mySigma%Ds)/2d0
!            IF (myid.eq.showid) write(*,*) "Spiel SG= ", mySigma%delta
!         CASE ("3_Gang_Standard")
!            mySigma%GANGZAHL= 3
!            IF (myid.eq.showid) write(*,*) "Elementdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Gang= ",mySigma%GANGZAHL
!            READ(string(iEq+1:),*) mySigma%Ds, mySigma%s
!            IF (myid.eq.showid) write(*,*) "Ds= ",mySigma%Ds
!            IF (myid.eq.showid) write(*,*) "Spiel SS= ", mySigma%s
!            mySigma%delta=(mySigma%Dz - mySigma%Ds)/2d0
!            IF (myid.eq.showid) write(*,*) "Spiel SG= ", mySigma%delta
!         CASE ("4_Gang_Standard")
!            mySigma%GANGZAHL= 4
!            IF (myid.eq.showid) write(*,*) "Elementdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Gang= ",mySigma%GANGZAHL
!            READ(string(iEq+1:),*) mySigma%Ds, mySigma%s
!            IF (myid.eq.showid) write(*,*) "Ds= ",mySigma%Ds
!            IF (myid.eq.showid) write(*,*) "Spiel SS= ", mySigma%s
!            mySigma%delta=(mySigma%Dz - mySigma%Ds)/2d0
!            IF (myid.eq.showid) write(*,*) "Spiel SG= ", mySigma%delta
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Foerderart1")
!         iSeg = 1
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("ZME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%ZME_N,mySigma%mySegment(iSeg)%ZME_DiscThick,&
!            mySigma%mySegment(iSeg)%ZME_gap_SS,mySigma%mySegment(iSeg)%ZME_gap_SG,mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + &
!            2d0*mySigma%mySegment(iSeg)%ZME_N * (mySigma%mySegment(iSeg)%ZME_DiscThick + mySigma%mySegment(iSeg)%ZME_gap_SS)
! 
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
! 
!            mySigma%mySegment(iSeg)%ART= "ZME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenanzahl= ", mySigma%mySegment(iSeg)%ZME_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenbreite= ", mySigma%mySegment(iSeg)%ZME_DiscThick
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SG = ", mySigma%mySegment(iSeg)%ZME_gap_SG
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SS = ", mySigma%mySegment(iSeg)%ZME_gap_SS
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("SME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "SME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("CONE_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%L, mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "CONE"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("FOERD_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "FOERD"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("KNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "KNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("SKNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "SKNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         END SELECT
! 
!       CASE ("Foerderart2")
!         iSeg = 2
!         mySigma%NumberOfSeg=2
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("ZME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%ZME_N,mySigma%mySegment(iSeg)%ZME_DiscThick,&
!            mySigma%mySegment(iSeg)%ZME_gap_SS,mySigma%mySegment(iSeg)%ZME_gap_SG,mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + &
!            2d0*mySigma%mySegment(iSeg)%ZME_N * (mySigma%mySegment(iSeg)%ZME_DiscThick + mySigma%mySegment(iSeg)%ZME_gap_SS)
! 
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
! 
!            mySigma%mySegment(iSeg)%ART= "ZME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenanzahl= ", mySigma%mySegment(iSeg)%ZME_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenbreite= ", mySigma%mySegment(iSeg)%ZME_DiscThick
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SG = ", mySigma%mySegment(iSeg)%ZME_gap_SG
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SS = ", mySigma%mySegment(iSeg)%ZME_gap_SS
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("SME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "SME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("CONE_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%L, mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "CONE"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("FOERD_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "FOERD"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("KNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "KNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("SKNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "SKNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         END SELECT
! 
!       CASE ("Foerderart3")
!         iSeg = 3
!         mySigma%NumberOfSeg=3
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("ZME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%ZME_N,mySigma%mySegment(iSeg)%ZME_DiscThick,&
!            mySigma%mySegment(iSeg)%ZME_gap_SS,mySigma%mySegment(iSeg)%ZME_gap_SG,mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + &
!            2d0*mySigma%mySegment(iSeg)%ZME_N * (mySigma%mySegment(iSeg)%ZME_DiscThick + mySigma%mySegment(iSeg)%ZME_gap_SS)
! 
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
! 
!            mySigma%mySegment(iSeg)%ART= "ZME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenanzahl= ", mySigma%mySegment(iSeg)%ZME_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Scheibenbreite= ", mySigma%mySegment(iSeg)%ZME_DiscThick
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SG = ", mySigma%mySegment(iSeg)%ZME_gap_SG
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Spiel SS = ", mySigma%mySegment(iSeg)%ZME_gap_SS
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("SME_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min,&
!            mySigma%mySegment(iSeg)%SecProf_N,mySigma%mySegment(iSeg)%SecProf_I,&
!            mySigma%mySegment(iSeg)%SecProf_D,mySigma%mySegment(iSeg)%SecProf_W,mySigma%mySegment(iSeg)%SecProf_L
! 
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "SME"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Num= "  , mySigma%mySegment(iSeg)%SecProf_N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Type= " , mySigma%mySegment(iSeg)%SecProf_I
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Width= ", mySigma%mySegment(iSeg)%SecProf_W
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Depth= ", mySigma%mySegment(iSeg)%SecProf_D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - CutOff-Pitch= ", mySigma%mySegment(iSeg)%SecProf_L
!         CASE ("CONE_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%L, mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "CONE"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("FOERD_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%t, mySigma%mySegment(iSeg)%L, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
!            mySigma%mySegment(iSeg)%ART= "FOERD"
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Steigung= ",mySigma%mySegment(iSeg)%t
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("KNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "KNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         CASE ("SKNET_Standard")
!            READ(string(iEq+1:),*) mySigma%mySegment(iSeg)%N,mySigma%mySegment(iSeg)%Alpha, mySigma%mySegment(iSeg)%D, &
!                 mySigma%mySegment(iSeg)%Min
!            mySigma%mySegment(iSeg)%ART= "SKNET"
!            mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
!            mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Foerderart= ",TRIM(ADJUSTL(cPar))
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Art= ",mySigma%mySegment(iSeg)%ART
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Versatzwinkel= ",mySigma%mySegment(iSeg)%Alpha
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenbreite= ", mySigma%mySegment(iSeg)%D
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Knetscheibenanzahl= ", mySigma%mySegment(iSeg)%N
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Laenge des Segmentes= ", mySigma%mySegment(iSeg)%L
!            IF (myid.eq.showid) write(*,*) "Segment ",iSeg," - Startposition= ", mySigma%mySegment(iSeg)%Min
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Verfahrensdaten")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Verfahren_Periodic")
!            READ(string(iEq+1:),*) myProcess%Umdr,myProcess%Ta,myProcess%Ti,myProcess%T0,&
!                 myProcess%dPress 
!            IF (myid.eq.showid) write(*,*) "Verfahrensdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Drehzahl= ", myProcess%Umdr 
!            IF (myProcess%Ta.LT.0d0) myProcess%Ta = 0d0/0d0
!            IF (myProcess%Ti.LT.0d0) myProcess%Ti = 0d0/0d0
! !            IF (myid.eq.1) WRITE(*,*) myProcess%Umdr,myProcess%Ta,myProcess%Ti,myProcess%T0,&
! !                 myProcess%Massestrom 
!            IF (myid.eq.showid) THEN
!               IF (isnan(myProcess%Ta)) THEN
!                  write(*,*) "Wand ist adiabat "
!               ELSE 
!                  write(*,*) "Wandtemperatur=  ", myProcess%Ta
!               END IF
!            END IF
!            IF (myid.eq.showid) THEN
!               IF (isnan(myProcess%Ti)) THEN
!                  write(*,*) "Schnecke ist adiabat "
!               ELSE 
!                  write(*,*) "Schneckentemperatur= ", myProcess%Ti
!               END IF
!            END IF
!            IF (myid.eq.showid) write(*,*)  "T0= ", myProcess%T0
!            IF (myid.eq.showid) write(*,*)  "DruckDiff= ", myProcess%dPress
!            myProcess%Massestrom = 0d0/0d0
! 
!         CASE ("Verfahren_Standard")
!            READ(string(iEq+1:),*) myProcess%Umdr,myProcess%Ta,myProcess%Ti,myProcess%T0,&
!                 myProcess%Massestrom 
!            IF (myid.eq.showid) write(*,*) "Verfahrensdaten= ",cPar
!            IF (myid.eq.showid) write(*,*) "Drehzahl= ", myProcess%Umdr 
!            IF (myProcess%Ta.LT.0d0) myProcess%Ta = 0d0/0d0
!            IF (myProcess%Ti.LT.0d0) myProcess%Ti = 0d0/0d0
! !            IF (myid.eq.1) WRITE(*,*) myProcess%Umdr,myProcess%Ta,myProcess%Ti,myProcess%T0,&
! !                 myProcess%Massestrom 
!            IF (myid.eq.showid) THEN
!               IF (isnan(myProcess%Ta)) THEN
!                  write(*,*) "Wand ist adiabat "
!               ELSE 
!                  write(*,*) "Wandtemperatur=  ", myProcess%Ta
!               END IF
!            END IF
!            IF (myid.eq.showid) THEN
!               IF (isnan(myProcess%Ti)) THEN
!                  write(*,*) "Schnecke ist adiabat "
!               ELSE 
!                  write(*,*) "Schneckentemperatur= ", myProcess%Ti
!               END IF
!            END IF
!            IF (myid.eq.showid) write(*,*)  "T0= ", myProcess%T0
!            IF (myid.eq.showid) write(*,*)  "Durchsatz= ", myProcess%Massestrom 
!            myProcess%dPress = 0d0/0d0
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Berechnungsart_Visko")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Carreau") 
!            IF (myid.eq.showid) write(*,*) "Berechnungsart_Visko= ",cPar
!            myRheology%Equation = 1
!            READ(string(iEq+1:),*) myRheology%A, myRheology%B, myRheology%C
!            IF (myid.eq.showid) write(*,*) "Carreau A= ", myRheology%A
!            IF (myid.eq.showid) write(*,*) "Carreau B= ", myRheology%B
!            IF (myid.eq.showid) write(*,*) "Carreau C= ", myRheology%C
!         CASE ("Potenz")
!            IF (myid.eq.showid) write(*,*) "Berechnungsart_Visko= ",cPar
!            myRheology%Equation = 2
!            READ(string(iEq+1:),*) myRheology%n, myRheology%K
!            IF (myid.eq.showid) write(*,*) "PowerLaw n= ", myRheology%n
!            IF (myid.eq.showid) write(*,*) "PowerLaw K= ", myRheology%K
!         CASE ("Polyflow")
!            IF (myid.eq.showid) write(*,*) "Berechnungsart_Visko= ",cPar
!            myRheology%Equation = 3
!            READ(string(iEq+1:),*) myRheology%A, myRheology%B, myRheology%C
!            IF (myid.eq.showid) write(*,*) "Polyflow A= ", myRheology%A
!            IF (myid.eq.showid) write(*,*) "Polyflow B= ", myRheology%B
!            IF (myid.eq.showid) write(*,*) "Polyflow C= ", myRheology%C
!         CASE ("Ellis") 
!            IF (myid.eq.showid) write(*,*) "Berechnungsart_Visko= ",cPar
!            myRheology%Equation = 4
!            READ(string(iEq+1:),*) myRheology%A, myRheology%B, myRheology%C
!            IF (myid.eq.showid) write(*,*) "Ellis nu_0= ", myRheology%A
!            IF (myid.eq.showid) write(*,*) "Ellis gamma_0= ", myRheology%B
!            IF (myid.eq.showid) write(*,*) "Ellis n= ", myRheology%C
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Grenzen_Visko")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Grenzen")
!           READ(string(iEq+1:),*) myRheology%eta_min, myRheology%eta_max
!           IF (myid.eq.showid) write(*,*) "min Viscosoty= ", myRheology%eta_min
!           IF (myid.eq.showid) write(*,*) "max Viscosity= ", myRheology%eta_max 
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Temperaturansatz")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("WLF_TBTS")
!            IF (myid.eq.showid) write(*,*) "Temperaturansatz= ",cPar
!            myRheology%AtFunc = 3
!            READ(string(iEq+1:),*) myRheology%Ts, myRheology%Tb, myRheology%C1, myRheology%C2
!            IF (myid.eq.showid) write(*,*) "WLF-Standardtemperatur Ts= ", myRheology%Ts
!            IF (myid.eq.showid) write(*,*) "WLF-Bezugstemperatur Tb= ", myRheology%Tb
!            IF (myid.eq.showid) write(*,*) "WLF 1.Konstante C1= ", myRheology%C1
!            IF (myid.eq.showid) write(*,*) "WLF 2.Konstante C2= ", myRheology%C2
!         CASE ("WLF_C1C2")
!            IF (myid.eq.showid) write(*,*) "Temperaturansatz= ",cPar
!            myRheology%AtFunc = 2
!            READ(string(iEq+1:),*) myRheology%C1, myRheology%C2
!            myRheology%Ts = 165d0
!            IF (myid.eq.showid) write(*,*) "WLF 1.reg.Konstante C1= ", myRheology%C1
!            IF (myid.eq.showid) write(*,*) "WLF 2.reg.Konstante C2= ", myRheology%C2
!            IF (myid.eq.showid) write(*,*) "WLF-Standardtemperatur Ts= ", myRheology%Ts
!         CASE ("ISOTHERM")
!            IF (myid.eq.showid) write(*,*) "Temperaturansatz= ",cPar
!            myRheology%AtFunc = 1
!            IF (myid.eq.showid) write(*,*) "TEMPERATURKOEFF=1.0 KONSTANT "
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Thermodyn_Waerme")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Waerme_Standard")
!            IF (myid.eq.showid) write(*,*) "Thermodyn_Waerme= ",cPar
!            READ(string(iEq+1:),*) myThermodyn%lambda,myThermodyn%lambdaSteig,  myThermodyn%Cp, myThermodyn%CpSteig
!            IF (myid.eq.showid) write(*,*) "Waermeleitfaehigkei lambda= ", myThermodyn%lambda
!            IF (myid.eq.showid) write(*,*) "Steigung lambda= ", myThermodyn%lambdaSteig
!            IF (myid.eq.showid) write(*,*) "Waermekapazitaet Cp= ", myThermodyn%Cp  
!            IF (myid.eq.showid) write(*,*) "Steigung Cp= ", myThermodyn%CpSteig
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Thermodyn_Dichte")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Dichte_Standard")
!            IF (myid.eq.showid) write(*,*) "Thermodyn_Dichte= ",cPar
!            READ(string(iEq+1:),*) myThermodyn%density, myThermodyn%densitySteig
!            IF (myid.eq.showid) write(*,*) "Density= ", myThermodyn%density
!            IF (myid.eq.showid) write(*,*) "Steigung= ", myThermodyn%densitySteig
!         CASE ("Dichte_Volumen")
!            IF (myid.eq.showid) write(*,*) "Thermodyn_Dichte= ",cPar
!            READ(string(iEq+1:),*) myThermodyn%density, myThermodyn%densitySteig
!            myThermodyn%density = 1d0/myThermodyn%density
!            IF (myid.eq.showid) write(*,*) "Density= ", myThermodyn%density
!            myThermodyn%densitySteig = -myThermodyn%densitySteig
!            IF (myid.eq.showid) write(*,*) "Steigung Density= ", myThermodyn%densitySteig
!         END SELECT
! 
! !---------------------------------------------------------- 
!       CASE ("Einstellungen")
! 
!         SELECT CASE (TRIM(ADJUSTL(cPar)))
!         CASE ("Simulation_Standard")
!            IF (myid.eq.showid) write(*,*) "Einstellungen= ",cPar
!            READ(string(iEq+1:),*)  netz , myProcess%Dreh
!            IF (myid.eq.showid) write(*,*) "Netzqualitaet= ", netz
!            IF (myid.eq.showid) write(*,*) "Number of rotations = ", myProcess%Dreh
!            IF (netz.EQ.'fein')   MeshResolution = 1
!            IF (netz.EQ.'mittel') MeshResolution = 2
!            IF (netz.EQ.'grob')  MeshResolution = 3
!         END SELECT
! !---------------------------------------------------------- 
! !       CASE ("RTD")
! ! 
! !         SELECT CASE (TRIM(ADJUSTL(cPar)))
! !         CASE ("RTD_Standard")
! !            READ(string(iEq+1:),*) myRTD%nTimeLevels,myRTD%nParticles,myRTD%nRotation,myRTD%minFrac
! !            IF (myid.eq.showid) write(*,*) "RTD parameters= ",cPar
! !            IF (myid.eq.showid) write(*,*) "Zeitebene= ", myRTD%nTimeLevels
! !            IF (myid.eq.showid) write(*,*) "Partikels= ", myRTD%nParticles
! !            IF (myid.eq.showid) write(*,*) "Drehungen= ", myRTD%nRotation
! !            IF (myid.eq.showid) write(*,*) "MinFraktion= ",myRTD%minFrac
! !         END SELECT
! !---------------------------------------------------------- 
!       END SELECT
!    end if
! END DO
!    
! myPI = dATAN(1d0)*4d0
! 
! IF (myid.eq.showid) write(*,*) 'Number Of Segments = ',mySigma%NumberOfSeg
! mySigma%mySegment(1)%StartAlpha = 0.0d0
! 
! mySigma%SegmentLength = mySigma%mySegment(1)%L
! IF (mySigma%NumberOfSeg.GE.2) THEN
!  DO iSeg=2,mySigma%NumberOfSeg
!   IF (mySigma%mySegment(iSeg-1)%ART.EQ.'FOERD'.OR.mySigma%mySegment(iSeg-1)%ART.EQ.'SME') THEN
!    mySigma%mySegment(iSeg)%StartAlpha = mySigma%mySegment(iSeg-1)%StartAlpha + mySigma%mySegment(iSeg-1)%L/mySigma%mySegment(iSeg-1)%t*2d0*myPI
!   END IF
!   IF (mySigma%mySegment(iSeg-1)%ART.EQ.'KNET'.or.mySigma%mySegment(iSeg-1)%ART.EQ.'SKNET') THEN
!    mySigma%mySegment(iSeg)%StartAlpha =  mySigma%mySegment(iSeg-1)%StartAlpha + myPI*DBLE(mySigma%mySegment(iSeg-1)%N-1)*mySigma%mySegment(iSeg-1)%Alpha/1.8d2
!    IF (mySigma%mySegment(iSeg)%ART.EQ.'KNET'.or.mySigma%mySegment(iSeg)%ART.EQ.'SKNET') THEN
!     mySigma%mySegment(iSeg)%StartAlpha = mySigma%mySegment(iSeg-1)%StartAlpha + mySigma%mySegment(iSeg)%StartAlpha + myPI*mySigma%mySegment(iSeg)%Alpha/1.8d2
!    END IF
! !  WRITE(*,*) 'alpha: ',mySigma%mySegment(2)%StartAlpha
!   END IF
!   mySigma%SegmentLength = mySigma%SegmentLength + mySigma%mySegment(iSeg)%L
!  END DO
! END IF
! 
! CLOSE (myFile)
! 
! IF (.NOT.ISNAN(myProcess%dPress)) THEN
!  dZPeriodicLenght = mySigma%L
!  myProcess%dPress = 1d3*myProcess%dPress
! ELSE
!  dZPeriodicLenght = 1d5*mySigma%L
! END IF
! 
! IF (myid.eq.showid) write(*,*) 'periodic length = ', dZPeriodicLenght,myProcess%dPress
! CONTAINS!----------------------------------------------
! 
! SUBROUTINE StrStuct()
! IMPLICIT NONE
! INTEGER i,n
! 
! bOk=.FALSE.
! 
! n = len(string)
! iAt = 0
! iEq = 0
! iStrEnd = 0
! DO i=1,n
!  IF (string(i:i).EQ. '=') iAt = i
!  IF (string(i:i).EQ. '(') iEq = i
!  IF (string(i:i).EQ. ')') iStrEnd = i
! END DO
! 
! bOk=.TRUE.
! IF (iEq.eq.0)       bOk=.FALSE.
! IF (iAt.eq.0)       bOk=.FALSE.
! IF (iStrEnd.eq.0)   bOk=.FALSE.
! IF (iAt.gt.iEq)     bOk=.FALSE.
! IF (iAt.gt.iStrEnd) bOk=.FALSE.
! !WRITE(*,*) string, iEq,iAt, iStrEnd,bOk
! 
! END SUBROUTINE StrStuct
! 
! END SUBROUTINE GetSigmaParameters
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE DistanceToScrewShell(X,Y,Z,d1,d2,t)
! REAL*8 X,Y,Z,t,d1,d2
! REAL*8 dSeg1,dSeg2,tt
! INTEGER k,iaux
! 
! d1 = 1d2
! d2 = 1d2
! ! RETURN
! 
! tt = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)
! 
! !----------------------------------------------------------
! DO k=1, mySigma%NumberOfSeg
!  IF (mySigma%mySegment(k)%ART.EQ.'KNET' ) CALL KNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'SKNET') CALL SKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'FOERD') CALL FOERD_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'CONE' ) CALL CONE_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'SME'  ) CALL SME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  IF (mySigma%mySegment(k)%ART.EQ.'ZME'  ) CALL ZME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,iaux)
!  d1 = MIN(dSeg1,d1)
!  d2 = MIN(dSeg2,d2)
! END DO
! !-------------------------------------------------------
! 
! ! CALL Shell_dist(X,Y,Z,d3)
! 
! return
! 
! END SUBROUTINE DistanceToScrewShell
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE GetMixerKnpr(X,Y,Z,iBndr,inpr,D,t)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t,d
! REAL*8 :: dKnetMin,dKnetMax
! INTEGER :: inpr,iBndr, l, k, nmbr,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 dScale, dist1, dist2,dEps,dCut
! REAL*8 dSeg1,dSeg2,tt
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! tt = (myProcess%Angle/360d0)/(myProcess%Umdr/60d0)
! d = 1d2
! inpr = 0
! ! RETURN
! ! return
! 
! !----------------------------------------------------------
! DO k=1, mySigma%NumberOfSeg
!  IF (mySigma%mySegment(k)%ART.EQ.'KNET' ) CALL KNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  IF (mySigma%mySegment(k)%ART.EQ.'SKNET') CALL SKNET_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  IF (mySigma%mySegment(k)%ART.EQ.'CONE' ) CALL CONE_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  IF (mySigma%mySegment(k)%ART.EQ.'FOERD') CALL FOERD_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  IF (mySigma%mySegment(k)%ART.EQ.'SME'  ) CALL SME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  IF (mySigma%mySegment(k)%ART.EQ.'ZME'  ) CALL ZME_elem(X,Y,Z,tt,k,dSeg1,dSeg2,inpr)
!  d = MIN(dSeg1,dSeg2,d)
! END DO
! !-------------------------------------------------------
! 
! return
! 
! END SUBROUTINE GetMixerKnpr
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE Shell_dist(X,Y,Z,d)
! REAL*8 X,Y,Z,d
! REAL*8 PX1,PX2,dZ_max,dd1,dd2,dd3
! 
! dZ_max = 0.5d0*mySigma%Dz
! 
! IF (Y.GT.0d0) THEN
! 
! PX1 = +SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
! PX2 = -SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
! 
! dd1 = dZ_max - SQRT(X*X + (Y-mySigma%a/2d0)*(Y-mySigma%a/2d0))
! dd2 =          SQRT((X-PX1)*(X-PX1) + Y*Y)
! dd3 =          SQRT((X-PX2)*(X-PX2) + Y*Y)
! 
! d = min(dd1,dd2,dd3)
! ELSE
! 
! PX1 = +SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
! PX2 = -SQRT(dZ_max*dZ_max - (mySigma%a/2d0)*(mySigma%a/2d0))
! 
! dd1 = dZ_max - SQRT(X*X + (Y+mySigma%a/2d0)*(Y+mySigma%a/2d0))
! dd2 =          SQRT((X-PX1)*(X-PX1) + Y*Y)
! dd3 =          SQRT((X-PX2)*(X-PX2) + Y*Y)
! 
! d = min(dd1,dd2,dd3)
! 
! END IF
! 
! END SUBROUTINE Shell_dist
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE PipeShell_dist(X,Y,Z,d)
! REAL*8 X,Y,Z,d
! REAL*8 PX1,PX2,dZ_max,dd1,dd2,dd3
! 
! dZ_max = 0.5d0*mySigma%Dz
! 
! d = dZ_max - SQRT(X*X + Y*Y)
! 
! END SUBROUTINE PipeShell_dist
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE ZME_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 :: ZMIN,ZMAX
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut,InnerDist,dInnerRadius
! REAL*8 myPI
! REAL*8 :: ZME_SegThick
! 
! dEps = mySigma%a/1d5
! 
! IF (mySigma%GANGZAHL.NE.0) THEN
!  dInnerRadius = 1d0*(mySigma%A - 0.5d0*mySigma%Ds - 1d0*mySigma%s)
! ELSE
!  dInnerRadius = 1d0*0.5d0*mySigma%Ds
! END IF
! 
! ZME_SegThick = 2d0*(mySigma%mySegment(iSeg)%ZME_DiscThick+mySigma%mySegment(iSeg)%ZME_gap_SS)
! 
! myPI = dATAN(1d0)*4d0
! 
! IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
!  DO l=1, mySigma%mySegment(iSeg)%ZME_N
!   dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick-dEps
!   dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*ZME_SegThick+dEps
!   IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
!    lKnet = l
!    EXIT
!   END IF
!  END DO
! ELSE
!  IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
!  IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%ZME_N
! END IF
! 
! ! First screw
! dist1 = 5d0
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = 0d0 - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! daux = -5d0
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%ZME_N,lKnet+1)
! ! l=1
! 
!  
! ZMIN = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 0.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 0d0*mySigma%mySegment(iSeg)%ZME_DiscThick
! ZMAX = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 0.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 1d0*mySigma%mySegment(iSeg)%ZME_DiscThick
! 
! CALL ZME_disc(XT,YT,ZT,ZMIN,ZMAX,iSeg,daux) 
! InnerDist = sqrt(xt*xt+yt*yt)-dInnerRadius
! ! daux = sqrt(xt*xt+yt*yt) - 37d0
! dist1 =MIN(daux,dist1,InnerDist)
! 
! END DO
! 
! CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
! dist1 = max(dist1,-dCut)
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! Second screw
! dist2 = 5d0
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = 0d0 - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! daux = -5d0
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%ZME_N,lKnet+1)
! ! l=1
! 
! ZMIN = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 1.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 1d0*mySigma%mySegment(iSeg)%ZME_DiscThick
! ZMAX = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*ZME_SegThick + 1.5d0*mySigma%mySegment(iSeg)%ZME_gap_SS + 2d0*mySigma%mySegment(iSeg)%ZME_DiscThick
! 
! CALL ZME_disc(XT,YT,ZT,ZMIN,ZMAX,iSeg,daux) 
! InnerDist = sqrt(xt*xt+yt*yt)-dInnerRadius
! 
! dist2 =MIN(daux,dist2,InnerDist)
! 
! END DO
! 
! CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
! dist2 = max(dist2,-dCut)
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE ZME_elem
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE SKNET_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax,dKnetMed
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta1,dBeta2,XTT,YTT,ZTT,dist
! REAL*8 :: dZ
! REAL*8 daux1,daux2
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
!  DO l=1, mySigma%mySegment(iSeg)%N
!   dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
!   dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
!   IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
!    lKnet = l
!    EXIT
!   END IF
!  END DO
! ELSE
!  IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
!  IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%N
! END IF
! 
! ! First screw
! dist1 = 5d0
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)
! 
! dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
! dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
! dKnetMed = 0.5d0*(dKnetMin+dKnetMax)
! dBeta1 = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! dBeta2 = myPI*DBLE(l-2)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta1) - YT*sin(dBeta1)
! YTT = XT*sin(dBeta1) + YT*cos(dBeta1)
! 
! CALL Schnecke_nGang(XTT,YTT,daux1)
! 
! XTT = XT*cos(dBeta2) - YT*sin(dBeta2)
! YTT = XT*sin(dBeta2) + YT*cos(dBeta2)
! 
! CALL Schnecke_nGang(XTT,YTT,daux2)
! 
! IF (ZT.GE.dKnetMed.and.ZT.LE.dKnetMax) THEN
!  daux=daux1
! END IF
! IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMed) THEN
!  daux=MAX(daux1,daux2)
!  dZ = ABS(ZT-dKnetMed)
!  IF (daux.le.0d0) THEN
!   daux=daux
!  ELSE
!   IF (daux1.lt.0d0) THEN
!    daux=min(dZ,daux)
!   ELSE
!    daux=min(daux,SQRT(daux1**2d0+dZ**2d0))
!   END IF
!  END IF
! 
! END IF
! 
! IF (ZT.GT.dKnetMax) THEN
!  dZ = ABS(ZT-dKnetMax)
!  daux=daux1
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! IF (ZT.LT.dKnetMin) THEN
!  dZ = ABS(ZT-dKnetMin)
!  daux=MAX(daux1,daux2)
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist1 =MIN(daux,dist1)
! 
! END DO
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! 
! ! Second screw
! dist2 = 5d0
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
! 
! ! First the point needs to be transformed back to time = 0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)
! 
! dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
! dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
! dKnetMed = 0.5d0*(dKnetMin+dKnetMax)
! dBeta1 = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! dBeta2 = myPI*DBLE(l-2)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta1) - YT*sin(dBeta1)
! YTT = XT*sin(dBeta1) + YT*cos(dBeta1)
! 
! CALL Schnecke_nGang(XTT,YTT,daux1)
! 
! XTT = XT*cos(dBeta2) - YT*sin(dBeta2)
! YTT = XT*sin(dBeta2) + YT*cos(dBeta2)
! 
! CALL Schnecke_nGang(XTT,YTT,daux2)
! 
! IF (ZT.GE.dKnetMed.and.ZT.LE.dKnetMax) THEN
!  daux=daux1
! END IF
! IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMed) THEN
!  daux=MAX(daux1,daux2)
!  dZ = ABS(ZT-dKnetMed)
!  IF (daux.le.0d0) THEN
!   daux=daux
!  ELSE
!   IF (daux1.lt.0d0) THEN
!    daux=min(dZ,daux)
!   ELSE
!    daux=min(daux,SQRT(daux1**2d0+dZ**2d0))
!   END IF
!  END IF
! 
! END IF
! 
! IF (ZT.GT.dKnetMax) THEN
!  dZ = ABS(ZT-dKnetMax)
!  daux=daux1
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! IF (ZT.LT.dKnetMin) THEN
!  dZ = ABS(ZT-dKnetMin)
!  daux=MAX(daux1,daux2)
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist2 =MIN(daux,dist2)
! 
! END DO
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE SKNET_elem
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE SME_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 :: dZ
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
!  dBeta=2d0*(Z-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
! ELSE
!  IF (Z.lt.mySigma%mySegment(iSeg)%Min) dBeta=2d0*(mySigma%mySegment(iSeg)%min-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
!  IF (Z.gt.mySigma%mySegment(iSeg)%Max) dBeta=2d0*(mySigma%mySegment(iSeg)%max-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
! END IF
! 
! ! First screw
! dist1 = 5d0
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
! daux = max(daux,-dCut)
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist1 =MIN(daux,dist1)
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! First screw
! dist2 = 5d0
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
! 
! ! First the point needs to be transformed back to time = 0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! CALL SecondaryProfile(XT,YT,ZT,dCut,iSeg)
! daux = max(daux,-dCut)
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist2 =MIN(daux,dist2)
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE SME_elem
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE FOERD_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 :: dZ
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
!  dBeta=2d0*(Z-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
! ELSE
!  IF (Z.lt.mySigma%mySegment(iSeg)%Min) dBeta=2d0*(mySigma%mySegment(iSeg)%min-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
!  IF (Z.gt.mySigma%mySegment(iSeg)%Max) dBeta=2d0*(mySigma%mySegment(iSeg)%max-mySigma%mySegment(iSeg)%min)*myPI/mySigma%mySegment(iSeg)%t 
! END IF
! 
! ! First screw
! dist1 = 1d2
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! 
! dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
! !  daux=daux
!  IF (daux.le.0d0) THEN
!   daux=max(-dZ,daux)
!  ELSE
!   daux=daux
!  END IF
! ELSE
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist1 =MIN(daux,dist1)
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! First screw
! dist2 = 1d2
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
! 
! ! First the point needs to be transformed back to time = 0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! 
! dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
! !  daux=daux
!  IF (daux.le.0d0) THEN
!   daux=max(-dZ,daux)
!  ELSE
!   daux=daux
!  END IF
! ELSE
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist2 =MIN(daux,dist2)
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE FOERD_elem
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE CONE_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax,dH
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 :: dZ
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! dH = mySigma%mySegment(iSeg)%Max-mySigma%mySegment(iSeg)%Min
! dBeta = 0d0
! 
! ! First screw
! dist1 = 1d2
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = 0d0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = (mySigma%mySegment(iSeg)%Max-ZT)/dH
! 
! CALL Schnecke_Cone(XTT,YTT,ZTT,dH,daux)
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist1 =MIN(daux,dist1)
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! First screw
! dist2 = 1d2
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! dAlpha = 0d0
! 
! ! First the point needs to be transformed back to time = 0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = (mySigma%mySegment(iSeg)%Max-ZT)/dH
! 
! CALL Schnecke_Cone(XTT,YTT,ZTT,dH,daux)
! 
! IF (ZT.GE.mySigma%mySegment(iSeg)%Min.and.ZT.LE.mySigma%mySegment(iSeg)%Max) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-mySigma%mySegment(iSeg)%Min),ABS(ZT-mySigma%mySegment(iSeg)%Max))
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist2 =MIN(daux,dist2)
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE CONE_elem
! !
! !********************GEOMETRY****************************
! !
! SUBROUTINE KNET_elem(X,Y,Z,t,iSeg,dist1,dist2,inpr)
! IMPLICIT NONE
! REAL*8 X,Y,Z,t
! REAL*8 :: dKnetMin,dKnetMax
! INTEGER :: inpr,iBndr, l, k, iSeg,lKnet
! REAL*8 :: dAlpha,XT,YT,ZT,XB,YB,ZB
! REAL*8 :: dBeta,XTT,YTT,ZTT,dist
! REAL*8 :: dZ
! REAL*8 dScale, daux,dist1, dist2,dEps,dCut
! REAL*8 myPI
! 
! dEps = mySigma%a/1d5
! 
! myPI = dATAN(1d0)*4d0
! 
! IF (Z.ge.mySigma%mySegment(iSeg)%Min.and.Z.le.mySigma%mySegment(iSeg)%Max) then
!  DO l=1, mySigma%mySegment(iSeg)%N
!   dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
!   dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
!   IF (Z.GE.dKnetMin.AND.Z.LT.dKnetMax) THEN
!    lKnet = l
!    EXIT
!   END IF
!  END DO
! ELSE
!  IF (Z.lt.mySigma%mySegment(iSeg)%Min) lknet = 1
!  IF (Z.gt.mySigma%mySegment(iSeg)%Max) lknet = mySigma%mySegment(iSeg)%N
! END IF
! 
! ! First screw
! dist1 = 5d0
! XB = X
! YB = Y-mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! dAlpha = mySigma%mySegment(iSeg)%StartAlpha - t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)
! 
! dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
! dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
! dBeta = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! 
! IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
! !  WRITE(*,*) zt
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist1 =MIN(daux,dist1)
! 
! END DO
! 
! IF (dist1.LT.0d0) THEN
!  inpr = 101
! END IF
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ! Second screw
! dist2 = 5d0
! XB = X
! YB = Y+mySigma%a/2d0
! ZB = Z
! 
! ! First the point needs to be transformed back to time = 0
! IF (mySigma%GANGZAHL .EQ. 1) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 2) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/2d0)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 3) dAlpha = mySigma%mySegment(iSeg)%StartAlpha -t*myPI*(myProcess%Umdr/3d1)*myProcess%ind
! IF (mySigma%GANGZAHL .EQ. 4) dAlpha = mySigma%mySegment(iSeg)%StartAlpha + (-t*myPI*(myProcess%Umdr/3d1)+myPI/4d0)*myProcess%ind
! 
! ! First the point needs to be transformed back to time = 0
! XT = XB*cos(dAlpha) - YB*sin(dAlpha)
! YT = XB*sin(dAlpha) + YB*cos(dAlpha)
! ZT = ZB
! 
! DO l=MAX(1,lKnet-1),MIN(mySigma%mySegment(iSeg)%N,lKnet+1)
! 
! dKnetMin = mySigma%mySegment(iSeg)%Min + DBLE(l-1)*mySigma%mySegment(iSeg)%D-dEps
! dKnetMax = mySigma%mySegment(iSeg)%Min + DBLE(l  )*mySigma%mySegment(iSeg)%D+dEps
! dBeta = myPI*DBLE(l-1)*mySigma%mySegment(iSeg)%Alpha/1.8d2
! 
! ! Next the point needs to be transformed back to Z = 0 level
! XTT = XT*cos(dBeta) - YT*sin(dBeta)
! YTT = XT*sin(dBeta) + YT*cos(dBeta)
! ZTT = ZT
! 
! CALL Schnecke_nGang(XTT,YTT,daux)
! 
! IF (ZT.GE.dKnetMin.and.ZT.LE.dKnetMax) THEN
!  daux=daux
! ELSE
!  dZ = MIN(ABS(ZT-dKnetMin),ABS(ZT-dKnetMax))
! !  WRITE(*,*) zt
!  IF (daux.le.0d0) THEN
!   daux=dZ
!  ELSE
!   daux=SQRT(daux**2d0+dZ**2d0)
!  END IF
! END IF
! 
! dist2 =MIN(daux,dist2)
! 
! END DO
! 
! IF (dist2.LT.0d0) THEN
!  inpr = 102
! END IF
! 
! END SUBROUTINE KNET_elem
! !
! !-----------------------------------------------------------
! !
! SUBROUTINE ZME_disc(XP,YP,ZP,ZPMIN,ZPMAX,iSeg,Dif)
! REAL*8 XP,YP,ZP, Dif,ZPMIN,ZPMAX
! REAL*8  Dz, a, delta, R_out
! REAL*8 d1,d2,d3 !,dShift
! INTEGER iSeg
! ! INTEGER iSide
! 
! Dz=mySigma%Dz
! a=mySigma%a
! delta=mySigma%mySegment(iSeg)%ZME_gap_SG
! 
! R_out = Dz/2D0 - delta
! 
! ! IF (iSide.EQ.1) dShift = -10.0d0
! ! IF (iSide.EQ.2) dShift = +10.0d0
! 
! d1 = dSQRT(XP*XP + YP*YP)-R_out
! d2 = ZP - (ZPMAX)
! d3 = (ZPMIN) - ZP
! 
! Dif = MAX(d1,d2,d3)
! 
! END SUBROUTINE ZME_disc
! !
! !-----------------------------------------------------------
! !
! SUBROUTINE Schnecke_nGang(XP,YP,Dif)
! !Dz= Zylinderdurchmesser, a=Achsabstand
! !delta= Spiel Schneckenkamm-Gehaeusewand
! !s= Spiel Schnecke-Schnecke (from Sigma!)
! REAL*8  Dz, a, delta, s, Ds    
! REAL*8 dTheta, XP, YP, Dif, XB, YB
! REAL*8 R1, R2, R3, dAlpha,dGamma,dCutAngle,dRadius,dN
! REAL*8  :: myPI = dATAN(1d0)*4d0
! INTEGER iAux
! 
! dN=DBLE(mySigma%GANGZAHL)
! Dz=mySigma%Dz
! a=mySigma%a
! delta=mySigma%delta
! s=mySigma%s
! 
! Ds= Dz - 2d0*delta
! 
! R2 = Dz/2D0 - delta
! R1 = a-s-R2
! R3 = a-s
! 
! dRadius = SQRT(XP*XP + YP*YP)
! IF (XP.LT.0d0) THEN 
!  dGamma  = 0d0
! ELSE
!  dGamma = myPi
! END IF
! 
! dGamma = dGamma + dATAN(YP/XP)
! dCutAngle = 2d0*myPi/(dN)
! iAux = FLOOR(dGamma/dCutAngle)
! dGamma = dGamma - DBLE(iAux)*dCutAngle - 0.5d0*dCutAngle
! dGamma = ABS(dGamma)
! XB = (dRadius*COS(dGamma))
! YB = (dRadius*SIN(dGamma))
! 
! dTheta = 0.5d0*(myPI/dN - 2d0*ACOS(R3/(2d0*R2)))
! IF (dTheta.lt.0d0) THEN
!  WRITE(*,*) "the screw could not be created (wrong theta)", dTheta
!  stop
! ENDIF
! 
! Dif = dSQRT(XB*XB+YB*YB)-R2
! 
! IF ( dGamma.LT.dTheta .AND. dGamma.GT.0d0 ) THEN
!   Dif = dSQRT(XB*XB+YB*YB)-R1
! END IF
! 
! IF ( dGamma.GT.dTheta .AND. dGamma.LT.(2d0*myPI/(dN*2d0)-dTheta)) THEN
!   Dif = dSQRT((XB + R2*dCOS(dTheta))**2d0 +(YB + R2*dSIN(dTheta))**2d0)-R3
! END IF
! 
! ! IF ( dAlpha.LT.-dTheta .AND. dAlpha.GT.(-myPI/2d0+dTheta)) THEN
! !   Dif = dSQRT((XB + dSIGN(1d0,XB)*R2*dCOS(dTheta))**2d0 +(YB - dSIGN(1d0,XB)*R2*dSIN(dTheta))**2d0)-R3
! ! END IF
! 
! 
! END SUBROUTINE Schnecke_nGang
! !-----------------------------------------------------------
! SUBROUTINE Schnecke_CONE(XP,YP,ZP,H,Dif)
! !Dz= Zylinderdurchmesser, a=Achsabstand
! !delta= Spiel Schneckenkamm-Gehaeusewand
! !s= Spiel Schnecke-Schnecke (from Sigma!)
! REAL*8  Dz, a, delta, s, Ds,H    
! REAL*8 dTheta, XP, YP, ZP, Dif, XB, YB
! REAL*8 R1, R2, R3, dCosAlpha,dGamma,dCutAngle,dRadius,dN
! REAL*8  :: myPI = dATAN(1d0)*4d0
! INTEGER iAux
! 
! Dz=mySigma%Dz
! a=mySigma%a
! delta=mySigma%delta
! s=mySigma%s
! 
! Ds= Dz - 2d0*delta
! 
! R2 = Dz/2D0 - delta
! R1 = a-s-R2
! R3 = a-s
! 
! dCosAlpha = R1/H
! 
! Dif = dSQRT(XP*XP + YP*YP) - ZP*R1
! dif = dCosAlpha*Dif
! 
! END SUBROUTINE Schnecke_CONE
! !------------------------------------------------------------
! SUBROUTINE SecondaryProfile(XB,YB,ZB,d,k)
! REAL*8 XB,YB,ZB,d, XT,YT,XP,YP
! INTEGER k
! REAL*8 dGamma,dRadius,dCutAngle,dCutRadius,dBeta,dStartAngle
! ! REAL*8  :: dCutDepth=15d0,dCutWidth=28d0,dCutPitch !,dCutAngleF=0.33333333333333d0 !2d0*0.66666666667d0
! ! INTEGER :: nCut = 4,SSE_Type = 3
! REAL*8  :: dCutDepth,dCutWidth,dCutPitch !,dCutAngleF=0.33333333333333d0 !2d0*0.66666666667d0
! INTEGER :: nCut,SSE_Type
! INTEGER iAux
! REAL*8  :: myPI = dATAN(1d0)*4d0
! 
! nCut      = mySigma%mySegment(k)%SecProf_N
! SSE_Type  = mySigma%mySegment(k)%SecProf_I
! dCutDepth = mySigma%mySegment(k)%SecProf_D
! dCutWidth = mySigma%mySegment(k)%SecProf_W
! dCutPitch = mySigma%mySegment(k)%SecProf_L
! 
! dCutRadius = mySigma%Dz/2D0 !- mySigma%delta
! 
! ! if (k.eq.1) dStartAngle =  myPi* 0d0/180d0
! ! if (k.eq.2) dStartAngle = 2d0*dCutAngleF*myPI/mySigma%mySegment(1)%t !+myPi*45d0/180d0
! 
! ! XT = XB
! ! YT = YB
! ! dBeta= -2d0*dCutAngleF*(ZB-mySigma%mySegment(k)%Min)*myPI/40d0 !mySigma%mySegment(k)%t
! 
! dBeta=2d0*(ZB-mySigma%mySegment(k)%min)*myPI/dCutPitch
! XT = XB*cos(dBeta) - YB*sin(dBeta)
! YT = XB*sin(dBeta) + YB*cos(dBeta)
! 
! dRadius = SQRT(XT*XT + YT*YT)
! IF (XT.LT.0d0) THEN 
!  dGamma  = 0d0
! ELSE
!  dGamma = myPi
! END IF
! 
! dGamma = dGamma + dATAN(YT/XT)
! dCutAngle = 2d0*myPi/DBLE(nCut)
! iAux = FLOOR(dGamma/dCutAngle)
! dGamma = dGamma - DBLE(iAux)*dCutAngle - 0.5d0*dCutAngle
! XB = (dRadius*COS(dGamma))
! YB = (dRadius*SIN(dGamma))
! 
! IF (SSE_Type.eq.1) THEN
!  d = SQRT((XB-dCutRadius)*(XB-dCutRadius) + YB*YB) - dCutDepth
! END IF
! 
! IF (SSE_Type.eq.2) THEN
!  d = MAX(-dCutDepth-(XB-dCutRadius),-0.5d0*dCutWidth-YB,-0.5d0*dCutWidth+YB) !XB-dCutRadius+dCutDepth !MAX(XB-dCutRadius+dCutDepth, YB-0.5d0*dCutWidth,0.5d0*dCutWidth-YB)
! ENDIF
! 
! IF (SSE_Type.eq.3) THEN
!  d = FindDist(YB,XB-dCutRadius)
! ENDIF
! 
! ! IF (SSE_Type.eq.2) THEN
! !  d = SQRT(dCutWidth*dCutWidth*(XB-dCutRadius)*(XB-dCutRadius) + dCutDepth*dCutDepth*YB*YB)-dCutWidth*dCutDepth
! !  d=0.1d0*d
! ! END IF
! 
!  CONTAINS
! 
!  FUNCTION FindDist(ppx, ppy)
!  REAL*8  ppx, ppy
!  REAL*8  FindDist
!  REAL*4  px, py, pz, ax, ay, az, bx, by, bz, tol
!  INTEGER np,noptfd
!  REAL*4  dpmin, fdmin, cx, cy, cz
!  INTEGER nlim, itrun, nerr
!  INTEGER i
!  REAL*8 ::  A,B,C,dSGN,dF,dDist
!  REAL*8 :: memory(2,33)
! 
!  C=-dCutDepth
!  A= -C/(dCutWidth*dCutWidth)*4d0
!  B= 0d0
! 
!  tol = 1e-6
!  np =1
!  noptfd = 0
! 
!  px = REAL(ppx)
!  py = REAL(ppy)
!  pz = 0e0
! 
!  dF = A*px*px + B*px + C- py
! ! dpmin = dF
!  dSGN = (dF)/ABS(dF)
! 
!  ax =  0e0
!  ay =  REAL(C)
!  az =  0e0
! 
!  bx = REAL(ppx)
!  by = REAL(A*px*px + B*px + C)
!  bz = 0e0
! 
!  DO i=1,33
!   cx = 0.5d0*(ax + bx)
!   cy = REAL(A*cx*cx + B*cx + C)
! 
!   memory(:,i) = [cx,cy]
! 
!   IF (ABS(ay).lt.ABS(by)) THEN
!    bx = cx
!    by = cy
!   ELSE
!    ax = cx
!    ay = cy
!   END IF
! 
!   IF (((ax-bx)**2d0+(ay-by)**2d0).lt. 1d-4*mySigma%Dz) EXIT
! 
!  END DO
! 
!  dDist = SQRT((px-cx)**2d0+(py-cy)**2d0)
! 
!  FindDist = dSGN*ABS(dDist) !max(-dF,ppx-mySigma%Dz/2D0)
! 
!  END FUNCTION FindDist
! 
! END SUBROUTINE SecondaryProfile
! !
!!!  End of geometry
!-----------------------------------------------------------
SUBROUTINE GetVeloMixerVal(X,Y,Z,ValU,ValV,ValW,iP,t)
IMPLICIT NONE
REAL*8 :: myPI
INTEGER iP
REAL*8 X,Y,Z,ValU,ValV,ValW,t

myPI = dATAN(1d0)*4d0

SELECT CASE(iP)
 CASE (101) ! Y positiv
  ValU =  -myPI*(Y-mySigma%a/2d0)*(myProcess%Umdr/3d1)*myProcess%ind
 CASE (102) ! Y negativ
  ValU =  -myPI*(Y+mySigma%a/2d0)*(myProcess%Umdr/3d1)*myProcess%ind
 END SELECT


ValV = myPI*X*(myProcess%Umdr/3d1)*myProcess%ind
ValW = 0d0


END SUBROUTINE GetVeloMixerVal
!-----------------------------------------------------------
SUBROUTINE Sigma_AdjustTimeParameters(dt,tmax,dtOut)

REAL*8 dt,tmax,dtOut
REAL*8 :: tFact = 10d0

dt    = 1d0/(myProcess%Umdr/6d1)/(tFact*1d3)
dtOut = tFact*24.9999d0*dt
tmax  = myProcess%Dreh/(myProcess%Umdr/6d1)

! IF (bGeoTest) dt = 125*dt

END SUBROUTINE Sigma_AdjustTimeParameters

END MODULE Sigma_User
