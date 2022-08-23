SUBROUTINE SolToFile(iOutput)
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
use solution_io, only: write_pres_sol,write_vel_sol,write_time_sol,write_q2_sol
use var_QuadScalar, only: myDump,istep_ns,fieldPtr
use, intrinsic :: iso_c_binding

IMPLICIT NONE
INTEGER iOutput
integer :: out_lev
integer :: ndof
character(60) :: fieldName

type(fieldPtr), dimension(3) :: packed

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER ifilen,iOut,nn
DATA ifilen/0/

! filen is a global counter of the output file number
IF (iOutput.LT.0) THEN
 ifilen=ifilen+1
 iOut=MOD(ifilen+insavn-1,insavn)+1
ELSE
 iOut = iOutput
END IF

!call WriteSol_Velo(iOut,0,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
nn = knel(nlmax)

ndof = KNVT(NLMAX) + KNAT(NLMAX) + KNET(NLMAX) + KNEL(NLMAX)

call write_vel_sol(iOut,0,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
                   QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

!call WriteSol_Pres(iOut,0,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
!     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)

call write_pres_sol(iOut,0,nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Elements,LinSc%ValP(NLMAX)%x)

!! CALL WriteSol_Coor(iOut,0,DWORK(L(KLCVG(NLMAX))),QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)

call write_time_sol(iOut,istep_ns, timens)


! This is how to output a custom field with 3 components:
fieldName = "coordinates"

packed(1)%p => QuadSc%auxU
packed(2)%p => QuadSc%auxV
packed(3)%p => QuadSc%auxW

call write_q2_sol(fieldName, iOut,0,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
                  3, packed)


END SUBROUTINE SolToFile
!
! ----------------------------------------------
!
SUBROUTINE SolFromFileRepart(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
use solution_io
use var_QuadScalar, only: myDump,istep_ns,fieldPtr

IMPLICIT NONE
INTEGER mfile,iLevel,nn
character(60) :: cInFile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

character(60) :: FileA
character(60) :: FileB
character(60) :: fieldName

integer :: ndof

type(fieldPtr), dimension(3) :: packed

nn = knel(nlmax)

ndof = KNVT(NLMAX) + KNAT(NLMAX) + KNET(NLMAX) + KNEL(NLMAX)

!subroutine read_vel_sol_single(startFrom, iiLev,nn, nmin, nmax,elemmap,edofs, u, v, w)
! read in the velocity solution
call read_vel_sol_single(cInFile,iLevel-1,nn,NLMIN,NLMAX,&
                         coarse%myELEMLINK,myDump%Vertices,&
                         QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

!FileA='single_v'
!call write_vel_test(FileA, nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)
!
! read in the pressure solution
call read_pres_sol_single(cInFile,iLevel-1,nn,NLMIN,NLMAX,&
                          coarse%myELEMLINK,myDump%Elements,&
                          LinSc%ValP(NLMAX)%x)

!FileB='single_p'
!call write_pres_test(FileB, nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Elements,LinSc%ValP(NLMAX)%x)
!
!FileA='time.dmp'
call read_time_sol_single(cInFile, istep_ns, timens)

fieldName = "coordinates"

packed(1)%p => QuadSc%auxU
packed(2)%p => QuadSc%auxV
packed(3)%p => QuadSc%auxW

!subroutine read_q2_sol_single(fieldName, startFrom, iiLev,nn, nmin, nmax,elemmap,edofs, icomp, field_pack)
call read_q2_sol_single(fieldName,cInFile,iLevel-1,nn,NLMIN,NLMAX,&
                        coarse%myELEMLINK,myDump%Vertices,&
                        3,packed)
! 
! call read_q2_sol(fieldName, cInFile,ilevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
!                  3, packed)

END SUBROUTINE SolFromFileRepart
!
! ----------------------------------------------
!
SUBROUTINE DumpAllValues(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
use solution_io
use var_QuadScalar, only: myDump,istep_ns,fieldPtr

IMPLICIT NONE
INTEGER mfile,iLevel,nn
character(60) :: cInFile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

character(60) :: FileA
character(60) :: FileB

integer :: ndof

nn = knel(nlmax)

ndof = KNVT(NLMAX) + KNAT(NLMAX) + KNET(NLMAX) + KNEL(NLMAX)


FileA='dump_v'
call write_vel_test(FileA, nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

FileB='dump_p'
call write_pres_test(FileB, nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Elements,LinSc%ValP(NLMAX)%x)

END SUBROUTINE DumpAllValues
!
! ----------------------------------------------
!
SUBROUTINE SolToFile_Compact(iOutput)
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
use var_QuadScalar, only: myDump,mg_mesh

IMPLICIT NONE
INTEGER iOutput

INTEGER ifilen,iOut,nn
DATA ifilen/0/

! filen is a global counter of the output file number
IF (iOutput.LT.0) THEN
 ifilen=ifilen+1
 iOut=MOD(ifilen+insavn-1,insavn)+1
ELSE
 iOut = iOutput
END IF

call WriteSol_Velo(iOut,0,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

nn = knel(nlmax)
call WriteSol_Pres(iOut,0,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
    LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)

CALL WriteSol_Coor(iOut,0,mg_mesh%level(nlmax+1)%dcorvg,QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)

call WriteSol_Time(iOut)


END SUBROUTINE SolToFile_Compact
!
! ----------------------------------------------
!
SUBROUTINE SolFromFile_Compact(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic,Temperature,MaterialDistribution
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
use var_QuadScalar, only: myDump,mg_mesh

IMPLICIT NONE
INTEGER mfile,iLevel,nn
character(60) :: cInFile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

nn = knel(nlmax)

CALL ReadSol_Velo(cInFile,iLevel,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

CALL ReadSol_Pres(cInFile,iLevel,LinSc%ValP(NLMAX)%x,LinSc%AuxP(NLMAX)%x(1),&
     LinSc%AuxP(NLMAX)%x(nn+1),LinSc%AuxP(NLMAX)%x(2*nn+1),LinSc%AuxP(NLMAX)%x(3*nn+1),nn)
    
CALL ReadSol_Coor(cInFile,iLevel,mg_mesh%level(nlmax+1)%dcorvg,QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof)

! This part here is responsible for creation of structures enabling the mesh coordinate 
! transfer to the master node so that it can create the corresponding matrices
IF (myid.EQ.0) THEN
  CALL CreateDumpStructures(0)
ELSE
  LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
  CALL CreateDumpStructures(LevDif)
END IF

ilev = LinSc%prm%MGprmIn%MedLev

nLengthV = (2**(ilev-1)+1)**3
nLengthE = mg_mesh%level(NLMIN)%nel

ALLOCATE(SendVect(3,nLengthV,nLengthE))

CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                            mg_mesh%level(ilev)%kvert,&
                            nLengthV,&
                            nLengthE,&
                            mg_mesh%level(ilev)%nel,&
                            mg_mesh%level(ilev)%nvt)
DEALLOCATE(SendVect)


END SUBROUTINE SolFromFile_Compact
!
! ----------------------------------------------
!
SUBROUTINE SolFromFile(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic,Temperature,MaterialDistribution
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
use solution_io
use var_QuadScalar, only: myDump,istep_ns,fieldPtr
use var_QuadScalar, only: GenLinScalar

IMPLICIT NONE
INTEGER mfile,iLevel,nn
character(60) :: cInFile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

character(60) :: FileA
character(60) :: FileB
character(60) :: fieldName

integer :: ndof,iFld

type(fieldPtr), dimension(3) :: packed

nn = knel(nlmax)

ndof = KNVT(NLMAX) + KNAT(NLMAX) + KNET(NLMAX) + KNEL(NLMAX)

call read_vel_sol(cInFile,iLevel-1,nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

! read in the pressure solution
call read_pres_sol(cInFile,iLevel-1,nn,NLMIN,NLMAX,coarse%myELEMLINK,&
                   myDump%Elements,LinSc%ValP(NLMAX)%x)

!FileB='new_p'
!call write_pres_test(FileB, nn,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Elements,LinSc%ValP(NLMAX)%x)

IF (allocated(Temperature)) then
 fieldName = "temperature"
 packed(1)%p => QuadSc%auxU
 call read_q2_sol(fieldName,cInFile,iLevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,1, packed)
 Temperature = QuadSc%auxU
END IF                  

IF (allocated(MaterialDistribution)) then
 fieldName = "MaterialDistribution"
 QuadSc%auxU = 0
 packed(1)%p => QuadSc%auxU
 call read_q2_sol(fieldName,cInFile,iLevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,1, packed)
 MaterialDistribution(NLMAX+iLevel-1)%x(1:knel(NLMAX+iLevel-1)) = QuadSc%auxU((knvt(NLMAX+iLevel-1) + knat(NLMAX+iLevel-1) + knet(NLMAX+iLevel-1))+1:) 
END IF                  

if (allocated(GenLinScalar%Fld)) then
 DO iFld=1,GenLinScalar%nOfFields
  fieldName = adjustl(trim(GenLinScalar%prm%cField(iFld)))
  QuadSc%auxU = 0
  packed(1)%p => QuadSc%auxU
  call read_q2_sol(fieldName,cInFile,iLevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,1, packed)
  GenLinScalar%fld(iFld)%Val = QuadSc%auxU
 END DO
end if

call read_time_sol(cInFile, istep_ns, timens)

fieldName = "coordinates"

packed(1)%p => QuadSc%auxU
packed(2)%p => QuadSc%auxV
packed(3)%p => QuadSc%auxW

call read_q2_sol(fieldName, cInFile,ilevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
                 3, packed)


if(bViscoElastic)then
  CALL ReadSol_Visco(cInFile, iLevel)
end if

END SUBROUTINE SolFromFile
!
! ----------------------------------------------
!
SUBROUTINE SolFromFile_heat(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,bViscoElastic,Temperature,MaterialDistribution
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:Tracer
use solution_io
use var_QuadScalar, only: myDump,istep_ns,fieldPtr

IMPLICIT NONE
INTEGER mfile,iLevel,nn
character(60) :: cInFile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------
INTEGER nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

character(60) :: FileA
character(60) :: FileB
character(60) :: fieldName

integer :: ndof

type(fieldPtr), dimension(3) :: packed

nn = knel(nlmax)

ndof = KNVT(NLMAX) + KNAT(NLMAX) + KNET(NLMAX) + KNEL(NLMAX)

IF (allocated(Temperature)) then
 fieldName = "temperature"
 packed(1)%p => QuadSc%auxU
 call read_q2_sol(fieldName,cInFile,iLevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,1, packed)
 Temperature = QuadSc%auxU
END IF                  

call read_time_sol(cInFile, istep_ns, timens)

fieldName = "coordinates"

packed(1)%p => QuadSc%auxU
packed(2)%p => QuadSc%auxV
packed(3)%p => QuadSc%auxW

call read_q2_sol(fieldName, cInFile,ilevel-1,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
                 3, packed)
                 
call read_heatingStatus(cInFile)
                 

END SUBROUTINE SolFromFile_heat
!
!-------------------------------------------------------------------------------
!
SUBROUTINE read_heatingStatus(startFrom)
USE Sigma_User,ONLY:mySigma
USE PP3D_MPI, ONLY:myid
character(60), intent(in) :: startFrom
integer iSeg

open(file='_dump/heatstatus_'//adjustl(trim(startFrom))//'.dmp',unit=472)

do iSeg=1,mySigma%NumberOfSeg
 IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
  read(472,*) mySigma%mySegment(iSeg)%TemperatureSensor%HeatingStatus 
  if (mySigma%mySegment(iSeg)%TemperatureSensor%HeatingStatus) then
   mySigma%mySegment(iSeg)%UseHeatSource = mySigma%mySegment(iSeg)%HeatSourceMax
  else
   mySigma%mySegment(iSeg)%UseHeatSource = mySigma%mySegment(iSeg)%HeatSourceMin
  end if
  if (myid.eq.1) WRITE(*,'(A,I0,L,ES12.4)') 'iSeg, HeatingStatus, HeatSource[kW] :',&
   iSeg,mySigma%mySegment(iSeg)%TemperatureSensor%HeatingStatus,mySigma%mySegment(iSeg)%UseHeatSource
 END IF
end do

close(472)

END SUBROUTINE read_heatingStatus
!
!-------------------------------------------------------------------------------
!
SUBROUTINE write_heatingStatus(iO)
USE Sigma_User,ONLY:mySigma
USE PP3D_MPI, ONLY:myid

integer :: iO
integer iSeg
character*(200) cFile

write(cFile,'(A,I0,A)') '_dump/heatstatus_',iO,'.dmp'
open(file=ADJUSTL(TRIM(cFile)),unit=472)

do iSeg=1,mySigma%NumberOfSeg
 IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
  WRITE(472,'(L)') mySigma%mySegment(iSeg)%TemperatureSensor%HeatingStatus 
 END IF
end do

close(472)

END SUBROUTINE write_heatingStatus
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Time(iOut)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
INTEGER iInd
INTEGER pID
character cFile*(40)

 IF (myid.eq.showid) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_Time.dmp'
  WRITE(MTERM,*) 'Releasing Time level into: "', TRIM(ADJUSTL(cFile)),'"'
 END IF

 IF (myid.eq.0) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_Time.dmp'
  OPEN(321,FILE=TRIM(ADJUSTL(cFile)))
  WRITE(321,*) TIMENS
  CLOSE(321)
 END IF

END SUBROUTINE WriteSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Velo(iInd,iiLev,Field1,Field2,Field3)
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iInd,iiLev
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'

CALL WriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE WriteSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Visco(iInd,iiLev)
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
INTEGER        iInd,iiLev
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'

CALL WriteSol(iInd,iType,iiLev,cFF)

END SUBROUTINE WriteSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Pres(iInd,iiLev,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*),Field4(*),Field5(*)
INTEGER        i,iInd,nn,iiLev
INTEGER        :: iType= 2
CHARACTER*(20) :: cFF='Pressure'

IF (myid.ne.0) THEN
 DO i=1,nn
  Field2(i) = Field1(4*(i-1)+1) ! mean
  Field3(i) = Field1(4*(i-1)+2) ! d/dx
  Field4(i) = Field1(4*(i-1)+3) ! d/dy 
  Field5(i) = Field1(4*(i-1)+4) ! d/dz
 END DO
!  WRITE(*,*) Field2(1:nn)
END IF


CALL WriteSol(iInd,iType,iiLev,cFF,Field2,Field3,Field4,Field5)

END SUBROUTINE WriteSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WriteSol_Coor(iInd,iiLev,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iOut,iType,iiLev
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20)
  END SUBROUTINE WriteSol
end interface
REAL*8         Field(3,*),Field1(*),Field2(*),Field3(*)
INTEGER        i,iInd,nn,iiLev
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Coordinates'

IF (myid.ne.0) THEN
 DO i=1,nn
  Field1(i) = Field(1,i)
  Field2(i) = Field(2,i)
  Field3(i) = Field(3,i)
 END DO
END IF

CALL WriteSol(iInd,iType,iiLev,cFF,Field1,Field2,Field3)

END SUBROUTINE WriteSol_Coor

SUBROUTINE WriteSol(iOut,iType,iiLev,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE var_QuadScalar,ONLY:myDump, ViscoSc
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
INTEGER iOut,iType,iiLev
REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
CHARACTER cField*(20)
!----------------------------------------------------
INTEGER i,ivt,jvt,jel,kel,iP,nLengthE,nLengthV
REAL*8,ALLOCATABLE :: Field(:,:),auxField(:,:)
REAL*8 dMaxNel
INTEGER pnel,pID
CHARACTER cFile*(40)

IF (myid.NE.0) THEN

 ILEV = NLMIN

 nLengthE = 8**((NLMAX+iiLev)-1)

 ! (2^(NLMAX+iiLev)+1)^3
 ! dofs on a cube
 nLengthV = (2**((NLMAX+iiLev))+1)**3

!  WRITE(*,*)  'WRITE :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(NLMIN))
 CALL COMM_Maximum(dMaxNel)

 IF (myid.eq.showid) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
  WRITE(MTERM,*) 'Releasing current '//TRIM(ADJUSTL(cField))//' solution into: "'//ADJUSTL(TRIM(cFile)),'"'
 END IF

 IF (myid.eq.0) THEN
  WRITE(cFile(1:),'(A,I2.2,A)') '_dump/',iOut,'_'//TRIM(ADJUSTL(cField))//'.dmp'
  OPEN(321,FILE=TRIM(ADJUSTL(cFile)))
 END IF

 IF (iType.eq.1) THEN
  IF (Present(Field1)) CALL CollectVertField(Field1)
  IF (Present(Field2)) CALL CollectVertField(Field2)
  IF (Present(Field3)) CALL CollectVertField(Field3)
 END IF

 IF (iType.eq.2) THEN
  IF (Present(Field1)) CALL CollectElemField(Field1)
  IF (Present(Field2)) CALL CollectElemField(Field2)
  IF (Present(Field3)) CALL CollectElemField(Field3)
  IF (Present(Field4)) CALL CollectElemField(Field4)
 END IF

 IF (iType.eq.3) THEN
   CALL CollectVertField(ViscoSc%val11)
   CALL CollectVertField(ViscoSc%val22)
   CALL CollectVertField(ViscoSc%val33)
   CALL CollectVertField(ViscoSc%val12)
   CALL CollectVertField(ViscoSc%val13)
   CALL CollectVertField(ViscoSc%val23)
 END IF

 IF (myid.eq.0) THEN
  CLOSE(321)
 END IF

 CONTAINS
!
! -----------------------------------------------------------------
!
SUBROUTINE CollectVertField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthV,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 

  ! row-column based output
  ! myDump%Vertices stores the global indices (on the output level)
  ! of vertices inside a coarse grid element
  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthV*KNEL(NLMIN),0)
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthV,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 
    ALLOCATE(auxField(nLengthV,pnel)) 
   END IF
   CALL RECVI_myMPI(pnel,pID)

   CALL RECVD_myMPI(auxField,pnel*nLengthV,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthV
     Field(ivt,iel) = auxField(ivt,I)
    END DO
   END DO
  END DO
 
  DEALLOCATE(auxField) 
  DO iel=1,KNEL(NLMIN)
   !WRITE(321,'(<nLengthV>ES14.6)') Field(1:nLengthV,iel)
   WRITE(321,*) Field(1:nLengthV,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectVertField
! -----------------------------------------------------------------
SUBROUTINE CollectElemField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 

  ! build a structure for each coarse mesh
  ! element that stores the values
  ! of the fine mesh dofs in it
  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthE
    ! jvt is the ivt-th local P1 dof
    ! of local element IEL
    jvt = myDump%Elements(IEL,ivt)
    Field(ivt,iel) = xField(jvt)
   END DO
  END DO
 
  CALL SENDD_myMPI(Field,nLengthE*KNEL(NLMIN),0)
 
 ELSE

  DO pID =1,subnodes

   ! number of sub-elements of a coarse mesh
   ! element: nLengthE
   CALL RECVI_myMPI(nLengthE,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)

    ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 

    ! an array that stores the values
    ! of the field from the partitions
    ALLOCATE(auxField(nLengthE,pnel)) 
   END IF

   ! set the pnel to the number of coarse
   ! grid elements for process pID
   CALL RECVI_myMPI(pnel,pID)

   CALL RECVD_myMPI(auxField,pnel*nLengthE,pID)

   DO I=1,pnel

   ! map from local coarse grid element number I from
   ! partition pID to global coarse grid elmenet number IEL
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthE
     Field(ivt,iel) = auxField(ivt,I)
    END DO
   END DO
  END DO
 
  DEALLOCATE(auxField) 
!  WRITE(321,'(A)') '-- - ---'
  DO iel=1,KNEL(NLMIN)
   WRITE(321,*) Field(1:nLengthE,iel)
   !WRITE(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
  END DO
 END IF

DEALLOCATE(Field) 

END SUBROUTINE CollectElemField
! -----------------------------------------------------------------

END SUBROUTINE WriteSol
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Coor(cInFile,iLevel,Field,Field1,Field2,Field3,nn)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field(3,*),Field1(*),Field2(*),Field3(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Coordinates'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field(1,i) = Field1(i) 
  Field(2,i) = Field2(i) 
  Field(3,i) = Field3(i) 
 END DO
END IF

END SUBROUTINE ReadSol_Coor
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Velo(cInFile,iLevel,Field1,Field2,Field3)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*)
INTEGER        iLevel
INTEGER        :: iType= 1
CHARACTER*(20) :: cFF='Velocity'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field1,Field2,Field3)

END SUBROUTINE ReadSol_Velo
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Visco(cInFile,iLevel)
USE PP3D_MPI, ONLY:myid
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
INTEGER        iLevel
INTEGER        :: iType= 3
CHARACTER*(20) :: cFF='Stress'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF)

END SUBROUTINE ReadSol_Visco
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Pres(cInFile,iLevel,Field1,Field2,Field3,Field4,Field5,nn)
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
REAL*8         Field1(*),Field2(*),Field3(*),Field4(*),Field5(*)
INTEGER        iLevel,i,nn
INTEGER        :: iType= 2
CHARACTER*(20) :: cFF='Pressure'
CHARACTER*(60) :: cInFile

CALL ReadSol(cInFile,iLevel,iType,cFF,Field2,Field3,Field4,Field5)

IF (myid.ne.0) THEN
!  WRITE(*,*) Field2(1:nn)
 DO i=1,nn
  Field1(4*(i-1)+1) = Field2(i) 
  Field1(4*(i-1)+2) = Field3(i) 
  Field1(4*(i-1)+3) = Field4(i) 
  Field1(4*(i-1)+4) = Field5(i) 
 END DO
END IF

END SUBROUTINE ReadSol_Pres
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol_Time(cInFile)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,subnodes,RECVDD_myMPI,SENDDD_myMPI
interface
  SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
  USE def_FEAT
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  USE var_QuadScalar,ONLY:myDump
  IMPLICIT NONE
  INTEGER iType,iLevel
  REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
  CHARACTER cField*(20),cInFile*(60)
  END SUBROUTINE ReadSol
end interface
CHARACTER cInFile*(60)
INTEGER pID

 IF (myid.eq.showid) THEN
  WRITE(MTERM,*) 'Loading Time level from: "'//ADJUSTL(TRIM(cInFile))//'_Time.dmp','"'
 END IF

 IF (myid.eq.0) THEN
  OPEN(321,FILE=TRIM(ADJUSTL(cInFile))//'_Time.dmp')
  READ(321,*) TIMENS
  CLOSE(321)
 END IF

 IF (myid.NE.0) THEN
  CALL RECVDD_myMPI(timens,0)
 ELSE
  DO pID =1,subnodes
   CALL SENDDD_myMPI(timens,pID)
  END DO
 END IF

END SUBROUTINE ReadSol_Time
!
!-------------------------------------------------------------------------------
!
SUBROUTINE ReadSol(cInFile,iLevel,iType,cField,Field1,Field2,Field3,Field4)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
    RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
USE var_QuadScalar,ONLY:myDump,ViscoSc
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
INTEGER iInd,iType,iLevel
REAL*8, OPTIONAL :: Field1(*),Field2(*),Field3(*),Field4(*)
CHARACTER cField*(20),cInFile*(60)
!----------------------------------------------------
INTEGER i,ivt,jvt,jel,kel,iP,nLengthE,nLengthV
REAL*8,ALLOCATABLE :: Field(:,:),auxField(:,:)
REAL*8 dMaxNel
INTEGER pnel,pID
CHARACTER cFile*(40)

IF (myid.NE.0) THEN

 ILEV = NLMIN

 nLengthE = 8**(NLMAX+(iLevel-1)-1)
 nLengthV = (2**(NLMAX+(iLevel-1))+1)**3

!  WRITE(*,*)  'READ  :::',nLengthE,nLengthV
END IF

 IF (myid.ne.0) dMaxNel = DBLE(KNEL(NLMIN))
 CALL COMM_Maximum(dMaxNel)

 IF (myid.eq.showid) THEN
  WRITE(MTERM,*) 'Loading dumped '//TRIM(ADJUSTL(cField))//' solution from: "'//ADJUSTL(TRIM(cInFile))//'_'//TRIM(ADJUSTL(cField))//'.dmp','"'
 END IF

 IF (myid.eq.0) THEN
  OPEN(321,FILE=TRIM(ADJUSTL(cInFile))//'_'//TRIM(ADJUSTL(cField))//'.dmp')
 END IF

 IF (iType.eq.1) THEN
  IF (Present(Field1)) CALL DistributeVertField(Field1)
  IF (Present(Field2)) CALL DistributeVertField(Field2)
  IF (Present(Field3)) CALL DistributeVertField(Field3)
 END IF

 IF (iType.eq.2) THEN
  IF (Present(Field1)) CALL DistributeElemField(Field1)
  IF (Present(Field2)) CALL DistributeElemField(Field2)
  IF (Present(Field3)) CALL DistributeElemField(Field3)
  IF (Present(Field4)) CALL DistributeElemField(Field4)
 END IF

 IF (iType.eq.3) THEN
   CALL DistributeVertField(ViscoSc%val11)
   CALL DistributeVertField(ViscoSc%val22)
   CALL DistributeVertField(ViscoSc%val33)
   CALL DistributeVertField(ViscoSc%val12)
   CALL DistributeVertField(ViscoSc%val13)
   CALL DistributeVertField(ViscoSc%val23)
 END IF

 IF (myid.eq.0) THEN
  CLOSE(321)
 END IF

 CONTAINS
! -----------------------------------------------------------------
SUBROUTINE DistributeVertField(xField)
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthV,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 

  CALL RECVD_myMPI(Field,nLengthV*KNEL(NLMIN),0)

  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthV
    jvt = myDump%Vertices(IEL,ivt)
    xField(jvt) = Field(ivt,iel)
   END DO
  END DO
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthV,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthV,KNEL(NLMIN))) 
    ALLOCATE(auxField(nLengthV,pnel)) 

    DO iel=1,KNEL(NLMIN)
!    READ(321,'(<nLengthV>ES14.6)') Field(1:nLengthV,iel)
     READ(321,*) Field(1:nLengthV,iel)
    END DO
   END IF
   CALL RECVI_myMPI(pnel,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthV
     auxField(ivt,I) = Field(ivt,iel)
    END DO
   END DO
 
   CALL SENDD_myMPI(auxField,pnel*nLengthV,pID)

  END DO

  DEALLOCATE(auxField) 

 END IF

DEALLOCATE(Field) 

END SUBROUTINE DistributeVertField
! -----------------------------------------------------------------
SUBROUTINE DistributeElemField(xField)
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
REAL*8 xField(*)

 IF (myid.ne.0) THEN

  CALL SENDI_myMPI(nLengthE,0)
  CALL SENDI_myMPI(KNEL(NLMIN),0)
  
  ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 

  CALL RECVD_myMPI(Field,nLengthE*KNEL(NLMIN),0)

  DO iel = 1,KNEL(NLMIN)
   DO ivt=1,nLengthE
    jvt = myDump%Elements(IEL,ivt)
    xField(jvt) = Field(ivt,iel)
   END DO
  END DO
 
 ELSE

  DO pID =1,subnodes

   CALL RECVI_myMPI(nLengthE,pID)
   IF (pID.EQ.1) THEN
    pnel = INT(dMaxNel)
    ALLOCATE(Field(nLengthE,KNEL(NLMIN))) 
    ALLOCATE(auxField(nLengthE,pnel)) 
    DO iel=1,KNEL(NLMIN)
!    READ(321,'(<nLengthE>ES14.6)') Field(1:nLengthE,iel)
     READ(321,*) Field(1:nLengthE,iel)
    END DO
   END IF
   CALL RECVI_myMPI(pnel,pID)

   DO I=1,pnel
   IEL = coarse%pELEMLINK(pID,I)
    DO ivt=1,nLengthE
     auxField(ivt,I) = Field(ivt,iel)
    END DO
   END DO

   CALL SENDD_myMPI(auxField,pnel*nLengthE,pID)
  END DO
 
  DEALLOCATE(auxField) 

 END IF

DEALLOCATE(Field) 

END SUBROUTINE DistributeElemField
! -----------------------------------------------------------------

END SUBROUTINE ReadSol
!
!-------------------------------------------------------------------------------
!
SUBROUTINE Output_Profiles(iOutput)
USE def_FEAT
USE var_QuadScalar,ONLY:QuadSc,LinSc,&
    Viscosity,Distance,Distamce,mgNormShearStress
USE var_QuadScalar,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,showid,Comm_Summ,myMPI_Barrier
USE var_QuadScalar,ONLY:myExport,myFBM,mg_mesh
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
! USE PLinScalar,ONLY:PLinScP1toQ1,OutputInterphase,PLinLS,&
!                dNorm,IntPhaseElem,FracFieldQ1
use cinterface, only: outputRigidBodies 

IMPLICIT NONE
INTEGER iOutput,mfile

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

IF (myExport%Format.EQ."GMV") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_GMV_fields(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 END IF

ELSEIF (myExport%Format.EQ."VTK") THEN

 IF (myid.NE.0) THEN
  NLMAX = NLMAX + 1
  ILEV = myExport%Level
  CALL SETLEV(2)
  CALL Output_VTK_piece(iOutput,&
    mg_mesh%level(ILEV)%dcorvg,&
    mg_mesh%level(ILEV)%kvert)
  NLMAX = NLMAX - 1
 ELSE
  CALL Output_VTK_main(iOutput)
!   CALL Output_VTK_piece(iOutput,&
!     mg_mesh%level(ILEV)%dcorvg,&
!     mg_mesh%level(ILEV)%kvert)

  ILEV = NLMAX
  CALL SETLEV(2)
  CALL Output_Mesh(iOutput,"_vtk")
 END IF

END IF

IF (myid.eq.1) THEN
  if(outputRigidBodies())then
    call write_rigid_bodies(iOutput)
  end if
end if

END
!
! ----------------------------------------------
!
SUBROUTINE Output_Mesh(iO,cFolder)
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:mg_mesh,ilev,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
INTEGER i,j,iO,iBnds
CHARACTER cf*(24)
CHARACTER cFolder*(*)

WRITE(cf,'(A11,I3.3,A4)') ADJUSTL(TRIM(cFolder))//'/cMESH_',iO, '.tri'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') mg_mesh%level(ilev)%NEL,mg_mesh%level(ilev)%NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1,mg_mesh%level(ilev)%nvt
 WRITE(1,'(3ES13.5)') mg_mesh%level(ilev)%dcorvg(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1,mg_mesh%level(ilev)%nel
 WRITE(1,'(8I8)') mg_mesh%level(ILEV)%kvert(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
IF (allocated(myBoundary%bWall)) then
 DO i = 1,mg_mesh%level(ilev)%nvt
  IF (myBoundary%bWall(i)) THEN
   WRITE(1,'(I8)') 1
  ELSE
   WRITE(1,'(I8)') 0
  END IF
 END DO
ELSE
 DO i = 1,mg_mesh%level(ilev)%nvt
  WRITE(1,'(I8)') 0
 END DO
END IF

CLOSE(1)

OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cFolder))//'/file.prj')
WRITE(cf,'(A11,I3.3,A4)') 'cMESH_',iO, '.tri'
WRITE(2,'(A)') ADJUSTL(TRIM(cf))
 
DO iBnds = 1, nBnds
 cf = ' '
 WRITE(cf,'(A)') ADJUSTL(TRIM(cFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
 WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
 j=0
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   j = j + 1
  END IF
 END DO
 WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
 WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   WRITE(1,'(I8,A)') i
  END IF
 END DO
 CLOSE(1)
END DO
CLOSE(2)

END SUBROUTINE Output_Mesh
!
! ----------------------------------------------
!
SUBROUTINE Output_Profiles_sub(U,V,W,P,T,kv,dx,&
           iO,nvt,nel,time)
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
REAL*8 dx(3,*),U(*),V(*),W(*),P(*),T(*),time
REAL*8 DGP,DAUX
INTEGER kv(8,*),nvt,nel,IP
CHARACTER cf*(20),cmm*(20),cm*(20)
INTEGER :: mf=120,ivt,iel
INTEGER iO


  WRITE(cf(13:15),'(A3)') "000"
  WRITE(cf(16:20),'(I1,A4)') iO,".gmv"


 cf="_gmv/               "

 IF(myid.lt.10) THEN
  WRITE(cf(6:12),'(A4,I1,I1,A1)') 'res_',0,myid,'_'
 ELSE
  WRITE(cf(6:12),'(A4,I2,A1)') 'res_',myid,'_'
 END IF

 IF     ((iO.GE.0   ).AND.(iO.LT.10   )) THEN
  WRITE(cf(13:15),'(A3)') "000"
  WRITE(cf(16:20),'(I1,A4)') iO,".gmv"
 ELSEIF ((iO.GE.10  ).AND.(iO.LT.100  )) THEN
  WRITE(cf(13:14),'(A2)') "00"
  WRITE(cf(15:20),'(I2,A4)') iO,".gmv"
 ELSEIF ((iO.GE.100 ).AND.(iO.LT.1000 )) THEN
  WRITE(cf(13:13),'(A1)') "0"
  WRITE(cf(14:20),'(I3,A4)') iO,".gmv"
 ELSEIF ((iO.GE.1000).AND.(iO.LT.10000)) THEN
  WRITE(cf(13:20),'(I4,A4)') iO,".gmv"
 ELSEIF (iO.GE.10000)                       THEN
  STOP
 END IF

  IF(myid.eq.showid) WRITE(*,*) "Outputting gmv file into ",cf

 cmm="msh                 "
 IF(myid.lt.10) WRITE(cmm(4:10),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
 IF(myid.ge.10) WRITE(cmm(4:10),'(A1,I2,A4)') '_',myid,".gmv"

 cm="_gmv/msh             "
 IF(myid.lt.10) WRITE(cm(9:15),'(A1,I1,I1,A4)') '_',0,myid,".gmv"
 IF(myid.ge.10) WRITE(cm(9:15),'(A1,I2,A4)') '_',myid,".gmv"

 IF (iO.EQ.0) CALL XGMVMS(mf,cm,nel,nvt,kv,dx)

 OPEN (UNIT=mf,FILE=cf)
 WRITE(mf,'(A)')'gmvinput ascii'
 WRITE(UNIT=mf,FMT=*) 'nodes fromfile "',TRIM(cmm),'"'
 WRITE(UNIT=mf,FMT=*) 'cells fromfile "',TRIM(cmm),'"'

 WRITE(mf,*)  'velocity 1'
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(U(IVT))
 END DO

 DO IVT=1,NVT
  WRITE(mf,1000) REAL(V(IVT))
 END DO
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(W(IVT))
 END DO

 WRITE(mf,*)  'variable'
 WRITE(mf,*)  'pressure 1'
 DO IVT=1,NVT
  WRITE(mf,1000) REAL(P(IVT))
 END DO

 WRITE(mf,*)  'temperature 1'
 DO IVT=1,NVT


  WRITE(mf,1000) REAL(T(IVT))
 END DO

 WRITE(mf,*)  'endvars'
 WRITE(mf,*)  'probtime ',REAL(time)
 WRITE(mf,*)  'endgmv'
 CLOSE(mf)

1000  FORMAT(E12.5)
1100  FORMAT(8I8)

END
!
! ----------------------------------------------
!
SUBROUTINE Output_DUMPProfiles()
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE var_QuadScalar,ONLY:Tracer
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
CHARACTER COFile*15
INTEGER itwx,ifilen,i
DATA ifilen/0/

IF (myid.ne.0) THEN
 ILEV=NLMAX
 CALL SETLEV(2)

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='_ns/QL        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1,A1)') 0,itwx,'_'
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2,A1)') itwx,'_'
 IF (myid.eq.showid) &
 WRITE(*,*) "Outputting dump profiles into ",COFile
 IF     (myid.LT.10 ) THEN
  WRITE(COFile(10:12),'(A,I1)') "00",myid
 ELSEIF (myid.LT.100) THEN
  WRITE(COFile(10:12),'(A,I2)') "0",myid
 ELSE
  WRITE(COFile(10:12),'(I3)') myid
 END IF

 OPEN (UNIT=2,FILE=COFile)

 DO I=1,SIZE(QuadSc%ValU)
  WRITE(2,*) QuadSc%ValU(i)
 END DO
 DO I=1,SIZE(QuadSc%ValV)
  WRITE(2,*) QuadSc%ValV(i)
 END DO
 DO I=1,SIZE(QuadSc%ValW)
  WRITE(2,*) QuadSc%ValW(i)
 END DO

 DO I=1,SIZE(LinSc%ValP(NLMAX)%x)
  WRITE(2,*) LinSc%ValP(NLMAX)%x(i)
 END DO

 IF (bTracer) THEN
  DO I=1,SIZE(Tracer%Val(NLMAX+1)%x)
   WRITE(2,*) Tracer%Val(NLMAX+1)%x(i)
  END DO
 END IF
 WRITE(2,*) timens

 CLOSE(2)
END IF

END
!
! ----------------------------------------------
!
SUBROUTINE Load_DUMPProfiles(cdump)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,ShareValueD_myMPI
USE var_QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE var_QuadScalar,ONLY:Tracer
IMPLICIT NONE
CHARACTER CIFile*(30),cdump*(*)
INTEGER IFl,I
REAL*8 ddd
LOGICAL bExist
REAL*8 daux(1)

 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.ne.0) THEN
  IF (myid.eq.showid) &
  WRITE(*,*) "Loading initial dump profiles from ",CIFile

  INQUIRE(FILE=CIFile,EXIST=bExist)
  IF (.NOT.bExist) THEN
   WRITE(*,*) "File ",CIFile," could not be read ...",myid
   RETURN
  END IF

  OPEN (UNIT=1,FILE=CIFile) 

  DO I=1,SIZE(QuadSc%ValU)
   READ(1,*) QuadSc%ValU(i)
  END DO
  DO I=1,SIZE(QuadSc%ValV)
   READ(1,*) QuadSc%ValV(i)
  END DO
  DO I=1,SIZE(QuadSc%ValW)
   READ(1,*) QuadSc%ValW(i)
  END DO

  DO I=1,SIZE(LinSc%ValP(NLMAX)%x)
   READ(1,*) LinSc%ValP(NLMAX)%x(i)
  END DO

  IF (bTracer) THEN
   DO I=1,SIZE(Tracer%Val(NLMAX+1)%x)
    READ(1,*) Tracer%Val(NLMAX+1)%x(i)
   END DO
  END IF

  READ(1,*) timens
  CLOSE(1)
 END IF

 daux(1) = timens
 CALL ShareValueD_myMPI(daux,1,1)
 timens = daux(1)

END
!
! ----------------------------------------------
!
SUBROUTINE Load_LowDUMPProfiles(cdump)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,ShareValueD_myMPI
USE var_QuadScalar,ONLY:QuadSc,LinSc,bTracer
USE var_QuadScalar,ONLY:Tracer
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
CHARACTER CIFile*(30),cdump*(*)
INTEGER IFl,I,ndof
REAL*8 ddd
LOGICAL bExist
REAL*8 daux(1)

 CIFile=TRIM(ADJUSTL(cdump))
 IF (myid.ne.0) THEN
  ILEV = NLMAX-1
  ndof = KNVT(ILEV) + KNAT(ILEV) + KNET(ILEV) + KNEL(ILEV)
  IF (myid.eq.showid) &
  WRITE(*,*) "Loading initial dump profiles from ",CIFile

  INQUIRE(FILE=CIFile,EXIST=bExist)
  IF (.NOT.bExist) THEN
   WRITE(*,*) "File ",CIFile," could not be read ...",myid
   RETURN
  END IF

  OPEN (UNIT=1,FILE=CIFile) 

  DO I=1,ndof
   READ(1,*) QuadSc%ValU(i)
  END DO
  DO I=1,ndof
   READ(1,*) QuadSc%ValV(i)
  END DO
  DO I=1,ndof
   READ(1,*) QuadSc%ValW(i)
  END DO

  DO I=1,SIZE(LinSc%ValP(ILEV)%x)
   READ(1,*) LinSc%ValP(ILEV)%x(i)
  END DO

  READ(1,*) timens
  CLOSE(1)
 END IF

 daux(1) = timens
 CALL ShareValueD_myMPI(daux,1,1)
 timens = daux(1)

END
!
! ----------------------------------------------
!
SUBROUTINE FBM_ToFile()
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
CHARACTER COFile*15
INTEGER itwx,ifilen,i
INTEGER iP
DATA ifilen/0/

IF (myid.eq.showid) THEN

 ifilen=ifilen+1
 itwx=MOD(ifilen+insavn-1,insavn)+1
 COFile='_ns/PT        '
 IF (itwx.lt.10) WRITE(COFile(7:9),'(I1,I1)') 0,itwx
 IF (itwx.ge.10) WRITE(COFile(7:9),'(I2)') itwx
 WRITE(*,*) "Outputting particle data into ",COFile

 OPEN (UNIT=2,FILE=COFile)

 WRITE(2,'(I10)') myFBM%nParticles
 DO iP = 1,myFBM%nParticles
  WRITE(2,'(A)')    myFBM%ParticleNew(iP)%cType
  WRITE(2,'(D12.4)')  myFBM%ParticleNew(iP)%density
  WRITE(2,'(D12.4)')  myFBM%ParticleNew(iP)%sizes(1)
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Position
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Velocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Angle
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%AngularVelocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%Acceleration
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%FrameVelocity
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%ResistanceForce
  WRITE(2,'(3D12.4)') myFBM%ParticleNew(iP)%TorqueForce
 END DO

 CLOSE(2)

END IF

END SUBROUTINE FBM_ToFile
!
! ----------------------------------------------
!
SUBROUTINE FBM_FromFile(CIFile)
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myFBM
 INTEGER iP
 CHARACTER*(*) CIFile

 IF (myid.eq.showid) &
 WRITE(*,*) "Loading particle data from ",CIFile

 OPEN(1,FILE=TRIM(ADJUSTL(CIFile)))

 READ(1,'(I10)') myFBM%nParticles
 ALLOCATE(myFBM%ParticleNew(myFBM%nParticles),myFBM%ParticleOld(myFBM%nParticles))
 ALLOCATE(myFBM%Force(6*myFBM%nParticles))
 DO iP = 1,myFBM%nParticles
  READ(1,'(A)')    myFBM%ParticleOld(iP)%cType
  myFBM%ParticleOld(iP)%cType = TRIM(ADJUSTL(myFBM%ParticleOld(iP)%cType))
  READ(1,'(D12.4)')  myFBM%ParticleOld(iP)%density
  READ(1,'(D12.4)')  myFBM%ParticleOld(iP)%sizes(1)
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Position
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Velocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Angle
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%AngularVelocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%Acceleration
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%FrameVelocity
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%ResistanceForce
  READ(1,'(3D12.4)') myFBM%ParticleOld(iP)%TorqueForce
 END DO

 myFBM%ParticleNew = myFBM%ParticleOld
 CLOSE(1)

END SUBROUTINE FBM_FromFile
!
!---------------------------------------------------------------------------
!
SUBROUTINE CreateDumpStructures(iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myDump,mg_mesh

IMPLICIT NONE

INTEGER iLevel
INTEGER JEL,KEL,ivt,jvt,nLengthE,nLengthV
INTEGER iaux,jaux,jj,kv(8),II
LOGICAL ,ALLOCATABLE :: bGot(:)
! -------------- workspace -------------------

! IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN
 NEL  = mg_mesh%level(ilev)%nel
 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3
 IF(ALLOCATED(myDump%Elements)) DEALLOCATE(myDump%Elements)
 IF(ALLOCATED(myDump%Vertices)) DEALLOCATE(myDump%Vertices)
 ALLOCATE(myDump%Elements(NEL,nLengthE))
 ALLOCATE(myDump%Vertices(NEL,nLengthV))

 DO IEL = 1,mg_mesh%level(NLMIN)%nel

  myDump%Elements(IEL,1) = IEL
  iaux = 1
  DO II=1+1,NLMAX
   jaux = iaux
   DO jel=1,jaux
    kel = myDump%Elements(iel,jel)
    CALL Get8Elem(mg_mesh%level(II)%kadj,&
      kv,kel)
    DO jj = 2,8
     iaux = iaux + 1
     myDump%Elements(iel,iaux) = kv(jj)
    END DO
   END DO
  END DO
 END DO 

 ALLOCATE(bGot(mg_mesh%level(NLMAX)%nvt))
 DO IEL = 1,mg_mesh%level(NLMIN)%nel
  bGot = .FALSE.
  iaux = 0
  DO JEL = 1,nLengthE
   KEL = myDump%Elements(IEL,JEL)
   
   CALL getVert(mg_mesh%level(NLMAX)%kvert,&
                kv,KEL)

   DO IVT = 1,8
    JVT = kv(IVT)
    IF (.NOT.bGot(JVT)) THEN
     iaux = iaux + 1
     myDump%Vertices(IEL,iaux) = JVT
     bGot(JVT) = .TRUE.
    END IF
   END DO

  END DO

!  IF (iaux.ne.729) WRITE(*,*) myid,iel,iaux
  
 END DO
 DEALLOCATE(bGot)

 NLMAX = NLMAX - iLevel
!   write(*,*) size(myDump%Vertices),"asdas dasd sad sa",myid
!   pause

! END IF


 CONTAINS

 SUBROUTINE Get8Elem(KADJ,k,el)
 INTEGER KADJ(6,*),k(8),el

 k(1) = el
 k(2) = KADJ(3,k(1))
 k(3) = KADJ(3,k(2))
 k(4) = KADJ(3,k(3))
 k(5) = KADJ(6,k(1))
 k(6) = KADJ(3,k(5))
 k(7) = KADJ(3,k(6))
 k(8) = KADJ(3,k(7))

 END SUBROUTINE Get8Elem

 SUBROUTINE getVert(BigKv,SmallKv,elem)
 INTEGER BigKv(8,*),SmallKv(8),elem

 SmallKv(:) = BigKv(:,elem)

 END SUBROUTINE getVert

END SUBROUTINE CreateDumpStructures
!
!---------------------------------------------------------------------------
!
SUBROUTINE ReleaseSmartDumpFiles(iOutO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier
USE var_QuadScalar,ONLY:QuadSc,LinSc,myDump
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
CHARACTER cOutFile*20,command*100
INTEGER iOutO
INTEGER iOut,myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lCmnd
INTEGER ifilen,i
DATA ifilen/0/

IF (myid.NE.0) THEN

 NLMAX = NLMAX+1

 IF (iOutO.LT.0) THEN
  ifilen=ifilen+1
  iOut=MOD(ifilen+insavn-1,insavn)+1
 ELSE
  iOut = iOuto
 END IF

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 IF (myid.eq.1) THEN

  command = 'mkdir _dump'
  WRITE(*,*) 'Executing command:"',TRIM(ADJUSTL(command)),'"'
  CALL system(TRIM(ADJUSTL(command)))
  lCmnd = LEN(TRIM(ADJUSTL(command)))
  IF     (iOut.LT.10    ) THEN
   WRITE(command(lCmnd+1:lCmnd+3),'(A2,I1)') '/0',iOut
  ELSEIF (iOut.LT.100   ) THEN
   WRITE(command(lCmnd+1:lCmnd+3),'(A1,I2)') '/',iOut
  ELSE
    WRITE(*,*) "Decrease the output index!"
    STOP
  END IF
  WRITE(*,*) 'Executing command:"',TRIM(ADJUSTL(command)),'"'
  CALL system(TRIM(ADJUSTL(command)))

 END IF

! pause
 CALL myMPI_Barrier()

 cOutFile = '_dump/00/******.prf'

 IF     (iOut.LT.10    ) THEN
  WRITE(cOutFile(7:8),'(A1,I1)') '0',iOut
 ELSEIF (iOut.LT.100   ) THEN
   WRITE(cOutFile(7:8),'(I2)') iOut
 ELSE
   WRITE(*,*) "Decrease the output index!"
   STOP
 END IF

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Storing the ",TRIM(ADJUSTL(cOutFile))," series for solution backup ..."
 END IF

 DO IEL = 1,KNEL(NLMIN)
  myCoarseElem = coarse%myELEMLINK(IEL)
  IF     (myCoarseElem.LT.10    ) THEN
   WRITE(cOutFile(10:15),'(A5,I1)') '00000',myCoarseElem
  ELSEIF (myCoarseElem.LT.100   ) THEN
   WRITE(cOutFile(10:15),'(A4,I2)') '0000',myCoarseElem
  ELSEIF (myCoarseElem.LT.1000  ) THEN
   WRITE(cOutFile(10:15),'(A3,I3)') '000',myCoarseElem
  ELSEIF (myCoarseElem.LT.10000 ) THEN
   WRITE(cOutFile(10:15),'(A2,I4)') '00',myCoarseElem
  ELSEIF (myCoarseElem.LT.100000) THEN
   WRITE(cOutFile(10:15),'(A1,I5)') '0',myCoarseElem
  ELSE
   WRITE(*,*) "Decrease the problem size!"
   STOP
  END IF

  ! Velocities
  OPEN (FILE=TRIM(ADJUSTL(cOutFile)),UNIT = 547)
  WRITE(547,*)  "Velocities"
  DO ivt=1,nLengthV
   jvt = myDump%Vertices(IEL,ivt)
   WRITE(547,'(3D16.8)') QuadSc%valU(jvt),QuadSc%valV(jvt),QuadSc%valW(jvt)
  END DO

  ! Pressure
  WRITE(547,*)  "Pressure"
  DO jel=1,nLengthE
   kel = myDump%Elements(IEL,jel)
   IF (kel.le.knel(NLMAX-1)) THEN
    iP = 4*(kel-1)
   WRITE(547,'(4D16.8)') LinSc%valP(NLMAX-1)%x(iP+1:iP+4)
   END IF   
  END DO

  WRITE(547,*)  "EOF"
  CLOSE(547)

 END DO 

 IF (myid.eq.1) THEN
  cOutFile = '_dump/00/time.prf'
  IF     (iOut.LT.10    ) THEN
   WRITE(cOutFile(7:8),'(A1,I1)') '0',iOut
  ELSEIF (iOut.LT.100   ) THEN
    WRITE(cOutFile(7:8),'(I2)') iOut
  ELSE
    WRITE(*,*) "Decrease the output index!"
    STOP
  END IF

  OPEN (FILE=TRIM(ADJUSTL(cOutFile)),UNIT = 547)
  WRITE (547,*) "timens"
  WRITE (547,'(D16.8)') timens
  CLOSE(547)
 END IF

 NLMAX = NLMAX-1

END IF

END SUBROUTINE ReleaseSmartDumpFiles
!
!---------------------------------------------------------------------------
!
SUBROUTINE LoadSmartDumpFiles(cFldrInFile,iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse
USE var_QuadScalar,ONLY:QuadSc,LinSc,myDump
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
IMPLICIT NONE
INTEGER iLevel
CHARACTER cInFile*99,cFldrInFile*(*)
INTEGER myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lStr

IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cFldrInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Reading the ",TRIM(ADJUSTL(cInFile))," series for initialization of profiles ..."
 END IF

 DO IEL = 1,KNEL(NLMIN)
  myCoarseElem = coarse%myELEMLINK(IEL)
  IF     (myCoarseElem.LT.10    ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A6,I1,A4)') '/00000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100   ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A5,I2,A4)') '/0000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.1000  ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A4,I3,A4)') '/000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.10000 ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A3,I4,A4)') '/00',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100000) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A2,I5,A4)') '/0',myCoarseElem,".prf"
  ELSE
   WRITE(*,*) "Decrease the problem size!"
   STOP
  END IF

  ! Velocities
  OPEN (FILE=TRIM(ADJUSTL(cInFile)),UNIT = 547)
  READ(547,*)  
  DO ivt=1,nLengthV
   jvt = myDump%Vertices(IEL,ivt)
   READ(547,*) QuadSc%valU(jvt),QuadSc%valV(jvt),QuadSc%valW(jvt)
  END DO

  ! Pressure
  READ(547,*)  
  DO jel=1,nLengthE
   kel = myDump%Elements(IEL,jel)
   IF (kel.le.knel(NLMAX-1)) THEN
    iP = 4*(kel-1)
   READ(547,*) LinSc%valP(NLMAX-1)%x(iP+1:iP+4)
   END IF   
  END DO

  CLOSE(547)

 END DO 
END IF

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cFldrInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))
 WRITE(cInFile(lStr+1:lStr+9),'(A)') '/time.prf'

 OPEN (FILE=TRIM(ADJUSTL(cInFile)),UNIT = 547)
 READ (547,*) 
 READ (547,*) timens
 CLOSE(547)

IF (myid.NE.0) THEN

 NLMAX = NLMAX-iLevel

END IF

END SUBROUTINE LoadSmartDumpFiles
!
!---------------------------------------------------------------------------
!
SUBROUTINE LoadSmartAdaptedMeshFile(dcorvg,cMeshInFile,iLevel)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,coarse
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:myDump
IMPLICIT NONE
REAL*8 dcorvg(3,*)
INTEGER iLevel
CHARACTER cInFile*99,cMeshInFile*(*)
INTEGER myCoarseElem,ivt,jvt,jel,kel,nLengthE,nLengthV
INTEGER iP,lStr

IF (myid.NE.0) THEN

 NLMAX = NLMAX+iLevel

 ILEV = NLMIN

 nLengthE = 8**(NLMAX-1)
 nLengthV = (2**(NLMAX-1)+1)**3

 cInFile = ' '
 WRITE(cInFile(1:),'(A)') TRIM(ADJUSTL(cMeshInFile))
 lStr = LEN(TRIM(ADJUSTL(cInFile)))

 IF (myid.EQ.1) THEN
  WRITE(*,'(3A)') "Reading the ",TRIM(ADJUSTL(cInFile))," series for preadaptation of the mesh ..."
 END IF

 DO IEL = 1,KNEL(NLMIN)
  myCoarseElem = coarse%myELEMLINK(IEL)
  IF     (myCoarseElem.LT.10    ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A6,I1,A4)') '/00000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100   ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A5,I2,A4)') '/0000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.1000  ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A4,I3,A4)') '/000',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.10000 ) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A3,I4,A4)') '/00',myCoarseElem,".prf"
  ELSEIF (myCoarseElem.LT.100000) THEN
   WRITE(cInFile(lStr+1:lStr+11),'(A2,I5,A4)') '/0',myCoarseElem,".prf"
  ELSE
   WRITE(*,*) "Decrease the problem size!"
   STOP
  END IF

  ! Velocities
  OPEN (FILE=TRIM(ADJUSTL(cInFile)),UNIT = 547)
  READ(547,*)
  DO ivt=1,nLengthV
   jvt = myDump%Vertices(IEL,ivt)
   READ(547,*) dcorvg(:,jvt)
  END DO

  CLOSE(547)

 END DO 

 NLMAX = NLMAX-iLevel

ELSE
  OPEN (FILE=TRIM(ADJUSTL(cMeshInFile))//"/coarse.prf",UNIT = 547)
  DO ivt = 1, NVT
   READ(547,*) dcorvg(:,ivt)
  END DO
  CLOSE (547)
END IF

END SUBROUTINE LoadSmartAdaptedMeshFile


SUBROUTINE Output_VTK_piece(iO,dcoor,kvert)
USE def_FEAT
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY: QuadSc,LinSc,Viscosity,Distance,Distamce,mgNormShearStress,myALE
USE var_QuadScalar,ONLY: MixerKnpr,FictKNPR,ViscoSc,myBoundary
USE var_QuadScalar,ONLY: iTemperature_AVG,Temperature_AVG,Temperature
USE var_QuadScalar,ONLY: Tracer
USE var_QuadScalar,ONLY: myExport, Properties, bViscoElastic,myFBM,mg_mesh,Shearrate,myHeatObjects,MaterialDistribution
USE var_QuadScalar,ONLY: myFBM,knvt,knet,knat,knel,ElemSizeDist,BoundaryNormal
USE var_QuadScalar,ONLY: GenLinScalar
USE def_LinScalar, ONLY: mg_RhoCoeff,mg_CpCoeff,mg_LambdaCoeff

IMPLICIT NONE
REAL*8 dcoor(3,*)
INTEGER kvert(8,*),iO,ioffset,ive,ivt,iField,i,istat
CHARACTER fileid*(5),filename*(27),procid*(3),cGenScalar*(50)
INTEGER NoOfElem,NoOfVert
REAL*8,ALLOCATABLE ::  tau(:,:)
REAL*8 psi(6)
REAL*8 Temp,dMF
integer :: iunit = 908070
integer iX, ifbm,iFld

NoOfElem = KNEL(ILEV)
NoOfVert = KNVT(ILEV)

filename=" "
WRITE(filename(1:),'(A,I5.5,A4)') '_vtk/res_node_***.',iO,".vtu"

IF(myid.eq.showid) WRITE(*,'(104("="))') 
IF(myid.eq.showid) WRITE(*,*) "Outputting vtk file into ",filename
WRITE(filename(15:17),'(I3.3)') myid

OPEN (UNIT=iunit,FILE=filename,action='write',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open file for writing. "
  stop          
end if

write(iunit, *)"<VTKFile type=""UnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(iunit, *)"  <UnstructuredGrid>"
write(iunit, *)"    <Piece NumberOfPoints=""",KNVT(ILEV),""" NumberOfCells=""",NoOfElem,""">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <PointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValU(ivt)),&
     REAL(QuadSc%ValV(ivt)),&
     REAL(QuadSc%ValW(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  
 CASE('GradVelocity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity_x",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValUx(ivt)),REAL(QuadSc%ValUy(ivt)),REAL(QuadSc%ValUz(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity_y",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValVx(ivt)),REAL(QuadSc%ValVy(ivt)),REAL(QuadSc%ValVz(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Velocity_z",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(QuadSc%ValWx(ivt)),REAL(QuadSc%ValWy(ivt)),REAL(QuadSc%ValWz(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  
 CASE('PartForce')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","PartForce",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   if((FictKNPR(ivt) .ne. 0).and.(FictKNPR(ivt).le.myFBM%nParticles))then
     ifbm = 1
     write(iunit, '(A,3E16.7)')"        ",&
     REAL(real(ifbm) * myFBM%particleNew(FictKNPR(ivt))%ResistanceForce(1)),&
     REAL(real(ifbm) * myFBM%particleNew(FictKNPR(ivt))%ResistanceForce(2)),&
     REAL(real(ifbm) * myFBM%particleNew(FictKNPR(ivt))%ResistanceForce(3))
   else
     ifbm = 0
     write(iunit, '(A,3E16.7)')"        ",&
     REAL(0),&
     REAL(0),&
     REAL(0)
   end if
  end do
  write(iunit, *)"        </DataArray>"
   
 CASE('Psi')
  if(bViscoElastic)then

  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Psi",""" NumberOfComponents=""6"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Val11(ivt)),REAL(ViscoSc%Val22(ivt)),REAL(ViscoSc%Val33(ivt)),&
                                        REAL(ViscoSc%Val12(ivt)),REAL(ViscoSc%Val13(ivt)),REAL(ViscoSc%Val23(ivt)) 
  end do

  write(iunit, *)"        </DataArray>"
  end if

 CASE('Stress')
  if(bViscoElastic)then

  ALLOCATE(tau(6,NoOfVert))
  DO i=1,NoOfVert
   psi = [ViscoSc%Val11(i),ViscoSc%Val22(i),ViscoSc%Val33(i),&
          ViscoSc%Val12(i),ViscoSc%Val13(i),ViscoSc%Val23(i)]
   CALL ConvertPsiToTau(psi,tau(:,i))   
!      tau(1,i) = (tau(1,i) - 1d0)/Properties%ViscoLambda
!      tau(2,i) = (tau(2,i) - 1d0)/Properties%ViscoLambda
!      tau(3,i) = (tau(3,i) - 1d0)/Properties%ViscoLambda
!      tau(4,i) = (tau(4,i) - 0d0)/Properties%ViscoLambda
!      tau(5,i) = (tau(5,i) - 0d0)/Properties%ViscoLambda
!      tau(6,i) = (tau(6,i) - 0d0)/Properties%ViscoLambda
  END DO

!  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6"" format=""ascii"">"
!   do ivt=1,NoOfVert
!    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Val11(ivt)),REAL(ViscoSc%Val22(ivt)),REAL(ViscoSc%Val33(ivt)),REAL(ViscoSc%Val12(ivt)),REAL(ViscoSc%Val13(ivt)),REAL(ViscoSc%Val23(ivt))
!   end do

  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,6E16.7)')"        ",REAL(tau(1,ivt)),REAL(tau(2,ivt)),REAL(tau(3,ivt)),REAL(tau(4,ivt)),REAL(tau(5,ivt)),REAL(tau(6,ivt))
  end do

  DEALLOCATE(tau)
  write(iunit, *)"        </DataArray>"
  end if
 
 CASE('GradStress')
  if(bViscoElastic)then
    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_11",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad11%x(ivt)),REAL(ViscoSc%Grad11%y(ivt)),REAL(ViscoSc%Grad11%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"

    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_22",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad22%x(ivt)),REAL(ViscoSc%Grad22%y(ivt)),REAL(ViscoSc%Grad22%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"

    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_33",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad33%x(ivt)),REAL(ViscoSc%Grad33%y(ivt)),REAL(ViscoSc%Grad33%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"

    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_12",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad12%x(ivt)),REAL(ViscoSc%Grad12%y(ivt)),REAL(ViscoSc%Grad12%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"

    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_13",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad13%x(ivt)),REAL(ViscoSc%Grad13%y(ivt)),REAL(ViscoSc%Grad13%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"

    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","GradStress_23",""" NumberOfComponents=""3"" format=""ascii"">"
    do ivt=1,NoOfVert
    write(iunit, '(A,6E16.7)')"        ",REAL(ViscoSc%Grad23%x(ivt)),REAL(ViscoSc%Grad23%y(ivt)),REAL(ViscoSc%Grad23%z(ivt))
    end do
    write(iunit, *)"        </DataArray>"
  end if
  
 CASE('GenScalar')
  DO iFld=1,GenLinScalar%nOfFields
  WRITE(cGenScalar,'(A)') TRIM(GenLinScalar%Fld(iFld)%cName)
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""",ADJUSTL(TRIM(cGenScalar)),""" format=""ascii"">"
   do ivt=1,NoOfVert
    write(iunit, '(A,E16.7)')"        ",REAL(GenLinScalar%Fld(iFld)%Val(ivt))
   end do
   write(iunit, *)"        </DataArray>"
 
  END DO

 CASE('MeshVelo')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(myALE%MeshVelo(1,ivt)),REAL(myALE%MeshVelo(2,ivt)),REAL(myALE%MeshVelo(3,ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('BoundaryNormal')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","BoundaryNormal",""" NumberOfComponents=""3"" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,3E16.7)')"        ",REAL(BoundaryNormal(1,ivt)),REAL(BoundaryNormal(2,ivt)),REAL(BoundaryNormal(3,ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Pressure_V')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_V",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(LinSc%Q2(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Temperature')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Temperature(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Temperature_AVG')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Temperature_AVG",""" format=""ascii"">"
  if (DBLE(iTemperature_AVG).gt.0) then
   do ivt=1,NoOfVert
    write(iunit, '(A,E16.7)')"        ",REAL(Temperature_AVG(ivt)/DBLE(iTemperature_AVG))
   end do
  else
   do ivt=1,NoOfVert
    write(iunit, '(A,E16.7)')"        ",0d0
   end do
  end if
  write(iunit, *)"        </DataArray>"
  
 CASE('Shearrate')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Shearrate",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(10d0**Shearrate(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Wall')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Wall",""" format=""ascii"">"
  do ivt=1,NoOfVert
   IF (myBoundary%bWall(ivt)) THEN
    write(iunit, '(A,E16.7)')"        ",1d0
   ELSE
    write(iunit, '(A,E16.7)')"        ",0d0
   END IF
  end do
  write(iunit, *)"        </DataArray>"

 CASE('LogShearrate')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","LogShearrate",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Shearrate(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('HeatObjects')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Block",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myHeatObjects%Block(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Wire",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myHeatObjects%Wire(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Melt",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",myHeatObjects%Channel(ivt)
  end do
  write(iunit, *)"        </DataArray>"

CASE('ElemSizeDist')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","ElemSizeDist",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(ElemSizeDist(ivt))
  end do
  write(iunit, *)"        </DataArray>"
  
 CASE('Distamce')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Shell",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Distance(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Mixer')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Mixer",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(FictKNPR(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Viscosity')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Viscosity(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Monitor')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Monitor",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(myALE%monitor(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('Distance')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Distance",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(Distance(ivt))
  end do
  write(iunit, *)"        </DataArray>"
 CASE('Shell')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Shell",""" format=""ascii"">"
  do ivt=1,NoOfVert
   write(iunit, '(A,E16.7)')"        ",REAL(QuadSc%AuxV(ivt))
  end do
  write(iunit, *)"        </DataArray>"

 CASE('BndryType')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","BndryType",""" format=""ascii"">"
  do ivt=1,NoOfVert
   iX = 5
   IF (mg_mesh%BndryNodes(ivt)%ParamTypes(1)) iX = Min(1,iX)
   IF (mg_mesh%BndryNodes(ivt)%ParamTypes(2)) iX = Min(2,iX)
   IF (mg_mesh%BndryNodes(ivt)%ParamTypes(3)) iX = Min(3,iX)
   IF (mg_mesh%BndryNodes(ivt)%ParamTypes(4)) iX = Min(4,iX)
   iF (iX.eq.5) iX = 0
   write(iunit, '(A,E16.7)')"        ",REAL(iX)
  end do
  write(iunit, *)"        </DataArray>"

 CASE('OuterPoint')
  write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","OuterPoint",""" format=""ascii"">"
  do ivt=1,NoOfVert
   iX = 0
   IF (mg_mesh%BndryNodes(ivt)%bOuterPoint) iX=1
   write(iunit, '(A,E16.7)')"        ",REAL(iX)
  end do
  write(iunit, *)"        </DataArray>"
 END SELECT 

END DO

write(iunit, '(A)')"    </PointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"    <CellData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
!  CASE('Alpha_E')
!   IF (ILEV.EQ.NLMAX-1) THEN
!    write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Alpha_E",""" format=""ascii"">"
!    do ivt=1,NoOfElem
!     write(iunit, '(A,E16.7)')"        ",REAL(mgDiffCoeff(NLMAX)%x(ivt))
!    end do
!    write(iunit, *)"        </DataArray>"
!   END IF
  
 CASE('Pressure_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Pressure_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    ive = 4*(ivt-1)+1
    write(iunit, '(A,E16.7)')"        ",REAL(LinSc%ValP(NLMAX-1)%x(ive))
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 CASE('MatProps_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MF_E_[kg/m3]",""" format=""ascii"">"
   do ivt=KNVT(ILEV+1)-NoOfElem+1,KNVT(ILEV+1)
    Temp = GenLinScalar%Fld(1)%Val(ivt)
    CALL MeltFunction_MF(dMF,Temp)
    write(iunit, '(A,E16.7)')"        ",REAL(dMF)
   end do
   write(iunit, *)"        </DataArray>"
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Rho_E_[kg/m3]",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(mg_RhoCoeff(NLMAX-1)%x(ivt))/1d-9
   end do
   write(iunit, *)"        </DataArray>"
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Cp_E_[J/kg/K]",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(mg_CpCoeff(NLMAX-1)%x(ivt))/1d+4
   end do
   write(iunit, *)"        </DataArray>"
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Lambda_E_[W/m/K]",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(mg_LambdaCoeff(NLMAX-1)%x(ivt))/1d-1
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 CASE('Viscosity_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Viscosity_E",""" format=""ascii"">"
   do ivt=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(Viscosity((nvt+net+nat+ivt)))
   end do
   write(iunit, *)"        </DataArray>"
  END IF

 CASE('Material_E')
  IF (ALLOCATED(MaterialDistribution)) THEN
  IF (ILEV.LE.NLMAX-1.and.ALLOCATED(MaterialDistribution(ILEV)%x)) THEN
   write(iunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","Material_E",""" format=""ascii"">"
   do iel=1,NoOfElem
    write(iunit, '(A,E16.7)')"        ",REAL(MaterialDistribution(ILEV)%x(iel))
   end do
   write(iunit, *)"        </DataArray>"
  END IF
  END IF
 END SELECT

END DO

write(iunit, '(A)')"    </CellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(iunit, '(A)')"      <Points>"
write(iunit, '(A)')"        <DataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"" RangeMin=""0"" RangeMax=""1.0"">"
do ivt=1,NoOfVert
 write(iunit,'(A10,3E16.7)')"          ",REAL(dcoor(1,ivt)),REAL(dcoor(2,ivt)),REAL(dcoor(3,ivt))
!  write(iunit,'(A10,3E16.7)')"          ",REAL(myALE%OldCoor(1,ivt)),REAL(myALE%OldCoor(2,ivt)),REAL(myALE%OldCoor(3,ivt))
end do
write(iunit, *)"        </DataArray>"
write(iunit, *)"      </Points>"

write(iunit, *)"      <Cells>"
write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""connectivity"" format=""ascii"" RangeMin=""0"" RangeMax=""",NoOfElem-1,""">"
do ive=1,NoOfElem   
 write(iunit, '(8I10)')kvert(1,ive)-1,kvert(2,ive)-1,kvert(3,ive)-1,kvert(4,ive)-1,&
                       kvert(5,ive)-1,kvert(6,ive)-1,kvert(7,ive)-1,kvert(8,ive)-1
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A,I10,A)')"        <DataArray type=""Int32"" Name=""offsets"" format=""ascii"" RangeMin=""8"" RangeMax=""",8*NoOfElem,""">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')ive*8,(ive+1)*8,(ive+2)*8,(ive+3)*8,(ive+4)*8,(ive+5)*8,(ive+6)*8,(ive+7)*8
end do

do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')ive*8
end do
write(iunit, '(A)')"        </DataArray>"

write(iunit, '(A)')"        <DataArray type=""UInt8"" Name=""types"" format=""ascii"" RangeMin=""12"" RangeMax=""12"">"
ioffset=NoOfElem/8
ioffset=ioffset*8
do ive=1,ioffset,8
 write(iunit, '(8I10)')12,12,12,12,12,12,12,12
end do
do ive=ioffset+1,NoOfElem
 write(iunit, '(I10)')12
end do
write(iunit, '(A)')"        </DataArray>"
 
write(iunit, *)"      </Cells>"
write(iunit, *)"    </Piece>"
   
write(iunit, *)"  </UnstructuredGrid>"
write(iunit, *)"</VTKFile>"
close(iunit)

END SUBROUTINE Output_VTK_piece



SUBROUTINE Output_VTK_main(iO)
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY:myExport,bViscoElastic,MaterialDistribution
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel
USE var_QuadScalar,ONLY:GenLinScalar
USE def_FEAT

IMPLICIT NONE
INTEGER iO,iproc,iField,iFld
INTEGER :: iMainUnit=555
CHARACTER mainname*(20) 
CHARACTER filename*(26)
CHARACTER cGenScalar*(50)

integer :: istat

! generate the file name
mainname=' '
WRITE(mainname(1:),'(A,I5.5,A5)') '_vtk/main.',iO,'.pvtu'

OPEN (UNIT=imainunit,FILE=mainname,action='write',iostat=istat)
if(istat .ne. 0)then
  write(*,*)"Could not open file for writing. "
  stop          
end if

write(imainunit, *)"<VTKFile type=""PUnstructuredGrid"" version=""0.1"" byte_order=""LittleEndian"">"
write(imainunit, *)"  <PUnstructuredGrid GhostLevel=""0"">"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the node field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PPointData>"

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity",""" NumberOfComponents=""3""/>"
 CASE('GradVelocity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity_x",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity_y",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Velocity_z",""" NumberOfComponents=""3""/>"
 CASE('PartForce')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","PartForce",""" NumberOfComponents=""3""/>"
 CASE('Psi')
  if(bViscoElastic)then
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Psi",""" NumberOfComponents=""6""/>"
  end if
 CASE('Stress')
  if(bViscoElastic)then
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Stress",""" NumberOfComponents=""6""/>"
  end if
 CASE('GradStress')
  if(bViscoElastic)then
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_11",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_22",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_33",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_12",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_13",""" NumberOfComponents=""3""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","GradStress_23",""" NumberOfComponents=""3""/>"
  end if
 CASE('GenScalar')
  DO iFld=1,GenLinScalar%nOfFields
  WRITE(cGenScalar,'(A)') TRIM(GenLinScalar%Fld(iFld)%cName)
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""",ADJUSTL(TRIM(cGenScalar)),"""/>"
 END DO
 CASE('MeshVelo')
  write(imainunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","MeshVelocity",""" NumberOfComponents=""3""/>"
 CASE('BoundaryNormal')
  write(imainunit, '(A,A,A)')"        <DataArray type=""Float32"" Name=""","BoundaryNormal",""" NumberOfComponents=""3""/>"
 CASE('Pressure_V')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_V","""/>"
 CASE('Temperature')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature","""/>"
 CASE('Temperature_AVG')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Temperature_AVG","""/>"
 CASE('Shearrate')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Shearrate","""/>"
 CASE('Wall')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Wall","""/>"
 CASE('LogShearrate')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","LogShearrate","""/>"
 CASE('HeatObjects')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Block","""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Wire","""/>"
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Melt","""/>"
 CASE('ElemSizeDist')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","ElemSizeDist","""/>"
 CASE('Distamce')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Shell","""/>"
 CASE('Mixer')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Mixer","""/>"
 CASE('Viscosity')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity","""/>"
 CASE('Monitor')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Monitor","""/>"
 CASE('Distance')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Distance","""/>"
 CASE('Shell')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Shell","""/>"
 CASE('BndryType')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","BndryType","""/>"
 CASE('OuterPoint')
  write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","OuterPoint","""/>"

 END SELECT
END DO

write(imainunit, '(A)')"    </PPointData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the cell field data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, '(A)')"    <PCellData>"
DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Pressure_E","""/>"
  END IF

CASE('MatProps_E')
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","MF_E_[-]","""/>"
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Rho_E_[kg/m3]","""/>"
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Cp_E_[J/kg/K]","""/>"
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Lambda_E_[W/m/K]","""/>"
  END IF

CASE('Viscosity_E')
!  WRITE(*,*) myExport%Level,myExport%LevelMax,myExport%Level.EQ.myExport%LevelMax
  IF (myExport%Level.EQ.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Viscosity_E","""/>"
  END IF

 CASE('Material_E')
  IF (ALLOCATED(MaterialDistribution)) THEN
  WRITE(*,*) "myExport%Level.LE.myExport%LevelMax  ",myExport%Level,myExport%LevelMax
  IF (myExport%Level.LE.myExport%LevelMax) THEN
   write(imainunit, '(A,A,A)')"       <PDataArray type=""Float32"" Name=""","Material_E","""/>"
  END IF
  END IF

 END SELECT
END DO
write(imainunit, '(A)')"    </PCellData>"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Here comes the mesh data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(imainunit, *)"    <PPoints>"
write(imainunit, *)"      <PDataArray type=""Float32"" Name=""Points"" NumberOfComponents=""3""/>"
write(imainunit, *)"    </PPoints>"

do iproc=1,subnodes
 filename=" "
 WRITE(filename(1:),'(A9,I3.3,A1,I5.5,A4)') 'res_node_',iproc,'.',iO,".vtu"
 write(imainunit, '(A,A,A)')"      <Piece Source=""",trim(adjustl(filename)),"""/>"  
end do
write(imainunit, *)"  </PUnstructuredGrid>"
write(imainunit, *)"  </VTKFile>"
close(imainunit)

END SUBROUTINE Output_VTK_main


SUBROUTINE Output_GMV_fields(iO,dcoor,kvert)
USE def_FEAT
USE  PP3D_MPI, ONLY:myid,showid,subnodes
USE var_QuadScalar,ONLY:QuadSc,LinSc,Viscosity,Distance,Distamce,mgNormShearStress,myALE
USE var_QuadScalar,ONLY:Tracer
USE var_QuadScalar,ONLY:myExport
USE var_QuadScalar,ONLY:myFBM,knvt,knet,knat,knel

IMPLICIT NONE
REAL*8 dcoor(3,*)
INTEGER kvert(8,*),iO
INTEGER NoOfVert,NoOfElem
CHARACTER cf*(30),cmm*(30),cm*(30)
INTEGER i,j,iField
INTEGER :: iOutUnit=555

 NoOfElem = KNEL(ILEV)
 NoOfVert = KNVT(ILEV)
 cf  = ' '
 cm  = ' '
 cmm = ' '

 WRITE(cf(1:),'(A9,A2,A1,I4.4,A4)') '_gmv/res_','**','_',iO,".gmv"
 WRITE(cmm(1:),'(A4,I2.2,A4)') 'msh_',myid,".gmv"
 WRITE(cm (1:),'(A9,I2.2,A4)') '_gmv/msh_',myid,".gmv"
!  WRITE(cf(1:),'(A9,A3,A1,I5.5,A4)') '_gmv/res_','***','_',iO,".gmv"
!  WRITE(cmm(1:),'(A4,I3.3,A4)') 'msh_',myid,".gmv"
!  WRITE(cm (1:),'(A9,I3.3,A4)') '_gmv/msh_',myid,".gmv"

 IF(myid.eq.showid) WRITE(*,'(104("="))') 
 IF(myid.eq.showid) WRITE(*,*) "Outputting gmv file into ","[",ADJUSTL(TRIM(cf)),"]"
 WRITE(cf(10:11),'(I2.2)') myid
!  WRITE(cf(10:12),'(I3.3)') myid

 IF (iO.EQ.0) THEN
  OPEN (UNIT=iOutUnit,FILE=ADJUSTL(TRIM(cm)))

  WRITE(iOutUnit,'(A)')'gmvinput ascii'
  WRITE(iOutUnit,*)'nodes ',NoOfVert

  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(1,i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(2,i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1200) REAL(dcoor(3,i))
  END DO

  WRITE(iOutUnit,*)'cells ',NoOfElem
  DO i=1,NoOfElem
   WRITE(iOutUnit,*)'hex 8'
   WRITE(iOutUnit,1300) (kvert(j,i),j=1,8)
  END DO

  WRITE(iOutUnit,*)  'endgmv'

  CLOSE  (iOutUnit)
 
 END IF
! CALL XGMVMS(iOutUnit,cm,NoOfElem,NoOfVert,kvert,dcoor)

 OPEN (UNIT=iOutUnit,FILE=cf)
 WRITE(iOutUnit,'(A)')'gmvinput ascii'
 WRITE(UNIT=iOutUnit,FMT=*) 'nodes fromfile "',ADJUSTL(TRIM(cmm)),'"'
 WRITE(UNIT=iOutUnit,FMT=*) 'cells fromfile "',ADJUSTL(TRIM(cmm)),'"'

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Velocity')
  WRITE(iOutUnit,*)  'velocity 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValU(i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValV(i))
  END DO
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(QuadSc%ValW(i))
  END DO

 END SELECT
END DO
 
WRITE(iOutUnit,*)  'variable'

DO iField=1,SIZE(myExport%Fields)

 SELECT CASE(ADJUSTL(TRIM(myExport%Fields(iField))))
 CASE('Pressure_V')
  WRITE(iOutUnit,*)  'pressure_V 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(LinSc%Q2(i))
  END DO

 CASE('Temperature')
  WRITE(iOutUnit,*)  'temperature 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Tracer%Val(NLMAX)%x(i))
  END DO

 CASE('Mixer')
  WRITE(iOutUnit,*)  'mixer 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Distamce(i))
  END DO

 CASE('Monitor')
  WRITE(iOutUnit,*)  'monitor 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(myALE%Monitor(i))
  END DO

 CASE('Viscosity')
  WRITE(iOutUnit,*)  'viscosity 1'
  DO i=1,NoOfVert
   WRITE(iOutUnit,1000) REAL(Viscosity(i))
  END DO

 CASE('Pressure_E')
  IF (ILEV.EQ.NLMAX-1) THEN
   WRITE(iOutUnit,*)  'pressure_E 0'
   DO i=1,NoOfElem
    j = 4*(i-1) + 1
    WRITE(iOutUnit,1000) REAL(LinSc%ValP(ILEV)%x(j))
   END DO
  END IF

 END SELECT
END DO

 WRITE(iOutUnit,*)  'endvars'
 WRITE(iOutUnit,*)  'probtime ',REAL(timens)
 WRITE(iOutUnit,*)  'endgmv'
 CLOSE(iOutUnit)

1000  FORMAT(E12.5)
1100  FORMAT(8I8)
1200  FORMAT(E12.5)
1300  FORMAT(8I8)

END SUBROUTINE Output_GMV_fields

SUBROUTINE OutputTriMesh(dcorvg,kvert,knpr,nvt,nel,iO)
!USE QuadScalar,ONLY:myQ2Coor
USE var_QuadScalar,ONLY:ilev
USE Parametrization, ONLY : myParBndr,nBnds
IMPLICIT NONE
REAL*8 dcorvg(3,*)
INTEGER nvt,nel,kvert(8,*),knpr(*),i,iO,j,iIt,iBnds
CHARACTER cf*(40)

WRITE(cf,'(A)') '_vtk/mesh.tri'
WRITE(*,*) "Outputting actual Coarse mesh into: '"//ADJUSTL(TRIM(cf))//"'"
OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
WRITE(1,*) 'Coarse mesh exported by DeViSoR TRI3D exporter'
WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
WRITE(1,'(2I8,A)') NEL,NVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

WRITE(1,'(A)') 'DCORVG'
DO i = 1 ,nvt
 WRITE(1,'(3ES13.5)') dcorvg(:,i)
END DO

WRITE(1,'(A)') 'KVERT'
DO i = 1 ,nel
 WRITE(1,'(8I8)') kvert(:,i)
END DO

WRITE(1,'(A)') 'KNPR'
DO i = 1 ,nvt
 WRITE(1,'(I8)') knpr(i)
END DO

CLOSE(1)

 DO iBnds = 1, nBnds
  cf = ' '
  WRITE(cf,'(A)') "_vtk/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
  WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
  OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
   j=0
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     j = j + 1
    END IF
   END DO
   WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
   WRITE(1,'(A)')    "'"//myParBndr(iBnds)%Parameters//"'"
   DO i=1,NVT
    IF (myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
     WRITE(1,'(I8,A)') i
    END IF
   END DO
   CLOSE(1)
 END DO
 
! DO iIt=1,Properties%nInterface
!  IF (allocated(myTSurf(iIT)%T)) then
!  WRITE(cf,'(A11,I3.3,A1,I3.3,A4)') '_vtk/cSURF_',iO,'_',iIT,'.csv'
!  write(*,*) "writing surface data to: '", ADJUSTL(trim(cf)),"'"
!  OPEN(FILE=ADJUSTL(TRIM(cf)),UNIT=171)
!  WRITE(171,'(A)') 'X, Y, Z'
!  DO i=1,myTSurf(iIT)%nT
!   DO j=1,9
!    WRITE(171,'(2(ES12.4,A1),ES12.4)') myTSurf(iIT)%T(i)%C(1,j),",",myTSurf(iIT)%T(i)%C(2,j),',',myTSurf(iIT)%T(i)%C(3,j)
!   END DO
!  END DO
!  CLOSE(171)
!  end if
! END DO

END SUBROUTINE OutputTriMesh
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Release_ListFiles_SSE(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
! call Release_ListFile('v',iO)
! call Release_ListFile('p',iO)
! call Release_ListFile('d',iO)
! call Release_ListFile('t',iO)
! call Release_ListFile('s',iO)
! 
!  CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
! ! CALL SetCoor(DWORK(L(KLCVG(NLMAX))))
!  call Release_ListFile('x',iO)
! 
!  CONTAINS
!  SUBROUTINE SetCoor(dc)
!  real*8 dc(3,*)
!  integer i
!  
!  do i=1,QuadSc%ndof
!   QuadSc%AuxU(i) = dc(1,i)
!   QuadSc%AuxV(i) = dc(2,i)
!   QuadSc%AuxW(i) = dc(3,i)
!  end do
!  
!  END SUBROUTINE SetCoor
!  
! END SUBROUTINE Release_ListFiles_SSE
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Release_ListFiles_SSE_Q1_Scalar(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
! call Release_ListFile('v',iO)
! call Release_ListFile('p',iO)
! call Release_ListFile('q',iO)
! ! call Release_ListFile('d',iO)
! call Release_ListFile('t',iO)
! ! call Release_ListFile('s',iO)
! 
!  CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
!  call Release_ListFile('x',iO)
! 
!  CONTAINS
!  SUBROUTINE SetCoor(dc)
!  real*8 dc(3,*)
!  integer i
!  
!  do i=1,QuadSc%ndof
!   QuadSc%AuxU(i) = dc(1,i)
!   QuadSc%AuxV(i) = dc(2,i)
!   QuadSc%AuxW(i) = dc(3,i)
!  end do
!  
!  END SUBROUTINE SetCoor
!  
! END SUBROUTINE Release_ListFiles_SSE_Q1_Scalar
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Release_ListFiles_SSE_temp(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
! call Release_ListFile('t',iO)
! 
! END SUBROUTINE Release_ListFiles_SSE_temp
!
!---------------------------------------------------------------------------
!
SUBROUTINE Load_ListFiles_PRT_Tracer(iO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,master
USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
implicit none
integer iO
integer nLengthV,nLengthE,LevDif
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

! -------------- workspace -------------------
INTEGER  NNWORK
PARAMETER (NNWORK=1)
INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)

INTEGER            :: KWORK(1)
REAL               :: VWORK(1)
DOUBLE PRECISION   :: DWORK(NNWORK)

COMMON       NWORK,IWORK,IWMAX,L,DWORK
EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! -------------- workspace -------------------

 call Load_ListFile('v',iO)
 call Load_ListFile('x',iO)
 call Load_ListFile('s',iO)
 
 if (myid.ne.master) CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
!  CALL SetCoor(DWORK(L(KLCVG(NLMAX))))
 
 IF (myid.EQ.0) THEN
   CALL CreateDumpStructures(0)
 ELSE
   LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
   CALL CreateDumpStructures(LevDif)
 END IF
 
 ILEV = LinSc%prm%MGprmIn%MedLev

 nLengthV = (2**(ILEV-1)+1)**3
 nLengthE = mg_mesh%level(NLMIN)%nel

 ALLOCATE(SendVect(3,nLengthV,nLengthE))

 CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                             mg_mesh%level(ILEV)%kvert,&
                             nLengthV,&
                             nLengthE,&
                             mg_mesh%level(ILEV)%nel,&
                             mg_mesh%level(ILEV)%nvt)
 DEALLOCATE(SendVect)
 
 CONTAINS
 SUBROUTINE SetCoor(dc)
 real*8 dc(3,*)
 integer i
 
 do i=1,QuadSc%ndof
  dc(1,i) = QuadSc%AuxU(i)
  dc(2,i) = QuadSc%AuxV(i)
  dc(3,i) = QuadSc%AuxW(i)
 end do
 
 END SUBROUTINE SetCoor
 
END SUBROUTINE Load_ListFiles_PRT_Tracer
!
!---------------------------------------------------------------------------
!
SUBROUTINE Load_ListFiles_General(iO,cLIST)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,master
USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
implicit none
integer,intent(in) :: iO
character(len=*),intent(in) :: cLIST
! local variables
character*(1),allocatable :: cFLD(:)
integer nLengthV,nLengthE,LevDif,iFld,iERR,nFld
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

IF (myid.eq.1) write(*,*) 'General List Reader call for:['//ADJUSTL(TRIM(cLIST))//']'

nFLD = (LEN(cLIST)+1)/2

allocate(cFLD(nFLD))
read(cLIST,*,IOSTAT=iERR) cFLD

if (ierr.ne.0) then
 write(*,*) 'Wrong field identifier defined for dump reading!'
 STOP 65
end if

DO iFld = 1,nFLD
 call Load_ListFile(cFLD(iFLD),iO)
END DO

DO iFld = 1,nFLD
 if (cFLD(iFLD).eq.'X'.or.cFLD(iFLD).eq.'x') then
  if (myid.ne.master) CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
  
  IF (myid.EQ.0) THEN
    CALL CreateDumpStructures(0)
  ELSE
    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
    CALL CreateDumpStructures(LevDif)
  END IF
  
  ILEV = LinSc%prm%MGprmIn%MedLev

  nLengthV = (2**(ILEV-1)+1)**3
  nLengthE = mg_mesh%level(NLMIN)%nel

  ALLOCATE(SendVect(3,nLengthV,nLengthE))

  CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
                              mg_mesh%level(ILEV)%kvert,&
                              nLengthV,&
                              nLengthE,&
                              mg_mesh%level(ILEV)%nel,&
                              mg_mesh%level(ILEV)%nvt)
  DEALLOCATE(SendVect)
 end if
END DO

 CONTAINS
 SUBROUTINE SetCoor(dc)
 real*8 dc(3,*)
 integer i
 
 do i=1,QuadSc%ndof
  dc(1,i) = QuadSc%AuxU(i)
  dc(2,i) = QuadSc%AuxV(i)
  dc(3,i) = QuadSc%AuxW(i)
 end do
 
 END SUBROUTINE SetCoor
 
END SUBROUTINE Load_ListFiles_General
!
!---------------------------------------------------------------------------
!
SUBROUTINE Release_ListFiles_General(iO,cLIST)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid,master
USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
implicit none
integer,intent(in) :: iO
character(len=*),intent(in) :: cLIST
! local variables
character*(1),allocatable :: cFLD(:)
integer nLengthV,nLengthE,LevDif,iFld,iERR,nFld
REAL*8 , ALLOCATABLE :: SendVect(:,:,:)

IF (myid.eq.1) write(*,*) 'General List Writer call for:['//ADJUSTL(TRIM(cLIST))//']'

nFLD = (LEN(cLIST)+1)/2

allocate(cFLD(nFLD))
read(cLIST,*,IOSTAT=iERR) cFLD

if (ierr.ne.0) then
 write(*,*) 'Wrong field identifier defined for dump reading!'
 STOP 65
end if

DO iFld = 1,nFLD
 if (cFLD(iFLD).eq.'X'.or.cFLD(iFLD).eq.'x') then
  CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
 END IF
 call Release_ListFile(cFLD(iFLD),iO)
END DO

 CONTAINS
 SUBROUTINE SetCoor(dc)
 real*8 dc(3,*)
 integer i
 
 do i=1,QuadSc%ndof
  QuadSc%AuxU(i) = dc(1,i)
  QuadSc%AuxV(i) = dc(2,i)
  QuadSc%AuxW(i) = dc(3,i)
 end do
 
 END SUBROUTINE SetCoor
 
END SUBROUTINE Release_ListFiles_General
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Load_ListFiles_Q1_Scalar(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid,master
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! integer nLengthV,nLengthE,LevDif
! REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
! 
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
!  call Load_ListFile('s',iO)
!  call Load_ListFile('p',iO)
!  call Load_ListFile('v',iO)
!  call Load_ListFile('x',iO)
!  call Load_ListFile('q',iO)
!  
!  if (myid.ne.master) CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
!  
!  IF (myid.EQ.0) THEN
!    CALL CreateDumpStructures(0)
!  ELSE
!    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
!    CALL CreateDumpStructures(LevDif)
!  END IF
!  
!  ILEV = LinSc%prm%MGprmIn%MedLev
! 
!  nLengthV = (2**(ILEV-1)+1)**3
!  nLengthE = mg_mesh%level(NLMIN)%nel
! 
!  ALLOCATE(SendVect(3,nLengthV,nLengthE))
! 
!  CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
!                              mg_mesh%level(ILEV)%kvert,&
!                              nLengthV,&
!                              nLengthE,&
!                              mg_mesh%level(ILEV)%nel,&
!                              mg_mesh%level(ILEV)%nvt)
!  DEALLOCATE(SendVect)
!  
!  CONTAINS
!  SUBROUTINE SetCoor(dc)
!  real*8 dc(3,*)
!  integer i
!  
!  do i=1,QuadSc%ndof
!   dc(1,i) = QuadSc%AuxU(i)
!   dc(2,i) = QuadSc%AuxV(i)
!   dc(3,i) = QuadSc%AuxW(i)
!  end do
!  
!  END SUBROUTINE SetCoor
!  
! END SUBROUTINE Load_ListFiles_Q1_Scalar
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Load_ListFiles_SSE(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid,master
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! integer nLengthV,nLengthE,LevDif
! REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
! 
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
! ! call Load_ListFile('t',iO)
!  call Load_ListFile('p',iO)
!  call Load_ListFile('v',iO)
!  call Load_ListFile('d',iO)
!  call Load_ListFile('x',iO)
!  call Load_ListFile('t',iO)
!  call Load_ListFile('q',iO)
! ! call Load_ListFile('s',iO)
!  
!  if (myid.ne.master) CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
! !  CALL SetCoor(DWORK(L(KLCVG(NLMAX))))
!  
!  IF (myid.EQ.0) THEN
!    CALL CreateDumpStructures(0)
!  ELSE
!    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
!    CALL CreateDumpStructures(LevDif)
!  END IF
!  
!  ILEV = LinSc%prm%MGprmIn%MedLev
! 
!  nLengthV = (2**(ILEV-1)+1)**3
!  nLengthE = mg_mesh%level(NLMIN)%nel
! 
!  ALLOCATE(SendVect(3,nLengthV,nLengthE))
! 
!  CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
!                              mg_mesh%level(ILEV)%kvert,&
!                              nLengthV,&
!                              nLengthE,&
!                              mg_mesh%level(ILEV)%nel,&
!                              mg_mesh%level(ILEV)%nvt)
!  DEALLOCATE(SendVect)
!  
!  CONTAINS
!  SUBROUTINE SetCoor(dc)
!  real*8 dc(3,*)
!  integer i
!  
!  do i=1,QuadSc%ndof
!   dc(1,i) = QuadSc%AuxU(i)
!   dc(2,i) = QuadSc%AuxV(i)
!   dc(3,i) = QuadSc%AuxW(i)
!  end do
!  
!  END SUBROUTINE SetCoor
!  
! END SUBROUTINE Load_ListFiles_SSE
!
!---------------------------------------------------------------------------
!
! SUBROUTINE Load_ListFiles_SSE_temp(iO)
! USE def_FEAT
! USE PP3D_MPI, ONLY:myid,showid,master
! USE var_QuadScalar,ONLY: LinSc,QuadSc,mg_mesh
! implicit none
! integer iO
! integer nLengthV,nLengthE,LevDif
! REAL*8 , ALLOCATABLE :: SendVect(:,:,:)
! 
! ! -------------- workspace -------------------
! INTEGER  NNWORK
! PARAMETER (NNWORK=1)
! INTEGER            :: NWORK,IWORK,IWMAX,L(NNARR)
! 
! INTEGER            :: KWORK(1)
! REAL               :: VWORK(1)
! DOUBLE PRECISION   :: DWORK(NNWORK)
! 
! COMMON       NWORK,IWORK,IWMAX,L,DWORK
! EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
! ! -------------- workspace -------------------
! 
!  call Load_ListFile('v',iO)
!  call Load_ListFile('d',iO)
!  call Load_ListFile('x',iO)
!  call Load_ListFile('t',iO)
!  call Load_ListFile('s',iO)
!  
!  if (myid.ne.master) CALL SetCoor(mg_mesh%level(NLMAX+1)%dcorvg)
!  
!  IF (myid.EQ.0) THEN
!    CALL CreateDumpStructures(0)
!  ELSE
!    LevDif = LinSc%prm%MGprmIn%MedLev - NLMAX
!    CALL CreateDumpStructures(LevDif)
!  END IF
!  
!  ILEV = LinSc%prm%MGprmIn%MedLev
! 
!  nLengthV = (2**(ILEV-1)+1)**3
!  nLengthE = mg_mesh%level(NLMIN)%nel
! 
!  ALLOCATE(SendVect(3,nLengthV,nLengthE))
! 
!  CALL SendNodeValuesToCoarse(SendVect,mg_mesh%level(NLMAX)%dcorvg,&
!                              mg_mesh%level(ILEV)%kvert,&
!                              nLengthV,&
!                              nLengthE,&
!                              mg_mesh%level(ILEV)%nel,&
!                              mg_mesh%level(ILEV)%nvt)
!  DEALLOCATE(SendVect)
!  
!  CONTAINS
!  SUBROUTINE SetCoor(dc)
!  real*8 dc(3,*)
!  integer i
!  
!  do i=1,QuadSc%ndof
!   dc(1,i) = QuadSc%AuxU(i)
!   dc(2,i) = QuadSc%AuxV(i)
!   dc(3,i) = QuadSc%AuxW(i)
!  end do
!  
!  END SUBROUTINE SetCoor
!  
! END SUBROUTINE Load_ListFiles_SSE_temp
!
!---------------------------------------------------------------------------
!
SUBROUTINE Load_ListFile(cF,iO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myDump,LinSc,QuadSc,Screw,temperature,mySegmentIndicator
USE var_QuadScalar,ONLY:GenLinScalar
USE Sigma_User, ONLY: myProcess
USE iniparser, ONLY : inip_openFileForReading
implicit none
integer :: iO
character :: cf*(*)
CHARACTER*(512) :: cFile,cFileBU
integer ilen,ndofs,ivt,ifile,iFld,ierr
logical :: bExists

if (myid.eq.0) return
 
 cfile = '_dump/processor_'
 ilen = len(trim(adjustl(cfile)))
 write(cfile(ilen+1:),'(I0)') myid
 
 ilen = len(trim(adjustl(cfile)))
 write(cfile(ilen+1:),'(A,I0)') "/",iO

 ilen = len(trim(adjustl(cfile)))
 
 if (adjustl(trim(cf)).eq.'q'.or.adjustl(trim(cf)).eq.'Q') THEN
  if (allocated(GenLinScalar%Fld)) then
   cFileBU = cFile
   DO iFld=1,GenLinScalar%nOfFields
    cFile = cFileBU
    write(cfile(ilen+1:),'(A)') "/"//adjustl(trim(GenLinScalar%prm%cField(iFld)))//".lst"
    call inip_openFileForReading(cfile, ifile, .true.)
    IF (ifile.gt.0) then
     ndofs = QuadSc%ndof
     read(ifile,*) 
     DO ivt=1,ndofs
      read(ifile,*,iostat=ierr) GenLinScalar%fld(iFld)%Val(ivt)
     END DO 
     write(*,*) 'list file loaded: "'//trim(adjustl(cfile))//'"'
     close(ifile)
    else
     if (adjustl(trim(GenLinScalar%prm%cField(iFld))).eq.'temp') THEN
      GenLinScalar%fld(iFld)%Val = myProcess%T0
     else
      GenLinScalar%fld(iFld)%Val = 0.0d0
     end if
     write(*,*) 'the LST file '//adjustl(trim(cfile))//' is NOT available! skipping the read sequence ... '
    end if
   END DO ! iFld
  end if
 END IF ! q,Q
 if (adjustl(trim(cf)).eq.'t'.or.adjustl(trim(cf)).eq.'T') THEN
   write(cfile(ilen+1:),'(A)') "/temperature.lst"

   call inip_openFileForReading(cfile, ifile, .true.)
   IF (ifile.gt.0) then
    ndofs = QuadSc%ndof
    read(ifile,*)
    DO ivt=1,ndofs
     read(ifile,*) temperature(ivt)
    END DO 
   else
    write(*,*) 'the LST file '//adjustl(trim(cfile))//' is NOT available! skipping the read sequence ... '
    
    temperature = myProcess%T0
   end if
   
  END IF
 if (adjustl(trim(cf)).eq.'p'.or.adjustl(trim(cf)).eq.'P') THEN
   write(cfile(ilen+1:),'(A)') "/pressure.lst"
   
   call inip_openFileForReading(cfile, ifile, .true.)
   
   ndofs = LinSc%ndof
   read(ifile,*)

   DO ivt=1,ndofs
    read(ifile,*) LinSc%ValP(NLMAX)%x(ivt)
   END DO 
  END IF
 if (adjustl(trim(cf)).eq.'v'.or.adjustl(trim(cf)).eq.'V') THEN
   write(cfile(ilen+1:),'(A)') "/velocity.lst"
   
   call inip_openFileForReading(cfile, ifile, .true.)
   
   ndofs = QuadSc%ndof
   read(ifile,*)

   DO ivt=1,ndofs
    read(ifile,*) QuadSc%ValU(ivt)
   END DO 
   DO ivt=1,ndofs
    read(ifile,*) QuadSc%ValV(ivt)
   END DO 
   DO ivt=1,ndofs
    read(ifile,*) QuadSc%ValW(ivt)
   END DO 
  END IF
 if (adjustl(trim(cf)).eq.'x'.or.adjustl(trim(cf)).eq.'X') THEN
   write(cfile(ilen+1:),'(A)') "/coordinates.lst"
   
   call inip_openFileForReading(cfile, ifile, .true.)
   
   ndofs = QuadSc%ndof
   read(ifile,*)

   DO ivt=1,ndofs
    read(ifile,*) QuadSc%AuxU(ivt)
   END DO 
   DO ivt=1,ndofs
    read(ifile,*) QuadSc%AuxV(ivt)
   END DO 
   DO ivt=1,ndofs
    read(ifile,*) QuadSc%AuxW(ivt)
   END DO 
  END IF
  if (adjustl(trim(cf)).eq.'d'.or.adjustl(trim(cf)).eq.'D') THEN
   write(cfile(ilen+1:),'(A)') "/distance.lst"
   
   call inip_openFileForReading(cfile, ifile, .true.)
   
   ndofs = QuadSc%ndof
   read(ifile,*)

   DO ivt=1,ndofs
    read(ifile,*) Screw(ivt)
   END DO 
  END IF
  
  if (adjustl(trim(cf)).eq.'s'.or.adjustl(trim(cf)).eq.'S') THEN
   IF (allocated(mySegmentIndicator)) then
    write(cfile(ilen+1:),'(A)') "/segment.lst"
    
    call inip_openFileForReading(cfile, ifile, .true.)
    
    ndofs = QuadSc%ndof
    read(ifile,*)

    DO ivt=1,ndofs
     read(ifile,*) mySegmentIndicator(2,ivt)
    END DO 
   END IF
  END IF

 if (adjustl(trim(cf)).eq.'q'.or.adjustl(trim(cf)).eq.'Q') THEN
 else
  write(*,*) 'list file loaded: "'//trim(adjustl(cfile))//'"'
  close(ifile)
 end if

END SUBROUTINE Load_ListFile
!
!---------------------------------------------------------------------------
!
SUBROUTINE Release_ListFile(cF,iO)
USE def_FEAT
USE PP3D_MPI, ONLY:myid,showid
USE var_QuadScalar,ONLY:myDump,LinSc,QuadSc,Screw,temperature,mySegmentIndicator
USE var_QuadScalar,ONLY:GenLinScalar
USE iniparser, ONLY : INIP_REPLACE,inip_openFileForWriting
implicit none

integer :: iO
character :: cf*(*)
CHARACTER*(512) :: cFile,cFileBU
integer ilen,ndofs,ivt,ifile,iFld,iRet,iNum
logical :: bExists

if (myid.eq.0) return
 
 cfile = '_dump/processor_'
 ilen = len(trim(adjustl(cfile)))
 write(cfile(ilen+1:),'(I0)') myid

!  inquire(=trim(adjustl(cfile)),Exist=bExists)
!  if (.not.bExists) call system('mkdir '//trim(adjustl(cfile)))
 
 ilen = len(trim(adjustl(cfile)))
 
 write(cfile(ilen+1:),'(A,I0)') "/",iO
!  inquire(file=trim(adjustl(cfile)),Exist=bExists)
 
!  if (.not.bExists) call system('mkdir '//trim(adjustl(cfile)))
 ilen = len(trim(adjustl(cfile)))
 if (adjustl(trim(cf)).eq.'q'.or.adjustl(trim(cf)).eq.'Q') THEN
  
   if (allocated(GenLinScalar%Fld)) then
    cFileBU = cFile
    DO iFld=1,GenLinScalar%nOfFields
    
     cFile = cFileBU
     write(cfile(ilen+1:),'(A)') "/"//adjustl(trim(GenLinScalar%prm%cField(iFld)))//".lst"
     call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
     
     ndofs = QuadSc%ndof
     write(ifile,'(A,I0)') "DofsTotal:",ndofs
     DO ivt=1,ndofs
      write(ifile,'(ES14.6)') GenLinScalar%fld(iFld)%Val(ivt)
     END DO 
     write(*,*) 'list file released: "'//trim(adjustl(cfile))//'"'
     call closeFile(ifile)
    END DO
   end if

 END IF
 if (adjustl(trim(cf)).eq.'t'.or.adjustl(trim(cf)).eq.'T') THEN
   write(cfile(ilen+1:),'(A)') "/temperature.lst"

   call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
   
   ndofs = QuadSc%ndof
   write(ifile,'(A,I0)') "DofsTotal:",ndofs

   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') temperature(ivt)
   END DO 
 END IF
 if (adjustl(trim(cf)).eq.'p'.or.adjustl(trim(cf)).eq.'P') THEN
   write(cfile(ilen+1:),'(A)') "/pressure.lst"

   call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
   ndofs = LinSc%ndof
   write(ifile,'(A,I0)') "DofsTotal:",ndofs

   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') LinSc%ValP(NLMAX)%x(ivt)
   END DO 
  END IF
 if (adjustl(trim(cf)).eq.'v'.or.adjustl(trim(cf)).eq.'V') THEN
   write(cfile(ilen+1:),'(A)') "/velocity.lst"

   call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
   ndofs = QuadSc%ndof
   write(ifile,'(A,I0)') "DofsTotal:",ndofs

   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%ValU(ivt)
   END DO 
   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%ValV(ivt)
   END DO 
   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%ValW(ivt)
   END DO 
  END IF
 if (adjustl(trim(cf)).eq.'x'.or.adjustl(trim(cf)).eq.'X') THEN
   write(cfile(ilen+1:),'(A)') "/coordinates.lst"

   call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
   ndofs = QuadSc%ndof
   write(ifile,'(A,I0)') "DofsTotal:",ndofs

   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%AuxU(ivt)
   END DO 
   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%AuxV(ivt)
   END DO 
   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') QuadSc%AuxW(ivt)
   END DO 
  END IF
  if (adjustl(trim(cf)).eq.'d'.or.adjustl(trim(cf)).eq.'D') THEN
   write(cfile(ilen+1:),'(A)') "/distance.lst"

   call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
   ndofs = QuadSc%ndof
   write(ifile,'(A,I0)') "DofsTotal:",ndofs

   DO ivt=1,ndofs
    write(ifile,'(ES14.6)') Screw(ivt)
   END DO 
  END IF
  if (adjustl(trim(cf)).eq.'s'.or.adjustl(trim(cf)).eq.'S') THEN
   IF (allocated(mySegmentIndicator)) then
    write(cfile(ilen+1:),'(A)') "/segment.lst"

    call inip_openFileForWriting(cfile, ifile, INIP_REPLACE, bExists, .true.)
    ndofs = QuadSc%ndof
    write(ifile,'(A,I0)') "DofsTotal:",ndofs

    DO ivt=1,ndofs
     write(ifile,'(ES14.6)') mySegmentIndicator(2,ivt)
    END DO 
   END IF
  END IF

  if (adjustl(trim(cf)).eq.'q'.or.adjustl(trim(cf)).eq.'Q') THEN
  else
   write(*,*) 'list file released: "'//trim(adjustl(cfile))//'"'
   call closeFile(ifile)
  end if

END SUBROUTINE Release_ListFile
!
!---------------------------------------------------------------------------
!

SUBROUTINE GetMeshVelocity()
  use var_QuadScalar, only:myALE, QuadSc, myQ2Coor 
  USE PP3D_MPI, ONLY:myid,showid,coarse,myMPI_Barrier,subnodes,&
      RECVI_myMPI,SENDI_myMPI,RECVD_myMPI,SENDD_myMPI,COMM_Maximum
  REAL*8 dMaxVelo,daux
  INTEGER i

  dMaxVelo = 0d0
  IF (myid.ne.0) then
    DO i=1,QuadSc%ndof
    myALE%MeshVelo(:,i) = (myQ2Coor(:,i) -  myALE%Q2Coor_old(:,i))/tstep
    daux = myALE%MeshVelo(1,i)**2d0+myALE%MeshVelo(2,i)**2d0+myALE%MeshVelo(3,i)**2d0
    IF (dMaxVelo.lt.daux) dMaxVelo = daux
    END DO
  END IF

  CALL COMM_Maximum(dMaxVelo)

  IF (myid.eq.1) THEN
    WRITE(*,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
  END IF

END SUBROUTINE GetMeshVelocity
!
! ----------------------------------------------
!
SUBROUTINE ResetTimer()
 use var_QuadScalar, only:myStat 
 myStat%iNonLin=0
 myStat%iLinUVW=0
 myStat%iLinP=0
 myStat%tMGUVW=0d0
 myStat%tMGP=0d0
 myStat%tDefUVW=0d0
 myStat%tDefP=0d0
 myStat%tCorrUVWP=0d0
 myStat%tGMVOut=0d0
 myStat%tDumpOut=0d0
 myStat%tSmat=0d0
 myStat%tKmat=0d0
 myStat%tDmat=0d0
 myStat%tMmat=0d0
 myStat%tCmat=0d0
 myStat%tRestUVW=0d0
 myStat%tProlUVW=0d0
 myStat%tSmthUVW=0d0
 myStat%tSolvUVW=0d0
 myStat%tRestP=0d0
 myStat%tProlP=0d0
 myStat%tSmthP=0d0
 myStat%tSolvP=0d0
 myStat%tCommP=0d0
 myStat%tCommV=0d0
 myStat%tCommS=0d0

END SUBROUTINE ResetTimer
!
! ----------------------------------------------
!
subroutine ReduceMesh_sse_mesh()
USE def_FEAT
USE var_QuadScalar,ONLY:LinSc,mg_Mesh,myBoundary
USE Parametrization, ONLY : myParBndr,nBnds
USE PP3D_MPI, ONLY:myid

implicit none
integer, allocatable :: ElementIndices(:),VertexIndices(:),FaceIndices(:)
integer, allocatable :: knpr(:)
integer i,j,k,ivt,iat,iadj,ibnds,nnel,nnvt,nnat
integer ivt1,ivt2,ivt3,ivt4
character cTRIFolder*(256),cf*(256)
INTEGER NeighA(4,6)
logical bBndry
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

if (myid.eq.0) then
 ilev=nlmin
 CALL SETLEV(2)
 
 nel = mg_Mesh%level(ILEV)%nel
 nvt = mg_Mesh%level(ILEV)%nvt
 nat = mg_Mesh%level(ILEV)%nat
 
 allocate(ElementIndices(nel),VertexIndices(nvt),FaceIndices(nat))
 allocate(knpr(nvt))
 
 knpr = 0
 ElementIndices = -1
 VertexIndices = -1
 FaceIndices = -1
 
 nnel = 0
 nnat = 0
 nnvt = 0
 do i=1,nel
  if (LinSc%knprP(ilev)%x(i).eq.0) THEN
   nnel = nnel + 1
   ElementIndices(i) = nnel
   do j=1,8
    ivt = mg_Mesh%level(ILEV)%kvert(j,i)
    if (VertexIndices(ivt).lt.0) then
     nnvt = nnvt + 1
     VertexIndices(ivt) = nnvt
    end if
   end do
  end if
 end do 
 
 DO i=1,nel
  do j=1,6
   if (ElementIndices(i).gt.0) then
    iadj = mg_Mesh%level(ILEV)%kadj(j,i)
    if (iadj.eq.0) then
     bBndry = .True.
    elseif (ElementIndices(iadj).lt.0) then
     bBndry = .True.
    else
     bBndry = .False.
    end if
    
    IF (bBndry) THEN
      FaceIndices(mg_Mesh%level(ILEV)%karea(j,i)) = 1
      knpr(mg_Mesh%level(ILEV)%kvert(neighA(:,j),i)) = 1
    end if
    
   end if
  end do
 end do
 
 cTRIFolder = 'ReducedMeshDir'
 CALL CheckIfFolderIsThereCreateIfNot(cTRIFolder,-1)
 
 write(*,'(2(A,I0,A,I0))') 'Reduction of NEL : ',nel,',',nnel, ' NVT : ',nvt,',',nnvt
 OPEN(FILE=ADJUSTL(TRIM(cTRIFolder))//'/ReducedMesh.tri',unit=1)
 
 WRITE(1,*) 'Coarse mesh exported by sse_mesh'
 WRITE(1,*) 'Parametrisierung PARXC, PARYC, TMAXC'
 WRITE(1,'(2I8,A)') nNEL,nNVT, " 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE"

 nnvt = 0
 VertexIndices = -1
 
 WRITE(1,*) 'DCORVG'
 do i=1,nel
  if (LinSc%knprP(ilev)%x(i).eq.0) THEN
   do j=1,8
    ivt = mg_Mesh%level(ILEV)%kvert(j,i)
    if (VertexIndices(ivt).lt.0) then
     nnvt = nnvt + 1
     VertexIndices(ivt) = nnvt
     write(1,'(3ES13.5)') mg_Mesh%level(ILEV)%dcorvg(:,ivt)
    end if
   end do
  end if
 end do 

 WRITE(1,*) 'KVERT'
 DO i=1,nel
  if (ElementIndices(i).gt.0) write(1,'(8I8)') VertexIndices(mg_Mesh%level(ILEV)%kvert(:,i))
 end do
 
 nnvt = 0
 VertexIndices = -1
 WRITE(1,*) 'KNPR'
 do i=1,nel
  if (LinSc%knprP(ilev)%x(i).eq.0) THEN
   do j=1,8
    ivt = mg_Mesh%level(ILEV)%kvert(j,i)
    if (VertexIndices(ivt).lt.0) then
     nnvt = nnvt + 1
     VertexIndices(ivt) = nnvt
     write(1,'(I0)') knpr(ivt)
    end if
   end do
  end if
 end do 

 CLOSE(1)

 OPEN(UNIT=2,FILE=ADJUSTL(TRIM(cTRIFolder))//'/file.prj')
 WRITE(cf,'(A)') 'ReducedMesh.tri' 
 WRITE(2,'(A)') ADJUSTL(TRIM(cf))
  
 DO iBnds = 1, nBnds
  cf = ' '
  WRITE(cf,'(A)') ADJUSTL(TRIM(cTRIFolder))//"/"//ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
  WRITE(2,'(A)') ADJUSTL(TRIM(myParBndr(iBnds)%Names))//".par"
  WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
  OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
  
  k=1
  DO i=1,nel
   do j=1,6
    iat = mg_Mesh%level(ILEV)%karea(j,i)
    IF (k.eq.iat) THEN
     if (FaceIndices(iat).eq.1) then
      ivt1 = mg_Mesh%level(ILEV)%kvert(NeighA(1,j),i)
      ivt2 = mg_Mesh%level(ILEV)%kvert(NeighA(2,j),i)
      ivt3 = mg_Mesh%level(ILEV)%kvert(NeighA(3,j),i)
      ivt4 = mg_Mesh%level(ILEV)%kvert(NeighA(4,j),i)
      if (myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt1).and.myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt2).and.&
          myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt3).and.myParBndr(iBnds)%Bndr(ILEV)%Vert(ivt4)) then
          FaceIndices(iat) = FaceIndices(iat) + 1
      end if
     end if
     k = k + 1
    end if 
   End do
  End do
  
  j=0
  DO i=1,mg_mesh%level(ilev)%nvt
   IF (VertexIndices(i).gt.0.and.myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
    j = j + 1
   END IF
  END DO
  WRITE(1,'(I8,A)') j," "//myParBndr(iBnds)%Types
  WRITE(1,'(A)')    "'"//ADJUSTL(TRIM(myParBndr(iBnds)%Parameters))//"'"
  DO i=1,mg_mesh%level(ilev)%nvt
   IF (VertexIndices(i).gt.0.and.myParBndr(iBnds)%Bndr(ILEV)%Vert(i)) THEN
   WRITE(1,'(I0,A,I0)') VertexIndices(i)!,", ",i
   END IF
  END DO
  CLOSE(1)
 END DO
 
 knpr=0
 k=1
 DO i=1,nel
  do j=1,6
   iat = mg_Mesh%level(ILEV)%karea(j,i)
   IF (k.eq.iat) THEN
    if (FaceIndices(iat).eq.1) then
     knpr(mg_Mesh%level(ILEV)%kvert(NeighA(1,j),i)) = 1
     knpr(mg_Mesh%level(ILEV)%kvert(NeighA(2,j),i)) = 1
     knpr(mg_Mesh%level(ILEV)%kvert(NeighA(3,j),i)) = 1
     knpr(mg_Mesh%level(ILEV)%kvert(NeighA(4,j),i)) = 1
!     write(*,*) ivt1,ivt2,ivt3,ivt4
    end if
    k = k + 1
   end if 
  End do
 End do
 
 cf = ' '
 WRITE(cf,'(A)') ADJUSTL(TRIM(cTRIFolder))//"/"//"NewWall.par"
 WRITE(2,'(A)') "NewWall.par"
 WRITE(*,*) "Outputting actual parametrization into: '"//ADJUSTL(TRIM(cf))//"'"
 OPEN(UNIT=1,FILE=ADJUSTL(TRIM(cf)))
 
 j=0
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (knpr(i).eq.1) THEN
   j = j + 1
  END IF
 END DO
 WRITE(1,'(I8,A)') j," Wall"
 WRITE(1,'(A)')    "' '"
 DO i=1,mg_mesh%level(ilev)%nvt
  IF (knpr(i).eq.1) THEN
   WRITE(1,'(I0,A,I0)') VertexIndices(i)!,", ",i
  END IF
 END DO
 CLOSE(1)

 CLOSE(2)

end if

end subroutine ReduceMesh_sse_mesh
!
! ----------------------------------------------
!

