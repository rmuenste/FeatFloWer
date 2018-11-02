PROGRAM AutoParam
USE MeshProcDef
!USE Parametrization, ONLY : InitParametrization
IMPLICIT NONE
INTEGER :: nLoops
INTEGER,ALLOCATABLE :: nTol(:)
REAL*8,ALLOCATABLE :: aTol(:)
integer Reason,iFile
CHARACTER sCommand*(200)
LOGICAL THERE

CALL GetParameters()

CALL ReadMEsh()

!CALL BuildKedge()
CALL BuildOctTree()

CALL BuildKarea()

CALL BuildFaceNeighborhood()

CALL ReadParametrization()

DO i=1,nLoops
 CALL BuildFaceLists(nTol(i),aTol(i))
 write(*,*) i,'--> NumberOfSurfaces=',NumberOfSurfaces
END DO

CALL ExportParFiles()

CALL Output_SingleSurfToVTK()

CALL Output_SurfToVTK()

CALL Output_VTK()

! IF (nStitchFiles.gt.0) THEN
!  sCommand = 'rm -fr '//ADJUSTL(TRIM(cProjectFolder))//'/STITCH/*.par'
!  CALL system(ADJUSTL(TRIM(sCommand)))
! 
!  iFile = 0
!  DO 
!   iFile = iFile + 1
!   write(sFile,'(A,I1,A)')
!   INQUIRE(FILE=ADJUSTL(TRIM(cProjectFolder))//"/PAR/stitch",iFile,".par", EXIST=THERE) 
!   
!   sCommand = 'mv '//ADJUSTL(TRIM(cProjectFolder))//'01.par '//ADJUSTL(TRIM(cProjectFolder))//'/PAR/stitch.par'
!   CALL system(ADJUSTL(TRIM(sCommand)))
!  END DO
! END IF

DEALLOCATE(dcorvg)
DEALLOCATE(kvert)

 CONTAINS
SUBROUTINE GetParameters

OPEN(1,file='param.cfg')

READ(1,*) cProjectFolder
READ(1,*) nLoops
allocate(nTol(nLoops),aTol(nLoops))

DO i=1,nLoops
 READ(1,*) nTol(i),aTol(i)
END DO

allocate(cParFile(1000))
DO 
 READ(1,*,IOSTAT=Reason) cParFile(nParFiles+1)
 IF (TRIM(ADJUSTL(cParFile(nParFiles+1))).eq.'STITCH') exit
 if (Reason.lt.0) EXIT
 nParFiles = nParFiles + 1
END DO

allocate(cStitchFile(1000))
DO 
 READ(1,*,IOSTAT=Reason) cStitchFile(nStitchFiles+1)
 if (Reason.lt.0) EXIT
 WRITE(*,*) ADJUSTL(TRIM(cStitchFile(nStitchFiles+1)))
 nStitchFiles = nStitchFiles + 1
END DO
CLOSE(1)

END SUBROUTINE GetParameters

END PROGRAM AutoParam
