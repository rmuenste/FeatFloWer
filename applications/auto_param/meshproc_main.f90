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

write(*,*)'param'
CALL GetParameters()

write(*,*)'readmesh'
CALL ReadMEsh()

write(*,*)'octtree'
!CALL BuildKedge()
CALL BuildOctTree()

write(*,*)'buildkarea'
CALL BuildKarea()

write(*,*)'buildfaceneigh'
CALL BuildFaceNeighborhood()

write(*,*)'param'
CALL ReadParametrization()

DO i=1,nLoops
 CALL BuildFaceLists(nTol(i),aTol(i))
 write(*,*) i,'--> NumberOfSurfaces=',NumberOfSurfaces
END DO

CALL ExportParFiles()

CALL Output_SingleSurfToVTK()

CALL Output_SurfToVTK()

CALL Output_VTK()

DEALLOCATE(dcorvg)
DEALLOCATE(kvert)

 CONTAINS
SUBROUTINE GetParameters

OPEN(1,file='param_vtu_mesher.cfg')

READ(1,*) cProjectFolder
write(*,*)cProjectFolder
READ(1,*) cShortProjectFile
write(*,*)cShortProjectFile
READ(1,*) nLoops
write(*,*)nLoops
allocate(nTol(nLoops),aTol(nLoops))

DO i=1,nLoops
 READ(1,*) nTol(i),aTol(i)
 write(*,*) nTol(i),aTol(i)
END DO

#if defined __INTEL_COMPILER
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
#endif

CLOSE(1)

END SUBROUTINE GetParameters

END PROGRAM AutoParam
