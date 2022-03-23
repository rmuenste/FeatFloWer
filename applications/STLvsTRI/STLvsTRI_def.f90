MODULE mSTLvsTRI
USE types, ONLY: tMultiMesh,tMesh
USE var_QuadScalar
USE PP3D_MPI, ONLY:myid,master
USE MESH_Structures
USE OcttreeSearch

implicit none

type tOFFmesh
 REAL*8,allocatable :: coor(:,:)
 integer, allocatable :: kvert(:,:)
 integer nel,nvt
end type tOFFmesh
type (tOFFmesh) :: OFFmesh 
CHARACTER*(200) :: cOutputFolder,cShortProjectFile,cOFFMeshFile
integer iOK,iTriang,iRecLevel

 CONTAINS
 
SUBROUTINE readOFFMesh(cF)
CHARACTER*(*) cF
integer i,iaux

OPEN(File=trim(cF),unit=1)
read(1,*)
read(1,*)
read(1,*) OFFmesh%nvt,OFFmesh%nel

allocate(OFFmesh%coor(OFFmesh%nvt,3))
allocate(OFFmesh%kvert(OFFmesh%nel,3))

DO i=1,OFFmesh%nvt
 read(1,*) OFFmesh%coor(i,:)
END DO

DO i=1,OFFmesh%nel
 read(1,*) iaux,OFFmesh%kvert(i,:)
END DO
OFFmesh%kvert =  OFFmesh%kvert + 1

close(1)

END SUBROUTINE readOFFMesh
!
! -------------------------------------------------------------------------------
!
SUBROUTINE CheckForIntersection
integer i

ILEV=1

iOK = 0
DO iTriang=1,OFFmesh%nel

 iRecLevel = 0
 CALL AssignTriangle(0.1d0*OFFmesh%coor(OFFmesh%kvert(i,1),:),&
                     0.1d0*OFFmesh%coor(OFFmesh%kvert(i,2),:),&
                     0.1d0*OFFmesh%coor(OFFmesh%kvert(i,3),:))
END DO

write(*,*) 'nOK',iOK
 
END SUBROUTINE CheckForIntersection
!
! -------------------------------------------------------------------------------
!
RECURSIVE SUBROUTINE AssignTriangle(P1,P2,P3)
REAL*8 P1(3),P2(3),P3(3)
REAL*8 P12(3),P23(3),P31(3)
REAL*8 dist
integer ivt1,ivt2,ivt3

CALL FindInOctTreeFORTRI(mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt,P1,ivt1,DIST)
CALL FindInOctTreeFORTRI(mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt,P2,ivt2,DIST)
CALL FindInOctTreeFORTRI(mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt,P3,ivt3,DIST)

if (iRecLevel.gt.10) then
 write(*,*) iTriang,ivt1,ivt2,ivt3
 return
end if

if (ivt1.lt.0.or.ivt2.lt.0.or.ivt3.lt.0) then
 if (ivt1.lt.0.and.ivt2.lt.0.and.ivt3.lt.0) then
  return
 else
  return
 end if
else
  IF (ivt1.eq.ivt2.and.ivt2.eq.ivt3.and.ivt1.eq.ivt3) then
   iOK = iOK + 1
   return
  ELSE
  
   P12 = 0.5d0*(P1+P2)
   P23 = 0.5d0*(P2+P3)
   P31 = 0.5d0*(P1+P3)
   
   iRecLevel = iRecLevel + 1

   CALL AssignTriangle(P1,P12,P31)
   CALL AssignTriangle(P2,P12,P23)
   CALL AssignTriangle(P3,P23,P31)
   CALL AssignTriangle(P12,P23,P31)
  END IF
end if

END SUBROUTINE AssignTriangle

END 