Program STLvsTRI
USE mSTLvsTRI
CHARACTER cInputFile*(256),cVal*(256),cKey*(256)


!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  cProjectFolder="e3d_input"
  cShortProjectFile = "file.prj"
  cProjectGridFile="Mesh.tri"
  cOFFMeshFile="surface.off"
  cAreaIntensityFile="area.txt"
  cInputFile = ADJUSTL(TRIM(cProjectFolder))//'/'//'param.txt'
!  geometryStart  = -50.6059641065686 -79.20000010999999 -74.55000011
!  geometryLength  = 101.2119282131452 158.40000022 74.55000021000002
!  voxelSize = 20.24238564262904 19.8000000275 14.910000042000004
!  voxelStart = -50.6059641065686 -79.20000010999999 -74.55000011
!  voxelAmount = 5 8 5
!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   BoxMesh%Extent(:,1) =[-5.06,-7.92,-7.46]
!   BoxMesh%Extent(:,2) = BoxMesh%Extent(:,1) + [10.12,15.84,+7.46]
!   BoxMesh%Division(:) =[+10,+16,+10]

 cKey='geometryStart'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) BoxMesh%Extent(:,1)

 cKey='geometryLength'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) BoxMesh%Extent(:,2)
 BoxMesh%Extent(:,2) = BoxMesh%Extent(:,1) + BoxMesh%Extent(:,2)

 cKey='voxelAmount'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) BoxMesh%Division(:)

call readOFFMesh(adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cOFFMeshFile)))

CALL CheckForIntersection()

CALL Output_VTK()

CALL Output_TriMesh()

CALL Output_AreaIntenisty()

END Program STLvsTRI

