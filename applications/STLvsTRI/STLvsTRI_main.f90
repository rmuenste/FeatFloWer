Program STLvsTRI
USE mSTLvsTRI
use f90getopt
use Sigma_User, only : mySigma,myProcess

CHARACTER cInputFile*(256),cVal*(256),cKey*(256),cE3Dfile*(256)
REAL*8 dSurfInt
integer nOfCritical
type(option_s)              :: opts(4)

opts(1) = option_s('inputfolder', .true.,  'f')
opts(2) = option_s('NumberOfElementsInit', .true.,  'n')
opts(3) = option_s('dSurfIntCrit', .true.,  'c')
opts(4) = option_s('help',  .false., 'h')

 cProjectFolder=" "
! check the command line arguments
do
    select case (getopt('f:n:c:h', opts))
        case (char(0))
            exit
        case ('f')
            read(optarg,*) cProjectFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cProjectFolder))
        case ('-f')
            read(optarg,*) cProjectFolder
            write(*,*)'Input Folder: ', ADJUSTL(TRIM(cProjectFolder))
        case ('n')
            read(optarg,*) NumberOfElementsInit
            write(*,*)'NumberOfElementsInit: ', NumberOfElementsInit
        case ('-n')
            read(optarg,*) NumberOfElementsInit
            write(*,*)'NumberOfElementsInit: ', NumberOfElementsInit
        case ('c')
            read(optarg,*) dSurfIntCrit
            write(*,*)'dSurfIntCrit: ', dSurfIntCrit
        case ('-c')
            read(optarg,*) dSurfIntCrit
            write(*,*)'dSurfIntCrit: ', dSurfIntCrit
        case ('h')
          call print_help()
          call exit(0)
        case default
    end select
end do

! INTEGER :: NumberOfElementsInit=1500
! REAL*8 :: dSurfIntCrit = 2.5d0

IF (cProjectFolder.eq." ") THEN
 WRITE(*,*) 'Input folder is an empty string!'
 call print_help()
 call exit(0)
END IF

write(cE3Dfile,'(A)') ADJUSTL(TRIM(cProjectFolder))//'/'//'setup.e3d'
CALL ReadS3Dfile(cE3Dfile)

!!Scale back to mm!!!
BoxMesh%Extent(:,1) = 10d0*mySigma%DIE_Start  
BoxMesh%dsize(:)    = 10d0*mySigma%DIE_Length
BoxMesh%Extent(:,2) = BoxMesh%Extent(:,1) + BoxMesh%dsize(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 cShortProjectFile = "file.prj"
 cProjectGridFile="Mesh.tri"
 cOFFMeshFile="surface.off"
 cAreaIntensityFile="area.txt"
 
 DomainVolume = BoxMesh%dsize(1)*BoxMesh%dsize(2)*BoxMesh%dsize(3)

 call readOFFMesh(adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cOFFMeshFile)))
 
 NumberOfElements = NumberOfElementsInit

 VoxelVolume = DomainVolume/DBLE(NumberOfElements)
 VoxelSize = 8d0*myProcess%ExtrusionGapSize*myProcess%ExtrusionGapFactor !VoxelVolume**(1d0/3d0)
 
 BoxMesh%Division(1) = 2*NINT(BoxMesh%dsize(1)/VoxelSize)+1
 BoxMesh%Division(2) = 2*NINT(BoxMesh%dsize(2)/VoxelSize)+1
 BoxMesh%Division(3) = 2*NINT(BoxMesh%dsize(3)/VoxelSize)+1

 NumberOfElements = (BoxMesh%Division(1)-1)*(BoxMesh%Division(2)-1)*(BoxMesh%Division(3)-1)
 VoxelVolume = DomainVolume/DBLE(NumberOfElements)
 VoxelSize = VoxelVolume**(1d0/3d0)
 UnityArea   = VoxelVolume**(2d0/3d0)
 
 WRITE(*,*) "Resolution, x,y,z: ",BoxMesh%Division
 WRITE(*,*) "Number of Elements: ",NumberOfElements

 CALL CheckForIntersection()

 ! Construct Surface Intensity
 TriMesh%I = TriMesh%d/UnityArea
 
 CALL CoarseSurfIntensity()
 
 CALL QualityCheck(dSurfInt,nOfCritical)
 WRITE(*,'(A,ES12.4,I0)') "MaxSurfaceIntensity & nOf Critical elements: ",dSurfInt,nOfCritical
 
 CALL Output_VTK()
 CALL Output_CoarseMeshVTK()

 CALL Output_TriMesh()

 CALL Output_AreaIntenisty()

 CONTAINS

 subroutine print_help()
     print '(a, /)', 'command-line options:'
     print '(a)',    '  -f      Input Folder'
     print '(a)',    '  -n      Number of Initial Hex Elements'
     print '(a)',    '  -c      Value of Surface Intensity Criterion '
     print '(a)',    '  -h      Print usage information and exit'
 end subroutine print_help  
 
END Program STLvsTRI

