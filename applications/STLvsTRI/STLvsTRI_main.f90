Program STLvsTRI
USE mSTLvsTRI


!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 cProjectFolder="IN"
 cProjectGridFile="mesh.tri"
 cOFFMeshFile="surface.off"
 mg_Mesh%nlmax = 1
 mg_Mesh%nlmin = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!! INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NLMAX = mg_Mesh%nlmax 
NLMIN = mg_Mesh%nlmin 
mg_Mesh%maxlevel = mg_Mesh%nlmax+1
allocate(mg_mesh%level(mg_Mesh%maxlevel))
myid = 1
master = 0

call readTriCoarse(adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cProjectGridFile)), mg_mesh)

call refineMesh(mg_mesh, mg_mesh%maxlevel)

DO ILEV=NLMIN,NLMAX
if(myid.eq.1)then

write(*,'(8A10)')"MESH:","NVT","NAT","NET","NEL"
write(*,'(A10,8I10)')"L1:",mg_mesh%level(ILEV)%nvt,&
                           mg_mesh%level(ILEV)%nat,&
                           mg_mesh%level(ILEV)%net,&
                           mg_mesh%level(ILEV)%nel

end if
end do

call readOFFMesh(adjustl(trim(cProjectFolder))//'/'//adjustl(trim(cOFFMeshFile)))

ILEV=1
bWrite=.true.
CALL InitOctTree(mg_mesh%level(ILEV)%dcorvg,mg_mesh%level(ILEV)%nvt)

CALL CheckForIntersection()

END Program STLvsTRI

