Module Mesh_Structures

  DIMENSION KIV(2,12)
  DATA KIV /1,2, 2,3, 3,4, 4,1, 1,5, 2,6, &
    3,7, 4,8, 5,6, 6,7, 7,8, 8,5/

  DIMENSION KIAD(4,6)
  DATA KIAD/1,2,3,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/

  type t_connector3D
    integer, dimension(6) :: I_conData
  end type t_connector3D

contains

!================================================================================================
!                                   Sub: writeTriFile 
!================================================================================================
subroutine writeTriFile(mesh)
  USE var_QuadScalar
  type(tMesh) :: mesh
  integer :: cunit = 123
  integer :: ivert,ielem,istat 

  open (unit=cunit,file='mesh2.tri',action="write",iostat=istat)
  if(istat .ne. 0)then
    write(*,*)"Could not open file for writing. "
    stop          
  end if    

  write(cunit,'(A)')'Coarse mesh 3D'
  write(cunit,'(A)')'modified by tr2to3'
  write(cunit,'(6I6,A)')mesh%nel,mesh%nvt,mesh%nbct,mesh%nve,&
    mesh%nee,mesh%nae,' NEL,NVT,NBCT,NVE,NEE,NAE'
  write(cunit,'(A)')' DCORVG'

  do ivert=1,mesh%nvt
  write(cunit, '(A,3E16.7)')"  ",REAL(mesh%dcorvg(1,ivert)),&
    REAL(mesh%dcorvg(2,ivert)),&
    REAL(mesh%dcorvg(3,ivert))
  end do

  write(cunit,'(A)')' KVERT'

  do ielem=1,mesh%nel
  write(cunit, '(A,8I6)')" ",(mesh%kvert(1,ielem)),&
    (mesh%kvert(2,ielem)),&
    (mesh%kvert(3,ielem)),&
    (mesh%kvert(4,ielem)),&
    (mesh%kvert(5,ielem)),&
    (mesh%kvert(6,ielem)),&
    (mesh%kvert(7,ielem)),&
    (mesh%kvert(8,ielem))
  end do

  write(cunit,'(A)')' KNPR'

  do ivert=1,mesh%nvt
  write(cunit, '(A,I1)')" ",mesh%knpr(ivert)
  end do

  close(cunit)

end subroutine

!================================================================================================
!                                 Sub: writeTriArrays  
!================================================================================================
subroutine writeTriArrays(mydcorvg, mykvert, mykedge, mykadj,&
    mykarea, mnvt, mnel, mnet, cfile)
  USE var_QuadScalar

  REAL*8  mydcorvg(3,*)

  integer mykvert(8,*)

  integer mykedge(12,*)

  integer mykadj(6,*)

  integer mykarea(6,*)

  CHARACTER (len = 60) :: cfile 

  integer :: mnvt
  integer :: mnel 
  integer :: mnet 

  integer :: cunit = 123
  integer :: ivert,ielem,iedge 

  open (unit=cunit,file=cfile)

  write(cunit,'(A)')'Coarse mesh 3D'
  write(cunit,'(A)')'modified by tr2to3'
  write(cunit,'(2I6,A)')mnel,mnvt,&
    ' 0 8 12 6 NEL,NVT,NBCT,NVE,NEE,NAE'
  write(cunit,'(A)')' DCORVG'

  do ivert=1,mnvt
  write(cunit, '(A,3E16.7)')"  ",REAL(mydcorvg(1,ivert)),&
    REAL(mydcorvg(2,ivert)),&
    REAL(mydcorvg(3,ivert))
  end do

  write(cunit,'(A)')' KVERT'

  do ielem=1,mnel
  write(cunit, '(A,8I6)')" ",(mykvert(1,ielem)),&
    (mykvert(2,ielem)),&
    (mykvert(3,ielem)),&
    (mykvert(4,ielem)),&
    (mykvert(5,ielem)),&
    (mykvert(6,ielem)),&
    (mykvert(7,ielem)),&
    (mykvert(8,ielem))
  end do

  write(cunit,'(A)')' KNPR'

  do ivert=1,mnvt
  write(cunit, '(A)')" 0"
  end do

  close(cunit)

end subroutine

!================================================================================================
!                                 Sub: writeArrays  
!================================================================================================
subroutine writeArrays(mydcorvg, mykvert, mykedge, mykadj,&
    mykarea, mnvt, mnel, mnet,cfile)
  USE var_QuadScalar

  REAL*8  mydcorvg(3,*)

  integer mykvert(8,*)

  integer mykedge(12,*)

  integer mykadj(6,*)

  integer mykarea(6,*)

  CHARACTER (len = 60) :: cfile 

  integer :: mnvt
  integer :: mnel 
  integer :: mnet 

  integer :: cunit = 123
  integer :: ivert,ielem,iedge 

  open (unit=cunit,file=cfile)

  write(cunit,'(3I6,A)')mnel,mnvt,mnet,' NEL,NVT,NET'

  do ivert=1,mnvt
  write(cunit, '(A,3E16.7)')"  ",REAL(mydcorvg(1,ivert)),&
    REAL(mydcorvg(2,ivert)),&
    REAL(mydcorvg(3,ivert))
  end do

  write(cunit,'(A)')'KVERT'
  do ielem=1,mnel
  write(cunit, '(A,8I6)')" ",(mykvert(1,ielem)),&
    (mykvert(2,ielem)),&
    (mykvert(3,ielem)),&
    (mykvert(4,ielem)),&
    (mykvert(5,ielem)),&
    (mykvert(6,ielem)),&
    (mykvert(7,ielem)),&
    (mykvert(8,ielem))
  end do


  write(cunit,'(A)')'KADJ'
  do ielem=1,mnel
  write(cunit, '(A,6I6)')" ",(mykadj(1,ielem)),&
    (mykadj(2,ielem)),&
    (mykadj(3,ielem)),&
    (mykadj(4,ielem)),&
    (mykadj(5,ielem)),&
    (mykadj(6,ielem))
  end do


  write(cunit,'(A)')'KEDGE'
  do ielem=1,mnel
  write(cunit, '(A,12I6)')" ",(mykedge(1,ielem)),&
    (mykedge(2,ielem)),&
    (mykedge(3,ielem)),&
    (mykedge(4,ielem)),&
    (mykedge(5,ielem)),&
    (mykedge(6,ielem)),&
    (mykedge(7,ielem)),&
    (mykedge(8,ielem)),&
    (mykedge(9,ielem)),&
    (mykedge(10,ielem)),&
    (mykedge(11,ielem)),&
    (mykedge(12,ielem))
  end do

  close(cunit)

end subroutine

!================================================================================================
!                                 Sub: readTriCoarse  
!================================================================================================
subroutine readTriCoarse(CFILE, mgMesh)
  USE var_QuadScalar
  USE PP3D_MPI
  CHARACTER (len = *) :: CFILE
  type(tMultiMesh) :: mgMesh
  logical :: bExist
  integer :: cunit = 18
  integer :: i,istat

  inquire(file=adjustl(trim(cfile)),exist=bExist)
  if (.not. bExist) then
    write(*,*) "File '"//adjustl(trim(cfile)),"' could not be read ..."
    return
  end if

  open(cunit,file=adjustl(trim(cfile)),action="read",iostat=istat)
  if(istat .ne. 0)then
    write(*,*)'Could not open file: ',adjustl(trim(cfile))
  end if  

  read(cunit, *)
  read(cunit, *)
  read(cunit,*)mgMesh%level(1)%nel, mgMesh%level(1)%nvt, &
    mgMesh%level(1)%nbct, mgMesh%level(1)%nve, &
    mgMesh%level(1)%nee , mgMesh%level(1)%nae

  read(cunit, *)

  if(.not.associated(mgMesh%level(1)%dcorvg))then
    allocate(mgMesh%level(1)%dcorvg(3,mgMesh%level(1)%nvt))
  end if

  if(.not.allocated(mgMesh%level(1)%kvert))then
    allocate(mgMesh%level(1)%kvert(8,mgMesh%level(1)%nel))
  end if

  if(.not.allocated(mgMesh%level(1)%knpr))then
    allocate(mgMesh%level(1)%knpr(mgMesh%level(1)%nvt))
  end if

  do i=1,mgMesh%level(1)%nvt
  read(cunit,*)mgMesh%level(1)%dcorvg(1:3,i)
  !    if(myid.eq.0)then
  !      write(*,*),mgMesh%level(1)%dcorvg(1:3,i)
  !    end if
  end do
  read(cunit, *)

  do i=1,mgMesh%level(1)%nel
  read(cunit,*)mgMesh%level(1)%kvert(1:8,i)
  !    if(myid.eq.0)then
  !      write(*,*),mgMesh%level(1)%kvert(1:8,i)
  !    end if
  end do

  read(cunit, *)
  do i=1,mgMesh%level(1)%nvt
  read(cunit,*),mgMesh%level(1)%knpr(i)
  !    if(myid.eq.0)then
  !      write(*,*),mgMesh%level(1)%knpr(i)
  !    end if
  end do

  close(cunit)


end subroutine

!================================================================================================
!                                 Sub: refineMesh  
!================================================================================================
subroutine refineMesh(mgMesh,maxlevel, extended)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  implicit none

  type(tMultiMesh), intent(inout) :: mgMesh
  integer, intent(in) :: maxlevel
  logical, optional, intent(in) :: extended
  integer :: nfine
  integer :: ilevel
  integer :: icurr

  logical :: calculateExtendedConnectivity 

  if(present(extended))then
    calculateExtendedConnectivity = extended
  else
    calculateExtendedConnectivity = .false.
  end if

  nfine = maxlevel-1

  icurr = 1

  call genMeshStructures(mgMesh%level(icurr), calculateExtendedConnectivity)

  if(myid.eq.1)then
    write(*,*)'Refining to level: ',maxlevel
    write(*,*)'Number of refinement steps needed: ',nfine
  end if

  do ilevel=1,nfine

    icurr = icurr + 1 
    call refineMeshLevel(mgMesh%level(icurr-1), mgMesh%level(icurr))
    call genMeshStructures(mgMesh%level(icurr), calculateExtendedConnectivity)

    if(myid.eq.1)then
      write(*,*)'-----RefinementLevelFinished-----'
    end if

  end do

  do ilevel=1,maxlevel-1 
    deallocate(mgMesh%level(ilevel)%dcorvg)
    mgMesh%level(ilevel)%dcorvg => mgMesh%level(maxlevel)%dcorvg
  end do

end subroutine

!================================================================================================
!                                 Sub: refineMeshLevel  
!================================================================================================
subroutine refineMeshLevel(mesh0, mesh1)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  implicit none

  type(tMesh) :: mesh0
  type(tMesh) :: mesh1

  integer :: nnvt
  integer :: ielOffset
  integer :: iivt,ivt1,ivt2
  integer :: iedge,ielem,j,iface
  integer :: iistart,iend
  integer :: iel1
  integer :: iel2
  integer :: iel3
  integer :: iel4
  integer :: iel5
  integer :: iel6
  integer :: iel7
  integer :: midpointA
  integer :: midpointB
  integer :: midpointC
  integer :: midpointD
  integer :: midpointE
  integer :: midpointF
  integer :: ielmid

  integer, allocatable, dimension(:) :: midpointsAtFace

  mesh1%nvt = mesh0%nvt + mesh0%net + mesh0%nat + mesh0%nel

  if(.not.associated(mesh1%dcorvg))then
    allocate(mesh1%dcorvg(3,mesh1%nvt))
  end if

  do iivt=1,mesh0%nvt 
  mesh1%dcorvg(1:3,iivt) = mesh0%dcorvg(1:3,iivt)
  end do

  iedge = 1
  do iivt=mesh0%nvt+1,mesh0%nvt+mesh0%net
  ivt1 = mesh0%kved(1,iedge)
  ivt2 = mesh0%kved(2,iedge)

  mesh1%dcorvg(1:3,iivt) = 0.5d0 * (mesh0%dcorvg(1:3,ivt1) + &
    mesh0%dcorvg(1:3,ivt2))

!        write(*,'(A,I3,3E16.7,3I3)')'vertex edge : ',iivt,&
!                                   mesh1%dcorvg(1:3,iivt),&
!                                   ivt1,&
!                                   ivt2,iedge
  iedge = iedge + 1
  end do

  mesh1%net = iedge-1

  allocate(midpointsAtFace(mesh0%nat))

  iistart = mesh0%nvt+mesh0%net+1 
  iend    = mesh0%nvt+mesh0%net+mesh0%nat 
  iface   = 1
  do iivt=iistart,iend
  mesh1%dcorvg(1:3,iivt) = 0d0
  do j=1,4
  mesh1%dcorvg(1:3,iivt) = mesh1%dcorvg(1:3,iivt) + &
    mesh0%dcorvg(1:3,mesh0%kvar(j,iface))                  
  end do

  mesh1%dcorvg(1:3,iivt) = 0.25d0 * mesh1%dcorvg(1:3,iivt)

!        write(*,*)'vertex face : ',iivt,mesh1%dcorvg(1:3,iivt)
  midpointsAtFace(iface) = iivt
  iface = iface + 1
  end do

  mesh1%nat = iface-1

  iistart = mesh0%nvt+mesh0%net+mesh0%nat+1 
  iend   = mesh0%nvt+mesh0%net+mesh0%nat+mesh0%nel 

  ielem = 1
  do iivt=iistart,iend
  mesh1%dcorvg(1:3,iivt) = 0d0

  do j=1,8
  mesh1%dcorvg(1:3,iivt) = mesh1%dcorvg(1:3,iivt) + &
    mesh0%dcorvg(1:3,mesh0%kvert(j,ielem))       
  end do

  mesh1%dcorvg(1:3,iivt) = 0.125d0 * mesh1%dcorvg(1:3,iivt)

!        write(*,*)'vertex elem: ',iivt,mesh1%dcorvg(1:3,iivt)

  ielem = ielem + 1
  end do

  mesh1%nel = mesh0%nel * 8

  if(.not.allocated(mesh1%kvert))then
    allocate(mesh1%kvert(8,mesh1%nel))
  end if

  ielem = 1
  ielOffset = mesh0%nel
  do ielem=1,mesh0%nel

  iel1 = ielOffset + 1
  iel2 = ielOffset + 2
  iel3 = ielOffset + 3
  iel4 = ielOffset + 4
  iel5 = ielOffset + 5
  iel6 = ielOffset + 6
  iel7 = ielOffset + 7

  midpointA = midpointsAtFace(mesh0%KAREA(1,ielem))
  midpointB = midpointsAtFace(mesh0%KAREA(2,ielem))
  midpointC = midpointsAtFace(mesh0%KAREA(3,ielem))
  midpointD = midpointsAtFace(mesh0%KAREA(4,ielem))
  midpointE = midpointsAtFace(mesh0%KAREA(5,ielem))
  midpointF = midpointsAtFace(mesh0%KAREA(6,ielem))

  ielmid = mesh0%nvt + mesh0%net + mesh0%nat + ielem

  ielOffset = ielOffset + 7 

  ! 1st new element

  mesh1%kvert(1,ielem) = mesh0%kvert(1,ielem)
  mesh1%kvert(2,ielem) = mesh0%kedge(1,ielem) + mesh0%nvt
  mesh1%kvert(3,ielem) = midpointA
  mesh1%kvert(4,ielem) = mesh0%kedge(4,ielem) + mesh0%nvt
  mesh1%kvert(5,ielem) = mesh0%kedge(5,ielem) + mesh0%nvt
  mesh1%kvert(6,ielem) = midpointB
  mesh1%kvert(7,ielem) = ielmid
  mesh1%kvert(8,ielem) = midpointE

  ! 2nd new element

  mesh1%kvert(1,iel1) = mesh0%kvert(2,ielem)
  mesh1%kvert(2,iel1) = mesh0%kedge(2,ielem) + mesh0%nvt
  mesh1%kvert(3,iel1) = midpointA
  mesh1%kvert(4,iel1) = mesh0%kedge(1,ielem) + mesh0%nvt
  mesh1%kvert(5,iel1) = mesh0%kedge(6,ielem) + mesh0%nvt
  mesh1%kvert(6,iel1) = midpointC
  mesh1%kvert(7,iel1) = ielmid
  mesh1%kvert(8,iel1) = midpointB

  ! 3rd new element

  mesh1%kvert(1,iel2) = mesh0%kvert(3,ielem)
  mesh1%kvert(2,iel2) = mesh0%kedge(3,ielem) + mesh0%nvt
  mesh1%kvert(3,iel2) = midpointA
  mesh1%kvert(4,iel2) = mesh0%kedge(2,ielem) + mesh0%nvt
  mesh1%kvert(5,iel2) = mesh0%kedge(7,ielem) + mesh0%nvt
  mesh1%kvert(6,iel2) = midpointD
  mesh1%kvert(7,iel2) = ielmid
  mesh1%kvert(8,iel2) = midpointC

  ! 4th new element

  mesh1%kvert(1,iel3) = mesh0%kvert(4,ielem)
  mesh1%kvert(2,iel3) = mesh0%kedge(4,ielem) + mesh0%nvt
  mesh1%kvert(3,iel3) = midpointA
  mesh1%kvert(4,iel3) = mesh0%kedge(3,ielem) + mesh0%nvt

  mesh1%kvert(5,iel3) = mesh0%kedge(8,ielem) + mesh0%nvt
  mesh1%kvert(6,iel3) = midpointE
  mesh1%kvert(7,iel3) = ielmid
  mesh1%kvert(8,iel3) = midpointD

  ! 5th new element

  mesh1%kvert(1,iel4) = mesh0%kvert(5,ielem)
  mesh1%kvert(4,iel4) = mesh0%kedge(12,ielem) + mesh0%nvt
  mesh1%kvert(3,iel4) = midpointF
  mesh1%kvert(2,iel4) = mesh0%kedge(9,ielem) + mesh0%nvt

  mesh1%kvert(5,iel4) = mesh0%kedge(5,ielem) + mesh0%nvt
  mesh1%kvert(8,iel4) = midpointE
  mesh1%kvert(7,iel4) = ielmid
  mesh1%kvert(6,iel4) = midpointB

  ! 6th new element

  mesh1%kvert(1,iel5) = mesh0%kvert(6,ielem)
  mesh1%kvert(4,iel5) = mesh0%kedge(9,ielem) + mesh0%nvt
  mesh1%kvert(3,iel5) = midpointF
  mesh1%kvert(2,iel5) = mesh0%kedge(10,ielem) + mesh0%nvt

  mesh1%kvert(5,iel5) = mesh0%kedge(6,ielem) + mesh0%nvt
  mesh1%kvert(8,iel5) = midpointB
  mesh1%kvert(7,iel5) = ielmid
  mesh1%kvert(6,iel5) = midpointC

  ! 7th new element

  mesh1%kvert(1,iel6) = mesh0%kvert(7,ielem)
  mesh1%kvert(4,iel6) = mesh0%kedge(10,ielem) + mesh0%nvt
  mesh1%kvert(3,iel6) = midpointF
  mesh1%kvert(2,iel6) = mesh0%kedge(11,ielem) + mesh0%nvt

  mesh1%kvert(5,iel6) = mesh0%kedge(7,ielem) + mesh0%nvt
  mesh1%kvert(8,iel6) = midpointC
  mesh1%kvert(7,iel6) = ielmid
  mesh1%kvert(6,iel6) = midpointD

  ! 8th new element

  mesh1%kvert(1,iel7) = mesh0%kvert(8,ielem)
  mesh1%kvert(4,iel7) = mesh0%kedge(11,ielem) + mesh0%nvt
  mesh1%kvert(3,iel7) = midpointF
  mesh1%kvert(2,iel7) = mesh0%kedge(12,ielem) + mesh0%nvt

  mesh1%kvert(5,iel7) = mesh0%kedge(8,ielem) + mesh0%nvt
  mesh1%kvert(8,iel7) = midpointD
  mesh1%kvert(7,iel7) = ielmid
  mesh1%kvert(6,iel7) = midpointE

  end do

  !---------------------------------
  if(.not.allocated(mesh1%knpr))then
    allocate(mesh1%knpr(mesh1%nvt))
  end if
  mesh1%knpr=0

end subroutine

!================================================================================================
!                                 Sub: genMeshStructures  
!================================================================================================
subroutine genMeshStructures(mesh, extended)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh), intent(inout) :: mesh
  logical, optional, intent(in) :: extended 

  ! local variables
  real :: ttt0=0.0,ttt1=0.0
  logical :: calculateExtendedConnectivity 

  if(present(extended))then
    calculateExtendedConnectivity = extended
  else
    calculateExtendedConnectivity = .false.
  end if

  CALL ZTIME(TTT0)

  call genKVEL(mesh) 

  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for KVEL : ',TTGRID
  END IF

  CALL ZTIME(TTT0)

  call genKEDGE2(mesh)

  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for KEDGE : ',TTGRID
  END IF

  !call genKADJ(mesh)

  CALL ZTIME(TTT0)

  call genKADJ2(mesh)

  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for KADJ2 : ',TTGRID
  END IF

  CALL ZTIME(TTT0)

  call genKAREA(mesh)

  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for KAREA : ',TTGRID
  END IF

  CALL ZTIME(TTT0)

  call genKVAR(mesh)

  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for KVAR : ',TTGRID
  END IF

  CALL ZTIME(TTT0)

  call genDCORAG(mesh)
  CALL ZTIME(TTT1)
  TTGRID=TTT1-TTT0

  IF (myid.eq.showid) THEN
    WRITE(*,*) 'time for DCORAG : ',TTGRID
  END IF

  if(calculateExtendedConnectivity)then
    CALL ZTIME(TTT0)
    call tria_genElementsAtVertex3D(mesh)
    CALL ZTIME(TTT1)
    TTGRID=TTT1-TTT0

    IF (myid.eq.showid) THEN
      WRITE(*,*) 'time for elementsAtVertex : ',TTGRID
    END IF
  end if

end subroutine

!================================================================================================
!                                 Sub: genKVEL  
!================================================================================================
subroutine genKVEL(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: i,j,k
  integer :: ive
  integer :: iglobal, itemp
  integer, allocatable, dimension(:) :: nvel_temp

  mesh%nvel = 0

  allocate(nvel_temp(mesh%nvt))

  nvel_temp = 0

  do i=1,mesh%nel
    do ive=1,mesh%nve
      iglobal = mesh%kvert(ive,i) 
      nvel_temp(iglobal) = nvel_temp(iglobal) + 1
      mesh%nvel=max(mesh%nvel,nvel_temp(iglobal))
    end do
  end do

  if(.not.allocated(mesh%kvel))then
    allocate(mesh%kvel(mesh%nvel,mesh%nvt))
    mesh%kvel = 0
  end if

  do i=1,mesh%nel
    do ive=1,mesh%nve
      iglobal = mesh%kvert(ive,i) 
      do j=1,mesh%nvel
        if(mesh%kvel(j,iglobal).eq.0)then
          mesh%kvel(j,iglobal) = i
          exit
        end if
      end do
    end do
  end do

  ! Sort kvel
  do i = 1,mesh%nvt
    do j = 1,mesh%nvel
      do k = j+1,mesh%nvel 
        if(mesh%kvel(j,i).gt.mesh%kvel(k,i).and.&
          (.not.mesh%kvel(k,i).eq.0))then
          itemp = mesh%kvel(j,i)
          mesh%kvel(j,i) = mesh%kvel(k,i)
          mesh%kvel(k,i) = itemp
        end if
      end do
    end do 
  end do

end subroutine

!================================================================================================
!                                 Sub: genKADJ  
!================================================================================================
subroutine genKADJ(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: iiel,iiar
  integer :: jiel,jiar
  integer :: ive

  integer :: ive1
  integer :: ive2
  integer :: ive3
  integer :: ive4

  integer :: ivt1
  integer :: ivt2
  integer :: ivt3
  integer :: ivt4

  integer :: jve1
  integer :: jve2
  integer :: jve3
  integer :: jve4

  integer :: jvt1
  integer :: jvt2
  integer :: jvt3
  integer :: jvt4

  logical :: bfound = .false.

  if(.not.allocated(mesh%kadj))then
    allocate(mesh%kadj(mesh%nae,mesh%nel))
  end if

  mesh%kadj = -1
  mesh%nat  =  0

  do iiel=1,mesh%nel
  do iiar=1,mesh%nae

  if(mesh%kadj(iiar,iiel).ge.0)cycle

  mesh%nat  =  mesh%nat + 1
  ive1=kiad(1,iiar)
  ive2=kiad(2,iiar)
  ive3=kiad(3,iiar)
  ive4=kiad(4,iiar)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  bfound = .false.
  do jiel=1,mesh%nel

  if(bfound)exit
  if(jiel.eq.iiel)cycle
  do jiar=1,mesh%nae

  jve1=kiad(1,jiar)
  jve2=kiad(2,jiar)
  jve3=kiad(3,jiar)
  jve4=kiad(4,jiar)

  jvt1=mesh%kvert(jve1,jiel)
  jvt2=mesh%kvert(jve2,jiel)
  jvt3=mesh%kvert(jve3,jiel)
  jvt4=mesh%kvert(jve4,jiel)

  if (((jvt1.eq.ivt1).or.(jvt2.eq.ivt1).or.&
    (jvt3.eq.ivt1).or.(jvt4.eq.ivt1)).and.&
    ((jvt1.eq.ivt2).or.(jvt2.eq.ivt2).or.&
    (jvt3.eq.ivt2).or.(jvt4.eq.ivt2)).and.&
    ((jvt1.eq.ivt3).or.(jvt2.eq.ivt3).or.&
    (jvt3.eq.ivt3).or.(jvt4.eq.ivt3)).and.&
    ((jvt1.eq.ivt4).or.(jvt2.eq.ivt4).or.&
    (jvt3.eq.ivt4).or.(jvt4.eq.ivt4))) then

    ! We found a matching face
    mesh%kadj(iiar,iiel)=jiel
    mesh%kadj(jiar,jiel)=iiel
    bfound = .true.

  end if

  end do ! jiar
  end do ! jiel

  ! if there is no matching neighbour then
  ! we have a boundary face
  if(.not.bfound)mesh%kadj(iiar,iiel)=0

  end do
  end do

end subroutine

!================================================================================================
!                                    Sub: genKADJ2  
!================================================================================================
subroutine genKADJ2(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: iiel,iiar
  integer :: j,k
  integer :: iElements

  ! the list of connectors
  type(t_connector3D), dimension(:), pointer :: p_IConnectList

  if(.not.allocated(mesh%kadj))then
    allocate(mesh%kadj(mesh%nae,mesh%nel))
  end if

  mesh%kadj = 0

  iElements = mesh%nel * 6
  allocate(p_IConnectList(iElements))

  call buildConnectorList(p_IConnectList, mesh)

  ! ConnectorList is build, now sort it
  call tria_sortElements3DInt(p_IConnectList, iElements)
  call tria_sortElements3D(p_IConnectList, iElements, 5)

  ! assign the neighbours at elements
  ! traverse the connector list
  do iiel = 2, iElements

  ! check for equivalent connectors... that means:
  ! check if all 4 vertices that define the face are equal.
  ! For mixed triangulations the fourth vertex may be zero but
  ! this is the case for both items of the connector list
  j = 0
  do while(p_IConnectList(iiel-1)%I_conData(j+1) .eq. &
    p_IConnectList(iiel)%I_conData(j+1) )
  ! increment counter
  j = j+1
  end do

  ! assign information
  if(j .eq. 4) then

    mesh%kadj(p_IConnectList(iiel-1)%I_conData(6),&
      p_IConnectList(iiel-1)%I_conData(5)) = & 
      p_IConnectList(iiel)%I_conData(5)

    mesh%kadj(p_IConnectList(iiel)%I_conData(6), &
      p_IConnectList(iiel)%I_conData(5)) = &
      p_IConnectList(iiel-1)%I_conData(5)
  end if

  end do

  ! free list of connectors
  deallocate(p_IConnectList)

  !if(myid.eq.1)then

!        k = 0 
!        do iiel=1,mesh%nel
!          do iiar=1,mesh%nae
!            if(mesh%kadj(iiar,iiel) .ne. mesh%kadj2(iiar,iiel))then
!              write(*,*)'not equal'
!              write(*,'(I5,A,I5)')mesh%kadj(iiar,iiel),':',mesh%kadj2(iiar,iiel)
!              k=k+1
!            end if
!          end do
!        end do
!        write(*,*)'not equal elements,myid: ',k,iElements,myid

  !end if

end subroutine genKADJ2

!================================================================================================
!                                    Sub: genKEDGE  
!================================================================================================
subroutine genKEDGE(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  implicit none
  type(tMesh) :: mesh
  integer :: i,j,k
  integer :: iv1,iv2,ivt1,ivt2,iedge
  integer :: smiel
  integer :: lnet

  if(.not.allocated(mesh%kedge))then
    allocate(mesh%kedge(mesh%nee,mesh%nel))
  end if

  mesh%net = 0
  do i=1,mesh%nel
  do j=1,mesh%nee
  iv1 = KIV(1,j)
  iv2 = KIV(2,j)

  ivt1 = mesh%kvert(iv1,i)
  ivt2 = mesh%kvert(iv2,i)

  call findSmallestIEL(ivt1,ivt2,mesh,i,smiel)

  if(i.gt.smiel)then
    cycle
  end if

  mesh%net   = mesh%net + 1

  end do
  end do

  if(.not.allocated(mesh%kved))then
    allocate(mesh%kved(2,mesh%net))
  end if

  mesh%kved=-1

  lnet = 0
  do i=1,mesh%nel
  do j=1,mesh%nee
  iv1 = KIV(1,j)
  iv2 = KIV(2,j)

  ivt1 = mesh%kvert(iv1,i)
  ivt2 = mesh%kvert(iv2,i)

  call findSmallestIEL(ivt1,ivt2,mesh,i,smiel)

  if(i.ne.smiel)then
    do k=1,mesh%nee 
    iedge = mesh%kedge(k,smiel)
    if(((mesh%kved(1,iedge).eq.ivt1).and.(mesh%kved(2,iedge).eq.ivt2)).or.&
      ((mesh%kved(1,iedge).eq.ivt2).and.(mesh%kved(2,iedge).eq.ivt1)))then

      ! we have found the edge
      mesh%kedge(j,i)=iedge

    end if
    end do

  else

    lnet   = lnet + 1
    do k=1,lnet
    if(((mesh%kved(1,k).eq.ivt1).and.(mesh%kved(2,k).eq.ivt2)).or.&
      ((mesh%kved(1,k).eq.ivt2).and.(mesh%kved(2,k).eq.ivt1)))then
      stop
    end if
    end do
    mesh%kved(1,lnet) = ivt1
    mesh%kved(2,lnet) = ivt2
    mesh%kedge(j,i)=lnet
  end if

  end do
  end do

contains

subroutine findSmallestIEL(i1,i2,mesh,ciel,siel)
  implicit none
  integer :: i1,i2
  type(tMesh) :: mesh
  integer :: ciel
  integer :: siel
  integer :: iel1,iel2
  integer :: i,j

  ! set the smallest iel to the current
  siel = ciel
  ! loop over elements at vertex i1
  do i=1,mesh%nvel
  iel1 = mesh%kvel(i,i1)
!        write(*,*)'element at vertex1:',i1,iel1
  ! if there are no more elements at the vertex: stop
  if(iel1.eq.0)return

  ! loop over the elements at vertex i2
  do j=1,mesh%nvel
  iel2 = mesh%kvel(j,i2)
!          write(*,*)'element at vertex2:',i2,iel2
  ! if there are no more elements at the vertex: exit
  if(iel2.eq.0)exit
  if((iel1.eq.iel2).and.(iel1.lt.siel))then
    siel = iel1
  end if

  end do
  end do

end subroutine

end subroutine

!================================================================================================
!                                    Sub: genKEDGE2  
!================================================================================================
subroutine genKEDGE2(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  implicit none
  type(tMesh) :: mesh
  integer :: i,j,k,smiel_edge
  integer :: iv1,iv2,ivt1,ivt2
  integer :: smiel
  integer :: iedge
  integer, allocatable, dimension(:,:,:) :: verticesAtEdge

  if(.not.allocated(mesh%kedge))then
    allocate(mesh%kedge(mesh%nee,mesh%nel))
  end if

  if(.not.allocated(verticesAtEdge))then
    allocate(verticesAtEdge(mesh%nee,mesh%nel,4))
  end if

  verticesAtEdge = 0

  mesh%net = 0
  iedge = 0

  do i=1,mesh%nel
  do j=1,mesh%nee
  iv1 = KIV(1,j)
  iv2 = KIV(2,j)

  ivt1 = mesh%kvert(iv1,i)
  ivt2 = mesh%kvert(iv2,i)

  call findSmallestIEL(ivt1,ivt2,mesh,i,smiel)

  verticesAtEdge(j,i,1)=ivt1
  verticesAtEdge(j,i,2)=ivt2
  verticesAtEdge(j,i,3)=smiel
  verticesAtEdge(j,i,4)=-1

  if(i.gt.smiel)then
    cycle
  end if

  iedge   = iedge + 1

  verticesAtEdge(j,i,4)=iedge

  end do
  end do

  mesh%net = iedge

  if(.not.allocated(mesh%kved))then
    allocate(mesh%kved(2,mesh%net))
  end if

  mesh%kved=-1

  do i=1,mesh%nel
  do j=1,mesh%nee

  if(verticesAtEdge(j,i,3) .eq. i)then
    iedge = verticesAtEdge(j,i,4)

    mesh%kved(1,iedge) = verticesAtEdge(j,i,1)
    mesh%kved(2,iedge) = verticesAtEdge(j,i,2)
    mesh%kedge(j,i)    = iedge
  else

    iv1 = KIV(1,j)
    iv2 = KIV(2,j)

    ivt1 = mesh%kvert(iv1,i)
    ivt2 = mesh%kvert(iv2,i)
    smiel = verticesAtEdge(j,i,3)

    do k=1,mesh%nee

    smiel_edge = mesh%kedge(k,smiel)

    if(( (mesh%kved(1,smiel_edge).eq.ivt1).and.&
      (mesh%kved(2,smiel_edge).eq.ivt2)).or.&
      ( (mesh%kved(1,smiel_edge).eq.ivt2).and.&
      (mesh%kved(2,smiel_edge).eq.ivt1)))then

      mesh%kedge(j,i)=smiel_edge
      exit              
    end if

    end do

!            if((myid.eq.0).and.(i.eq.7).and.(j.eq.7))then
!              write(*,*)'else..',k,smiel
!              call findSmallestIEL(ivt1,ivt2,mesh,i,smiel)
!              write(*,*)'else..',k,smiel
!            end if

  end if

  end do
  end do

!      do i=1,mesh%nel
!        do j=1,mesh%nee
!        if((mesh%kedge(j,i) .lt. 0).or.(mesh%kedge(j,i) .gt. mesh%net))then
!          write(*,*)'error: kedge,myid',mesh%kedge(j,i),i,j,mesh%net,myid
!          stop
!        end if
!
!        end do
!      end do

contains

subroutine findSmallestIEL(i1,i2,mesh,ciel,siel)
  implicit none
  integer :: i1,i2
  type(tMesh) :: mesh
  integer :: ciel
  integer :: siel
  integer :: iel1,iel2
  integer :: i,j

  ! set the smallest iel to the current
  siel = ciel
  ! loop over elements at vertex i1
  do i=1,mesh%nvel
  iel1 = mesh%kvel(i,i1)

  ! if there are no more elements at the vertex: stop
  if(iel1.eq.0)return

  ! loop over the elements at vertex i2
  do j=1,mesh%nvel
  iel2 = mesh%kvel(j,i2)

  ! if there are no more elements at the vertex: exit
  if(iel2.eq.0)exit
  if((iel1.eq.iel2).and.(iel1.lt.siel))then
    siel = iel1
  end if

  end do
  end do

end subroutine

end subroutine genKEDGE2

!================================================================================================
!                                    Sub: genKVAR  
!================================================================================================
subroutine genKVAR(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: iiel,iface
  integer :: jiel,jface
  integer :: ive

  integer :: ive1
  integer :: ive2
  integer :: ive3
  integer :: ive4

  integer :: ivt1
  integer :: ivt2
  integer :: ivt3
  integer :: ivt4

  integer :: jve1
  integer :: jve2
  integer :: jve3
  integer :: jve4

  integer :: jvt1
  integer :: jvt2
  integer :: jvt3
  integer :: jvt4
  integer :: ifaceGlobal
  integer :: ifaceNumber


  if(.not.allocated(mesh%kvar))then
    allocate(mesh%kvar(4,mesh%nat))
  end if

  ifaceGlobal =  0

  do iiel=1,mesh%nel
  do iface=1,mesh%nae

  if((mesh%kadj(iface,iiel).eq.0).or.&
    (mesh%kadj(iface,iiel) > iiel))then

    ifaceGlobal = ifaceGlobal + 1 

    do ive=1,4
    mesh%kvar(ive,ifaceGlobal) = &
      mesh%kvert(kiad(ive,iface),iiel)
    end do

  else

    jiel = mesh%kadj(iface,iiel) 

    jface = 1
    do while(iiel.ne.mesh%kadj(jface,jiel))
    jface = jface + 1
    end do

    ifaceNumber = mesh%karea(jface,jiel)

    do ive=1,4
    mesh%kvar(ive,ifaceNumber) = &
      mesh%kvert(kiad(ive,iface),iiel)
    end do

  end if

  end do
  end do

end subroutine

!================================================================================================
!                                    Sub: genDCORAG  
!================================================================================================
subroutine genDCORAG(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: iface
  integer :: ive
  integer :: ivt

  if(.not.allocated(mesh%dcorag))then
    allocate(mesh%dcorag(3,mesh%nat))
  end if

  do iface=1,mesh%nat

  mesh%dcorag(1:3,iface) = 0d0
  do ive=1,4
  ivt = mesh%kvar(ive,iface)
  mesh%dcorag(1:3,iface) = mesh%dcorag(1:3,iface) +&
    mesh%dcorvg(1:3,ivt)   
  end do
  mesh%dcorag(1:3,iface) = 0.25d0 * mesh%dcorag(1:3,iface)

  end do

end subroutine

!================================================================================================
!                                    Sub: genKAREA  
!================================================================================================
subroutine genKAREA(mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE
  type(tMesh) :: mesh
  integer :: iiel,iface
  integer :: jiel,jface
  integer :: ive

  integer :: ive1
  integer :: ive2
  integer :: ive3
  integer :: ive4

  integer :: ivt1
  integer :: ivt2
  integer :: ivt3
  integer :: ivt4

  integer :: jve1
  integer :: jve2
  integer :: jve3
  integer :: jve4

  integer :: jvt1
  integer :: jvt2
  integer :: jvt3
  integer :: jvt4
  integer :: ifaceGlobal

  logical :: bfound = .false.


  if(.not.allocated(mesh%karea))then
    allocate(mesh%karea(mesh%nae,mesh%nel))
  end if


  ifaceGlobal =  0
  mesh%karea=0

  do iiel=1,mesh%nel
    do iface=1,mesh%nae

      if((mesh%kadj(iface,iiel).eq.0).or.&
        (mesh%kadj(iface,iiel) > iiel))then

        ifaceGlobal = ifaceGlobal + 1 

        mesh%karea(iface, iiel) = ifaceGlobal
      else
        jiel = mesh%kadj(iface,iiel) 

        do jface=1,6

        if(iiel.eq.mesh%kadj(jface,jiel))then
          mesh%karea(iface,iiel) = mesh%karea(jface,jiel)
          exit
        end if
        end do

      end if

    end do
  end do

  mesh%nat = ifaceGlobal

end subroutine

subroutine buildConnectorList(IConnectList, mesh)
  use PP3D_MPI, only:myid,showid
  use var_QuadScalar
  IMPLICIT NONE

  type(tMesh) :: mesh

  ! the list of connectors this routine is supposed to build
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList

  ! local variables
  integer :: iiel,k,nfaces

  integer :: ive1
  integer :: ive2
  integer :: ive3
  integer :: ive4

  integer :: ivt1
  integer :: ivt2
  integer :: ivt3
  integer :: ivt4

  ! function body

  ! initialise the number of faces
  nfaces = 0

  ! loop through all elements
  do iiel = 1, mesh%NEL

  ! build connectors for each hexahedron

  !=========================================================
  ! first face
  nfaces = nfaces+1

  ive1=kiad(1,1)
  ive2=kiad(2,1)
  ive3=kiad(3,1)
  ive4=kiad(4,1)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4

  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 1

  !=========================================================
  ! sixth face
  nfaces = nfaces+1

  ive1=kiad(1,6)
  ive2=kiad(2,6)
  ive3=kiad(3,6)
  ive4=kiad(4,6)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4

  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 6

  !=========================================================
  ! second face
  nfaces = nfaces+1

  ive1=kiad(1,2)
  ive2=kiad(2,2)
  ive3=kiad(3,2)
  ive4=kiad(4,2)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4

  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 2

  !=========================================================
  ! fourth face
  nfaces = nfaces+1

  ive1=kiad(1,4)
  ive2=kiad(2,4)
  ive3=kiad(3,4)
  ive4=kiad(4,4)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4

  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 4

  !=========================================================
  ! third face
  nfaces = nfaces+1

  ive1=kiad(1,3)
  ive2=kiad(2,3)
  ive3=kiad(3,3)
  ive4=kiad(4,3)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4

  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 3

  !=========================================================
  ! fifth face
  nfaces = nfaces+1

  ive1=kiad(1,5)
  ive2=kiad(2,5)
  ive3=kiad(3,5)
  ive4=kiad(4,5)

  ivt1=mesh%kvert(ive1,iiel)
  ivt2=mesh%kvert(ive2,iiel)
  ivt3=mesh%kvert(ive3,iiel)
  ivt4=mesh%kvert(ive4,iiel)

  IConnectList(nfaces)%I_conData(1) = ivt1
  IConnectList(nfaces)%I_conData(2) = ivt2
  IConnectList(nfaces)%I_conData(3) = ivt3
  IConnectList(nfaces)%I_conData(4) = ivt4


  ! save the number of the element this face was found from
  IConnectList(nfaces)%I_conData(5) = iiel

  ! assign the local face number
  IConnectList(nfaces)%I_conData(6) = 5

  !=========================================================

  end do

end subroutine buildConnectorList

!================================================================================================
!                                    Sub: tria_sortElements3D  
!================================================================================================
subroutine tria_genElementsAtVertex3D(mesh)
  use types
  implicit none
  ! The triangulation structure to be updated.
  type(tMesh), intent(inout) :: mesh

  ! local variables
  integer :: iiel
  integer :: iive
  integer :: iivt, isize
  integer, allocatable, dimension(:) :: auxArray

  allocate(mesh%elementsAtVertexIdx(mesh%NVT+1)) 

  mesh%elementsAtVertexIdx = 0

  ! first we calculate the number of elements at each vertex simply
  ! by counting; thus, loop over all elements
  do iiel = 1, mesh%NEL

    ! loop over all vertices at the element
    do iive = 1, mesh%NVE

      ! ivt is the ive-th vertex at element iel
      iivt = mesh%kvert(iive,iiel)

      ! increase the number of elements by one
      mesh%elementsAtVertexIdx(iivt+1) = mesh%elementsAtVertexIdx(iivt+1) + 1

    end do ! end ive
  end do ! end iel

  mesh%NNelAtVertex = 0

  ! In the next step we sum up the number of elements at two
  ! successive vertices to create the index array and find the 
  ! length of elementsAtVertex. Calculate NNelAtVertex.
  do iivt = 2, mesh%NVT+1
    mesh%NNelAtVertex = max(mesh%NNelAtVertex, &
                            mesh%elementsAtVertexIdx(iivt))

    mesh%elementsAtVertexIdx(iivt) = mesh%elementsAtVertexIdx(iivt) + &
                                    mesh%elementsAtVertexIdx(iivt-1)
  end do

  ! set the size
  isize = mesh%elementsAtVertexIdx(mesh%NVT+1)
  mesh%elementsAtVertexIdx(1) = 1

  ! Isize contains now the length of the array where we store the
  ! adjacency information (IelementsAtVertex).  Do we have (enough)
  ! memory for that array?
  allocate(mesh%elementsAtVertex(isize))


  ! Duplicate the elementsAtVertexIdx array. We use that as
  ! pointer and index if new elements at a vertex are found.
  allocate(auxArray(mesh%NVT+1)) 
  auxArray = mesh%elementsAtVertexIdx 

  ! loop over all elements
  do iiel = 1, mesh%NEL

    ! loop over all vertices of the element
    do iive = 1, mesh%NVE

      ! ivt is the ive-th vertex at element iel
      iivt = mesh%kvert(iive,iiel)

      ! store the adjacency information at position p_Iaux1(ivt)
      mesh%elementsAtVertex( auxArray(iivt) ) = iiel

      ! increase the position of the next element in p_Iaux1(ivt)
      auxArray(iivt) = auxArray(iivt) + 1

    end do ! end iive
  end do ! end iiel 

end subroutine tria_genElementsAtVertex3D

!================================================================================================
!                                    Sub: tria_sortElements3D  
!================================================================================================
subroutine tria_sortElements3D(IConnectList, iElements, components)

  ! This subroutine establishes the lexicographic
  ! ordering on the list of connectors in 3D

  integer, intent(in) :: iElements

  integer, intent(in) :: components

  type(t_connector3D), dimension(:), intent(inout) :: IConnectList

  ! local
  integer :: j

  do j = components, 1, -1
  call tria_mergesort(IConnectList, 1, iElements, j)
  end do

end subroutine tria_sortElements3D

!================================================================================================
!                                    Sub: tria_sortElements3DInt  
!================================================================================================
subroutine tria_sortElements3DInt(IConnectList, iElements)

  ! This subroutine establishes the internal sorted numbering
  ! on the entries of the connector list in 3D

  ! parameter values

  integer, intent(in) :: iElements

  type(t_connector3D), dimension(:), intent(inout) :: IConnectList

  ! local variables
  integer :: i

  ! create a sorted numbering in all connectors
  do i = 1, iElements
  call sort(IConnectList(i)%I_conData(1:4))
  end do

contains

! ---------------------------------------------------------------

pure subroutine sort(Idata)
  integer, intent(inout), dimension(4) :: Idata

  if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
  if (Idata(3) < Idata(2)) call swap(Idata(3), Idata(2))
  if (Idata(4) < Idata(3)) call swap(Idata(4), Idata(3))
  if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
  if (Idata(3) < Idata(2)) call swap(Idata(3), Idata(2))
  if (Idata(2) < Idata(1)) call swap(Idata(2), Idata(1))
end subroutine sort

! ---------------------------------------------------------------

elemental pure subroutine swap(a,b)
  integer, intent(inout) :: a,b

  ! local variables
  integer :: c

  c = a
  a = b
  b = c
end subroutine swap

end subroutine tria_sortElements3DInt

!************************************************************************

recursive subroutine tria_mergesort(IConnectList, l, r, pos)

  ! This routine sorts a connector list it is used as an
  ! auxilliary routine during the Neighbours at elements routine

  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer, intent(in) :: l,r,pos

  ! the list of connectors
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList

  ! local variables
  integer :: m

  if(l < r) then

    m = l + (r-l)/2

    call tria_mergesort(IConnectList, l,   m, pos)
    call tria_mergesort(IConnectList, m+1, r, pos)
    call tria_merge(IConnectList, l, m, r, pos)

  end if

end subroutine tria_mergesort

subroutine tria_merge(IConnectList, l, m, r, pos)

  ! the array positions l...r will be sorted
  ! the sorting key is element 'pos' of the connector
  integer, intent(in) :: l,r,m,pos

  ! the list of connectors
  type(t_connector3D), dimension(:), intent(inout) :: IConnectList


  ! local variables
  integer :: i,j,n1,n2,k
  type(t_connector3D), dimension(:), pointer :: p_L, p_R

  ! init counters
  n1 = m - l + 1

  n2 = r - m

  k = l

  ! allocate memory for merging
  allocate(p_L(n1))
  allocate(p_R(n2))

  ! fill left array
  do i = 1, n1
  p_L(i) = IConnectList(l+i-1)
  end do

  ! fill right array
  do j = 1, n2
  p_R(j) = IConnectList(m+j)
  end do

  i = 1
  j = 1

  ! merge
  do
  if( (i > n1 ) .or. (j > n2) ) exit

  ! if the current element of the left array is smaller
  ! copy it to p_ConnectorList
  ! else
  ! copy the element from the right array
  if(p_L(i)%I_conData(pos) .le. p_R(j)%I_conData(pos)) then
    IConnectList(k) = p_L(i)
    i = i + 1
    k = k + 1
  else
    IConnectList(k) = p_R(j)
    j = j + 1
    k = k + 1
  end if

  end do

  ! copy the remaining entries of p_L (if present)
  do
  if(i > n1) exit

  IConnectList(k) = p_L(i)
  ! increment counters
  k = k + 1
  i = i + 1

  end do

  ! copy the remaining entries of p_R (if present)
  do
  if(j > n2) exit

  IConnectList(k) = p_R(j)
  ! increment counters
  k = k + 1
  j = j + 1

  end do

  ! done merging

  ! free p_L and p_R
  deallocate(p_L)
  deallocate(p_R)

end subroutine tria_merge

!================================================================================================
!                                  Sub: release_mesh 
!================================================================================================
subroutine release_mesh(mgMesh) 

use var_QuadScalar, only: tMultiMesh 
implicit none

type(tMultiMesh) :: mgMesh
integer :: maxlevel

maxlevel = mgMesh%maxlevel

if(associated(mgMesh%level(maxlevel)%dcorvg))then
  deallocate(mgMesh%level(maxlevel)%dcorvg)
  mgMesh%level(maxlevel)%dcorvg => null()
end if

end subroutine release_mesh

!================================================================================================
!                                  Sub: release_mesh 
!================================================================================================
subroutine testElementsAtVertex(mesh) 

use types 
implicit none

type(tMesh) :: mesh

integer :: ive

do ive = 1, mesh%NVT
  write(*,*)'Vertex: ', ive,'Number of elements:', mesh%elementsAtVertexIdx(ive+1) - mesh%elementsAtVertexIdx(ive)
end do

end subroutine testElementsAtVertex
!================================================================================================
!                                  Sub: release_mesh 
!================================================================================================
subroutine SETARE(DVOL,NEL,KVERT,DCORVG)
implicit none
!***********************************************************************
!
!   Purpose: - writes on  AVOL(IEL)  the VOLUME of the element IEL,
!              IEL=1,...,NEL
!            - writes on  AVOL(NEL+1) the sum of all  AVOL(IEL)
!            - KVERT,DCORVG are the usual FEAT arrays
!
!***********************************************************************
!=======================================================================
!     Declarations
!=======================================================================
integer, parameter :: nnve=8
real*8,parameter :: a1=1d0/6d0
!      
integer, intent(in) :: nel
real*8, intent(inout), dimension(:) :: dvol
real*8, intent(in), dimension(:,:) :: dcorvg
integer, intent(in), dimension(:,:) :: kvert

! locals    
!=======================================================================
real*8 :: sum,aaa
integer :: i1,i2,i3,i4,i5,i6,i7,i8,iel
real*8  :: x1,x2,x3,x4,x5,x6,x7,x8
real*8  :: y1,y2,y3,y4,y5,y6,y7,y8
real*8  :: z1,z2,z3,z4,z5,z6,z7,z8

sum = 0.0

DO IEL=1,NEL
 
 I1=KVERT(1,IEL)
 I2=KVERT(2,IEL)
 I3=KVERT(3,IEL)
 I4=KVERT(4,IEL)
 I5=KVERT(5,IEL)
 I6=KVERT(6,IEL)
 I7=KVERT(7,IEL)
 I8=KVERT(8,IEL)
 
 X1=DCORVG(1,I1)
 X2=DCORVG(1,I2)
 X3=DCORVG(1,I3)
 X4=DCORVG(1,I4)
 X5=DCORVG(1,I5)
 X6=DCORVG(1,I6)
 X7=DCORVG(1,I7)
 X8=DCORVG(1,I8)
 
 Y1=DCORVG(2,I1)
 Y2=DCORVG(2,I2)
 Y3=DCORVG(2,I3)
 Y4=DCORVG(2,I4)
 Y5=DCORVG(2,I5)
 Y6=DCORVG(2,I6)
 Y7=DCORVG(2,I7)
 Y8=DCORVG(2,I8)
 
 Z1=DCORVG(3,I1)
 Z2=DCORVG(3,I2)
 Z3=DCORVG(3,I3)
 Z4=DCORVG(3,I4)
 Z5=DCORVG(3,I5)
 Z6=DCORVG(3,I6)
 Z7=DCORVG(3,I7)
 Z8=DCORVG(3,I8)
 
 AAA=A1*((DABS((X4-X1)*(Y4-Y3)*(Z4-Z8)+(Y4-Y1)*  &
        (Z4-Z3)*(X4-X8)+(Z4-Z1)*(X4-X3)*(Y4-Y8)- &
        (X4-X8)*(Y4-Y3)*(Z4-Z1)-(Y4-Y8)*(Z4-Z3)* &
        (X4-X1)-(Z4-Z8)*(X4-X3)*(Y4-Y1)))+       &
        (DABS((X2-X3)*(Y2-Y1)*(Z2-Z6)+(Y2-Y3)*   &
        (Z2-Z1)*(X2-X6)+(Z2-Z3)*(X2-X1)*(Y2-Y6)- &
        (X2-X6)*(Y2-Y1)*(Z2-Z3)-(Y2-Y6)*(Z2-Z1)* &
        (X2-X3)-(Z2-Z6)*(X2-X1)*(Y2-Y3)))+       &
        (DABS((X5-X8)*(Y5-Y6)*(Z5-Z1)+(Y5-Y8)*   &
        (Z5-Z6)*(X5-X1)+(Z5-Z8)*(X5-X6)*(Y5-Y1)- &
        (X5-X1)*(Y5-Y6)*(Z5-Z8)-(Y5-Y1)*(Z5-Z6)* &
        (X5-X8)-(Z5-Z1)*(X5-X6)*(Y5-Y8)))+       &
        (DABS((X7-X6)*(Y7-Y8)*(Z7-Z3)+(Y7-Y6)*   &
        (Z7-Z8)*(X7-X3)+(Z7-Z6)*(X7-X8)*(Y7-Y3)- &
        (X7-X3)*(Y7-Y8)*(Z7-Z6)-(Y7-Y3)*(Z7-Z8)* &
        (X7-X6)-(Z7-Z3)*(X7-X8)*(Y7-Y6)))+       &
        (DABS((X1-X3)*(Y1-Y8)*(Z1-Z6)+(Y1-Y3)*   &
        (Z1-Z8)*(X1-X6)+(Z1-Z3)*(X1-X8)*(Y1-Y6)- &
        (X1-X6)*(Y1-Y8)*(Z1-Z3)-(Y1-Y6)*(Z1-Z8)* &
        (X1-X3)-(Z1-Z6)*(X1-X8)*(Y1-Y3))))
 DVOL(IEL)=REAL(AAA)
 SUM=SUM+AAA
END DO

DVOL(NEL+1)=SUM

End subroutine SETARE

end Module
