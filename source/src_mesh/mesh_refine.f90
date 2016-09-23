      Module Mesh_Structures

      DIMENSION KIV(2,12)
      DATA KIV /1,2, 2,3, 3,4, 4,1, 1,5, 2,6, &
                3,7, 4,8, 5,6, 6,7, 7,8, 8,5/

      DIMENSION KIAD(4,6)
      DATA KIAD/1,2,3,4, 1,2,6,5, 2,3,7,6, 3,4,8,7, 4,1,5,8, 5,6,7,8/
      !DATA KIAD/1,2,3,4, 1,5,6,2, 2,6,7,3, 3,7,8,4, 1,4,8,5, 5,8,7,6/

      contains

      subroutine writeTriFile(mesh)
      USE var_QuadScalar
      type(tMesh) :: mesh
      integer :: cunit = 123
      integer :: ivert,ielem 

      open (unit=cunit,file='mesh2.tri')

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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
      subroutine writeArrays(mydcorvg, mykvert, mykedge, mykadj,&
      mykarea, mnvt, mnel, mnet,cfile)
      USE var_QuadScalar

      !REAL*8  mydcorvg(3,*)
      !REAL*8, intent(in) ::  mydcorvg(3,:)
      REAL*8  mydcorvg(3,*)

      integer mykvert(8,*)

      integer mykedge(12,*)

      integer mykadj(6,*)

      integer mykarea(6,*)

      CHARACTER (len = 60) :: cfile 

!REAL*8  DCORAG(3,*),DCORVG(3,*)
!INTEGER KVERT(8,*),KEDGE(12,*),KAREA(6,*)

      integer :: mnvt
      integer :: mnel 
      integer :: mnet 

      integer :: cunit = 123
      integer :: ivert,ielem,iedge 

      open (unit=cunit,file=cfile)

      !write(cunit,'(A)')'modified by tr2to3'
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
      subroutine readTriCoarse(CFILE, mgMesh)
        USE var_QuadScalar
        USE PP3D_MPI
        CHARACTER (len = 60) :: CFILE 
        type(tMultiMesh) :: mgMesh
        logical :: bExist
        integer :: cunit = 18
        integer :: i

        INQUIRE(FILE=CFILE,EXIST=bExist)
        IF (.NOT.bExist) THEN
         WRITE(*,*) "File ",CFILE," could not be read ..."
         RETURN
        END IF
        
      !  WRITE(*,*) "Found File ",CFile
        open(cunit,file=cfile)
        read(cunit, *)
        read(cunit, *)
        read(cunit,*),mgMesh%level(1)%nel, mgMesh%level(1)%nvt, &
                      mgMesh%level(1)%nbct, mgMesh%level(1)%nve, &
                      mgMesh%level(1)%nee , mgMesh%level(1)%nae

      !  if(myid.eq.0)then
      !    write(*,*),mgMesh%level(1)%nel, mgMesh%level(1)%nvt, &
      !               mgMesh%level(1)%nbct, mgMesh%level(1)%nve, &
      !               mgMesh%level(1)%nee , mgMesh%level(1)%nae
      !  end if

        read(cunit, *)

        if(.not.allocated(mgMesh%level(1)%dcorvg))then
          allocate(mgMesh%level(1)%dcorvg(3,mgMesh%level(1)%nvt))
        end if

        if(.not.allocated(mgMesh%level(1)%kvert))then
          allocate(mgMesh%level(1)%kvert(8,mgMesh%level(1)%nel))
        end if

        if(.not.allocated(mgMesh%level(1)%knpr))then
          allocate(mgMesh%level(1)%knpr(mgMesh%level(1)%nvt))
        end if

        do i=1,mgMesh%level(1)%nvt
          read(cunit,*),mgMesh%level(1)%dcorvg(1:3,i)
      !    if(myid.eq.0)then
      !      write(*,*),mgMesh%level(1)%dcorvg(1:3,i)
      !    end if
        end do
        read(cunit, *)

        do i=1,mgMesh%level(1)%nel
          read(cunit,*),mgMesh%level(1)%kvert(1:8,i)
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
      subroutine refineMesh(mgMesh,maxlevel)
      use PP3D_MPI, only:myid,showid
      use var_QuadScalar

      type(tMultiMesh) :: mgMesh
      integer :: maxlevel
      integer :: nfine
      integer :: ilevel
      integer :: icurr

      nfine = maxlevel-1
       
      icurr = 1

      call genMeshStructures(mgMesh%level(icurr))

      write(*,*)'Refining to level: ',maxlevel
      write(*,*)'Number of refinement steps needed: ',nfine
      do ilevel=1,nfine
        write(*,*)'Refining to level: ',ilevel + mgMesh%nlmin
        icurr = icurr + 1 
        call refineMeshLevel(mgMesh%level(icurr-1), mgMesh%level(icurr))
        call genMeshStructures(mgMesh%level(icurr))
      end do

      end subroutine
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
      subroutine refineMeshLevel(mesh0, mesh1)
      use PP3D_MPI, only:myid,showid
      use var_QuadScalar

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

      if(.not.allocated(mesh1%dcorvg))then
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
      subroutine genMeshStructures(mesh)
      use PP3D_MPI, only:myid,showid
      use var_QuadScalar
      IMPLICIT NONE
      type(tMesh) :: mesh

        call genKVEL(mesh) 

        call genKEDGE(mesh)

        call genKADJ(mesh)

        call genKAREA(mesh)

        call genKVAR(mesh)

        call genDCORAG(mesh)

      end subroutine
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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
!+––––––––––––––––––––––––––––––––––––––––––––––
!| 
!|      
!|      
!+––––––––––––––––––––––––––––––––––––––––––––––
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

end Module
