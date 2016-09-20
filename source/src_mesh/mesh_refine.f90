      Module Mesh_Structures

      DIMENSION KIV(2,12)
      DATA KIV /1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/

      DIMENSION KIAD(4,6)
      DATA KIAD/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

      contains

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

      if(myid.ne.0)then
        write(*,*)'Refining to level: ',maxlevel
        write(*,*)'Number of refinement steps needed: ',nfine
        do ilevel=1,nfine
          icurr = icurr + 1 
          if(myid.eq.1)then
            !call refineMeshLevel(mgMesh%level(icurr-1), mgMesh%level(icurr))
            !call genMeshStructures(mgMesh%level(icurr))
          end if
          return
        end do
      end if

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
      integer :: iivt,ivt1,ivt2
      integer :: iedge,ielem,j,iface
      integer :: iistart,iend

      mesh1%nvt = mesh0%nvt + mesh0%net + mesh0%nat + mesh0%nel

      write(*,*)'New NVT: ',mesh1%nvt
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
        iedge = iedge + 1
      end do

      iistart = mesh0%nvt+mesh0%net+1 
      iend    = mesh0%nvt+mesh0%net+mesh0%nat 
      iface   = 1
      do iivt=iistart,iend
        mesh1%dcorvg(1:3,iivt) = 0d0
        do j=1,4
          mesh1%dcorvg(1:3,iivt) = mesh1%dcorvg(1:3,iivt) + &
                                   mesh0%kvar(j,iface)
        end do

        mesh1%dcorvg(1:3,iivt) = 0.25d0 * mesh1%dcorvg(1:3,iivt)
        iface = iface + 1
      end do

      iistart = mesh0%nvt+mesh0%net+mesh0%nat+1 
      iend   = mesh0%nvt+mesh0%net+mesh0%nat+mesh0%nel 
      ielem = 1
      do iivt=iistart,iend
        mesh1%dcorvg(1:3,iivt) = 0d0
        do j=1,8
          mesh1%dcorvg(1:3,iivt) = mesh1%dcorvg(1:3,iivt) + &
                                   mesh0%kvert(j,ielem)
        end do

        mesh1%dcorvg(1:3,iivt) = 0.125d0 * mesh1%dcorvg(1:3,iivt)
        ielem = ielem + 1
      end do
      
      mesh1%nel = mesh0%nel * 8

      write(*,*)'New NEL: ',mesh1%nel
      allocate(mesh1%kvert(8,mesh1%nel))

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

        call genKVAR(mesh)

        call genKAREA(mesh)

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

      write(*,*)'max nvel:',mesh%nvel

      if(.not.allocated(mesh%kvel))then
        allocate(mesh%kvel(mesh%nvel,mesh%nvt))
        mesh%kvel = 0
      end if

      write(*,*)'Allocated kvel1',size(mesh%kvel,1)
      write(*,*)'Allocated kvel2',size(mesh%kvel,2)

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
            if(mesh%kvel(j,i).gt.mesh%kvel(k,i))then
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

          mesh%nat=mesh%nat+1
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

      write(*,*)'Number of faces: ', mesh%nat

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
      write(*,*)'Number of edges: ',mesh%net

      lnet = 0
      do i=1,mesh%nel
        do j=1,mesh%nee
          iv1 = KIV(1,j)
          iv2 = KIV(2,j)

          ivt1 = mesh%kvert(iv1,i)
          ivt2 = mesh%kvert(iv2,i)

          call findSmallestIEL(ivt1,ivt2,mesh,i,smiel)

          if(i.gt.smiel)then
            do k=1,mesh%nee 
              iedge = mesh%kedge(k,smiel)
              if(((mesh%kved(1,iedge).eq.ivt1).or.(mesh%kved(2,iedge).eq.ivt1)).and.&
                 ((mesh%kved(1,iedge).eq.ivt2).or.(mesh%kved(2,iedge).eq.ivt2)))then

                ! we have found the edge
                mesh%kedge(j,i)=iedge

              end if
            end do
            cycle
          end if

          lnet   = lnet + 1
          mesh%kved(1,lnet) = ivt1
          mesh%kved(2,lnet) = ivt2
          mesh%kedge(j,i)=lnet

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

      siel = ciel
      do i=1,mesh%nvel
        iel1 = mesh%kvel(i,i1)
        if(iel1.eq.0)return
        do j=1,mesh%nvel
          iel2 = mesh%kvel(j,i2)
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
      integer :: nnat
      integer, allocatable, dimension(:,:) :: lkadj

      logical :: bfound = .false.

      if(.not.allocated(mesh%kvar))then
        allocate(mesh%kvar(4,mesh%nat))
      end if

      allocate(lkadj(mesh%nae,mesh%nel))

      lkadj = -1
      nnat  =  0

      do iiel=1,mesh%nel
        do iiar=1,mesh%nae

          if(lkadj(iiar,iiel).ge.0)cycle

          nnat=nnat+1
          ive1=kiad(1,iiar)
          ive2=kiad(2,iiar)
          ive3=kiad(3,iiar)
          ive4=kiad(4,iiar)

          ivt1=mesh%kvert(ive1,iiel)
          ivt2=mesh%kvert(ive2,iiel)
          ivt3=mesh%kvert(ive3,iiel)
          ivt4=mesh%kvert(ive4,iiel)

          mesh%kvar(1,nnat) = ivt1
          mesh%kvar(2,nnat) = ivt2
          mesh%kvar(3,nnat) = ivt3
          mesh%kvar(4,nnat) = ivt4

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
                 lkadj(iiar,iiel)=jiel
                 lkadj(jiar,jiel)=iiel
                 bfound = .true.

             end if

            end do ! jiar
          end do ! jiel

          ! if there is no matching neighbour then
          ! we have a boundary face
          if(.not.bfound)lkadj(iiar,iiel)=0

        end do
      end do

      write(*,*)'Kvar array built, number of faces: ', nnat

      end subroutine

      subroutine genKAREA(mesh)
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
      integer :: ifaceGlobal

      logical :: bfound = .false.

      if(.not.allocated(mesh%karea))then
        allocate(mesh%karea(mesh%nae,mesh%nat))
      end if

      ifaceGlobal =  0

      do iiel=1,mesh%nel
        do iiar=1,mesh%nae

          if((mesh%kadj(iiar,iiel).eq.0).or.&
             (mesh%kadj(iiar,iiel) > iiel))then

            ifaceGlobal = ifaceGlobal + 1 

            mesh%karea(iiar,iiel) = ifaceGlobal

          else

            jiel = mesh%kadj(iiar,iiel)
            jiar = 1
            do while(iiel.ne.mesh%kadj(jiar,jiel))
              jiar = jiar + 1
            end do

            mesh%karea(iiar,iiel) = mesh%kadj(jiar,jiel)

          end if

        end do
      end do

      write(*,*)'Karea array built, number of faces: ', mesh%nat

      end subroutine

end Module
