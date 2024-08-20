 Subroutine FindNodes()
 
  use mpi
  USE PP3D_MPI, ONLY : myid,numnodes,subnodes,MPI_COMM_SUBS,MPI_COMM_SUBGROUP
  use var_QuadScalar, ONLY : myRecComm
  implicit none

  integer :: ierr, hostname_len, max_hostname_len
  integer :: i, myNodeGroup,world_group,new_group
  INTEGER pID,pJD
  character(len=256) :: cFMT

  call MPI_COMM_GROUP(MPI_COMM_WORLD, world_group, ierr)
  
  if (myid.eq.0) return
  
  ! Get the hostname of the current process
  call MPI_Get_processor_name(myRecComm%hostname, hostname_len, ierr)

  ! Determine the maximum hostname length across all processes
  max_hostname_len = 256

  ! Adjust hostname length to max_hostname_len
  myRecComm%hostname = adjustl(myRecComm%hostname)

  ! Allocate myRecComm%all_hostnames array
  allocate(myRecComm%all_hostnames(subnodes))
  allocate(myRecComm%groupIDs(subnodes))

  ! Collect all hostnames at all processes
  call MPI_Allgather(myRecComm%hostname, max_hostname_len, MPI_CHARACTER, myRecComm%all_hostnames, max_hostname_len, MPI_CHARACTER, MPI_COMM_SUBS, ierr)

  ! Determine the unique hostnames
  call unique(myRecComm%all_hostnames, myRecComm%unique_hostnames, myRecComm%hostleaders,myRecComm%hostgroup)
  myRecComm%NumHosts = size(myRecComm%unique_hostnames)
  
  do i=1,myRecComm%NumHosts
   if (adjustl(trim(myRecComm%hostname)).eq.adjustl(trim(myRecComm%unique_hostnames(i)))) then
    myRecComm%myNodeGroup = i
    write(cFMT,*) '(A,I0,A,A,A,I0,A,',size(myRecComm%hostgroup)-1,'(I0,(",")),I0)'
    write(*,cFMT) "myid: ",myid," My Node is :",adjustl(trim(myRecComm%hostname))," My GroupLeader is :",myRecComm%hostleaders(myNodeGroup), " My group is: ",myRecComm%hostgroup
   end if
  end do

  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
  
  !----------------------------------------------------------------------------------------------
  call MPI_GROUP_INCL(world_group, size(myRecComm%hostgroup), myRecComm%hostgroup, new_group, ierr)
  
  ! Create a new communicator for the subgroup
  call MPI_COMM_CREATE(MPI_COMM_SUBS, new_group, MPI_COMM_SUBGROUP, ierr)
  !----------------------------------------------------------------------------------------------

  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
  
  call MPI_Comm_rank(MPI_COMM_SUBGROUP, myRecComm%myid, ierr)
  call MPI_Comm_size(MPI_COMM_SUBGROUP, myRecComm%numnodes, ierr)  
  
  DO pID=1,subnodes
   DO pJD=1,myRecComm%NumHosts
    IF (adjustl(trim(myRecComm%all_hostnames(pID))).eq.adjustl(trim(myRecComm%unique_hostnames(pJD)))) then
     myRecComm%groupIDs(pID) = pJD
    END IF
   END DO
  END DO
  
  if (myid == 1) then
    ! Output the total number of unique host-nodes
    print *, "Total number of host-nodes: ", myRecComm%NumHosts
    print *, "Hostnames: "
    do i = 1, myRecComm%NumHosts
      print *, i,trim(myRecComm%unique_hostnames(i)),myRecComm%hostleaders(i)
    end do
    cFMT=' '
    write(cFMT,*) '(A,',subnodes-1,'(I0,(",")),I0)'
    write(*,cFMT) "GroupIDs: ",myRecComm%groupIDs
  end if

  CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
!   pause
  
 contains

  subroutine unique(input, output, leaders, group)
    character(len=256), allocatable, intent(in) :: input(:)
    character(len=256), allocatable, intent(out) :: output(:)
    integer, allocatable, intent(out) :: leaders(:),group(:)
    integer :: i, j, unique_count,nGroup
    logical :: is_unique
    character(len=256), allocatable :: temp(:)

    unique_count = 0

    do i = 1, size(input)
      is_unique = .true.
      do j = 1, unique_count
        if (trim(input(i)) == trim(output(j))) then
          is_unique = .false.
          exit
        end if
      end do
      if (is_unique) then
        ! Append to the temporary array
        if (unique_count == 0) then
          allocate(output(1))
          output(1) = input(i)
        else
          allocate(temp(unique_count))
          temp = output
          deallocate(output)
          allocate(output(unique_count + 1))
          output(1:unique_count) = temp
          output(unique_count + 1) = input(i)
          deallocate(temp)
        end if
        unique_count = unique_count + 1
      end if
    end do
    
    allocate(leaders(unique_count))
    
    leaders = size(input)
    
    do i = 1, size(input)
      do j = 1, unique_count
        if (trim(input(i)) == trim(output(j))) then
          if (leaders(j).gt.i) leaders(j) = i
        end if
      end do
    end do
    
    nGroup = 0
    do i = 1, size(input)
     if (trim(input(i)) == trim(input(myid))) then
      nGroup = nGroup + 1
     end if
    end do
    
    allocate(group(nGroup))
    nGroup = 0
    do i = 1, size(input)
     if (trim(input(i)) == trim(input(myid))) then
      nGroup = nGroup + 1
      group(nGroup) = i
     end if
    end do
    
  end subroutine unique

end subroutine FindNodes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
subroutine MemoryPrint(iDo,cGroup,cMSG)
 USE PP3D_MPI
 
 integer, intent(in) :: ido
 integer id,mem1,mem2
 integer, allocatable :: ddd(:),sendcounts1(:),sendcounts2(:)
 character*256 cFormat,filename1,filename2
 integer :: status
 character, intent(in):: cGroup*(*)
 character, intent(in):: cMSG*(*)
 
 CALL Get_PID(id,myid,mem1,mem2)
 
 allocate(sendcounts1(0:subnodes))
 allocate(sendcounts2(0:subnodes))
 
 if (cGroup.eq.'w') THEN ! with master
  call MPI_allgather(mem1, 1, MPI_INTEGER, sendcounts1, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
  call MPI_allgather(mem2, 1, MPI_INTEGER, sendcounts2, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
 END IF
 
 if (cGroup.eq.'s') THEN ! without master
  if (myid.ne.0) then
   call MPI_allgather(mem1, 1, MPI_INTEGER, sendcounts1(1:), 1, MPI_INTEGER, MPI_COMM_SUBS, ierr)
   sendcounts1(0) = 0
   call MPI_allgather(mem2, 1, MPI_INTEGER, sendcounts2(1:), 1, MPI_INTEGER, MPI_COMM_SUBS, ierr)
   sendcounts2(0) = 0
  end if
 END IF
 
 if (myid.eq.showid) then
 
  filename1 = '_data/memory_consumption_max.txt'
  filename2 = '_data/memory_consumption.txt'
  
  ! Open the file in append mode
  if (iDo.eq.0) THEN
    open(unit=4712, file=filename1)
    open(unit=4713, file=filename2)
  ELSE
  
   open(unit=4712, file=filename1, status='old', action='write', position='append', iostat=status)
  ! Check if the file opened successfully
   if (status /= 0) then
     open(unit=4712, file=filename1, status='new', action='write', iostat=status)
   endif
   
   open(unit=4713, file=filename2, status='old', action='write', position='append', iostat=status)
  ! Check if the file opened successfully
   if (status /= 0) then
     open(unit=4713, file=filename2, status='new', action='write', iostat=status)
   endif
  
  end if
 
!  OPEN (file=filename,unit=4712)
  write(cFormat,'(A,I0,A)') '(A,',subnodes,'(I0,(",")),I0)'
  
  write(4712,TRIM(cFormat)) "MemConsumptionIn_MB_["//trim(TRIM(cMSG))//"]:",sendcounts1/1024
  close(4712)
  
  write(4713,TRIM(cFormat)) "MemConsumptionIn_MB_["//trim(TRIM(cMSG))//"]:",sendcounts2/1024
  close(4713)
 END IF
 
 deallocate(sendcounts1)
 deallocate(sendcounts2)
 
END subroutine MemoryPrint
!
!----------------------------------------------
!
SUBROUTINE Get_PID(pid,id,mem1,mem2)

use iso_c_binding
implicit none
integer(C_INT) :: cpid
integer :: pid,id,mem1,mem2
character*256 cF

interface
    ! Declare the C function using Fortran-C interoperability
   function get_C_pid() bind(c, name="get_C_pid") result(pid)
       import :: C_INT
       integer(C_INT) :: pid
   end function get_C_pid
end interface

cpid = get_C_pid()

pid = int(cpid)


write(cF,'(A,I0,A)') "/proc/",pid,"/status"

CALL KeywordSearch(cF,mem1,mem2)

!     print *, "Rank", id, "PID:", pid, "Mem:",mem/1024, "MB"


END SUBROUTINE Get_PID


subroutine KeywordSearch(filename,mem1,mem2)

implicit none
character*256 :: line
character*7 :: keyword1 = "VmPeak:"
character*7 :: keyword2 = "VmSize:"
character*256 :: filename
integer :: iunit, i, pos,mem1,mem2
character cunit*8

mem1 = 0
mem2 = 0

! Open the input file
open(unit=iunit, file=filename, status='old', action='read', iostat=i)
if (i /= 0) then
    write(*,*) "Error opening the "//TRIM(filename)//" file"
    stop
end if

! Read and search line by line
do
    read(iunit, '(A)', iostat=i) line
    if (i /= 0) then
        exit
    end if

    ! Search for the keyword
    pos = index(line, keyword1)
    if (pos /= 0) then
        read(line(pos+len(keyword1):),*) mem1,cunit
        if (trim(cunit).eq."kB") THEN
!            write(*, *) "memory consumed: ", mem/1024, 'MB'
        ELSE
         mem1 = -1
!             write(*, *) "unusual unit ... "q
        END IF
    end if

    pos = index(line, keyword2)
    if (pos /= 0) then
        read(line(pos+len(keyword2):),*) mem2,cunit
        if (trim(cunit).eq."kB") THEN
!            write(*, *) "memory consumed: ", mem/1024, 'MB'
        ELSE
         mem2 = -1
!             write(*, *) "unusual unit ... "q
        END IF
    end if
    
end do

! Close the file
close(iunit)

end subroutine KeywordSearch
    