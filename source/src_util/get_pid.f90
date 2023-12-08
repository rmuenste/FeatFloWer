subroutine MemoryPrint(iDo,cGroup,cMSG)
 USE PP3D_MPI
 
 integer, intent(in) :: ido
 integer id,mem
 integer, allocatable :: ddd(:),sendcounts(:)
 character*256 cFormat,filename
 integer :: status
 character, intent(in):: cGroup*(*)
 character, intent(in):: cMSG*(*)
 
 CALL Get_PID(id,myid,mem)
 
 allocate(sendcounts(0:subnodes))
 
 if (cGroup.eq.'w') THEN ! with master
  call MPI_allgather(mem, 1, MPI_INTEGER, sendcounts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
 END IF
 
 if (cGroup.eq.'s') THEN ! without master
  if (myid.ne.0) then
   call MPI_allgather(mem, 1, MPI_INTEGER, sendcounts(1:), 1, MPI_INTEGER, MPI_COMM_SUBS, ierr)
   sendcounts(0) = 0
  end if
 END IF
 
 if (myid.eq.showid) then
 
  filename = '_data/memory_consumption.txt'
  
  ! Open the file in append mode
  if (iDo.eq.0) THEN
    open(unit=4712, file=filename)
  ELSE
  
   open(unit=4712, file=filename, status='old', action='write', position='append', iostat=status)
  
  ! Check if the file opened successfully
   if (status /= 0) then
     open(unit=4712, file=filename, status='new', action='write', iostat=status)
   endif
  end if
 
!  OPEN (file=filename,unit=4712)
  write(cFormat,'(A,I0,A)') '(A,',subnodes,'(I0,(",")),I0)'
  write(4712,TRIM(cFormat)) "MemConsumptionIn_MB_["//trim(TRIM(cMSG))//"]:",sendcounts/1024
  close(4712)
 END IF
 deallocate(sendcounts)
 
END subroutine MemoryPrint
!
!----------------------------------------------
!
SUBROUTINE Get_PID(pid,id,mem)

use iso_c_binding
implicit none
integer(C_INT) :: cpid
integer :: pid,id,mem
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

CALL KeywordSearch(cF,mem)

!     print *, "Rank", id, "PID:", pid, "Mem:",mem/1024, "MB"


END SUBROUTINE Get_PID


subroutine KeywordSearch(filename,mem)

implicit none
character*256 :: line
character*7 :: keyword = "VmPeak:"
character*256 :: filename
integer :: iunit, i, pos,mem
character cunit*8

mem = 0

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
    pos = index(line, keyword)
    if (pos /= 0) then
        read(line(pos+len(keyword):),*) mem,cunit
        if (trim(cunit).eq."kB") THEN
!            write(*, *) "memory consumed: ", mem/1024, 'MB'
        ELSE
         mem = -1
!             write(*, *) "unusual unit ... "q
        END IF
    end if
end do

! Close the file
close(iunit)

end subroutine KeywordSearch
    