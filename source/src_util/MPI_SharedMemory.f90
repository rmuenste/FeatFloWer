module shared_memory_module
    use mpi
    use, intrinsic :: ISO_C_BINDING, only: c_ptr, c_f_pointer
    implicit none
    
    TYPE mg_dV
     REAL*8  , DIMENSION(:)  , pointer :: x
    END TYPE mg_dV

    TYPE mg_kV
     INTEGER  , DIMENSION(:) , pointer :: x
    END TYPE mg_kV

    TYPE tRecComm
     logical :: bPrepared = .false.
     TYPE(mg_dV), DIMENSION(:), pointer :: s,r
     TYPE(mg_kV), DIMENSION(:), pointer :: CODECs,CODECr
     integer, allocatable :: CODECs_win(:),CODECr_win(:)
     integer :: SIZEs_win,SIZEr_win
     type(c_ptr), allocatable :: CODECsptr(:),CODECrptr(:)
     type(c_ptr) :: SIZEsptr,SIZErptr
     integer, allocatable :: StartOfAllRecords(:,:)
     integer , pointer  :: sendSIZE(:),recvSIZE(:)
     
     integer, allocatable :: winS(:),winR(:)
     type(c_ptr), allocatable :: ptrS(:),ptrR(:)
    END TYPE tRecComm
    TYPE (tRecComm), allocatable :: myRC(:),myPC(:)
    
    contains

    subroutine get_shared_memory_INT(comm, shared_mem, shared_size, win, baseptr, ierr)
        integer, intent(in) :: comm
        integer, pointer :: shared_mem(:)
        integer, intent(in) :: shared_size
        integer, intent(out) :: win
        type(c_ptr), intent(out) :: baseptr
        integer, intent(out) :: ierr
        integer(kind=MPI_ADDRESS_KIND) :: win_size
        integer :: disp_unit
        INTEGER :: info
        
        ! Set the size and displacement unit
        win_size = shared_size * 4  ! Size in bytes (4 bytes per INT)
        disp_unit = 4               ! Size of each element (INT)

        
        call MPI_Info_create(info, ierr)
        call MPI_Info_set(info, "no_locks", "true", ierr)
        
        ! Allocate shared memory window
!         call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, comm, baseptr, win, ierr)
        call MPI_Win_allocate_shared(win_size, disp_unit, info, comm, baseptr, win, ierr)
        if (ierr /= MPI_SUCCESS) then
            print *, "Error in MPI_Win_allocate_shared"
            return
        end if

        ! Query the base pointer to the shared memory
        call MPI_Win_shared_query(win, MPI_PROC_NULL, win_size, disp_unit, baseptr, ierr)
        if (ierr /= MPI_SUCCESS) then
            print *, "Error in MPI_Win_shared_query"
            return
        end if

        ! Convert the base pointer to a Fortran pointer
        call c_f_pointer(baseptr, shared_mem, [shared_size])
    end subroutine get_shared_memory_INT

    subroutine get_shared_memory_DBL(comm, shared_mem, shared_size, win, baseptr, ierr)
        integer, intent(in) :: comm
        real*8, pointer :: shared_mem(:)
        integer, intent(in) :: shared_size
        integer, intent(out) :: win
        type(c_ptr), intent(out) :: baseptr
        integer, intent(out) :: ierr
        integer(kind=MPI_ADDRESS_KIND) :: win_size
        integer :: disp_unit
        INTEGER :: info

        ! Set the size and displacement unit
        win_size = shared_size * 8  ! Size in bytes (4 bytes per DOUBLE)
        disp_unit = 8               ! Size of each element (DOUBLE)

        call MPI_Info_create(info, ierr)
        call MPI_Info_set(info, "no_locks", "true", ierr)
        
        ! Allocate shared memory window
!        call MPI_Win_allocate_shared(win_size, disp_unit, MPI_INFO_NULL, comm, baseptr, win, ierr)
        call MPI_Win_allocate_shared(win_size, disp_unit, info, comm, baseptr, win, ierr)
        if (ierr /= MPI_SUCCESS) then
            print *, "Error in MPI_Win_allocate_shared"
            return
        end if

        ! Query the base pointer to the shared memory
        call MPI_Win_shared_query(win, MPI_PROC_NULL, win_size, disp_unit, baseptr, ierr)
        if (ierr /= MPI_SUCCESS) then
            print *, "Error in MPI_Win_shared_query"
            return
        end if

        ! Convert the base pointer to a Fortran pointer
        call c_f_pointer(baseptr, shared_mem, [shared_size])
    end subroutine get_shared_memory_DBL
    
end module shared_memory_module
