SUBROUTINE CheckIfFolderIsThereCreateIfNot(cFolder,iBARR)
use iniparser
USE PP3D_MPI, ONLY:myid,MPI_COMM_SUBS,MPI_COMM_WORLD

IMPLICIT NONE
INTEGER IBARR
CHARACTER cFolder*(256)
INTEGER IERR

if (myid.eq.1) then
 if (.not. inip_isDirectory(cFolder)) then
   write(*,*) 'Creating Directory: '//adjustl(trim(cFolder))
   call inip_makeDirectory(cFolder)
 end if
end if

if (iBarr.eq.1) then! without master
 CALL MPI_BARRIER(MPI_COMM_SUBS,IERR)
else            ! with master
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
end if

END SUBROUTINE CheckIfFolderIsThereCreateIfNot
