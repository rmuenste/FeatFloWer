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


subroutine closeFile(iFile)
#ifdef __INTEL_COMPILER 
  USE IFCORE
#endif  

#ifdef __GFORTRAN__
  interface
    function fsync (fd) bind(c,name="fsync")
    use iso_c_binding, only: c_int
      integer(c_int), value :: fd
      integer(c_int) :: fsync
    end function fsync
  end interface
  integer :: iRet,iNum
#endif  
  integer, intent(in) :: ifile
  
#ifdef __GFORTRAN__
  flush(iFile)
  iNum = fnum(iFile)
  iRet = fsync(iNum)
  ! Handle possible error
  if (iRet.ne.0) WRITE (*,*) 'Closing(GFORTRAN) of data file has failed ... '
  close(iFile)
#elif __INTEL_COMPILER 
  IF (.NOT. COMMITQQ(iFile)) WRITE (*,*) 'Closing(INTEL) of data file has failed ... '
  close(iFile)
#elif
  close(iFile)
#endif  

end subroutine closeFile
