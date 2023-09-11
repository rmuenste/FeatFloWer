! CALL OUTPUT_BurgersMatrix_slave(mumps_par%N,mumps_par%NZ,mumps_par%IRN,mumps_par%JCN,mumps_par%A_loc)
! CALL OUTPUT_BurgersMatrix_master(mumps_par%NZ_loc,mumps_par%IRN_loc,mumps_par%JCN_loc,mumps_par%RHS)

SUBROUTINE OUTPUT_BurgersMatrix_master(N,NZ,IRN,JCN,RHS)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes
implicit none
integer N,NZ
integer :: IRN(*),JCN(*)
REAL*8 :: RHS(*)
!!!!
character :: cF*(256)
integer  :: myFile=5489
integer i,istat

RETURN

write(cF,'(A,I2.2,A,I2.2,A)')  'input_',subnodes,'/input_',MYID,'.txt'
open(file=adjustl(trim(cF)),unit=myFile, action="write",iostat=istat)

WRITE(myFile,'(I0)') N
WRITE(myFile,'(I0)') NZ

DO i=1,NZ
 WRITE(myFile,'(I0,A,I0)') IRN(i),' ',JCN(i)
END DO

DO i=1,N
 WRITE(myFile,'(ES12.4)') RHS(i)
END DO

close(myFile)

write(*,*) 'here we go'
pause

END

SUBROUTINE OUTPUT_BurgersMatrix_slave(NZ_loc,IRN_loc,JCN_loc,A_loc)
USE PP3D_MPI, ONLY : myid,master,showid,subnodes
implicit none
integer NZ_loc
integer :: IRN_loc(*),JCN_loc(*)
REAL*8 :: A_loc(*)
!!!!
character :: cF*(256)
integer  :: myFile=5489
integer i,istat

RETURN

write(cF,'(A,I2.2,A,I2.2,A)')  'input_',subnodes,'/input_',MYID,'.txt'
open(file=adjustl(trim(cF)),unit=myFile, action="write",iostat=istat)

WRITE(myFile,'(I0)') NZ_loc

DO i=1,NZ_loc
 WRITE(myFile,'(I0,A,I0,A,ES12.4)') IRN_loc(i),' ',JCN_loc(i),' ',A_loc(i)
END DO


close(myFile)

write(*,*) 'here we go'
pause
END 