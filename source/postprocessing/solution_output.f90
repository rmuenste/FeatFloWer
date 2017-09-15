module sol_out
!-------------------------------------------------------------------------------------------------
! A module for saving the solution values to 
! a file. This output(dump) is mainly done
! to a binary file.
!-------------------------------------------------------------------------------------------------
use var_QuadScalar, only: knel
contains

!------------------------------------------------------------------------------------------------- 
! Write the pressure solution to file
!-------------------------------------------------------------------------------------------------
! write_pres_sol: The structure of the pressure solution array is:
! pres(1:4*nn)
! the entries pres(4*(iel-1)+1) to pres(4*(iel-1)+4) contain
! the mean value and the dx,dy,dz derivatives in the element iel.
!
! The file format of the dump is:
! Header until the '\n' character
! After the header:
! Global coarse element idx
! #dofsInCoarseElement entries of the mean pressure values 
! #dofsInCoarseElement entries of the d/dx values 
! #dofsInCoarseElement entries of the d/dy values 
! #dofsInCoarseElement entries of the d/dz values 
! ...
!  
! 
! @param iInd 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 
subroutine write_pres_sol(iInd,iiLev,nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

integer, intent(in) :: iInd
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt

elemCoarse = KNEL(nmin)

! DO i=1,nn
!  Field2(i) = Field1(4*(i-1)+1)
!  Field3(i) = Field1(4*(i-1)+2)
!  Field4(i) = Field1(4*(i-1)+3)
!  Field5(i) = Field1(4*(i-1)+4)
! END DO

! the subdivision level of an element on the 
! output level, NLMAX = 2, iiLev = 0
! 8**(2-1) = 8
! meaning on level 2 a coarse grid
! element is divided into 8 elements
dofsInCoarseElement = 8**((nmax+iiLev)-1)

if(myid.ne.0)then

  call write_sol_pres(iInd, iiLev, nn ,elemCoarse, dofsInCoarseElement, elemmap, edofs, pres);

end if

end subroutine write_pres_sol
!------------------------------------------------------------------------------------------------- 
! Read the pressure solution from a file
!-------------------------------------------------------------------------------------------------
! read_pres_sol: The structure of the pressure solution array is:
! pres(1:4*nn)
! the entries pres(4*(iel-1)+1) to pres(4*(iel-1)+4) contain
! the mean value and the dx,dy,dz derivatives in the element iel.
! @param startFrom 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 
subroutine read_pres_sol(startFrom,iiLev,nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

character(60), intent(in) :: startFrom
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt

elemCoarse = KNEL(nmin)

! the subdivision level of an element on the 
! output level, NLMAX = 2, iiLev = 0
! 8**(2-1) = 8
! meaning on level 2 a coarse grid
! element is divided into 8 elements
dofsInCoarseElement = 8**((nmax+iiLev)-1)

if(myid.ne.0)then

  call read_sol_pres(startFrom, iiLev, nn ,elemCoarse, dofsInCoarseElement, elemmap, edofs, pres);

end if

end subroutine read_pres_sol
!
!-------------------------------------------------------------------------------------------------
! Write the velocity solution to file
!-------------------------------------------------------------------------------------------------
! @param iInd 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param u the array of u-velocity component 
! @param v the array of v-velocity component 
! @param w the array of w-velocity component 
subroutine write_vel_sol(iInd,iiLev,nn, nmin, nmax,elemmap,edofs, u, v, w)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
implicit none

integer, intent(in) :: iInd
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:), target :: u
real*8, dimension(:), target :: v
real*8, dimension(:), target :: w

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt
integer :: comp
type(fieldPtr), dimension(3) :: packed

comp = 3

elemCoarse = KNEL(nmin)

packed(1)%p => u
packed(2)%p => v
packed(3)%p => w

! the subdivision level of an element on the 
! output level, i.e. lvl = 1, iiLev = 0
! (2**(1) + 1)[#dofs on an edge] * 3[#edges in y] * 3[#layers in z]
! = (2**(1)+1)**3 = 27
!
! Q2 dofs on a cube on level NLMAX+iiLev
dofsInCoarseElement = (2**((nmax+iiLev))+1)**3

if(myid.ne.0)then

  call write_sol_vel(iInd, iiLev, comp, nn,& 
                     elemCoarse, dofsInCoarseElement,&
                     elemmap, edofs, u, v, w)

end if

end subroutine write_vel_sol
!
!-------------------------------------------------------------------------------------------------
! Read the velocity solution from file
!-------------------------------------------------------------------------------------------------
! @param startFrom A string representation of the start directory 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param u the array of u-velocity component 
! @param v the array of v-velocity component 
! @param w the array of w-velocity component 
subroutine read_vel_sol(startFrom, iiLev,nn, nmin, nmax,elemmap,edofs, u, v, w)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
implicit none

character(60), intent(in) :: startFrom
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:), target :: u
real*8, dimension(:), target :: v
real*8, dimension(:), target :: w

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: ivt
integer :: jvt
integer :: comp

comp = 3

elemCoarse = KNEL(nmin)

! the subdivision level of an element on the 
! output level, i.e. lvl = 1, iiLev = 0
! (2**(1) + 1)[#dofs on an edge] * 3[#edges in y] * 3[#layers in z]
! = (2**(1)+1)**3 = 27
!
! Q2 dofs on a cube on level NLMAX+iiLev
dofsInCoarseElement = (2**((nmax+iiLev))+1)**3

if(myid.ne.0)then

  call read_sol_vel(startFrom, iiLev, comp, nn,& 
                    elemCoarse, dofsInCoarseElement,&
                    elemmap, edofs, u, v, w)

end if

end subroutine read_vel_sol
!
!-------------------------------------------------------------------------------------------------
! Unit test function for P1 dump output 
!-------------------------------------------------------------------------------------------------
! @param fileName Name of the output file
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param pres the array of pressure values on lvl NLMAX 

subroutine write_pres_test(fileName, nn, nmin, nmax,elemmap,edofs,pres)
use pp3d_mpi, only:myid,coarse
implicit none

CHARACTER*(60) :: fileName
CHARACTER*(20) :: idName
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: pres

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: i
integer :: jvt
integer :: iiLev = 0

elemCoarse = KNEL(nmin)

if(myid.ne.0)then
  write(idName,'(I0)')myid
  open(unit=iunit, file=trim(adjustl(fileName))//trim(adjustl(idName)), iostat=istatus, action="write")

 do i=1,elemCoarse
  write(iunit,*)pres(4*(i-1)+1)
  write(iunit,*)pres(4*(i-1)+2)
  write(iunit,*)pres(4*(i-1)+3)
  write(iunit,*)pres(4*(i-1)+4)
 end do

 close(iunit)

end if

end subroutine write_pres_test
!
!-------------------------------------------------------------------------------------------------
! Unit test function for Q2 dump output 
!-------------------------------------------------------------------------------------------------
! @param fileName Name of the output file
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param u u-component velocity 
! @param v v-component velocity 
! @param w w-component velocity 
subroutine write_vel_test(fileName, nn, nmin, nmax,elemmap,edofs,u, v, w)
use pp3d_mpi, only:myid,coarse
implicit none

character(60) :: fileName
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs
real*8, dimension(:) :: u
real*8, dimension(:) :: v
real*8, dimension(:) :: w

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel
integer :: i
integer :: jvt
integer :: iiLev = 0
character(20) :: idName

elemCoarse = KNEL(nmin)

if(myid.ne.0)then
  write(idName,'(I0)')myid
  open(unit=iunit, file=trim(adjustl(fileName))//trim(adjustl(idName)), iostat=istatus, action="write")

 do i=1,elemCoarse
  write(iunit,*)u(i)
  write(iunit,*)v(i)
  write(iunit,*)w(i)
 end do

 close(iunit)

end if

end subroutine write_vel_test
!
!-------------------------------------------------------------------------------------------------
! Write the output time to file
!-------------------------------------------------------------------------------------------------
! @param iInd number of the output
! @param istep number of the discrete time step
! @param simTime current simulation time
subroutine write_time_sol(iInd, istep, simTime)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
implicit none

integer, intent(in) :: iInd
integer, intent(in) :: istep
real*8 :: simTime

! locals
integer :: iunit = 321

if(myid.ne.0)then

  call write_sol_time(iInd, istep, simTime)

end if

end subroutine write_time_sol
!
!-------------------------------------------------------------------------------------------------
! Read the time from file
!-------------------------------------------------------------------------------------------------
! @param iInd number of the output
! @param istep number of the discrete time step
! @param simTime current simulation time
subroutine read_time_sol(startFrom, istep, simTime)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
implicit none

character(60), intent(in) :: startFrom

integer :: istep
real*8  :: simTime

! locals
integer :: iunit = 321

if(myid.ne.0)then

  call read_sol_time(startFrom, istep, simTime)
  istep = istep + 1
end if

end subroutine read_time_sol

end module sol_out

