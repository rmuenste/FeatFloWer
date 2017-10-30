module sol_out
USE var_QuadScalar,ONLY:knvt,knet,knat,knel
!-------------------------------------------------------------------------------------------------
! A module for saving the solution values to 
! a file. This output(dump) is mainly done
! to a binary file.
!-------------------------------------------------------------------------------------------------

! a variable for counting the outputs
integer :: ifile = 0

contains
!
!-------------------------------------------------------------------------------------------------
! Wrapper routine for writing the solution to a file 
!-------------------------------------------------------------------------------------------------
! @param iInd number of the output
! @param istep number of the discrete time step
! @param simTime current simulation time
subroutine write_sol_to_file(imax_out, time_ns, output_idx)
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
use var_QuadScalar, only: myDump,istep_ns,myFBM,fieldPtr
USE Transport_Q1,ONLY:Tracer
USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier

implicit none
integer, intent(in) :: imax_out
real*8, intent(in) :: time_ns

integer, optional :: output_idx

! locals
integer :: iout
integer :: ndof
integer :: nelem
character(60) :: fieldName

type(fieldPtr), dimension(3) :: packed

if(.not.present(output_idx))then
  ifile = ifile+1
  iout=mod(ifile+imax_out-1,imax_out)+1
else
  iout = output_idx
end if

nelem = knel(nlmax)

ndof = knvt(NLMAX) + knat(NLMAX) + knet(NLMAX) + knel(NLMAX)

call write_vel_sol(iout,0,ndof,NLMIN,NLMAX,&
                   coarse%myELEMLINK,myDump%Vertices,&
                   QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

call write_pres_sol(iout,0,nelem,NLMIN,NLMAX,&
                    coarse%myELEMLINK,myDump%Elements,LinSc%ValP(NLMAX)%x)

call write_time_sol(iout,istep_ns, time_ns)

fieldName = "myvel"

packed(1)%p => QuadSc%ValU
packed(2)%p => QuadSc%ValV
packed(3)%p => QuadSc%ValW

call write_q2_sol(fieldName, iOut,0,ndof,NLMIN,NLMAX,coarse%myELEMLINK,myDump%Vertices,&
                  3, packed)

end subroutine write_sol_to_file
!
!-------------------------------------------------------------------------------------------------
! Wrapper routine for reading the solution from a file 
!-------------------------------------------------------------------------------------------------
! @param startFrom character string containing the start folder
! @param iLevel level adjustment for reading a solution 
! @param time_ns simulation time
subroutine read_sol_from_file(startFrom, iLevel, time_ns)

USE PP3D_MPI, ONLY:myid,coarse,myMPI_Barrier
USE def_FEAT
USE Transport_Q2P1,ONLY:QuadSc,LinSc,SetUp_myQ2Coor,bViscoElastic
USE var_QuadScalar,ONLY:myFBM,myDump,istep_ns
USE Transport_Q1,ONLY:Tracer

implicit none

character(60), intent(in) :: startFrom
integer, intent(in) :: iLevel
real*8, intent(inout) :: time_ns

integer :: nelem

nelem = knel(nlmax)

call read_vel_sol(startFrom,iLevel-1,nelem,NLMIN,NLMAX,&
                  coarse%myELEMLINK,myDump%Vertices,&
                  QuadSc%ValU,QuadSc%ValV,QuadSc%ValW)

call read_pres_sol(startFrom,iLevel-1,nelem,NLMIN,NLMAX,&
                   coarse%myELEMLINK,&
                   myDump%Elements,LinSc%ValP(NLMAX)%x)

call read_time_sol(startFrom, istep_ns, time_ns)

end subroutine read_sol_from_file
!
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
! Read the velocity solution from a single file
!-------------------------------------------------------------------------------------------------
! @param startFrom A string representation of the start directory 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param u the array of u-velocity component 
! @param v the array of v-velocity component 
! @param w the array of w-velocity component 
subroutine read_vel_sol_single(startFrom, iiLev,nn, nmin, nmax,elemmap,edofs, u, v, w)
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
real*8, dimension(:), intent(inout) :: u
real*8, dimension(:), intent(inout) :: v
real*8, dimension(:), intent(inout) :: w

! locals
integer :: iunit = 321
integer :: istatus
integer :: dofsInCoarseElement
integer :: elemCoarse
integer :: iel

integer :: ivt
integer :: jvt

integer :: global_idx
integer :: i_local
integer :: comp

real*8, dimension(:), allocatable :: buf

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

  iunit = iunit + myid

!  call read_sol_vel(startFrom, iiLev, comp, nn,& 
!                    elemCoarse, dofsInCoarseElement,&
!                    elemmap, edofs, u, v, w)

! if(myid.eq.1)then
!    DO i=1,elemCoarse
!    Write(*,*)i, " :  ", elemmap(i)
!    END DO
! end if

 i_local = 1
 open(unit=iunit, file="velocity.dmp", iostat=istatus, action="read")

 allocate(buf(dofsInCoarseElement)) 

 ! skip header
 READ(iunit,*) 

 DO iel=1,130

  READ(iunit,"(I3)") global_idx 

  if(elemmap(i_local) .eq. iel)then

    ! skip 3 lines
    read(iunit,*) buf(1:dofsInCoarseElement)

    do ivt=1,dofsInCoarseElement
     jvt = edofs(ivt,iel)
     !write(*,*)'jvt: ', jvt 

     if(myid.eq.3)then 
       u(1) = buf(ivt)
     end if

    end do

    read(iunit,*)
!    read(iunit,*) buf(1:dofsInCoarseElement)
!
!    do ivt=1,dofsInCoarseElement
!     jvt = edofs(iel,ivt)
!     v(jvt) = buf(ivt)
!    end do

    read(iunit,*)
!    read(iunit,*) buf(1:dofsInCoarseElement)
!
!    do ivt=1,dofsInCoarseElement
!     jvt = edofs(iel,ivt)
!     w(jvt) = buf(ivt)
!    end do

    if(myid.eq.1)write(*,*)iel, " : ", i_local, " : ", elemmap(i_local)
!    if(myid.eq.1)write(*,*)buf(1:dofsInCoarseElement)
    i_local = i_local + 1

  else
    READ(iunit,*) 
    READ(iunit,*) 
    READ(iunit,*) 
  end if

 END DO

 close(iunit)

end if

end subroutine read_vel_sol_single
!
!-------------------------------------------------------------------------------------------------
! Write a custom q2 field to file 
!-------------------------------------------------------------------------------------------------
! @param fieldName Name of the user-defined field 
! @param idx index of the output file 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param icomp Number of components of the output field 
! @param field_pack An array of structures that contain pointers to the 
!                   components of the output field
subroutine write_q2_sol(fieldName, idx, iiLev,nn, nmin, nmax,elemmap,edofs, icomp, field_pack)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
implicit none

character(60) :: fieldName
integer, intent(in) :: idx
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs

integer, intent(in) :: icomp

type(fieldPtr), dimension(:) :: field_pack

! locals
integer :: iunit = 321
integer :: i
integer :: dofsInCoarseElement
integer :: elemCoarse

elemCoarse = KNEL(nmin)

!packed(1)%p => u
!packed(2)%p => v
!packed(3)%p => w

! the subdivision level of an element on the 
! output level, i.e. lvl = 1, iiLev = 0
! (2**(1) + 1)[#dofs on an edge] * 3[#edges in y] * 3[#layers in z]
! = (2**(1)+1)**3 = 27
!
! Q2 dofs on a cube on level NLMAX+iiLev
dofsInCoarseElement = (2**((nmax+iiLev))+1)**3


if(myid.ne.0)then

  call clean_output_array(); 

  do i=1,icomp

    call wrap_pointer(idx, iiLev, nn,& 
                      nmin, nmax,&
                      elemmap, edofs, icomp, field_pack(i)%p)
  end do

  call write_sol_q2(TRIM(ADJUSTL(fieldName))//CHAR(0), idx, iiLev, icomp, nn,& 
                    elemCoarse, dofsInCoarseElement,&
                    elemmap, edofs)

end if

contains

subroutine wrap_pointer(idx, iiLev,nn, nmin, nmax,elemmap,edofs, icomp, p)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
implicit none

integer, intent(in) :: idx
integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs

integer, intent(in) :: icomp

real*8, dimension(:) :: p

integer :: dofsInCoarseElement
integer :: elemCoarse

elemCoarse = KNEL(nmin)

!packed(1)%p => u
!packed(2)%p => v
!packed(3)%p => w

! the subdivision level of an element on the 
! output level, i.e. lvl = 1, iiLev = 0
! (2**(1) + 1)[#dofs on an edge] * 3[#edges in y] * 3[#layers in z]
! = (2**(1)+1)**3 = 27
!
! Q2 dofs on a cube on level NLMAX+iiLev
dofsInCoarseElement = (2**((nmax+iiLev))+1)**3

 call add_output_array(p)

end subroutine wrap_pointer

end subroutine write_q2_sol
!
!-------------------------------------------------------------------------------------------------
! Read a custom q2 field from file 
!-------------------------------------------------------------------------------------------------
! @param fieldName Name of the user-defined field 
! @param idx index of the output file 
! @param iiLev the solution is written out on lvl: NLMAX+iiLev 
! @param nn the number of mesh elements on level NLMAX 
! @param nmin NLMIN 
! @param nmax NLMAX 
! @param elemmap a map from local to global element index 
! @param edofs a map from local to global element index 
! @param edofs an array of the fine level dofs in a coarse mesh element 
! @param icomp Number of components of the output field 
! @param field_pack An array of structures that contain pointers to the 
!                   components of the output field
subroutine read_q2_sol(fieldName, startFrom, iiLev,nn, nmin, nmax,elemmap,edofs, icomp, field_pack)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
implicit none

character(60) :: fieldName
character(60) :: startFrom

integer, intent(in) :: iiLev
integer, intent(in) :: nn
integer, intent(in) :: nmin
integer, intent(in) :: nmax

integer, dimension(:) :: elemmap
integer, dimension(:,:) :: edofs

integer, intent(in) :: icomp

type(fieldPtr), dimension(:) :: field_pack

! locals
integer :: iunit = 321
integer :: i
integer :: idx = 0
integer :: dofsInCoarseElement
integer :: elemCoarse

elemCoarse = KNEL(nmin)

!packed(1)%p => u
!packed(2)%p => v
!packed(3)%p => w

! the subdivision level of an element on the 
! output level, i.e. lvl = 1, iiLev = 0
! (2**(1) + 1)[#dofs on an edge] * 3[#edges in y] * 3[#layers in z]
! = (2**(1)+1)**3 = 27
!
! Q2 dofs on a cube on level NLMAX+iiLev
dofsInCoarseElement = (2**((nmax+iiLev))+1)**3


if(myid.ne.0)then

  call clean_output_array(); 

  do i=1,icomp

    call wrap_pointer(field_pack(i)%p)
  end do

  call read_sol_q2(TRIM(ADJUSTL(fieldName))//CHAR(0),TRIM(ADJUSTL(startFrom))//CHAR(0),&
                   idx, iiLev, icomp, nn,& 
                   elemCoarse, dofsInCoarseElement,&
                   elemmap, edofs)

end if

contains

subroutine wrap_pointer(p)
use pp3d_mpi, only:myid,coarse
use var_QuadScalar, only: fieldPtr
USE Transport_Q2P1,ONLY:QuadSc,LinSc,bViscoElastic
implicit none

real*8, dimension(:) :: p

 call add_output_array(p)

end subroutine wrap_pointer

end subroutine read_q2_sol
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
!
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
!
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
!
!-------------------------------------------------------------------------------------------------
! A general postprocessing for a Feat_FloWer application
!-------------------------------------------------------------------------------------------------
! @param dout Output interval
! @param iogmv Output index of the current file 
! @param istep number of the discrete time step
! @param inlU   
! @param inlT 
! @param filehandle Unit of the output file
!
subroutine postprocessing_app(dout, iogmv, inlU,inlT,filehandle)

include 'defs_include.h'
use var_QuadScalar, only: istep_ns
use post_utils, only: TimeStepCtrl

implicit none

integer, intent(in) :: filehandle

integer, intent(inout) :: iogmv
real, intent(inout) :: dout

INTEGER :: inlU,inlT,MFILE

! Output the solution in GMV or GiD format
IF (itns.eq.1) THEN
  CALL ZTIME(myStat%t0)
  CALL Output_Profiles(0)
  CALL ZTIME(myStat%t1)
  myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
END IF

IF(dout.LE.(timens+1e-10)) THEN

  iOGMV = istep_ns
  IF (itns.ne.1) THEN
    iogmv = iogmv - 1
    CALL ZTIME(myStat%t0)
    CALL Output_Profiles(iOGMV)
    CALL ZTIME(myStat%t1)
    myStat%tGMVOut = myStat%tGMVOut + (myStat%t1-myStat%t0)
  END IF
  dout=dout+dtgmv

  ! Save intermediate solution to a dump file
  IF (insav.NE.0.AND.itns.NE.1) THEN
    IF (MOD(iOGMV,insav).EQ.0) THEN
      CALL ZTIME(myStat%t0)
      call write_sol_to_file(insavn, timens)
      CALL ZTIME(myStat%t1)
      myStat%tDumpOut = myStat%tDumpOut + (myStat%t1-myStat%t0)
    END IF
  END IF

END IF

! Timestep control
CALL TimeStepCtrl(tstep,inlU,inlT,filehandle)

! Interaction from user
CALL ProcessControl(filehandle,mterm)

end subroutine postprocessing_app

end module sol_out

