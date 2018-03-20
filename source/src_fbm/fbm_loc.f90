!##############################################################################
!# ****************************************************************************
!# <name> Point location module loc</name>
!# ****************************************************************************
!# loc_LocatePoints
!# loc_QueryPointElement_noeval
!# loc_QueryPointElement_eval
!# loc_QueryRemotePoint
!# loc_QueryRemotePoint_func
!# loc_QueryPointElement_func
!# loc_RemoteLocation_noeval
!# loc_RemoteLocation_eval
!#
!##############################################################################
MODULE loc 

use var_QuadScalar

include 'mpif.h'

! 
integer, parameter, public :: LOC_NOTFOUND = -1
! 
integer, parameter, public :: LOC_DUMMY    = -2
! 
integer, parameter, public :: LOC_FOUND    = 1

type myPoint
  Real*8, dimension(3) :: xyz
  Real*8, dimension(3) :: res
  integer :: status
end type

type myPointGD
  Real*8, dimension(3) :: xyz
  Real*8, dimension(5) :: res
  integer :: status
end type

interface loc_QueryPointElement
  module procedure loc_QueryPointElement_noeval
  module procedure loc_QueryPointElement_eval
  module procedure loc_QueryRemotePoint    
end interface

interface loc_RemoteLocation
  module procedure loc_RemoteLocation_eval
  module procedure loc_RemoteLocation_noeval
  module procedure loc_RemoteLocationGD
end interface

interface loc_LocatePoints
  module procedure loc_LocatePointsVel
  module procedure loc_LocatePointsGD
end interface

contains
!
!****************************************************************************  
!
subroutine loc_LocatePointsVel(dcorvg,kvert,kedge,karea,dpoints,npoints,dres,u,v,w)
!**@Description
  ! Wrapper function for point location in different meshes
  ! dpoints: coordinates that should be located and evaluated
  ! dres: stores result values that are evaluated in the input coordinates
!**@EndDescription  
use pp3d_mpi
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_DomainCandidates,fbmaux_evalE013_sing
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoints(3,*)
Real*8  dres(3,*)
Real*8, dimension(:) ::u,v,w
Integer kvert(8,*),kedge(12,*),karea(6,*)
integer :: npoints

! local variables
logical, dimension(:), allocatable :: found
integer, dimension(:), allocatable :: index_arr
integer :: i,ires,imissed,offset
Real*8  ref(3)

! bookkeeping array for found points
allocate(found(npoints))

found=.false.
imissed=0

write(*,*)'Number of points: ',npoints,myid

! search for points in domain
do i=1,npoints

  ref(:)=0d0
  ires=loc_QueryPointElement(dcorvg,kvert,kedge,karea,dpoints(:,i),ref)
      
  if(ires.ne.LOC_NOTFOUND)then
    found(i)=.true.
    ! evaluate the field in the element
    call fbmaux_evalE013_sing(dcorvg,kvert,kedge,karea,u,v,w,ref,dres(:,i),ires)      
  else
    found(i)=.false.  
    imissed=imissed+1
  end if
end do

allocate(index_arr(imissed))

offset=1
do i=1,npoints
  if(.not.found(i))then
    index_arr(offset)=i
    offset=offset+1
  end if
end do

! free memory
deallocate(found)

write(*,*)'Number of points not found in domain: ',imissed,myid
  
call loc_RemoteLocation(dcorvg,kvert,kedge,karea,dpoints,dres,index_arr,imissed,npoints,u,v,w)

end subroutine loc_LocatePointsVel
!
!****************************************************************************  
!
integer function loc_QueryPointElement_noeval(dcorvg,kvert,kedge,karea,dpoint,ref) result(iel)
!**@Description
  ! Wrapper function for the PointInElement query
!**@EndDescription  
use pp3d_mpi, only: myid
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use fbmaux, only: fbmaux_TestDomainBdry
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoint(3)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
real*8 ref(3)

! local variables
real*8 vertices(3,8)
integer :: ibody,i,ibodyc,nnel
real*8 point(3),dvalues(3)
integer, dimension(:), allocatable :: elements


! query the domain bounding box
if(.not.fbmaux_TestDomainBdry(dpoint,myid))then
  iel = LOC_NOTFOUND
  return
end if

! query the uniform grid
call ug_pointquery(dpoint,nnel)

! the body is not in the domain -> skip
if(nnel .eq. 0)then
  iel = LOC_NOTFOUND
  return
end if

! allocate memory for element array
if(allocated(elements))deallocate(elements)
allocate(elements(nnel))

! get the element ids
call ug_getelements(elements)


! start element determination
iel = loc_QueryPointElement_func(dcorvg,kvert,kedge,&
                                                   karea,elements,&
                                                   nnel,dpoint,ref) 
                                                                                                                                 
end function loc_QueryPointElement_noeval
!
!****************************************************************************  
!
integer function loc_QueryPointElement_eval(dcorvg,kvert,kedge,karea,dpoint,u,v,w) result(iel)
!**@Description
  ! Wrapper function for the direct PointInElement query
!**@EndDescription  
use pp3d_mpi, only: myid
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_evalE013_sing
implicit none
REAL*8  dcorvg(3,*)
REAL*8  u(*),v(*),w(*)
Real*8  dpoint(3)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)

! local variables
real*8 vertices(3,8)
integer :: ibody,i,ibodyc,nnel
real*8 point(3),dvalues(3)
integer, dimension(:), allocatable :: elements
Real*8  ref(3)

! query the domain bounding box
if(.not.fbmaux_TestDomainBdry(dpoint,myid))then
  iel = LOC_NOTFOUND
  return
end if

! query the uniform grid
call ug_pointquery(dpoint,nnel)

! the body is not in the domain -> skip
if(nnel .eq. 0)then
  iel = LOC_NOTFOUND
  return
end if

! allocate memory for element array
if(allocated(elements))deallocate(elements)
allocate(elements(nnel))

! get the element ids
call ug_getelements(elements)


! start element determination
iel = loc_QueryPointElement_func(dcorvg,kvert,kedge,&
                                                   karea,elements,&
                                                   nnel,dpoint,ref) 
         
! directly evaluate at the reference coordinate in the element
call fbmaux_evalE013_sing(dcorvg,kvert,kedge,karea,u,v,w,ref,dvalues,iel)               
                                    
end function loc_QueryPointElement_eval
!
!****************************************************************************  
!
integer function loc_QueryRemotePoint(dcorvg,kvert,kedge,karea,dpoint,ref) result(iel)
!**@Description
  ! Wrapper function for the PointInElement query
!**@EndDescription  
use pp3d_mpi, only: myid
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use fbmaux, only: fbmaux_TestDomainBdry
implicit none
REAL*8  dcorvg(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
type(myPoint) :: dpoint
real*8 ref(3)

! local variables
real*8 vertices(3,8)
integer :: nnel
real*8 point(3)
integer, dimension(:), allocatable :: elements

if(.not.fbmaux_TestDomainBdry(dpoint%xyz,myid))then
  iel=LOC_NOTFOUND
  dpoint%status=LOC_NOTFOUND
  return
end if

! query the uniform grid
call ug_pointquery(dpoint%xyz,nnel)

! the body is not in the domain -> skip
if(nnel .eq. 0)then
  dpoint%status=LOC_NOTFOUND
  iel=LOC_NOTFOUND
  return
end if

! allocate memory for element array
if(allocated(elements))deallocate(elements)
allocate(elements(nnel))

! get the element ids
call ug_getelements(elements)

! start element determination
iel = loc_QueryRemotePoint_func(dcorvg,kvert,kedge,&
                                                 karea,elements,&
                                                 nnel,dpoint,ref) 
                                    
! set the found marker
if(iel.ne.LOC_NOTFOUND)then
  dpoint%status=LOC_FOUND
end if
                   
end function loc_QueryRemotePoint
!
!****************************************************************************  
!
integer function loc_QueryRemotePoint_func(dcorvg,kvert,kedge,karea,elements,nnel,dpoint,ref) 
!**@Description
  ! The function tests a number of elements
  ! if they contain the query point by using a transformation
  ! to the reference element
!**@EndDescription     
use def_FEAT
use var_QuadScalar,only:myFBM 
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use fbmaux, only: fbmaux_IsInBoundingBox,fbmaux_PointInHex
implicit none
      
REAL*8  dcorvg(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ibody,nnel
integer, dimension(*) :: elements
type(myPoint) :: dpoint      
Real*8, dimension(3) :: ref

! local variables
logical :: found
integer :: i,ielem,iconv


found = .false.

loc_QueryRemotePoint_func = LOC_NOTFOUND

do i=1,nnel

  ref(:) = 0d0

  ielem = elements(i)
 
  if(.not.fbmaux_IsInBoundingBox(dcorvg,kvert,ielem,dpoint%xyz))cycle

  found = fbmaux_PointInHex(dpoint%xyz(1),dpoint%xyz(2),dpoint%xyz(3),&
                            ref(1),ref(2),ref(3),&
                            ielem,iconv,dcorvg,kvert)
 
  if(found .eqv. .true.)then
    loc_QueryRemotePoint_func = ielem
    return
  end if

end do

end function loc_QueryRemotePoint_func
!
!****************************************************************************  
!
integer function loc_QueryPointElement_func(dcorvg,kvert,kedge,karea,elements,nnel,dpoint,ref) 
!**@Description
  ! The function tests a number of elements
  ! if they contain the query point by using a transformation
  ! to the reference element
!**@EndDescription     
use def_FEAT
use var_QuadScalar,only:myFBM 
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use fbmaux, only: fbmaux_IsInBoundingBox,fbmaux_PointInHex
implicit none
      
REAL*8  dcorvg(3,*)
Real*8  dpoint(3)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ibody,nnel
integer, dimension(*) :: elements
Real*8, dimension(3) :: ref

! local variables
logical :: found
integer :: i,ielem,iconv


found = .false.

loc_QueryPointElement_func = LOC_NOTFOUND

do i=1,nnel

  ref(:) = 0d0

  ielem = elements(i)
 
  if(.not.fbmaux_IsInBoundingBox(dcorvg,kvert,ielem,dpoint))cycle

  ref(:)=0d0
  found = fbmaux_PointInHex(dpoint(1),dpoint(2),dpoint(3),&
                            ref(1),ref(2),ref(3),&
                            ielem,iconv,dcorvg,kvert)
 
  if(found .eqv. .true.)then
    loc_QueryPointElement_func = ielem
    return
  end if

end do

end function loc_QueryPointElement_func
!
!****************************************************************************  
!
subroutine loc_RemoteLocation_noeval(dcorvg,kvert,kedge,karea,dpoints,npoints)
!**@Description
  ! Wrapper function for the distributed point location
  ! the function only locates the element that contain the point
!**@EndDescription  
!use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use pp3d_mpi
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_DomainCandidates
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoints(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: npoints,ierror
REAL*8  ref(3)

! local variables
integer, dimension(:,:), allocatable :: icandidates
integer, dimension(:), allocatable :: inumcand,ijobs,iremotejobs,ierrors
integer :: i,j,ielementsInGroup,imax,ireduce,irank,itotal
integer, dimension(:), allocatable :: igroup,iallmax
type(myPoint), dimension(:), allocatable :: allpoints
type(myPoint), dimension(:), allocatable :: recpoints

integer :: pointtype
integer, dimension(2) :: oldtypes(0:1),blockcounts(0:1),offsets(0:1)
integer :: extent,ielement


! make a copy of dcorvg
allocate(icandidates(subnodes,npoints))
allocate(inumcand(npoints))
allocate(ijobs(subnodes))
allocate(iremotejobs(subnodes))
allocate(igroup(subnodes))
allocate(iallmax(subnodes))
allocate(ierrors(npoints))

ireduce=0

call MPI_Reduce(npoints,ireduce,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_SUBS,IERROR)

if(myid.eq.1)then
  write(*,*)'Maximum of points: ',ireduce
end if

call MPI_BCast(ireduce,1,MPI_INTEGER,0,MPI_COMM_SUBS,ierror)

! first block
offsets(0)=0
oldtypes(0)=MPI_REAL8
blockcounts(0)=6

! next block
call MPI_TYPE_EXTENT(MPI_REAL8,extent,ierror)
offsets(1)=6*extent
oldtypes(1)=MPI_INTEGER
blockcounts(1)=1

call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,pointtype,ierror)
call MPI_TYPE_COMMIT(pointtype,ierror)

allocate(allpoints(subnodes*ireduce))
allocate(recpoints(subnodes*ireduce))

do i=1,subnodes
 do j=1,ireduce
   if(j.le.npoints)then
   allpoints((i-1)*ireduce+j)%xyz(:)=dpoints(:,j)
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=-1
   else
   allpoints((i-1)*ireduce+j)%xyz(:)=0d0
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=-2
   end if
 end do
end do

! communicate array to all processes
CALL MPI_ALLTOALL(allpoints,ireduce,pointtype,recpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   if(recpoints((i-1)*ireduce+j)%status.eq.-2)cycle
   ielement = loc_QueryRemotePoint(dcorvg,kvert,kedge,karea,recpoints((i-1)*ireduce+j),ref)   
   if(recpoints((i-1)*ireduce+j)%status.eq.1)then
!     write(*,'(A,I8,A,I8)')'found point from: ',i,' in domain: ',myid
   end if
 end do
end do

! communicate the result
CALL MPI_ALLTOALL(recpoints,ireduce,pointtype,allpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

ierrors=-1
itotal=0

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   if(allpoints((i-1)*ireduce+j)%status.eq.1)then
     if(ierrors(j).lt.0)then
       ierrors(j)=0
     end if
!     write(*,'(A,I8,A,I8)')'Process: ',myid,' my point was found in domain: ',i
   end if
 end do
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    if(fbmaux_TestDomainBdry(dpoints(:,j),myid))then
      ierrors(j)=0
    end if
  end if
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    itotal=itotal+1
    !write(*,*)'Missed point:  ',dpoints(:,j)
  end if
end do
write(*,'(A,I8,A,I8)')'Process: ',myid,' total number of missed points: ',itotal

end subroutine loc_RemoteLocation_noeval
!
!****************************************************************************  
!
subroutine loc_RemoteLocation_eval(dcorvg,kvert,kedge,karea,dpoints,dres,index_arr,imissed,ipoints,u,v,w)
!**@Description
  ! Wrapper function for the distributed point location
  ! that directly evaluate a function in the point and store the result
!**@EndDescription  
!use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use pp3d_mpi
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_DomainCandidates,fbmaux_evalE013_sing
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoints(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ipoints,ierror,imissed
integer, dimension(:) :: index_arr
Real*8, dimension(:) ::u,v,w
Real*8  dres(3,*)
Real*8  ref(3)

! local variables
integer, dimension(:,:), allocatable :: icandidates
integer, dimension(:), allocatable :: inumcand,ijobs,iremotejobs,ierrors
integer :: i,j,ielementsInGroup,imax,ireduce,irank,itotal
integer, dimension(:), allocatable :: igroup,iallmax
type(myPoint), dimension(:), allocatable :: allpoints
type(myPoint), dimension(:), allocatable :: recpoints

integer :: pointtype
integer, dimension(2) :: oldtypes(0:1),blockcounts(0:1),offsets(0:1)
integer :: extent,npoints,ires
real*8, dimension(3) :: center,vec
real*8 :: length

npoints=imissed

! make a copy of dcorvg
allocate(icandidates(subnodes,npoints))
allocate(inumcand(npoints))
allocate(ijobs(subnodes))
allocate(iremotejobs(subnodes))
allocate(igroup(subnodes))
allocate(iallmax(subnodes))
allocate(ierrors(npoints))

! do i=1,npoints
!   call fbmaux_DomainCandidates(dpoints(:,i),icandidates(:,i),inumcand(i))
! !  write(*,*)'Number of candidates: ',inumcand(i),myid
! end do

! ijobs=0

! do i=1,npoints
!   do j=1,inumcand(i)
!     ! if there is a job increase the job number by 1
!     ijobs(icandidates(j,i)) = ijobs(icandidates(j,i)) + 1
!   end do  
! end do

! imax=0

! do j=1,subnodes
!   if(ijobs(j).gt.imax)then
!     imax=ijobs(j)
!   end if
! !  write(*,'(A,3I8)')'Number of jobs for process/myid: ',ijobs(j),j,myid
! end do

ireduce=0

call MPI_Reduce(npoints,ireduce,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_SUBS,IERROR)

! if(myid.eq.1)then
!   write(*,*)'Maximum of points: ',ireduce
! end if

call MPI_BCast(ireduce,1,MPI_INTEGER,0,MPI_COMM_SUBS,ierror)

! make point array of maximum size

! first block
offsets(0)=0
oldtypes(0)=MPI_REAL8
blockcounts(0)=6

! next block
call MPI_TYPE_EXTENT(MPI_REAL8,extent,ierror)
offsets(1)=6*extent
oldtypes(1)=MPI_INTEGER
blockcounts(1)=1

call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,pointtype,ierror)
call MPI_TYPE_COMMIT(pointtype,ierror)

allocate(allpoints(subnodes*ireduce))
allocate(recpoints(subnodes*ireduce))

do i=1,subnodes
 do j=1,ireduce
   if(j.le.npoints)then
   allpoints((i-1)*ireduce+j)%xyz(:)=dpoints(:,index_arr(j))
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=LOC_NOTFOUND
   else
   allpoints((i-1)*ireduce+j)%xyz(:)=0d0
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=LOC_DUMMY
   end if
 end do
end do

! communicate array to all processes
CALL MPI_ALLTOALL(allpoints,ireduce,pointtype,recpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   ref(:)=0d0
   if(recpoints((i-1)*ireduce+j)%status.eq.LOC_DUMMY)cycle
   ires = loc_QueryPointElement(dcorvg,kvert,kedge,karea,recpoints((i-1)*ireduce+j),ref)  

   if(ires.ne.LOC_NOTFOUND)then
     recpoints((i-1)*ireduce+j)%status=LOC_FOUND

     ! evaluate the field in the element
     call fbmaux_evalE013_sing(dcorvg,kvert,kedge,karea,u,v,w,&
                          ref,&
                          recpoints((i-1)*ireduce+j)%res,ires)         
                        
   end if
 end do
end do

! communicate the result
CALL MPI_ALLTOALL(recpoints,ireduce,pointtype,allpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

ierrors=-1
itotal=0

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   if(allpoints((i-1)*ireduce+j)%status.eq.LOC_FOUND)then
     
     ! store the result of the point
     dres(:,index_arr(j))=allpoints((i-1)*ireduce+j)%res(:)
   
     if(ierrors(j).lt.0)then
       ierrors(j)=0
     end if
   end if
 end do
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    center(1)=0.5
    center(2)=0.2
    center(3)=dpoints(3,index_arr(j))  
    vec(:)=dpoints(:,index_arr(j))  -center(:)
    length=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
    
    if(fbmaux_TestDomainBdry(dpoints(:,j),myid))then
      ierrors(j)=0
    elseif(length.le.0.05d0)then
      ierrors(j)=0    
    end if
 
  end if
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    itotal=itotal+1
  end if
end do
write(*,'(A,I8,A,I8)')'Process: ',myid,' number of missed points after distributed search: ',itotal

end subroutine loc_RemoteLocation_eval
!
!****************************************************************************  
!
subroutine loc_LocatePointsGD(dcorvg,dpoints,npoints,kvert,kedge,karea,dResult,Ux,Uy,Uz)
!**@Description
  ! Wrapper function for point location in different meshes
  ! dpoints: coordinates that should be located and evaluated
  ! dres: stores result values that are evaluated in the input coordinates
!**@EndDescription  
use pp3d_mpi
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_DomainCandidates,fbmaux_evalPhi_sing
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoints(3,*)
Real*8, dimension(:) ::Ux,Uy,Uz
Integer kvert(8,*),kedge(12,*),karea(6,*)
integer :: npoints
type(tGDefFunc), dimension(:) :: dResult

! local variables
logical, dimension(:), allocatable :: found
integer, dimension(:), allocatable :: index_arr
integer :: i,ires,imissed,offset
Real*8  ref(3)

! bookkeeping array for found points
allocate(found(npoints))

found=.false.
imissed=0

!write(*,*)'Number of points: ',npoints,myid

! search for points in domain
do i=1,npoints

  ref(:)=0d0
  ires=loc_QueryPointElement(dcorvg,kvert,kedge,karea,dpoints(:,i),ref)
      
  if(ires.ne.LOC_NOTFOUND)then
    found(i)=.true.
    ! evaluate the field in the element
    call fbmaux_evalPhi_sing(dcorvg,kvert,kedge,karea,Ux,Uy,Uz,ref,dResult(i)%dvalues,ires)      
  else
    found(i)=.false.  
    imissed=imissed+1
  end if
end do

allocate(index_arr(imissed))

offset=1
do i=1,npoints
  if(.not.found(i))then
    index_arr(offset)=i
    offset=offset+1
  end if
end do

! free memory
deallocate(found)

! if(imissed .gt. 0)then
!   write(*,*)'Number of points not found in domain: ',imissed,myid
! end if
  
call loc_RemoteLocation(dcorvg,kvert,kedge,karea,dpoints,dResult,index_arr,imissed,npoints,Ux,Uy,Uz)

end subroutine loc_LocatePointsGD
!
!****************************************************************************  
!
subroutine loc_RemoteLocationGD(dcorvg,kvert,kedge,karea,dpoints,dres,index_arr,imissed,ipoints,Ux,Uy,Uz)
!**@Description
  ! Wrapper function for distributed point location
  ! that directly evaluates the function values for the GD problem in the points
!**@EndDescription  
!use pp3d_mpi, only: myid,mgBoundingBox,subnodes
use pp3d_mpi
use fbmaux, only: fbmaux_TestDomainBdry,fbmaux_DomainCandidates,fbmaux_evalPhi_sing
implicit none
REAL*8  dcorvg(3,*)
Real*8  dpoints(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ipoints,ierror,imissed
integer, dimension(:) :: index_arr
Real*8, dimension(:) ::Ux,Uy,Uz
type(tGDefFunc), dimension(:) :: dres
Real*8  ref(3)

! local variables
integer, dimension(:,:), allocatable :: icandidates
integer, dimension(:), allocatable :: inumcand,ijobs,iremotejobs,ierrors
integer :: i,j,ielementsInGroup,imax,ireduce,irank,itotal
integer, dimension(:), allocatable :: igroup,iallmax
type(myPointGD), dimension(:), allocatable :: allpoints
type(myPointGD), dimension(:), allocatable :: recpoints
type(myPoint) :: dpoint

integer :: pointtype
integer, dimension(2) :: oldtypes(0:1),blockcounts(0:1),offsets(0:1)
integer :: extent,npoints,ires
real*8, dimension(3) :: center,vec
real*8 :: length

npoints=imissed

! make a copy of dcorvg
allocate(icandidates(subnodes,npoints))
allocate(inumcand(npoints))
allocate(ijobs(subnodes))
allocate(iremotejobs(subnodes))
allocate(igroup(subnodes))
allocate(iallmax(subnodes))
allocate(ierrors(npoints))

ireduce=0

! calculate the maximum number of points
call MPI_Reduce(npoints,ireduce,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_SUBS,IERROR)

! communicate the maximum number of points
call MPI_BCast(ireduce,1,MPI_INTEGER,0,MPI_COMM_SUBS,ierror)

! first block
offsets(0)=0
oldtypes(0)=MPI_REAL8
blockcounts(0)=8

! next block
call MPI_TYPE_EXTENT(MPI_REAL8,extent,ierror)
offsets(1)=8*extent
oldtypes(1)=MPI_INTEGER
blockcounts(1)=1

call MPI_TYPE_STRUCT(2,blockcounts,offsets,oldtypes,pointtype,ierror)
call MPI_TYPE_COMMIT(pointtype,ierror)

allocate(allpoints(subnodes*ireduce))
allocate(recpoints(subnodes*ireduce))

do i=1,subnodes
 do j=1,ireduce
   if(j.le.npoints)then
   allpoints((i-1)*ireduce+j)%xyz(:)=dpoints(:,index_arr(j))
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=LOC_NOTFOUND
   else
   allpoints((i-1)*ireduce+j)%xyz(:)=0d0
   allpoints((i-1)*ireduce+j)%res(:)=0d0
   allpoints((i-1)*ireduce+j)%status=LOC_DUMMY
   end if
 end do
end do

! communicate array to all processes
CALL MPI_ALLTOALL(allpoints,ireduce,pointtype,recpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   ref(:)=0d0
   if(recpoints((i-1)*ireduce+j)%status.eq.LOC_DUMMY)cycle
   dpoint%xyz(:)= recpoints((i-1)*ireduce+j)%xyz(:)
   ires = loc_QueryPointElement(dcorvg,kvert,kedge,karea,dpoint,ref)  

   if(ires.ne.LOC_NOTFOUND)then
     recpoints((i-1)*ireduce+j)%status=LOC_FOUND

     ! evaluate the field in the element
     call fbmaux_evalPhi_sing(dcorvg,kvert,kedge,karea,Ux,Uy,Uz,ref,&
                          recpoints((i-1)*ireduce+j)%res,&
                          ires)         
                        
   end if
 end do
end do

! communicate the result
CALL MPI_ALLTOALL(recpoints,ireduce,pointtype,allpoints,ireduce,pointtype,MPI_COMM_SUBS,ierror)

ierrors=-1
itotal=0

do i=1,subnodes
 if(i.eq.myid)cycle
 do j=1,ireduce
   if(allpoints((i-1)*ireduce+j)%status.eq.LOC_FOUND)then
     
     ! store the result of the point    
     dres(index_arr(j))%dvalues=allpoints((i-1)*ireduce+j)%res(:)
   
     if(ierrors(j).lt.0)then
       ierrors(j)=0
     end if
   end if
 end do
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    center(1)=0.5
    center(2)=0.2
    center(3)=dpoints(3,index_arr(j))  
    vec(:)=dpoints(:,index_arr(j))  -center(:)
    length=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
    
    if(fbmaux_TestDomainBdry(dpoints(:,j),myid))then
      ierrors(j)=0
    elseif(length.le.0.05d0)then
      ierrors(j)=0    
    end if
 
  end if
end do

do j=1,npoints
  if(ierrors(j).ne.0)then
    itotal=itotal+1
  end if
end do

!write(*,'(A,I8,A,I8)')'Process: ',myid,' number of missed points after distributed search: ',itotal

end subroutine loc_RemoteLocationGD
!
!****************************************************************************  
!
end module
