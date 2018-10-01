!##############################################################################
!# ****************************************************************************
!# <name> fbm </name>
!# ****************************************************************************
!# In this module we define functions that related to the handling of
!# fictitious boundary components.
!#
!# 
!##############################################################################
MODULE fbm 

use var_QuadScalar

integer, dimension(:), allocatable :: mykvel

! 
integer, parameter, public :: SRCH_UNDEF = 0
! 
integer, parameter, public :: SRCH_OK = 1
! 
integer, parameter, public :: SRCH_FAIL = 2
! 
integer, parameter, public :: SRCH_EPIC_FAIL = 3
! 
integer, parameter, public :: SRCH_LEFTDOMAIN = 4
! 
integer, parameter, public :: SRCH_WRONGELEMENT = 5
! 
integer, parameter, public :: SRCH_MAXITER = 6
! 
integer, parameter, public :: SRCH_RES7 = 7

contains
!=========================================================================
! 
!=========================================================================
Subroutine fbm_updateFBM(DensityL,dTime,simTime,Gravity,mfile,myid,u,v,w,p,usr_updateFBM)
use PP3D_MPI, only: myMPI_Barrier
use cinterface

include 'fbm_up_include.h'

procedure(update_fbm_handler) :: usr_updateFBM 

! Liquid density, time step and simulation time
real*8, intent(in) :: DensityL,dTime,simTime

! The gravity vector
real*8, intent(inout) :: Gravity(3)

! The process id
integer, intent(in) :: myid

! The protocol file
integer, intent(in) :: mfile

! U/V/W velocity components
REAL*8, dimension(:), intent(inout) :: u,v,w

! Pressure, viscosity
Real*8, dimension(:), intent(inout) :: p

external E013

! local variables
integer :: ilevel 

if (calculateDynamics()) then

 if(myid.eq.1) write(*,*)'> Dynamics update'

 ilevel=mg_mesh%nlmax
 CALL SETLEV(2)
 CALL GetForces(Properties%ForceScale,u,&
                v,w,p,&
                FictKNPR,Viscosity,&
                mg_mesh%level(ilevel)%kvert,&
                mg_mesh%level(ilevel)%karea,&
                mg_mesh%level(ilevel)%kedge,&
                mg_mesh%level(ilevel)%dcorvg,&
                E013)

 call usr_updateFBM(DensityL,dTime,simTime,Gravity,mfile,myid)

end if

end subroutine fbm_updateFBM
!=========================================================================
! 
!=========================================================================
subroutine fbm_updateFBMGeom(px, py, pz, bndryId, fictId, dist, usr_geomFBM)
use PP3D_MPI, only: myMPI_Barrier
use cinterface

! Coordinates of the query point 
real*8, intent(in) :: px, py, pz 

! Id of the boundary component
integer, intent(inout) :: bndryId

! fictId
integer, intent(inout) :: fictId

! Distance solution in the query point 
real*8, intent(inout) :: dist 

include 'fbm_geom_include.h'
procedure(fbm_geom_handler) :: usr_geomFBM 

! Begin function

 call usr_geomFBM(px, py, pz, bndryId, fictId, dist)

end subroutine fbm_updateFBMGeom
!=========================================================================
! 
!=========================================================================
subroutine fbm_getSoftKnpr(x,y,z, bndryId, inpr, dist)
! 
!   This subroutine handles the FBM geometric computations
!   for a single(!) elastic or soft body.
!
use var_QuadScalar, only : myFBM
implicit none

! Coordinates of the query point 
real*8, intent(in) :: x, y, z 

! Id of the boundary component
integer, intent(inout) :: bndryId

! fictId
integer, intent(inout) :: inpr

! Distance solution in the query point 
real*8, intent(inout) :: dist 

! local variables
integer :: IP,ipc,isin

  inpr = 0
  dist = 1000.0d0

  ! We are dealing with a single soft body
  ! so set IP to 1
  IP = 1

  ipc=ip-1
  isin = 0

  ! Here we call the interface function to
  ! the geometry library to handle the FBM
  ! computation for a soft body.
  call insidesoftbody(x,y,z,ipc,isin)

  if(isin .gt. 0)then
  inpr = isin
  end if

end subroutine fbm_getSoftKnpr
!=========================================================================
! 
!=========================================================================
subroutine fbm_updateSoftBodyDynamics(DensityL,dTime,simTime,Gravity,mfile,myid)
use var_QuadScalar, only : myFBM
use PP3D_MPI, only: myMPI_Barrier
use cinterface

real*8, intent(in) :: DensityL ! fluid density
real*8, intent(in) :: dTime    ! time step
real*8, intent(in) :: simTime  ! simulation time
real*8, dimension(3), intent(in) :: Gravity  ! simulation time

integer, intent(in) :: mfile   ! prot file unit
integer, intent(in) :: myid    ! process id

integer :: IP,ipc

real*8 :: volume,mass,massR,radius,dimomir,dSubStep
real*8,parameter :: PI = 3.1415926535897931D0
real*8 :: RForce(3),dVelocity(3),dOmega(3),timecoll
integer :: iSubSteps

 iSubSteps = 1

 call settimestep(dTime)

 if(myid.eq.1) write(*,*)'> Dynamics update'

 ! After the rigid body solver has computed a 
 ! step, we have to get the new particle state
 ! values from the rigid body solver
 call FBM_GetParticleStateUpdate()

 ! Communicate new force + torque
 IP = 1
 ipc = ip-1
 ! Communicate the force
 call setforce(myFBM%particleOld(IP)%ResistanceForce(1),&
               myFBM%particleOld(IP)%ResistanceForce(2),&
               myFBM%particleOld(IP)%ResistanceForce(3),ipc)

 ! Communicate the torque
 call settorque(myFBM%particleOld(IP)%TorqueForce(1),&
 myFBM%particleOld(IP)%TorqueForce(2),&
 myFBM%particleOld(IP)%TorqueForce(3),ipc)
     
 ! in the first loop we updated the external force and
 ! torque, with this we can start the collision handling
 call settimestep(dTime)      

 ! update velocities by the force determined in the time step
 ! This function will transfer the current hydrodynamic forces
 ! to the rigid body solver, so it can add those values in
 ! the next rigid body solver step
 call velocityupdate_soft()     

 call settimestep(dTime/real(iSubSteps))
 call starttiming()      

 ! Invoke a rigid body solver step 
 !     do iStep=1,iSubSteps
 !       call startcollisionpipeline()
 !     end do

 ! Integrate the equations of motion for the soft body
 call stepsoftbody(myFBM%particleOld(IP)%ResistanceForce(1),&
                   myFBM%particleOld(IP)%ResistanceForce(1),&
                   myFBM%particleOld(IP)%ResistanceForce(1),&
                   dTime)

 call gettiming(timecoll)

 ! After the rigid body solver has computed a 
 ! step, we have to get the new particle state
 ! values from the rigid body solver
 call FBM_GetParticleStateUpdate()

 !-------------------------------------------------------------------------------------
 ! Backup the particle parameters so the old values are available in the next timestep
 !-------------------------------------------------------------------------------------
 DO IP = 1,myFBM%nParticles
 ipc = ip-1

 ! Backup the forces
 myFBM%particleOld(IP)%ResistanceForce = &
 myFBM%particleNew(IP)%ResistanceForce
    
 ! Backup the torque
 myFBM%particleOld(IP)%TorqueForce = &
 myFBM%particleNew(IP)%TorqueForce

 ! Backup the positions
 myFBM%particleOld(IP)%Position = &
 myFBM%particleNew(IP)%Position

 ! Backup the velocities
 myFBM%particleOld(IP)%Velocity = &
 myFBM%particleNew(IP)%Velocity

 ! Backup the angles
 myFBM%particleOld(IP)%Angle = &
 myFBM%particleNew(IP)%Angle

 ! Backup the angular velocity
 myFBM%particleOld(IP)%AngularVelocity = &
 myFBM%particleNew(IP)%AngularVelocity

 IF((myid.eq.1) .and. (IP .lt. 2)) THEN
   if(myid.eq.1)then
   write(*,*)'After collision handling'
   write(mfile,'(A19,3ES14.4)') "ResistanceForce: ",&
       myFBM%particleNew(IP)%ResistanceForce
   write(*,'(A19,3ES14.4)') "ResistanceForce: ",&
       myFBM%particleNew(IP)%ResistanceForce
   write(mfile,'(A19,3ES14.4)') "TorqueForce: ",&
       myFBM%particleNew(IP)%TorqueForce        
   write(*,'(A19,3ES14.4)') "TorqueForce: ",&
       myFBM%particleNew(IP)%TorqueForce
   write(mfile,'(A19,4ES14.4)') "Position: ",&
       myFBM%particleNew(IP)%Position,simTime        
   write(*,'(A19,4ES14.4)') "Position: ",&
       myFBM%particleNew(IP)%Position,simTime
   write(mfile,'(A19,4ES14.4)') "PartVel: ",&
       myFBM%particleNew(IP)%Velocity,simTime        
   write(*,'(A19,4ES14.4)') "PartVel: ",&
       myFBM%particleNew(IP)%Velocity,simTime
   write(mfile,'(A19,3ES14.4)') "Angle: ",&
       myFBM%particleNew(IP)%Angle        
   write(*,'(A19,3ES14.4)') "Angle: ",&
       myFBM%particleNew(IP)%Angle
   write(mfile,'(A19,3ES14.4)') "AngularVelocity: ",&
       myFBM%particleNew(IP)%AngularVelocity        
   write(*,'(A19,3ES14.4)') "AngularVelocity: ",&
       myFBM%particleNew(IP)%AngularVelocity

   end if
 end if
 end do ! all particles

end subroutine fbm_updateSoftBodyDynamics
!=========================================================================
! 
!=========================================================================
Subroutine fbm_velBCUpdate(x,y,z,valu,valv,valw,ip,t,usr_velBCUpdate)
use PP3D_MPI, only: myMPI_Barrier
use cinterface

include 'fbm_vel_bc_include.h'

procedure(fbm_velBC_handler) :: usr_velBCUpdate 

! The FBM value of the boundary vertex
integer, intent(in) :: ip

! The coordinates of the boundary vertex
! and the simulation time
real*8 , intent(in) :: x,y,z,t

! The velocitiy values of the boundary vertex
real*8 , intent(inout) :: valu,valv,valw

call usr_velBCUpdate(x,y,z,valu,valv,valw,ip,t)

end subroutine fbm_velBCUpdate
!=========================================================================
! 
!=========================================================================
subroutine fbm_velBC(x,y,z,valu,valv,valw,ip,t)
use var_QuadScalar, only : myFBM,bRefFrame
implicit none

! Parameters
integer, intent(in) :: ip
real*8 , intent(in) :: x,y,z,t
real*8 , intent(inout) :: valu,valv,valw

! local variables
REAL*8 :: dvelz_x,dvelz_y,dvely_z,dvely_x,dvelx_y,dvelx_z
real*8 :: Velo(3),Pos(3),Omega(3)
integer :: ipc

valu = 0d0
valv = 0d0
valw = 0d0

if(ip .le. myFBM%nParticles)then
  Velo  = myFBM%particleNEW(IP)%Velocity
  Pos   = myFBM%particleNEW(IP)%Position
  Omega = myFBM%particleNEW(IP)%AngularVelocity

  dvelz_x = -(y-Pos(2))*Omega(3)
  dvelz_y = +(x-Pos(1))*Omega(3)

  dvely_z = -(x-Pos(1))*Omega(2)
  dvely_x = +(z-Pos(3))*Omega(2)

  dvelx_y = -(z-Pos(3))*Omega(1)
  dvelx_z = +(y-Pos(2))*Omega(1)

  valu = velo(1) + dvelz_x + dvely_x
  valv = velo(2) + dvelz_y + dvelx_y
  valw = velo(3) + dvelx_z + dvely_z
end if

end subroutine fbm_velBC
!=========================================================================
! 
!=========================================================================
SUBROUTINE GetFictKnpr_DIE(X,Y,Z,iBndr,inpr,Dist)
! 
!   This subroutine handles the FBM geometric computations
!   for a single(!) elastic or soft body.
!
use var_QuadScalar, only : dCGALtoRealFactor,activeFBM_Z_Position
implicit none

! Coordinates of the query point 
real*8, intent(in) :: x, y, z 

! Id of the boundary component
integer, intent(inout) :: iBndr

! fictId
integer, intent(inout) :: inpr

! Distance solution in the query point 
real*8, intent(inout) :: dist 

! local variables
integer :: IP,ipc,isin, idynType
real*8 :: dist_sign, cpx,cpy,cpz, d_temp

 inpr = 0
 dist_sign = 1
 Dist = 1d8
  
 DO IP = 1,myFBM%nParticles
  
  ipc=ip-1
  isin = 0
  call get_dynamics_type(ipc, idynType) 
  dist_sign = +1d0
  if (idynType.eq.2) dist_sign = -1d0
!  if (IP.eq.1) dist_sign = -1d0

  call isinelementid(dCGALtoRealFactor*x,dCGALtoRealFactor*y,dCGALtoRealFactor*z,ipc,isin)
  if(isin .gt. 0)then
   dist_sign = -1d0*dist_sign
   call getclosestpointid(dCGALtoRealFactor*x,dCGALtoRealFactor*y,dCGALtoRealFactor*z,cpx,cpy,cpz,d_temp,ipc);        
  else
   dist_sign = +1d0*dist_sign
   call getclosestpointid(dCGALtoRealFactor*x,dCGALtoRealFactor*y,dCGALtoRealFactor*z,cpx,cpy,cpz,d_temp,ipc);        
  end if
  
  dist = min(dist,dist_sign * d_temp)
  
 end do
 
 if (dist.lt.0d0) THEN
  inpr = 100
  if (z.lt.activeFBM_Z_Position/dCGALtoRealFactor) THEN
   inpr = 0
  end if
 end if

END SUBROUTINE GetFictKnpr_DIE
!=========================================================================
! 
!=========================================================================
subroutine fbm_getFictKnpr(x,y,z,bndryId,inpr,dist)
! 
!   This subroutine handles the FBM geometric computations
!
use var_QuadScalar, only : myFBM
implicit none

! Coordinates of the query point 
real*8, intent(in) :: x, y, z 

! Id of the boundary component
integer, intent(inout) :: bndryId

! fictId
integer, intent(inout) :: inpr

! Distance solution in the query point 
real*8, intent(inout) :: dist 

! local variables
integer :: IP,ipc,isin

 inpr = 0
 dist = 1000.0d0

 DO IP = 1,myFBM%nParticles
  ipc=ip-1
  isin = 0
  call isinelementid(x,y,z,ipc,isin)
  if(isin .gt. 0)then
   inpr = isin+1
   inpr = IP
  end if
 end do

end subroutine fbm_getFictKnpr

!=========================================================================
! 
!=========================================================================
!
!****************************************************************************  
!
!subroutine fbm_DomainCandidates(dpoint,icandidates)
!  !**@Description
!  !    Computes the domains that a body can potentially be located in
!  !**End@Description
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!
!  implicit none
!  real*8, dimension(3)   :: dpoint
!  integer, dimension(:)  :: icandidates
!
!  ! local variables
!  integer :: i,inumcand
!
!  inumcand=1
!
!  do i=1,subnodes
!  if(i.eq.myid)cycle
!
!  ! is body intersects domain
!  if(fbm_IsInBndryBox(dpoint,i))then
!    icandidates(inumcand)=i
!    inumcand = inumcand + 1
!  end if
!
!  end do
!
!end subroutine fbm_DomainCandidates
!!
!!****************************************************************************  
!!
!logical function fbm_IsInBndryBox(dpoint,pid)
!  !**@Description
!
!  !**End@Description
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!
!  implicit none
!  real*8, dimension(3)  :: dpoint
!  integer :: pid
!
!  ! local variables
!  integer :: i
!
!  if(dpoint(1) .lt. mgBoundingBox(pid)%vertices(1,1))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  if(dpoint(2) .lt. mgBoundingBox(pid)%vertices(2,1))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  if(dpoint(3) .lt. mgBoundingBox(pid)%vertices(3,1))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  if(dpoint(1) .gt. mgBoundingBox(pid)%vertices(1,2))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  if(dpoint(2) .gt. mgBoundingBox(pid)%vertices(2,2))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  if(dpoint(3) .gt. mgBoundingBox(pid)%vertices(3,2))then
!    fbm_IsInBndryBox = .false.
!    return
!  end if
!
!  fbm_IsInBndryBox = .true.
!
!end function fbm_IsInBndryBox
!!
!! ----------------------------------------------
!!
!SUBROUTINE fbm_FictPre(dcorvg,kvert,kedge,karea)
!  USE var_QuadScalar,ONLY:myExport,myFBM 
!  REAL*8  dcorvg(3,*)
!  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  REAL*8 PX,PY,PZ,DIST
!  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
!  INTEGER NeighE(2,12),NeighA(4,6)
!  real*8 nvertices(3,8)
!  real*8 point(3)
!  integer :: iel,idofs
!  integer :: iae,ineighbour
!  integer :: ivert,iedge,iarea,ielem,ibodyc,in,myelement,ibody
!  logical :: found
!  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!
!  do ibody=1,myFBM%nParticles
!
!  myelement = 0
!  ibodyc = ibody - 1
!  found  = .false.
!  call setelement(myelement,ibodyc)
!
!  do i=1,NEL
!  nvertices(:,1) = dcorvg(:,kvert(1,i))
!  nvertices(:,2) = dcorvg(:,kvert(2,i))
!  nvertices(:,3) = dcorvg(:,kvert(3,i))
!  nvertices(:,4) = dcorvg(:,kvert(4,i))
!  nvertices(:,5) = dcorvg(:,kvert(5,i))
!  nvertices(:,6) = dcorvg(:,kvert(6,i))
!  nvertices(:,7) = dcorvg(:,kvert(7,i))
!  nvertices(:,8) = dcorvg(:,kvert(8,i))
!
!  ! calculate the dofs inside
!  do ivert=1,8
!  in=0
!  call isinelementid(nvertices(1,ivert),nvertices(2,ivert),nvertices(3,ivert),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!    exit
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found)exit
!
!  ! calculate the edge dofs inside
!  do iedge=1,12
!  in=0
!  ivt1 = kvert(NeighE(1,iedge),i)
!  ivt2 = kvert(NeighE(2,iedge),i)
!  point(:)= 0.5 * (dcorvg(:,ivt1) + dcorvg(:,ivt2))
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found)exit
!
!  ! calculate the face dofs inside
!  do iarea=1,6
!  in=0
!  ivt1 = kvert(NeighE(1,iarea),i)
!  ivt2 = kvert(NeighE(2,iarea),i)
!  ivt3 = kvert(NeighA(3,iarea),i)
!  ivt4 = kvert(NeighA(4,iarea),i)
!  point(:)= 0.25 * (dcorvg(:,ivt1) + dcorvg(:,ivt2) + dcorvg(:,ivt3) + dcorvg(:,ivt4))
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found)exit
!
!  ! calculate the element dof
!  in=0
!  point(:)= 0.125 * (nvertices(:,1) + nvertices(:,2) + nvertices(:,3) + nvertices(:,4) + &
!    nvertices(:,5) + nvertices(:,6) + nvertices(:,7) + nvertices(:,8))                   
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!
!  end do ! end for all elements
!
!  end do ! end for all particles
!
!END SUBROUTINE fbm_FictPre
!!
!! ----------------------------------------------
!!
!! subroutine fbm_FictPre2(dcorvg,kvert,kedge,karea)
!! USE var_QuadScalar,ONLY:myExport,myFBM 
!! use pproc
!! use def_feat
!! implicit none
!! REAL*8  dcorvg(3,*)
!! INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!! REAL*8 PX,PY,PZ,DIST
!! INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,ibody,ibodyc,in
!! INTEGER NeighE(2,12),NeighA(4,6)
!! integer :: myelement
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! do ibody=1,myFBM%nParticles
!! 
!!   call fbm_FictPreBodyAll(dcorvg,kvert,kedge,karea,ibody)
!! 
!! end do
!! 
!! end subroutine fbm_FictPre2
!!
!! ----------------------------------------------
!!
!! subroutine fbm_FictPre2Body(dcorvg,kvert,kedge,karea,ibody)
!! USE var_QuadScalar,ONLY:myExport,myFBM 
!! use pp3d_mpi, only: myid
!! use pproc
!! use def_feat
!! implicit none
!! REAL*8  dcorvg(3,*)
!! INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!! REAL*8 PX,PY,PZ,DIST
!! integer :: ibody
!! INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,ibodyc,in
!! INTEGER NeighE(2,12),NeighA(4,6)
!! integer :: myelement
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! myelement = 0
!! ibodyc = ibody - 1
!! in = 0
!! call setelement(myelement,ibodyc)
!! 
!! do i=1,nvt
!!    PX = dcorvg(1,I)
!!    PY = dcorvg(2,I)
!!    PZ = dcorvg(3,I)
!!    call isinelementid(px,py,pz,ibodyc,in)
!!    if(in .gt. 0)then
!!       ! we found an element set the data and go on to the next body
!!       call setelement(mykvel(i),ibodyc)
!!       return
!!    end if
!! end do
!! 
!! k=1
!! do i=1,nel
!!    do j=1,12
!!       if (k.eq.kedge(j,i)) then
!!          ivt1 = kvert(NeighE(1,j),i)
!!          ivt2 = kvert(NeighE(2,j),i)
!!          PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
!!          PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
!!          PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
!!          call isinelementid(px,py,pz,ibodyc,in)
!!          if(in .gt. 0)then
!!             ! we found an element set the data and go on to the next body
!!             ! if(myid.eq.5)write(*,*)'incremental: setting edge dof ',k+nvt,i
!!             call setelement(i,ibodyc)
!!             return
!!          end if
!!          k = k + 1
!!       end if
!!    end do
!! end do
!! 
!! k=1
!! do i=1,nel
!!    do j=1,6
!!       if (k.eq.karea(j,i)) then
!!          ivt1 = kvert(NeighA(1,j),i)
!!          ivt2 = kvert(NeighA(2,j),i)
!!          ivt3 = kvert(NeighA(3,j),i)
!!          ivt4 = kvert(NeighA(4,j),i)
!!          PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
!!          PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
!!          PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
!!          call isinelementid(px,py,pz,ibodyc,in)
!!          if(in .gt. 0)then
!!             ! we found an element set the data and go on to the next body
!!             ! if(myid.eq.5)write(*,*)'incremental: setting area dof ',nvt+net+k,i
!!             call setelement(i,ibodyc)
!!             return
!!          end if
!!          k = k + 1
!!       end if
!!    end do
!! end do
!! 
!! do i=1,nel
!!    PX = 0d0
!!    PY = 0d0
!!    PZ = 0d0
!!    do j=1,8
!!       PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
!!       PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
!!       PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
!!    end do
!!    call isinelementid(px,py,pz,ibodyc,in)
!!    if(in .gt. 0)then
!!       ! we found an element set the data and go on to the next body
!!       ! if(myid.eq.5)write(*,*)'incremental: setting elem dof ',i
!!       call setelement(i,ibodyc)
!!       return
!!    end if
!! end do
!! 
!! end subroutine fbm_FictPre2Body
!!
!! ----------------------------------------------
!!
!! SUBROUTINE fbm_Knpr_Perf(dcorvg,kvert,kedge,karea)
!! USE var_QuadScalar,ONLY:myExportExt,elementsvisited 
!! use pp3d_mpi, only: myid
!! REAL*8  dcorvg(3,*)
!! INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!! INTEGER NeighE(2,12),NeighA(4,6)
!! real*8 vertices(3,8)
!! logical, dimension(:), allocatable :: visited
!! integer :: ibody,i,isboundary
!! real*8 point(3)
!! integer :: iel,idofs,dofcount,dof,testdof
!! integer :: ivert,iedge,iarea,ielem,ibodyc,in
!! integer :: ivt1,ivt2,ivt3,ivt4,nnel,itotal
!! integer, dimension(:), allocatable :: elements
!! integer, dimension(:), allocatable :: idofselement
!! logical :: startFound
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! FictKNPR=0
!! allocate(visited(NEL))
!! 
!! ! initialize the init array
!! 
!! ! loop over all bodies and compute the elements that intersect the 
!! ! bodies
!! DO ibody=1,myFBM%nParticles
!!  elementsvisited=0
!!  isboundary = 0
!!  ibodyc=ibody-1
!!  idofs = 0
!!  dofcount = 0 
!!  dof = 0
!!  do i=1,NEL
!!   visited(i)=.false.
!!  end do
!! 
!! 
!!  ! if the body is a boundary then skip
!!  call isboundarybody(isboundary,ibodyc)
!!  if(isboundary .eq. 1) cycle
!! 
!!   if(fbm_ContainsBdryElement(ibodyc))then
!!    call fbm_Knpr_PerfBodyBdr(dcorvg,kvert,kedge,karea,ibody)
!!   end if
!! 
!!  ! traverse the element list until we found
!!  ! an element that inside has dofs inside
!!  call getelement(iel,ibodyc)
!! 
!!  ! if the body has no element in the domain
!!  if(iel .eq. 0)cycle 
!! 
!!  ! check if this is a special case
!!  ! if the element had boundary elements inside that are
!!  ! on a subdomain boundary we call a special routine  
!! 
!!  ! get the vertices of the element found to be inside
!!  vertices(:,1) = dcorvg(:,kvert(1,iel))
!!  vertices(:,2) = dcorvg(:,kvert(2,iel))
!!  vertices(:,3) = dcorvg(:,kvert(3,iel))
!!  vertices(:,4) = dcorvg(:,kvert(4,iel))
!!  vertices(:,5) = dcorvg(:,kvert(5,iel))
!!  vertices(:,6) = dcorvg(:,kvert(6,iel))
!!  vertices(:,7) = dcorvg(:,kvert(7,iel))
!!  vertices(:,8) = dcorvg(:,kvert(8,iel))
!! 
!!  ! check each dof of the current element
!!  call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!! 
!!   ! if there are dofs inside start the recursive search from there
!!   if(idofs.eq.27)then
!!     call clearelementlists(ibodyc)
!!     !write(*,*)'start recursive search: '
!!     call addelement2list(iel,ibodyc)
!!     call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else if(idofs.gt.0)then
!!     call clearelementlists(ibodyc)
!!     !write(*,*)'dofs inside start element: ',idofs
!!     !write(*,*)'start recursive search: '
!!     call addelement2bndlist(iel,dof,ibodyc)
!!     myExportExt%p_DataScalarCell(3)%pData(iel)=real(ibody)
!!     call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else
!!     !write(*,*)'no dof inside old start element, trying to determine a new start element',myid
!!     ! loop through the element list and test until we found an element that is inside
!!     ! *** Loop over the boundary elements
!!     ! get the number of boundary elements and loop 
!!     nnel = 0
!!     startFound = .false.
!!     call getelementsbndry(nnel,ibodyc)
!!     if(allocated(elements))deallocate(elements)
!!     if(allocated(idofselement))deallocate(idofselement)
!!     allocate(elements(nnel))
!!     allocate(idofselement(nnel))
!!     call getelementarray(elements,idofselement,ibodyc)
!!     do i=1,nnel
!!      iel = elements(i)
!!      dof=0
!!      idofs=0
!!      vertices(:,1) = dcorvg(:,kvert(1,iel))
!!      vertices(:,2) = dcorvg(:,kvert(2,iel))
!!      vertices(:,3) = dcorvg(:,kvert(3,iel))
!!      vertices(:,4) = dcorvg(:,kvert(4,iel))
!!      vertices(:,5) = dcorvg(:,kvert(5,iel))
!!      vertices(:,6) = dcorvg(:,kvert(6,iel))
!!      vertices(:,7) = dcorvg(:,kvert(7,iel))
!!      vertices(:,8) = dcorvg(:,kvert(8,iel))
!!      call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!!      if(idofs.eq.27)then
!!       call setelement(iel,ibodyc)
!!       call clearelementlists(ibodyc)
!!       !write(*,*)'ok,found...start recursive search: ',myid
!!       call addelement2list(iel,ibodyc)
!!       call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!       startFound = .true.
!!       exit
!!      else if(idofs.gt.0)then
!!       call setelement(iel,ibodyc)
!!       call clearelementlists(ibodyc)
!!        !write(*,*)'dofs inside new start element: ',idofs,myid
!!        !write(*,*)'ok,... start recursive search: ',myid
!!       call addelement2bndlist(iel,dof,ibodyc)
!!       myExportExt%p_DataScalarCell(3)%pData(iel)=real(ibody)
!!       call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!       startFound = .true.
!!       exit
!!      end if
!!     end do ! end do i
!!     if(startFound)then   
!!       !write(*,*)'found new start element...'
!!     else
!!      ! we could not find a suitable start element,
!!      ! which can only mean that the body has left the domain
!!      iel = 0
!!      call setelement(iel,ibodyc)
!!      call clearelementlists(ibodyc)
!!     end if
!!   end if
!!  !if(myid.eq.1)write(*,*)'end recursive search... elementsvisited: ', elementsvisited,NEL
!!  call getelementsinside(elementsvisited,ibodyc)
!!  !write(*,*)'Elements inside: ',elementsvisited,' for body: ' ,ibody,myid
!!  call getelementsbndry(elementsvisited,ibodyc)
!!  !write(*,*)'Elements on boundary: ',elementsvisited,myid
!!  call gettotalelements(itotal, ibodyc)
!!  !write(*,*)'itotal: ',itotal,myid
!!  if(itotal.eq.0)then
!!    iel = 0
!!    call setelement(iel,ibodyc)
!!    call clearelementlists(ibodyc)
!!  end if
!! END DO ! end for all bodies
!! 
!! END SUBROUTINE fbm_Knpr_Perf
!!
!! ----------------------------------------------
!!
!! RECURSIVE SUBROUTINE fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,kadj)
!! USE var_QuadScalar,ONLY:myExportExt
!! use pp3d_mpi, only:myid 
!! REAL*8  dcorvg(3,*)
!! INTEGER kvert(8,*),kedge(12,*),karea(6,*),kadj(6,*)
!! REAL*8 PX,PY,PZ,DIST
!! INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,ibody
!! INTEGER NeighE(2,12),NeighA(4,6)
!! real*8 nvertices(3,8)
!! real*8 point(3)
!! integer :: iel,idofs,dofcount,dof,testdof,testidofs
!! integer :: iae,ineighbour
!! integer :: ivert,iedge,iarea,ielem,ibodyc,in
!! logical, dimension(*) :: visited
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! ! TODO: still need to consider the cases when the real boundary or
!! ! the subdomain boundary is hit
!!  
!!  visited(iel)=.true.
!! 
!!  ibodyc = ibody - 1
!! 
!!  elementsvisited = elementsvisited + 1
!! 
!! do iae=1,6
!!   idofs  = 0
!!   dofcount = 0 
!!   dof = 0
!! 
!!   ineighbour = KWORK(L(KLADJ(NLMAX))+(iel-1)*NNAE+iae)
!!   ineighbour = kadj(iae,iel)
!! 
!!   if(ineighbour .eq. 0)then
!!    !write(*,*)'skip 0 ineighbour: ',ineighbour
!!    cycle
!!   end if
!! 
!!   ! check if the neighbour was already visited
!!   if(visited(ineighbour))then
!!    !write(*,*)'skip visited ineighbour: ',ineighbour
!!     cycle
!!   end if
!! 
!!   ! calculate the number of vertices that are inside
!!   ! add them to the partices element list
!!   ! and set the corresponding degrees of freedom to the 
!!   ! particle id     
!!   nvertices(:,1) = dcorvg(:,kvert(1,ineighbour))
!!   nvertices(:,2) = dcorvg(:,kvert(2,ineighbour))
!!   nvertices(:,3) = dcorvg(:,kvert(3,ineighbour))
!!   nvertices(:,4) = dcorvg(:,kvert(4,ineighbour))
!!   nvertices(:,5) = dcorvg(:,kvert(5,ineighbour))
!!   nvertices(:,6) = dcorvg(:,kvert(6,ineighbour))
!!   nvertices(:,7) = dcorvg(:,kvert(7,ineighbour))
!!   nvertices(:,8) = dcorvg(:,kvert(8,ineighbour))
!!   
!! !  write(*,*)'--------------------------------Knpr_Element---------------------------------------'
!!    call fbm_CheckKnpr_Element(nvertices,ineighbour,ibody,dof,idofs)
!! 
!!   ! if there are dofs inside start the recursive search from there
!!   if(idofs.eq.27)then
!!     call addelement2list(ineighbour,ibodyc)
!!     !write(*,*)'completely inside',   ineighbour
!!     call fbm_FBM_Element(ibody,ineighbour,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else if(idofs.gt.0)then
!!     myExportExt%p_DataScalarCell(3)%pData(ineighbour)=1
!!     !write(*,*)'add2 bndryelement list...',ineighbour
!!     call addelement2bndlist(ineighbour,dof,ibodyc)
!!     call fbm_FBM_Element(ibody,ineighbour,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else
!!    !write(*,*)'no dof inside...',ineighbour
!!   end if
!! 
!! end do ! end loop over neighbours
!! 
!! END SUBROUTINE fbm_FBM_Element
!!
!! ----------------------------------------------
!!
!SUBROUTINE fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!  use pp3d_mpi, only: myid
!  implicit none
!  real *8, dimension(3,8) :: vertices
!  integer :: iel
!  integer :: ibody
!  integer :: dof
!  integer :: idofs
!  ! local variables
!  real*8  :: point(3)
!  integer :: dofcount
!  integer :: iae,ineighbour,indexdof
!  integer :: ivert,iedge,iarea,ielem,ibodyc,in
!  !REAL*8, dimension(3,*) :: dcorvg
!  !INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
!  integer, dimension(2,12) ::  NeighE
!  integer, dimension(4,6)  ::  NeighA
!  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!
!  ibodyc = ibody-1
!  dof = 0
!  idofs = 0
!  dofcount = 0
!
!  ! calculate the dofs inside
!  do ivert=1,8
!  in=0
!  call isinelementid(vertices(1,ivert),vertices(2,ivert),vertices(3,ivert),ibodyc,in)
!  ! dof is inside, set the corresponding degree of freedom
!  if(in .gt. 0)then 
!    indexdof = kwork(l(lvert)+(iel-1)*8+ivert-1)
!    FictKNPR(indexdof)=ibody
!    idofs = idofs + 1
!    dof = ior(dof,2**dofcount)
!  else
!    indexdof = kwork(l(lvert)+(iel-1)*8+ivert-1)
!  end if
!  dofcount = dofcount + 1
!  end do 
!
!  ! calculate the edge dofs inside
!  do iedge=1,12
!  in=0
!  point(:)= 0.5 * (vertices(:,NeighE(1,iedge)) + vertices(:,NeighE(2,iedge)))  
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    indexdof = nvt + kwork(l(ledge)+(iel-1)*12+iedge-1)
!    FictKNPR(indexdof)=ibody
!    idofs = idofs + 1
!    dof = ior(dof,2**dofcount)
!  else
!    indexdof = nvt + kwork(l(ledge)+(iel-1)*12+iedge-1)
!  end if
!  dofcount = dofcount + 1
!  end do 
!
!  ! calculate the face dofs inside
!  do iarea=1,6
!  in=0
!  point(:)= 0.25 * (vertices(:,NeighA(1,iarea)) + vertices(:,NeighA(2,iarea)) + vertices(:,NeighA(3,iarea)) + vertices(:,NeighA(4,iarea)))
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    indexdof = nvt+net+kwork(l(larea)+(iel-1)*6+iarea-1)
!    FictKNPR(indexdof)=ibody
!    idofs = idofs + 1
!    dof = ior(dof,2**dofcount)
!  else
!    indexdof = nvt+net+kwork(l(larea)+(iel-1)*6+iarea-1)
!  end if
!  dofcount = dofcount + 1
!  end do 
!
!  ! calculate the element dof
!  in=0
!  point(:)= 0.125 * (vertices(:,1) + vertices(:,2) + vertices(:,3) + vertices(:,4) + &
!    vertices(:,5) + vertices(:,6) + vertices(:,7) + vertices(:,8))                   
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    indexdof=nvt+net+nat+iel
!    FictKNPR(nvt+net+nat+iel)=ibody
!    idofs = idofs + 1
!    dof = ior(dof,2**dofcount)
!  else
!    indexdof=nvt+net+nat+iel
!  end if
!
!END SUBROUTINE fbm_CheckKnpr_Element
!!
!! ----------------------------------------------
!!
!! subroutine fbm_RemoteDomains()
!! USE var_QuadScalar,ONLY:myExportExt,elementsvisited 
!! use pp3d_mpi, only:subnodes,myid,mgBoundingBox
!! implicit none
!! INTEGER NeighE(2,12),NeighA(4,6)
!! integer :: ibody,i,isboundary,inter,idom,itotal,iprev,ibodyc
!! integer, dimension(:), allocatable :: elements
!! integer, dimension(:), allocatable :: idofselement
!! integer, dimension(subnodes) :: icandidates
!! logical :: bIntersection
!! integer :: ichecks
!! 
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! icandidates=-1
!! ichecks    = 0
!! 
!! ! loop over all bodies and compute the elements that intersect the 
!! ! bodies
!! DO ibody=1,myFBM%nParticles
!!   elementsvisited=0
!!   isboundary = 0
!!   ibodyc=ibody-1
!! 
!!   ! get the local info on the body
!!   call gettotalelements(itotal, ibodyc)
!!   call getelementsprev(iprev, ibodyc)
!! 
!! !  if(iprev.gt.0).and.(itotal.eq.0)then
!! !    call clearelementlists(ibodyc)
!! !  end if
!! 
!!   ! if the body is known in the domain
!!   if(iprev.gt.0)then
!!     call updateelementsprev(ibodyc)
!!     cycle
!!   end if 
!! 
!!   ! We want to determine whether the current body could have entered the domain.  
!!   bIntersection =  fbm_IntersectCoarseGrid(knvt(nlmin),knel(nlmin),dwork(l(klcvg(nlmin))),kwork(l(klvert(nlmin))),ibodyc) 
!! 
!!   ! check if body new in domain
!!   if(bIntersection)then
!!      ! write(*,*)'checking candidate: body/myid',ibody,myid
!!      ! if new in domain start FictPre for it
!!      call fbm_FictPreBodyAll(dwork(l(lcorvg)),kwork(l(lvert)),kwork(l(ledge)),kwork(l(larea)),ibody)
!!      ! after we ran FictPre we can start Knpr_PerfBody
!!      call fbm_Knpr_PerfBodyBdr(dwork(l(lcorvg)),kwork(l(lvert)),kwork(l(ledge)),kwork(l(larea)),ibody)
!!      ichecks = ichecks + 1
!!      call updateelementsprev(ibodyc)
!!   else
!!      ! has left the domain and the domain BVH
!!      call updateelementsprev(ibodyc)
!!   end if
!!  
!! END DO ! end for all bodies
!! !write(*,*)'Total checks: ',ichecks,myid
!! end subroutine
!!
!! ----------------------------------------------
!!
!subroutine fbm_FictPreBody(dcorvg,kvert,kedge,karea,ibody)
!  use var_QuadScalar,only:myExport,myFBM 
!  REAL*8  dcorvg(3,*)
!  integer kvert(8,*),kedge(12,*),karea(6,*)
!  integer :: ibody
!  REAL*8 PX,PY,PZ,DIST
!  integer i,j,k,ivt1,ivt2,ivt3,ivt4
!  integer NeighE(2,12),NeighA(4,6)
!  real*8 nvertices(3,8)
!  real*8 point(3)
!  integer :: iel,idofs
!  integer :: iae,ineighbour
!  integer :: ivert,iedge,iarea,ielem,ibodyc,in,myelement
!  logical :: found
!  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!
!
!  myelement = 0
!  ibodyc = ibody - 1
!  found  = .false.
!  call setelement(myelement,ibodyc)
!
!  do i=1,NEL
!  nvertices(:,1) = dcorvg(:,kvert(1,i))
!  nvertices(:,2) = dcorvg(:,kvert(2,i))
!  nvertices(:,3) = dcorvg(:,kvert(3,i))
!  nvertices(:,4) = dcorvg(:,kvert(4,i))
!  nvertices(:,5) = dcorvg(:,kvert(5,i))
!  nvertices(:,6) = dcorvg(:,kvert(6,i))
!  nvertices(:,7) = dcorvg(:,kvert(7,i))
!  nvertices(:,8) = dcorvg(:,kvert(8,i))
!
!  ! calculate the dofs inside
!  do ivert=1,8
!  in=0
!  call isinelementid(nvertices(1,ivert),nvertices(2,ivert),nvertices(3,ivert),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!    exit
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found)exit
!
!  ! calculate the edge dofs inside
!  do iedge=1,12
!  in=0
!  ivt1 = kvert(NeighE(1,iedge),i)
!  ivt2 = kvert(NeighE(2,iedge),i)
!  point(:)= 0.5 * (dcorvg(:,ivt1) + dcorvg(:,ivt2))
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found)exit
!
!  ! calculate the face dofs inside
!  do iarea=1,6
!  in=0
!  ivt1 = kvert(NeighE(1,iarea),i)
!  ivt2 = kvert(NeighE(2,iarea),i)
!  ivt3 = kvert(NeighA(3,iarea),i)
!  ivt4 = kvert(NeighA(4,iarea),i)
!  point(:)= 0.25 * (dcorvg(:,ivt1) + dcorvg(:,ivt2) + dcorvg(:,ivt3) + dcorvg(:,ivt4))
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!  end do 
!
!  ! if we found an element exit the loop
!  if(found .eqv. .true.)exit
!
!  ! calculate the element dof
!  in=0
!  point(:)= 0.125 * (nvertices(:,1) + nvertices(:,2) + nvertices(:,3) + nvertices(:,4) + &
!    nvertices(:,5) + nvertices(:,6) + nvertices(:,7) + nvertices(:,8))                   
!  call isinelementid(point(1),point(2),point(3),ibodyc,in)
!  if(in .gt. 0)then
!    ! we found an element set the data and go on to the next body
!    found = .true.
!    myelement = i
!    !write(*,*)'found element: ',i,'for body: ',ibodyc,' in sub: ',myid
!    call setelement(myelement,ibodyc)
!  end if
!
!  end do ! end for all elements
!
!end subroutine fbm_FictPreBody
!!
!!****************************************************************************  
!!
!subroutine fbm_elementsAtVertex(kvert)
!  !**@Description
!  ! In this function we calculate the number of elements
!  ! that are adjacent to a vertex 
!  !**@EndDescription
!
!  USE def_FEAT
!
!  implicit none
!
!  INTEGER kvert(8,*)
!  INTEGER i,j,ivt
!
!  allocate(mykvel(NVT))
!
!  ! for all elements
!  do i=1,nel
!  ! for all vertices at the element
!  do j=1,8
!  ivt = kvert(j,i)
!  mykvel(ivt) = i
!  end do 
!  end do
!
!end subroutine
!!
!!****************************************************************************  
!!
!logical function fbm_IntersectCoarseGrid(invt,inel,dcorvgcoarse,kvertcoarse,ibody) 
!  !**@Description
!  !  Checks for an intersection of a body with the coarse grid
!  !**@EndDescription
!
!  implicit none
!  real*8  dcorvgcoarse(3,*)
!  integer kvertcoarse(8,*)
!  integer :: invt
!  integer :: inel
!  integer :: ibody
!
!  ! local variables
!  integer :: iel,i,intersection
!  real*8, dimension(3,2) :: dMinMax
!
!  do iel=1,inel
!  ! compute bounding box of element
!  dMinMax(:,1)=dcorvgcoarse(:,kvertcoarse(1,iel))
!  dMinMax(:,2)=dcorvgcoarse(:,kvertcoarse(1,iel))
!  do i=2,8
!  if(dMinMax(1,1).gt.dcorvgcoarse(1,kvertcoarse(i,iel)))then
!    dMinMax(1,1)=dcorvgcoarse(1,kvertcoarse(i,iel))
!  end if
!
!  if(dMinMax(2,1).gt.dcorvgcoarse(2,kvertcoarse(i,iel)))then
!    dMinMax(2,1)=dcorvgcoarse(2,kvertcoarse(i,iel))
!  end if
!
!  if(dMinMax(3,1).gt.dcorvgcoarse(3,kvertcoarse(i,iel)))then
!    dMinMax(3,1)=dcorvgcoarse(3,kvertcoarse(i,iel))
!  end if
!
!  if(dMinMax(1,2).lt.dcorvgcoarse(1,kvertcoarse(i,iel)))then
!    dMinMax(1,2)=dcorvgcoarse(1,kvertcoarse(i,iel))
!  end if
!
!  if(dMinMax(2,2).lt.dcorvgcoarse(2,kvertcoarse(i,iel)))then
!    dMinMax(2,2)=dcorvgcoarse(2,kvertcoarse(i,iel))
!  end if
!
!  if(dMinMax(3,2).lt.dcorvgcoarse(3,kvertcoarse(i,iel)))then
!    dMinMax(3,2)=dcorvgcoarse(3,kvertcoarse(i,iel))
!  end if
!  end do
!
!  ! check for intersection with coarse grid element   
!  call intersecthexbody(dMinMax,ibody,intersection)
!  if(intersection .eq. 1)then
!    fbm_IntersectCoarseGrid = .true.
!    return
!  end if
!  end do
!
!  fbm_IntersectCoarseGrid = .false.
!
!end function
!!
!!----------------------------------------------
!!
!! SUBROUTINE fbm_Knpr_PerfBody(dcorvg,kvert,kedge,karea,ibody)
!! !**@Description
!! !  Checks for an intersection of a body with the coarse grid
!! !**@EndDescription
!! USE var_QuadScalar,ONLY:myExportExt,elementsvisited 
!! use pp3d_mpi, only: myid
!! REAL*8  dcorvg(3,*)
!! INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!! integer :: ibody
!! INTEGER NeighE(2,12),NeighA(4,6)
!! 
!! real*8 vertices(3,8)
!! logical, dimension(:), allocatable :: visited
!! integer :: i,isboundary
!! real*8 point(3)
!! integer :: iel,idofs,dofcount,dof,testdof
!! integer :: ivert,iedge,iarea,ielem,ibodyc,in
!! integer :: ivt1,ivt2,ivt3,ivt4,nnel
!! integer, dimension(:), allocatable :: elements
!! integer, dimension(:), allocatable :: idofselement
!! logical :: startFound
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! ! allocate the visited array for this body
!! allocate(visited(NEL))
!! 
!!  elementsvisited=0
!!  isboundary = 0
!!  ibodyc=ibody-1
!!  idofs = 0
!!  dofcount = 0 
!!  dof = 0
!!  do i=1,NEL
!!   visited(i)=.false.
!!  end do
!! 
!! 
!!  ! if the body is a boundary then skip
!!  call isboundarybody(isboundary,ibodyc)
!!  if(isboundary .eq. 1) return
!! 
!!  ! traverse the element list until we found
!!  ! an element that is inside
!! 
!!  ! erase the old element list
!!  ! TODO: alternative ~~ check the old elements that are
!!  ! inside, mark them as visited and run the
!!  ! algorithm
!!  call getelement(iel,ibodyc)
!! 
!!  if(iel .eq. 0)return 
!! 
!!  ! write(*,*)'start element :',iel,myid
!!  
!!  ! get the vertices of the element found to be inside
!!  vertices(:,1) = dcorvg(:,kvert(1,iel))
!!  vertices(:,2) = dcorvg(:,kvert(2,iel))
!!  vertices(:,3) = dcorvg(:,kvert(3,iel))
!!  vertices(:,4) = dcorvg(:,kvert(4,iel))
!!  vertices(:,5) = dcorvg(:,kvert(5,iel))
!!  vertices(:,6) = dcorvg(:,kvert(6,iel))
!!  vertices(:,7) = dcorvg(:,kvert(7,iel))
!!  vertices(:,8) = dcorvg(:,kvert(8,iel))
!! 
!!  call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!! 
!!  ! write(*,*)'---------------------------------------------------------------'
!! 
!!   ! if there are dofs inside start the recursive search from there
!!   if(idofs.eq.27)then
!!     call clearelementlists(ibodyc)
!!     !write(*,*)'start recursive search: '
!!     call addelement2list(iel,ibodyc)
!!     call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else if(idofs.gt.0)then
!!     call clearelementlists(ibodyc)
!!     !write(*,*)'dofs inside start element: ',idofs
!!     !write(*,*)'start recursive search: '
!!     call addelement2bndlist(iel,dof,ibodyc)
!!     myExportExt%p_DataScalarCell(3)%pData(iel)=real(ibody)
!!     call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!   else
!!     ! write(*,*)'no dof inside old start element, trying to determine a new start element',myid
!!     ! loop through the element list and test until we found an element that is inside
!!     ! *** Loop over the boundary elements
!!     ! get the number of boundary elements and loop 
!!     nnel = 0
!!     startFound = .false.
!!     call getelementsbndry(nnel,ibodyc)
!!     if(allocated(elements))deallocate(elements)
!!     if(allocated(idofselement))deallocate(idofselement)
!!     allocate(elements(nnel))
!!     allocate(idofselement(nnel))
!!     call getelementarray(elements,idofselement,ibodyc)
!!     do i=1,nnel
!!      iel = elements(i)
!!      dof=0
!!      idofs=0
!!      vertices(:,1) = dcorvg(:,kvert(1,iel))
!!      vertices(:,2) = dcorvg(:,kvert(2,iel))
!!      vertices(:,3) = dcorvg(:,kvert(3,iel))
!!      vertices(:,4) = dcorvg(:,kvert(4,iel))
!!      vertices(:,5) = dcorvg(:,kvert(5,iel))
!!      vertices(:,6) = dcorvg(:,kvert(6,iel))
!!      vertices(:,7) = dcorvg(:,kvert(7,iel))
!!      vertices(:,8) = dcorvg(:,kvert(8,iel))
!!      call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!!      if(idofs.eq.27)then
!!       call setelement(iel,ibodyc)
!!       call clearelementlists(ibodyc)
!!       ! write(*,*)'ok,found...start recursive search: ',myid
!!       call addelement2list(iel,ibodyc)
!!       call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!       startFound = .true.
!!       exit
!!      else if(idofs.gt.0)then
!!       call setelement(iel,ibodyc)
!!       call clearelementlists(ibodyc)
!!        ! write(*,*)'dofs inside new start element: ',idofs,myid
!!        ! write(*,*)'ok,... start recursive search: ',myid
!!       call addelement2bndlist(iel,dof,ibodyc)
!!       myExportExt%p_DataScalarCell(3)%pData(iel)=real(ibody)
!!       call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!       startFound = .true.
!!      exit
!!      end if
!!     end do ! end do i
!!     if(startFound)then   
!!       ! write(*,*)'found new start element...'
!!     else
!!      ! we could not find a suitable start element,
!!      ! which can only mean that the body has left the domain
!!      ! write(*,*)'could not find start element...'
!!      iel = 0
!!      call setelement(iel,ibodyc)
!!      call clearelementlists(ibodyc)
!!     end if
!!   end if
!! 
!!  !if(myid.eq.1)write(*,*)'end recursive search... elementsvisited: ', elementsvisited,NEL
!!  call getelementsinside(elementsvisited,ibodyc)
!!  ! write(*,*)'Elements inside: ',elementsvisited,' for body: ' ,ibody,myid
!!  call getelementsbndry(elementsvisited,ibodyc)
!!  ! write(*,*)'Elements on boundary: ',elementsvisited,myid
!!  call gettotalelements(itotal, ibodyc)
!!  ! write(*,*)'itotal: ',itotal,myid
!!  if(itotal.eq.0)then
!!    iel = 0
!!    call setelement(iel,ibodyc)
!!    call clearelementlists(ibodyc)
!!  end if
!! 
!! END SUBROUTINE fbm_Knpr_PerfBody
!!
!!****************************************************************************  
!!
!! subroutine fbm_Knpr_PerfBodyBdr(dcorvg,kvert,kedge,karea,ibody)
!! !**@Description
!! !  Checks for an intersection of a body with the coarse grid
!! !**@EndDescription
!! use var_QuadScalar,only:myExportExt,elementsvisited 
!! use pp3d_mpi, only: myid
!! implicit none
!! REAL*8  dcorvg(3,*)
!! integer kvert(8,*),kedge(12,*),karea(6,*)
!! integer :: ibody
!! integer NeighE(2,12),NeighA(4,6)
!! 
!! real*8 vertices(3,8)
!! logical, dimension(:), allocatable :: visited
!! integer :: i,isboundary
!! real*8 point(3)
!! integer :: iel,idofs,dofcount,dof,testdof
!! integer :: ivert,iedge,iarea,ielem,ibodyc,in
!! integer :: ivt1,ivt2,ivt3,ivt4,nnel,itotal
!! integer, dimension(:), allocatable :: elements
!! logical :: startFound
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! ! allocate the visited array for this body
!! allocate(visited(NEL))
!! 
!!  elementsvisited=0
!!  isboundary = 0
!!  ibodyc=ibody-1
!!  idofs = 0
!!  dofcount = 0 
!!  dof = 0
!!  do i=1,NEL
!!   visited(i)=.false.
!!  end do
!! 
!!  ! if the body is a boundary then skip
!!  call isboundarybody(isboundary,ibodyc)
!!  if(isboundary .eq. 1) return
!! 
!!  call clearelementlists(ibodyc)
!!  call fbm_FictPreBodyAll(dwork(l(lcorvg)),kwork(l(lvert)),kwork(l(ledge)),kwork(l(larea)),ibody)
!! 
!!  ! get a start element
!!  call getelement(iel,ibodyc)
!! 
!!  ! if the body is not in the domain we exit the routine
!!  if(iel .eq. 0)return 
!! 
!!  ! get the vertices of the element found to be inside
!!  vertices(:,1) = dcorvg(:,kvert(1,iel))
!!  vertices(:,2) = dcorvg(:,kvert(2,iel))
!!  vertices(:,3) = dcorvg(:,kvert(3,iel))
!!  vertices(:,4) = dcorvg(:,kvert(4,iel))
!!  vertices(:,5) = dcorvg(:,kvert(5,iel))
!!  vertices(:,6) = dcorvg(:,kvert(6,iel))
!!  vertices(:,7) = dcorvg(:,kvert(7,iel))
!!  vertices(:,8) = dcorvg(:,kvert(8,iel))
!! 
!!  call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!! 
!!  ! get the total number of elements
!!  call gettotalelements(itotal, ibodyc)
!! 
!!  ! allocate the element array
!!  nnel = itotal
!!  startFound = .false.
!!  if(allocated(elements))deallocate(elements)
!!  allocate(elements(nnel))! max (inner,bndry)
!! 
!!  ! get the element array
!!  call getallelements(elements,ibodyc)
!!  call clearelementlists(ibodyc)
!! 
!!  ! start the PerfSearch from all the old
!!  ! bndry elements
!!  do i=1,nnel
!!     iel = elements(i)
!!   
!!     ! check if this element was already visited
!!     if(visited(iel))cycle
!! 
!!     ! reset the dofs
!!     dof=0
!!     idofs=0
!! 
!!     ! get the element vertices
!!     vertices(:,1) = dcorvg(:,kvert(1,iel))
!!     vertices(:,2) = dcorvg(:,kvert(2,iel))
!!     vertices(:,3) = dcorvg(:,kvert(3,iel))
!!     vertices(:,4) = dcorvg(:,kvert(4,iel))
!!     vertices(:,5) = dcorvg(:,kvert(5,iel))
!!     vertices(:,6) = dcorvg(:,kvert(6,iel))
!!     vertices(:,7) = dcorvg(:,kvert(7,iel))
!!     vertices(:,8) = dcorvg(:,kvert(8,iel))
!!     call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!!     if(idofs.eq.27)then
!!        call setelement(iel,ibodyc)
!!        ! write(*,*)'ok,found...start recursive search: ',myid
!!        call addelement2list(iel,ibodyc)
!!        call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!        startFound = .true.
!!     else if(idofs.gt.0)then
!!        ! write(*,*)'dofs inside new start element: ',idofs,myid
!!        ! write(*,*)'ok,... start recursive search: ',myid
!!        call addelement2bndlist(iel,dof,ibodyc)
!!        myExportExt%p_DataScalarCell(3)%pData(iel)=real(ibody)
!!        call fbm_FBM_Element(ibody,iel,dcorvg,kvert,visited,kedge,karea,KWORK(L(KLADJ(NLMAX))))
!!        startFound = .true.
!!     end if
!!  end do ! end do i
!! 
!!  ! if(startFound)then   
!!  !    ! write(*,*)'found new start element...'
!!  ! else
!!  !    ! we could not find a suitable start element,
!!  !    ! which can only mean that the body has left the domain
!!  !    ! write(*,*)'could not find start element...'
!!  !    iel = 0
!!  !    call setelement(iel,ibodyc)
!!  !    call clearelementlists(ibodyc)
!!  ! end if
!! 
!!  ! check the new number of elements
!!  call gettotalelements(itotal, ibodyc)
!!  ! write(*,*)'itotal: ',itotal,myid
!!  if(itotal.eq.0)then
!!    iel = 0
!!    call setelement(iel,ibodyc)
!!    call clearelementlists(ibodyc)
!!  end if
!! 
!! end subroutine fbm_Knpr_PerfBodyBdr
!!
!!****************************************************************************  
!!
!logical function fbm_ContainsBdryElement(ibody)
!  !**@Description
!  !  Checks for an intersection of a body with the coarse grid
!  !**@EndDescription
!  implicit none
!
!  integer :: ibody
!
!  ! local variables
!  integer :: nnel,i
!  integer, dimension(:), allocatable :: elements
!
!
!  fbm_ContainsBdryElement = .true.
!  return
!
!  ! allocate the element array
!  call gettotalelements(nnel, ibody)
!  if(allocated(elements))deallocate(elements)
!  allocate(elements(nnel))
!
!  ! get the element array
!  call getallelements(elements,ibody)
!
!  do i=1,nnel
!  if(elements(i))then
!    ! 
!    fbm_ContainsBdryElement = .true.
!    return
!  end if
!  end do
!
!  fbm_ContainsBdryElement = .false.
!
!end function fbm_ContainsBdryElement
!!
!! ----------------------------------------------
!!
!! subroutine fbm_FictPreBodyAll(dcorvg,kvert,kedge,karea,ibody)
!! use var_QuadScalar,only:myExport,myFBM 
!! REAL*8  dcorvg(3,*)
!! integer kvert(8,*),kedge(12,*),karea(6,*)
!! integer :: ibody
!! REAL*8 PX,PY,PZ,dist
!! integer i,j,k,ivt1,ivt2,ivt3,ivt4
!! integer NeighE(2,12),NeighA(4,6)
!! real*8 nvertices(3,8)
!! real*8 point(3)
!! real*8  :: Center(3)
!! integer :: iel,idofs
!! integer :: iae,ineighbour
!! integer :: ivert,iedge,iarea,ielem,ibodyc,in,myelement,dof
!! logical :: found
!! DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!! DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!! 
!! 
!!   myelement = 0
!!   ibodyc = ibody - 1
!!   found  = .false.
!!   ! when we call this routine we want to invalidate the
!!   ! old settings for start element and element list
!!   call setelement(myelement,ibodyc)
!!   call clearelementlists(ibodyc)
!!   
!!   do i=1,NEL
!!     idofs=0
!!     dof=0
!!     nvertices(:,1) = dcorvg(:,kvert(1,i))
!!     nvertices(:,2) = dcorvg(:,kvert(2,i))
!!     nvertices(:,3) = dcorvg(:,kvert(3,i))
!!     nvertices(:,4) = dcorvg(:,kvert(4,i))
!!     nvertices(:,5) = dcorvg(:,kvert(5,i))
!!     nvertices(:,6) = dcorvg(:,kvert(6,i))
!!     nvertices(:,7) = dcorvg(:,kvert(7,i))
!!     nvertices(:,8) = dcorvg(:,kvert(8,i))
!! 
!!     point(:)= 0.125 * (nvertices(:,1) + nvertices(:,2) + nvertices(:,3) + nvertices(:,4) + &
!!                        nvertices(:,5) + nvertices(:,6) + nvertices(:,7) + nvertices(:,8))                   
!!    
!!     Center = myFBM%particleNew(ibody)%Position
!! 
!!     point = point - Center
!!     dist  = sqrt(point(1)**2+point(2)**2+point(3)**2)
!! 
!!     ! if the distance is larger than 3*diameter then skip the element
!!     if(dist .gt. 6d0*myFBM%ParticleNew(ibody)%sizes(1))cycle
!!  
!!     call fbm_CheckKnpr_Element(nvertices,i,ibody,dof,idofs)
!! 
!!     ! if there are dofs inside start the recursive search from there
!!     if(idofs.eq.27)then
!!       call addelement2list(i,ibodyc)
!!       if(.not.found)then
!!         call setelement(i,ibodyc)
!!         found=.true.
!!       end if
!!     else if(idofs.gt.0)then
!!       if(.not.found)then
!!         call setelement(i,ibodyc)
!!         found=.true.
!!       end if
!!       call addelement2bndlist(i,dof,ibodyc)
!!     end if
!! 
!!   end do ! end for all elements
!! 
!! end subroutine fbm_FictPreBodyAll
!!
!! ----------------------------------------------
!!
!SUBROUTINE fbm_Knpr_UniformGrid(dcorvg,kvert,kedge,karea)
!  !USE var_QuadScalar,ONLY:myExportExt
!  use pp3d_mpi, only: myid
!  REAL*8  dcorvg(3,*)
!  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  INTEGER NeighE(2,12),NeighA(4,6)
!  real*8 vertices(3,8)
!  integer :: ibody,i
!  real*8 point(3)
!  integer :: iel,idofs,dofcount,dof
!  integer :: ielem,ibodyc,indomain,iid
!  integer :: nnel,itotal
!  integer, dimension(:), allocatable :: elements
!  integer, dimension(:), allocatable :: idofselement
!  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
!  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
!
!  FictKNPR=0
!
!  ! loop over all bodies and compute the elements that intersect the 
!  ! bodies
!  do ibody=1,myFBM%nParticles
!  ibodyc=ibody-1
!  idofs = 0
!  dofcount = 0 
!  dof = 0
!
!  call queryuniformgrid(ibodyc)
!
!  end do
!
!
!
!  !   ! check whether the domain box is intersected
!  !   ! we can skip the following tests if the box is
!  !   ! not intersected
!  !   indomain = 0
!  !   iid = myid
!
!  !   ! clear the elementlist 
!  !   call clearelementlists(ibodyc)  
!
!  !   ! check for intersection
!  !   call intersectdomainbody(ibodyc, iid, indomain)
!  !   if(indomain .eq. 0)cycle
!
!  !   ! determine the elements that need to be checked
!  !   call queryuniformgrid(ibodyc)
!
!  !   ! get the number of elements
!  !   call getelementsinside(nnel,ibodyc);
!
!  !   ! get the element array
!  !   if(allocated(elements))deallocate(elements)
!  !   allocate(elements(nnel))
!  !   call getelements(elements,ibodyc)
!
!  !   ! check each element
!  !   do i=1,nnel
!  !     ! check each dof of the current element
!  !     iel = elements(i)
!  !     vertices(:,1) = dcorvg(:,kvert(1,iel))
!  !     vertices(:,2) = dcorvg(:,kvert(2,iel))
!  !     vertices(:,3) = dcorvg(:,kvert(3,iel))
!  !     vertices(:,4) = dcorvg(:,kvert(4,iel))
!  !     vertices(:,5) = dcorvg(:,kvert(5,iel))
!  !     vertices(:,6) = dcorvg(:,kvert(6,iel))
!  !     vertices(:,7) = dcorvg(:,kvert(7,iel))
!  !     vertices(:,8) = dcorvg(:,kvert(8,iel))    
!  !     call fbm_CheckKnpr_Element(vertices,iel,ibody,dof,idofs)
!  !     if((idofs.gt.0).and.(idofs.lt.27))then
!  !       call addelement2bndlist(iel,dof,ibodyc)
!  !     end if    
!  !   end do
!
!  ! end do 
!
!end subroutine fbm_Knpr_UniformGrid
!!
!! ----------------------------------------------
!!
!subroutine fbm_BuildUniformGrid(dcorvg,kvert,kedge,karea)
!  use var_QuadScalar,only:myExport,myFBM 
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!  REAL*8  dcorvg(3,*)
!  integer kvert(8,*),kedge(12,*),karea(6,*)
!  REAL*8 PX,PY,PZ,DIST
!  integer iel,j,NEL2
!  real*8 vertices(3,8)
!  real*8 center(3),minsize,maxsize,avgsize,avgsize2
!  real*8, dimension(:), allocatable :: size
!  integer, dimension(5) :: idist
!
!  avgsize=0d0
!  idist=0
!  allocate(size(NEL))
!
!  vertices(:,1) = dcorvg(:,kvert(1,1))
!  vertices(:,2) = dcorvg(:,kvert(2,1))
!  vertices(:,3) = dcorvg(:,kvert(3,1))
!  vertices(:,4) = dcorvg(:,kvert(4,1))
!  vertices(:,5) = dcorvg(:,kvert(5,1))
!  vertices(:,6) = dcorvg(:,kvert(6,1))
!  vertices(:,7) = dcorvg(:,kvert(7,1))
!  vertices(:,8) = dcorvg(:,kvert(8,1))
!  call elementsize(vertices, size(1))
!  minsize=size(1)
!  maxsize=size(1)
!
!  do iel=2,NEL
!  ! we expect an evenly distributed mesh, so we get a sample element
!  vertices(:,1) = dcorvg(:,kvert(1,iel))
!  vertices(:,2) = dcorvg(:,kvert(2,iel))
!  vertices(:,3) = dcorvg(:,kvert(3,iel))
!  vertices(:,4) = dcorvg(:,kvert(4,iel))
!  vertices(:,5) = dcorvg(:,kvert(5,iel))
!  vertices(:,6) = dcorvg(:,kvert(6,iel))
!  vertices(:,7) = dcorvg(:,kvert(7,iel))
!  vertices(:,8) = dcorvg(:,kvert(8,iel))
!  call elementsize(vertices, size(iel))
!  avgsize=avgsize+size(iel)
!  if(minsize .gt. size(iel))then
!    minsize=size(iel)
!  end if
!
!  if(maxsize .lt. size(iel))then
!    maxsize=size(iel)
!  end if
!  end do
!
!  call setelementarray(size,NEL)
!
!  !  if(myid.eq.1)then
!  !    write(*,*)'minsize,myid: ',minsize,myid
!  !    write(*,*)'maxsize,myid: ',maxsize,myid
!  !    write(*,*)'maxsize/minsize,myid: ',maxsize/minsize,myid
!  !  end if
!
!  do iel=1,NEL
!  center(:)=0d0
!  ! compute the mid point of the element
!  do j=1,8
!  center(:) = center(:) + dcorvg(:,kvert(j,iel))
!  end do
!  center(:) = center(:) * 0.125d0
!  ! insert the element into the uniform grid
!  call ug_insertelement(iel,center,size(iel))
!  end do
!
!  call ug_querystatus()
!
!end subroutine fbm_BuildUniformGrid
!!
!! ----------------------------------------------
!!
!SUBROUTINE fbm_QueryBruteForce(dcorvg,kvert,kedge,karea)
!  use def_FEAT
!  use pp3d_mpi, only: myid
!  use var_QuadScalar,only:myFBM 
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!  use fbmaux, only: fbmaux_IsInBoundingBox,fbmaux_PointInHex
!  implicit none
!  REAL*8  dcorvg(3,*)
!  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  real*8 vertices(3,8)
!  integer :: ibody,i,ibodyc,nnel,j,iconv
!  real*8 point(3)
!  integer, dimension(:), allocatable :: elements
!  logical :: found
!  real*8, dimension(3) :: center,ref
!  integer :: ichecks
!
!  ichecks=0
!
!  found = .false.
!  !write(*,*)'calling brute force with particles: ',myFBM%nParticles,nel
!
!  do ibody=1,myFBM%nParticles
!  myFBM%iel_ug(ibody)=0
!  center = myFBM%particleNew(ibody)%Position
!  ref(:) = 0d0
!
!  do j=1,nel
!
!  if(.not.fbmaux_IsInBoundingBox(dcorvg,kvert,j,center))cycle
!
!  found = fbmaux_PointInHex(center(1),center(2),center(3),&
!    ref(1),ref(2),ref(3),&
!    j,iconv,dcorvg,kvert)
!  ichecks=ichecks+1
!
!  if(found .eqv. .true.)then
!    myFBM%iel_ug(ibody)=j
!    exit
!  end if
!
!  end do
!
!  end do
!
!  !write(*,*)'total checks=',ichecks
!
!end subroutine fbm_QueryBruteForce
!!
!! ----------------------------------------------
!!
!SUBROUTINE fbm_QueryPointElement(dcorvg,kvert,kedge,karea)
!  use pp3d_mpi, only: myid
!  use var_QuadScalar,only:myFBM 
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!  implicit none
!  REAL*8  dcorvg(3,*)
!  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  real*8 vertices(3,8)
!  integer :: ibody,i,ibodyc,nnel,iel
!  real*8 point(3)
!  integer, dimension(:), allocatable :: elements
!
!  myFBM%iel_ug=0
!
!  do ibody=1,myFBM%nParticles
!
!  ibodyc=ibody-1
!
!  ! clear the list for the body
!  call clearelementlists(ibodyc)    
!
!  ! query the uniform grid
!  call queryuniformgrid(ibodyc)
!
!  ! get the number of elements
!  call getelementsinside(nnel,ibodyc);
!
!  ! the body is not in the domain -> skip
!  if(nnel .eq. 0)cycle
!
!  ! get the element array
!  if(allocated(elements))deallocate(elements)
!  allocate(elements(nnel))
!  call getelements(elements,ibodyc)
!
!  ! start element determination
!  iel = fbm_QueryPointElement_sub(dcorvg,kvert,kedge,&
!    karea,elements,&
!    nnel,ibody) 
!
!  myFBM%iel_ug(ibody)=iel
!
!  end do
!
!end subroutine fbm_QueryPointElement
!!
!!****************************************************************************  
!!
!integer function fbm_QueryPointElement_sub(dcorvg,kvert,kedge,karea,elements,nnel,ibody) 
!  !**@Description
!  ! 
!  !**@EndDescription     
!  use def_FEAT
!  use var_QuadScalar,only:myFBM 
!  use pp3d_mpi, only: myid,mgBoundingBox,subnodes
!  use fbmaux, only: fbmaux_IsInBoundingBox,fbmaux_PointInHex
!  implicit none
!
!  REAL*8  dcorvg(3,*)
!  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
!  integer :: ibody,nnel
!  integer, dimension(*) :: elements
!
!  ! local variables
!  logical :: found
!  integer :: i,ielem,iconv
!  Real*8, dimension(3) :: center,ref
!
!  found = .false.
!
!  fbm_QueryPointElement_sub = 0
!
!  ! get the center
!  center = myFBM%particleNew(ibody)%Position
!  ref(:) = 0d0
!
!  do i=1,nnel
!
!  ielem = elements(i)
!
!  if(.not.fbmaux_IsInBoundingBox(dcorvg,kvert,ielem,center))cycle
!
!  found = fbmaux_PointInHex(center(1),center(2),center(3),&
!    ref(1),ref(2),ref(3),&
!    ielem,iconv,dcorvg,kvert)
!
!  if(found .eqv. .true.)then
!    fbm_QueryPointElement_sub = ielem
!    return
!  end if
!
!  end do
!
!end function
!!
!!****************************************************************************  
!!  
!subroutine fbm_InitPointLocation(dcorvg,kvert,kedge,karea,nverts)
!  !**@Description
!  ! This function initializes the point location module
!  ! This function should be called when only the point location
!  ! functionality is required
!  !**@EndDescription  
!  use pp3d_mpi, only: myid, mgBoundingBox,updatebvh
!  implicit none
!  REAL*8  dcorvg(3,*)
!  Integer kvert(8,*),kedge(12,*),karea(6,*)
!  integer :: nverts
!
!  ! local variables
!
!  call initpointlocation()
!
!  call setmyid(myid)
!
!  call updatebvh(nverts,dcorvg)
!
!  if(myid.ne.0)then
!    call setdomainbox(mgBoundingBox(myid)%vertices(:,1),&
!      mgBoundingBox(myid)%vertices(:,2))
!  end if
!
!  if(myid.ne.0)then
!    call fbm_BuildUniformGrid(dcorvg,&
!      kvert,kedge,karea)
!  end if
!
!end subroutine fbm_InitPointLocation
!
!****************************************************************************  
!  
end module
