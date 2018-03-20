!##############################################################################
!# ****************************************************************************
!# <name> ptrace </name>
!# ****************************************************************************
!#
!# 
!##############################################################################
MODULE ptrace 

use var_QuadScalar

contains
!
!****************************************************************************  
!
subroutine ptrace_TraceParticles(U,V,W)
!**@Description
!    Starts the particle tracing module
!    
!**End@Description
use fbmaux, only:fbmaux_evalE013_mult
implicit none

Real*8, dimension(:) :: U,V,W

! local variables

 ! query the the uniform grid to find the elements
 ! containing the particles
 call ptrace_QueryPointElement(DWORK(L(LCORVG)),&
               KWORK(L(LVERT)),KWORK(L(LEDGE)),KWORK(L(LAREA)))

 ! evaluate the velocity for each particle in the element
 ! in which it is located
 call fbmaux_evalE013_mult(DWORK(L(LCORVG)),KWORK(L(LVERT)),&
               KWORK(L(LEDGE)),KWORK(L(LAREA)),&
               U,V,W)

 ! The new velocities have been calculated,
 ! to keep the global information up-to-date
 ! we have to synchronize it
 call ptrace_SyncVelocities()

 ! start the particle solver
 call startcollisionpipeline()

 ! get the values from the particle solver
 call FBM_SetParticles()
 
end subroutine ptrace_TraceParticles
!
!****************************************************************************  
!
subroutine ptrace_SyncVelocities()
!**@Description
!    Synchronizes the velocities and positions of all particles
!    between all processes
!**End@Description
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
implicit none

! local variables
integer :: i,errors,ip,ipointer

do ip = 1,myFBM%nParticles
  ! synchronize velocity
  iPointer = 6*(IP-1)
  myFBM%Force(iPointer+1) = myFBM%ParticleNew(IP)%Position(1)
  myFBM%Force(iPointer+2) = myFBM%ParticleNew(IP)%Position(2)
  myFBM%Force(iPointer+3) = myFBM%ParticleNew(IP)%Position(3)
  myFBM%Force(iPointer+4) = myFBM%ParticleNew(iP)%Velocity(1)
  myFBM%Force(iPointer+5) = myFBM%ParticleNew(iP)%Velocity(2)
  myFBM%Force(iPointer+6) = myFBM%ParticleNew(iP)%Velocity(3)
end do ! nParticles

CALL COMM_SUMMN(myFBM%Force,6*myFBM%nParticles)

DO IP = 1,myFBM%nParticles
  iPointer = 6*(IP-1)+1
  myFBM%ParticleNew(IP)%Position(1)= &
  myFBM%Force(iPointer)
  myFBM%ParticleNew(IP)%Position(2) = &
  myFBM%Force(iPointer+1)
  myFBM%ParticleNew(IP)%Position(3) = &
  myFBM%Force(iPointer+2)
  myFBM%ParticleNew(IP)%Velocity(1) = &
  myFBM%Force(iPointer+3)
  myFBM%ParticleNew(IP)%Velocity(2) = &
  myFBM%Force(iPointer+4)
  myFBM%ParticleNew(IP)%Velocity(3) = &
  myFBM%Force(iPointer+5)
END DO

DO IP = 1,myFBM%nParticles
  iPointer = IP-1
  call setvelocityid(myFBM%ParticleNew(ip)%Velocity(1),&
                     myFBM%ParticleNew(ip)%Velocity(2),&
                     myFBM%ParticleNew(ip)%Velocity(3),ipointer)
END DO

end subroutine ptrace_SyncVelocities
!
!****************************************************************************  
!
subroutine ptrace_QueryPointElement(dcorvg,kvert,kedge,karea)
!**@Description
  ! Wrapper function for the PointInElement query
!**@EndDescription  
use pp3d_mpi, only: myid
use var_QuadScalar,only:myFBM 
use pp3d_mpi, only: myid,mgBoundingBox,subnodes
implicit none
REAL*8  dcorvg(3,*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)

! local variables
real*8 vertices(3,8)
integer :: ibody,i,ibodyc,nnel,iel
real*8 point(3)
integer, dimension(:), allocatable :: elements

myFBM%iel_ug=0

do ibody=1,myFBM%nParticles

  ibodyc=ibody-1

  ! clear the list for the body
  call clearelementlists(ibodyc)    

  ! query the uniform grid
  call queryuniformgrid(ibodyc)

  ! get the number of elements
  call getelementsinside(nnel,ibodyc);
  
  ! the body is not in the domain -> skip
  if(nnel .eq. 0)cycle
  
  ! get the element array
  if(allocated(elements))deallocate(elements)
  allocate(elements(nnel))
  call getelements(elements,ibodyc)

  ! start element determination
  iel = ptrace_QueryPointElement_func(dcorvg,kvert,kedge,&
                                      karea,elements,&
                                      nnel,ibody) 

  myFBM%iel_ug(ibody)=iel

end do
    
end subroutine ptrace_QueryPointElement
!
!****************************************************************************  
!
integer function ptrace_QueryPointElement_func(dcorvg,kvert,kedge,karea,elements,nnel,ibody) 
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
      
! local variables
logical :: found
integer :: i,ielem,iconv
Real*8, dimension(3) :: center,ref

found = .false.

ptrace_QueryPointElement_func = 0

! get the center
center = myFBM%particleNew(ibody)%Position
ref(:) = 0d0

do i=1,nnel

  ielem = elements(i)

  if(.not.fbmaux_IsInBoundingBox(dcorvg,kvert,ielem,center))cycle

  found = fbmaux_PointInHex(center(1),center(2),center(3),&
                            ref(1),ref(2),ref(3),&
                            ielem,iconv,dcorvg,kvert)
 
  if(found .eqv. .true.)then
   ptrace_QueryPointElement_func = ielem
   return
  end if

end do

end function
!
!****************************************************************************  
!
end module