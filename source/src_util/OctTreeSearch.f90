Module OcttreeSearch
Implicit none

TYPE tCube
 INTEGER :: n=0
 INTEGER, ALLOCATABLE :: L(:)
END TYPE tCube

TYPE tOctTree
 integer :: nn=0
 integer nx,ny,nz
 real*8 xmax,xmin,ymax,ymin,zmax,zmin,dx,dy,dz,dL
 TYPE(tCube), ALLOCATABLE :: E(:,:,:)
END TYPE tOctTree
TYPE (tOctTree) OctTree 

Logical :: bWrite=.false.

CONTAINS
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE InitOctTree(dcorvg,nvt)
integer :: nvt
real*8 :: dcorvg(3,nvt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer i,i1,iii,iix,iiy,iiz
real*8 dmax,dmin
real*8 p(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 :: dSize

OctTree%xmin = +1d30
OctTree%xmax = -1d30

OctTree%ymin = +1d30
OctTree%ymax = -1d30

OctTree%zmin = +1d30
OctTree%zmax = -1d30

DO i=1,nvt
 if (OctTree%xmin.gt.dcorvg(1,i)) OctTree%xmin=dcorvg(1,i)
 if (OctTree%ymin.gt.dcorvg(2,i)) OctTree%ymin=dcorvg(2,i)
 if (OctTree%zmin.gt.dcorvg(3,i)) OctTree%zmin=dcorvg(3,i)
 
 if (OctTree%xmax.lt.dcorvg(1,i)) OctTree%xmax=dcorvg(1,i)
 if (OctTree%ymax.lt.dcorvg(2,i)) OctTree%ymax=dcorvg(2,i)
 if (OctTree%zmax.lt.dcorvg(3,i)) OctTree%zmax=dcorvg(3,i)
end do

OctTree%dx = OctTree%xmax - OctTree%xmin
OctTree%dy = OctTree%ymax - OctTree%ymin
OctTree%dz = OctTree%zmax - OctTree%zmin

dSize = SQRT(OctTree%dx**2d0 + OctTree%dy**2d0 + OctTree%dz**2d0)*0.025d0

OctTree%xmin = OctTree%xmin - 1e-3*OctTree%dx      
OctTree%xmax = OctTree%xmax + 1e-3*OctTree%dx

OctTree%ymin = OctTree%ymin - 1e-3*OctTree%dy
OctTree%ymax = OctTree%ymax + 1e-3*OctTree%dy

OctTree%zmin = OctTree%zmin - 1e-3*OctTree%dz
OctTree%zmax = OctTree%zmax + 1e-3*OctTree%dz
               
OctTree%dx = OctTree%xmax - OctTree%xmin
OctTree%dy = OctTree%ymax - OctTree%ymin
OctTree%dz = OctTree%zmax - OctTree%zmin

if (bWrite) then
write(*,*) 'OctTree%xmin,OctTree%xmax: ',OctTree%xmin,OctTree%xmax,OctTree%dx
write(*,*) 'OctTree%ymin,OctTree%ymax: ',OctTree%ymin,OctTree%ymax,OctTree%dy
write(*,*) 'OctTree%zmin,OctTree%zmax: ',OctTree%zmin,OctTree%zmax,OctTree%dz
end if

if (OctTree%dx.ge.OctTree%dy.and.OctTree%dx.ge.OctTree%dz) OctTree%dL = OctTree%dx
if (OctTree%dy.ge.OctTree%dx.and.OctTree%dy.ge.OctTree%dz) OctTree%dL = OctTree%dy
if (OctTree%dz.ge.OctTree%dx.and.OctTree%dz.ge.OctTree%dy) OctTree%dL = OctTree%dz

OctTree%nn = MAX(1,NINT(OctTree%dL / dSize))
if (bWrite) WRITE(*,*) 'max cube size:', dSize, OctTree%nn

OctTree%nx = max(1,nint(DBLE(OctTree%nn)*OctTree%dx/OctTree%dL))
OctTree%ny = max(1,nint(DBLE(OctTree%nn)*OctTree%dy/OctTree%dL))
OctTree%nz = max(1,nint(DBLE(OctTree%nn)*OctTree%dz/OctTree%dL))
if (bWrite) write(*,*) 'Limiting size: ',OctTree%dL,OctTree%nx,OctTree%ny,OctTree%nz

ALLOCATE(OctTree%E(0:OctTree%nx+1,0:OctTree%ny+1,0:OctTree%nz+1))

DO i=1,nvt

 p = dcorvg(:,i)
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
end do

do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   ALLOCATE(OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n))
   OctTree%E(iix,iiy,iiz)%n = 0
  end do
 end do
end do


DO i=1,nvt

 p = dcorvg(:,i)
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
 OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n) = i
 
end do

iii = 0
do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   iii = iii + OctTree%E(iix,iiy,iiz)%n
  end do
 end do
end do

!if (bWrite) write(*,'(A,2I)') 'Check if all entities were assigned to OctTree Structure nn(1,2) = ', iii, nvt
if (iii.ne.nvt) then
 write(*,'(A,2I0)') 'problem with initialization of the OctTree Structure nn(1,2) = ', iii, nvt
 stop
end if

END SUBROUTINE InitOctTree
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE InitOctTreeC(dcorvgX,dcorvgY,dcorvgZ,nvt)
integer :: nvt
real*8 :: dcorvgX(nvt),dcorvgY(nvt),dcorvgZ(nvt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer i,i1,iii,iix,iiy,iiz
real*8 dmax,dmin
real*8 p(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 :: dSize

OctTree%xmin = +1d30
OctTree%xmax = -1d30

OctTree%ymin = +1d30
OctTree%ymax = -1d30

OctTree%zmin = +1d30
OctTree%zmax = -1d30

DO i=1,nvt
 if (OctTree%xmin.gt.dcorvgX(i)) OctTree%xmin=dcorvgX(i)
 if (OctTree%ymin.gt.dcorvgY(i)) OctTree%ymin=dcorvgY(i)
 if (OctTree%zmin.gt.dcorvgZ(i)) OctTree%zmin=dcorvgZ(i)
 
 if (OctTree%xmax.lt.dcorvgX(i)) OctTree%xmax=dcorvgX(i)
 if (OctTree%ymax.lt.dcorvgY(i)) OctTree%ymax=dcorvgY(i)
 if (OctTree%zmax.lt.dcorvgZ(i)) OctTree%zmax=dcorvgZ(i)
end do

OctTree%dx = OctTree%xmax - OctTree%xmin
OctTree%dy = OctTree%ymax - OctTree%ymin
OctTree%dz = OctTree%zmax - OctTree%zmin

dSize = SQRT(OctTree%dx**2d0 + OctTree%dy**2d0 + OctTree%dz**2d0)*0.025d0

OctTree%xmin = OctTree%xmin - 1e-3*OctTree%dx      
OctTree%xmax = OctTree%xmax + 1e-3*OctTree%dx

OctTree%ymin = OctTree%ymin - 1e-3*OctTree%dy
OctTree%ymax = OctTree%ymax + 1e-3*OctTree%dy

OctTree%zmin = OctTree%zmin - 1e-3*OctTree%dz
OctTree%zmax = OctTree%zmax + 1e-3*OctTree%dz
               
OctTree%dx = OctTree%xmax - OctTree%xmin
OctTree%dy = OctTree%ymax - OctTree%ymin
OctTree%dz = OctTree%zmax - OctTree%zmin

if (bWrite) then
write(*,*) 'OctTree%xmin,OctTree%xmax: ',OctTree%xmin,OctTree%xmax,OctTree%dx
write(*,*) 'OctTree%ymin,OctTree%ymax: ',OctTree%ymin,OctTree%ymax,OctTree%dy
write(*,*) 'OctTree%zmin,OctTree%zmax: ',OctTree%zmin,OctTree%zmax,OctTree%dz
end if

if (OctTree%dx.ge.OctTree%dy.and.OctTree%dx.ge.OctTree%dz) OctTree%dL = OctTree%dx
if (OctTree%dy.ge.OctTree%dx.and.OctTree%dy.ge.OctTree%dz) OctTree%dL = OctTree%dy
if (OctTree%dz.ge.OctTree%dx.and.OctTree%dz.ge.OctTree%dy) OctTree%dL = OctTree%dz

OctTree%nn = MAX(1,NINT(OctTree%dL / dSize))
if (bWrite) WRITE(*,*) 'max cube size:', dSize, OctTree%nn

OctTree%nx = max(1,nint(DBLE(OctTree%nn)*OctTree%dx/OctTree%dL))
OctTree%ny = max(1,nint(DBLE(OctTree%nn)*OctTree%dy/OctTree%dL))
OctTree%nz = max(1,nint(DBLE(OctTree%nn)*OctTree%dz/OctTree%dL))
if (bWrite) write(*,*) 'Limiting size: ',OctTree%dL,OctTree%nx,OctTree%ny,OctTree%nz

ALLOCATE(OctTree%E(0:OctTree%nx+1,0:OctTree%ny+1,0:OctTree%nz+1))

DO i=1,nvt

 p = [dcorvgX(i),dcorvgY(i),dcorvgZ(i)]
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
end do

do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   ALLOCATE(OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n))
   OctTree%E(iix,iiy,iiz)%n = 0
  end do
 end do
end do


DO i=1,nvt

 p = [dcorvgX(i),dcorvgY(i),dcorvgZ(i)]
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 OctTree%E(iix,iiy,iiz)%n = OctTree%E(iix,iiy,iiz)%n + 1
 OctTree%E(iix,iiy,iiz)%L(OctTree%E(iix,iiy,iiz)%n) = i
 
end do

iii = 0
do iix=1,OctTree%nx
 do iiy=1,OctTree%ny
  do iiz=1,OctTree%nz
   iii = iii + OctTree%E(iix,iiy,iiz)%n
  end do
 end do
end do

!if (bWrite) write(*,'(A,2I)') 'Check if all entities were assigned to OctTree Structure nn(1,2) = ', iii, nvt
if (iii.ne.nvt) then
 write(*,'(A,2I0)') 'problem with initialization of the OctTree Structure nn(1,2) = ', iii, nvt
 stop
end if

END SUBROUTINE InitOctTreeC
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FindInOctTreeC(dcorvgX,dcorvgY,dcorvgZ,nvt,p,iP,dist)
integer :: nvt,iP
real*8 :: dcorvgX(nvt),dcorvgY(nvt),dcorvgZ(nvt),p(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer i,i1,iii,iix,iiy,iiz,iiix,iiiy,iiiz
real*8 dmax,dmin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 :: q(3),dist,daux
integer k
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 if ((iiX.lt.1).or.(iiX.gt.OctTree%nx).or.&
     (iiY.lt.1).or.(iiY.gt.OctTree%ny).or.&
     (iiZ.lt.1).or.(iiZ.gt.OctTree%nz)) then
      dist = 1d30
      iP = -1
      return
!      write(*,*) 'fatal problem in OctTreeSearch'
!      pause
 end if

 dist = 1d30
 iP = -1
 
 do iiix=iix-1,iix+1
  do iiiy=iiy-1,iiy+1
   do iiiz=iiz-1,iiz+1
    do iii = 1,OctTree%E(iiix,iiiy,iiiz)%N
     k = OctTree%E(iiix,iiiy,iiiz)%L(iii)
     q = [dcorvgX(k),dcorvgY(k),dcorvgZ(k)]
     daux = sqrt((p(1)-q(1))**2d0 + (p(2)-q(2))**2d0 + (p(3)-q(3))**2d0)
     IF (daux.lt.dist) then
      dist = daux
      iP = k
     END IF
!      bFound = .true.
!      GOTO 1
    end do
   end do
  end do
 end do

END SUBROUTINE FindInOctTreeC
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FindInOctTree(dcorvg,nvt,p,iP,dist)
integer :: nvt,iP
real*8 :: dcorvg(3,nvt),p(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer i,i1,iii,iix,iiy,iiz,iiix,iiiy,iiiz
real*8 dmax,dmin
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 :: q(3),dist,daux
integer k
 
 iiX = 0
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 iiY = 0
 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 iiZ = 0
 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 if ((iiX.lt.1).or.(iiX.gt.OctTree%nx).or.&
     (iiY.lt.1).or.(iiY.gt.OctTree%ny).or.&
     (iiZ.lt.1).or.(iiZ.gt.OctTree%nz)) then
      dist = 1d30
      iP = -1
!       return
     write(*,*) 'fatal problem in OctTreeSearch'
     pause
 end if

 dist = 1d30
 iP = -1
 
 do iiix=iix-1,iix+1
  do iiiy=iiy-1,iiy+1
   do iiiz=iiz-1,iiz+1
    do iii = 1,OctTree%E(iiix,iiiy,iiiz)%N
     k = OctTree%E(iiix,iiiy,iiiz)%L(iii)
     q = dcorvg(:,k)
     daux = sqrt((p(1)-q(1))**2d0 + (p(2)-q(2))**2d0 + (p(3)-q(3))**2d0)
     IF (daux.lt.dist) then
      dist = daux
      iP = k
     END IF
!      bFound = .true.
!      GOTO 1
    end do
   end do
  end do
 end do

END SUBROUTINE FindInOctTree
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FindInPeriodicOctTree(dcorvg,nvt,pp,iP,dist,dPerio)
integer :: nvt,iP
real*8 :: dcorvg(3,nvt),pp(3)
real*8, intent(in), optional :: dPerio(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer i,i1,iii,iix,iiy,iiz,iiix,iiiy,iiiz
real*8 dmax,dmin,p(3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL*8 :: q(3),dist,daux
REAL*8 :: dPSign(2)=[-1d0,1d0]
integer k,iDirection,iPerio
 
 dist = 1d30
 iP = -1
 
 iF (.not.present(dPerio)) return
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 DO iDirection = 1,3
 
 DO iPerio=1,2
 
 p = pp
 
 p(iDirection) = p(iDirection) + dPSign(iPerio)*dPerio(iDirection)
 
 if ((p(1).ge. OctTree%xmin.and.p(1).le. OctTree%xmax).and.&
     (p(2).ge. OctTree%ymin.and.p(2).le. OctTree%ymax).and.&
     (p(3).ge. OctTree%zmin.and.p(3).le. OctTree%zmax)) then
 
 do i1=1,OctTree%nx
  dmin = OctTree%xmin + dble(i1-1)*OctTree%dx/dble(OctTree%nx)
  dmax = OctTree%xmin + dble(i1-0)*OctTree%dx/dble(OctTree%nx)
  if (p(1).ge.dmin.and.p(1).le.dmax) then
   iiX = i1
  end if  
 end do

 do i1=1,OctTree%ny
  dmin = OctTree%ymin + dble(i1-1)*OctTree%dy/dble(OctTree%ny)
  dmax = OctTree%ymin + dble(i1-0)*OctTree%dy/dble(OctTree%ny)
  if (p(2).ge.dmin.and.p(2).le.dmax) then
   iiY = i1
  end if  
 end do

 do i1=1,OctTree%nz
  dmin = OctTree%zmin + dble(i1-1)*OctTree%dz/dble(OctTree%nz)
  dmax = OctTree%zmin + dble(i1-0)*OctTree%dz/dble(OctTree%nz)
  if (p(3).ge.dmin.and.p(3).le.dmax) then
   iiZ = i1
  end if  
 end do

 if ((iiX.lt.1).or.(iiX.gt.OctTree%nx).or.&
     (iiY.lt.1).or.(iiY.gt.OctTree%ny).or.&
     (iiZ.lt.1).or.(iiZ.gt.OctTree%nz)) then
      dist = 1d30
      iP = -1
      return
!      write(*,*) 'fatal problem in OctTreeSearch'
!      pause
 end if

 do iiix=iix-1,iix+1
  do iiiy=iiy-1,iiy+1
   do iiiz=iiz-1,iiz+1
    do iii = 1,OctTree%E(iiix,iiiy,iiiz)%N
     k = OctTree%E(iiix,iiiy,iiiz)%L(iii)
     q = dcorvg(:,k)
     daux = sqrt((p(1)-q(1))**2d0 + (p(2)-q(2))**2d0 + (p(3)-q(3))**2d0)
     IF (daux.lt.dist) then
      dist = daux
      iP = k
!       write(*,*) dist
!       pause
     END IF
!      bFound = .true.
!      GOTO 1
    end do
   end do
  end do
 end do
 
 END IF !point is in octtree structures .... 
 
 END DO
 
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE FindInPeriodicOctTree
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FreeOctTree

deallocate(OctTree%E)
OctTree%nn = 0

END SUBROUTINE FreeOctTree
!
!-------------------------------------------------------------------------------------
!

End Module OcttreeSearch
