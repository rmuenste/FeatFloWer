subroutine Initfield3(MarkerE,kvert,dcorvg,nel,dEps)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
real*8 dEps
!
integer iel,i
real*8 dc(3),dist1,dist2,dS
integer :: isin, ipc = 0

MarkerE(1:nel) = 0

do iel=1,nel

 dc = 0d0
 do i=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(i,iel))
 end do

!  isin = 0
!  call isinelementid(dc(1),dc(2),dc(3),0,isin)
! 
!  if(isin .gt. 0)then
!    dS = -1d0
!  else
!    dS = +1d0
!  end if
!  
!  call getdistanceid(dc(1),dc(2),dc(3),dist1,ipc)
!  dist1 = dist1*dS
!  if (abs(dist1).lt.dEps) then
!   MarkerE(iel) = 1
!  end if

 dist1 = sqrt((dc(1)+0d0)**2d0 + (dc(2)-0.5d0)**2d0 + (dc(3)-1.5d0)**2d0 )
 dist2 = sqrt((dc(1)+0d0)**2d0 + (dc(2)-0.5d0)**2d0 + (dc(3)-1.5d0)**2d0 )
 
 if (abs(dist1).lt.3.00.and.abs(dist2).gt.2.7) then
  MarkerE(iel) = 1
 end if


end do

end subroutine Initfield3
!
!-----------------------------------------------------------
!
subroutine Initfield2(MarkerE,kvert,dcorvg,nel,dEps)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
real*8 dEps,dist1
!
integer iel,i
real*8 dc(3)

MarkerE(1:nel) = 1

do iel=1,nel

 dc = 0d0
 do i=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(i,iel))
 end do

 if (abs(dc(1)).lt.1.5d0.and.abs(dc(2)).lt.1.5d0.and.abs(dc(3)).lt.1.5d0) then
  MarkerE(iel) = 0
 end if

 if (abs(dc(1)).lt.0.5d0.and.abs(dc(2)).lt.0.5d0.and.abs(dc(3)).lt.0.5d0) then
  MarkerE(iel) = 1
 end if

end do

end subroutine Initfield2
!
!-----------------------------------------------------------
!
subroutine Initfield1(MarkerE,kvert,dcorvg,nel,dEps)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
real*8 dEps,dist1
!
integer iel,i
real*8 dc(3)

MarkerE(1:nel) = 0

do iel=1,nel

 dc = 0d0
 do i=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(i,iel))
 end do

 dist1 = sqrt((dc(1)-2d0)**2d0 + (dc(2)-2.0d0)**2d0 + (dc(3)-2.0d0)**2d0 )
 
 if (dist1.lt.3.6d0) MarkerE(iel) = 1
 
end do

end subroutine Initfield1
!
!-----------------------------------------------------------
!
subroutine Initfield0(MarkerE,kvert,dcorvg,nel)
USE MeshRefVar, only : cIntputFolder,AreaIntensity
USE MeshRefDef, only : GetValueFromFile
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
real*8 dEps,dist1,d1,d2,d3,dAreaCrit,dTotalArea,dCritArea
!
integer iel,jel,i
CHARACTER cInputFile*(256),cVal*(256),cKey*(256)
real*8 dc(3)
real*8 xbox(3),nbox(3),dAreaThreshold,dVolume,dSize,dzmax
integer ii,nnel,ivt

 cInputFile = ADJUSTL(TRIM(cIntputFolder))//'/'//'param.txt'

!  cKey='RefinementFraction'
!  CALL GetValueFromFile(cInputFile,cVal,cKey)
!  read(cVal,*) dEps
! 
 cKey='geometryLength'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) xbox

 cKey='voxelAmount'
 CALL GetValueFromFile(cInputFile,cVal,cKey)
 read(cVal,*) nBox

 dVolume = xbox(1)*xbox(2)*xbox(3)/(nbox(1)*nbox(2)*nbox(3)) ! volume of one voxel
 dSize   = dVolume**(1d0/3d0)
 dAreaThreshold = dVolume**(2d0/3d0)
 WRITE(*,*) dAreaThreshold
 
MarkerE(1:nel) = 0
allocate(AreaIntensity(2,nel))

dzmax = -1d8
do iel=1,nel
 dc = 0d0
 do ivt=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(ivt,iel))
 end do
 if (dzmax.lt.dc(3)) dzmax = dc(3)
end do

open(file=ADJUSTL(TRIM(cIntputFolder))//'/'//'area.txt',unit=3)
do iel=1,nel
 AreaIntensity(1,iel) = iel
 read(3,*) AreaIntensity(2,iel)
end do

close(3)
call Sortmy2D(AreaIntensity(2,:),AreaIntensity(1,:),nel)

call CreateHistogram(AreaIntensity(2,:),nel)

dTotalArea = 0d0
do iel=1,nel
 if (AreaIntensity(2,iel).lt.1d-8) exit
 dTotalArea = dTotalArea + AreaIntensity(2,iel)
end do

! dCritArea = dEps*dTotalArea
! write(*,*) dCritArea,"/",dTotalArea, ' || ', iel,' / ',nel

do iel=1,nel
! if (dTotalArea.gt.dCritArea) exit
 if (AreaIntensity(2,iel).lt.2d0*dAreaThreshold) exit
end do

write(*,*) iel,"/",nel
nnel=iel

ii = 0
do iel=1,nnel
  ii = ii + 1
  jel = AreaIntensity(1,iel)
  MarkerE(jel) = 1
end do

do iel=1,nel
 jel = AreaIntensity(1,iel)
 dc = 0d0
 do ivt=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(ivt,jel))
 end do
 if (dc(3).gt.dzmax-dSize.and.AreaIntensity(2,iel).gt.0d0) then
  MarkerE(jel) = 1
  ii = ii + 1
 end if
end do

WRITE(*,*) 'Number of refined elements: ', ii,dzmax,dSize

!!!! Reload the file for output 
open(file=ADJUSTL(TRIM(cIntputFolder))//'/'//'area.txt',unit=3)

do iel=1,nel
 AreaIntensity(1,iel) = iel
 read(3,*) AreaIntensity(2,iel)
end do

close(3)

! do iel=1,nel
! 
!  dc = 0d0
!  do i=1,8
!   dc = dc + 0.125d0*dcorvg(:,kvert(i,iel))
!  end do
! 
!  dist1 = sqrt((dc(1)-0d0)**2d0 + (dc(2)-0.0d0)**2d0 + (dc(3)-0.0d0)**2d0 )
!  
!  if (dist1.lt.80d0) MarkerE(iel) = 1
!  
! end do
 CONTAINs

SUBROUTINE SORTmy2D(LW,KW,N)
  INTEGER N
  REAL*8 LW(N),KW(N),LWA,KWA
  INTEGER I,J

  DO I=2,N
  DO J=N,I,-1
  IF (LW(J).GT.LW(J-1)) THEN
    LWA     = LW(J)
    KWA     = KW(J)
    LW(J)   = LW(J-1)
    KW(J)   = KW(J-1)
    LW(J-1) = LWA
    KW(J-1) = KWA
  END IF
  END DO
  END DO

END SUBROUTINE SORTmy2D

SUBROUTINE CreateHistogram(LW,N)
 implicit none
  INTEGER N
  REAL*8 LW(N),KW(N),LWA,KWA
  INTEGER I,J,ii
  INTEGER :: nC=16
  real*8, allocatable :: dHist(:),dPivot(:)
  real*8 :: dW, dM, dV,dX,dS,yRef,nRef
  
  dM = LW(1)
  dW = dM/DBLE(nC)
  
  allocate(dHist(nC)); dHist = 0d0
  allocate(dPivot(nC+1)); dPivot = 0d0
  
  dPivot(1) = 0d0
  do j=2,nC+1
   dPivot(j) = dPivot(j-1) + dW
  end do
  dPivot(nC+1) = dPivot(nC+1) + dW

  ii = 0
  do i=1,n
   dV = LW(i)
   do j=1,nC
    if (dV.ge.dPivot(j).and.dV.lt.dPivot(j+1).and.dV.gt.0d0) then
     dHist(j) = dHist(j) + 1
     ii = ii + 1
    end if
   end do
  end do
  
  open(file=ADJUSTL(TRIM(cIntputFolder))//'/'//'hist.txt',unit=4)
  write(4,*) "p h"
  yRef= 0
  nRef= 0
  
  do j=1,nC
   dX = (dPivot(j)+0.5d0*dW)
   if (dX.lt.2d0*dAreaThreshold) then
    nRef = nRef + dHist(j)
    dS = -1d0
   else
    yRef = yRef + dHist(j)
    dS = +1d0
   end if
   write(4,*) dX,dS*(dHist(j)/dble(ii))
  end do
  close(4)
  
  write(*,*) 'nRef/yRef',nRef/dble(ii),yRef/dble(ii)
  
END 

end subroutine Initfield0
