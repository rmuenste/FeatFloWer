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
subroutine Initfield0(MarkerE,kvert,dcorvg,nel,dEps)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
real*8 dEps,dist1,d1,d2,d3,dAreaCrit
real*8 , allocatable :: dArea(:,:)
!
integer iel,jel,i
real*8 dc(3)


MarkerE(1:nel) = 0
allocate(dArea(2,nel))

open(file='Out/area.txt',unit=3)
do iel=1,nel
 dArea(1,iel) = iel
 read(3,*) dArea(2,iel)
end do

call Sortmy2D(dArea(2,:),dArea(1,:),nel)

! do iel=1,nel
!  write(*,*) int(dArea(1,iel)),dArea(2,iel)
! end do
close(3)

do iel=1,nel
 if (dArea(2,iel).lt.1d-8) exit
end do

dAreaCrit = dEps*dArea(2,1)
write(*,*) iel,"/",nel, dAreaCrit

do iel=1,nel
 if (darea(2,iel).gt.dAreaCrit) then
  jel = darea(1,iel)
  MarkerE(jel) = 1
 end if
end do


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
  REAL*8 LW(N),KW(N),LWA,KWA
  INTEGER I,J,N

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

end subroutine Initfield0
