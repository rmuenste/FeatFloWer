subroutine Initfield3(MarkerE,kvert,dcorvg,nel)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
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

 isin = 0
 call isinelementid(dc(1),dc(2),dc(3),0,isin)

 if(isin .gt. 0)then
   dS = -1d0
 else
   dS = +1d0
 end if
 
 call getdistanceid(dc(1),dc(2),dc(3),dist1,ipc)
 dist1 = dist1*dS
 if (abs(dist1).lt.0.07d0) then
  MarkerE(iel) = 1
 end if

!  dist1 = sqrt((dc(1)+1.68d0)**2d0 + (dc(2)+1.68d0)**2d0 + (dc(3)-0d0)**2d0 )
!  dist2 = sqrt((dc(1)+1.68d0)**2d0 + (dc(2)+1.68d0)**2d0 + (dc(3)-0d0)**2d0 )
!  
!  if (abs(dist1).lt.3.00.and.abs(dist2).gt.2.7) then
!   MarkerE(iel) = 1
!  end if


end do

end subroutine Initfield3
!
!-----------------------------------------------------------
!
subroutine Initfield2(MarkerE,kvert,dcorvg,nel)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)
!
integer iel,i
real*8 dc(3)

MarkerE(1:nel) = 1

do iel=1,nel

 dc = 0d0
 do i=1,8
  dc = dc + 0.125d0*dcorvg(:,kvert(i,iel))
 end do

 if (abs(dc(1)).gt.0.69d0.or.abs(dc(2)).gt.0.69d0.or.abs(dc(3)).gt.0.69d0) then
  MarkerE(iel) = 0
 end if

 if (abs(dc(1)).lt.0.59d0.and.abs(dc(2)).lt.0.59d0.and.abs(dc(3)).lt.0.59d0) then
  MarkerE(iel) = 0
 end if
 
 if (dc(1).gt.0.59d0) then
  MarkerE(iel) = 0
 end if

end do

end subroutine Initfield2
!
!-----------------------------------------------------------
!
subroutine Initfield1(MarkerE,kvert,dcorvg,nel)
implicit none
integer nel
integer MarkerE(*),kvert(8,*)
real*8 dcorvg(3,*)

! markerE(1) = 1
! markerE(5) = 1
! markerE(25) = 1
! markerE(913) = 1
! markerE(914) = 1
! markerE(543) = 1
! markerE(793) = 1
! markerE(17) = 1

markerE(1043)=1
markerE(1067)=1
markerE(1072)=1
markerE(1074)=1
markerE(1297)=1
markerE(1327)=1
markerE(1340)=1
markerE(1472)=1
markerE(1482)=1
markerE(1558)=1
markerE(1621)=1
markerE(1644)=1
markerE(1757)=1
markerE(1776)=1
markerE(1887)=1
markerE(1892)=1

! markerE(520)=1
! markerE(524)=1
! markerE(542)=1
! markerE(617)=1
! markerE(622)=1
! markerE(638)=1
! markerE(639)=1
! markerE(640)=1
! markerE(644)=1
! markerE(685)=1
! markerE(705)=1
! markerE(721)=1
! markerE(727)=1
! markerE(748)=1
! markerE(762)=1
! markerE(766)=1
! markerE(822)=1
! markerE(836)=1
! markerE(885)=1
! markerE(944)=1
! markerE(977)=1
! markerE(980)=1
! markerE(996)=1
! markerE(1009)=1

! markerE(296)=1
! markerE(334)=1
! markerE(335)=1
! markerE(406)=1
! markerE(592)=1
! markerE(656)=1
! markerE(669)=1
! markerE(670)=1
! markerE(1107)=1
! markerE(1187)=1
! markerE(1711)=1
! markerE(1848)=1
! 
! 
! markerE(1067)=1
! markerE(1096)=1
! markerE(1241)=1
! markerE(1324)=1
! markerE(1396)=1
! markerE(1492)=1
! markerE(1510)=1
! markerE(2361)=1
! markerE(2428)=1
! markerE(2745)=1
! markerE(3083)=1
! markerE(3088)=1
! markerE(3120)=1
! markerE(3542)=1
! markerE(2230)=1
! markerE(2260)=1
! markerE(2290)=1
! markerE(2487)=1
! markerE(2960)=1
! markerE(2991)=1
! markerE(3073)=1
! markerE(3126)=1
! markerE(3180)=1
! markerE(3230)=1
! markerE(3242)=1
! markerE(3439)=1
! markerE(3702)=1
! markerE(3977)=1
! markerE(4001)=1
! markerE(4045)=1

end subroutine Initfield1
