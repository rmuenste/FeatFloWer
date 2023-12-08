Subroutine ALStructExtractorSub(kvert,karea,kedge,dcorvg,&
           LdA,ColA,nvt,net,nat,nel)
USE PP3D_MPI, ONLY:myid
implicit none
real*8 dcorvg(3,*)
integer kvert(8,*),karea(6,*),kedge(12,*)
integer LdA(*),ColA(*)
integer nvt,net,nat,nel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer iel, ivt,jvt, jel, i,ii, iCorner,indice
integer iSet(8,8),xlist(100),nlist

DATA iSet / 1,12, 9,13,21,25,22,27,&
            2, 9,10,14,21,22,23,27,& 
            3,10,11,15,21,23,24,27,&
            4,11,12,16,21,24,25,27,&
            5,20,17,13,26,25,22,27,&
            6,17,18,14,26,22,23,27,&
            7,18,19,15,26,23,24,27,&
            8,19,20,16,26,24,25,27/
 
DO ivt = 1, nvt
 xlist = 0
 nlist = 0
 do ii = LdA(ivt),LdA(ivt+1)-1
  iel = ColA(ii)
  iCorner = -1
  do i=1,8
   if (kvert(i,iel).eq.ivt) then
    iCorner = i
    exit
   end if
  end do
    
  indice = ivt
  call addDOF(indice)
  
  indice = nvt + kedge(iSet(2,iCorner) - 8,iel)
  call addDOF(indice)
  indice = nvt + kedge(iSet(3,iCorner) - 8,iel)
  call addDOF(indice)
  indice = nvt + kedge(iSet(4,iCorner) - 8,iel)
  call addDOF(indice)
  
  indice = nvt + net + karea(iSet(5,iCorner) - 8 - 12,iel)
  call addDOF(indice)
  indice = nvt + net + karea(iSet(6,iCorner) - 8 - 12,iel)
  call addDOF(indice)
  indice = nvt + net + karea(iSet(7,iCorner) - 8 - 12,iel)
  call addDOF(indice)

  indice = nvt + net + nat + iel
  call addDOF(indice)
  
 end do
 
 if (myid.eq.1) write(*,'(3(I0,(" ")),A,100(I0,(" ")))') ivt,nList,(LdA(ivt+1)-LdA(ivt))," : ",xList(1:nList)
end do
pause

 contains

 SUBROUTINE addDOF(iNode)
 integer iNode
 integer i
 
 if (myid.eq.1.and.ivt.eq.23) write(*,*) 'indice =',indice
 
 DO i = 1, nList
  if (xList(i).eq.iNode) RETURN
 END DO
  
 nList = nList + 1 
 xList(nList) = iNode
 
 END SUBROUTINE addDOF

END Subroutine ALStructExtractorSub