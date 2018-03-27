!##############################################################################
!# ****************************************************************************
!# <name> fbmaux </name>
!# ****************************************************************************
!#
!# 
!##############################################################################
MODULE fbmaux 

use var_QuadScalar

contains
!
!****************************************************************************  
!
subroutine fbmaux_getCoordinates(kvert,karea,kedge,dcorvg,dpointsCub,ele)
!**@Description
!   Function returns the cubature points in real coordinates
!   for each element
!**End@Description
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
            NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8  DCORVG(NNDIM,*),DVAL(3),dpointsCub(3,*)

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks
INTEGER  VIPARM 
DIMENSION VIPARM(100)
EQUIVALENCE (IAUSAV,VIPARM)
COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
              IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
              ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
              INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
              ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
              IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

EXTERNAL ELE

! local variables
integer iel1,i,j,ive,ieltyp,jp,icub,icubp,iel,ipoints
real*8 dxpos,dypos,dzpos
real*8 x,y,z

IF (myid.ne.0)then

ipoints=0

do i = 1,nnder
  bder(i)=.false.
enddo

! we only need the function
! we do not need jacobians
bder(1)=.true.

! Ask the element about its type:
ieltyp=-1
call ele(0d0,0d0,0d0,ieltyp)

! how many basis functions=DOF's do we have on each element?
idfl=ndfl(ieltyp)

! select the cubature formula
icub=9

! initialize the cubature formula
call cb3h(icub)
if (ier.lt.0) return

icubp=icub
CALL ELE (0d0,0d0,0d0,-2)

! loop over all elements
do iel=1,nel

  ! Calculate the local and global DOF's on our current element.
  ! We later have to loop about them...
  call ndfgl(iel,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)
  if (ier.lt.0) return

  ! Evaluation of vertex coordinates
  do ive=1,nve
    jp=kvert(ive,iel)
    kve(ive)=jp
    dx(ive)=dcorvg(1,jp)
    dy(ive)=dcorvg(2,jp)
    dz(ive)=dcorvg(3,jp)
  end do

  ! Compute coefficients of the jacobian
  dj11=( dx(1)+dx(2)+dx(3)+dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
  dj12=( dy(1)+dy(2)+dy(3)+dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
  dj13=( dz(1)+dz(2)+dz(3)+dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
  dj21=(-dx(1)+dx(2)+dx(3)-dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
  dj22=(-dy(1)+dy(2)+dy(3)-dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
  dj23=(-dz(1)+dz(2)+dz(3)-dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
  dj31=(-dx(1)-dx(2)+dx(3)+dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
  dj32=(-dy(1)-dy(2)+dy(3)+dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
  dj33=(-dz(1)-dz(2)+dz(3)+dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
  dj41=(-dx(1)-dx(2)-dx(3)-dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
  dj42=(-dy(1)-dy(2)-dy(3)-dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
  dj43=(-dz(1)-dz(2)-dz(3)-dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
  dj51=( dx(1)-dx(2)+dx(3)-dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
  dj52=( dy(1)-dy(2)+dy(3)-dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
  dj53=( dz(1)-dz(2)+dz(3)-dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
  dj61=( dx(1)-dx(2)-dx(3)+dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
  dj62=( dy(1)-dy(2)-dy(3)+dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
  dj63=( dz(1)-dz(2)-dz(3)+dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
  dj71=( dx(1)+dx(2)-dx(3)-dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
  dj72=( dy(1)+dy(2)-dy(3)-dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
  dj73=( dz(1)+dz(2)-dz(3)-dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
  dj81=(-dx(1)+dx(2)-dx(3)+dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
  dj82=(-dy(1)+dy(2)-dy(3)+dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
  dj83=(-dz(1)+dz(2)-dz(3)+dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8  
  
  ! loop over the cubature points
  do icubp=1,ncubp
    
    xi1=dxi(icubp,1)
    xi2=dxi(icubp,2)
    xi3=dxi(icubp,3)
  
    ! build the jacobian of bilinear mapping onto the reference element 
    djac(1,1)=dj21+dj51*xi2+dj61*xi3+dj81*xi2*xi3
    djac(1,2)=dj31+dj51*xi1+dj71*xi3+dj81*xi1*xi3
    djac(1,3)=dj41+dj61*xi1+dj71*xi2+dj81*xi1*xi2
    djac(2,1)=dj22+dj52*xi2+dj62*xi3+dj82*xi2*xi3
    djac(2,2)=dj32+dj52*xi1+dj72*xi3+dj82*xi1*xi3
    djac(2,3)=dj42+dj62*xi1+dj72*xi2+dj82*xi1*xi2
    djac(3,1)=dj23+dj53*xi2+dj63*xi3+dj83*xi2*xi3
    djac(3,2)=dj33+dj53*xi1+dj73*xi3+dj83*xi1*xi3
    djac(3,3)=dj43+dj63*xi1+dj73*xi2+dj83*xi1*xi2
    detj= djac(1,1)*(djac(2,2)*djac(3,3)-djac(3,2)*djac(2,3))&
        -djac(2,1)*(djac(1,2)*djac(3,3)-djac(3,2)*djac(1,3))&
        +djac(3,1)*(djac(1,2)*djac(2,3)-djac(2,2)*djac(1,3))
     
    ! Computation of the cartesian coordiante of the cubature point
    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2

    ipoints=ipoints+1        
    dpointsCub(1,ipoints)=XX
    dpointsCub(2,ipoints)=YY
    dpointsCub(3,ipoints)=ZZ  
    

  end do ! icubp
end do ! iel

write(*,*)'Number of cub points: ',ipoints

end if

end subroutine fbmaux_getCoordinates
!
!****************************************************************************  
!
subroutine fbmaux_L2Project(kvert,karea,kedge,dcorvg,dfield,rhs,ele)
!**@Description
!   Function returns the cubature points in real coordinates
!   for each element
!**End@Description
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
            NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8  DCORVG(NNDIM,*),DVAL(3),dfield(3,*),rhs(3,*)

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks
INTEGER  VIPARM 
DIMENSION VIPARM(100)
EQUIVALENCE (IAUSAV,VIPARM)
COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
              IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
              ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
              INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
              ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
              IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

EXTERNAL ELE

! local variables
integer iel1,i,j,ive,ieltyp,jp,icub,icubp,iel,ipoints
real*8 dxpos,dypos,dzpos
real*8 x,y,z

IF (myid.ne.0)then

ipoints=0

do i = 1,nnder
  bder(i)=.false.
enddo

! we only need the function
! we do not need jacobians
bder(1)=.true.

! Ask the element about its type:
ieltyp=-1
call ele(0d0,0d0,0d0,ieltyp)

! how many basis functions=DOF's do we have on each element?
idfl=ndfl(ieltyp)

! select the cubature formula
icub=9

! initialize the cubature formula
call cb3h(icub)
if (ier.lt.0) return

icubp=icub
CALL ELE (0d0,0d0,0d0,-2)

! loop over all elements
do iel=1,nel

  ! Calculate the local and global DOF's on our current element.
  ! We later have to loop about them...
  call ndfgl(iel,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)
  if (ier.lt.0) return

  ! Evaluation of vertex coordinates
  do ive=1,nve
    jp=kvert(ive,iel)
    kve(ive)=jp
    dx(ive)=dcorvg(1,jp)
    dy(ive)=dcorvg(2,jp)
    dz(ive)=dcorvg(3,jp)
  end do

  ! Compute coefficients of the jacobian
  dj11=( dx(1)+dx(2)+dx(3)+dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
  dj12=( dy(1)+dy(2)+dy(3)+dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
  dj13=( dz(1)+dz(2)+dz(3)+dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
  dj21=(-dx(1)+dx(2)+dx(3)-dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
  dj22=(-dy(1)+dy(2)+dy(3)-dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
  dj23=(-dz(1)+dz(2)+dz(3)-dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
  dj31=(-dx(1)-dx(2)+dx(3)+dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
  dj32=(-dy(1)-dy(2)+dy(3)+dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
  dj33=(-dz(1)-dz(2)+dz(3)+dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
  dj41=(-dx(1)-dx(2)-dx(3)-dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
  dj42=(-dy(1)-dy(2)-dy(3)-dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
  dj43=(-dz(1)-dz(2)-dz(3)-dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
  dj51=( dx(1)-dx(2)+dx(3)-dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
  dj52=( dy(1)-dy(2)+dy(3)-dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
  dj53=( dz(1)-dz(2)+dz(3)-dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
  dj61=( dx(1)-dx(2)-dx(3)+dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
  dj62=( dy(1)-dy(2)-dy(3)+dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
  dj63=( dz(1)-dz(2)-dz(3)+dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
  dj71=( dx(1)+dx(2)-dx(3)-dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
  dj72=( dy(1)+dy(2)-dy(3)-dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
  dj73=( dz(1)+dz(2)-dz(3)-dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
  dj81=(-dx(1)+dx(2)-dx(3)+dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
  dj82=(-dy(1)+dy(2)-dy(3)+dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
  dj83=(-dz(1)+dz(2)-dz(3)+dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8  
  
  ! loop over the cubature points
  do icubp=1,ncubp
    
    xi1=dxi(icubp,1)
    xi2=dxi(icubp,2)
    xi3=dxi(icubp,3)
  
    ! build the jacobian of bilinear mapping onto the reference element 
    djac(1,1)=dj21+dj51*xi2+dj61*xi3+dj81*xi2*xi3
    djac(1,2)=dj31+dj51*xi1+dj71*xi3+dj81*xi1*xi3
    djac(1,3)=dj41+dj61*xi1+dj71*xi2+dj81*xi1*xi2
    djac(2,1)=dj22+dj52*xi2+dj62*xi3+dj82*xi2*xi3
    djac(2,2)=dj32+dj52*xi1+dj72*xi3+dj82*xi1*xi3
    djac(2,3)=dj42+dj62*xi1+dj72*xi2+dj82*xi1*xi2
    djac(3,1)=dj23+dj53*xi2+dj63*xi3+dj83*xi2*xi3
    djac(3,2)=dj33+dj53*xi1+dj73*xi3+dj83*xi1*xi3
    djac(3,3)=dj43+dj63*xi1+dj73*xi2+dj83*xi1*xi2
    detj= djac(1,1)*(djac(2,2)*djac(3,3)-djac(3,2)*djac(2,3))&
        -djac(2,1)*(djac(1,2)*djac(3,3)-djac(3,2)*djac(1,3))&
        +djac(3,1)*(djac(1,2)*djac(2,3)-djac(2,2)*djac(1,3))
     
    ipoints=ipoints+1
    
    XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
    YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
    ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2 
    
!     if(myid.eq.1)then
!     write(*,*)'',abs(XX-dfield(1,ipoints))
!     end if
     
    OM=DOMEGA(ICUBP)*ABS(DETJ)

    CALL ELE(XI1,XI2,XI3,-3)
    IF (IER.LT.0) return

! *** Summing up over all pairs of multiindices
      DO JDFL=1,IDFL
        IG=KDFG(JDFL)
        IL=KDFL(JDFL)
        HBAS=DBAS(1,IL,1)
        rhs(1,IG)=rhs(1,IG) + OM*1d0*HBAS
        rhs(2,IG)=rhs(2,IG) + OM*dfield(2,ipoints)*HBAS
        rhs(3,IG)=rhs(3,IG) + OM*dfield(3,ipoints)*HBAS
     end do
    
  end do ! icubp
end do ! iel

!write(*,*)'Number of cub points: ',ipoints

end if

end subroutine fbmaux_L2Project
!
!****************************************************************************  
!
subroutine fbmaux_evalE013_mult(dcorvg,kvert,kedge,karea,u,v,w)
!**@Description
!  Wrapper function that starts the evaluatation
!  for all particles
!**End@Description
use var_QuadScalar,only:myFBM 
implicit none
REAL*8  dcorvg(3,*),u(*),v(*),w(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)

! local variables
real*8, dimension(3) :: pos
real*8, dimension(3) :: dvalues
logical :: bRes

external e013

integer :: i,j,k,ielem

do i=1,myFBM%nParticles

  dvalues(:)=0d0
  pos(:)=myFBM%ParticleNew(i)%Position(:)
  ielem = myFBM%iel_ug(i)
  myFBM%ParticleNew(i)%Velocity(:)=0d0

  if(ielem.ne.0)then
    ! we have found the element, now we can evaluate
    call fbmaux_evalVelAux(pos(1),pos(2),pos(3),dvalues,ielem,&
                               kvert,kedge,karea,dcorvg,&
                               u,v,w,&
                               e013)

    myFBM%ParticleNew(i)%Velocity(:)=dvalues(:)
  end if

end do

end subroutine fbmaux_evalE013_mult
!
!****************************************************************************  
!
subroutine fbmaux_evalE013_sing(dcorvg,kvert,kedge,karea,u,v,w,dpos,dvalues,ielem)
!**@Description
!  Wrapper function that starts the evaluatation
!  for a certain point in a certain element
!**End@Description
use var_QuadScalar,only:myFBM 
implicit none
REAL*8  dcorvg(3,*),u(*),v(*),w(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ielem
real*8, dimension(3) :: dpos
real*8, dimension(3) :: dvalues
logical :: bRes

! local variables
external e013

dvalues(:)=0d0

! we have found the element, now we can evaluate
call fbmaux_evalVelAux(dpos(1),dpos(2),dpos(3),dvalues,ielem,&
                            kvert,kedge,karea,dcorvg,&
                            u,v,w,&
                            e013)

end subroutine fbmaux_evalE013_sing
!
!****************************************************************************  
!
subroutine fbmaux_evalVelAux(dxpos,dypos,dzpos,dval,ielem,&
                           kvert,kedge,karea,dcorvg,&
                           u,v,w,&
                           ele)
!**@Description
!    Evaluates the velocity in the cell in which
!    the particle is located
!**End@Description
USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myFBM,myExport
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
            NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8  DCORVG(NNDIM,*),DVAL(3),u(*),v(*),w(*)

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks
INTEGER  VIPARM 
DIMENSION VIPARM(100)
EQUIVALENCE (IAUSAV,VIPARM)
COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
              IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
              ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
              INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
              ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
              IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

EXTERNAL ELE

! local variables
integer iel1,i,j,ive,ieltyp,jp,icub
real*8 dxpos,dypos,dzpos
real*8 x,y,z

IF (myid.ne.0)then

! initialize the result array
dval(:) = 0.0d0

! set the element
IEL1=IELEM

! set the evaluation point
X = DXPOS
Y = DYPOS
Z = DZPOS

do i = 1,nnder
  bder(i)=.false.
enddo

! we only need the function
! we do not need jacobians
bder(1)=.true.

! Ask the element about its type:
ieltyp=-1
call ele(0d0,0d0,0d0,ieltyp)

! how many basis functions=DOF's do we have on each element?
idfl=ndfl(ieltyp)

icub=9

! initialize the cubature formula
call cb3h(icub)
if (ier.lt.0) return

! Calculate the local and global DOF's on our current element.
! We later have to loop about them...
call ndfgl(iel1,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)

if (ier.lt.0) return

! Initialise the COMMON block of the element to inform the it
! about our current element
do ive = 1, nve
  jp=kvert(ive,iel1)
  kve(ive)=jp
  dx(ive)=dcorvg(1,jp)
  dy(ive)=dcorvg(2,jp)
  dz(ive)=dcorvg(3,jp)
end do

! Call the element routine to calculate the weights of the basis 
! functions in the current point:
CALL ELE (X,Y,Z,0)

! Is there at all something to do, or is the solution vector empty?!?
if (ier.lt.0) return

! Loop about the local DOF's and sum up the values there together with
! the just calculated weights to obtain the function value in (X,Y,Z):
DO I=1,IDFL
  IG=KDFG(I)
  IL=KDFL(I)
  
  ! multiply the basis function by the FE-coefficients and sum up
  dval(1) = dval(1) + u(IG)  * DBAS(1,IL,1)
  dval(2) = dval(2) + v(IG)  * DBAS(1,IL,1)
  dval(3) = dval(3) + w(IG)  * DBAS(1,IL,1)
!  write(*,'(A,4E16.4)')'global sol: ',u(IG),v(IG),w(IG),dval(3)
END DO

end if

end subroutine fbmaux_evalVelAux
!
!****************************************************************************  
!
subroutine fbmaux_evalPhi_sing(dcorvg,kvert,kedge,karea,u,v,w,dpos,dvalues,ielem)
!**@Description
!  Wrapper function that starts the evaluatation
!  for a certain point in a certain element
!**End@Description
use var_QuadScalar,only:myFBM 
implicit none
REAL*8  dcorvg(3,*),u(*),v(*),w(*)
INTEGER kvert(8,*),kedge(12,*),karea(6,*)
integer :: ielem
real*8, dimension(3) :: dpos
real*8, dimension(5) :: dvalues
logical :: bRes

! local variables
external e013

dvalues(:)=0d0

! we have found the element, now we can evaluate
call fbmaux_evalPhiAux(dpos(1),dpos(2),dpos(3),dvalues,ielem,&
                            kvert,kedge,karea,dcorvg,&
                            u,v,w,&
                            e013)

end subroutine fbmaux_evalPhi_sing
!
!****************************************************************************  
!
subroutine fbmaux_evalPhiAux(dxpos,dypos,dzpos,dval,ielem,&
                           kvert,kedge,karea,dcorvg,&
                           ux,uy,uz,&
                           ele)

USE PP3D_MPI, ONLY:myid,showID,COMM_SUMMN
USE var_QuadScalar, ONLY : myFBM !,myExport
IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)

PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,&
            NNDIM=3,NNCOF=10)
PARAMETER (Q2=0.5D0,Q8=0.125D0)
!
INTEGER KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
INTEGER KDFG(NNBAS),KDFL(NNBAS)
REAL*8  DCORVG(NNDIM,*),DVAL(5),ux(*),uy(*),uz(*)

COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
COMMON /ERRCTL/ IER,ICHECK
COMMON /CHAR/   SUB,FMT(3),CPARAM
COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,&
                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),&
                IEL,NDIM
COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,&
                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
COMMON /COAUX1/ KDFG,KDFL,IDFL
!
! *** user COMMON blocks
INTEGER  VIPARM 
DIMENSION VIPARM(100)
EQUIVALENCE (IAUSAV,VIPARM)
COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,&
              IMASS,IMASSL,IUPW,IPRECA,IPRECB,&
              ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,&
              INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,&
              ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,&
              IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA

EXTERNAL ELE

! local variables
integer iel1,i,j,ive,ieltyp,jp,icub
real*8 dxpos,dypos,dzpos
real*8 x,y,z

IF (myid.ne.0)then

! initialize the result array
dval(:) = 0.0d0

! set the element
IEL1=IELEM

! set the evaluation point
X = DXPOS
Y = DYPOS
Z = DZPOS

do i = 1,nnder
  bder(i)=.false.
enddo

! we only need the function
! we do not need jacobians
bder(1)=.true.

! Ask the element about its type:
ieltyp=-1
call ele(0d0,0d0,0d0,ieltyp)

! how many basis functions=DOF's do we have on each element?
idfl=ndfl(ieltyp)

icub=9

! initialize the cubature formula
call cb3h(icub)
if (ier.lt.0) return

! Calculate the local and global DOF's on our current element.
! We later have to loop about them...
call ndfgl(iel1,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)

if (ier.lt.0) return

! Initialise the COMMON block of the element to inform the it
! about our current element
do ive = 1, nve
  jp=kvert(ive,iel1)
  kve(ive)=jp
  dx(ive)=dcorvg(1,jp)
  dy(ive)=dcorvg(2,jp)
  dz(ive)=dcorvg(3,jp)
end do

! Call the element routine to calculate the weights of the basis 
! functions in the current point:
CALL ELE (X,Y,Z,0)

! Is there at all something to do, or is the solution vector empty?!?
if (ier.lt.0) return

! Loop about the local DOF's and sum up the values there together with
! the just calculated weights to obtain the function value in (X,Y,Z):
DO I=1,IDFL
  IG=KDFG(I)
  IL=KDFL(I)
  
  ! multiply the basis function by the FE-coefficients and sum up
  dval(1) = dval(1) + ux(IG)     * DBAS(1,IL,1)
  dval(2) = dval(2) + uy(IG)     * DBAS(1,IL,1)
  dval(3) = dval(3) + uz(IG)     * DBAS(1,IL,1)
!  dval(4) = dval(4) + GDefInfo%monitorFunc(IG)      * DBAS(1,IL,1)
!  dval(5) = dval(5) + GDefInfo%volumesAtVertex(IG)  * DBAS(1,IL,1)
END DO

dval(4) = 1.0d0/dval(4)
dval(5) = 1.0d0/dval(5)

end if

end subroutine fbmaux_evalPhiAux
!
!****************************************************************************  
!
logical function fbmaux_IsInBoundingBox(dcorvg,kvert,iel,dPoint) 
!**@Description
  ! Check if dPoint is in the bounding box of the element iel
!**@EndDescription
implicit none
      
real*8  :: dcorvg(3,*)
integer :: kvert(8,*)
integer :: iel
real*8, dimension(3) :: dPoint
      
! local variables
logical :: found
integer :: i
real*8, dimension(3,2) :: dMinMax

dMinMax(:,1)=dcorvg(:,kvert(1,iel))
dMinMax(:,2)=dcorvg(:,kvert(1,iel))
do i=2,8

  if(dMinMax(1,1).gt.dcorvg(1,kvert(i,iel)))then
    dMinMax(1,1)=dcorvg(1,kvert(i,iel))
  end if

  if(dMinMax(2,1).gt.dcorvg(2,kvert(i,iel)))then
    dMinMax(2,1)=dcorvg(2,kvert(i,iel))
  end if

  if(dMinMax(3,1).gt.dcorvg(3,kvert(i,iel)))then
    dMinMax(3,1)=dcorvg(3,kvert(i,iel))
  end if

  if(dMinMax(1,2).lt.dcorvg(1,kvert(i,iel)))then
    dMinMax(1,2)=dcorvg(1,kvert(i,iel))
  end if

  if(dMinMax(2,2).lt.dcorvg(2,kvert(i,iel)))then
    dMinMax(2,2)=dcorvg(2,kvert(i,iel))
  end if

  if(dMinMax(3,2).lt.dcorvg(3,kvert(i,iel)))then
    dMinMax(3,2)=dcorvg(3,kvert(i,iel))
  end if

end do

if(dPoint(1).lt.dMinMax(1,1))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

if(dPoint(2).lt.dMinMax(2,1))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

if(dPoint(3).lt.dMinMax(3,1))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

if(dPoint(1).gt.dMinMax(1,2))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

if(dPoint(2).gt.dMinMax(2,2))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

if(dPoint(3).gt.dMinMax(3,2))then
  fbmaux_IsInBoundingBox = .false.
  return
end if

fbmaux_IsInBoundingBox = .true.
     
end function
!
!****************************************************************************  
!
!logical function fbmaux_PointInHex(X0,Y0,Z0,DPARX,DPARY,DPARZ,IEL,ICONV,DCORVG,KVERT)
logical function fbmaux_PointInHex(X0,Y0,Z0,x,y,z,DPARX,DPARY,DPARZ,IEL)
!**@Description
 ! The subroutine transformes a point from a deformed element back to
 ! the reference element by computing the inverse transformation and
 ! applying it to the point. The point on the reference element in returned
 ! in dparx,dpary,dparz
!**End@Description  
implicit none

integer iel
integer iconv
real*8  :: x0,y0,z0,dparx,dpary,dparz

real*8, dimension(8) :: x,y,z

integer ive,ivt,i
real*8 eps 


!real*8 dcorvg(3,*)
!integer kvert(8,*)

! a tolerance that is too high can 
! cause the algorithm to fail, especially when 
! the nodes are located on the boundary of two subdomains.
! These nodes need to have the same coordinates, so when they 
! differ only a little bit this can be problematic
eps = 0.1d-6

! assign the points of the hexahedron
!do i=1,8
!  x(i)=dcorvg(1,kvert(i,iel))
!  y(i)=dcorvg(2,kvert(i,iel))
!  z(i)=dcorvg(3,kvert(i,iel))
!end do

! transform the point on the reference element
call fbmaux_InvTrans(X,Y,Z,X0,Y0,Z0,DPARX,DPARY,DPARZ,iel,iconv)

! check if the point is in the reference element
if(.not.((dparx.ge.-1.0d0-eps).and.(dparx.le.1.0d0+eps)))then
  fbmaux_PointInHex=.false.
  return
end if

if(.not.((dpary.ge.-1.0d0-eps).and.(dpary.le.1.0d0+eps)))then
  fbmaux_PointInHex=.false.
  return
end if

if(.not.((dparz.ge.-1.0d0-eps).and.(dparz.le.1.0d0+eps)))then
  fbmaux_PointInHex=.false.
  return
end if
 
fbmaux_PointInHex = .true.

end function
!
!****************************************************************************  
!
subroutine fbmaux_InvTrans(x,y,z, x0,y0,z0,dparx,dpary,dparz,iel,iconv)
!**@Description
 ! The subroutine transformes a point from a deformed element back to
 ! the reference element by computing the inverse transformation and
 ! applying it to the point. The point on the reference element in returned
 ! in dparx,dpary,dparz
!**End@Description
use pp3d_mpi, only:myid
      
implicit none
      
! parameters
      
! coordinates of the evaluation point
real*8 :: dxreal,dyreal,dzreal

! coordinates of the element vertices
real*8, dimension(3,8) :: coord

! coordinates of (x,y,z) on the reference element
real*8 :: dparx,dpary,dparz

integer :: iconv

! coordinates on the element vertices
real*8 x(8),y(8),z(8)

! local variables
real*8 :: doldx,doldy,doldz

!      double precision q8,factor,f(3),fd(3,3)
real*8 q8,factor,f(24),jg(3,3),g(3),fd(3,3),&
                 det,detx,dety,detz,x0,y0,z0

integer          maxit,i, iel

q8 = 0.125d0

f( 1)=( x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8))*q8
f( 2)=( y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8))*q8
f( 3)=( z(1)+z(2)+z(3)+z(4)+z(5)+z(6)+z(7)+z(8))*q8

f( 4)=(-x(1)+x(2)+x(3)-x(4)-x(5)+x(6)+x(7)-x(8))*q8
f( 5)=(-y(1)+y(2)+y(3)-y(4)-y(5)+y(6)+y(7)-y(8))*q8
f( 6)=(-z(1)+z(2)+z(3)-z(4)-z(5)+z(6)+z(7)-z(8))*q8

f( 7)=(-x(1)-x(2)+x(3)+x(4)-x(5)-x(6)+x(7)+x(8))*q8
f( 8)=(-y(1)-y(2)+y(3)+y(4)-y(5)-y(6)+y(7)+y(8))*q8
f( 9)=(-z(1)-z(2)+z(3)+z(4)-z(5)-z(6)+z(7)+z(8))*q8  

f(10)=(-x(1)-x(2)-x(3)-x(4)+x(5)+x(6)+x(7)+x(8))*q8
f(11)=(-y(1)-y(2)-y(3)-y(4)+y(5)+y(6)+y(7)+y(8))*q8
f(12)=(-z(1)-z(2)-z(3)-z(4)+z(5)+z(6)+z(7)+z(8))*q8

f(13)=( x(1)-x(2)+x(3)-x(4)+x(5)-x(6)+x(7)-x(8))*q8
f(14)=( y(1)-y(2)+y(3)-y(4)+y(5)-y(6)+y(7)-y(8))*q8
f(15)=( z(1)-z(2)+z(3)-z(4)+z(5)-z(6)+z(7)-z(8))*q8

f(16)=( x(1)-x(2)-x(3)+x(4)-x(5)+x(6)+x(7)-x(8))*q8
f(17)=( y(1)-y(2)-y(3)+y(4)-y(5)+y(6)+y(7)-y(8))*q8
f(18)=( z(1)-z(2)-z(3)+z(4)-z(5)+z(6)+z(7)-z(8))*q8

f(19)=( x(1)+x(2)-x(3)-x(4)-x(5)-x(6)+x(7)+x(8))*q8
f(20)=( y(1)+y(2)-y(3)-y(4)-y(5)-y(6)+y(7)+y(8))*q8
f(21)=( z(1)+z(2)-z(3)-z(4)-z(5)-z(6)+z(7)+z(8))*q8

f(22)=(-x(1)+x(2)-x(3)+x(4)+x(5)-x(6)+x(7)-x(8))*q8
f(23)=(-y(1)+y(2)-y(3)+y(4)+y(5)-y(6)+y(7)-y(8))*q8
f(24)=(-z(1)+z(2)-z(3)+z(4)+z(5)-z(6)+z(7)-z(8))*q8

! set an epsilon for zero testing
factor= 0.1e-7

! set the start solution to 0
dparx = 0.0d0
dpary = 0.0d0
dparz = 0.0d0

maxit = 100

do i=1,maxit
  
  ! save the old solution vector
  doldx = dparx
  doldy = dpary
  doldz = dparz

  ! set up the jacobian of the mapping
  jg(1,1)=f(4)+f(13)*dpary+f(16)*dparz+f(22)*dpary*dparz

  jg(1,2)=f(7)+f(13)*dparx+f(19)*dparz+f(22)*dparx*dparz

  jg(1,3)=f(10)+f(16)*dparx+f(19)*dpary+f(22)*dparx*dpary

  jg(2,1)=f(5)+f(14)*dpary+f(17)*dparz+f(23)*dpary*dparz

  jg(2,2)=f(8)+f(14)*dparx+f(20)*dparz+f(23)*dparx*dparz

  jg(2,3)=f(11)+f(17)*dparx+f(20)*dpary+f(23)*dparx*dpary

  jg(3,1)=f(6)+f(15)*dpary+f(18)*dparz+f(24)*dpary*dparz

  jg(3,2)=f(9)+f(15)*dparx+f(21)*dparz+f(24)*dparx*dparz

  jg(3,3)=f(12)+f(18)*dparx+f(21)*dpary+f(24)*dparx*dpary


  ! compute the RHS of the linear system
  g(1)=f(1)+jg(1,1)*dparx+f(7)*dpary+f(10)*dparz+f(19)*dpary*dparz-x0

  g(2)=f(2)+f(5)*dparx+jg(2,2)*dpary+f(11)*dparz+f(17)*dparx*dparz-y0

  g(3)=f(3)+f(6)*dparx+f(9)*dpary+jg(3,3)*dparz+ f(15)*dparx*dpary-z0

  ! solve the linear system using cramer's rule
  det  =  jg(1,1)*(jg(2,2)*jg(3,3)-jg(3,2)*jg(2,3))- &
          jg(1,2)*(jg(2,1)*jg(3,3)-jg(3,1)*jg(2,3))+ &
          jg(1,3)*(jg(2,1)*jg(3,2)-jg(3,1)*jg(2,2))

  detx =  g(1)   *(jg(2,2)*jg(3,3)-jg(3,2) *jg(2,3))- &
          jg(1,2)*( g(2)  *jg(3,3)- g(3)   *jg(2,3))+ &
          jg(1,3)*( g(2)  *jg(3,2)- g(3)   *jg(2,2))

  dety =  jg(1,1)*( g(2)  *jg(3,3)- g(3)   *jg(2,3))- &
          g(1)   *(jg(2,1)*jg(3,3)-jg(3,1) *jg(2,3))+ &
          jg(1,3)*(jg(2,1)* g(3)  -jg(3,1) * g(2)  )

  detz =  jg(1,1)*(jg(2,2)* g(3)  -jg(3,2) * g(2)  )- &
          jg(1,2)*(jg(2,1)* g(3)  -jg(3,1) * g(2)  )+ &
          g(1)   *(jg(2,1)*jg(3,2)-jg(3,1) *jg(2,2))

  dparx = detx/det
  dpary = dety/det
  dparz = detz/det

  dparx = doldx-dparx
  dpary = doldy-dpary
  dparz = doldz-dparz

  ! check for divergence
  if(i.eq.10)then
  if(dparx.lt.-2.0d0.or.dparx.gt.2.0d0)then
    iconv=1
    return
  else if(dpary.lt.-2.0d0.or.dpary.gt.2.0d0)then
    iconv=1
    return
  else if(dparz.lt.-2.0d0.or.dparz.gt.2.0d0)then
    iconv=1
    return
  endif
  endif

  ! check convergence criterion
  if ((abs(dparx-doldx).le.factor).and.(abs(dpary-doldy) &
      .le.factor).and.(abs(dparz-doldz).le.factor))then
      return
  endif

end do
 
end subroutine fbmaux_InvTrans
!
!****************************************************************************  
!
logical function fbmaux_TestDomainBdry(dpoint,pid)
!**@Description
!    Tests if a point is inside the domain BdryBox
!**End@Description
use pp3d_mpi, only: myid,subnodes
use var_QuadScalar, only: mgBoundingBox


implicit none
real*8, dimension(3)  :: dpoint
integer :: pid

! local variables
integer :: i

if(dpoint(1) .lt. mgBoundingBox(pid)%vertices(1,1))then
  fbmaux_TestDomainBdry = .false.
  return
end if

if(dpoint(2) .lt. mgBoundingBox(pid)%vertices(2,1))then
  fbmaux_TestDomainBdry = .false.
  return
end if

if(dpoint(3) .lt. mgBoundingBox(pid)%vertices(3,1))then
  fbmaux_TestDomainBdry = .false.
  return
end if

if(dpoint(1) .gt. mgBoundingBox(pid)%vertices(1,2))then
  fbmaux_TestDomainBdry = .false.
  return
end if

if(dpoint(2) .gt. mgBoundingBox(pid)%vertices(2,2))then
  fbmaux_TestDomainBdry = .false.
  return
end if

if(dpoint(3) .gt. mgBoundingBox(pid)%vertices(3,2))then
  fbmaux_TestDomainBdry = .false.
  return
end if

fbmaux_TestDomainBdry = .true.

end function fbmaux_TestDomainBdry
!
!****************************************************************************  
!
subroutine fbmaux_DomainCandidates(dpoint,icandidates,inumcand)
!**@Description
!   Determines in which domains a point could be located
!   based on their bounding box
!**End@Description
use pp3d_mpi, only: myid,subnodes

implicit none
real*8, dimension(3)  :: dpoint
integer, dimension(:)  :: icandidates
integer :: inumcand

! local variables
integer :: i

inumcand = 0

do i=1,subnodes
  if(i.eq.myid)cycle

  if(fbmaux_TestDomainBdry(dpoint,i))then
    inumcand = inumcand + 1  
    icandidates(inumcand)=i
  end if

end do

end subroutine fbmaux_DomainCandidates

end module
