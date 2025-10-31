!=========================================================================
! QuadSc_boundary.f90
!
! Boundary condition functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE QuadScalar_Knpr()
  INTEGER i,ndof

  ndof  = mg_mesh%level(ilev)%nvt+&
  mg_mesh%level(ilev)%net+&
  mg_mesh%level(ilev)%nat+&
  mg_mesh%level(ilev)%nel

  QuadSc%knprU(ILEV)%x = 0
  QuadSc%knprV(ILEV)%x = 0
  QuadSc%knprW(ILEV)%x = 0

  DO i=1,ndof

  IF (myBoundary%bWall(i).OR.myBoundary%iInflow(i).NE.0) THEN
    QuadSc%knprU(ILEV)%x(i) = 1
    QuadSc%knprV(ILEV)%x(i) = 1
    QuadSc%knprW(ILEV)%x(i) = 1
  END IF

  IF (myBoundary%bSymmetry(1,i)) THEN
    QuadSc%knprU(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry u"
  END IF
  IF (myBoundary%bSymmetry(2,i)) THEN
    QuadSc%knprV(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry v"
  END IF
  IF (myBoundary%bSymmetry(3,i)) THEN
    QuadSc%knprW(ILEV)%x(i) = 1
    !  WRITE(*,*) myid,"Symmetry w"
  END IF

  END DO

  IF (myid.eq.0) THEN
    ! bNoOutFlow = .TRUE.
    ! bNoOutFlow = .FALSE.
    ! DO i=1,nvt+net+nat+nel
    !  iBndr = QuadScBoundary(i)
    !  inprU  = QuadSc%knprU(ILEV)%x(i)
    !  inprV  = QuadSc%knprV(ILEV)%x(i)
    !  inprW  = QuadSc%knprW(ILEV)%x(i)
    !  IF (iBndr.EQ.1.AND.inprU.EQ.0) bNoOutFlow = .FALSE.
    ! END DO
  END IF

END SUBROUTINE QuadScalar_Knpr
!=========================================================================
!
!=========================================================================
SUBROUTINE QuadScalar_FictKnpr(dcorvg,dcorag,kvert,kedge,karea, silent)
  use fbm, only: fbm_updateFBMGeom
#ifdef HAVE_PE
  use dem_query, only: getTotalParticles
#endif
  include 'mpif.h'

  ! Function parameters
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  logical, intent(in), optional :: silent

  ! Local variables
  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4,totalInside, reducedVal, ierr, totalP, reducedP
  REAL*8 PX,PY,PZ,DIST
  real*8 :: dofsPerParticle
  REAL tttx0,tttx1
  logical :: isSilent

  ! Constants
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  isSilent = .true.
  if(present(silent)) isSilent = silent

  totalInside = 0

  if (myid /= 0) then

#ifdef HAVE_PE
    call clear_fbm_maps()
#endif

    CALL myMPI_Barrier()
    call ztime(tttx0)

    DO i=1,nvt
    PX = dcorvg(1,I)
    PY = dcorvg(2,I)
    PZ = dcorvg(3,I)
    call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(i),FictKNPR(i),Distance(i), i, FictKNPR_uint64(i), fbm_geom_handler_ptr)
    END DO

    k=1
    DO i=1,nel
    DO j=1,12
    IF (k.eq.kedge(j,i)) THEN
      ivt1 = kvert(NeighE(1,j),i)
      ivt2 = kvert(NeighE(2,j),i)
      PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
      PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
      PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
      call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+k),FictKNPR(nvt+k),Distance(nvt+k), nvt+k, FictKNPR_uint64(nvt+k),fbm_geom_handler_ptr)
      k = k + 1
    END IF
    END DO
    END DO

    k=1
    DO i=1,nel
    DO j=1,6
    IF (k.eq.karea(j,i)) THEN
      ivt1 = kvert(NeighA(1,j),i)
      ivt2 = kvert(NeighA(2,j),i)
      ivt3 = kvert(NeighA(3,j),i)
      ivt4 = kvert(NeighA(4,j),i)
      PX = 0.25d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2)+dcorvg(1,ivt3)+dcorvg(1,ivt4))
      PY = 0.25d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2)+dcorvg(2,ivt3)+dcorvg(2,ivt4))
      PZ = 0.25d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2)+dcorvg(3,ivt3)+dcorvg(3,ivt4))
      call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+net+k),FictKNPR(nvt+net+k),Distance(nvt+net+k), nvt+net+k, FictKNPR_uint64(nvt+net+k),fbm_geom_handler_ptr)
      k = k + 1
    END IF
    END DO
    END DO

    ! DO i=1,nat
    !  PX = dcorag(1,I)
    !  PY = dcorag(2,I)
    !  PZ = dcorag(3,I)
    !  CALL GetFictKnpr(PX,PY,PZ,QuadScBoundary(nvt+net+i),FictKNPR(nvt+net+i),Distance(nvt+net+i))
    ! END DO

    DO i=1,nel
    PX = 0d0
    PY = 0d0
    PZ = 0d0
    DO j=1,8
    PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
    PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
    PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
    END DO
    call fbm_updateFBMGeom(PX,PY,PZ,QuadScBoundary(nvt+net+nat+i),FictKNPR(nvt+net+nat+i),Distance(nvt+net+nat+i), nvt+net+nat+i, FictKNPR_uint64(nvt+net+nat+i),fbm_geom_handler_ptr)
    END DO

    CALL myMPI_Barrier()
    call ztime(tttx1)
    if (myid.eq.1) WRITE(*,*) 'FBM time : ', tttx1-tttt0, ' s'

    do i=1,nvt+net+nat+nel
      myALE%Monitor(i)=distance(i)

      if(FictKnpr(i) .ne. 0)then
        totalInside = totalInside + 1
      end if

    end do

  end if ! myid /= 0

#ifdef HAVE_PE
  call MPI_Reduce(totalInside, reducedVal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (myid /= 0) then
    totalP = getTotalParticles()
  end if

  call MPI_Reduce(totalP, reducedP, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  IF (myid.eq.0) THEN
    write(*,'(A,I0)') '> Total dofs inside: ', reducedVal
    dofsPerParticle = real(reducedVal) / reducedP
    write(*,'(A,I0)') '> Dofs per Particle: ', NINT(dofsPerParticle)
  end if
#endif

END SUBROUTINE QuadScalar_FictKnpr
!=========================================================================
!
!=========================================================================
SUBROUTINE QuadScalar_FictKnpr_Wangen(dcorvg,dcorag,kvert,kedge,karea)
  use fbm, only: fbm_updateFBMGeom,myFBM
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST,P(3),dR,distX,distY
  REAL tttx0,tttx1
  INTEGER i,j,k,ivt,ndof,IP,minpr,IPP
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
  LOGICAL, allocatable :: bMark(:)
  REAL*8 :: iCount(3)

  iCount = 0d0
  ndof = nvt+nat+net+nel

  if (myid.eq.0) goto 1

  CALL myMPI_Barrier()
  call ztime(tttx0)

  ALLOCATE(bMark(ndof))
  bMark = .false.
  Distance = 1d8

  DO IPP = 1,myFBM%nParticles
  IP = IPP

  DO i=1,nel
   PX = 0d0
   PY = 0d0
   PZ = 0d0
   DO j=1,8
    PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
    PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
    PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
   END DO
   call fbm_updateFBMGeom(PX,PY,PZ,IP,minpr,distX,i, FictKNPR_uint64(i),fbm_geom_handler_ptr)

   if (distX.lt.Distance(nvt+net+nat+i)) then
    FictKNPR(nvt+net+nat+i) = minpr
    Distance(nvt+net+nat+i) = distX
    bMark(nvt+net+nat+i) = .true.
   end if
   iCount(1) = iCount(1) + 1d0

   dR = 0d0
   do j = 1,8
    ivt = kvert(j,i)
    P = dcorvg(:,ivt)
    dR = max(dR, SQRT((P(1)-PX)**2d0 + (P(2)-PY)**2d0 + (P(3)-PZ)**2d0))
   end do

   IF (1.4d0*dR.gt.distX) then
    ! nvt
    do j = 1,8
     ivt = kvert(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,i, FictKNPR_uint64(i),fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(2) = iCount(2) + 1d0
     end if
    end do

    ! net
    do j = 1,12
     ivt = nvt + kedge(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,ivt, FictKNPR_uint64(ivt),fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(2) = iCount(2) + 1d0
     end if
    end do

    ! nat
    do j = 1,6
     ivt = nvt + net + karea(j,i)
     P = dcorvg(:,ivt)
     if (.not.bMark(ivt).or.DistX.lt.Distance(ivt)) then
      call fbm_updateFBMGeom(P(1),P(2),P(3),IP,minpr,distY,ivt, FictKNPR_uint64(ivt),fbm_geom_handler_ptr)
      if (distY.lt.Distance(ivt)) then
       FictKNPR(ivt) = minpr
       Distance(ivt) = distY
       bMark(ivt) = .True.
      end if
      iCount(3) = iCount(3) + 1d0
     end if
    end do

   ELSE

    ! nvt
    do j = 1,8
     ivt = kvert(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do

    ! net
    do j = 1,12
     ivt = nvt + kedge(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do

    ! nat
    do j = 1,6
     ivt = nvt + net + karea(j,i)
     if ((.not.bMark(ivt)).and.(distX.lt.Distance(ivt))) then
      Distance(ivt) = distX
      FictKNPR(ivt) = minpr
     end if
    end do

   END IF

  END DO

  END DO ! IP

1 continue

  iCount(3) = 4*ndof
  call Comm_SummN(iCount,3)
  call ztime(tttx1)
  if (myid.eq.1) WRITE(*,'(A,ES12.4,A,I0,A,I0,A,2(ES12.4,A))') &
  'FBM time : ', tttx1-tttt0, ' s | ',int(iCount(1))+int(iCount(2)),' / ',int(iCount(3)) ,'=> ',1d2*iCount(1)/iCount(3),'%, ',1d2*iCount(2)/iCount(3),'%'
!  if (myid.eq.1) WRITE(*,'(3I10,2ES12.4)') int(iCount(1:3)),1d2*iCount(1)/iCount(3),1d2*iCount(2)/iCount(3)

  do i=1,nvt+net+nat+nel
   myALE%Monitor(i)=distance(i)
  end do
  if (myid.ne.0) then
!    do i=1,nvt+net+nat+nel
!     myALE%Monitor(i)=0d0
!     if (bMark(i))  myALE%Monitor(i)=1d0
!    end do

   DEALLOCATE(bMark)
  end if

END SUBROUTINE QuadScalar_FictKnpr_Wangen
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_QuadScalar_Def()
  INTEGER i
  REAL*8 daux

  DO i=1,QuadSc%ndof

    IF (QuadSc%knprU(ILEV)%x(i).eq.1) QuadSc%defU(i) = 0d0
    IF (QuadSc%knprV(ILEV)%x(i).eq.1) QuadSc%defV(i) = 0d0
    IF (QuadSc%knprW(ILEV)%x(i).eq.1) QuadSc%defW(i) = 0d0

    IF (FictKNPR(i).ne.0) THEN
      QuadSc%defU(i) = 0d0
      QuadSc%defV(i) = 0d0
      QuadSc%defW(i) = 0d0
    END IF
    IF (MixerKNPR(i).ne.0) THEN
      QuadSc%defU(i) = 0d0
      QuadSc%defV(i) = 0d0
      QuadSc%defW(i) = 0d0
    END IF

    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then

     DAUX = QuadSc%defU(i) * BoundaryNormal(1,i) + &
            QuadSc%defV(i) * BoundaryNormal(2,i) + &
            QuadSc%defW(i) * BoundaryNormal(3,i)

     QuadSc%defU(i) = QuadSc%defU(i) - DAUX*BoundaryNormal(1,i)
     QuadSc%defV(i) = QuadSc%defV(i) - DAUX*BoundaryNormal(2,i)
     QuadSc%defW(i) = QuadSc%defW(i) - DAUX*BoundaryNormal(3,i)

    END IF


  END DO

END SUBROUTINE Boundary_QuadScalar_Def
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_QuadScalar_Val()
  use fbm, only: fbm_velBCUpdate
  implicit none
  REAL*8 PX,PY,PZ,DAUX
  INTEGER i,inpr,finpr,minpr,inprU,inprV,inprW,ndof,iType

  ilev = NLMAX
  ndof = mg_mesh%level(ilev)%nvt + mg_mesh%level(ilev)%net +&
    mg_mesh%level(ilev)%nat + mg_mesh%level(ilev)%nel


  DO i=1,ndof
    PX = myQ2Coor(1,i);  PY = myQ2Coor(2,i);  PZ = myQ2Coor(3,i)
    inpr = 0

    IF (QuadSc%knprU(ilev)%x(i).EQ.1) QuadSc%valU(i)=0d0
    IF (QuadSc%knprV(ilev)%x(i).EQ.1) QuadSc%valV(i)=0d0
    IF (QuadSc%knprW(ilev)%x(i).EQ.1) QuadSc%valW(i)=0d0

    IF (myBoundary%iInflow(i).NE.0) THEN
      inpr = 1
      iType = myBoundary%iInflow(i)
      CALL GetVeloBCVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),iType,timens)
    END IF
    finpr = FictKNPR(i)
    minpr = MixerKNPR(i)
    IF (finpr.ne.0.and.inpr.eq.0) THEN
      CALL fbm_velBCUpdate(PX,PY,PZ,&
                           QuadSc%valU(i),QuadSc%valV(i),&
                           QuadSc%valW(i),finpr,timens,i,&
                           fbm_vel_bc_handler_ptr)
    END IF
    IF (minpr.ne.0.and.inpr.eq.0) THEN
      CALL GetVeloMixerVal(PX,PY,PZ,QuadSc%valU(i),QuadSc%valV(i),QuadSc%valW(i),minpr,timens)
    END IF
  END DO

  DO i=1,ndof
    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then
     DAUX = QuadSc%ValU(i) * BoundaryNormal(1,i) + &
            QuadSc%ValV(i) * BoundaryNormal(2,i) + &
            QuadSc%ValW(i) * BoundaryNormal(3,i)

     QuadSc%ValU(i) = QuadSc%ValU(i) - DAUX*BoundaryNormal(1,i)
     QuadSc%ValV(i) = QuadSc%ValV(i) - DAUX*BoundaryNormal(2,i)
     QuadSc%ValW(i) = QuadSc%ValW(i) - DAUX*BoundaryNormal(3,i)
    END IF
  END DO

END SUBROUTINE Boundary_QuadScalar_Val
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_LinScalar_Mat(DA,KLD,KNPRP,NEL,iS)
  REAL*8  DA(*)
  INTEGER KLD(*),KNPRP(*),ICOL,I,NEL,J,JJ,iS

  DO I=1,NEL
   IF (KNPRP(I).EQ.1) THEN
    DO J=1,4
     JJ = 4*(I-1) + J
     IF (iS.eq.1) DA(KLD(JJ)) = 1d-8
     DO ICOL=KLD(JJ)+iS,KLD(JJ+1)-1
      DA(ICOL) = 0d0
     END DO
    END DO
   END IF
  END DO

END SUBROUTINE Boundary_LinScalar_Mat
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_LinScalar_Def(DD,KNPRP,NEL)
  REAL*8  DD(*)
  INTEGER KNPRP(*),ICOL,I,NEL,J,JJ

  DO I=1,NEL
   IF (KNPRP(I).EQ.1) THEN
    DO J=1,4
     JJ = 4*(I-1) + J
     DD(JJ) = 0d0
    END DO
   END IF
  END DO

END SUBROUTINE Boundary_LinScalar_Def
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_QuadScalar_Mat(DA11,DA22,DA33,KLD,&
    KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8  DA11(*),DA22(*),DA33(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF
  REAL*8 DAUX

  DO I=1,NDOF
  IF (KNPRU(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRV(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA22(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRW(I).EQ.1) THEN
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA33(ICOL) = 0d0
    END DO
  END IF
  END DO

  DO I=1,NDOF
  IF (FictKNPR(I).NE.0.OR.MixerKNPR(I).NE.0) THEN
    !    DA(KLD(I))=1d-8
    DO ICOL=KLD(I)+1,KLD(I+1)-1
      DA11(ICOL)=0E0
      DA22(ICOL)=0E0
      DA33(ICOL)=0E0
    END DO
  END IF
  END DO

  DO I=1,NDOF
    IF (myBoundary%bSlip(i).and.(.not.(myBoundary%bWall(i).or.myBoundary%iInflow(i).gt.0))) then
    ICOL = KLD(I)
    DO ICOL=KLD(I)+1,KLD(I+1)-1
      DA11(ICOL)=0E0
      DA22(ICOL)=0E0
      DA33(ICOL)=0E0
    END DO

    END IF
  END DO



END SUBROUTINE Boundary_QuadScalar_Mat
!=========================================================================
!
!=========================================================================
SUBROUTINE Boundary_QuadScalar_Mat_9(DA11,DA22,DA33,DA12,DA13,DA23,DA21,DA31,DA32,&
    KLD,KNPRU,KNPRV,KNPRW,NDOF)
  REAL*8 DA11(*),DA22(*),DA33(*),DA12(*),DA13(*),DA23(*),DA21(*),DA31(*),DA32(*)
  INTEGER KLD(*),KNPRU(*),KNPRV(*),KNPRW(*),ICOL,I,NDOF

  DO I=1,NDOF
  IF (KNPRU(I).EQ.1) THEN
    ICOL = KLD(I)
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRV(I).EQ.1) THEN
    ICOL = KLD(I)
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA22(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    END DO
  END IF
  IF (KNPRW(I).EQ.1) THEN
    ICOL = KLD(I)
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA33(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    END DO
  END IF
  END DO

  DO I=1,NDOF
  IF (FictKNPR(I).NE.0.OR.MixerKNPR(I).NE.0) THEN
    ICOL = KLD(I)
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    DO ICOL=KLD(I)+1,KLD(I+1)-1
    DA11(ICOL) = 0d0
    DA22(ICOL) = 0d0
    DA33(ICOL) = 0d0
    DA12(ICOL) = 0d0
    DA13(ICOL) = 0d0
    DA23(ICOL) = 0d0
    DA21(ICOL) = 0d0
    DA31(ICOL) = 0d0
    DA32(ICOL) = 0d0
    END DO
  END IF
  END DO

END SUBROUTINE Boundary_QuadScalar_Mat_9
!=========================================================================
!
!=========================================================================
SUBROUTINE ProlongateSolution()

 IF (allocated(Temperature)) then
  CALL ProlongateSolutionSub(QuadSc,LinSc,Boundary_QuadScalar_Val,Temperature)
  if (myid.ne.0) Tracer%Val(NLMAX+1)%x = Temperature
 else
  CALL ProlongateSolutionSub(QuadSc,LinSc,Boundary_QuadScalar_Val)
 end if
 CALL QuadScP1toQ2(LinSc,QuadSc)

END SUBROUTINE ProlongateSolution
!=========================================================================
!
!=========================================================================
SUBROUTINE InitBoundaryList(knpr,kvert,kedge,karea)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*),knpr(*)
  INTEGER i,j,k,ivt1,ivt2
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  DO i=1,nvt
  QuadScBoundary(i) = mg_mesh%level(ilev)%knpr(i)
  !  IF (QuadScBoundary(i).eq.1) write(*,*) "type 1"
  END DO

  k=1
  DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    !IF (knpr(ivt1).EQ.1.AND.knpr(ivt2).EQ.1) THEN
    QuadScBoundary(nvt+k) = 0
    !END IF
    !    IF (QuadScBoundary(nvt+k).eq.1) write(*,*) "type 2"
    k = k + 1
  END IF
  END DO
  END DO

  DO i=1,nat
  !QuadScBoundary(nvt+net+i) = knpr(nvt+i)
  QuadScBoundary(nvt+net+i) = 0! knpr(nvt+i)
  !  IF (QuadScBoundary(nvt+net+i).eq.1) write(*,*) "type 3"
  END DO

  DO i=1,nel
  QuadScBoundary(nvt+net+nat+i) = 0
  END DO

END SUBROUTINE InitBoundaryList
!=========================================================================
!
!=========================================================================
SUBROUTINE GetMG_KNPRP(mgMesh)
type(tMultiMesh), intent(inout) :: mgMesh
integer iel,jel(8),i,jjj

jjj = 0
do iel = 1,nel/8
 if (LinSc%knprP(ilev-1)%x(iel).eq.1) then
  JEL(1)  = iel
  JEL(2)  = mgMesh%level(ilev)%kadj(3,JEL(1))
  JEL(3)  = mgMesh%level(ilev)%kadj(3,JEL(2))
  JEL(4)  = mgMesh%level(ilev)%kadj(3,JEL(3))
  JEL(5)  = mgMesh%level(ilev)%kadj(6,JEL(1))
  JEL(6)  = mgMesh%level(ilev)%kadj(6,JEL(2))
  JEL(7)  = mgMesh%level(ilev)%kadj(6,JEL(3))
  JEL(8)  = mgMesh%level(ilev)%kadj(6,JEL(4))
  do i=1,8
    LinSc%knprP(ilev)%x(jel(i)) = 1
    jjj = jjj + 1
  end do
 end if
end do

if (myid.eq.1) write(*,'(A,I0,A,I0)') 'KNPRP nodes on level ',ilev, ' are :', jjj

END SUBROUTINE GetMG_KNPRP
!=========================================================================
!
!=========================================================================
SUBROUTINE IncludeFBM_BCs(mgMesh)

type(tMultiMesh), intent(inout) :: mgMesh
integer iel,i

do iel = 1,nel
 if (LinSc%knprP(ilev)%x(iel).eq.1) then
  do i=1,8
   myBoundary%bWall(mgMesh%level(ilev)%kvert(i,iel)) = .true.
  end do

  do i=1,12
   myBoundary%bWall(nvt + mgMesh%level(ilev)%kedge(i,iel)) = .true.
  end do

  do i=1,6
   myBoundary%bWall(nvt + net + mgMesh%level(ilev)%karea(i,iel)) = .true.
  end do

  myBoundary%bWall(nvt + net + nat + iel) = .true.
 end if
end do

QuadSc%auxU = 0d0
DO i=1,QuadSc%ndof
 if (myBoundary%bWall(i)) QuadSc%auxU(i) = 1d0
END DO

CALL E013Sum(QuadSc%auxU)

DO i=1,QuadSc%ndof
 if (QuadSc%auxU(i).gt.0d0) myBoundary%bWall(i) = .true.
END DO

!!!! COMMUNICATION needed !!!!!


END SUBROUTINE IncludeFBM_BCs
!=========================================================================
!
!=========================================================================
SUBROUTINE  Setup_PeriodicVelocityRHS()
INTEGER i

IF (.NOT.(ALLOCATED(dPeriodicVector))) ALLOCATE(dPeriodicVector(QuadSc%ndof))

IF (ieee_is_finite(myProcess%dPress)) THEN
!IF (.NOT.ISNAN(myProcess%dPress)) THEN

DO i=1,SIZE(LinSc%AuxP(NLMAX)%x)
 IF (MOD(i,4).EQ.1) then
  LinSc%AuxP(NLMAX)%x(i) = myProcess%dPress
 ELSE
  LinSc%AuxP(NLMAX)%x(i) = 0d0
 END IF
END DO

CALL B_Mul_U(qlMat%ColA,qlMAt%LdA,BXMat,BYMat,BZMat,LinSc%auxP(NLMAX)%x,&
     QuadSc%AuxU,QuadSc%AuxV,QuadSc%AuxW,QuadSc%ndof,+1d0,0d0)

DO I=1,QuadSc%ndof
 IF ((myQ2Coor(3,i)     .LT.+1e-3)) THEN
  dPeriodicVector(i) = QuadSc%auxW(i)
 ELSE
  dPeriodicVector(i) = 0d0
 END IF
END DO
ELSE
  dPeriodicVector = 0d0
END IF

END SUBROUTINE  Setup_PeriodicVelocityRHS
!=========================================================================
!
!=========================================================================
SUBROUTINE RestrictWallBC()
INTEGER ndof
INTEGER i,j,k
integer iat,ivt1,ivt2,ivt3,ivt4
INTEGER NeighA(4,6),NeighU(4,6)
DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/
DATA NeighU/1,2,3,4, 1,6,9,5, 2,7,10,6, 3,8,11,7, 4,5,12,8, 9,10,11,12/

 ndof  = mg_mesh%level(ilev)%nvt+&
 mg_mesh%level(ilev)%net+&
 mg_mesh%level(ilev)%nat+&
 mg_mesh%level(ilev)%nel

 nvt  = mg_mesh%level(ilev)%nvt
 net  = mg_mesh%level(ilev)%net
 nat  = mg_mesh%level(ilev)%nat
 nel  = mg_mesh%level(ilev)%nel

 ! overlap outper points with wall ==> inner 'wall' markers will disappear
 DO i=1,ndof
  IF (myBoundary%bWall(i)) THEN
   IF (.not.(mg_mesh%BndryNodes(i)%bOuterPoint)) myBoundary%bWall(i) = .FALSE.
  END IF
 END DO

!  remove all outflow / inflow dofs
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%iInflow(nvt+net+k).ne.0.or.myBoundary%bOutflow(nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .false.
      myBoundary%bWall(nvt+net+k) = .false.
     end if
     k = k + 1
    END IF
   END DO
  END DO

!  remove all symmetry dofs
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%bSymmetry(1,nvt+net+k).or.myBoundary%bSymmetry(2,nvt+net+k).or.myBoundary%bSymmetry(3,nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .false.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .false.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .false.
      myBoundary%bWall(nvt+net+k) = .false.
     end if
     k = k + 1
    END IF
   END DO
  END DO

  ! refresh Wall
  k=1
  DO i=1,nel
   DO j=1,6
    IF (k.eq.mg_mesh%level(ilev)%karea(j,i)) THEN
     if (myBoundary%bWall(nvt+net+k)) then
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(1,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(2,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(3,j),i)) = .true.
      myBoundary%bWall(mg_mesh%level(ilev)%kvert(NeighA(4,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(1,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(2,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(3,j),i)) = .true.
      myBoundary%bWall(nvt+mg_mesh%level(ilev)%kedge(NeighU(4,j),i)) = .true.
     end if
     k = k + 1
    END IF
   END DO
  END DO

END SUBROUTINE RestrictWallBC
