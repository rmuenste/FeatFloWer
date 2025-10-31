!=========================================================================
! QuadSc_mesh_operations.f90
!
! Mesh deformation and ALE (Arbitrary Lagrangian-Eulerian) operations
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE  STORE_OLD_MESH(dcoor)
  REAL*8 dcoor(3,*)
  INTEGER i

  !write(*,*)'ndof:',QuadSc%ndof,mg_mesh%level(NLMAX+1)%nvt
  DO i=1,QuadSc%ndof
  myALE%OldCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  STORE_OLD_MESH
!=========================================================================
!
!=========================================================================
SUBROUTINE  STORE_NEW_MESH(dcoor)
  REAL*8 dcoor(3,*)
  INTEGER i

  DO i=1,QuadSc%ndof
  myALE%NewCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  STORE_NEW_MESH
!=========================================================================
!
!=========================================================================
SUBROUTINE  GET_MESH_VELO()
  INTEGER i
  REAL*8 dmax,daux

  dmax = 0d0
  DO i=1,QuadSc%ndof
  myALE%MeshVelo(:,i) = (myALE%NewCoor(:,i)-myALE%OldCoor(:,i))/tstep
  ! myALE%MeshVelo(:,i) = (myALE%NewCoor(:,i)-myALE%OldCoor(:,i))/tstep
  daux = SQRT(myALE%MeshVelo(1,i)**2d0 + myALE%MeshVelo(2,i)**2d0 +&
    myALE%MeshVelo(3,i)**2d0)
  IF (daux.gt.dmax) dmax = daux
  END DO

  CALL COMM_Maximum(dmax)
  IF (myid.eq.showid) WRITE(*,*) 'max mesh velocity: ', dmax

END SUBROUTINE  GET_MESH_VELO
!=========================================================================
!
!=========================================================================
SUBROUTINE  RotateMyMesh(dcoor)
  REAL*8 dcoor(3,*)
  REAL*8 X,Y,Z,R,A,dAlpha
  REAL*8 :: PI = 3.141592654d0
  INTEGER i,nnn

  dAlpha = 1d0*(250d0/60d0)*2d0*timens*PI
  ! dAlpha = 2d0*(250d0/60d0)*2d0*tstep*PI

  IF (myid.eq.0) THEN
    nnn = KNVT(NLMAX)
  ELSE
    nnn = KNVT(NLMAX+1)
  END IF

  DO i=1,nnn
  !  X = dcoor(1,i)
  !  Y = dcoor(2,i)
  X = myALE%OrigCoor(1,i)
  Y = myALE%OrigCoor(2,i)
  R = DSQRT(X*X + Y*Y)
  A = DATAN(Y/X)
  IF (X.LT.0d0) A = A + PI
  A = A + dAlpha
  dcoor(1,i) = R*DCOS(A)
  dcoor(2,i) = R*DSIN(A)
  !IF (myid.eq.1) WRITE(*,'(2(2E12.4,A))') X,Y,' : ', dcoor(1:2,i)
  END DO

  ILEV = NLMIN
  CALL SETLEV(2)

  CALL ExchangeNodeValuesOnCoarseLevel(DWORK(L(LCORVG)),KWORK(L(LVERT)),NVT,NEL)

  ILEV = NLMAX
  CALL SETLEV(2)

END SUBROUTINE  RotateMyMesh
!=========================================================================
!
!=========================================================================
SUBROUTINE  StoreOrigCoor(dcoor)
  INTEGER i,nnn
  REAL*8 dcoor(3,*)

  IF (myid.eq.0) THEN
    nnn = KNVT(NLMAX)
  ELSE
    nnn = KNVT(NLMAX+1)
  END IF

  DO i=1,nnn
  myALE%OrigCoor(:,i) = dcoor(:,i)
  END DO

END SUBROUTINE  StoreOrigCoor
!=========================================================================
!
!=========================================================================
SUBROUTINE GetMeshVelocity2(mfile)
  integer mfile
  REAL*8 dMaxVelo,daux
  INTEGER i

  dMaxVelo = 0d0
  IF (myid.ne.0) then
    DO i=1,QuadSc%ndof
    myALE%MeshVelo(:,i) = (myQ2Coor(:,i) -  myALE%Q2Coor_old(:,i))/tstep
    daux = myALE%MeshVelo(1,i)**2d0+myALE%MeshVelo(2,i)**2d0+myALE%MeshVelo(3,i)**2d0
    IF (dMaxVelo.lt.daux) dMaxVelo = daux
    END DO
  END IF

  CALL COMM_Maximum(dMaxVelo)

  IF (myid.eq.1) THEN
    WRITE(mfile,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
    WRITE(mterm,*)  "Maximum Mesh Velocity: ", SQRT(dMaxVelo)
  END IF

END SUBROUTINE GetMeshVelocity2
!=========================================================================
!
!=========================================================================
SUBROUTINE StaticMeshAdaptation()
  INTEGER iAdaptMeshLevel,idL

  IF (.NOT.bMeshAdaptation) RETURN

  OPEN(474,FILE=ADJUSTL(TRIM(cAdaptedMeshFile))//"/level.prf")
  READ(474,*) iAdaptMeshLevel
  CLOSE(474)

  idL = iAdaptMeshLevel - NLMAX

  CALL CreateDumpStructures(idL)
  CALL LoadSmartAdaptedMeshFile(DWORK(L(KLCVG(1))),cAdaptedMeshFile,idL)

  ! ---------------------------- -  -    - -- -- - - - - - - -   ----------------

  DO ILEV = iAdaptMeshLevel,NLMAX
  CALL SETLEV(2)
  WRITE(*,*) 'mesh levels', ilev,iAdaptMeshLevel,NLMAX
  CALL RefreshCoordinates(DWORK(L(KLCVG(ILEV+1))),DWORK(L(KLCAG(ILEV))),&
    KWORK(L(KLVERT(ILEV))),KWORK(L(KLEDGE(ILEV))),KWORK(L(KLAREA(ILEV))))
  END DO

  ! ---------------------------- -  -    - -- -- - - - - - - -   ----------------

END SUBROUTINE StaticMeshAdaptation
!=========================================================================
!
!=========================================================================
SUBROUTINE CoorWriter(dcorvg,nvt,cF)
  CHARACTER*(*) cF
  REAL*8 dcorvg(3,*)
  INTEGER i,nvt

  OPEN(547,FILE=TRIM(ADJUSTL(cF)))
  DO i=1,nvt
  WRITE(547,*) dcorvg(:,i)
  END DO
  CLOSE(547)

END SUBROUTINE CoorWriter
!=========================================================================
!
!=========================================================================
SUBROUTINE RefreshCoordinates(dcorvg,dcorag,kvert,kedge,karea)
  REAL*8  dcorvg(3,*),dcorag(3,*)
  INTEGER kvert(8,*),kedge(12,*),karea(6,*)
  REAL*8 PX,PY,PZ,DIST
  INTEGER i,j,k,ivt1,ivt2,ivt3,ivt4
  INTEGER NeighE(2,12),NeighA(4,6)
  DATA NeighE/1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5/
  DATA NeighA/1,2,3,4,1,2,6,5,2,3,7,6,3,4,8,7,4,1,5,8,5,6,7,8/

  k=1
  DO i=1,nel
  DO j=1,12
  IF (k.eq.kedge(j,i)) THEN
    ivt1 = kvert(NeighE(1,j),i)
    ivt2 = kvert(NeighE(2,j),i)
    PX = 0.5d0*(dcorvg(1,ivt1)+dcorvg(1,ivt2))
    PY = 0.5d0*(dcorvg(2,ivt1)+dcorvg(2,ivt2))
    PZ = 0.5d0*(dcorvg(3,ivt1)+dcorvg(3,ivt2))
    dcorvg(1,nvt+k) = PX
    dcorvg(2,nvt+k) = PY
    dcorvg(3,nvt+k) = PZ
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
    dcorag(1,k) = PX
    dcorag(2,k) = PY
    dcorag(3,k) = PZ
    dcorvg(1,nvt+net+k) = PX
    dcorvg(2,nvt+net+k) = PY
    dcorvg(3,nvt+net+k) = PZ
    k = k + 1
  END IF
  END DO
  END DO

  DO i=1,nel
  PX = 0d0
  PY = 0d0
  PZ = 0d0
  DO j=1,8
  PX = PX + 0.125d0*(dcorvg(1,kvert(j,i)))
  PY = PY + 0.125d0*(dcorvg(2,kvert(j,i)))
  PZ = PZ + 0.125d0*(dcorvg(3,kvert(j,i)))
  END DO
  dcorvg(1,nvt+net+nat+i) = PX
  dcorvg(2,nvt+net+nat+i) = PY
  dcorvg(3,nvt+net+nat+i) = PZ
  END DO

END SUBROUTINE RefreshCoordinates
