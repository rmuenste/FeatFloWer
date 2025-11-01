MODULE fbm_particle_reynolds

USE PP3D_MPI,     ONLY : myid, COMM_SUMMN
USE var_QuadScalar, ONLY : myFBM, mg_mesh

IMPLICIT NONE

CONTAINS

SUBROUTINE fbm_compute_particle_reynolds(U1,U2,U3,DVISC,RHOFLUID,mfile)
  USE fbmaux, ONLY : fbmaux_PointInHex
  INTEGER, PARAMETER :: NNBAS=27, NNDER=10, NNVE=8, NNDIM=3

  REAL*8, DIMENSION(:), INTENT(IN) :: U1,U2,U3
  REAL*8, DIMENSION(:), INTENT(IN) :: DVISC
  REAL*8,               INTENT(IN) :: RHOFLUID
  INTEGER,              INTENT(IN) :: mfile

  INTEGER :: ip,iel,ive,ig,il,ieltTyp,ilev,nel_local
  REAL*8 :: xi1,xi2,xi3
  REAL*8 :: vel_sample(3), slip(3), speed, diameter, mu_loc
  REAL*8 :: xverts(8), yverts(8), zverts(8)
  REAL*8 :: xmin,xmax,ymin,ymax,zmin,zmax,eps_box
  LOGICAL :: found

  REAL*8, ALLOCATABLE :: re_local(:),re_weight(:)
!
!  INTEGER KDFG(NNBAS),KDFL(NNBAS)
!  REAL*8  DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ
!  REAL*8  DBAS(NNDIM,NNBAS,NNDER)
!  LOGICAL BDER(NNDER)
!  INTEGER KVE(NNVE),NDIM,IEL
!
!  COMMON /ELEM/   DX,DY,DZ,DJAC,DETJ,DBAS,BDER,KVE,IEL,NDIM
!  COMMON /COAUX1/ KDFG,KDFL,IDFL
!  INTEGER :: IDFL
!
!  EXTERNAL :: ELE, NDFGL, SETLEV
!  INTEGER :: NDFL
!
!  IF (.NOT.ALLOCATED(myFBM%ParticleRe)) RETURN
!  IF (myFBM%nParticles.EQ.0) THEN
!    myFBM%ParticleRe = 0D0
!    RETURN
!  END IF
!
!  ilev = mg_mesh%nlmax
!  CALL SETLEV(2)
!  nel_local = mg_mesh%level(ilev)%nel
!  IF (SIZE(DVISC).LT.nel_local) THEN
!    myFBM%ParticleRe = 0D0
!    RETURN
!  END IF
!
!  ALLOCATE(re_local(myFBM%nParticles))
!  ALLOCATE(re_weight(myFBM%nParticles))
!  re_local  = 0D0
!  re_weight = 0D0
!
!  DO ive = 1,NNDER
!    BDER(ive) = .FALSE.
!  END DO
!  BDER(1) = .TRUE.
!
!  IELTYP = -1
!  CALL ELE(0D0,0D0,0D0,IELTYP)
!  IDFL = NDFL(IELTYP)
!
!  eps_box = 1D-10
!
!  IF (myid.NE.0) THEN
!    DO ip = 1,myFBM%nParticles
!      vel_sample = 0D0
!      found = .FALSE.
!      xi1 = 0D0
!      xi2 = 0D0
!      xi3 = 0D0
!
!      DO iel = 1,nel_local
!        DO ive = 1,NNVE
!          ig = mg_mesh%level(ilev)%kvert(ive,iel)
!          xverts(ive) = mg_mesh%level(ilev)%dcorvg(1,ig)
!          yverts(ive) = mg_mesh%level(ilev)%dcorvg(2,ig)
!          zverts(ive) = mg_mesh%level(ilev)%dcorvg(3,ig)
!        END DO
!
!        xmin = MINVAL(xverts)
!        xmax = MAXVAL(xverts)
!        ymin = MINVAL(yverts)
!        ymax = MAXVAL(yverts)
!        zmin = MINVAL(zverts)
!        zmax = MAXVAL(zverts)
!
!        IF (myFBM%ParticleNew(ip)%Position(1).LT.xmin-eps_box) CYCLE
!        IF (myFBM%ParticleNew(ip)%Position(1).GT.xmax+eps_box) CYCLE
!        IF (myFBM%ParticleNew(ip)%Position(2).LT.ymin-eps_box) CYCLE
!        IF (myFBM%ParticleNew(ip)%Position(2).GT.ymax+eps_box) CYCLE
!        IF (myFBM%ParticleNew(ip)%Position(3).LT.zmin-eps_box) CYCLE
!        IF (myFBM%ParticleNew(ip)%Position(3).GT.zmax+eps_box) CYCLE
!
!        xi1 = 0D0
!        xi2 = 0D0
!        xi3 = 0D0
!        found = fbmaux_PointInHex(myFBM%ParticleNew(ip)%Position(1), &
!                                  myFBM%ParticleNew(ip)%Position(2), &
!                                  myFBM%ParticleNew(ip)%Position(3), &
!                                  xverts,yverts,zverts,xi1,xi2,xi3,iel)
!        IF (.NOT.found) CYCLE
!
!        DO ive = 1,NNDER
!          BDER(ive) = .FALSE.
!        END DO
!        BDER(1) = .TRUE.
!
!        CALL NDFGL(iel,1,IELTYP,mg_mesh%level(ilev)%kvert, &
!                   mg_mesh%level(ilev)%kedge,mg_mesh%level(ilev)%karea, &
!                   KDFG,KDFL)
!
!        DO ive = 1,NNVE
!          ig = mg_mesh%level(ilev)%kvert(ive,iel)
!          KVE(ive) = ig
!          DX(ive) = xverts(ive)
!          DY(ive) = yverts(ive)
!          DZ(ive) = zverts(ive)
!        END DO
!
!        CALL ELE(xi1,xi2,xi3,0)
!
!        vel_sample = 0D0
!        DO il = 1,IDFL
!          ig = KDFG(il)
!          ive = KDFL(il)
!          vel_sample(1) = vel_sample(1) + U1(ig)*DBAS(1,ive,1)
!          vel_sample(2) = vel_sample(2) + U2(ig)*DBAS(1,ive,1)
!          vel_sample(3) = vel_sample(3) + U3(ig)*DBAS(1,ive,1)
!        END DO
!
!        mu_loc = DVISC(iel)
!        slip = vel_sample - myFBM%ParticleNew(ip)%Velocity
!        speed = SQRT(slip(1)**2 + slip(2)**2 + slip(3)**2)
!        diameter = 2D0*myFBM%ParticleNew(ip)%sizes(1)
!
!        IF (mu_loc.GT.0D0.AND.diameter.GT.0D0) THEN
!          re_local(ip) = RHOFLUID*speed*diameter/mu_loc
!        ELSE
!          re_local(ip) = 0D0
!        END IF
!        re_weight(ip) = 1D0
!        IF (ALLOCATED(myFBM%iel_ug)) THEN
!          IF (ip.LE.SIZE(myFBM%iel_ug)) myFBM%iel_ug(ip) = iel
!        END IF
!        EXIT
!      END DO
!    END DO
!  END IF
!
!  CALL COMM_SUMMN(re_local,myFBM%nParticles)
!  CALL COMM_SUMMN(re_weight,myFBM%nParticles)
!
!  DO ip = 1,myFBM%nParticles
!    IF (re_weight(ip).GT.0D0) THEN
!      myFBM%ParticleRe(ip) = re_local(ip)/re_weight(ip)
!    ELSE
!      myFBM%ParticleRe(ip) = 0D0
!    END IF
!  END DO
!
!  IF (myid.EQ.1) THEN
!    IF (myFBM%nParticles.GT.0) THEN
!      WRITE(mfile,'(A,2E16.8)') 'Particle Reynolds number range: ', &
!           MINVAL(myFBM%ParticleRe), MAXVAL(myFBM%ParticleRe)
!    END IF
!  END IF
!
!  DEALLOCATE(re_local)
!  DEALLOCATE(re_weight)

END SUBROUTINE fbm_compute_particle_reynolds

END MODULE fbm_particle_reynolds
