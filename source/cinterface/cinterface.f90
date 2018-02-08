module cinterface

  implicit none

  contains
  !
  !----------------------------------------------
  !
  logical function calculateDynamics()
    use var_QuadScalar, only:myFBM

    implicit none 

    integer :: dodyn

      call get_do_dynamics(dodyn)

      if ((myFBM%nParticles.gt.0).and.(dodyn.gt.0)) then
       calculateDynamics = .true.
      else
       calculateDynamics = .false.
      end if

  end function calculateDynamics
  !
  !----------------------------------------------
  !
  logical function calculateFBM()
    use var_QuadScalar, only:myFBM
    implicit none 

    integer :: dofbm

      call get_do_fbm(dofbm)

      if((myFBM%nParticles.gt.0).and.(dofbm.gt.0)) then
       calculateFBM  = .true.
      else
       calculateFBM  = .false.
      end if

  end function calculateFBM
  !
  !----------------------------------------------
  !
  SUBROUTINE FBM_GetParticleStateUpdate()
    USE PP3D_MPI, ONLY:myid,showid
    USE var_QuadScalar,ONLY:myFBM
    INTEGER iP,iParticles,ipc
    REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
    REAL*8 posx,posy,posz

    DO iP = 1,myFBM%nParticles
    ipc=iP-1

    call update_particle_state(posx,posy,posz,&
      velx,vely,velz,&
      angx,angy,angz,&
      angvelx,angvely,angvelz,&
      ipc)

    myFBM%ParticleNew(iP)%Position(1)=posx
    myFBM%ParticleNew(iP)%Position(2)=posy
    myFBM%ParticleNew(iP)%Position(3)=posz

    myFBM%ParticleNew(iP)%Velocity(1)=velx
    myFBM%ParticleNew(iP)%Velocity(2)=vely
    myFBM%ParticleNew(iP)%Velocity(3)=velz 

    myFBM%ParticleNew(iP)%Angle(1)=angx
    myFBM%ParticleNew(iP)%Angle(2)=angy
    myFBM%ParticleNew(iP)%Angle(3)=angz

    myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
    myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
    myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

    END DO

  END SUBROUTINE FBM_GetParticleStateUpdate
  !
  !----------------------------------------------
  !
  SUBROUTINE FBM_SetParticles()
    USE PP3D_MPI, ONLY:myid,showid
    USE var_QuadScalar,ONLY:myFBM
    INTEGER iP,iParticles,ipc
    REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
    REAL*8 posx,posy,posz

    DO iP = 1,myFBM%nParticles
    ipc=iP-1

    call getpos(posx,posy,posz,ipc)
    call getvel(velx,vely,velz,ipc)
    call getangle(angx,angy,angz,ipc)
    call getangvel(angvelx,angvely,angvelz,ipc)

    myFBM%ParticleNew(iP)%Position(1)=posx
    myFBM%ParticleNew(iP)%Position(2)=posy
    myFBM%ParticleNew(iP)%Position(3)=posz

    myFBM%ParticleNew(iP)%Velocity(1)=velx
    myFBM%ParticleNew(iP)%Velocity(2)=vely
    myFBM%ParticleNew(iP)%Velocity(3)=velz 

    myFBM%ParticleNew(iP)%Angle(1)=angx
    myFBM%ParticleNew(iP)%Angle(2)=angy
    myFBM%ParticleNew(iP)%Angle(3)=angz

    myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
    myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
    myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

    END DO

  END SUBROUTINE FBM_SetParticles
  !
  !----------------------------------------------
  !
  SUBROUTINE FBM_GetSoftParticles
    USE PP3D_MPI, ONLY:myid,showid
    USE var_QuadScalar,ONLY:myFBM
    INTEGER iP,iParticles,ipc
    REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
    REAL*8 posx,posy,posz

    DO iP = 1,myFBM%nParticles
    ipc=iP-1

    call getpos(posx,posy,posz,ipc)
    call getvel(velx,vely,velz,ipc)
    call getangle(angx,angy,angz,ipc)
    call getangvel(angvelx,angvely,angvelz,ipc)

    myFBM%ParticleNew(iP)%Position(1)=posx
    myFBM%ParticleNew(iP)%Position(2)=posy
    myFBM%ParticleNew(iP)%Position(3)=posz

    myFBM%ParticleNew(iP)%Velocity(1)=velx
    myFBM%ParticleNew(iP)%Velocity(2)=vely
    myFBM%ParticleNew(iP)%Velocity(3)=velz 

    myFBM%ParticleNew(iP)%Angle(1)=angx
    myFBM%ParticleNew(iP)%Angle(2)=angy
    myFBM%ParticleNew(iP)%Angle(3)=angz

    myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
    myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
    myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

    END DO

  END SUBROUTINE FBM_GetSoftParticles
  !
  !----------------------------------------------
  !
  SUBROUTINE FBM_GetParticles(icorr)
    USE PP3D_MPI, ONLY:myid,showid
    USE var_QuadScalar,ONLY:myFBM
    integer, optional :: icorr
    INTEGER iP,iParticles,ipc,index_correction
    REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
    REAL*8 posx,posy,posz,drad,density
    REAL*8 forcex,forcey,forcez
    REAL*8 torquex,torquey,torquez

    if(present(icorr))then
      index_correction = icorr
    else
      index_correction = 1
    end if

    call getnumparticles(iParticles)

    myFBM%nParticles=iParticles-index_correction

    ALLOCATE(myFBM%ParticleNew(myFBM%nParticles),myFBM%ParticleOld(myFBM%nParticles))
    ALLOCATE(myFBM%iel_ug(myFBM%nParticles))

    myFBM%iel_ug=0

    ALLOCATE(myFBM%Force(6*myFBM%nParticles))

    DO iP = 1,myFBM%nParticles
    ipc=iP-1
    myFBM%ParticleOld(iP)%cType = 'Sphere'

    ! get the density
    call getdensity(density,ipc)
    myFBM%ParticleOld(iP)%density = density
    myFBM%ParticleNew(iP)%density = density

    ! get the radius
    call getradius(drad,ipc)
    myFBM%ParticleOld(iP)%sizes(1)=drad
    myFBM%ParticleNew(iP)%sizes(1)=drad

    call getpos(posx,posy,posz,ipc)
    call getvel(velx,vely,velz,ipc)
    call getangle(angx,angy,angz,ipc)
    call getangvel(angvelx,angvely,angvelz,ipc)

    call getforce(forcex,forcey,forcez,ipc)
    call gettorque(torquex,torquey,torquez,ipc)

    myFBM%ParticleOld(iP)%Position(1)=posx
    myFBM%ParticleOld(iP)%Position(2)=posy
    myFBM%ParticleOld(iP)%Position(3)=posz

    myFBM%ParticleNew(iP)%Position(1)=posx
    myFBM%ParticleNew(iP)%Position(2)=posy
    myFBM%ParticleNew(iP)%Position(3)=posz

    myFBM%ParticleOld(iP)%Velocity(1)=velx
    myFBM%ParticleOld(iP)%Velocity(2)=vely
    myFBM%ParticleOld(iP)%Velocity(3)=velz 

    myFBM%ParticleNew(iP)%Velocity(1)=velx
    myFBM%ParticleNew(iP)%Velocity(2)=vely
    myFBM%ParticleNew(iP)%Velocity(3)=velz 

    myFBM%ParticleOld(iP)%Angle(1)=angx
    myFBM%ParticleOld(iP)%Angle(2)=angy
    myFBM%ParticleOld(iP)%Angle(3)=angz

    myFBM%ParticleNew(iP)%Angle(1)=angx
    myFBM%ParticleNew(iP)%Angle(2)=angy
    myFBM%ParticleNew(iP)%Angle(3)=angz

    myFBM%ParticleOld(iP)%AngularVelocity(1)=angvelx
    myFBM%ParticleOld(iP)%AngularVelocity(2)=angvely
    myFBM%ParticleOld(iP)%AngularVelocity(3)=angvelz

    myFBM%ParticleNew(iP)%AngularVelocity(1)=angvelx
    myFBM%ParticleNew(iP)%AngularVelocity(2)=angvely
    myFBM%ParticleNew(iP)%AngularVelocity(3)=angvelz

    myFBM%ParticleOld(iP)%Acceleration(1)=0
    myFBM%ParticleOld(iP)%Acceleration(2)=0
    myFBM%ParticleOld(iP)%Acceleration(3)=0

    myFBM%ParticleNew(iP)%Acceleration(1)=0
    myFBM%ParticleNew(iP)%Acceleration(2)=0
    myFBM%ParticleNew(iP)%Acceleration(3)=0

    myFBM%ParticleOld(iP)%FrameVelocity(1)=0
    myFBM%ParticleOld(iP)%FrameVelocity(2)=0
    myFBM%ParticleOld(iP)%FrameVelocity(3)=0

    myFBM%ParticleOld(iP)%ResistanceForce(1)=forcex
    myFBM%ParticleOld(iP)%ResistanceForce(2)=forcey
    myFBM%ParticleOld(iP)%ResistanceForce(3)=forcez

    myFBM%ParticleNew(iP)%ResistanceForce(1)=forcex
    myFBM%ParticleNew(iP)%ResistanceForce(2)=forcey
    myFBM%ParticleNew(iP)%ResistanceForce(3)=forcez

    myFBM%ParticleOld(iP)%TorqueForce(1)=torquex
    myFBM%ParticleOld(iP)%TorqueForce(2)=torquey
    myFBM%ParticleOld(iP)%TorqueForce(3)=torquez

    myFBM%ParticleNew(iP)%TorqueForce(1)=torquex
    myFBM%ParticleNew(iP)%TorqueForce(2)=torquey
    myFBM%ParticleNew(iP)%TorqueForce(3)=torquez

    END DO

  END SUBROUTINE FBM_GetParticles
  !
  !----------------------------------------------
  !
  SUBROUTINE FBM_ScatterParticles
    USE PP3D_MPI, ONLY:myid,showid,SENDD_myMPI,RECVD_myMPI,subnodes
    USE var_QuadScalar,ONLY:myFBM
    IMPLICIT NONE
    INTEGER iP,iParticles,ipc,pID
    REAL*8 velx,vely,velz,angx,angy,angz,angvelx,angvely,angvelz
    REAL*8 posx,posy,posz
    real*8, allocatable :: particleArray(:)

    if(myid.eq.0)then

      ! master copies particle data from gpu
      DO iP = 1,myFBM%nParticles
      ipc=iP-1
      call getpos(posx,posy,posz,ipc)
      call getvel(velx,vely,velz,ipc)
      call getangle(angx,angy,angz,ipc)
      call getangvel(angvelx,angvely,angvelz,ipc)

      myFBM%ParticleNew(iP)%Position(1)=posx
      myFBM%ParticleNew(iP)%Position(2)=posy
      myFBM%ParticleNew(iP)%Position(3)=posz
      myFBM%ParticleNew(iP)%Velocity(1)=velx
      myFBM%ParticleNew(iP)%Velocity(2)=vely
      myFBM%ParticleNew(iP)%Velocity(3)=velz 
      END DO

    end if

    allocate(particleArray(6*myFBM%nParticles))
    particleArray=0d0
    ! master scatters data to non-master processes
    IF(myid.eq.0)THEN

      ! copy particle date into an array for communication
      DO iP = 0,myFBM%nParticles-1
      particleArray(iP*6+1)=myFBM%ParticleNew(iP+1)%Position(1)
      particleArray(iP*6+2)=myFBM%ParticleNew(iP+1)%Position(2)
      particleArray(iP*6+3)=myFBM%ParticleNew(iP+1)%Position(3)
      particleArray(iP*6+4)=myFBM%ParticleNew(iP+1)%Velocity(1)
      particleArray(iP*6+5)=myFBM%ParticleNew(iP+1)%Velocity(2)
      particleArray(iP*6+6)=myFBM%ParticleNew(iP+1)%Velocity(3)
      END DO

      ! send to the other processes
      DO pID=1,subnodes
      CALL SENDD_myMPI(particleArray,myFBM%nParticles*6,pID)  
      END DO
    ELSE
      ! non-master processes receive data from master
      CALL RECVD_myMPI(particleArray,myFBM%nParticles*6,0)
      DO iP = 0,myFBM%nParticles-1
      myFBM%ParticleNew(iP+1)%Position(1)=particleArray(iP*6+1)
      myFBM%ParticleNew(iP+1)%Position(2)=particleArray(iP*6+2)
      myFBM%ParticleNew(iP+1)%Position(3)=particleArray(iP*6+3)
      myFBM%ParticleNew(iP+1)%Velocity(1)=particleArray(iP*6+4)
      myFBM%ParticleNew(iP+1)%Velocity(2)=particleArray(iP*6+5)
      myFBM%ParticleNew(iP+1)%Velocity(3)=particleArray(iP*6+6)
      END DO 
    END IF
    ! end scatter  

    ! executed on the non-master processes
    ! copy particle data to c++ lib
    if(myid.ne.0)then
      DO iP = 1,myFBM%nParticles
      ipc=iP-1
      posx=myFBM%ParticleNew(iP)%Position(1)
      posy=myFBM%ParticleNew(iP)%Position(2)
      posz=myFBM%ParticleNew(iP)%Position(3)

      velx=myFBM%ParticleNew(iP)%Velocity(1)
      vely=myFBM%ParticleNew(iP)%Velocity(2)
      velz=myFBM%ParticleNew(iP)%Velocity(3)

      call setpositionid(posx,posy,posz,ipc)
      call setvelocityid(velx,vely,velz,ipc)
      END DO
    end if

  END SUBROUTINE FBM_ScatterParticles

  SUBROUTINE ExitError(msg,line)
    implicit none
    character(*) :: msg  
    !character(*) :: myfile  
    integer      :: line

    print *,msg
    write(*,*)'Fatal error in file:  Line: ',line
    stop

  END SUBROUTINE ExitError

end module cinterface
