module particles_input

use iniparser
use types
USE PP3D_MPI, ONLY:myid,showid

  implicit none
  
  ! Top-level file that will be read. Path is relative
  ! to the binary
  character(len=256), parameter, private :: configfile_particles = "_data/config_particles.dat"

contains
  subroutine prt_read_config(ParticleParam, mfile, mterm)
    type(tParticleParam), intent(inout) :: ParticleParam
    ! Output units for file-output and terminal output
    integer, intent(in) :: mfile
    integer, intent(in) :: mterm
    character(len=INIP_STRLEN) :: particleDatafile, machineDatafile

    ! Local variables
    ! A parameterlist
    type(t_parlist) :: parameterlist

    character(len=256) :: tmpstring,cSeg,cInitType,cInitTypeUpper
    integer :: i, numberOfElements,ierr

    ! Declare default parameters
    real*8 :: minFrac, dEps1, dEps2,Z_seed,Epsilon,Z1,Z2,hSize
    minFrac = 0.95
    dEps1 = 0.9995
    dEps2 = 0.9990
    Z_seed = 0.00
    Z1 = 0.00
    Z2 = 0.00
    hSize = 0.1
    Epsilon = 5d-2
    
    ! Init the output from the iniparser
    call inip_output_init(myid,showid,mfile,mterm)

    ! Init the parameterlist
    call inip_init(parameterlist)

    ! Now read the configfile for the particles
    call inip_readfromfile(parameterlist,trim(adjustl(configfile_particles)))

    ! Now lets search for parameters
    ! First find out the other datafiles which we are going to read as well
    call inip_getvalue_string(parameterlist,"GeneralSettings","particle_datafile",particleDatafile,bdequote=.TRUE.)
    call inip_getvalue_string(parameterlist,"GeneralSettings","machine_datafile",machineDatafile,bdequote=.TRUE.)

    ! Now lets read these files as well
    call inip_readfromfile(parameterlist,trim(adjustl(particleDatafile)))
    call inip_readfromfile(parameterlist,trim(adjustl(machineDatafile)))

    ! Control-Output to the terminal and in a file
    call inip_info(parameterlist)
    call inip_dumpToFile(parameterlist,"_data/particleparameters_controloutput.txt",INIP_REPLACE)

    ! First thing (as most important one) is to get the units
    ! Example: If Z_Seed = 2.0 cm, we can leave it as it is.
    ! If it is 2.0mm, we have to multiply by 0.1 because internally the code works in cm!
    call inip_getvalue_string(parameterlist,"GeneralSettings","unitIN",tmpstring,bdequote=.TRUE.)
    call inip_toupper_replace(tmpstring)
    if (trim(adjustl(tmpstring)) .eq. "CM") then
      ParticleParam%dFacUnitIn = 1.0d0
    else if (trim(adjustl(tmpstring)) .eq. "MM") then
      ParticleParam%dFacUnitIn = 0.1d0
    else
      write(*,*) "Unknown unit: ", tmpstring, "Stopping"
      stop
    end if

    ! Example:
    ! If the particle is at z = 7cm (cm is the internal unit of the code!)
    ! and we want to output in cm, we can leave it as it is.
    ! If we want to output it in mm, we have to multiply by 10
    call inip_getvalue_string(parameterlist,"GeneralSettings","unitOUT",tmpstring,bdequote=.TRUE.)
    call inip_toupper_replace(tmpstring)
    if (trim(adjustl(tmpstring)) .eq. "CM") then
      ParticleParam%dFacUnitOut = 1.0d0
    else if (trim(adjustl(tmpstring)) .eq. "MM") then
      ParticleParam%dFacUnitOut = 10.0d0
    else
      write(*,*) "Unknown unit: ", tmpstring, "Stopping"
      stop
    end if

    ParticleParam%bRotationalMovement=.false.
    call inip_getvalue_string(parameterlist,"GeneralSettings","RotationalMovement",tmpstring,"YES")
    call inip_toupper_replace(tmpstring)
    IF (tmpstring.eq."NO") THEN
     ParticleParam%bRotationalMovement = .FALSE.
    ELSE
     ParticleParam%bRotationalMovement = .TRUE.
    END IF

    IF (ParticleParam%bRotationalMovement) THEN
     call inip_getvalue_double(parameterlist, "GeneralSettings","f",ParticleParam%f)
    ELSE
     call inip_getvalue_double(parameterlist, "GeneralSettings","dt",ParticleParam%f,0d0)
     IF (ParticleParam%f.eq.0d0) THEN
      write(*,*) "No TimeStep size defined: ", ParticleParam%f
      stop
     ELSE
      ParticleParam%f = 1d0/ParticleParam%f
     END IF
    END IF
    
    call inip_getvalue_string(parameterlist,"GeneralSettings","DumpFormat",tmpstring,"DMP")
    call inip_toupper_replace(tmpstring)
    IF (tmpstring.eq."DMP") ParticleParam%DumpFormat = 1
    IF (tmpstring.eq."LST") ParticleParam%DumpFormat = 2
    IF (tmpstring.eq."REPART") ParticleParam%DumpFormat = 3

    ! Get the number of Inflow regions for Particle Backtracing
    call inip_getvalue_int(parameterlist,"GeneralSettings","NumberOfInflowRegions",ParticleParam%NumberOfInflowRegions,0)
    IF (ParticleParam%NumberOfInflowRegions.eq.0) THEN
     ParticleParam%NumberOfInflowRegions = 1
     GOTO 56
    ELSE
     allocate(ParticleParam%InflowRegion(ParticleParam%NumberOfInflowRegions))
     DO i=1,ParticleParam%NumberOfInflowRegions
       WRITE(cSeg,'(A,I1.1)') 'InflowRegion_',i
       call inip_getvalue_string(parameterlist,"GeneralSettings",cSeg,tmpstring,'unknown')
       if (myid.eq.1) write(*,*) "'",adjustl(trim(tmpstring)),"'"
       read(tmpstring,*,err=55) ParticleParam%InflowRegion(i)%Center,ParticleParam%InflowRegion(i)%Radius
     END DO
     GOTO 56
    END IF
    
55     write(*,*) 'WRONGLY DEFINED parameters for Inflow ',i,' !!'
56  CONTINUE

    ! Get the starting procedure
    call inip_getvalue_string(parameterlist,"GeneralSettings","startingprocedure",cInitType,bdequote=.TRUE.)
    cInitTypeUpper = cInitType
    call inip_toupper_replace(cInitTypeUpper)
    IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'CSV').ne.0) THEN
      ParticleParam%inittype  = ParticleSeed_CSVFILE
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'FILE=').ne.0) THEN
       ParticleParam%Plane = 1
       READ(cInitType(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'FILE=')+5:),*) ParticleParam%sourcefile
       if (myid.eq.1) WRITE(*,*) "Particle sourcefile: '",ADJUSTL(TRIM(ParticleParam%sourcefile)),"'"
      END IF
    ELSEIF (ADJUSTL(TRIM(cInitTypeUpper)).eq.'ELEM') THEN
      ParticleParam%inittype  = ParticleSeed_ELEMCENTER
    ELSEIF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'PLANE').ne.0) THEN
      ParticleParam%inittype  = ParticleSeed_PLANE
      ParticleParam%Plane = 0
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'X=').ne.0) THEN
       ParticleParam%Plane = 1
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'X=')+2:),*) ParticleParam%PlaneOffset
      END IF
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'Y=').ne.0) THEN
       ParticleParam%Plane = 2
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'Y=')+2:),*) ParticleParam%PlaneOffset
      END IF
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'Z=').ne.0) THEN
       ParticleParam%Plane = 3
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'Z=')+2:),*) ParticleParam%PlaneOffset
      END IF
      IF (ParticleParam%Plane.eq.0) THEN
       WRITE(*,*) 'no plane for planar seedoing is defined'
       STOP
      END IF
      
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=').ne.0) THEN
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=')+2:),*) ParticleParam%PlaneParticles
      ELSE
       WRITE(*,*) 'no numbner of particle is defined ...'
       STOP
      END IF
      
      WRITE(*,*) ParticleSeed_PLANE,ParticleParam%PlaneParticles, ParticleParam%Plane ,ParticleParam%PlaneOffset
    ELSEIF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'VOLUME').ne.0) THEN
      ParticleParam%inittype  = ParticleSeed_VOLUME
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=').ne.0) THEN
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=')+2:),*) ParticleParam%VolumeParticles
      ELSE
       WRITE(*,*) 'no numbner of particle is defined ...'
       STOP
      END IF
      WRITE(*,*) ParticleSeed_VOLUME,ParticleParam%VolumeParticles
      
    ELSEIF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'E3D').ne.0) THEN
      ParticleParam%inittype  = ParticleSeed_E3D
      
      IF (INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=').ne.0) THEN
       READ(cInitTypeUpper(INDEX(ADJUSTL(TRIM(cInitTypeUpper)),'N=')+2:),*) ParticleParam%PlaneParticles
      ELSE
       WRITE(*,*) 'no numbner of particle is defined ...'
       STOP
      END IF
      WRITE(*,*) ParticleSeed_E3D,ParticleParam%PlaneParticles

    ELSE
    
       WRITE(*,*) 'no valid particle seed is defined ...'
       STOP
     
    END IF
    
    
!    call inip_getvalue_int(parameterlist,"GeneralSettings","startingprocedure",ParticleParam%inittype)

    ! If we have to read something then lets do it
    if (ParticleParam%inittype .eq. ParticleSeed_CSVFILE .or. &
        ParticleParam%inittype .eq. ParticleSeed_OUTPUTFILE) then
!       call inip_getvalue_string(parameterlist,"GeneralSettings","sourcefile",ParticleParam%sourcefile,bdequote=.TRUE.)

      ! Get the unit of the sourcefile
      ! Example: If we read that the particle has the position 2.5cm , we can leave it as it is.
      ! If it is 2.5mm, we have to multiply by 0.1 because internally the code works in cm!
      call inip_getvalue_string(parameterlist,"GeneralSettings","unitSourcefile",tmpstring,bdequote=.TRUE.)
      call inip_toupper_replace(tmpstring)
      if (trim(adjustl(tmpstring)) .eq. "CM") then
        ParticleParam%dFacUnitSourcefile = 1.0d0
      else if (trim(adjustl(tmpstring)) .eq. "MM") then
        ParticleParam%dFacUnitSourcefile = 0.1d0
      else
        write(*,*) "Unknown unit: ", tmpstring, "Stopping"
        stop
      end if
    end if

    call inip_getvalue_string(parameterlist,"GeneralSettings","backtrace",tmpstring,"NO")
    call inip_toupper_replace(tmpstring)
    IF (tmpstring.eq."NO") THEN
     ParticleParam%bBacktrace = .FALSE.
    ELSE
     ParticleParam%bBacktrace = .TRUE.
    END IF

    call inip_getvalue_string(parameterlist,"GeneralSettings","PhysicalParticles",tmpstring,"NO")
    call inip_toupper_replace(tmpstring)
    IF (tmpstring.eq."NO") THEN
     ParticleParam%bPhysParticles = .FALSE.
    ELSE
     ParticleParam%bPhysParticles = .TRUE.
    END IF
    
    IF (ParticleParam%bPhysParticles) THEN
     call inip_getvalue_double(parameterlist,"GeneralSettings/PhysicalParticles","rho_l_gcm3",ParticleParam%PhysParticles%rho_l,-1d0)
     if (ParticleParam%PhysParticles%rho_l.eq.-1d0) GOTO 65
     call inip_getvalue_double(parameterlist,"GeneralSettings/PhysicalParticles","mu_l_Pas",ParticleParam%PhysParticles%mu_l,-1d0)
     if (ParticleParam%PhysParticles%mu_l.eq.-1d0) GOTO 65
     ParticleParam%PhysParticles%mu_l = 1d1*ParticleParam%PhysParticles%mu_l
     call inip_getvalue_double(parameterlist,"GeneralSettings/PhysicalParticles","rho_p_gcm3",ParticleParam%PhysParticles%rho_p,-1d0)
     if (ParticleParam%PhysParticles%rho_p.eq.-1d0) GOTO 65
     call inip_getvalue_double(parameterlist,"GeneralSettings/PhysicalParticles","d_p_mm",ParticleParam%PhysParticles%d_p,-1d0)
     if (ParticleParam%PhysParticles%d_p.eq.-1d0) GOTO 65
     ParticleParam%PhysParticles%d_p = 0.1d0*ParticleParam%PhysParticles%d_p
     call inip_getvalue_string(parameterlist,"GeneralSettings/PhysicalParticles","gravity_ms2",tmpstring,"0d0,0d0,0d0")
     READ(tmpstring,*,err=65) ParticleParam%PhysParticles%gravity
     ParticleParam%PhysParticles%gravity = 100d0*ParticleParam%PhysParticles%gravity
     GOTO 66
65   write(*,*) 'WRONGLY DEFINED parameter for physical particles model !!' 
     STOP
66   CONTINUE
    END IF
    
    
    call inip_getvalue_int(parameterlist,"GeneralSettings","dump_in_file",ParticleParam%dump_in_file,-1)
    
    ! TimeLevels. Default is 72
    call inip_getvalue_int(parameterlist,"GeneralSettings","TimeLevels",ParticleParam%nTimeLevels,72)

    call inip_getvalue_int(parameterlist,"GeneralSettings","OutputFreq",ParticleParam%OutputFreq,1)
    
    call inip_getvalue_int(parameterlist,"GeneralSettings","Periodicity",ParticleParam%nPeriodicity,1)
    

    ! The number of particles. We only need them if we seed from the parameterfile
    ! Default is 32000
    if (ParticleParam%inittype .eq. ParticleSeed_Parameterfile) then
      call inip_getvalue_int(parameterlist, "GeneralSettings","nElements",numberOfElements)
      ParticleParam%nParticles = 4*numberOfElements
      call inip_getvalue_double(parameterlist,"GeneralSettings","nSegmentsColoring",ParticleParam%numberSegments)
    end if

    ! Maximum number of rotations that has to be performed
    call inip_getvalue_int(parameterlist, "GeneralSettings","nRotation",ParticleParam%nRotation)

    call inip_getvalue_double(parameterlist, "GeneralSettings","d_CorrDist",ParticleParam%d_CorrDist,1d-2)

    call inip_getvalue_double(parameterlist, "GeneralSettings","minFrac",ParticleParam%minFrac,minFrac)
    call inip_getvalue_int(parameterlist, "GeneralSettings","Raster",ParticleParam%Raster,50)
    call inip_getvalue_double(parameterlist, "GeneralSettings","dEps1",ParticleParam%dEps1,dEps1)
    call inip_getvalue_double(parameterlist, "GeneralSettings","dEps2",ParticleParam%dEps2,dEps2)

    ! Outer diameter for the particle seed
    call inip_getvalue_double(parameterlist, "GeneralSettings","D_Out",ParticleParam%D_Out)
    ! Scale it with the unit factor
    ParticleParam%D_Out=ParticleParam%D_Out*ParticleParam%dFacUnitIn
    call inip_getvalue_double(parameterlist, "GeneralSettings","D_In",ParticleParam%D_In)
    ! Scale it with the unit factor
    ParticleParam%D_In=ParticleParam%D_In*ParticleParam%dFacUnitIn
    call inip_getvalue_double(parameterlist, "GeneralSettings","Z_Seed",ParticleParam%Z_seed,Z_seed)
    ! Scale it with the unit factor
    ParticleParam%Z_seed=ParticleParam%Z_seed*ParticleParam%dFacUnitIn

    call inip_getvalue_double(parameterlist, "GeneralSettings","Epsilon",ParticleParam%Epsilon,Epsilon)
    call inip_getvalue_double(parameterlist, "GeneralSettings","hSize",ParticleParam%hSize,hSize)

    call inip_getvalue_double(parameterlist, "GeneralSettings","OutflowZPosition",ParticleParam%OutflowZPos,1d8)

    ! Particles at positions - new implementation
    ParticleParam%nZposCutplanes = INIP_querysubstrings(parameterlist,"GeneralSettings","zCutplanePositions")
    if (ParticleParam%nZposCutplanes .gt. 0 .and. ParticleParam%nZposCutplanes .lt. 100 ) then
      allocate(ParticleParam%cutplanePositions(ParticleParam%nZposCutplanes))
      do i=1,ParticleParam%nZposCutplanes
        call inip_getvalue_double(parameterlist,"GeneralSettings","zCutplanePositions",&
                                  ParticleParam%cutplanePositions(i),iarrayindex=i)
        ParticleParam%cutplanePositions(i) = ParticleParam%cutplanePositions(i)*ParticleParam%dFacUnitIn
      end do
    ! Don't do more than 99 cutplanes - this will give inconsistent filenames
    ! and be super-slow
    else if (ParticleParam%nZposCutplanes .ge. 100) then
      write(*,*) "More than 99 Cutplanes is really something you do not want to go for. This is too slow"
      write(*,*) "Stopping"
      stop
    end if

    call inip_done(parameterlist)

  end subroutine

end module
