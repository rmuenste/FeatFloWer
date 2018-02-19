    Subroutine ReadS3Dfile(cE3Dfile)
    use iniparser
    use Sigma_User

    implicit none

    character(len=*), intent(in) :: cE3Dfile
    logical :: bReadError=.FALSE.
    integer :: i,iSeg,iFile,iaux

    real*8 :: myPI = dATAN(1d0)*4d0
    character(len=INIP_STRLEN) cCut,cElement_i,cElemType,cKindOfConveying,cTemperature
    character(len=INIP_STRLEN) cProcessType,cRotation,cRheology,CDensity,cMeshQuality,cKTP
    
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    call inip_output_init(myid,showid,unitProtfile,unitTerminal)

    ! Init the parameterlist
    call inip_init(parameterlist)

!     READ (myFile,*) myProcess%Angle,myProcess%Phase

    call inip_readfromfile(parameterlist,ADJUSTL(TRIM(cE3Dfile)))

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","RotationDirection",cRotation)
    call inip_toupper_replace(cRotation)
    IF (ADJUSTL(TRIM(cRotation)).eq."R".OR.ADJUSTL(TRIM(cRotation)).eq."RECHTS".OR.ADJUSTL(TRIM(cRotation)).eq."RIGHT") myProcess%ind=1
    IF (ADJUSTL(TRIM(cRotation)).eq."LINKS".OR.ADJUSTL(TRIM(cRotation)).eq."LEFT".OR.ADJUSTL(TRIM(cRotation)).eq."L") myProcess%ind=-1
    myProcess%Rotation = cRotation
    IF (myProcess%ind.eq.0) THEN
     WRITE(*,*) "rotation direction is not defined"
     WRITE(*,*) '"',TRIM(myProcess%Rotation),'"'
    END IF
! 
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelDiameter", mySigma%Dz_out ,0d0/0d0)
    DistTolerance = 1d0*mySigma%Dz_Out
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","InnerDiameter", mySigma%Dz_in ,0d0/0d0)
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelLength", mySigma%L ,0d0/0d0)

    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfElements", mySigma%NumberOfSeg ,-1)
    
    IF (mySigma%NumberOfSeg.ge.1.and.mySigma%NumberOfSeg.le.9) THEN
     ALLOCATE (mySigma%mySegment(mySigma%NumberOfSeg))
    ELSE
     WRITE(*,*) "not a valid number of segments"
     WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
    ENDIF
    
    DO iSeg=1,mySigma%NumberOfSeg
     WRITE(cElement_i,'(A,I1.1)') 'E3DGeometryData/Machine/Element_',iSeg

     call INIP_getvalue_string(parameterlist,cElement_i,"ObjectType",mySigma%mySegment(iSeg)%ObjectType)
     call inip_toupper_replace(mySigma%mySegment(iSeg)%ObjectType)
     IF (.NOT.(mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW'.OR.mySigma%mySegment(iSeg)%ObjectType.eq.'DIE'.OR.mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE')) THEN
       WRITE(*,*) "STL object type is invalid. Only screw, die, or obstacle types are allowed"
     END IF

     call INIP_getvalue_string(parameterlist,cElement_i,"Type",cElemType)
     mySigma%mySegment(iSeg)%ART = ' '
     call inip_toupper_replace(cElemType)

     IF (ADJUSTL(TRIM(cElemType)).eq."STL") THEN
      mySigma%mySegment(iSeg)%ART   = "STL"
      call INIP_getvalue_string(parameterlist,cElement_i,"Unit",mySigma%mySegment(iSeg)%Unit,'CM')
      call inip_toupper_replace(mySigma%mySegment(iSeg)%Unit)
      IF (.NOT.(mySigma%mySegment(iSeg)%Unit.eq.'MM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'CM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'DM')) THEN
        WRITE(*,*) "STL unit type is invalid. Only MM, CM, DM units are allowed",mySigma%mySegment(iSeg)%Unit
      END IF

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition_[cm]",mySigma%mySegment(iSeg)%Min,0d0)
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,0d0/0d0)
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,0d0/0d0)
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL geometry dscription files are missing"
       WRITE(*,*) 'screwOFF'
       bReadError=.TRUE.
       GOTO 10
      END IF
      do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
        call INIP_getvalue_string(parameterlist,cElement_i,"screwOFF",mySigma%mySegment(iSeg)%OFFfiles(iFile),isubstring=iFile)
      end do
      mySigma%mySegment(iSeg)%Max = mySigma%mySegment(iSeg)%L + mySigma%mySegment(iSeg)%Min
     END IF
     IF (mySigma%mySegment(iSeg)%ART.eq.' ') THEN
      WRITE(*,*) "not a valid ",iSeg, "-segment"
      WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
      bReadError=.TRUE.
      GOTO 10
     ENDIF
    END DO

    
    myProcess%pTYPE = " "
    myProcess%dPress=0d0/0d0
    myProcess%Massestrom=0d0/0d0
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","ProcessType", cProcessType)
    call inip_toupper_replace(cProcessType)
    IF (ADJUSTL(TRIM(cProcessType)).eq."PRESSUREDROP") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","deltaP", myProcess%dPress,0d0/0d0)
     myProcess%pTYPE = "PRESSUREDROP"
    END IF
    IF (ADJUSTL(TRIM(cProcessType)).eq."THROUGHPUT") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","massthroughput", myProcess%Massestrom,0d0/0d0)
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MinInflowDiameter", myProcess%MinInflowDiameter,mySigma%Dz_In)
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaxInflowDiameter", myProcess%MaxInflowDiameter,mySigma%Dz_Out)
     myProcess%pTYPE = "THROUGHPUT"
    END IF
    IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq." ") THEN
     WRITE(*,*) "no valid process type is defined"
     WRITE(*,*) '"',TRIM(cProcessType),'"'
     bReadError=.TRUE.
     GOTO 10
    END IF
    
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewSpeed", myProcess%umdr,0d0/0d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaterialTemperature",myProcess%T0,0d0/0d0)

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","ScrewTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewTemperature",myProcess%Ti,0d0/0d0)
    ELSE
     myProcess%Ti=0d0/0d0
    END IF
    cTemperature=" "
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","BarrelTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","BarrelTemperature",myProcess%Ta,0d0/0d0)
    ELSE
     myProcess%Ta=0d0/0d0
    END IF
   
    myRheology%Equation = 0
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/RheologicalData","CalcVisco", cRheology)
    call inip_toupper_replace(cRheology)
    IF (ADJUSTL(TRIM(cRheology)).eq."CARREAU") THEN
      myRheology%Equation = 1
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","ZeroViscosity",myRheology%A,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","RecipVelocity",myRheology%B,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","Exponent",myRheology%C,0d0/0d0)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POWERLAW".OR.ADJUSTL(TRIM(cRheology)).eq."POTENZ".OR.ADJUSTL(TRIM(cRheology)).eq."POWER") THEN
      myRheology%Equation = 2
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Power","Consistence", myRheology%K,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Power","Exponent",myRheology%n,0d0/0d0)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POLYFLOW") THEN
      myRheology%Equation = 3
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_A",myRheology%A,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_B",myRheology%B,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_C",myRheology%C,0d0/0d0)
    END IF

    IF (myRheology%Equation.eq.0) THEN
     WRITE(*,*) "no valid rheology is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
     GOTO 10
    END IF
   
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material","LimitViscoMin",myRheology%ViscoMin,1d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material","LimitViscoMax",myRheology%ViscoMax,1d5)

    myRheology%AtFunc = 0
    cRheology = ' '
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/RheologicalData","CalcTemp", cRheology)
    call inip_toupper_replace(cRheology)
    
    IF (ADJUSTL(TRIM(cRheology)).eq."ISOTHERM") THEN
      myRheology%AtFunc = 1
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."C1C2") THEN
      myRheology%AtFunc = 2
      myRheology%Ts = 165d0
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/C1C2","C1",myRheology%C1,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/C1C2","C2",myRheology%C2,0d0/0d0)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."TBTS") THEN
      myRheology%AtFunc = 3
      myRheology%C1 = 8.86d0
      myRheology%C2 = 101.6d0
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/TBTS","ReferenceTemperature",myRheology%Tb,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/TBTS","StandardTemperature",myRheology%Ts,0d0/0d0)
    END IF
    IF (myRheology%AtFunc.eq.0) THEN
     WRITE(*,*) "no temperature correction is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
     GOTO 10
    END IF

    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatConductivity",myThermodyn%lambda,0d0/0d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatConductivitySlope",myThermodyn%lambdaSteig,0d0/0d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatCapacity",myThermodyn%cp,0d0/0d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatCapacitySlope",myThermodyn%CpSteig,0d0/0d0)

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/ThermoData","DensityModel", cDensity)
    call inip_toupper_replace(cDensity)
    myThermodyn%density=0d0/0d0
    IF (ADJUSTL(TRIM(cDensity)).eq."DENSITY") THEN
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/Density","Density",myThermodyn%Density,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/Density","DensitySlope",myThermodyn%DensitySteig,0d0/0d0)
    END IF
    IF (ADJUSTL(TRIM(cDensity)).eq."SPECVOLUME") THEN
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/SpecVolume","SpecVolume",myThermodyn%Density,0d0/0d0)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/SpecVolume","SpecVolumeSlope",myThermodyn%DensitySteig,0d0/0d0)
      myThermodyn%Density = 1d0/myThermodyn%Density
    END IF
    IF (myRheology%AtFunc.eq.0) THEN
     WRITE(*,*) "density is not defined"
     WRITE(*,*) '"',TRIM(cDensity),'"'
     bReadError=.TRUE.
     GOTO 10
    END IF

    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output","NumberOf1DLayers",myOutput%nOf1DLayers,16)
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output","NumberOfHistRegions",myOutput%nRegion,4)
    
    cKTP=' '
    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","KTPRelease",cKTP,"YES")
    call inip_toupper_replace(cKTP)
    IF (ADJUSTL(TRIM(cKTP)).eq."NO") THEN
     bKTPRelease = .FALSE.
    END IF

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","HexMesher", mySetup%cMesher,"OFF")
    call inip_toupper_replace(mySetup%cMesher)

    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."OFF") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","MeshPath", mySetup%cMeshPath)
    END IF
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."HOLLOWCYLINDER") THEN
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Tangential",mySetup%m_nT,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Radial",mySetup%m_nR,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Axial",mySetup%m_nZ,0)
     IF (mySetup%m_nT.eq.0.or.mySetup%m_nR.eq.0.or.mySetup%m_nZ.eq.0) THEN
      WRITE(*,*) "mesh resolution is not correctly defined"
      WRITE(*,*) '"',mySetup%m_nT,mySetup%m_nR,mySetup%m_nZ,'"'
      bReadError=.TRUE.
      GOTO 10
     END IF
    END IF


!     cMeshQuality=' '
!     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","SendEmail",cMeshQuality,"YES")
!     call inip_toupper_replace(cMeshQuality)
!     IF (ADJUSTL(TRIM(cMeshQuality)).eq."NO") THEN
!      mySetup%bSendEmail = .FALSE.
!     END IF
!     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","Periodicity",myProcess%Periodicity,1)
!     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nSolutions",mySetup%nSolutions,1)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","dAlpha",myProcess%dAlpha,10d0)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","Angle",myProcess%Angle,0d0/0d0)
!     call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","Phase",myProcess%Phase,0d0/0d0)
    
    
    IF (myid.eq.1) then
    write(*,*) "=========================================================================="
    write(*,*) "mySigma%Dz_Out",'=',mySigma%Dz_out
    write(*,*) "mySigma%Dz_In",'=',mySigma%Dz_In
    write(*,*) "mySigma%L",'=',mySigma%L
    
    write(*,*) 
    DO iSeg=1,mySigma%NumberOfSeg
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Art=',mySigma%mySegment(iSeg)%ART
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%ObjectType=',mySigma%mySegment(iSeg)%ObjectType
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Unit=',mySigma%mySegment(iSeg)%Unit
     write(*,'(A,I1.1,A,f12.2)') " mySIGMA%Segment(",iSeg,')%Min=',mySigma%mySegment(iSeg)%Min
     write(*,'(A,I1.1,A,f12.2)') " mySIGMA%Segment(",iSeg,')%Max=',mySigma%mySegment(iSeg)%Max
     write(*,'(A,I1.1,A,f12.2)') " mySIGMA%Segment(",iSeg,')%L=',mySigma%mySegment(iSeg)%L
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN
      write(*,'(A,I1.1,A,f12.2)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I1.1,A,I)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     write(*,*) 
    END DO    

    write(*,*) "myProcess%Rotation",'=',myProcess%Rotation
    write(*,*) "myProcess%ind",'=',myProcess%ind
    write(*,*) "myProcess%deltaP",'=',myProcess%dPress
    write(*,*) "myProcess%Massestrom",'=',myProcess%Massestrom
    write(*,*) "myProcess%f",'=',myProcess%umdr
    write(*,*) "myProcess%Ti",'=',myProcess%Ti
    write(*,*) "myProcess%Ta",'=',myProcess%Ta
    write(*,*) "myProcess%T0",'=',myProcess%T0
    write(*,*) 
    write(*,*) "myRheology%ViscoMin",'=',myRheology%ViscoMin
    write(*,*) "myRheology%ViscoMax",'=',myRheology%ViscoMax
    IF (myRheology%Equation.eq.2) THEN
     write(*,*) "myRheology%model",'=','Powerlaw'
     write(*,*) "myRheology%K",'=',myRheology%K
     write(*,*) "myRheology%n",'=',myRheology%n
    END IF
    IF (myRheology%Equation.eq.1) THEN
     write(*,*) "myRheology%model",'=','Carreau'
     write(*,*) "myRheology%A",'=',myRheology%A
     write(*,*) "myRheology%B",'=',myRheology%B
     write(*,*) "myRheology%C",'=',myRheology%C
    END IF
    IF (myRheology%Equation.eq.3) THEN
     write(*,*) "myRheology%model",'=','Polyflow'
     write(*,*) "myRheology%A",'=',myRheology%A
     write(*,*) "myRheology%B",'=',myRheology%B
     write(*,*) "myRheology%C",'=',myRheology%C
    END IF
    write(*,*) 
    IF (myRheology%AtFunc.eq.1) THEN
     write(*,*) "myRheology%TempModel",'=','ISOTHERM'
     write(*,*) "myRheology%aT",'=',1.0
    END IF
    IF (myRheology%AtFunc.eq.2) THEN
     write(*,*) "myRheology%TempModel",'=','C1C2'
     write(*,*) "myRheology%C1",'=',myRheology%C1
     write(*,*) "myRheology%C2",'=',myRheology%C2
     write(*,*) "myRheology%Ts",'=',myRheology%Ts
    END IF
    IF (myRheology%AtFunc.eq.3) THEN
     write(*,*) "myRheology%TempModel",'=','TBTS'
     write(*,*) "myRheology%C1",'=',myRheology%C1
     write(*,*) "myRheology%C2",'=',myRheology%C2
     write(*,*) "myRheology%Tb",'=',myRheology%Tb
     write(*,*) "myRheology%Ts",'=',myRheology%Ts
    END IF
    write(*,*) 
    write(*,*) "myThermodyn%HeatConductivity",'=',myThermodyn%lambda
    write(*,*) "myThermodyn%HeatConductivitySlope",'=',myThermodyn%lambdaSteig
    write(*,*) "myThermodyn%HeatCapacity",'=',myThermodyn%cp
    write(*,*) "myThermodyn%HeatCapacitySlope",'=',myThermodyn%cpSteig
    write(*,*) "myThermodyn%Density",'=',myThermodyn%density
    write(*,*) "myThermodyn%DensitySlope",'=',myThermodyn%densitySteig
    write(*,*) 
!     write(*,*) "mySetup%MeshQuality",'=',mySetup%MeshResolution
    write(*,*) "muOutput%NumberOf1DLayers = ",myOutput%nOf1DLayers
    write(*,*) "muOutput%NumberOfHistRegions = ",myOutput%nRegion
    write(*,*) 
   IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."OFF") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,*) "mySetup%MeshPath",'=',ADJUSTL(TRIM(mySetup%cMeshPath))
    END IF
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."HOLLOWCYLINDER") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,*) "mySetup%Resolution[nR,nT,nZ]",'=',mySetup%m_nR,mySetup%m_nT,mySetup%m_nZ
    END IF
    write(*,*)
    write(*,*) "myProcess%dAlpha",'=',myProcess%dAlpha
!     write(*,*) "mySetup%bSendEmail",'=',mySetup%bSendEmail
    write(*,*) "myProcess%Angle",'=',myProcess%Angle
!     write(*,*) "myProcess%Phase",'=',myProcess%Phase
!     write(*,*) "=========================================================================="
    END IF
! 
!     myThermodyn%Alpha     = (1e6*myThermodyn%lambda)/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
!     myThermodyn%Beta      = 1d1 ! !myProcess%Cnst_Lam/(myProcess%Cnst_Dens*myProcess%Cnst_Cp)
!     myThermodyn%Gamma     = 1d0/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
! 
!     IF (myid.eq.1) then
!      WRITE(*,*) 
!      WRITE(*,*) 'myThermodyn%Alpha = ', myThermodyn%Alpha
!      WRITE(*,*) 'myThermodyn%Beta = ',  myThermodyn%Beta
!      WRITE(*,*) 'myThermodyn%Gamma = ', myThermodyn%Gamma
!     end if
! 
    IF (.NOT.ISNAN(myProcess%dPress)) THEN
    dZPeriodicLenght = mySigma%L
    myProcess%dPress = 1d3*myProcess%dPress
    bNoOutflow = .TRUE.
    ELSE
    bNoOutflow = .FALSE.
    dZPeriodicLenght = 1d5*mySigma%L
    END IF

    IF (myid.eq.1) then
     write(*,*) "myProcess%dZPeriodicLenght",'=',dZPeriodicLenght
     write(*,*) "myProcess%dPress",'=',myProcess%dPress
     write(*,*) "myProcess%NoOutFlow",'=',bNoOutFlow
     write(*,*) "=========================================================================="
    END IF
    10  CONTINUE

!    ! Make some output on the terminal
!     call inip_info(parameterlist)

!     ! Write it into a file
!      call inip_dumpToFile(parameterlist,"test_dump_parlst.dat",INIP_REPLACE)

    ! Clean up the parameterlist
    call inip_done(parameterlist)

end Subroutine ReadS3Dfile
