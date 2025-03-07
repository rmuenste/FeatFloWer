 Subroutine ReadS3Dfile(cE3Dfile)
        
    use iniparser
    USE PP3D_MPI, ONLY:myid,showid,subnodes,dZPeriodicLength,dPeriodicity

    use Sigma_User, only : mySigma,myThermodyn,mySetup,myOutput,myTransientSolution,&
        myProcess,myMultiMat,SoftwareRelease,bKTPRelease,DistTolerance
    USE var_QuadScalar

    use, intrinsic :: ieee_arithmetic

    implicit none

    character(len=*), intent(in) :: cE3Dfile
    logical :: bReadError=.FALSE.
    integer :: i,iSeg,iFile,iaux,iInflow,iInflowErr,iMat,ierr,iSubInflow,iTempBC

    real*8 :: myPI = dATAN(1d0)*4d0
    character(len=INIP_STRLEN) cCut,cElement_i,cElemType,cKindOfConveying,cTemperature,cPressureFBM
    character(len=INIP_STRLEN) cBCtype,cInflow_i,cCenter,cNormal,cauxD,cauxZ,cOnlyBarrelAdaptation,cVelo,cTempBC_i
    character(len=INIP_STRLEN) cMidpointA, cMidpointB 
    character(len=INIP_STRLEN) cParserString,cSCR,cALE,cDissip,cSensor,cTVals, cRobin

    character(len=INIP_STRLEN) cProcessType,cRotation,cRheology,cMeshQuality,cKTP,cUnit,cOFF_Files,cShearRateRest,cTXT
    
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    real*8 :: myInf,dSizeScale,daux,dElemSizeScale,dFlow,dArea

    real*8 dExtract_Val
    character(len=INIP_STRLEN) cText,sExtract_Dim


    if(ieee_support_inf(myInf))then
      myInf = ieee_value(myInf, ieee_negative_inf)
    endif


    call inip_output_init(myid,showid,unitProtfile,unitTerminal)

    ! Init the parameterlist
    call inip_init(parameterlist)

!     READ (myFile,*) myProcess%Angle,myProcess%Phase

    call inip_readfromfile(parameterlist,ADJUSTL(TRIM(cE3Dfile)))

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","RotationType",cRotation,'co')
    call inip_toupper_replace(cRotation)
    myProcess%iInd =+1
    IF (ADJUSTL(TRIM(cRotation)).eq."COUNTER".OR.ADJUSTL(TRIM(cRotation)).eq."COUNTERROTATING") myProcess%iInd=-1
    
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","RotationDirection",cRotation,'NoRotation')
    call inip_toupper_replace(cRotation)
    IF (ADJUSTL(TRIM(cRotation)).eq."R".OR.ADJUSTL(TRIM(cRotation)).eq."RECHTS".OR.ADJUSTL(TRIM(cRotation)).eq."RIGHT") myProcess%ind=1
    IF (ADJUSTL(TRIM(cRotation)).eq."LINKS".OR.ADJUSTL(TRIM(cRotation)).eq."LEFT".OR.ADJUSTL(TRIM(cRotation)).eq."L") myProcess%ind=-1
    myProcess%Rotation = cRotation
    IF (myProcess%ind.eq.0) THEN
     WRITE(*,*) "rotation direction is not defined"
     WRITE(*,*) '"',TRIM(myProcess%Rotation),'"'
    END IF
! 
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","Unit",cUnit,'MM')
    call inip_toupper_replace(cUnit)
    IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
      WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
      cUnit = 'MM'
    END IF
    if (TRIM(cUnit).eq.'MM') dSizeScale = 0.100d0
    if (TRIM(cUnit).eq.'CM') dSizeScale = 1.000d0
    if (TRIM(cUnit).eq.'DM') dSizeScale = 10.00d0
    if (TRIM(cUnit).eq.'M')  dSizeScale = 100.0d0

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","ScrewCylinderRendering",cSCR,'YES')
    call inip_toupper_replace(cSCR)
    IF (ADJUSTL(TRIM(cSCR)).EQ."YES".OR.ADJUSTL(TRIM(cSCR)).EQ."ON") THEN
     mySigma%ScrewCylinderRendering=.true.
    END IF
    IF (ADJUSTL(TRIM(cSCR)).EQ."NO".OR.ADJUSTL(TRIM(cSCR)).EQ."OFF") THEN
     mySigma%ScrewCylinderRendering=.false.
    END IF
    
    cSCR=" "
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","BarrelRendering",cSCR,'OFF')
    call inip_toupper_replace(cSCR)
    IF (ADJUSTL(TRIM(cSCR)).EQ."YES".OR.ADJUSTL(TRIM(cSCR)).EQ."ON") THEN
     mySigma%BarrelRendering=.true.
    END IF
    IF (ADJUSTL(TRIM(cSCR)).EQ."NO".OR.ADJUSTL(TRIM(cSCR)).EQ."OFF") THEN
     mySigma%BarrelRendering=.false.
    END IF

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","Type",mySigma%cType,'SSE')
    call inip_toupper_replace(mySigma%cType)
    IF (.NOT.(ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.&
              ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".OR.&
              ADJUSTL(TRIM(mySigma%cType)).EQ."DIE".OR.&
              ADJUSTL(TRIM(mySigma%cType)).EQ."XSE".OR.&
              ADJUSTL(TRIM(mySigma%cType)).EQ."NETZSCH")) THEN
     WRITE(*,*) "not a valid Extruder type:", ADJUSTL(TRIM(mySigma%cType))
     CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
    END IF
    
!     WRITE(*,*) "asdsadsa sad sa dsad as as :",ADJUSTL(TRIM(mySigma%cType))
    
    IF (SoftwareRelease.ne."XSE") THEN
    IF (SoftwareRelease.NE.ADJUSTL(TRIM(mySigma%cType))) THEN
     WRITE(*,*) "This is not a valid 'TSE' (Twinscrew Extrusion) Simulation setup."
     WRITE(*,*) "The configured '"//ADJUSTL(TRIM(mySigma%cType))//"' is not supported."
     WRITE(*,*) "Program stops!"
     CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
    END IF
    END IF

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","GeometryStart", mySigma%GeometryStart ,'_INVALID_')
    call inip_toupper_replace(mySigma%GeometryStart)
    IF (adjustl(trim(mySigma%GeometryStart)).ne.'_INVALID_') THEN
     Read(mySigma%GeometryStart,*,err=550) mySigma%DIE_Start
     mySigma%DIE_Start = dSizeScale*mySigma%DIE_Start
     GOTO 551
550  CONTINUE
     if (myid.eq.1) WRITE(*,*) "   Wrong DIE start has been set: ",adjustl(trim(mySigma%GeometryStart))
     CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
551  CONTINUE
    END IF
    
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","GeometryLength", mySigma%GeometryLength ,'_INVALID_')
    call inip_toupper_replace(mySigma%GeometryLength)
    IF (adjustl(trim(mySigma%GeometryLength)).ne.'_INVALID_') THEN
     Read(mySigma%GeometryLength,*,err=552) mySigma%DIE_Length
     mySigma%DIE_Length = dSizeScale*mySigma%DIE_Length
     GOTO 553
552  CONTINUE
     if (myid.eq.1) WRITE(*,*) "   Wrong DIE Length has been set: ",adjustl(trim(mySigma%GeometryLength))
     CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
553  CONTINUE
    END IF
   
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","GeometrySymmetryBC", mySigma%GeometrySymmetryBC ,'_INVALID_')
    call inip_toupper_replace(mySigma%GeometrySymmetryBC)
    IF (adjustl(trim(mySigma%GeometrySymmetryBC)).ne.'_INVALID_') THEN
     Read(mySigma%GeometrySymmetryBC,*,err=554) mySigma%DIE_SymmetryBC
     GOTO 555
554  CONTINUE
     if (myid.eq.1) WRITE(*,*) "   Wrong DIE Symmetry BC has been set: ",adjustl(trim(mySigma%GeometrySymmetryBC))
     CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
555  CONTINUE
    END IF

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelDiameter", mySigma%Dz_out ,myInf)
    mySigma%Dz_out = dSizeScale*mySigma%Dz_out
    DistTolerance = 1d0*mySigma%Dz_Out

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelStraightCut", mySigma%W ,myInf)
    mySigma%W = dSizeScale*mySigma%W
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelCurvedCut", mySigma%Dzz ,myInf)
    mySigma%Dzz = dSizeScale*mySigma%Dzz
    
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","BarrelCut", mySigma%cZwickel ,"STRAIGHT")
    call inip_toupper_replace(mySigma%cZwickel)
    IF (.NOT.(ADJUSTL(TRIM(mySigma%cZwickel)).EQ."STRAIGHT".OR.ADJUSTL(TRIM(mySigma%cZwickel)).EQ."ROUND").OR.ADJUSTL(TRIM(mySigma%cZwickel)).EQ."CURVED") THEN
     WRITE(*,*) "not a valid Zwickel region definition:", ADJUSTL(TRIM(mySigma%cZwickel))
    END IF
    IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."STRAIGHT") mySigma%Dzz = myInf
    IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."ROUND".OR.ADJUSTL(TRIM(mySigma%cZwickel)).EQ."CURVED") THEN
     mySigma%cZwickel = "ROUND"
     mySigma%W = myInf
    END IF

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","InnerDiameter", mySigma%Dz_in ,mySigma%Dz_Out/dSizeScale)
    mySigma%Dz_in = dSizeScale*mySigma%Dz_in
    
    call INIP_getvalue_Int(parameterlist,"E3DGeometryData/Machine","InnerDiameterNList", mySigma%InnerDiamNParam ,0)
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","InnerDiameterAxiList", cauxD ,"no")
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","InnerDiameterRadList", cauxZ ,"no")
    allocate(mySigma%InnerDiamDParam(mySigma%InnerDiamNParam))
    allocate(mySigma%InnerDiamZParam(mySigma%InnerDiamNParam))
    read(cauxD,*) mySigma%InnerDiamDParam
    read(cauxZ,*) mySigma%InnerDiamZParam
    mySigma%InnerDiamDParam = dSizeScale*mySigma%InnerDiamDParam
    mySigma%InnerDiamZParam = dSizeScale*mySigma%InnerDiamZParam

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelLength", mySigma%L ,myInf)
    mySigma%L = dSizeScale*mySigma%L

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelAxialStartPos", mySigma%L0 ,0d0)
    mySigma%L0 = dSizeScale*mySigma%L0
    
    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfElements", mySigma%NumberOfSeg ,0)
    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfFlights", mySigma%GANGZAHL ,0)
    

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","OnlyBarrelAdaptation", cOnlyBarrelAdaptation,"NO")
    call inip_toupper_replace(cOnlyBarrelAdaptation)
    IF (ADJUSTL(TRIM(cOnlyBarrelAdaptation)).EQ."YES") THEN
     mySigma%bOnlyBarrelAdaptation = .true.
    END IF

    IF (mySigma%NumberOfSeg.ge.1) THEN
     ALLOCATE (mySigma%mySegment(mySigma%NumberOfSeg))
    ELSE
     WRITE(*,*) "not a valid number of segments"
     WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
    ENDIF
    
    IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
     
     call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","RotationAxis", mySigma%RotationAxis ,"PARALLEL")
     call inip_toupper_replace(mySigma%RotationAxis)
     
     IF (.not.(ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."NONPARALLEL")) mySigma%RotationAxis = "PARALLEL"
     
     IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
      call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","CenterlineDistance", mySigma%a ,myInf)
      mySigma%a = dSizeScale*mySigma%a
     ELSE
      call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","RotAxisAngle", mySigma%RotAxisAngle ,myInf)
      mySigma%RotAxisAngle = 0.5d0*mySigma%RotAxisAngle*datan(1d0)/45d0
      call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","RotAxisCenter", mySigma%RotAxisCenter ,myInf)
      mySigma%RotAxisCenter = dSizeScale*mySigma%RotAxisCenter
     END IF
    END IF
    
    DO iSeg=1,mySigma%NumberOfSeg
     WRITE(cElement_i,'(A,I0)') 'E3DGeometryData/Machine/Element_',iSeg

     call INIP_getvalue_string(parameterlist,cElement_i,"ObjectType",mySigma%mySegment(iSeg)%ObjectType)
     call inip_toupper_replace(mySigma%mySegment(iSeg)%ObjectType)
     IF (.NOT.(mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW'.OR.mySigma%mySegment(iSeg)%ObjectType.eq.'MIXER'.OR.&
               mySigma%mySegment(iSeg)%ObjectType.eq.'DIE'.OR.mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE')) THEN
       WRITE(*,*) "STL object type is invalid. Only screw, mixer, die, or obstacle types are allowed"
     END IF

     call INIP_getvalue_string(parameterlist,cElement_i,"Unit",mySigma%mySegment(iSeg)%Unit,'MM')
     call inip_toupper_replace(mySigma%mySegment(iSeg)%Unit)
     IF (.NOT.(mySigma%mySegment(iSeg)%Unit.eq.'MM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'CM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'DM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'M')) THEN
       WRITE(*,*) "STL unit type is invalid. Only MM, CM, DM or 'M' units are allowed",mySigma%mySegment(iSeg)%Unit
       bReadError=.TRUE.
!       mySigma%mySegment(iSeg)%Unit = 'MM'
     END IF
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'MM') dElemSizeScale = 0.100d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'CM') dElemSizeScale = 1.000d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'DM') dElemSizeScale = 10.00d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'M')  dElemSizeScale = 100.0d0
      
!     WRITE(*,*) "'",TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)),"'", dElemSizeScale

     IF (mySigma%mySegment(iSeg)%ObjectType.eq.'MIXER') THEN
      call INIP_getvalue_double(parameterlist,cElement_i,"SegRotFreq", mySigma%mySegment(iSeg)%SegRotFreq ,myInf)
     END IF
     
     call INIP_getvalue_int(parameterlist,cElement_i,"NoOfFlights", mySigma%mySegment(iSeg)%GANGZAHL,-1)
     
     call INIP_getvalue_string(parameterlist,cElement_i,"Type",cElemType)
     mySigma%mySegment(iSeg)%ART = ' '
     call inip_toupper_replace(cElemType)

!!!==============================================     FOERD    =================================================================
!!!=============================================================================================================================
     IF (ADJUSTL(TRIM(cElemType)).eq."FOERD".or.ADJUSTL(TRIM(cElemType)).eq."THREADED") THEN
      mySigma%mySegment(iSeg)%ART   = "FOERD"
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Lead", mySigma%mySegment(iSeg)%t,myInf)
      mySigma%mySegment(iSeg)%t = dElemSizeScale*mySigma%mySegment(iSeg)%t*(myProcess%ind)


      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s

      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))

      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfConveying", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."CONVEYING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."CONVEYING") THEN
        mySigma%mySegment(iSeg)%t=1d0*mySigma%mySegment(iSeg)%t
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%t=-1d0*mySigma%mySegment(iSeg)%t
       END IF
      ELSE
       WRITE(*,*) "invalid kind of conveying segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
      END IF
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
     
!!!==============================================     KNET     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."KNET".or.ADJUSTL(TRIM(cElemType)).eq."KNEADING") THEN
      mySigma%mySegment(iSeg)%ART   = "KNET"
     
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
     
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out-mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))

      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dElemSizeScale*mySigma%mySegment(iSeg)%D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      mySigma%mySegment(iSeg)%alpha = mySigma%mySegment(iSeg)%alpha*(myProcess%ind)
      
      call INIP_getvalue_int(parameterlist,cElement_i,"KneadingDiscs", mySigma%mySegment(iSeg)%N,-1)
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfkneading", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."REKNEADING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING") THEN
        mySigma%mySegment(iSeg)%alpha=1d0*mySigma%mySegment(iSeg)%alpha
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%alpha=-1d0*mySigma%mySegment(iSeg)%alpha
       END IF
      ELSE
       WRITE(*,*) "invalid kind of kneading segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
      END IF
     
!!!==============================================    TKNET     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."TKNET".or.ADJUSTL(TRIM(cElemType)).eq."TRUEKNEADING") THEN
      mySigma%mySegment(iSeg)%ART   = "TKNET"
     
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
     
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out-mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))

      call INIP_getvalue_double(parameterlist,cElement_i,"DiscFrac", mySigma%mySegment(iSeg)%DiscFrac,0.05d0)
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dElemSizeScale*mySigma%mySegment(iSeg)%D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      mySigma%mySegment(iSeg)%alpha = mySigma%mySegment(iSeg)%alpha*(myProcess%ind)
      
      call INIP_getvalue_int(parameterlist,cElement_i,"KneadingDiscs", mySigma%mySegment(iSeg)%N,-1)
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfkneading", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."REKNEADING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING") THEN
        mySigma%mySegment(iSeg)%alpha=1d0*mySigma%mySegment(iSeg)%alpha
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%alpha=-1d0*mySigma%mySegment(iSeg)%alpha
       END IF
      ELSE
       WRITE(*,*) "invalid kind of kneading segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
      END IF
     
!!!==============================================     EKNET     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."EKNET".or.ADJUSTL(TRIM(cElemType)).eq."ECCENTRIC".or.ADJUSTL(TRIM(cElemType)).eq."EXCENTRIC") THEN
      mySigma%mySegment(iSeg)%ART   = "EKNET"
     
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds

      call INIP_getvalue_double(parameterlist,cElement_i,"Eccentricity", mySigma%mySegment(iSeg)%excentre, myInf)
      mySigma%mySegment(iSeg)%excentre = dElemSizeScale*mySigma%mySegment(iSeg)%excentre
     
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out-2d0*mySigma%mySegment(iSeg)%excentre-mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s - mySigma%mySegment(iSeg)%excentre))
    
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dElemSizeScale*mySigma%mySegment(iSeg)%D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      mySigma%mySegment(iSeg)%alpha = mySigma%mySegment(iSeg)%alpha*(myProcess%ind)
      
      call INIP_getvalue_int(parameterlist,cElement_i,"KneadingDiscs", mySigma%mySegment(iSeg)%N,-1)
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfkneading", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."REKNEADING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING") THEN
        mySigma%mySegment(iSeg)%alpha=1d0*mySigma%mySegment(iSeg)%alpha
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%alpha=-1d0*mySigma%mySegment(iSeg)%alpha
       END IF
      ELSE
       WRITE(*,*) "invalid kind of kneading segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
      END IF
     
!!!==============================================     SKNET    =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."SKNET".or.ADJUSTL(TRIM(cElemType)).eq."SHOULDEREDKNEADING") THEN
      mySigma%mySegment(iSeg)%ART   = "SKNET"
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds

      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dElemSizeScale*mySigma%mySegment(iSeg)%D

      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      mySigma%mySegment(iSeg)%alpha = mySigma%mySegment(iSeg)%alpha*(myProcess%ind)
      
      call INIP_getvalue_int(parameterlist,cElement_i,"KneadingDiscs", mySigma%mySegment(iSeg)%N,-1)
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%N*mySigma%mySegment(iSeg)%D
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfkneading", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."REKNEADING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."KNEADING") THEN
        mySigma%mySegment(iSeg)%alpha=1d0*mySigma%mySegment(iSeg)%alpha
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%alpha=-1d0*mySigma%mySegment(iSeg)%alpha
       END IF
      ELSE
       WRITE(*,*) "invalid kind of kneading segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
       bReadError=.TRUE.
!        GOTO 10
      END IF
     
!!!==============================================      SME     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."SME".or.ADJUSTL(TRIM(cElemType)).eq."SCREWMIXING") THEN
      mySigma%mySegment(iSeg)%ART   = "SME"
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Lead", mySigma%mySegment(iSeg)%t,myInf)
      mySigma%mySegment(iSeg)%t = dElemSizeScale*mySigma%mySegment(iSeg)%t*(myProcess%ind)
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))

      call INIP_getvalue_string(parameterlist,cElement_i,"KindOfConveying", cKindOfConveying," ")
      call inip_toupper_replace(cKindOfConveying)
      IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."CONVEYING".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq." ".OR.ADJUSTL(TRIM(cKindOfConveying)).eq."CONVEYING") THEN
        mySigma%mySegment(iSeg)%t=1d0*mySigma%mySegment(iSeg)%t
       END IF
       IF (ADJUSTL(TRIM(cKindOfConveying)).eq."RECONVEYING") THEN
        mySigma%mySegment(iSeg)%t=-1d0*mySigma%mySegment(iSeg)%t
       END IF
      ELSE
       WRITE(*,*) "invalid kind of screwmixing segment"
       WRITE(*,*) '"',ADJUSTL(TRIM(cKindOfConveying)),'"'
       bReadError=.TRUE.
!        GOTO 10
      END IF
      call INIP_getvalue_int(parameterlist,cElement_i,"NoOfGrooves",mySigma%mySegment(iSeg)%SecProf_N ,-1)
      
      call INIP_getvalue_int(parameterlist,cElement_i,"KindOfGrooves",mySigma%mySegment(iSeg)%SecProf_I  ,1)
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesDepth",mySigma%mySegment(iSeg)%SecProf_D  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_D = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesWidth",mySigma%mySegment(iSeg)%SecProf_W  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_W = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_W
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesLead",mySigma%mySegment(iSeg)%SecProf_L  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_L = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_L*(myProcess%ind)
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
     
!!!==============================================      ZME     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."ZME".or.ADJUSTL(TRIM(cElemType)).eq."TOOTHMIXING") THEN
      mySigma%mySegment(iSeg)%ART   = "ZME"
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_int(parameterlist,cElement_i,"NoOfRows",mySigma%mySegment(iSeg)%ZME_N ,-1)
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%ZME_DiscThick,myInf)
      mySigma%mySegment(iSeg)%ZME_DiscThick = dElemSizeScale*mySigma%mySegment(iSeg)%ZME_DiscThick
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscDiscGap", mySigma%mySegment(iSeg)%ZME_gap_SS,myInf)
      mySigma%mySegment(iSeg)%ZME_gap_SS = dElemSizeScale*mySigma%mySegment(iSeg)%ZME_gap_SS
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscShellGap", mySigma%mySegment(iSeg)%ZME_gap_SG,myInf)
      mySigma%mySegment(iSeg)%ZME_gap_SG = dElemSizeScale*mySigma%mySegment(iSeg)%ZME_gap_SG
!       call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = mySigma%Dz_Out - 2d0*mySigma%mySegment(iSeg)%ZME_gap_SG


      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,myInf)
      mySigma%mySegment(iSeg)%Dss = dElemSizeScale*mySigma%mySegment(iSeg)%Dss
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)

      call INIP_getvalue_int(parameterlist,cElement_i,"NoOfTeeth",mySigma%mySegment(iSeg)%SecProf_N ,-1)
      call INIP_getvalue_int(parameterlist,cElement_i,"KindOfGrooves",mySigma%mySegment(iSeg)%SecProf_I  ,1)
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesDepth",mySigma%mySegment(iSeg)%SecProf_D  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_D = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesWidth",mySigma%mySegment(iSeg)%SecProf_W  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_W = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_W
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesLead",mySigma%mySegment(iSeg)%SecProf_L  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_L = dElemSizeScale*mySigma%mySegment(iSeg)%SecProf_L*(myProcess%ind)

      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + &
      2d0*mySigma%mySegment(iSeg)%ZME_N * (mySigma%mySegment(iSeg)%ZME_DiscThick + mySigma%mySegment(iSeg)%ZME_gap_SS)
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
     
!!!==============================================      STL     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL".OR.ADJUSTL(TRIM(cElemType)).eq."OFF") THEN
      mySigma%mySegment(iSeg)%ART   = "STL"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
!       write(*,*) 'dElemSizeScale: ',dElemSizeScale
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dElemSizeScale)
      mySigma%mySegment(iSeg)%Dss = dElemSizeScale*mySigma%mySegment(iSeg)%Dss
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_string(parameterlist,cElement_i,"TranslationVelocity_cmps",cVelo,'unknown')
      call inip_toupper_replace(cVelo)
      if (adjustl(trim(cVelo)).ne.'UNKNOWN') then
       read(cVelo,*,err=166) mySigma%mySegment(iSeg)%FBMVeloBC
      else
       mySigma%mySegment(iSeg)%FBMVeloBC = [0d0,0d0,0d0]
      end if
      cVelo='unknown'
      GOTO 167     
166   write(*,*) 'WRONGLY DEFINED FBM velocity for the related segment!!'
      CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)
167   CONTINUE

      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      call INIP_getvalue_double(parameterlist,cElement_i,"OffsetAngle", mySigma%mySegment(iSeg)%OffsetAngle,myInf)
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0

      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
      
            
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")

      IF (mySigma%mySegment(iSeg)%nOFFfiles.ne.0) THEN
       
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
       
       do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
         call INIP_getvalue_string(parameterlist,cElement_i,"screwOFF",mySigma%mySegment(iSeg)%OFFfiles(iFile),isubstring=iFile)
       end do
       
      ELSE
       call INIP_getvalue_string(parameterlist,cElement_i,"OFF_FileList", cOFF_Files,' ')
       
       CALL ExtractNomOfCharFromString(cOFF_Files,mySigma%mySegment(iSeg)%nOFFfiles)
       
       IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
        ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
       ELSE
        WRITE(*,*) "STL geometry description files are missing"
        WRITE(*,*) 'screwOFF'
        bReadError=.TRUE.
        !GOTO 10
       END IF
       
       CALL ExtractCharArrayFromString(cOFF_Files,mySigma%mySegment(iSeg)%OFFfiles)
       
      END IF

      do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
        if (trim(adjustl(mySigma%mySegment(iSeg)%OFFfiles(iFile))) .eq. '') then
          write(*,'(A30,I3)') 'geometry description file ', iFile
          write(*,'(A15,I3,A20)') 'for element ', iSeg, ' is empty string!'
          bReadError = .TRUE.
        end if
      end do
      
      mySigma%mySegment(iSeg)%Max = mySigma%mySegment(iSeg)%L + mySigma%mySegment(iSeg)%Min
     
!!!==============================================      STL     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL_LR".OR.ADJUSTL(TRIM(cElemType)).eq."OFF_LR".or.&
      ADJUSTL(TRIM(cElemType)).eq."STL_RL".OR.ADJUSTL(TRIM(cElemType)).eq."OFF_RL") THEN
      mySigma%mySegment(iSeg)%ART   = "STL_LR"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
!       write(*,*) 'dElemSizeScale: ',dElemSizeScale
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dElemSizeScale)
      mySigma%mySegment(iSeg)%Dss = dElemSizeScale*mySigma%mySegment(iSeg)%Dss
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dElemSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dElemSizeScale*mySigma%mySegment(iSeg)%s
      
      call INIP_getvalue_double(parameterlist,cElement_i,"OffsetAngle", mySigma%mySegment(iSeg)%OffsetAngle,myInf)
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
            
      mySigma%mySegment(iSeg)%nOFFfilesL = INIP_querysubstrings(parameterlist,cElement_i,"screwOFFL")

      IF (mySigma%mySegment(iSeg)%nOFFfilesL.ne.0) THEN
       
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfilesL(mySigma%mySegment(iSeg)%nOFFfilesL))
       
       do iFile=1,mySigma%mySegment(iSeg)%nOFFfilesL
         call INIP_getvalue_string(parameterlist,cElement_i,"screwOFFL",mySigma%mySegment(iSeg)%OFFfilesL(iFile),isubstring=iFile)
       end do
       
      ELSE
       call INIP_getvalue_string(parameterlist,cElement_i,"OFF_FileListL", cOFF_Files,' ')
       
       CALL ExtractNomOfCharFromString(cOFF_Files,mySigma%mySegment(iSeg)%nOFFfilesL)
       
       IF (mySigma%mySegment(iSeg)%nOFFfilesL.gt.0) THEN
        ALLOCATE(mySigma%mySegment(iSeg)%OFFfilesL(mySigma%mySegment(iSeg)%nOFFfilesL))
       ELSE
        WRITE(*,*) "STL geometry description files are missing"
        WRITE(*,*) 'screwOFFL'
        bReadError=.TRUE.
        !GOTO 10
       END IF
       
       CALL ExtractCharArrayFromString(cOFF_Files,mySigma%mySegment(iSeg)%OFFfilesL)
       
      END IF

      do iFile=1,mySigma%mySegment(iSeg)%nOFFfilesL
        if (trim(adjustl(mySigma%mySegment(iSeg)%OFFfilesL(iFile))) .eq. '') then
          write(*,'(A30,I3)') 'geometry description file ', iFile
          write(*,'(A15,I3,A20)') 'for element ', iSeg, 'L is empty string!'
          bReadError = .TRUE.
        end if
      end do
      
      mySigma%mySegment(iSeg)%nOFFfilesR = INIP_querysubstrings(parameterlist,cElement_i,"screwOFFR")

      IF (mySigma%mySegment(iSeg)%nOFFfilesR.ne.0) THEN
       
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfilesR(mySigma%mySegment(iSeg)%nOFFfilesR))
       
       do iFile=1,mySigma%mySegment(iSeg)%nOFFfilesR
         call INIP_getvalue_string(parameterlist,cElement_i,"screwOFFR",mySigma%mySegment(iSeg)%OFFfilesR(iFile),isubstring=iFile)
       end do
       
      ELSE
       call INIP_getvalue_string(parameterlist,cElement_i,"OFF_FileListR", cOFF_Files,' ')
       
       CALL ExtractNomOfCharFromString(cOFF_Files,mySigma%mySegment(iSeg)%nOFFfilesR)
       
       IF (mySigma%mySegment(iSeg)%nOFFfilesR.gt.0) THEN
        ALLOCATE(mySigma%mySegment(iSeg)%OFFfilesR(mySigma%mySegment(iSeg)%nOFFfilesR))
       ELSE
        WRITE(*,*) "STL geometry description files are missing"
        WRITE(*,*) 'screwOFFR'
        bReadError=.TRUE.
        !GOTO 10
       END IF
       
       CALL ExtractCharArrayFromString(cOFF_Files,mySigma%mySegment(iSeg)%OFFfilesR)
       
      END IF

      do iFile=1,mySigma%mySegment(iSeg)%nOFFfilesR
        if (trim(adjustl(mySigma%mySegment(iSeg)%OFFfilesR(iFile))) .eq. '') then
          write(*,'(A30,I3)') 'geometry description file ', iFile
          write(*,'(A15,I3,A20)') 'for element ', iSeg, 'L is empty string!'
          bReadError = .TRUE.
        end if
      end do
      
      mySigma%mySegment(iSeg)%Max = mySigma%mySegment(iSeg)%L + mySigma%mySegment(iSeg)%Min
     
!!!==============================================      STL_R     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL_R") THEN
      mySigma%mySegment(iSeg)%ART   = "STL_R"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dElemSizeScale)
      mySigma%mySegment(iSeg)%Dss = dElemSizeScale*mySigma%mySegment(iSeg)%Dss
      
      call INIP_getvalue_double(parameterlist,cElement_i,"OffsetAngle", mySigma%mySegment(iSeg)%OffsetAngle,myInf)
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
      
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL_R geometry description files are missing"
       WRITE(*,*) 'screwOFF'
       bReadError=.TRUE.
       !GOTO 10
      END IF
      do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
        call INIP_getvalue_string(parameterlist,cElement_i,"screwOFF",mySigma%mySegment(iSeg)%OFFfiles(iFile),isubstring=iFile)
      end do
      mySigma%mySegment(iSeg)%Max = mySigma%mySegment(iSeg)%L + mySigma%mySegment(iSeg)%Min
     
!!!==============================================      STL_R     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL_L") THEN
      mySigma%mySegment(iSeg)%ART   = "STL_L"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
      mySigma%mySegment(iSeg)%Min = dElemSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dElemSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dElemSizeScale)
      mySigma%mySegment(iSeg)%Dss = dElemSizeScale*mySigma%mySegment(iSeg)%Dss
      
      call INIP_getvalue_double(parameterlist,cElement_i,"OffsetAngle", mySigma%mySegment(iSeg)%OffsetAngle,myInf)
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
      
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL_L geometry description files are missing"
       WRITE(*,*) 'screwOFF'
       bReadError=.TRUE.
       !GOTO 10
      END IF
      do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
        call INIP_getvalue_string(parameterlist,cElement_i,"screwOFF",mySigma%mySegment(iSeg)%OFFfiles(iFile),isubstring=iFile)
      end do
      mySigma%mySegment(iSeg)%Max = mySigma%mySegment(iSeg)%L + mySigma%mySegment(iSeg)%Min

     ELSE
       write(*,*) 'Unknown Elementtype ' // trim(adjustl(cElemType))
       bReadError = .TRUE.
     END IF

!!!=============================================================
!!! End reading Elementtyps
!!!=============================================================
     IF (mySigma%mySegment(iSeg)%ART.eq.' ') THEN
      WRITE(*,*) "not a valid ",iSeg, "-segment"
      WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
      bReadError=.TRUE.
      !GOTO 10
     ENDIF
    END DO

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","SegmentThermoPhysProps",cTXT,'OFF')
    call inip_toupper_replace(cTXT)
    IF (ADJUSTL(TRIM(cTXT)).eq."ON".or.ADJUSTL(TRIM(cTXT)).eq."YES") THEN
     myProcess%SegmentThermoPhysProps = .TRUE.
     allocate(myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)) 
     IF (myProcess%SegmentThermoPhysProps) THEN
      call INIP_getvalue_string(parameterlist,"E3DProcessParameters/SegmentThermoPhysProps","density",cTXT,' ')
      READ(cTXT,*) myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)%rho
      call INIP_getvalue_string(parameterlist,"E3DProcessParameters/SegmentThermoPhysProps","heatconductivity",cTXT,' ')
      READ(cTXT,*) myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)%lambda
      call INIP_getvalue_string(parameterlist,"E3DProcessParameters/SegmentThermoPhysProps","heatcapacity",cTXT,' ')
      READ(cTXT,*) myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)%cp
      call INIP_getvalue_string(parameterlist,"E3DProcessParameters/SegmentThermoPhysProps","isothermal",cTXT,' ')
      call inip_toupper_replace(cTXT)
      READ(cTXT,*) myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)%cConstTemp
      DO iSeg=1,mySigma%NumberOfSeg
       IF (TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."Y".or.&
           TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."YES".or.&
           TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."ON") THEN
         myProcess%SegThermoPhysProp(iSeg)%bConstTemp = .true.
       ELSE
         myProcess%SegThermoPhysProp(iSeg)%bConstTemp = .false.
       END IF
      END DO
      DO iSeg=1,mySigma%NumberOfSeg
       IF (TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."S".or.&
           TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."SRC".or.&
           TRIM(ADJUSTL(myProcess%SegThermoPhysProp(iSeg)%cConstTemp)).eq."SOURCE") THEN
         myProcess%SegThermoPhysProp(iSeg)%bHeatSource = .true.
       ELSE
         myProcess%SegThermoPhysProp(iSeg)%bHeatSource = .false.
       END IF
      END DO
      call INIP_getvalue_string(parameterlist,"E3DProcessParameters/SegmentThermoPhysProps","temperature",cTXT,' ')
      READ(cTXT,*) myProcess%SegThermoPhysProp(0:mySigma%NumberOfSeg)%T_Const
     END IF

    ELSE
     myProcess%SegmentThermoPhysProps = .FALSE.
    END IF

    myProcess%pTYPE = " "
    myProcess%dPress=myInf
    myProcess%Massestrom=myInf
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","ProcessType", cProcessType,'NoProcessType')
    call inip_toupper_replace(cProcessType)
    
    IF (ADJUSTL(TRIM(cProcessType)).eq."PRESSUREDROP") THEN
    
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","deltaP", cText,"_INVALID_")
    call inip_toupper_replace(cText)
    call ReadDoubleFromDimensionalString()
    IF (TRIM(adjustl(sExtract_Dim)).eq.'BAR') myProcess%dPress = 1d6*dExtract_Val
    IF (TRIM(adjustl(sExtract_Dim)).eq.'PA')  myProcess%dPress  = 1d1*dExtract_Val
    IF (TRIM(adjustl(sExtract_Dim)).eq.'KPA') myProcess%dPress = 1d4*dExtract_Val
    IF (TRIM(adjustl(sExtract_Dim)).eq.'MPA') myProcess%dPress = 1d7*dExtract_Val
!     
!      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","deltaP", myProcess%dPress,myInf)
     myProcess%pTYPE = "PRESSUREDROP"
    END IF
    IF (ADJUSTL(TRIM(cProcessType)).eq."THROUGHPUT") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","massthroughput", myProcess%Massestrom,myInf)
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MinInflowDiameter", myProcess%MinInflowDiameter,mySigma%Dz_In)
     myProcess%MinInflowDiameter = dSizeScale*myProcess%MinInflowDiameter
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaxInflowDiameter", myProcess%MaxInflowDiameter,mySigma%Dz_Out)
     myProcess%MaxInflowDiameter = dSizeScale*myProcess%MaxInflowDiameter
     myProcess%pTYPE = "THROUGHPUT"
    END IF
    IF (ADJUSTL(TRIM(myProcess%pTYPE)).eq." ") THEN
     WRITE(*,*) "no valid process type is defined"
     WRITE(*,*) '"',TRIM(cProcessType),'"'
     bReadError=.TRUE.
     !GOTO 10
    END IF
    
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","UseAirCooling", cText,"OFF")
    call inip_toupper_replace(cText)
    if (adjustl(trim(cText)).EQ.'ON'.or.adjustl(trim(cText)).EQ.'YES') then
     myProcess%UseAirCooling = .TRUE.
    else
     myProcess%UseAirCooling = .FALSE.
    end if
    
    IF (myProcess%UseAirCooling) THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","AirCoolingHeatTransCoeff", myProcess%AirCoolingHeatTransCoeff,20d0)
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","AirCoolingRoomTemperature", myProcess%AirCoolingRoomTemperature,25d0)
    END IF

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","UseHeatDissipationForQ1Scalar",cDissip,'OFF')
    call inip_toupper_replace(cDissip)
    if (adjustl(trim(cVelo)).EQ.'ON'.or.adjustl(trim(cVelo)).EQ.'YES') then
     myProcess%UseHeatDissipationForQ1Scalar = .TRUE.
    else
     myProcess%UseHeatDissipationForQ1Scalar = .FALSE.
    end if
    
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","FBMVeloBC",cVelo,'unknown')
    call inip_toupper_replace(cVelo)
    if (adjustl(trim(cVelo)).ne.'UNKNOWN') then
     read(cVelo,*,err=66) myProcess%FBMVeloBC
    else
     myProcess%FBMVeloBC = [0d0,0d0,0d0]
    end if
    GOTO 67     
66   write(*,*) 'WRONGLY DEFINED FBM velocity !!'
67  CONTINUE

    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewSpeed", myProcess%umdr,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","FillingDegree", myProcess%FillingDegree,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaterialTemperature",myProcess%T0,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaterialTemperatureSlope",myProcess%T0_slope,0d0)
    call INIP_getvalue_int(parameterlist,"E3DProcessParameters","MaterialTemperatureValuesN",myProcess%T0_N,0)
    IF (myProcess%T0_N.gt.0) THEN
     call INIP_getvalue_string(parameterlist,"E3DProcessParameters","MaterialTemperatureValuesT",cTVals,'')
     ALLOCATE(myProcess%T0_T(1:myProcess%T0_N))
     READ(cTVals,*) myProcess%T0_T(1:myProcess%T0_N)
    END IF

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","ScrewTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewTemperature",myProcess%Ti,myInf)
     call INIP_getvalue_string(parameterlist,"E3DProcessParameters","RobinBCScrew",cRobin,"NO")
     call inip_toupper_replace(cRobin)
     IF (ADJUSTL(TRIM(cRobin)).EQ."YES") THEN
      myProcess%bRobinBCScrew = .true.
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","RobinBCScrew_HTC",myProcess%RobinBCScrew_HTC,myInf)
     END IF
    ELSE
     myProcess%Ti=myInf
    END IF
    cTemperature=" "
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","BarrelTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","BarrelTemperature",myProcess%Ta,myInf)
     call INIP_getvalue_string(parameterlist,"E3DProcessParameters","RobinBCBarrel",cRobin,"NO")
     call inip_toupper_replace(cRobin)
     IF (ADJUSTL(TRIM(cRobin)).EQ."YES") THEN
      myProcess%bRobinBCBarrel = .true.
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","RobinBCBarrel_HTC",myProcess%RobinBCBarrel_HTC,myInf)
     END IF
    ELSE
     IF (ADJUSTL(TRIM(cTemperature)).EQ."FLUX".or.ADJUSTL(TRIM(cTemperature)).EQ."YES") THEN
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","HeatFluxThroughBarrelWall_kWm2",myProcess%HeatFluxThroughBarrelWall_kWm2,0d0)
     END IF
     myProcess%Ta=myInf
    END IF
       
    call INIP_getvalue_int(parameterlist,"E3DProcessParameters",   "nOfTempBCs"      ,myProcess%nOfTempBCs,0)
    ALLOCATE(myProcess%myTempBC(myProcess%nOfTempBCs))
    
    cParserString = "E3DProcessParameters"
    CALL FillUpTempBCs(myProcess%nOfTempBCs,cParserString)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material and Material-specific read section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    call INIP_getvalue_int(parameterlist,"E3DProcessParameters",   "nOfInflows"      ,myProcess%nOfInflows,0)
    ALLOCATE(myProcess%myInflow(myProcess%nOfInflows))
    cParserString = "E3DProcessParameters"
    CALL FillUpInflows(myProcess%nOfInflows,cParserString)

    if (myProcess%FillingDegree.eq.myInf .or. myProcess%FillingDegree .eq. 1d0) then
      myProcess%FillingDegree=myInf
      ! this applies for the use of the Q2P1 sse code !IF! there IS NO partial filling defined
      call INIP_getvalue_Int(parameterlist,"E3DMaterialParameters","NoOfMaterials", myMultiMat%nOfMaterials,0)
      IF (myMultiMat%nOfMaterials.ne.0) then

       call INIP_getvalue_Int(parameterlist,"E3DMaterialParameters","InitMaterial", myMultiMat%InitMaterial,1)
       if (myMultiMat%nOfMaterials.ge.2) bMultiMat = .true.
       if (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") bMultiMat = .true.

       ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
       
       DO iMat = 1, myMultiMat%nOfMaterials
       
       WRITE(cParserString,'(A,I0)') "E3DMaterialParameters/Mat_",iMat
  !      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cParserString)),"LimitViscoMin",myMultiMat%Mat(iMat)%Rheology%ViscoMin,1d0)
  !      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cParserString)),"LimitViscoMax",myMultiMat%Mat(iMat)%Rheology%ViscoMax,1d5)

        WRITE(cParserString,'(A,I0,A)') "E3DMaterialParameters/Mat_",iMat,"/RheologicalData"
        CALL FillUpRheoData(myMultiMat%Mat(iMat)%Rheology,cParserString)

        WRITE(cParserString,'(A,I0,A)') "E3DMaterialParameters/Mat_",iMat,"/ThermoData"
        CALL FillUpThermoData(myMultiMat%Mat(iMat)%Thermodyn,cParserString)
        
        myThermodyn = myMultiMat%Mat(myMultiMat%InitMaterial)%Thermodyn

       END DO
      
      ELSE
       myMultiMat%nOfMaterials = 1
       ALLOCATE(myMultiMat%Mat(1))
       myMultiMat%InitMaterial = 1
       
       cParserString = "E3DProcessParameters/Material/ThermoData"
       CALL FillUpThermoData(myThermodyn,cParserString)

       cParserString = "E3DProcessParameters/Material/RheologicalData"
       CALL FillUpRheoData(myMultiMat%Mat(1)%Rheology,cParserString)
       
       myMultiMat%Mat(1)%Thermodyn = myThermodyn
      END IF
      
    else
      ! this applies for the use of the Q2P1 sse code !IF! there IS partial filling defined
      myMultiMat%nOfMaterials = 2
      myMultiMat%InitMaterial = 1
      bMultiMat = .true.
      ALLOCATE(myMultiMat%Mat(myMultiMat%nOfMaterials))
      
      cParserString = "E3DProcessParameters/Material/ThermoData"
      CALL FillUpThermoData(myThermodyn,cParserString)

      cParserString = "E3DProcessParameters/Material/RheologicalData"
      CALL FillUpRheoData(myMultiMat%Mat(1)%Rheology,cParserString)
      
      myMultiMat%Mat(1)%Thermodyn = myThermodyn
      
      CALL FillUpRheoAndThermoDataForPseudoAir(myMultiMat%Mat(2)%Rheology,myMultiMat%Mat(2)%Thermodyn)
    
    end if
    
!     pause
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Material and Material-specific read section !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output",   "nOf1DLayers"      ,myOutput%nOf1DLayers,16)
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output",   "nOfHistogramBins" ,myOutput%nOfHistogramBins,16)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramShearMax",myOutput%HistogramShearMax,1d5)
!     call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramShearMin",myOutput%HistogramShearMin,1d-2)
!     call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramViscoMax",myOutput%HistogramViscoMax,1d6)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramViscoMin",myOutput%HistogramViscoMin,1d0)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","CutDtata_1D",myOutput%CutDtata_1D,0.001d0)
    
    cKTP=' '
    IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","RotationalFramOfReference",cALE,"NO")
     call inip_toupper_replace(cALE)
     IF (ADJUSTL(TRIM(cALE)).eq."YES".OR.ADJUSTL(TRIM(cALE)).eq."ON") THEN
      mySetup%bRotationalFramOfReference = .TRUE.
     ELSE
      mySetup%bRotationalFramOfReference = .FALSE.
     END IF
    END IF
    
    IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","AutomaticTimeStepControl",cKTP,"YES")
     call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","TimeStepEnlargmentFactor",dTimeStepEnlargmentFactor,5d0)
    ELSE
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","AutomaticTimeStepControl",cKTP,"NO")
     call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","TimeStepEnlargmentFactor",dTimeStepEnlargmentFactor,1d0)
    END IF
    
    call inip_toupper_replace(cKTP)
    IF (ADJUSTL(TRIM(cKTP)).eq."NO") THEN
     mySetup%bAutomaticTimeStepControl = .FALSE.
    END IF
    IF (ADJUSTL(TRIM(cKTP)).eq."YES") THEN
     mySetup%bAutomaticTimeStepControl = .TRUE.
    END IF
    
    IF (mySetup%bAutomaticTimeStepControl.and.ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") then
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ExtrusionSpeed_CMpS",myProcess%ExtrusionSpeed,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ExtrusionOutflowArea_MM2",dArea,myInf)
      
      
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ExtrusionGapSize_MM",myProcess%ExtrusionGapSize,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ExtrusionGapFactor",myProcess%ExtrusionGapFactor,1d0)
      if (myProcess%ExtrusionSpeed.eq.myInf) then
       if (dArea.eq.myInf) then
        if (myid.eq.1) WRITE(*,*) "   Extrusion Speed is not set 'E3DProcessParameters@ExtrusionSpeed_CMpS'"
        if (myid.eq.1) WRITE(*,*) "or Extrusion Area is not set 'E3DProcessParameters@ExtrusionOutflowArea_MM2'"
        CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

       else
        dFlow = 0d0
        do iInflow=1,myProcess%nOfInflows
         IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
          dFlow = dFlow + myProcess%myInflow(iInflow)%massflowrate  
         ELSE
          DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
           dFlow = dFlow + myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate  
          END DO
         END IF
        enddo
        myProcess%ExtrusionSpeed = 1d2*(dFlow/3600d0)/(dArea*1e-6)/(1d3*(myThermodyn%densityT0 - myProcess%T0 * myThermodyn%densitySteig))
       end if
      end if
      if (myProcess%ExtrusionGapSize.eq.myInf) then
       if (myid.eq.1) WRITE(*,*) "Extrusion GapSize is not set 'E3DProcessParameters@ExtrusionGapSize_MM'"
        CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

      end if
    End if
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","CharacteristicShearRate",mySetup%CharacteristicShearRate,1d0)
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","activeFBM_Z_Position",activeFBM_Z_Position,myInf)

    cKTP=' '
    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","KTPRelease",cKTP,"YES")
    call inip_toupper_replace(cKTP)
    IF (ADJUSTL(TRIM(cKTP)).eq."NO") THEN
     bKTPRelease = .FALSE.
    END IF

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","PressureFBM", cPressureFBM,"OFF")
    call inip_toupper_replace(cPressureFBM)
    IF (ADJUSTL(TRIM(cPressureFBM)).eq."ON".OR.ADJUSTL(TRIM(cPressureFBM)).eq."YES") THEN
     mySetup%bPressureFBM = .true.
    ENDIF
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","PressureConvergenceTolerance", mySetup%PressureConvergenceTolerance,5d-3)

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","HexMesher", mySetup%cMesher,"OFF")
    call inip_toupper_replace(mySetup%cMesher)

    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."OFF") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","MeshPath", mySetup%cMeshPath,'.')
    END IF
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."HOLLOWCYLINDER") THEN
    
     mySetup%MeshResolution = 0
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","MeshQuality",cMeshQuality,'NO')
     call inip_toupper_replace(cMeshQuality)
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."EXTREME".or.ADJUSTL(TRIM(cMeshQuality)).eq."EXTREM") THEN
      mySetup%MeshResolution = 5
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."SUPER") THEN
      mySetup%MeshResolution = 4
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."FINE".or.ADJUSTL(TRIM(cMeshQuality)).eq."FEIN") THEN
      mySetup%MeshResolution = 3
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."MEDIUM".or.ADJUSTL(TRIM(cMeshQuality)).eq."MITTEL") THEN
      mySetup%MeshResolution = 2
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."ROUGH".or.ADJUSTL(TRIM(cMeshQuality)).eq."GROB".or.ADJUSTL(TRIM(cMeshQuality)).eq."COARSE") THEN
      mySetup%MeshResolution = 1
     END IF
     
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Tangential",mySetup%m_nT,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Radial",mySetup%m_nR,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Axial",mySetup%m_nZ,0)
     IF ((mySetup%m_nT.eq.0.or.mySetup%m_nR.eq.0.or.mySetup%m_nZ.eq.0).and.mySetup%MeshResolution.eq.0) THEN
      WRITE(*,*) "mesh resolution is not correctly defined: nT,nR,nZ & nRes"
      WRITE(*,*) '"',mySetup%m_nT,mySetup%m_nR,mySetup%m_nZ,mySetup%MeshResolution,'"'
      bReadError=.TRUE.
      !GOTO 10
     END IF
    END IF

    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."FULLCYLINDER") THEN
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Repetition",mySetup%m_nP,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Core",mySetup%m_nT,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Radial",mySetup%m_nR,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Axial",mySetup%m_nZ,0)
     IF (mySetup%m_nT.eq.0.or.mySetup%m_nZ.eq.0.or.mySetup%m_nP.eq.0) THEN
      WRITE(*,*) "mesh resolution is not correctly defined"
      WRITE(*,*) '"',mySetup%m_nT,mySetup%m_nZ,mySetup%m_nP,'"'
      bReadError=.TRUE.
!       GOTO 10
     END IF
    END IF
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."TWINSCREW") THEN
     mySetup%MeshResolution = 0
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","MeshQuality",cMeshQuality,'GROB')
     call inip_toupper_replace(cMeshQuality)
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."EXTREME".or.ADJUSTL(TRIM(cMeshQuality)).eq."EXTREM") THEN
      mySetup%MeshResolution = 5
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."SUPER") THEN
      mySetup%MeshResolution = 4
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."FINE".or.ADJUSTL(TRIM(cMeshQuality)).eq."FEIN") THEN
      mySetup%MeshResolution = 3
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."MEDIUM".or.ADJUSTL(TRIM(cMeshQuality)).eq."MITTEL") THEN
      mySetup%MeshResolution = 2
     END IF
     IF (ADJUSTL(TRIM(cMeshQuality)).eq."ROUGH".or.ADJUSTL(TRIM(cMeshQuality)).eq."GROB".or.ADJUSTL(TRIM(cMeshQuality)).eq."COARSE") THEN
      mySetup%MeshResolution = 1
     END IF
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Tangential1",mySetup%m_nT1,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Tangential2",mySetup%m_nT2,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Radial",mySetup%m_nR,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nEl_Axial",mySetup%m_nZ,0)
     IF (mySetup%MeshResolution.eq.0) THEN
      WRITE(*,*) "mesh quality/resolution is not defined"
      WRITE(*,*) '"',TRIM(cMeshQuality),'"'
      bReadError=.TRUE.
  !      GOTO 10
     END IF
    END IF
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."BOX") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherELems",cMeshQuality,'0,0,0')
     read(cMeshQuality,*) mySetup%m_nX,mySetup%m_nY,mySetup%m_nZ
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","BoxMesherNumberOfELems",mySetup%nBoxElem,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","MeshResolution",mySetup%MeshResolution,0)
     
     IF ((mySetup%nBoxElem.le.0).and.&
         (mySetup%m_nX.le.0.and.mySetup%m_nY.le.0.and.mySetup%m_nZ.le.0).and.&
          mySetup%MeshResolution.le.0) THEN
          if (myid.eq.1) WRITE(*,*) 'No rules defined to create a mesh...'
         CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

     END IF
     
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherUnit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0
    
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherX",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(1,:)
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherY",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(2,:)
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherZ",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(3,:)
     
     !!!!!!!!!! Scaling the geo to its right dimension !!!!
     mySetup%m_box = daux*mySetup%m_box
    
    END IF
!     cMeshQuality=' '
!     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","SendEmail",cMeshQuality,"YES")
!     call inip_toupper_replace(cMeshQuality)
!     IF (ADJUSTL(TRIM(cMeshQuality)).eq."NO") THEN
!      mySetup%bSendEmail = .FALSE.
!     END IF
!     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nSolutions",mySetup%nSolutions,1)

    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","Periodicity",myProcess%Periodicity,1)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","dAlpha",myProcess%dAlpha,10d0)
    
    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","AnalyticalShearRateRestriction",cShearRateRest,'NO')
    call inip_toupper_replace(cShearRateRest)
    IF (adjustl(trim(cShearRateRest)).eq."YES") then
     mySigma%bAnalyticalShearRateRestriction = .true.
    end if

    myProcess%nTimeLevels = INT(360d0/myProcess%dAlpha)
    if (myProcess%nTimeLevels*INT(myProcess%dAlpha)-360.ne.0) THEN
     if (myid.eq.1) WRITE(*,*) 'Wrong [E3DSimulationSettings/dAlpha] value defined (results in non-integer time-level number!) '
    END IF
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","Angle",myProcess%Angle,myInf)
    if (SSE_HAS_ANGLE)then
      myProcess%Angle = extruder_angle
    else 
      extruder_angle = myProcess%Angle
    end if
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","Phase",myProcess%Phase,-1)
    
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Sensors",   "nOfDIESensors"      ,myProcess%nOfDIESensors,0)
    IF (myProcess%nOfDIESensors.gt.0) then
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings/Sensors","SensorUnit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0

     ALLOCATE(myProcess%mySensor(myProcess%nOfDIESensors)) 
     DO i=1,myProcess%nOfDIESensors
      cSensor = ' '
      if (i.gt.0.and.i.le.9) WRITE(cSensor,'(A,I1.1)') 'Sensor_',i
      if (i.gt.9.and.i.le.99) WRITE(cSensor,'(A,I2.2)') 'Sensor_',i
      if (i.gt.99.and.i.le.999) WRITE(cSensor,'(A,I3.3)') 'Sensor_',i
      
      call INIP_getvalue_string(parameterlist,"E3DSimulationSettings/Sensors",cSensor,cCenter,'0,0,0,0,OFF,0.0,0')
      read(cCenter,*,err=153,end=153) myProcess%mySensor(i)%Center,myProcess%mySensor(i)%Radius,myProcess%mySensor(i)%Type,myProcess%mySensor(i)%PID_Ctrl%T_set, myProcess%mySensor(i)%iSeg
      CALL inip_toupper_replace(myProcess%mySensor(i)%Type)
      myProcess%mySensor(i)%PID_Ctrl%Omega_P = 0.010d0
      myProcess%mySensor(i)%PID_Ctrl%Omega_I = 1d2
      myProcess%mySensor(i)%PID_Ctrl%Omega_D = 0.1d0
      mySigma%mySegment(i)%PID_Ctrl%SumI = 0d0
      mySigma%mySegment(i)%PID_Ctrl%e_old = 0d0

      IF (.not.(myProcess%mySensor(i)%Type.eq."PID".OR.myProcess%mySensor(i)%Type.eq."OFF".or.&
                myProcess%mySensor(i)%Type.eq."-STOP".OR.myProcess%mySensor(i)%Type.eq."+STOP")) THEN
       write(*,*) 'WRONGLY DEFINED Regulation type for Sensor',i,' !!',adjustl(trim(myProcess%mySensor(i)%Type))
       goto 153
      END IF
      IF (myProcess%mySensor(i)%iSeg.lt.1.OR.myProcess%mySensor(i)%iSeg.gt.mySigma%NumberOfSeg) THEN
       write(*,*) 'WRONGLY DEFINED Regulation pointer, no valid segment ID',i,' !!',myProcess%mySensor(i)%iSeg
       goto 153
      end if
      goto 156

153   CONTINUE
      myProcess%mySensor(i)%Type = "OFF"
      myProcess%mySensor(i)%PID_Ctrl%T_set = myInf
      read(cCenter,*,err=155) myProcess%mySensor(i)%Center,myProcess%mySensor(i)%Radius
      goto 156
155   write(*,*) 'WRONGLY DEFINED Center and/or Radius for Sensor',i,' !!'//"|",ADJUSTL(TRIM(cCenter)),"|"
156   continue 
      myProcess%mySensor(i)%Center = daux*myProcess%mySensor(i)%Center
      myProcess%mySensor(i)%Radius = daux*myProcess%mySensor(i)%Radius
     END DO
     
    END IF
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2021.09.30 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO iMat = 1, myMultiMat%nOfMaterials
     IF (ADJUSTL(TRIM(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)).eq."DENSITY") THEN
      myMultiMat%Mat(iMat)%Thermodyn%density = myMultiMat%Mat(iMat)%Thermodyn%densityT0 - myProcess%T0 * myMultiMat%Mat(iMat)%Thermodyn%densitySteig
     END IF
     IF (ADJUSTL(TRIM(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)).eq."SPECVOLUME") THEN
      myMultiMat%Mat(iMat)%Thermodyn%density = 1d0/(myMultiMat%Mat(iMat)%Thermodyn%densityT0 + myProcess%T0 * myMultiMat%Mat(iMat)%Thermodyn%densitySteig)
     END IF
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2021.09.30 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2020.11.20 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF (ADJUSTL(TRIM(myThermodyn%DensityModel)).eq."DENSITY") THEN
     myThermodyn%density = myThermodyn%densityT0 - myProcess%T0 * myThermodyn%densitySteig
    END IF
    IF (ADJUSTL(TRIM(myThermodyn%DensityModel)).eq."SPECVOLUME") THEN
     myThermodyn%density = 1d0/(myThermodyn%densityT0 + myProcess%T0 * myThermodyn%densitySteig)
    END IF
    
    myThermodyn%lambda = myThermodyn%lambda + myProcess%T0 * myThermodyn%lambdaSteig
    myThermodyn%cp = myThermodyn%cp + myProcess%T0 * myThermodyn%cpSteig
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2020.11.20 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    DO iSeg=1,mySigma%NumberOfSeg
     IF (mySigma%mySegment(iSeg)%GANGZAHL.eq.-1) mySigma%mySegment(iSeg)%GANGZAHL=mySigma%GANGZAHL
    END DO
    
    IF (myid.eq.1.or.subnodes.eq.0) then
    write(*,*) "=========================================================================="
    write(*,*) "mySigma%Type",'=',trim(mySigma%cType)
    write(*,*) "mySigma%ScrewCylinderRendering",'=',mySigma%ScrewCylinderRendering
    write(*,*) "mySigma%BarrelRendering",'=',mySigma%BarrelRendering
    write(*,*) "mySigma%Zwickel",'=',trim(mySigma%cZwickel)
    write(*,*) "mySigma%InnerDiamNParam",'=',mySigma%InnerDiamNParam
    write(*,'(A,A,100ES12.4)') "mySigma%InnerDiamDParam",'=',mySigma%InnerDiamDParam
    write(*,'(A,A,100ES12.4)') "mySigma%InnerDiamZParam",'=',mySigma%InnerDiamZParam
    
    IF (mySigma%Dzz.ne.myinf) write(*,*) "mySigma%Dzz",'=',mySigma%Dzz
    IF (mySigma%W.ne.myinf)   write(*,*) "mySigma%W",'=',mySigma%W
    if (mySigma%bOnlyBarrelAdaptation) then
     write(*,*) "mySigma%OnlyBarrelAdaptation",'= YES'
    else
     write(*,*) "mySigma%OnlyBarrelAdaptation",'= NO'
    end if
    
    write(*,*) "mySigma%ShearRateRest",'=',mySigma%bAnalyticalShearRateRestriction
    
    IF (adjustl(trim(mySigma%GeometryStart)).ne.'_INVALID_') THEN
     write(*,'(A,A,3f13.3)') "mySigma%GeometryStart",'=',mySigma%DIE_Start
    END IF
    IF (adjustl(trim(mySigma%GeometryLength)).ne.'_INVALID_') THEN
     write(*,'(A,A,3f13.3)') "mySigma%GeometryLength",'=',mySigma%DIE_Length
    END IF
    IF (adjustl(trim(mySigma%GeometrySymmetryBC)).ne.'_INVALID_') THEN
     write(*,'(A,A,6L2)') "mySigma%GeometrySymmetryBC",'=',mySigma%DIE_SymmetryBC
    END IF
    
    write(*,*) "mySigma%Dz_Out",'=',mySigma%Dz_out
    write(*,*) "mySigma%Dz_In",'=',mySigma%Dz_In
    write(*,*) "mySigma%L",'=',mySigma%L
    write(*,*) "mySigma%BarrelAxialStartPos",'=',mySigma%L0
    write(*,*) "mySigma%GANGZAHL",'=',mySigma%GANGZAHL
    if (myProcess%iInd.eq.-1) write(*,*) "mySigma%RotationType",'=','CounterRotating'
    if (myProcess%iInd.eq.+1) write(*,*) "mySigma%RotationType",'=','CoRotating'
    IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN
     
     write(*,*) "mySigma%RotationAxis",'=',mySigma%RotationAxis
     
     IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
      write(*,*) "mySigma%a",'=',mySigma%a
     ELSE
      write(*,*) "mySigma%RotAxisAngle",'=',mySigma%RotAxisAngle
      write(*,*) "mySigma%RotAxisCenter",'=',mySigma%RotAxisCenter
     END IF
    END IF
    
    
    write(*,*) 
    IF (myProcess%SegmentThermoPhysProps) THEN
     write(*,*) "Thermal properties has been assigned to the segments!"
     DO iSeg=0,mySigma%NumberOfSeg
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%density=',myProcess%SegThermoPhysProp(iSeg)%rho
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%lambda=' ,myProcess%SegThermoPhysProp(iSeg)%lambda
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%cp='     ,myProcess%SegThermoPhysProp(iSeg)%cp
      write(*,'(A,I0,A,L1)') " mySIGMA%Segment(",iSeg,')%isothermal='     ,myProcess%SegThermoPhysProp(iSeg)%bConstTemp
      write(*,'(A,I0,A,L1)') " mySIGMA%Segment(",iSeg,')%heatsource='     ,myProcess%SegThermoPhysProp(iSeg)%bHeatSource
      if (myProcess%SegThermoPhysProp(iSeg)%bConstTemp) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%temperature='     ,myProcess%SegThermoPhysProp(iSeg)%T_Const
      end if
      if (myProcess%SegThermoPhysProp(iSeg)%bHeatSource) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%HeatSource='     ,myProcess%SegThermoPhysProp(iSeg)%T_Const
      end if
     END DO  
    ELSE
     write(*,*) "No thermal properties has been assigned to the segments!"
    END IF
    

    write(*,*) 
    write(*,*) "mySigma%NumberOfSeg",'=',mySigma%NumberOfSeg
    write(*,*) 
    DO iSeg=1,mySigma%NumberOfSeg
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Art=',mySigma%mySegment(iSeg)%ART
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%ObjectType=',mySigma%mySegment(iSeg)%ObjectType
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Unit=',mySigma%mySegment(iSeg)%Unit
     write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Min=',mySigma%mySegment(iSeg)%Min
     write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Max=',mySigma%mySegment(iSeg)%Max
     write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%L=',mySigma%mySegment(iSeg)%L
     write(*,'(A,I0,A,I5)') " mySIGMA%Segment(",iSeg,')%nFl=',mySigma%mySegment(iSeg)%GANGZAHL
     IF (mySigma%mySegment(iSeg)%ObjectType.eq.'MIXER') THEN
      IF (ieee_is_finite(mySigma%mySegment(iSeg)%OffsetAngle)) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SegRotFreq=',mySigma%mySegment(iSeg)%SegRotFreq
      else
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SegRotFreq=',myProcess%Umdr
      end if
     END IF
     
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."FOERD") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%t=',abs(mySigma%mySegment(iSeg)%t)
      IF (mySigma%mySegment(iSeg)%t.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','RECONVEYING'
      END IF
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."KNET") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."TKNET") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%DiscFrac=',mySigma%mySegment(iSeg)%DiscFrac
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."EKNET") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%E=',mySigma%mySegment(iSeg)%excentre
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."SKNET") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta

      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."SME") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%t=',abs(mySigma%mySegment(iSeg)%t)
      IF (mySigma%mySegment(iSeg)%t.lt.0d0) THEN
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','RECONVEYING'
      END IF
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta

      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_N=',mySigma%mySegment(iSeg)%SecProf_N
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_I=',mySigma%mySegment(iSeg)%SecProf_I
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_D=',mySigma%mySegment(iSeg)%SecProf_d
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_W=',mySigma%mySegment(iSeg)%SecProf_w
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_L=',mySigma%mySegment(iSeg)%SecProf_l
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."ZME") THEN
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%ZME_N=',mySigma%mySegment(iSeg)%ZME_N
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_DiscThick=',mySigma%mySegment(iSeg)%ZME_DiscThick
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_gap_SG=',mySigma%mySegment(iSeg)%ZME_gap_SG
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_gap_SS=',mySigma%mySegment(iSeg)%ZME_gap_SS
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_N=',mySigma%mySegment(iSeg)%SecProf_N
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_I=',mySigma%mySegment(iSeg)%SecProf_I
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_D=',mySigma%mySegment(iSeg)%SecProf_d
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_W=',mySigma%mySegment(iSeg)%SecProf_w
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_L=',mySigma%mySegment(iSeg)%SecProf_l
     END IF
     
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
      IF (ieee_is_finite(mySigma%mySegment(iSeg)%OffsetAngle)) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%OffsetAngle=',mySigma%mySegment(iSeg)%OffsetAngle
      END IF
      
      write(*,'(A,I0,A,3f13.3)') " mySIGMA%Segment(",iSeg,')%FBMVeloBC=', mySigma%mySegment(iSeg)%FBMVeloBC
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_LR") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
      IF (ieee_is_finite(mySigma%mySegment(iSeg)%OffsetAngle)) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%OffsetAngle=',mySigma%mySegment(iSeg)%OffsetAngle
      END IF
      
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfilesL=",mySigma%mySegment(iSeg)%nOFFfilesL
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesL
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfilesL(iFile))),'"'
      END DO
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfilesR=",mySigma%mySegment(iSeg)%nOFFfilesR
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfilesR
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfilesR(iFile))),'"'
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_R") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      IF (ieee_is_finite(mySigma%mySegment(iSeg)%OffsetAngle)) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%OffsetAngle=',mySigma%mySegment(iSeg)%OffsetAngle
      END IF
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_L") THEN
      write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      IF (ieee_is_finite(mySigma%mySegment(iSeg)%OffsetAngle)) then
       write(*,'(A,I0,A,f13.3)') " mySIGMA%Segment(",iSeg,')%OffsetAngle=',mySigma%mySegment(iSeg)%OffsetAngle
      END IF
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     write(*,*) 
    END DO    

    write(*,'(A,A,L2)') "myProcess%UseHeatDissipationForQ1Scalar",'=',myProcess%UseHeatDissipationForQ1Scalar
    
    IF (myProcess%UseAirCooling) THEN
     write(*,'(A)') "AirCooling of the outer surface is 'ON'!"
     write(*,'(A,A,3ES12.4)') "myProcess%AirCoolingHeatTransCoeff"," = ", myProcess%AirCoolingHeatTransCoeff
     write(*,'(A,A,3ES12.4)') "myProcess%AirCoolingRoomTemperature"," = ", myProcess%AirCoolingRoomTemperature
     write(*,'(A)') " "
    END IF
    
    write(*,'(A,A,3ES12.4)') "myProcess%FBMVeloBC",'=',myProcess%FBMVeloBC
    write(*,*) "myProcess%Rotation",'=',myProcess%Rotation
    write(*,*) "myProcess%ind",'=',myProcess%ind
    write(*,*) "myProcess%deltaP",'=',myProcess%dPress
    write(*,*) "myProcess%Massestrom",'=',myProcess%Massestrom
    write(*,*) "myProcess%MinInflowDiameter",'=',myProcess%MinInflowDiameter
    write(*,*) "myProcess%MaxInflowDiameter",'=',myProcess%MaxInflowDiameter
    write(*,*) "myProcess%f",'=',myProcess%umdr
    IF (ieee_is_finite(myProcess%FillingDegree)) then
     write(*,*) "myProcess%FillingDegree",'=',myProcess%FillingDegree
    END IF

    write(*,*) "myProcess%Ti",'=',myProcess%Ti
    IF (myProcess%bRobinBCScrew) then
     write(*,*) "myProcess%SCREW_RobinHTC",'=',myProcess%RobinBCScrew_HTC
    end if
    write(*,*) "myProcess%Ta",'=',myProcess%Ta
    IF (myProcess%bRobinBCBarrel) then
     write(*,*) "myProcess%BARREL_RobinHTC",'=',myProcess%RobinBCBarrel_HTC
    END IF
    IF (myProcess%Ta.eq.myInf) THEN
     write(*,*) "myProcess%FLUX",'=',myProcess%HeatFluxThroughBarrelWall_kWm2
    END IF
    write(*,*) "myProcess%T0",'=',myProcess%T0
    IF (myProcess%T0_N.gt.0) THEN
     write(*,'(A,A,100ES12.4)')  "myProcess%T0_Zones",'=',myProcess%T0_T
    else
     write(*,*) "myProcess%T0_Slope",'=',myProcess%T0_Slope
    END IF

    write(*,*) 
    write(*,*) "myProcess%nOfTempBCs",'=',myProcess%nOfTempBCs
    DO iTempBC=1,myProcess%nOfTempBCs
     if (iTempBC.gt.0.and.iTempBC.le.9) WRITE(cTempBC_i,'(A,I1.1)') 'TempBC_',iTempBC
     if (iTempBC.gt.9.and.iTempBC.le.99) WRITE(cTempBC_i,'(A,I2.2)') 'TempBC_',iTempBC
     if (iTempBC.gt.99.and.iTempBC.le.999) WRITE(cTempBC_i,'(A,I3.3)') 'TempBC_',iTempBC
    
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cTempBC_i))//'_Temperature','=',myProcess%myTempBC(iTempBC)%temperature
     write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cTempBC_i))//'_center'//'=',myProcess%myTempBC(iTempBC)%center
     write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cTempBC_i))//'_gradient'//'=',myProcess%myTempBC(iTempBC)%gradient
    END DO
    
    write(*,*) 
    write(*,*) "myProcess%nOfInflows",'=',myProcess%nOfInflows
    DO iInflow=1,myProcess%nOfInflows
     if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'Inflow_',iInflow
     if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'Inflow_',iInflow
     if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'Inflow_',iInflow
     
     IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Type','=',myProcess%myInflow(iInflow)%iBCtype
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Material','=',myProcess%myInflow(iInflow)%Material
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Massflowrate','=',myProcess%myInflow(iInflow)%massflowrate
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_TemperatureType','=',myProcess%myInflow(iInflow)%temperatureType
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_TemperatureRange','=',myProcess%myInflow(iInflow)%temperatureRange
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Temperature','=',myProcess%myInflow(iInflow)%temperature
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_InnerRadius','=',myProcess%myInflow(iInflow)%InnerRadius
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_OuterRadius','=',myProcess%myInflow(iInflow)%OuterRadius
      write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_center'//'=',myProcess%myInflow(iInflow)%center
      write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_normal'//'=',myProcess%myInflow(iInflow)%normal
     ELSE
      write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'nSubInflows = ',myProcess%myInflow(iInflow)%nSubInflows
      DO iSubInflow=1,myProcess%myInflow(iInflow)%nSubInflows
       write(*,'(A,I0,A)') "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Sub_',iSubinflow
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Type','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Material','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Massflowrate','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Temperature','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_InnerRadius','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%InnerRadius
       write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_OuterRadius','=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%OuterRadius
       write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_center'//'=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%center
       write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_normal'//'=',myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%normal
      END DO
     END IF
     
    END DO

    write(*,*)
    
    if (myMultiMat%nOfMaterials.gt.1) then
    
    write(*,*) "myMultiMat%nOfMaterials",'=',myMultiMat%nOfMaterials
    write(*,*) "myMultiMat%InitMaterial",'=',myMultiMat%InitMaterial
    write(*,*)
    DO iMat=1,myMultiMat%nOfMaterials

!      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%ViscoMin",'=',myMultiMat%Mat(iMat)%Rheology%ViscoMin
!      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%ViscoMax",'=',myMultiMat%Mat(iMat)%Rheology%ViscoMax
     
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.2) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Powerlaw'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%K",'=',myMultiMat%Mat(iMat)%Rheology%K
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(iMat)%Rheology%n
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.1) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Carreau'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(iMat)%Rheology%C
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.3) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Polyflow'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(iMat)%Rheology%C
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.4) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Ellis'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(iMat)%Rheology%C
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.6) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Bingham'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(iMat)%Rheology%C
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.7) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','MAS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%mu_0",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%k",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(iMat)%Rheology%C
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%a",'=',myMultiMat%Mat(iMat)%Rheology%D
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%Equation.eq.8) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','YASUDA'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%mu_0",'=',myMultiMat%Mat(iMat)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Lambda",'=',myMultiMat%Mat(iMat)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(iMat)%Rheology%C
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%a",'=',myMultiMat%Mat(iMat)%Rheology%D
     END IF
     write(*,*) 

     write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%TemperatureMin",'=',myMultiMat%Mat(iMat)%Rheology%TemperatureMin
     write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%TemperatureMax",'=',myMultiMat%Mat(iMat)%Rheology%TemperatureMax

     IF (myMultiMat%Mat(iMat)%Rheology%AtFunc.eq.1) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','ISOTHERM'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%aT",'=',1.0
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%AtFunc.eq.2) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','C1C2'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(iMat)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(iMat)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(iMat)%Rheology%Tb
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%AtFunc.eq.3) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','TBTS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(iMat)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(iMat)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(iMat)%Rheology%Tb
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Ts",'=',myMultiMat%Mat(iMat)%Rheology%Ts
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%AtFunc.eq.4) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','Melt-TBTS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(iMat)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(iMat)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(iMat)%Rheology%Tb
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Ts",'=',myMultiMat%Mat(iMat)%Rheology%Ts
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%log_aT_Tilde_Max",'=',myMultiMat%Mat(iMat)%Rheology%log_aT_Tilde_Max
     END IF
     IF (myMultiMat%Mat(iMat)%Rheology%AtFunc.eq.5) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','ETB'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%E",'=',myMultiMat%Mat(iMat)%Rheology%E
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(iMat)%Rheology%Tb
     END IF
     
     IF (myMultiMat%Mat(iMat)%Rheology%bWallSlip) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%WallSlip",'=','ACTIVATED'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_d",'=',myMultiMat%Mat(iMat)%Rheology%WS_d
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_TauMin",'=',myMultiMat%Mat(iMat)%Rheology%WS_TauMin
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_TauMax",'=',myMultiMat%Mat(iMat)%Rheology%WS_TauMax
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_SlipFactor",'=',myMultiMat%Mat(iMat)%Rheology%WS_SlipFactor
     ELSE
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%WallSlip",'=','OFF'
     END IF
     
     write(*,*)
     write(*,*) "myThermodyn%DensityModel",'=',TRIM(ADJUSTL(myMultiMat%Mat(iMat)%Thermodyn%DensityModel))
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%HeatConductivity",'=',myMultiMat%Mat(iMat)%Thermodyn%lambda
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%HeatConductivitySlope",'=',myMultiMat%Mat(iMat)%Thermodyn%lambdaSteig
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%HeatCapacity",'=',myMultiMat%Mat(iMat)%Thermodyn%cp
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%HeatCapacitySlope",'=',myMultiMat%Mat(iMat)%Thermodyn%cpSteig
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%DensityT0",'=',myMultiMat%Mat(iMat)%Thermodyn%densityT0
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%DensitySlope",'=',myMultiMat%Mat(iMat)%Thermodyn%densitySteig
     write(*,'(A,I0,A,A,ES12.4)') " myThermodyn(",iMat,")%Density",'=',myMultiMat%Mat(iMat)%Thermodyn%density
     write(*,*)
    END DO
  
    else
     
     iMat = 1
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.2) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Powerlaw'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%K",'=',myMultiMat%Mat(1)%Rheology%K
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(1)%Rheology%n
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.1) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Carreau'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(1)%Rheology%C
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.3) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Polyflow'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(1)%Rheology%C
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.4) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Ellis'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(1)%Rheology%C
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.6) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','Bingham'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%A",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%B",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C",'=',myMultiMat%Mat(1)%Rheology%C
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.7) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','MAS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%mu_0",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%k",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(1)%Rheology%C
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%a",'=',myMultiMat%Mat(1)%Rheology%D
     END IF
     IF (myMultiMat%Mat(1)%Rheology%Equation.eq.8) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%model",'=','YASUDA'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%mu_0",'=',myMultiMat%Mat(1)%Rheology%A
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Lambda",'=',myMultiMat%Mat(1)%Rheology%B
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%n",'=',myMultiMat%Mat(1)%Rheology%C
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%a",'=',myMultiMat%Mat(1)%Rheology%D
     END IF
     write(*,*)

     write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%TemperatureMin",'=',myMultiMat%Mat(1)%Rheology%TemperatureMin
     write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%TemperatureMax",'=',myMultiMat%Mat(1)%Rheology%TemperatureMax

     IF (myMultiMat%Mat(1)%Rheology%AtFunc.eq.1) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','ISOTHERM'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%aT",'=',1.0
     END IF
     IF (myMultiMat%Mat(1)%Rheology%AtFunc.eq.2) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','C1C2'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(1)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(1)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(1)%Rheology%Tb
     END IF
     IF (myMultiMat%Mat(1)%Rheology%AtFunc.eq.3) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','TBTS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(1)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(1)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(1)%Rheology%Tb
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Ts",'=',myMultiMat%Mat(1)%Rheology%Ts
     END IF
     IF (myMultiMat%Mat(1)%Rheology%AtFunc.eq.4) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','Melt-TBTS'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C1",'=',myMultiMat%Mat(1)%Rheology%C1
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%C2",'=',myMultiMat%Mat(1)%Rheology%C2
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(1)%Rheology%Tb
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Ts",'=',myMultiMat%Mat(1)%Rheology%Ts
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%log_aT_Tilde_Max",'=',myMultiMat%Mat(1)%Rheology%log_aT_Tilde_Max
     END IF
     IF (myMultiMat%Mat(1)%Rheology%AtFunc.eq.5) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%TempModel",'=','ETB'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%E",'=',myMultiMat%Mat(1)%Rheology%E
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%Tb",'=',myMultiMat%Mat(1)%Rheology%Tb
     END IF
     write(*,*)

     IF (myMultiMat%Mat(iMat)%Rheology%bWallSlip) THEN
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%WallSlip",'=','ACTIVATED'
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_d",'=',myMultiMat%Mat(iMat)%Rheology%WS_d
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_TauMin",'=',myMultiMat%Mat(iMat)%Rheology%WS_TauMin
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_TauMax",'=',myMultiMat%Mat(iMat)%Rheology%WS_TauMax
      write(*,'(A,I0,A,A,ES12.4)') " myRheology(",iMat,")%WS_SlipFactor",'=',myMultiMat%Mat(iMat)%Rheology%WS_SlipFactor
     ELSE
      write(*,'(A,I0,A,A,A)') " myRheology(",iMat,")%WallSlip",'=','OFF'
     END IF
     
     write(*,*)
    write(*,*) "myThermodyn%DensityModel",'=',TRIM(ADJUSTL(myThermodyn%DensityModel))
    write(*,*) "myThermodyn%HeatConductivity",'=',myThermodyn%lambda
    write(*,*) "myThermodyn%HeatConductivitySlope",'=',myThermodyn%lambdaSteig
    write(*,*) "myThermodyn%HeatCapacity",'=',myThermodyn%cp
    write(*,*) "myThermodyn%HeatCapacitySlope",'=',myThermodyn%cpSteig
    write(*,*) "myThermodyn%DensityT0",'=',myThermodyn%densityT0
    write(*,*) "myThermodyn%DensitySlope",'=',myThermodyn%densitySteig
    write(*,*) "myThermodyn%Density",'=',myThermodyn%density
    write(*,*) 
    
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2020.11.20 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,'(A,F10.2,A)') "Material properties interpolated to Material temperature",myProcess%T0,"C"
    write(*,'(A,A,F12.4)') "Density_[g/cm3]",'=',myThermodyn%density
    write(*,'(A,A,F12.4)') "Lambda_[W//m/K]",'=',myThermodyn%lambda
    write(*,'(A,A,F12.4)') "Cp_[kJ/kg/K]",'=',myThermodyn%cp
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PATCH 2020.11.20 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    write(*,*) 
!     write(*,*) "mySetup%MeshQuality",'=',mySetup%MeshResolution
    write(*,*) "myOutput%nOf1DLayers = "      ,myOutput%nOf1DLayers
    write(*,*) "myOutput%nOfHistogramBins = " ,myOutput%nOfHistogramBins
    write(*,*) "myOutput%HistogramShearMax = ",myOutput%HistogramShearMax    
    write(*,*) "myOutput%HistogramShearMin = ",myOutput%HistogramShearMin    
    write(*,*) "myOutput%HistogramViscoMax = ",myOutput%HistogramViscoMax    
    write(*,*) "myOutput%HistogramViscoMin = ",myOutput%HistogramViscoMin    
    write(*,*) "myOutput%CutDtata_1D = ",myOutput%CutDtata_1D
    
    write(*,*) 
    
    write(*,*) "E3DProcessParameters@nOfDIESensors = ", myProcess%nOfDIESensors
    if (myProcess%nOfDIESensors.gt.0) then
     do i=1,myProcess%nOfDIESensors
      write(*,'(A,i0,A,3ES12.4)') "E3DProcessParameters@DIESensor[",i,"]Center = ", myProcess%mySensor(i)%Center
      write(*,'(A,i0,A,3ES12.4)') "E3DProcessParameters@DIESensor[",i,"]Radius = ", myProcess%mySensor(i)%Radius
      write(*,'(A,i0,A,A)') "E3DProcessParameters@DIESensor[",i,"]RegType = ", adjustl(trim(myProcess%mySensor(i)%Type))
      write(*,'(A,i0,A,ES12.4)') "E3DProcessParameters@DIESensor[",i,"]SetValue = ", myProcess%mySensor(i)%PID_Ctrl%T_set
      write(*,'(A,i0,A,I8)') "E3DProcessParameters@DIESensor[",i,"]SetSegment = ",myProcess%mySensor(i)%iSeg
      write(*,'(A,i0,A,3ES12.4)') "E3DProcessParameters@DIESensor[",i,"]PID = ", myProcess%mySensor(i)%PID_Ctrl%omega_P, myProcess%mySensor(i)%PID_Ctrl%omega_I, myProcess%mySensor(i)%PID_Ctrl%omega_D
     end do
    else
     write(*,*) "no sensor data to report"
    end if
    
    write(*,'(A,F12.4)') "mySetup%PressureConvergenceTolerance = ",mySetup%PressureConvergenceTolerance

    write(*,*) "mySetup%PressureFBM = ",mySetup%bPressureFBM
    write(*,*) "mySetup%AutomaticTimeStepControl = ",mySetup%bAutomaticTimeStepControl
    IF (mySetup%bAutomaticTimeStepControl.and.ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") then
     write(*,*) "E3DProcessParameters@ExtrusionSpeed_CMpS = ", myProcess%ExtrusionSpeed 
     write(*,*) "E3DProcessParameters@ExtrusionGapSize_MM = ", myProcess%ExtrusionGapSize
     write(*,*) "E3DProcessParameters@ExtrusionGapFactor = ", myProcess%ExtrusionGapFactor
    END IF

    IF (mySetup%bRotationalFramOfReference) then
     write(*,*) "RotationalFramOfReference for (noPin-SSE) Temperature Simulation is on!"
    END IF
    write(*,*) "mySetup%CharacteristicShearRate = ",mySetup%CharacteristicShearRate
    write(*,*) "activeFBM_Z_Position = ",activeFBM_Z_Position   
    write(*,*) "TimeStepEnlargmentFactor = ",dTimeStepEnlargmentFactor   
    write(*,*) "AnalyticalShearRateRestriction = ",  mySigma%bAnalyticalShearRateRestriction

    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."OFF") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,*) "mySetup%MeshPath",'=',ADJUSTL(TRIM(mySetup%cMeshPath))
    END IF
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."HOLLOWCYLINDER") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,'(A,A,4I6)') " mySetup%Resolution",'=',mySetup%MeshResolution
     write(*,*) "mySetup%Resolution[nR,nT,nZ]",'=',mySetup%m_nR,mySetup%m_nT,mySetup%m_nZ
    END IF
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."FULLCYLINDER") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,'(A,A,4I6)') " mySetup%Resolution[nP,nT,nR,nZ]",'=',mySetup%m_nP,mySetup%m_nR,mySetup%m_nT,mySetup%m_nZ
    END IF
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."TWINSCREW") THEN
     write(*,*) "mySetup%HexMesher",'=',ADJUSTL(TRIM(mySetup%cMesher))
     write(*,'(A,A,4I6)') " mySetup%Resolution",'=',mySetup%MeshResolution
     write(*,*) "mySetup%Resolution[nR,nT1,nT2,nZ]",'=',mySetup%m_nR,mySetup%m_nT1,mySetup%m_nT2,mySetup%m_nZ
    END IF
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."BOX") THEN
     IF (mySetup%nBoxElem.gt.0) write(*,'(A,I0)') "myMesh%nBoxElem=",mySetup%nBoxElem
     IF (mySetup%m_nX.gt.0.and.mySetup%m_nY.gt.0.and.mySetup%m_nZ.gt.0) write(*,'(A,3I0)') &
     "myMesh%nBoxElemXYZ=",mySetup%m_nX,mySetup%m_nY,mySetup%m_nZ
     IF (mySetup%MeshResolution.gt.0) write(*,'(A,I0)') "myMesh%MeshResolution=",mySetup%MeshResolution
     write(*,'(A,2ES12.4)') "myMesh%Extent_X=",mySetup%m_box(1,:)
     write(*,'(A,2ES12.4)') "myMesh%Extent_Y=",mySetup%m_box(2,:)
     write(*,'(A,2ES12.4)') "myMesh%Extent_Z=",mySetup%m_box(3,:)
    END IF
    write(*,*)
    write(*,*) "myProcess%dAlpha",'=',myProcess%dAlpha
    write(*,*) "myProcess%Periodicity",'=',myProcess%Periodicity
    write(*,*) "myProcess%nTimeLevels",'=',myProcess%nTimeLevels
!     write(*,*) "mySetup%bSendEmail",'=',mySetup%bSendEmail
    write(*,*) "myProcess%Angle",'=',myProcess%Angle
    write(*,*) "myProcess%Phase",'=',myProcess%Phase
!     write(*,*) "=========================================================================="
    END IF
! 
    if (mySigma%NumberOfSeg.gt.0) THEN
     mySigma%mySegment(1)%StartAlpha = 0d0 !myPI/0d0
     mySigma%SegmentLength = mySigma%mySegment(1)%L
     IF (mySigma%NumberOfSeg.GE.2) THEN
     DO iSeg=2,mySigma%NumberOfSeg
       IF (mySigma%mySegment(iSeg-1)%ART.EQ.'FOERD'.OR.mySigma%mySegment(iSeg-1)%ART.EQ.'SME') THEN
        mySigma%mySegment(iSeg)%StartAlpha = mySigma%mySegment(iSeg-1)%StartAlpha + mySigma%mySegment(iSeg-1)%L/mySigma%mySegment(iSeg-1)%t*2d0*myPI
       END IF
       IF (mySigma%mySegment(iSeg-1)%ART.EQ.'KNET'.or.mySigma%mySegment(iSeg-1)%ART.EQ.'SKNET'.or.mySigma%mySegment(iSeg-1)%ART.EQ.'EKNET') THEN
        mySigma%mySegment(iSeg)%StartAlpha =  mySigma%mySegment(iSeg-1)%StartAlpha + myPI*DBLE(mySigma%mySegment(iSeg-1)%N-1)*mySigma%mySegment(iSeg-1)%Alpha/1.8d2
        IF (mySigma%mySegment(iSeg)%ART.EQ.'KNET'.or.mySigma%mySegment(iSeg)%ART.EQ.'SKNET'.or.mySigma%mySegment(iSeg-1)%ART.EQ.'EKNET') THEN
         mySigma%mySegment(iSeg)%StartAlpha = mySigma%mySegment(iSeg)%StartAlpha + myPI*mySigma%mySegment(iSeg)%Alpha/1.8d2
        END IF
       END IF
       mySigma%SegmentLength = mySigma%SegmentLength + mySigma%mySegment(iSeg)%L
     END DO
     END IF
    END IF

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
    IF (ieee_is_finite(myProcess%dPress)) THEN
      dZPeriodicLength = mySigma%L
      dPeriodicity     = [1d9,1d9,mySigma%L]
!       myProcess%dPress = 1d6*myProcess%dPress
      bNoOutflow = .TRUE.
    END IF
    
    IF (ieee_is_finite(myProcess%Massestrom)) THEN
      bNoOutflow = .FALSE.
!      bNoOutflow = .FALSE.
      dZPeriodicLength = 1d5*mySigma%L
    END IF

    IF (myid.eq.1.or.subnodes.eq.0) then
     write(*,*) "myProcess%dZPeriodicLength",'=',dZPeriodicLength
     write(*,*) "myProcess%dPress",'=',myProcess%dPress
     write(*,*) "myProcess%NoOutFlow",'=',bNoOutFlow
     write(*,*) "=========================================================================="
    END IF

!    ! Make some output on the terminal
!     call inip_info(parameterlist)

!     ! Write it into a file
!      call inip_dumpToFile(parameterlist,"test_dump_parlst.dat",INIP_REPLACE)

    ! Clean up the parameterlist
    call inip_done(parameterlist)

    if (bReadError) then
      write(*,*) 'Error during reading the file ', trim(adjustl(ce3dfile)), '. Stopping. See output above.'
      CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

    end if

    CONTAINS
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE FillUpRheoAndThermoDataForPseudoAir(r,t)
    TYPE(tRheology) :: r
    TYPE(tThermodyn) :: t
    
     t%Lambda = 0.10D0
     t%LambdaSteig = 0.0D0
     t%Cp = 2.25D0
     t%CpSteig = 0.0D0

     t%DensityModel = "DENSITY"
     t%DensityT0    = 1d-2
     t%DensitySteig = 0d0
     
     r%Equation = 1
     r%A = 1.0d0
     r%B = 1.0d0
     r%C = 0.0d0
     
     r%AtFunc = 1
    
    END SUBROUTINE FillUpRheoAndThermoDataForPseudoAir
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE FillUpTempBCs(nT,cINI)
    implicit none
    INTEGER :: nT   
    character(len=INIP_STRLEN) cINI,cUnit
    real*8 daux
    
    DO iTempBC=1,myProcess%nOfTempBCs
    
     WRITE(cTempBC_i,'(A,A,I0)') ADJUSTL(TRIM(cINI)),"/TempBC_",iTempBC
     if (myid.eq.1) write(*,*) "|",ADJUSTL(TRIM(cTempBC_i)),"|"
    
     call INIP_getvalue_string(parameterlist,cTempBC_i,"Unit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0
     
     call INIP_getvalue_double(parameterlist,cTempBC_i,"ReferenceTemperature",myProcess%myTempBC(iTempBC)%temperature,myInf)
     call INIP_getvalue_string(parameterlist,cTempBC_i,"center",cCenter,'0d0, 0d0, 0d0')
     read(cCenter,*,err=55) myProcess%myTempBC(iTempBC)%Center
     myProcess%myTempBC(iTempBC)%Center = daux*myProcess%myTempBC(iTempBC)%Center
     call INIP_getvalue_string(parameterlist,cTempBC_i,"gradient",cNormal,'OFF')
     call inip_toupper_replace(cNormal) 
     IF (ADJUSTL(TRIM(cNormal)).eq."OFF") THEN
      myProcess%myTempBC(iTempBC)%bGradient = .false.
     ELSE
      read(cNormal,*,err=56) myProcess%myTempBC(iTempBC)%Gradient
      myProcess%myTempBC(iTempBC)%Gradient = myProcess%myTempBC(iTempBC)%Gradient/daux
      myProcess%myTempBC(iTempBC)%bGradient = .true.
     END IF

     call INIP_getvalue_string(parameterlist,cTempBC_i,"LinSphereGradient",cNormal,'OFF')
     call inip_toupper_replace(cNormal) 
     IF (ADJUSTL(TRIM(cNormal)).eq."OFF") THEN
      myProcess%myTempBC(iTempBC)%bLinSphereGradient = .false.
     ELSE
      read(cNormal,*,err=57) myProcess%myTempBC(iTempBC)%LinSphereGradient
      myProcess%myTempBC(iTempBC)%LinSphereGradient = myProcess%myTempBC(iTempBC)%LinSphereGradient/daux
      myProcess%myTempBC(iTempBC)%bLinSphereGradient = .true.
     END IF
     
     call INIP_getvalue_string(parameterlist,cTempBC_i,"QuadSphereGradient",cNormal,'OFF')
     call inip_toupper_replace(cNormal) 
     IF (ADJUSTL(TRIM(cNormal)).eq."OFF") THEN
      myProcess%myTempBC(iTempBC)%bQuadSphereGradient = .false.
     ELSE
      read(cNormal,*,err=58) myProcess%myTempBC(iTempBC)%QuadSphereGradient
      myProcess%myTempBC(iTempBC)%QuadSphereGradient = myProcess%myTempBC(iTempBC)%QuadSphereGradient/daux
      myProcess%myTempBC(iTempBC)%bQuadSphereGradient = .true.
     END IF

     GOTO 80     
55   write(*,*) 'WRONGLY DEFINED Center for TempBC',iTempBC,' !!'//"|",ADJUSTL(TRIM(cCenter)),"|"
     GOTO 80
56   write(*,*) 'WRONGLY DEFINED Gradient for TempBC',iTempBC,' !!'//"|",ADJUSTL(TRIM(cNormal)),"|"  
     GOTO 80
57   write(*,*) 'WRONGLY DEFINED LinSphereGradient for TempBC',iTempBC,' !!'//"|",ADJUSTL(TRIM(cNormal)),"|"  
     GOTO 80
58   write(*,*) 'WRONGLY DEFINED QuadSphereGradient for TempBC',iTempBC,' !!'//"|",ADJUSTL(TRIM(cNormal)),"|"  
     GOTO 80
80   CONTINUE

    END DO

    END SUBROUTINE FillUpTempBCs

    SUBROUTINE FillUpInflows(nT,cINI)
    implicit none
    INTEGER :: nT   
    character(len=INIP_STRLEN) cINI,cUnit
    real*8 daux
    
    DO iInflow=1,myProcess%nOfInflows
!      if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'E3DProcessParameters/Inflow_',iInflow
!      if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'E3DProcessParameters/Inflow_',iInflow
!      if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'E3DProcessParameters/Inflow_',iInflow
     
     WRITE(cInflow_i,'(A,A,I0)') ADJUSTL(TRIM(cINI)),"/Inflow_",iInflow
     
     if (myid.eq.1) write(*,*) "|",ADJUSTL(TRIM(cInflow_i)),"|"
     
     call INIP_getvalue_string(parameterlist,cInflow_i,"Unit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0

     call INIP_getvalue_int(parameterlist,cInflow_i,"nSubInflows",myProcess%myInflow(iInflow)%nSubInflows,0)
     
     IF (myProcess%myInflow(iInflow)%nSubInflows.eq.0) then
     
      call INIP_getvalue_int(parameterlist,cInflow_i,"Material",myProcess%myInflow(iInflow)%Material,0)
      if (myProcess%myInflow(iInflow)%Material.eq.0) write(*,*) 'UNDEFINED material from Inflow ',iInflow,' !!'

      call INIP_getvalue_string(parameterlist,cInflow_i,"Type",cBCtype,'unknown')
      call inip_toupper_replace(cBCtype)
      myProcess%myInflow(iInflow)%iBCtype = 0
      IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA1") THEN
       myProcess%myInflow(iInflow)%iBCtype = 1
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA2") THEN
       myProcess%myInflow(iInflow)%iBCtype = 2
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."FLAT") THEN
       myProcess%myInflow(iInflow)%iBCtype = 3
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDFLAT") THEN
       myProcess%myInflow(iInflow)%iBCtype = 4
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."RECTANGLE") THEN
       myProcess%myInflow(iInflow)%iBCtype = 5
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDRECTANGLE") THEN
       myProcess%myInflow(iInflow)%iBCtype = 6
      END IF
      if (myProcess%myInflow(iInflow)%iBCtype.eq.0) then
       write(*,*) 'UNDEFINED Inflow type!!'
      end if
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"massflowrate",myProcess%myInflow(iInflow)%massflowrate,myInf)
      if (myProcess%myInflow(iInflow)%massflowrate.eq.myInf) write(*,*) 'UNDEFINED massflowrate through Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"temperature",myProcess%myInflow(iInflow)%temperature,myInf)
      if (myProcess%myInflow(iInflow)%temperature.eq.myInf) write(*,*) 'UNDEFINED temperature for Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"TemperatureRange",myProcess%myInflow(iInflow)%TemperatureRange,0d0)

      call INIP_getvalue_string(parameterlist,cInflow_i,"TemperatureType",cBCtype,'CNST') ! Default is cnst temperature (0)
      call inip_toupper_replace(cBCtype)
      myProcess%myInflow(iInflow)%temperatureType = 0
      IF (ADJUSTL(TRIM(cBCtype)).eq."LINEAR") THEN
       myProcess%myInflow(iInflow)%temperatureType = 1
      END IF
      IF (ADJUSTL(TRIM(cBCtype)).eq."QUADRATIC") THEN
       myProcess%myInflow(iInflow)%temperatureType = 2
      END IF
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"innerradius",myProcess%myInflow(iInflow)%innerradius,myInf)
      myProcess%myInflow(iInflow)%innerradius = daux*myProcess%myInflow(iInflow)%innerradius
      if (myProcess%myInflow(iInflow)%innerradius.eq.myInf) write(*,*) 'UNDEFINED inner radius for Inflow',iInflow,' !!'
      
      call INIP_getvalue_double(parameterlist,cInflow_i,"outerradius",myProcess%myInflow(iInflow)%outerradius,myInf)
      myProcess%myInflow(iInflow)%outerradius = daux*myProcess%myInflow(iInflow)%outerradius
      if (myProcess%myInflow(iInflow)%outerradius.eq.myInf) write(*,*) 'UNDEFINED outer radius for Inflow',iInflow,' !!'
      
      call INIP_getvalue_string(parameterlist,cInflow_i,"center",cCenter,'unknown')
      call INIP_getvalue_string(parameterlist,cInflow_i,"normal",cNormal,'unknown')
      read(cCenter,*,err=55) myProcess%myInflow(iInflow)%Center
      myProcess%myInflow(iInflow)%Center = daux*myProcess%myInflow(iInflow)%Center
      read(cNormal,*,err=56) myProcess%myInflow(iInflow)%Normal

      call INIP_getvalue_string(parameterlist,cInflow_i,"midpointA",cMidpointA,'unknown')
      call INIP_getvalue_string(parameterlist,cInflow_i,"midpointB",cMidpointB,'unknown')
      read(cMidpointA,*,err=57) myProcess%myInflow(iInflow)%MidpointA
      myProcess%myInflow(iInflow)%MidpointA = daux*myProcess%myInflow(iInflow)%MidpointA
      read(cMidpointB,*,err=57) myProcess%myInflow(iInflow)%MidpointB
      myProcess%myInflow(iInflow)%MidpointB = daux*myProcess%myInflow(iInflow)%MidpointB
     ELSE
      
      ALLOCATE(myProcess%myInflow(iInflow)%mySubInflow(myProcess%myInflow(iInflow)%nSubInflows))
      
      DO iSubInflow = 1,myProcess%myInflow(iInflow)%nSubInflows
       WRITE(cInflow_i,'(A,A,I0,A,I0)') ADJUSTL(TRIM(cINI)),"/Inflow_",iInflow,"/Sub_",iSubInflow
!       WRITE(*,*)"'",ADJUSTL(TRIM(cInflow_i)),"'"
      
       call INIP_getvalue_int(parameterlist,cInflow_i,"Material",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material,0)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Material.eq.0) write(*,*) 'UNDEFINED material from Inflow ',iInflow,' !!'

       call INIP_getvalue_string(parameterlist,cInflow_i,"Type",cBCtype,'unknown')
       call inip_toupper_replace(cBCtype)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 0
       IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA1") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 1
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."ROTATEDPARABOLA2") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 2
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."FLAT") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 3
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDFLAT") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 4
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."RECTANGLE") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 5
       END IF
       IF (ADJUSTL(TRIM(cBCtype)).eq."CURVEDRECTANGLE") THEN
        myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype = 6
       END IF
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%iBCtype.eq.0) then
        write(*,*) 'UNDEFINED Inflow type!!'
       end if
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"massflowrate",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate,myInf)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%massflowrate.eq.myInf) write(*,*) 'UNDEFINED massflowrate through Inflow',iInflow,' !!'
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"temperature",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature,myInf)
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%temperature.eq.myInf) write(*,*) 'UNDEFINED temperature for Inflow',iInflow,' !!'

       call INIP_getvalue_double(parameterlist,cInflow_i,"innerradius",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius,myInf)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%innerradius.eq.myInf) write(*,*) 'UNDEFINED inner radius for Inflow',iInflow,' !!'
       
       call INIP_getvalue_double(parameterlist,cInflow_i,"outerradius",myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius,myInf)
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius
       if (myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%outerradius.eq.myInf) write(*,*) 'UNDEFINED outer radius for Inflow',iInflow,' !!'
       
       call INIP_getvalue_string(parameterlist,cInflow_i,"center",cCenter,'unknown')
       call INIP_getvalue_string(parameterlist,cInflow_i,"normal",cNormal,'unknown')
       read(cCenter,*,err=55) myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center
       myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center = daux*myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Center
       read(cNormal,*,err=56) myProcess%myInflow(iInflow)%mySubInflow(iSubInflow)%Normal
       
      END DO
      
     END IF
     GOTO 57     
55   write(*,*) 'WRONGLY DEFINED center for Inflow',iInflow,' !!'//"|",ADJUSTL(TRIM(cCenter)),"|"
     GOTO 57
56   write(*,*) 'WRONGLY DEFINED normal for Inflow',iInflow,' !!'//"|",ADJUSTL(TRIM(cNormal)),"|"  
     GOTO 57
57   CONTINUE
     
    END DO
    
    END SUBROUTINE FillUpInflows    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE FillUpThermoData(t,cINI)
    TYPE(tThermodyn) :: t
    character(len=INIP_STRLEN) cINI
    
    call INIP_getvalue_double(parameterlist,cINI,"HeatConductivity",t%lambda,myInf)
    call INIP_getvalue_double(parameterlist,cINI,"HeatConductivitySlope",t%lambdaSteig,0d0)
    call INIP_getvalue_double(parameterlist,cINI,"HeatCapacity",t%cp,myInf)
    call INIP_getvalue_double(parameterlist,cINI,"HeatCapacitySlope",t%CpSteig,0d0)

    call INIP_getvalue_string(parameterlist,cINI,"DensityModel", t%DensityModel,'no')
    call inip_toupper_replace(t%DensityModel)
    t%density=myInf
    IF (ADJUSTL(TRIM(t%DensityModel)).eq."DENSITY") THEN
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Density","Density",t%DensityT0,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Density","DensitySlope",t%DensitySteig,0d0)
    END IF
    IF (ADJUSTL(TRIM(t%DensityModel)).eq."SPECVOLUME") THEN
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/SpecVolume","SpecVolume",t%DensityT0,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/SpecVolume","SpecVolumeSlope",t%DensitySteig,0d0)
!       t%Density = 1d0/t%Density
!       t%DensitySteig = 1d0/t%DensitySteig
    END IF
    IF (myThermodyn%densityT0.eq.myinf) THEN
     WRITE(*,*) "density is not defined"
     WRITE(*,*) '"',TRIM(myThermodyn%DensityModel),'"'
     bReadError=.TRUE.
    END IF
   
    END SUBROUTINE FillUpThermoData
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SUBROUTINE FillUpRheoData(t,cINI)
    TYPE(tRheology) :: t
    character(len=INIP_STRLEN) cINI
    
    t%Equation = 0
    call INIP_getvalue_string(parameterlist,ADJUSTL(TRIM(cINI)),"CalcVisco", cRheology,'NoRheology')
    call inip_toupper_replace(cRheology)
    IF (ADJUSTL(TRIM(cRheology)).eq."BINGHAM") THEN
      t%Equation = 6
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Bingham","ZeroViscosity",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Bingham","Regularization",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Bingham","YieldStress",t%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."CARREAU") THEN
      t%Equation = 1
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Carreau","ZeroViscosity",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Carreau","RecipVelocity",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Carreau","Exponent",t%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POWERLAW".OR.ADJUSTL(TRIM(cRheology)).eq."POTENZ".OR.ADJUSTL(TRIM(cRheology)).eq."POWER") THEN
      t%Equation = 2
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"Consistence", t%K,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"Exponent",t%n,myInf)
!       if (t%K.eq.myInf) then
!        call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Power","Consistence", t%K,myInf)
!       end if
!       if (t%n.eq.myInf) then
!        call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Power","Exponent",t%n,myInf)
!       end if
      
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POLYFLOW") THEN
      t%Equation = 3
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Polyflow","Polyflow_A",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Polyflow","Polyflow_B",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Polyflow","Polyflow_C",t%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."ELLIS") THEN
      t%Equation = 4
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Ellis","ZeroViscosity",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Ellis","Gamma0",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/Ellis","Exponent",t%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."MAS") THEN
      t%Equation = 7
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/MAS","MAS_mu0",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/MAS","MAS_k",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/MAS","MAS_n",t%C,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/MAS","MAS_a",t%D,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."YASUDA") THEN
      t%Equation = 8
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/YASUDA","ZeroViscosity",t%A,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/YASUDA","RecipVelocity",t%B,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/YASUDA","Exponent",t%C,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/YASUDA","Plateau",t%D,myInf)
    END IF
    IF (t%Equation.eq.0) THEN
     WRITE(*,*) "no valid rheology is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
    END IF

    CAll INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI)),"TemperatureMin",t%TemperatureMin,-1d30)
    call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI)),"TemperatureMax",t%TemperatureMax,+1d30)

    t%AtFunc = 0
    cRheology = ' '
    call INIP_getvalue_string(parameterlist,ADJUSTL(TRIM(cINI)),"CalcTemp", cRheology,'TEMPERATUREINDEPENDENT')
    call inip_toupper_replace(cRheology)
    
    IF (ADJUSTL(TRIM(cRheology)).eq."ISOTHERM".or.ADJUSTL(TRIM(cRheology)).eq."TEMPERATUREINDEPENDENT") THEN
      t%AtFunc = 1
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."C1C2") THEN
      t%AtFunc = 2
!       t%Tb = 165d0
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/C1C2","ReferenceTemperature",t%TB,165d0)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/C1C2","C1",t%C1,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/C1C2","C2",t%C2,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."TBTS") THEN
      t%AtFunc = 3
!       t%C1 = 8.86d0
!       t%C2 = 101.6d0
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/TBTS","C1",t%C1,8.86d0)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/TBTS","C2",t%C2,101.6d0)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/TBTS","ReferenceTemperature",t%Tb,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/TBTS","StandardTemperature",t%Ts,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."MELTTBTS".or.ADJUSTL(TRIM(cRheology)).eq."MELT-TBTS") THEN
      t%AtFunc = 4
!       t%C1 = 8.86d0
!       t%C2 = 101.6d0
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"C1",t%C1,8.86d0)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"C2",t%C2,101.6d0)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"ReferenceTemperature",t%Tb,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"StandardTemperature",t%Ts,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/"//ADJUSTL(TRIM(cRheology)),"log_aT_Tilde_Max",t%log_aT_Tilde_Max,myInf)
    END IF
    !ETB myRheology%E is in J/mol
    IF (ADJUSTL(TRIM(cRheology)).eq."ETB") THEN
      t%AtFunc = 5
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/ETB","ActivatingEnergy",t%E,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/ETB","ReferenceTemperature",t%TB,myInf)
    END IF
    IF (t%AtFunc.eq.0) THEN
     WRITE(*,*) "no temperature correction is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
    END IF
    
    cRheology = ' '
    call INIP_getvalue_string(parameterlist,ADJUSTL(TRIM(cINI)),"WallSlip", cRheology,'NO')
    call inip_toupper_replace(cRheology)
    IF (ADJUSTL(TRIM(cRheology)).eq."NO") THEN
      t%bWallSlip = .FALSE.
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."YES".OR.ADJUSTL(TRIM(cRheology)).eq."ON".or.ADJUSTL(TRIM(cRheology)).eq."ACTIVATED") THEN
      t%bWallSlip = .TRUE.
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/WallSlip","WS_d",t%WS_d,1d0)
      t%WS_d = dSizeScale*t%WS_d
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/WallSlip","WS_Tau0",t%WS_Tau0,myInf)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/WallSlip","WS_TauD",t%WS_TauD,10.0d3)
      call INIP_getvalue_double(parameterlist,ADJUSTL(TRIM(cINI))//"/WallSlip","WS_SlipFactor",t%WS_SlipFactor,myInf)
      
      IF (t%WS_Tau0.eq.myinf) THEN
       WRITE(*,*) "WallShearStress for WallSlip is not defined"
       bReadError=.TRUE.
      ELSE
       t%WS_TauMin = t%WS_Tau0 - 0.5d0*t%WS_TauD
       t%WS_TauMax = t%WS_Tau0 + 0.5d0*t%WS_TauD
      END IF
      
      IF (t%WS_SlipFactor.eq.myinf) THEN
       WRITE(*,*) "SlipFactor for WallSlip is not defined"
       bReadError=.TRUE.
      END IF
    END IF
    
    END SUBROUTINE FillUpRheoData
    
    SUBROUTINE ReadDoubleFromDimensionalString()
    IMPLICIT NONE
    INTEGER i,n,i1,i2
    LOGICAL :: bOK

    n = len(cText)
    i1 = 0
    i2 = 0
    DO i=1,n
     IF (cText(i:i).EQ. '[') i1 = i
     IF (cText(i:i).EQ. ']') i2 = i
    END DO

    bOk=.FALSE.
    IF (i1.ne.0.AND.i2.ne.0) bOk=.TRUE.
    
    IF (bOK) THEN
     read(cText(1:i1-1),*)     dExtract_Val
     read(cText(i1+1:i2-1),*) sExtract_Dim
    ELSE
     read(cText(1:),*)     dExtract_Val
     sExtract_Dim = 'BAR'
    END IF
!     write(*,*) 'hoho',i1,i2,adjustl(trim(cText))
!     write(*,*) dExtract_Val,'"'//adjustl(trim(sExtract_Dim))//'"'
!     pause


    END SUBROUTINE ReadDoubleFromDimensionalString

    
    SUBROUTINE ExtractNomOfCharFromString(sL,j)
    CHARACTER*(*) :: sL
    CHARACTER*(256) saux
    INTEGER i,n,i1,i2,i3,j,nS
    
    n = len(sL)
    
    j=0
    i2 = 1
    DO i=1,n
     IF (sL(i:i).EQ. ',') THEN
      j = j + 1
      
      i3 = i
      READ(sL(i2:i3-1),'(A)') saux
      
      saux = trim(adjustl(saux))
      CALL inip_toupper_replace(saux)
      nS = LEN(trim(saux))
      IF (saux(nS-3:nS).ne.'.OFF') GOTO 10      
      
      i2 = i3+1
     END IF
    END DO
    
    j = j + 1
    i3 = i
    READ(sL(i2:i3-1),'(A)') saux
    saux = trim(adjustl(saux))
    CALL inip_toupper_replace(saux)
    nS = LEN(trim(saux))
    IF (saux(nS-3:nS).ne.'.OFF') GOTO 10      
    
10  CONTINUE
     
    RETURN
    
    WRITE(*,*) 'The file in the list is not a OFF file ... ',trim(saux)
    WRITE(*,*) 'Program stops ... '
    CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

    
    
    END SUBROUTINE ExtractNomOfCharFromString
    
    
    SUBROUTINE ExtractCharArrayFromString(sL,sS)
    IMPLICIT NONE
    CHARACTER*(*) :: sS(:),sL
    CHARACTER*(256) saux
    INTEGER i,n,i1,i2,i3,j,nS
    
    n = len(sL)
    
    j=0
    i2 = 1
    DO i=1,n
     IF (sL(i:i).EQ. ',') THEN
      j = j + 1
      
      i3 = i
      READ(sL(i2:i3-1),'(A)') saux
      sS(j) = ADJUSTL(TRIM(saux))
      if (myid.eq.1) write(*,*) "'"//trim(sS(j))//"'"
      
      saux = sS(j)
      CALL inip_toupper_replace(saux)
      nS = LEN(trim(saux))
      
      i2 = i3+1
     END IF
    END DO
    
    j = j + 1
    i3 = i
    READ(sL(i2:i3-1),'(A)') saux
    sS(j) = ADJUSTL(TRIM(saux))
    if (myid.eq.1) write(*,*) "'"//trim(sS(j))//"'"
    
    saux = sS(j)
    CALL inip_toupper_replace(saux)
    nS = LEN(trim(saux))
    
    return

    END SUBROUTINE ExtractCharArrayFromString
    
    
    end Subroutine ReadS3Dfile
!
! ---------------------------------------------------------------------------------------------
!
    Subroutine ReadEWIKONfile(cE3Dfile)
    use iniparser
    use var_QuadScalar
    use Sigma_User
    USE PP3D_MPI, ONLY:MPI_COMM_WORLD

    use, intrinsic :: ieee_arithmetic

    implicit none

    character(len=*), intent(in) :: cE3Dfile
    logical :: bReadError=.FALSE.
    integer :: i,iSeg,iFile,iaux,ierr

    real*8 :: myPI = dATAN(1d0)*4d0
    character(len=INIP_STRLEN) cCut,cElement_i,cElemType,cKindOfConveying,cTemperature,cConvergenceEstimator
    character(len=INIP_STRLEN) cProcessType,cRotation,cRheology,CDensity,cMeshQuality,cKTP,cUnit,sCoorString
    
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    real*8 :: myInf,dSizeScale,daux
    integer iMat

    if(ieee_support_inf(myInf))then
      myInf = ieee_value(myInf, ieee_negative_inf)
    endif


    call inip_output_init(myid,showid,unitProtfile,unitTerminal)

    ! Init the parameterlist
    call inip_init(parameterlist)

!     READ (myFile,*) myProcess%Angle,myProcess%Phase

    call inip_readfromfile(parameterlist,ADJUSTL(TRIM(cE3Dfile)))

    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","Unit",cUnit,'CM')
    call inip_toupper_replace(cUnit)
    IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
      WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
      cUnit = 'CM'
    END IF
    if (TRIM(cUnit).eq.'MM') dSizeScale = 0.100d0
    if (TRIM(cUnit).eq.'CM') dSizeScale = 1.000d0
    if (TRIM(cUnit).eq.'DM') dSizeScale = 10.00d0
    if (TRIM(cUnit).eq.'M')  dSizeScale = 100.0d0
    
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","Type",mySigma%cType,'HEAT')
    call inip_toupper_replace(mySigma%cType)
    IF (.NOT.(ADJUSTL(TRIM(mySigma%cType)).EQ."HEAT")) THEN
     WRITE(*,*) "not a valid Machine type:", ADJUSTL(TRIM(mySigma%cType))
    END IF

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelDiameter", mySigma%Dz_out ,myInf)
    mySigma%Dz_out = dSizeScale*mySigma%Dz_out
    DistTolerance = 1d0*mySigma%Dz_Out

    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfElements", mySigma%NumberOfSeg ,0)

    IF (mySigma%NumberOfSeg.ge.1) THEN
     ALLOCATE (mySigma%mySegment(mySigma%NumberOfSeg))
    ELSE
     WRITE(*,*) "not a valid number of segments"
     WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
    ENDIF

    DO iSeg=1,mySigma%NumberOfSeg
     WRITE(cElement_i,'(A,I0)') 'E3DGeometryData/Machine/Element_',iSeg

     call INIP_getvalue_string(parameterlist,cElement_i,"ObjectType",mySigma%mySegment(iSeg)%ObjectType)
     call inip_toupper_replace(mySigma%mySegment(iSeg)%ObjectType)
     IF (.NOT.(TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE'.OR.TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'BLOCK'.OR.&
               TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'MELT')) THEN
       WRITE(*,*) "STL object type is invalid. Only BLOCK, MELT or WIRE types are allowed: '"//TRIM(mySigma%mySegment(iSeg)%ObjectType)//"'"
     END IF

     call INIP_getvalue_int(parameterlist,cElement_i,"MaterialIndex", mySigma%mySegment(iSeg)%MatInd ,0)
     
     call INIP_getvalue_double(parameterlist,cElement_i,"InitialTemperature", mySigma%mySegment(iSeg)%InitTemp ,0d0)

     call INIP_getvalue_double(parameterlist,cElement_i,"VolumetricHeatSourceMax", mySigma%mySegment(iSeg)%HeatSourceMax ,0d0)
     call INIP_getvalue_double(parameterlist,cElement_i,"VolumetricHeatSourceMin", mySigma%mySegment(iSeg)%HeatSourceMin ,0d0)
     mySigma%mySegment(iSeg)%UseHeatSource = mySigma%mySegment(iSeg)%HeatSourceMax

     call INIP_getvalue_string(parameterlist,cElement_i,"TemperatureBC",mySigma%mySegment(iSeg)%TemperatureBC,'NO')
     call inip_toupper_replace(mySigma%mySegment(iSeg)%TemperatureBC)
     if (.NOT.(mySigma%mySegment(iSeg)%TemperatureBC.eq.'CONSTANT'.or.&
               mySigma%mySegment(iSeg)%TemperatureBC.eq.'FULLCONSTANT'.or.&
               mySigma%mySegment(iSeg)%TemperatureBC.eq.'FLUX'.or.&
               mySigma%mySegment(iSeg)%TemperatureBC.eq.'NO')) THEN
       WRITE(*,*) "Undefined thermal condition: ",mySigma%mySegment(iSeg)%TemperatureBC
     END IF
     

     call INIP_getvalue_string(parameterlist,cElement_i,"Unit",mySigma%mySegment(iSeg)%Unit,'CM')
     call inip_toupper_replace(mySigma%mySegment(iSeg)%Unit)
     IF (.NOT.(mySigma%mySegment(iSeg)%Unit.eq.'MM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'CM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'DM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'M')) THEN
       WRITE(*,*) "STL unit type is invalid. Only MM, CM, DM or 'M' units are allowed",mySigma%mySegment(iSeg)%Unit
       mySigma%mySegment(iSeg)%Unit = 'CM'
     END IF
     if (TRIM(cUnit).eq.'MM') dSizeScale = 0.100d0
     if (TRIM(cUnit).eq.'CM') dSizeScale = 1.000d0
     if (TRIM(cUnit).eq.'DM') dSizeScale = 10.00d0
     if (TRIM(cUnit).eq.'M')  dSizeScale = 100.0d0
      
     IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
     
       call INIP_getvalue_double(parameterlist,cElement_i,"ConvergenceCondition", mySigma%mySegment(iSeg)%ConvergenceDetector%Condition ,0.001d0)
       call INIP_getvalue_int(parameterlist,cElement_i,"ConvergenceLimit", mySigma%mySegment(iSeg)%ConvergenceDetector%Limit, 250)
     
       call INIP_getvalue_string(parameterlist,cElement_i,"SensorType", mySigma%mySegment(iSeg)%TemperatureSensor%SensorType,"COOR")
       CALL inip_toupper_replace(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)
!        IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'COOR') then
         call INIP_getvalue_string(parameterlist,cElement_i,"TemperatureSensorCoor", sCoorString ," 0d0, 0d0, 0d0")
         read(sCoorString,*) mySigma%mySegment(iSeg)%TemperatureSensor%Coor
         mySigma%mySegment(iSeg)%TemperatureSensor%Coor = dSizeScale*mySigma%mySegment(iSeg)%TemperatureSensor%Coor
         call INIP_getvalue_double(parameterlist,cElement_i,"TemperatureSensorRadius", mySigma%mySegment(iSeg)%TemperatureSensor%Radius,myInf)
         mySigma%mySegment(iSeg)%TemperatureSensor%Radius = dSizeScale*mySigma%mySegment(iSeg)%TemperatureSensor%Radius
!        end if
       IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'OFF'.OR.ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'STL') then
        mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"sensorOFF")
        IF (mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles.gt.0) THEN
         ALLOCATE(mySigma%mySegment(iSeg)%TemperatureSensor%OFFfiles(mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles))
        ELSE
         WRITE(*,*) "STL geometry description files are missing"
         WRITE(*,*) 'sensorOFF'
         bReadError=.TRUE.
         !GOTO 10
        END IF
        do iFile=1,mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles
          call INIP_getvalue_string(parameterlist,cElement_i,"sensorOFF",mySigma%mySegment(iSeg)%TemperatureSensor%OFFfiles(iFile),isubstring=iFile)
        end do
       end if
     

       call INIP_getvalue_string(parameterlist,cElement_i,"Regulation", mySigma%mySegment(iSeg)%Regulation,"SIMPLE")
       CALL inip_toupper_replace(mySigma%mySegment(iSeg)%Regulation)
       
       IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%Regulation)).eq."SIMPLE") THEN
        call INIP_getvalue_double(parameterlist,cElement_i,"TemperatureSensorMinRegValue", mySigma%mySegment(iSeg)%TemperatureSensor%MinRegValue,myInf)
        call INIP_getvalue_double(parameterlist,cElement_i,"TemperatureSensorMaxRegValue", mySigma%mySegment(iSeg)%TemperatureSensor%MaxRegValue,myInf)
        IF (mySigma%mySegment(iSeg)%TemperatureSensor%Radius.eq.myinf.or. &
            mySigma%mySegment(iSeg)%HeatSourceMax.eq.myinf.or.&
            mySigma%mySegment(iSeg)%HeatSourceMin.eq.myinf.or.&
            mySigma%mySegment(iSeg)%TemperatureSensor%MinRegValue.eq.myinf.or.&
            mySigma%mySegment(iSeg)%TemperatureSensor%MaxRegValue.eq.myinf) THEN
            
            WRITE(*,*) "Wrongly defined heating wire parameterlist: ",&
            mySigma%mySegment(iSeg)%TemperatureSensor%Radius,&
            mySigma%mySegment(iSeg)%HeatSourceMax,&
            mySigma%mySegment(iSeg)%HeatSourceMin,&
            mySigma%mySegment(iSeg)%TemperatureSensor%MinRegValue,&
            mySigma%mySegment(iSeg)%TemperatureSensor%MaxRegValue
        END IF
       ELSEIF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%Regulation)).eq."PID") then
        call INIP_getvalue_double(parameterlist,cElement_i,"PID_TemperatureSetValue", mySigma%mySegment(iSeg)%PID_Ctrl%T_set,myInf)
        call INIP_getvalue_double(parameterlist,cElement_i,"PID_P_CNST", mySigma%mySegment(iSeg)%PID_Ctrl%omega_P,myInf)
        call INIP_getvalue_double(parameterlist,cElement_i,"PID_I_CNST", mySigma%mySegment(iSeg)%PID_Ctrl%omega_I,myInf)
        call INIP_getvalue_double(parameterlist,cElement_i,"PID_D_CNST", mySigma%mySegment(iSeg)%PID_Ctrl%omega_D,myInf)
        mySigma%mySegment(iSeg)%PID_Ctrl%SumI = 0d0
        
        IF (mySigma%mySegment(iSeg)%HeatSourceMax.eq.myinf.or.&
            mySigma%mySegment(iSeg)%HeatSourceMin.eq.myinf.or.&
            mySigma%mySegment(iSeg)%PID_Ctrl%T_Set.eq.myinf.or.&
            mySigma%mySegment(iSeg)%PID_Ctrl%Omega_P.eq.myinf.or.&
            mySigma%mySegment(iSeg)%PID_Ctrl%Omega_I.eq.myinf.or.&
            mySigma%mySegment(iSeg)%PID_Ctrl%Omega_D.eq.myinf) THEN
            
            WRITE(*,*) "Wrongly defined heating wire parameterlist: ",&
            mySigma%mySegment(iSeg)%HeatSourceMax,&
            mySigma%mySegment(iSeg)%HeatSourceMin,&
            mySigma%mySegment(iSeg)%PID_Ctrl%T_Set,&
            mySigma%mySegment(iSeg)%PID_Ctrl%omega_P,&
            mySigma%mySegment(iSeg)%PID_Ctrl%omega_I,&
            mySigma%mySegment(iSeg)%PID_Ctrl%omega_D
        END IF
       ELSE
            WRITE(*,*) "Unknown regulation mechanism ... "
            CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

       END IF
     END IF
     
     call INIP_getvalue_string(parameterlist,cElement_i,"Type",cElemType)
     mySigma%mySegment(iSeg)%ART = ' '
     call inip_toupper_replace(cElemType)

!!!=============================================================================================================================
     IF (ADJUSTL(TRIM(cElemType)).eq."STL") THEN
      mySigma%mySegment(iSeg)%ART   = "STL"

      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL geometry description files are missing"
       WRITE(*,*) 'screwOFF'
       bReadError=.TRUE.
       !GOTO 10
      END IF
      do iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
        call INIP_getvalue_string(parameterlist,cElement_i,"screwOFF",mySigma%mySegment(iSeg)%OFFfiles(iFile),isubstring=iFile)
      end do
     END IF
     IF (mySigma%mySegment(iSeg)%ART.eq.' ') THEN
      WRITE(*,*) "not a valid ",iSeg, "-segment"
      WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
      bReadError=.TRUE.
      !GOTO 10
     ENDIF
    END DO

    call INIP_getvalue_int(parameterlist,"E3DMaterialParameters","NoOfMaterials", mySigma%NumberOfMat ,0)
    IF (mySigma%NumberOfMat.ge.1.and.mySigma%NumberOfMat.le.29) THEN
     ALLOCATE (myMaterials(0:mySigma%NumberOfMat-1))
    ELSE
     WRITE(*,*) "not a valid number of materials"
     WRITE(*,*) '"',mySigma%NumberOfMat,'"'
    ENDIF
    
    DO iMat=0,mySigma%NumberOfMat-1
      WRITE(cElement_i,'(A,I0)') 'E3DMaterialParameters/Mat_',iMat
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatConductivity",myMaterials(iMat)%lambda,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatConductivitySlope",myMaterials(iMat)%lambdaSteig,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatCapacity",myMaterials(iMat)%cp,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatCapacitySlope",myMaterials(iMat)%CpSteig,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"Density",myMaterials(iMat)%DensityT0,myInf)
      myMaterials(iMat)%Density = myMaterials(iMat)%DensityT0
      call INIP_getvalue_double(parameterlist,cElement_i,"DensitySlope",myMaterials(iMat)%DensitySteig,myInf)
!     myThermodyn%Alpha     = (1e6*myThermodyn%lambda)/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
!     myThermodyn%Beta      = 1d1 ! !myProcess%Cnst_Lam/(myProcess%Cnst_Dens*myProcess%Cnst_Cp)
!     myThermodyn%Gamma     = 1d0/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
      myMaterials(iMat)%Alpha = (1e5*myMaterials(iMat)%lambda)/((1e0*myMaterials(iMat)%DensityT0)*(myMaterials(iMat)%Cp*1e7))
    END DO
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","AmbientTemperature", myProcess%AmbientTemperature ,myInf)
    if (myProcess%AmbientTemperature.eq.myInf) then
     if (myid.eq.1) WRITE(*,*)  "Ambient temperature is undefined ==>", myProcess%AmbientTemperature
    end if
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","CoolingWaterTemperatureC", myProcess%CoolingWaterTemperature ,55d0)
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","WorkBenchThicknessCM", myProcess%WorkBenchThickness ,5d0)
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","MeltInflowTemperature", myProcess%MeltInflowTemperature ,290d0)

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","ConvergenceEstimator",cConvergenceEstimator,'OFF')
    call inip_toupper_replace(cConvergenceEstimator)
    IF (TRIM(ADJUSTL(cConvergenceEstimator)).eq."ON".or.TRIM(ADJUSTL(cConvergenceEstimator)).eq."YES") THEN
     mySetup%bConvergenceEstimator = .TRUE.
    END IF

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","SensorPositions", mySigma%cSensorPositions,"_INVALID_")
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","SensorRadius", mySigma%SensorRadius ,myInf)
    if (ADJUSTL(TRIM(mySigma%cSensorPositions)).NE."_INVALID_") THEN
     IF (mySigma%SensorRadius.eq.myInf) THEN
      if (myid.eq.1) WRITE(*,*) "SENSORRADIUS in  not defined for the sensorcandidates !!! "
     ELSE
      mySigma%SensorRadius = mySigma%SensorRadius*0.1d0 !! convert MM to cm
     END IF
    END IF
    
    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","HexMesher", mySetup%cMesher,"OFF")
    call inip_toupper_replace(mySetup%cMesher)
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."BOX") THEN
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherELems",cMeshQuality,'0,0,0')
     read(cMeshQuality,*) mySetup%m_nX,mySetup%m_nY,mySetup%m_nZ
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","BoxMesherNumberOfELems",mySetup%nBoxElem,0)
     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","MeshResolution",mySetup%MeshResolution,0)
     
     IF ((mySetup%nBoxElem.le.0).and.&
         (mySetup%m_nX.le.0.and.mySetup%m_nY.le.0.and.mySetup%m_nZ.le.0).and.&
          mySetup%MeshResolution.le.0) THEN
          if (myid.eq.1) WRITE(*,*) 'No rules defined to create a mesh...'
          CALL StopTheProgramFromReader(subnodes,myErrorCode%SIGMA_READER)

     END IF
     
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherUnit",cUnit,'cm')
     call inip_toupper_replace(cUnit)
     IF (.NOT.(TRIM(cUnit).eq.'MM'.OR.TRIM(cUnit).eq.'CM'.OR.TRIM(cUnit).eq.'DM'.OR.TRIM(cUnit).eq.'M')) THEN
       WRITE(*,*) "Unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
       cUnit = 'cm'
     END IF
     if (TRIM(cUnit).eq.'MM') daux = 0.100d0
     if (TRIM(cUnit).eq.'CM') daux = 1.000d0
     if (TRIM(cUnit).eq.'DM') daux = 10.00d0
     if (TRIM(cUnit).eq.'M')  daux = 100.0d0
    
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherX",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(1,:)
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherY",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(2,:)
     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","BoxMesherZ",cMeshQuality,'-1.0,1.0')
     read(cMeshQuality,*) mySetup%m_box(3,:)
     
     !!!!!!!!!! Scaling the geo to its right dimension !!!!
     mySetup%m_box = daux*mySetup%m_box
     
    END IF
    
    
    IF (myid.eq.1.or.subnodes.eq.0) then
    write(*,*) "=========================================================================="
    write(*,*) "mySigma%Type",'=',trim(mySigma%cType)
    write(*,*) "mySigma%Dz_Out",'=',mySigma%Dz_out
    write(*,*) "mySigma%NumberOfSeg",'=',mySigma%NumberOfSeg
    write(*,*) "mySigma%NumberOfMat",'=',mySigma%NumberOfMat
    
    write(*,*) 
    DO iSeg=1,mySigma%NumberOfSeg
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Art=',mySigma%mySegment(iSeg)%ART
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%ObjectType=',mySigma%mySegment(iSeg)%ObjectType
     write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Unit=',mySigma%mySegment(iSeg)%Unit
     write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%MaterialIndex=',mySigma%mySegment(iSeg)%MatInd
     write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%InitialTemperature=',mySigma%mySegment(iSeg)%InitTemp
     IF (TRIM(mySigma%mySegment(iSeg)%ObjectType).eq.'WIRE') THEN
      write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%VolumetricHeatSourceMax=',mySigma%mySegment(iSeg)%HeatSourceMax
      write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%VolumetricHeatSourceMin=',mySigma%mySegment(iSeg)%HeatSourceMin
!       IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'COOR') then
!        write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%SensorType=','COOR'
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%TemperatureSensorCoor=',mySigma%mySegment(iSeg)%TemperatureSensor%Coor
       write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%TemperatureSensorRadius=',mySigma%mySegment(iSeg)%TemperatureSensor%Radius
!       end if
      IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'STL'.or.ADJUSTL(TRIM(mySigma%mySegment(iSeg)%TemperatureSensor%SensorType)).eq.'OFF') then
       write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%SensorType=','OFF file'
       write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")%TemperatureSensor%nOFFfiles=",mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles
       DO iFile=1,mySigma%mySegment(iSeg)%TemperatureSensor%nOFFfiles
        write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%TemperatureSensor%OFFfiles(iFile))),'"'
       END DO
      end if

      write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%ConvergenceCondition=',mySigma%mySegment(iSeg)%ConvergenceDetector%Condition
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,')%ConvergenceLimit=',mySigma%mySegment(iSeg)%ConvergenceDetector%Limit

      write(*,'(A,I0,A,A)') " mySIGMA%Segment(",iSeg,')%Regulation=',ADJUSTL(TRIM(mySigma%mySegment(iSeg)%Regulation))
      IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%Regulation)).eq."SIMPLE") then
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%TemperatureSensorMinRegValue=',mySigma%mySegment(iSeg)%TemperatureSensor%MinRegValue
       write(*,'(A,I0,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%TemperatureSensorMaxRegValue=',mySigma%mySegment(iSeg)%TemperatureSensor%MaxRegValue
      ELSEIF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%Regulation)).eq."PID") then
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%PID_TemperatureSetValue=',mySigma%mySegment(iSeg)%PID_Ctrl%T_Set
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%PID_OmegaP=',mySigma%mySegment(iSeg)%PID_Ctrl%Omega_P
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%PID_OmegaI=',mySigma%mySegment(iSeg)%PID_Ctrl%Omega_I
       write(*,'(A,I0,A,3ES12.4)') " mySIGMA%Segment(",iSeg,')%PID_OmegaD=',mySigma%mySegment(iSeg)%PID_Ctrl%Omega_D
      END IF
     END IF
     
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN
      write(*,'(A,I0,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     write(*,*) 
    END DO    

    write(*,*) 
    DO iMat=0,mySigma%NumberOfMat-1
     write(*,'(A,I2.2,A)') "Material(",iMat,') properties:'
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%HeatConductivity(",iMat,')=',myMaterials(iMat)%lambda
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%HeatConductivitySlope(",iMat,')=',myMaterials(iMat)%lambdaSteig
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%HeatCapacity(",iMat,')=',myMaterials(iMat)%cp
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%HeatCapacitySlope(",iMat,')=',myMaterials(iMat)%cpSteig
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%DensityT0(",iMat,')=',myMaterials(iMat)%densityT0
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%Density(",iMat,')=',myMaterials(iMat)%density
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%DensitySlope(",iMat,')=',myMaterials(iMat)%densitySteig
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%ALPHA(",iMat,')=',myMaterials(iMat)%Alpha
     write(*,'(A)') 
    END DO
    write(*,*) 
    
    write(*,'(A,A)') "mySIGMA%SensorPositionsFile= ",adjustl(trim(mySigma%cSensorPositions))
    if (ADJUSTL(TRIM(mySigma%cSensorPositions)).NE."_INVALID_") THEN
      write(*,'(A,E12.2)') "mySIGMA%SensorRadius= ",mySigma%SensorRadius
    END IF
    
    IF (ADJUSTL(TRIM(mySetup%cMesher)).eq."BOX") THEN
     IF (mySetup%nBoxElem.gt.0) write(*,'(A,I0)') "myMesh%nBoxElem=",mySetup%nBoxElem
     IF (mySetup%m_nX.gt.0.and.mySetup%m_nY.gt.0.and.mySetup%m_nZ.gt.0) write(*,'(A,3I0)') &
     "myMesh%nBoxElemXYZ=",mySetup%m_nX,mySetup%m_nY,mySetup%m_nZ
     IF (mySetup%MeshResolution.gt.0) write(*,'(A,I0)') "myMesh%MeshResolution=",mySetup%MeshResolution
     write(*,'(A,2ES12.4)') "myMesh%Extent_X=",mySetup%m_box(1,:)
     write(*,'(A,2ES12.4)') "myMesh%Extent_Y=",mySetup%m_box(2,:)
     write(*,'(A,2ES12.4)') "myMesh%Extent_Z=",mySetup%m_box(3,:)
    END IF
    
    write(*,'(A,ES12.4)') "myProcess%AmbientTemperature=",myProcess%AmbientTemperature 
!     write(*,'(A,3ES12.4)') "myProcess%TemperatureSensorCoor=",myProcess%TemperatureSensorCoor
!     write(*,'(A,3ES12.4)') "myProcess%TemperatureSensorRadius=",myProcess%TemperatureSensorRadius
    write(*,'(A,ES12.4)') "myProcess%WorkBenchThicknessCM=",myProcess%WorkBenchThickness
    write(*,'(A,ES12.4)') "myProcess%CoolingWaterTemperatureC=",myProcess%CoolingWaterTemperature
    write(*,'(A,ES12.4)') "myProcess%MeltInflowTemperature=",myProcess%MeltInflowTemperature
    write(*,*) 
    write(*,'(A,L)') "mySetup%ConvergenceEstimator=",mySetup%bConvergenceEstimator

!     write(*,'(A,ES12.4)') "myProcess%ConductiveLambda=",myProcess%ConductiveLambda 
!     write(*,'(A,ES12.4)') "myProcess%ConductiveGradient=", myProcess%ConductiveGradient
    write(*,*) "=========================================================================="
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
    10  CONTINUE

!    ! Make some output on the terminal
!     call inip_info(parameterlist)

!     ! Write it into a file
!      call inip_dumpToFile(parameterlist,"test_dump_parlst.dat",INIP_REPLACE)

    ! Clean up the parameterlist
    
    call inip_done(parameterlist)

    end Subroutine ReadEWIKONfile
    
    SUBROUTINE StopTheProgramFromReader(n,ie)
    
     USE PP3D_MPI, ONLY:subnodes,MPI_COMM_WORLD
     USE var_QuadScalar

     implicit none
     integer n,ie,ierr
    
     if (n.eq.1) then
      error stop 40 !! cannot get an argument here ....
     else
      call MPI_Abort(MPI_COMM_WORLD, ie, ierr)
     end if
     
    END SUBROUTINE StopTheProgramFromReader
