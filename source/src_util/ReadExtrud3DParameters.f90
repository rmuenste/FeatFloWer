    Subroutine ReadS3Dfile(cE3Dfile)
    use iniparser
    use Sigma_User

    use, intrinsic :: ieee_arithmetic

    implicit none

    character(len=*), intent(in) :: cE3Dfile
    logical :: bReadError=.FALSE.
    integer :: i,iSeg,iFile,iaux,iInflow,iInflowErr

    real*8 :: myPI = dATAN(1d0)*4d0
    character(len=INIP_STRLEN) cCut,cElement_i,cElemType,cKindOfConveying,cTemperature,cPressureFBM
    character(len=INIP_STRLEN) cBCtype,cInflow_i,cCenter,cNormal

    character(len=INIP_STRLEN) cProcessType,cRotation,cRheology,cMeshQuality,cKTP,cUnit,cOFF_Files
    
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    real*8 :: myInf,dSizeScale

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

    
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Machine","Type",mySigma%cType,'SSE')
    call inip_toupper_replace(mySigma%cType)
    IF (.NOT.(ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."DIE".OR.ADJUSTL(TRIM(mySigma%cType)).EQ."NETZSCH")) THEN
     WRITE(*,*) "not a valid Extruder type:", ADJUSTL(TRIM(mySigma%cType))
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
    IF (ADJUSTL(TRIM(mySigma%cZwickel)).eq."ROUND".OR.ADJUSTL(TRIM(mySigma%cZwickel)).EQ."CURVED")    mySigma%W = myInf

    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","InnerDiameter", mySigma%Dz_in ,mySigma%Dz_Out/dSizeScale)
    mySigma%Dz_in = dSizeScale*mySigma%Dz_in
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Machine","BarrelLength", mySigma%L ,myInf)
    mySigma%L = dSizeScale*mySigma%L

    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfElements", mySigma%NumberOfSeg ,0)
    call INIP_getvalue_int(parameterlist,"E3DGeometryData/Machine","NoOfFlights", mySigma%GANGZAHL ,0)
    
    IF (mySigma%NumberOfSeg.ge.1.and.mySigma%NumberOfSeg.le.9) THEN
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
     WRITE(cElement_i,'(A,I1.1)') 'E3DGeometryData/Machine/Element_',iSeg

     call INIP_getvalue_string(parameterlist,cElement_i,"ObjectType",mySigma%mySegment(iSeg)%ObjectType)
     call inip_toupper_replace(mySigma%mySegment(iSeg)%ObjectType)
     IF (.NOT.(mySigma%mySegment(iSeg)%ObjectType.eq.'SCREW'.OR.&
               mySigma%mySegment(iSeg)%ObjectType.eq.'DIE'.OR.mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE')) THEN
       WRITE(*,*) "STL object type is invalid. Only screw, die, or obstacle types are allowed"
     END IF

     call INIP_getvalue_string(parameterlist,cElement_i,"Unit",mySigma%mySegment(iSeg)%Unit,'MM')
     call inip_toupper_replace(mySigma%mySegment(iSeg)%Unit)
     IF (.NOT.(mySigma%mySegment(iSeg)%Unit.eq.'MM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'CM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'DM'.OR.mySigma%mySegment(iSeg)%Unit.eq.'M')) THEN
       WRITE(*,*) "STL unit type is invalid. Only MM, CM, DM or 'M' units are allowed",mySigma%mySegment(iSeg)%Unit
       bReadError=.TRUE.
!       mySigma%mySegment(iSeg)%Unit = 'MM'
     END IF
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'MM') dSizeScale = 0.100d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'CM') dSizeScale = 1.000d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'DM') dSizeScale = 10.00d0
     if (TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)).eq.'M')  dSizeScale = 100.0d0
      
!     WRITE(*,*) "'",TRIM(ADJUSTL(mySigma%mySegment(iSeg)%Unit)),"'", dSizeScale

     call INIP_getvalue_string(parameterlist,cElement_i,"Type",cElemType)
     mySigma%mySegment(iSeg)%ART = ' '
     call inip_toupper_replace(cElemType)

!!!==============================================     FOERD    =================================================================
!!!=============================================================================================================================
     IF (ADJUSTL(TRIM(cElemType)).eq."FOERD".or.ADJUSTL(TRIM(cElemType)).eq."THREADED") THEN
      mySigma%mySegment(iSeg)%ART   = "FOERD"
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Lead", mySigma%mySegment(iSeg)%t,myInf)
      mySigma%mySegment(iSeg)%t = dSizeScale*mySigma%mySegment(iSeg)%t


      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s

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
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds
     
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out-mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))

      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dSizeScale*mySigma%mySegment(iSeg)%D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      
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
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds

      call INIP_getvalue_double(parameterlist,cElement_i,"Eccentricity", mySigma%mySegment(iSeg)%excentre, myInf)
      mySigma%mySegment(iSeg)%excentre = dSizeScale*mySigma%mySegment(iSeg)%excentre
     
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out-2d0*mySigma%mySegment(iSeg)%excentre-mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s - mySigma%mySegment(iSeg)%excentre))
    
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dSizeScale*mySigma%mySegment(iSeg)%D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      
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
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min

      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds

      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s
      
      mySigma%mySegment(iSeg)%delta=(mySigma%Dz_Out - mySigma%mySegment(iSeg)%Ds)/2d0
      mySigma%Dz_In = min(mySigma%Dz_In,2D0*(mySigma%a - 0.5d0*mySigma%mySegment(iSeg)%Ds - mySigma%mySegment(iSeg)%s))
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%D,myInf)
      mySigma%mySegment(iSeg)%D = dSizeScale*mySigma%mySegment(iSeg)%D

      call INIP_getvalue_double(parameterlist,cElement_i,"StaggeringAngle", mySigma%mySegment(iSeg)%alpha,myInf)
      
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
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Lead", mySigma%mySegment(iSeg)%t,myInf)
      mySigma%mySegment(iSeg)%t = dSizeScale*mySigma%mySegment(iSeg)%t
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s
      
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
      mySigma%mySegment(iSeg)%SecProf_D = dSizeScale*mySigma%mySegment(iSeg)%SecProf_D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesWidth",mySigma%mySegment(iSeg)%SecProf_W  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_W = dSizeScale*mySigma%mySegment(iSeg)%SecProf_W
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesLead",mySigma%mySegment(iSeg)%SecProf_L  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_L = dSizeScale*mySigma%mySegment(iSeg)%SecProf_L
      
      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + mySigma%mySegment(iSeg)%L
     
!!!==============================================      ZME     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."ZME".or.ADJUSTL(TRIM(cElemType)).eq."TOOTHMIXING") THEN
      mySigma%mySegment(iSeg)%ART   = "ZME"
      
      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,myInf)
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_int(parameterlist,cElement_i,"NoOfRows",mySigma%mySegment(iSeg)%ZME_N ,-1)
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscWidth", mySigma%mySegment(iSeg)%ZME_DiscThick,myInf)
      mySigma%mySegment(iSeg)%ZME_DiscThick = dSizeScale*mySigma%mySegment(iSeg)%ZME_DiscThick
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscDiscGap", mySigma%mySegment(iSeg)%ZME_gap_SS,myInf)
      mySigma%mySegment(iSeg)%ZME_gap_SS = dSizeScale*mySigma%mySegment(iSeg)%ZME_gap_SS
      
      call INIP_getvalue_double(parameterlist,cElement_i,"DiscShellGap", mySigma%mySegment(iSeg)%ZME_gap_SG,myInf)
      mySigma%mySegment(iSeg)%ZME_gap_SG = dSizeScale*mySigma%mySegment(iSeg)%ZME_gap_SG
!       call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = mySigma%Dz_Out - 2d0*mySigma%mySegment(iSeg)%ZME_gap_SG


      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,myInf)
      mySigma%mySegment(iSeg)%Dss = dSizeScale*mySigma%mySegment(iSeg)%Dss
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)

      call INIP_getvalue_int(parameterlist,cElement_i,"NoOfTeeth",mySigma%mySegment(iSeg)%SecProf_N ,-1)
      call INIP_getvalue_int(parameterlist,cElement_i,"KindOfGrooves",mySigma%mySegment(iSeg)%SecProf_I  ,1)
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesDepth",mySigma%mySegment(iSeg)%SecProf_D  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_D = dSizeScale*mySigma%mySegment(iSeg)%SecProf_D
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesWidth",mySigma%mySegment(iSeg)%SecProf_W  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_W = dSizeScale*mySigma%mySegment(iSeg)%SecProf_W
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GroovesLead",mySigma%mySegment(iSeg)%SecProf_L  ,myInf)
      mySigma%mySegment(iSeg)%SecProf_L = dSizeScale*mySigma%mySegment(iSeg)%SecProf_L

      mySigma%mySegment(iSeg)%Max= mySigma%mySegment(iSeg)%Min + &
      2d0*mySigma%mySegment(iSeg)%ZME_N * (mySigma%mySegment(iSeg)%ZME_DiscThick + mySigma%mySegment(iSeg)%ZME_gap_SS)
      mySigma%mySegment(iSeg)%L = mySigma%mySegment(iSeg)%Max - mySigma%mySegment(iSeg)%Min
     
!!!==============================================      STL     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL".OR.ADJUSTL(TRIM(cElemType)).eq."OFF") THEN
      mySigma%mySegment(iSeg)%ART   = "STL"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
!       write(*,*) 'dSizeScale: ',dSizeScale
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dSizeScale)
      mySigma%mySegment(iSeg)%Dss = dSizeScale*mySigma%mySegment(iSeg)%Dss
      
      call INIP_getvalue_double(parameterlist,cElement_i,"Diameter", mySigma%mySegment(iSeg)%Ds,myInf)
      mySigma%mySegment(iSeg)%Ds = dSizeScale*mySigma%mySegment(iSeg)%Ds
      
      call INIP_getvalue_double(parameterlist,cElement_i,"GapScrewScrew", mySigma%mySegment(iSeg)%s,myInf)
      mySigma%mySegment(iSeg)%s = dSizeScale*mySigma%mySegment(iSeg)%s
      
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
        WRITE(*,*) "STL geometry dscription files are missing"
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
     
!!!==============================================      STL_R     =================================================================
!!!=============================================================================================================================
     ELSE IF (ADJUSTL(TRIM(cElemType)).eq."STL_R") THEN
      mySigma%mySegment(iSeg)%ART   = "STL_R"

      call INIP_getvalue_double(parameterlist,cElement_i,"StartPosition",mySigma%mySegment(iSeg)%Min,0d0)
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dSizeScale)
      mySigma%mySegment(iSeg)%Dss = dSizeScale*mySigma%mySegment(iSeg)%Dss
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
      
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL_R geometry dscription files are missing"
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
      mySigma%mySegment(iSeg)%Min = dSizeScale*mySigma%mySegment(iSeg)%Min
      
      call INIP_getvalue_double(parameterlist,cElement_i,"ElementLength", mySigma%mySegment(iSeg)%L ,myInf)
      mySigma%mySegment(iSeg)%L = dSizeScale*mySigma%mySegment(iSeg)%L
      
      call INIP_getvalue_double(parameterlist,cElement_i,"InnerDiameter", mySigma%mySegment(iSeg)%Dss,mySigma%Dz_In/dSizeScale)
      mySigma%mySegment(iSeg)%Dss = dSizeScale*mySigma%mySegment(iSeg)%Dss
      
      mySigma%Dz_In = min(mySigma%Dz_In,mySigma%mySegment(iSeg)%Dss)
!       mySigma%mySegment(iSeg)%Ds = mySigma%Dz_In
      
      mySigma%mySegment(iSeg)%nOFFfiles = INIP_querysubstrings(parameterlist,cElement_i,"screwOFF")
      IF (mySigma%mySegment(iSeg)%nOFFfiles.gt.0) THEN
       ALLOCATE(mySigma%mySegment(iSeg)%OFFfiles(mySigma%mySegment(iSeg)%nOFFfiles))
      ELSE
       WRITE(*,*) "STL_L geometry dscription files are missing"
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
    IF (TRIM(adjustl(sExtract_Dim)).eq.'PA') myProcess%dPress  = 1d1*dExtract_Val
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
    
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewSpeed", myProcess%umdr,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters","MaterialTemperature",myProcess%T0,myInf)

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","ScrewTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","ScrewTemperature",myProcess%Ti,myInf)
    ELSE
     myProcess%Ti=myInf
    END IF
    cTemperature=" "
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters","BarrelTemperatureAdiabatic", cTemperature,"YES")
    call inip_toupper_replace(cTemperature)
    IF (ADJUSTL(TRIM(cTemperature)).EQ."NO") THEN
     call INIP_getvalue_double(parameterlist,"E3DProcessParameters","BarrelTemperature",myProcess%Ta,myInf)
    ELSE
     myProcess%Ta=myInf
    END IF
       
    call INIP_getvalue_int(parameterlist,"E3DProcessParameters",   "nOfInflows"      ,myProcess%nOfInflows,0)
    ALLOCATE(myProcess%myInflow(myProcess%nOfInflows))
    DO iInflow=1,myProcess%nOfInflows
     if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'E3DProcessParameters/Inflow_',iInflow
     if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'E3DProcessParameters/Inflow_',iInflow
     if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'E3DProcessParameters/Inflow_',iInflow
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
     if (myProcess%myInflow(iInflow)%iBCtype.eq.0) write(*,*) 'UNDEFINED Inflow type!!'
     call INIP_getvalue_double(parameterlist,cInflow_i,"massflowrate",myProcess%myInflow(iInflow)%massflowrate,myInf)
     if (myProcess%myInflow(iInflow)%massflowrate.eq.myInf) write(*,*) 'UNDEFINED massflowrate through Inflow',iInflow,' !!'
     call INIP_getvalue_double(parameterlist,cInflow_i,"density",myProcess%myInflow(iInflow)%density,myInf)
     if (myProcess%myInflow(iInflow)%density.eq.myInf) write(*,*) 'UNDEFINED density forInflow',iInflow,' !!'
     call INIP_getvalue_double(parameterlist,cInflow_i,"innerradius",myProcess%myInflow(iInflow)%innerradius,myInf)
     if (myProcess%myInflow(iInflow)%innerradius.eq.myInf) write(*,*) 'UNDEFINED inner radius for Inflow',iInflow,' !!'
     call INIP_getvalue_double(parameterlist,cInflow_i,"outerradius",myProcess%myInflow(iInflow)%outerradius,myInf)
     if (myProcess%myInflow(iInflow)%outerradius.eq.myInf) write(*,*) 'UNDEFINED outer radius for Inflow',iInflow,' !!'
     call INIP_getvalue_string(parameterlist,cInflow_i,"center",cCenter,'unknown')
     call INIP_getvalue_string(parameterlist,cInflow_i,"normal",cNormal,'unknown')
     read(cCenter,*,err=55) myProcess%myInflow(iInflow)%Center
     if (iInflowErr.ne.0) write(*,*) 'WRONGLY DEFINED center for Inflow',iInflow,' !!'
     read(cNormal,*,err=56) myProcess%myInflow(iInflow)%Normal
     if (iInflowErr.ne.0) write(*,*) 'WRONGLY DEFINED normal for Inflow',iInflow,' !!'
     GOTO 57     
55   write(*,*) 'WRONGLY DEFINED center for Inflow',iInflow,' !!'
     GOTO 57
56   write(*,*) 'WRONGLY DEFINED normal for Inflow',iInflow,' !!'  
     GOTO 57
57   CONTINUE
     
    END DO
    
    myRheology%Equation = 0
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/RheologicalData","CalcVisco", cRheology,'NoRheology')
    call inip_toupper_replace(cRheology)
    IF (ADJUSTL(TRIM(cRheology)).eq."BINGHAM") THEN
      myRheology%Equation = 6
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Bingham","ZeroViscosity",myRheology%A,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Bingham","Regularization",myRheology%B,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Bingham","YieldStress",myRheology%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."CARREAU") THEN
      myRheology%Equation = 1
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","ZeroViscosity",myRheology%A,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","RecipVelocity",myRheology%B,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Carreau","Exponent",myRheology%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POWERLAW".OR.ADJUSTL(TRIM(cRheology)).eq."POTENZ".OR.ADJUSTL(TRIM(cRheology)).eq."POWER") THEN
      myRheology%Equation = 2
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Powerlaw","Consistence", myRheology%K,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Powerlaw","Exponent",myRheology%n,myInf)
      if (myRheology%K.eq.myInf) then
       call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Power","Consistence", myRheology%K,myInf)
      end if
      if (myRheology%n.eq.myInf) then
       call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Power","Exponent",myRheology%n,myInf)
      end if
      
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."POLYFLOW") THEN
      myRheology%Equation = 3
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_A",myRheology%A,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_B",myRheology%B,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Polyflow","Polyflow_C",myRheology%C,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."ELLIS") THEN
      myRheology%Equation = 4
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Ellis","ZeroViscosity",myRheology%A,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Ellis","Gamma0",myRheology%B,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/Ellis","Exponent",myRheology%C,myInf)
    END IF

    IF (myRheology%Equation.eq.0) THEN
     WRITE(*,*) "no valid rheology is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
     !GOTO 10
    END IF
   
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material","LimitViscoMin",myRheology%ViscoMin,1d0)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material","LimitViscoMax",myRheology%ViscoMax,1d5)

    myRheology%AtFunc = 0
    cRheology = ' '
    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/RheologicalData","CalcTemp", cRheology,'ISOTHERM')
    call inip_toupper_replace(cRheology)
    
    IF (ADJUSTL(TRIM(cRheology)).eq."ISOTHERM") THEN
      myRheology%AtFunc = 1
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."C1C2") THEN
      myRheology%AtFunc = 2
      myRheology%Ts = 165d0
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/C1C2","C1",myRheology%C1,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/C1C2","C2",myRheology%C2,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."TBTS") THEN
      myRheology%AtFunc = 3
      myRheology%C1 = 8.86d0
      myRheology%C2 = 101.6d0
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/TBTS","ReferenceTemperature",myRheology%Tb,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/TBTS","StandardTemperature",myRheology%Ts,myInf)
    END IF
    IF (ADJUSTL(TRIM(cRheology)).eq."ETB") THEN
      myRheology%AtFunc = 4
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/ETB","ActivatingEnergy",myRheology%E,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/RheologicalData/ETB","ReferenceTemperature",myRheology%TB,myInf)
    END IF
    IF (myRheology%AtFunc.eq.0) THEN
     WRITE(*,*) "no temperature correction is defined"
     WRITE(*,*) '"',TRIM(cRheology),'"'
     bReadError=.TRUE.
     !GOTO 10
    END IF

    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatConductivity",myThermodyn%lambda,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatConductivitySlope",myThermodyn%lambdaSteig,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatCapacity",myThermodyn%cp,myInf)
    call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData","HeatCapacitySlope",myThermodyn%CpSteig,myInf)

    call INIP_getvalue_string(parameterlist,"E3DProcessParameters/Material/ThermoData","DensityModel", myThermodyn%DensityModel,'no')
    call inip_toupper_replace(myThermodyn%DensityModel)
    myThermodyn%density=myInf
    IF (ADJUSTL(TRIM(myThermodyn%DensityModel)).eq."DENSITY") THEN
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/Density","Density",myThermodyn%Density,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/Density","DensitySlope",myThermodyn%DensitySteig,myInf)
    END IF
    IF (ADJUSTL(TRIM(myThermodyn%DensityModel)).eq."SPECVOLUME") THEN
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/SpecVolume","SpecVolume",myThermodyn%Density,myInf)
      call INIP_getvalue_double(parameterlist,"E3DProcessParameters/Material/ThermoData/SpecVolume","SpecVolumeSlope",myThermodyn%DensitySteig,myInf)
      myThermodyn%Density = 1d0/myThermodyn%Density
    END IF
    IF (myThermodyn%density.eq.myinf) THEN
     WRITE(*,*) "density is not defined"
     WRITE(*,*) '"',TRIM(myThermodyn%DensityModel),'"'
     bReadError=.TRUE.
     !GOTO 10
    END IF

    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output",   "nOf1DLayers"      ,myOutput%nOf1DLayers,16)
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings/Output",   "nOfHistogramBins" ,myOutput%nOfHistogramBins,16)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramShearMax",myOutput%HistogramShearMax,1d5)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramShearMin",myOutput%HistogramShearMin,1d-2)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramViscoMax",myOutput%HistogramViscoMax,1d6)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","HistogramViscoMin",myOutput%HistogramViscoMin,1d0)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings/Output","CutDtata_1D",myOutput%CutDtata_1D,0.04d0)
    
    cKTP=' '
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
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","CharacteristicShearRate",mySetup%CharacteristicShearRate,1d0)
    
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","activeFBM_Z_Position",activeFBM_Z_Position,myInf)

    cKTP=' '
    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","KTPRelease",cKTP,"YES")
    call inip_toupper_replace(cKTP)
    IF (ADJUSTL(TRIM(cKTP)).eq."NO") THEN
     bKTPRelease = .FALSE.
    END IF

    call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","PressureFBM", cPressureFBM,"OFF")
    IF (ADJUSTL(TRIM(cPressureFBM)).eq."ON".OR.ADJUSTL(TRIM(cPressureFBM)).eq."YES") THEN
     mySetup%bPressureFBM = .true.
    ENDIF

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
    
!     cMeshQuality=' '
!     call INIP_getvalue_string(parameterlist,"E3DSimulationSettings","SendEmail",cMeshQuality,"YES")
!     call inip_toupper_replace(cMeshQuality)
!     IF (ADJUSTL(TRIM(cMeshQuality)).eq."NO") THEN
!      mySetup%bSendEmail = .FALSE.
!     END IF
!     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","Periodicity",myProcess%Periodicity,1)
!     call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","nSolutions",mySetup%nSolutions,1)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","dAlpha",myProcess%dAlpha,10d0)
    call INIP_getvalue_double(parameterlist,"E3DSimulationSettings","Angle",myProcess%Angle,myInf)
    call INIP_getvalue_int(parameterlist,"E3DSimulationSettings","Phase",myProcess%Phase,-1)
    
    IF (myid.eq.1.or.subnodes.eq.0) then
    write(*,*) "=========================================================================="
    write(*,*) "mySigma%Type",'=',trim(mySigma%cType)
    write(*,*) "mySigma%Zwickel",'=',trim(mySigma%cZwickel)
    IF (mySigma%Dzz.ne.myinf) write(*,*) "mySigma%Dzz",'=',mySigma%Dzz
    IF (mySigma%W.ne.myinf)   write(*,*) "mySigma%W",'=',mySigma%W
    write(*,*) "mySigma%Dz_Out",'=',mySigma%Dz_out
    write(*,*) "mySigma%Dz_In",'=',mySigma%Dz_In
    write(*,*) "mySigma%L",'=',mySigma%L
    write(*,*) "mySigma%GANGZAHL",'=',mySigma%GANGZAHL
    if (myProcess%iInd.eq.-1) write(*,*) "mySigma%RotationType",'=','CounterRotating'
    if (myProcess%iInd.eq.+1) write(*,*) "mySigma%RotationType",'=','CoRotating'
    IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
     
     write(*,*) "mySigma%RotationAxis",'=',mySigma%RotationAxis
     
     IF (ADJUSTL(TRIM(mySigma%RotationAxis)).EQ."PARALLEL") THEN
      write(*,*) "mySigma%a",'=',mySigma%a
     ELSE
      write(*,*) "mySigma%RotAxisAngle",'=',mySigma%RotAxisAngle
      write(*,*) "mySigma%RotAxisCenter",'=',mySigma%RotAxisCenter
     END IF
    END IF
    
    write(*,*) "mySigma%NumberOfSeg",'=',mySigma%NumberOfSeg
    
    write(*,*) 
    DO iSeg=1,mySigma%NumberOfSeg
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Art=',mySigma%mySegment(iSeg)%ART
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%ObjectType=',mySigma%mySegment(iSeg)%ObjectType
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Unit=',mySigma%mySegment(iSeg)%Unit
     write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Min=',mySigma%mySegment(iSeg)%Min
     write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Max=',mySigma%mySegment(iSeg)%Max
     write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%L=',mySigma%mySegment(iSeg)%L
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."FOERD") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%t=',abs(mySigma%mySegment(iSeg)%t)
      IF (mySigma%mySegment(iSeg)%t.lt.0d0) THEN
       write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','RECONVEYING'
      END IF
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."KNET") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."EKNET") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%E=',mySigma%mySegment(iSeg)%excentre
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."SKNET") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%alpha=',abs(mySigma%mySegment(iSeg)%alpha)
      IF (mySigma%mySegment(iSeg)%alpha.lt.0d0) THEN
       write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','REKNEADING'
      END IF
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta

      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%N=',mySigma%mySegment(iSeg)%N
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%D=',mySigma%mySegment(iSeg)%D
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."SME") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%t=',abs(mySigma%mySegment(iSeg)%t)
      IF (mySigma%mySegment(iSeg)%t.lt.0d0) THEN
       write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Kind=','RECONVEYING'
      END IF
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta

      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_N=',mySigma%mySegment(iSeg)%SecProf_N
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_I=',mySigma%mySegment(iSeg)%SecProf_I
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_D=',mySigma%mySegment(iSeg)%SecProf_d
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_W=',mySigma%mySegment(iSeg)%SecProf_w
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_L=',mySigma%mySegment(iSeg)%SecProf_l
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."ZME") THEN
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%ZME_N=',mySigma%mySegment(iSeg)%ZME_N
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_DiscThick=',mySigma%mySegment(iSeg)%ZME_DiscThick
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_gap_SG=',mySigma%mySegment(iSeg)%ZME_gap_SG
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%ZME_gap_SS=',mySigma%mySegment(iSeg)%ZME_gap_SS
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_N=',mySigma%mySegment(iSeg)%SecProf_N
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%SecProf_I=',mySigma%mySegment(iSeg)%SecProf_I
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_D=',mySigma%mySegment(iSeg)%SecProf_d
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_W=',mySigma%mySegment(iSeg)%SecProf_w
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%SecProf_L=',mySigma%mySegment(iSeg)%SecProf_l
     END IF
     
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Ds=',mySigma%mySegment(iSeg)%Ds
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%s=',mySigma%mySegment(iSeg)%s
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%delta=',mySigma%mySegment(iSeg)%delta
      
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_R") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
      DO iFile=1,mySigma%mySegment(iSeg)%nOFFfiles
       write(*,*) '"',adjustl(trim(mySigma%mySegment(iSeg)%OFFfiles(iFile))),'"'
      END DO
     END IF
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL_L") THEN
      write(*,'(A,I1.1,A,f13.3)') " mySIGMA%Segment(",iSeg,')%Dss=',mySigma%mySegment(iSeg)%Dss
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
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
    write(*,*) "myProcess%nOfInflows",'=',myProcess%nOfInflows
    DO iInflow=1,myProcess%nOfInflows
     if (iInflow.gt.0.and.iInflow.le.9) WRITE(cInflow_i,'(A,I1.1)') 'Inflow_',iInflow
     if (iInflow.gt.9.and.iInflow.le.99) WRITE(cInflow_i,'(A,I2.2)') 'Inflow_',iInflow
     if (iInflow.gt.99.and.iInflow.le.999) WRITE(cInflow_i,'(A,I3.3)') 'Inflow_',iInflow
     
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Type','=',myProcess%myInflow(iInflow)%iBCtype
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Massflowrate','=',myProcess%myInflow(iInflow)%massflowrate
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_Density','=',myProcess%myInflow(iInflow)%density
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_InnerRadius','=',myProcess%myInflow(iInflow)%InnerRadius
     write(*,*) "myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_OuterRadius','=',myProcess%myInflow(iInflow)%OuterRadius
     write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_center'//'=',myProcess%myInflow(iInflow)%center
     write(*,'(A,3ES12.4)') " myProcess%"//ADJUSTL(TRIM(cInflow_i))//'_normal'//'=',myProcess%myInflow(iInflow)%normal
    END DO

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
    IF (myRheology%Equation.eq.4) THEN
     write(*,*) "myRheology%model",'=','Ellis'
     write(*,*) "myRheology%A",'=',myRheology%A
     write(*,*) "myRheology%B",'=',myRheology%B
     write(*,*) "myRheology%C",'=',myRheology%C
    END IF
    IF (myRheology%Equation.eq.6) THEN
     write(*,*) "myRheology%model",'=','Bingham'
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
    IF (myRheology%AtFunc.eq.4) THEN
     write(*,*) "myRheology%TempModel",'=','ETB'
     write(*,*) "myRheology%E",'=',myRheology%E
     write(*,*) "myRheology%Tb",'=',myRheology%Tb
    END IF
    write(*,*)
    write(*,*) "myThermodyn%DensityModel",'=',TRIM(ADJUSTL(myThermodyn%DensityModel))
    write(*,*) "myThermodyn%HeatConductivity",'=',myThermodyn%lambda
    write(*,*) "myThermodyn%HeatConductivitySlope",'=',myThermodyn%lambdaSteig
    write(*,*) "myThermodyn%HeatCapacity",'=',myThermodyn%cp
    write(*,*) "myThermodyn%HeatCapacitySlope",'=',myThermodyn%cpSteig
    write(*,*) "myThermodyn%Density",'=',myThermodyn%density
    write(*,*) "myThermodyn%DensitySlope",'=',myThermodyn%densitySteig
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

    write(*,*) "mySetup%PressureFBM = ",mySetup%bPressureFBM
    write(*,*) "mySetup%AutomaticTimeStepControl = ",mySetup%bAutomaticTimeStepControl
    write(*,*) "mySetup%CharacteristicShearRate = ",mySetup%CharacteristicShearRate
    write(*,*) "activeFBM_Z_Position = ",activeFBM_Z_Position   
    write(*,*) "TimeStepEnlargmentFactor = ",dTimeStepEnlargmentFactor   

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
    write(*,*)
    write(*,*) "myProcess%dAlpha",'=',myProcess%dAlpha
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
      stop
    end if

    CONTAINS

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
     dExtract_Val = myInf
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
    DO i=i1,n
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
    STOP
    
    
    END SUBROUTINE ExtractNomOfCharFromString
    
    
    SUBROUTINE ExtractCharArrayFromString(sL,sS)
    IMPLICIT NONE
    CHARACTER*(*) :: sS(:),sL
    CHARACTER*(256) saux
    INTEGER i,n,i1,i2,i3,j,nS
    
    n = len(sL)
    
    j=0
    i2 = 1
    DO i=i1,n
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
    use Sigma_User

    use, intrinsic :: ieee_arithmetic

    implicit none

    character(len=*), intent(in) :: cE3Dfile
    logical :: bReadError=.FALSE.
    integer :: i,iSeg,iFile,iaux

    real*8 :: myPI = dATAN(1d0)*4d0
    character(len=INIP_STRLEN) cCut,cElement_i,cElemType,cKindOfConveying,cTemperature
    character(len=INIP_STRLEN) cProcessType,cRotation,cRheology,CDensity,cMeshQuality,cKTP,cUnit,sCoorString
    
    integer :: unitProtfile = -1 ! I guess you use mfile here
    integer :: unitTerminal = 6 ! I guess you use mterm here

    type(t_parlist) :: parameterlist

    real*8 :: myInf,dSizeScale
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
      WRITE(*,*) "STL unit type is invalid. Only MM, CM, DM or 'M' units are allowed ",TRIM(cUnit)
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

    IF (mySigma%NumberOfSeg.ge.1.and.mySigma%NumberOfSeg.le.29) THEN
     ALLOCATE (mySigma%mySegment(mySigma%NumberOfSeg))
    ELSE
     WRITE(*,*) "not a valid number of segments"
     WRITE(*,*) '"',mySigma%NumberOfSeg,'"'
    ENDIF

    DO iSeg=1,mySigma%NumberOfSeg
     WRITE(cElement_i,'(A,I1.1)') 'E3DGeometryData/Machine/Element_',iSeg

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
     if (.NOT.(mySigma%mySegment(iSeg)%TemperatureBC.eq.'CONSTANT'.or.mySigma%mySegment(iSeg)%TemperatureBC.eq.'NO')) THEN
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
       WRITE(*,*) "STL geometry dscription files are missing"
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
      WRITE(cElement_i,'(A,I1.1)') 'E3DMaterialParameters/Mat_',iMat
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatConductivity",myMaterials(iMat)%lambda,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatConductivitySlope",myMaterials(iMat)%lambdaSteig,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatCapacity",myMaterials(iMat)%cp,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"HeatCapacitySlope",myMaterials(iMat)%CpSteig,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"Density",myMaterials(iMat)%Density,myInf)
      call INIP_getvalue_double(parameterlist,cElement_i,"DensitySlope",myMaterials(iMat)%DensitySteig,myInf)
!     myThermodyn%Alpha     = (1e6*myThermodyn%lambda)/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
!     myThermodyn%Beta      = 1d1 ! !myProcess%Cnst_Lam/(myProcess%Cnst_Dens*myProcess%Cnst_Cp)
!     myThermodyn%Gamma     = 1d0/((1e-3*myThermodyn%Density)*(myThermodyn%Cp*1e9))
      myMaterials(iMat)%Alpha = (1e5*myMaterials(iMat)%lambda)/((1e0*myMaterials(iMat)%Density)*(myMaterials(iMat)%Cp*1e7))
    END DO
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","AmbientTemperature", myProcess%AmbientTemperature ,myInf)
    if (myProcess%AmbientTemperature.eq.myInf) then
     if (myid.eq.1) WRITE(*,*)  "Ambient temperature is undefined ==>", myProcess%AmbientTemperature
    end if
    call INIP_getvalue_string(parameterlist,"E3DGeometryData/Process","TemperatureSensorCoor", sCoorString ," 0d0, 0d0, 0d0")
    read(sCoorString,*) myProcess%TemperatureSensorCoor
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","TemperatureSensorRadius", myProcess%TemperatureSensorRadius,myInf)
    if (myProcess%TemperatureSensorRadius.eq.myInf) then
     myProcess%TemperatureSensorRadius = 0d0
     if (myid.eq.1) WRITE(*,*) "Temperature sensor is undefined, Radius=", myProcess%TemperatureSensorRadius
    end if
    
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","HeatTransferCoeff", myProcess%HeatTransferCoeff ,0d0)
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","ConductiveLambda", myProcess%ConductiveLambda ,0d0)
    call INIP_getvalue_double(parameterlist,"E3DGeometryData/Process","ConductiveGradient", myProcess%ConductiveGradient ,0d0)

    IF (myid.eq.1.or.subnodes.eq.0) then
    write(*,*) "=========================================================================="
    write(*,*) "mySigma%Type",'=',trim(mySigma%cType)
    write(*,*) "mySigma%Dz_Out",'=',mySigma%Dz_out
    write(*,*) "mySigma%NumberOfSeg",'=',mySigma%NumberOfSeg
    write(*,*) "mySigma%NumberOfMat",'=',mySigma%NumberOfMat
    
    write(*,*) 
    DO iSeg=1,mySigma%NumberOfSeg
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Art=',mySigma%mySegment(iSeg)%ART
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%ObjectType=',mySigma%mySegment(iSeg)%ObjectType
     write(*,'(A,I1.1,A,A)') " mySIGMA%Segment(",iSeg,')%Unit=',mySigma%mySegment(iSeg)%Unit
     write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,')%MaterialIndex=',mySigma%mySegment(iSeg)%MatInd
     write(*,'(A,I1.1,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%InitialTemperature=',mySigma%mySegment(iSeg)%InitTemp
     write(*,'(A,I1.1,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%VolumetricHeatSourceMax=',mySigma%mySegment(iSeg)%HeatSourceMax
     write(*,'(A,I1.1,A,ES12.4)') " mySIGMA%Segment(",iSeg,')%VolumetricHeatSourceMin=',mySigma%mySegment(iSeg)%HeatSourceMin
     
     IF (ADJUSTL(TRIM(mySigma%mySegment(iSeg)%ART)).eq."STL") THEN
      write(*,'(A,I1.1,A,I0)') " mySIGMA%Segment(",iSeg,")nOFFfiles=",mySigma%mySegment(iSeg)%nOFFfiles
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
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%Density(",iMat,')=',myMaterials(iMat)%density
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%DensitySlope(",iMat,')=',myMaterials(iMat)%densitySteig
     write(*,'(A,I2.2,A,ES12.4)') "myThermodyn%ALPHA(",iMat,')=',myMaterials(iMat)%Alpha
     write(*,'(A)') 
    END DO
    write(*,*) 
    write(*,'(A,ES12.4)') "myProcess%AmbientTemperature=",myProcess%AmbientTemperature 
    write(*,'(A,3ES12.4)') "myProcess%TemperatureSensorCoor=",myProcess%TemperatureSensorCoor
    write(*,'(A,3ES12.4)') "myProcess%TemperatureSensorRadius=",myProcess%TemperatureSensorRadius
    write(*,'(A,ES12.4)') "myProcess%HeatTransferCoeff=", myProcess%HeatTransferCoeff    
    write(*,'(A,ES12.4)') "myProcess%ConductiveLambda=",myProcess%ConductiveLambda 
    write(*,'(A,ES12.4)') "myProcess%ConductiveGradient=", myProcess%ConductiveGradient
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
