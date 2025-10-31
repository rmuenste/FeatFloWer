!=========================================================================
! QuadSc_force_torque_calc.f90
!
! Force and torque computation functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE FBM_GetForces()

  CALL EvaluateDragLift(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
    LinSc%valP(NLMAX)%x,Viscosity,FictKNPR)

END SUBROUTINE FBM_GetForces
!=========================================================================
!
!=========================================================================
SUBROUTINE FAC_GetForcesParT(mfile,iT)
  INTEGER mfile,iT
  REAL*8 :: Force2(3),ForceV(3),ForceP(3),Force(3),Factor,PI = 3.141592654d0
  REAL*8 :: Scale
  REAL*8 :: U_mean=1.0d0,H=0.20d0,D=1d0
  INTEGER i,nn
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
    CALL EvaluateDragLift9_ParT_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP,iT)
  ELSE
    CALL EvaluateDragLift_ParT_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP,iT)
  END IF

  if(bViscoElastic)then
    Force = ForceV + ForceP + ViscoElasticForce
    IF (b2DViscoBench) THEN
     Scale = 2d0/(U_mean*U_mean*D*H)
     Force = Scale*Force
     ViscoElasticForce = (ViscoElasticForce)*Scale
    END IF
    IF (b3DViscoBench) THEN
     Scale = 6d0*PI*postParams%Sc_Mu*postParams%Sc_U*postParams%Sc_a
     Force = (4d0*Force)/Scale
     ViscoElasticForce = (4d0*ViscoElasticForce)/Scale
    END IF
  else
    Factor = 2d0/(postParams%U_mean*postParams%U_mean*postParams%D*postParams%H)
    ForceP = Factor*ForceP
    ForceV = Factor*ForceV
  end if

  IF (myid.eq.showID) THEN
    if(bViscoElastic)then
      WRITE(MTERM,5)
      WRITE(MFILE,5)
      IF (b2DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
      IF (b3DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
    else
      WRITE(MTERM,5)
      WRITE(MFILE,5)
       write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
       write(*,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
    end if
      write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(*,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(mfile,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
      write(*,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
      WRITE(MTERM,5)
      WRITE(MFILE,5)
  END IF

  5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForcesParT
!=========================================================================
!
!=========================================================================
SUBROUTINE FAC_GetSurfForces(mfile)

 INTEGER mfile
 EXTERNAL E013

 ILEV=NLMAX
 CALL SETLEV(2)

 CALL GetForcesOnSubmeshX(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                          LinSc%P_new,&
                          mg_Mesh%level(ILEV)%kvert,&
                          mg_Mesh%level(ILEV)%karea,&
                          mg_Mesh%level(ILEV)%kedge,&
                          BndrForce,&
                          mg_Mesh%level(ILEV)%dcorvg,&
                          Properties%Viscosity(1),mfile,&
                          E013)

END SUBROUTINE FAC_GetSurfForces
!=========================================================================
!
!=========================================================================
SUBROUTINE FAC_GetForces(mfile)
  INTEGER mfile
  REAL*8 :: Force2(3),ForceV(3),ForceP(3),Force(3),Factor,PI = 3.141592654d0
  REAL*8 :: Scale
  REAL*8 :: U_mean=1.0d0,H=0.20d0,D=1d0
  INTEGER i,nn
  EXTERNAL E013

  ILEV=NLMAX
  CALL SETLEV(2)

  IF (bNonNewtonian.AND.myMatrixRenewal%S.NE.0) THEN
    CALL EvaluateDragLift9_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  ELSE
    CALL EvaluateDragLift_old(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
      LinSc%P_new,BndrForce,ForceV,ForceP)
  END IF

  if(bViscoElastic)then
    Force = ForceV + ForceP + ViscoElasticForce
    IF (b2DViscoBench) THEN
     Scale = 2d0/(U_mean*U_mean*D*H)
     Force = Scale*Force
     ViscoElasticForce = (ViscoElasticForce)*Scale
    END IF
    IF (b3DViscoBench) THEN
     Scale = 6d0*PI*postParams%Sc_Mu*postParams%Sc_U*postParams%Sc_a
     Force = (4d0*Force)/Scale
     ViscoElasticForce = (4d0*ViscoElasticForce)/Scale
    END IF
  else
    Factor = 2d0/(postParams%U_mean*postParams%U_mean*postParams%D*postParams%H)
    ForceP = Factor*ForceP
    ForceV = Factor*ForceV
  end if

  IF (myid.eq.showID) THEN

    if(bViscoElastic)then
      WRITE(MTERM,5)
      WRITE(MFILE,5)
      IF (b2DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(1),(Force(1)-ViscoElasticForce(1)),Force(1)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
      IF (b3DViscoBench) THEN
       write(mfile,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       write(mterm,'(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):",&
         timens,ViscoElasticForce(3),(Force(3)-ViscoElasticForce(3)),Force(3)
       WRITE(666,'(10ES13.5)')timens,ViscoElasticForce,&
         (Force-ViscoElasticForce),Force
      END IF
    else
      WRITE(MTERM,5)
      WRITE(MFILE,5)
       write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
       write(*,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
    end if
      write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(*,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
      write(mfile,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
      write(*,'(A12,3ES15.7E2)') "BenchForce: ",timens,ForceV(1:2)+forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "FBMForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2)

!      write(mfile,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)

      WRITE(666,'(10ES16.8)') Timens,ForceV+forceP,ForceV,forceP
      WRITE(MTERM,5)
      WRITE(MFILE,5)
  END IF

  5  FORMAT(104('-'))

END SUBROUTINE FAC_GetForces
!=========================================================================
!
!=========================================================================
SUBROUTINE Calculate_Torque(mfile)
implicit none
INTEGER mfile,i,iSeg
REAL*8 Torque1(3), Torque2(3),dVolFlow1,dVolFlow2,myPI,daux
REAL*8 dHeat,Ml_i,Shear,Visco,dVol,dArea1,dArea2,dArea3
REAL*8 dIntPres1,dIntPres2,dIntPres3,dPressureDifference,zMin, zMax,dS
REAL*8, allocatable :: SegmentForce(:,:)

integer :: ilevel

EXTERNAL E013

ilevel = mg_mesh%nlmax

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,101)

 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque2, E013,102)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,103)
END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 call GetTorqueMixer(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                     LinSc%ValP(NLMAX)%x,MixerKNPR,& !How separate????
                     mg_mesh%level(ilevel)%kvert,&
                     mg_mesh%level(ilevel)%karea,&
                     mg_mesh%level(ilevel)%kedge,&
                     mg_mesh%level(ilevel)%dcorvg,&
                     Viscosity,Torque1, E013,103)

 if (.not.allocated(SegmentForce) )allocate(SegmentForce(3,mySigma%NumberOfSeg))
 DO iSeg=1,mySigma%NumberOfSeg
  call GetSegmentForce(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                      LinSc%ValP(NLMAX)%x,mySegmentIndicator,&
                      mg_mesh%level(ilevel)%kvert,&
                      mg_mesh%level(ilevel)%karea,&
                      mg_mesh%level(ilevel)%kedge,&
                      mg_mesh%level(ilevel)%dcorvg,&
                      Viscosity,SegmentForce(:,iSeg), E013,iSeg)

 END DO

END IF

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
 IF (myid.ne.0) then
  call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow1,0)

  call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow2,1)
 END IF
ELSE
 IF (myid.ne.0) then
  call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow1,mySigma%L0)

  call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dVolFlow2,mySigma%L0+mySigma%L)
 END IF
END IF

IF (myid.ne.0) then
 call IntegratePressureAtInflow(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dIntPres3,dArea3)

 call IntegratePressure(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dIntPres1,dArea1,mySigma%L0)

 call IntegratePressure(mg_mesh%level(ilevel)%dcorvg,&
                        mg_mesh%level(ilevel)%karea,&
                        mg_mesh%level(ilevel)%kvert,&
                        mg_mesh%level(ilevel)%nel,&
                        dIntPres2,dArea2,mySigma%L0+mySigma%L)
END IF

dHeat = 0d0
dVol  = 0d0

if (myid.ne.0) then
 DO i=1,QuadSc%ndof
  IF (MixerKNPR(i).eq.0) THEN
   Shear = Shearrate(i)
   Visco = 0.1d0*Viscosity(i)
   Ml_i = mg_MlRhoMat(NLMAX)%a(i)*1e-6
   dHeat = dHeat + Ml_i*Shear*Shear*Visco
   dVol = dVol + mg_MlRhoMat(NLMAX)%a(i)*1e-3
  END IF
 END DO
END IF

CALL COMM_SUMM(dVolFlow1)
CALL COMM_SUMM(dVolFlow2)
CALL COMM_SUMM(dArea1)
CALL COMM_SUMM(dArea2)
CALL COMM_SUMM(dArea3)
CALL COMM_SUMM(dIntPres1)
CALL COMM_SUMM(dIntPres2)
CALL COMM_SUMM(dIntPres3)
CALL COMM_SUMM(dHeat)
CALL COMM_SUMM(dVol)


dPressureDifference = ((dIntPres2/dArea2) - (dIntPres1/dArea1))*1e-6 !!! [Bar]

myPI = dATAN(1d0)*4d0
daux = 1D0*1e-7*myPI*(myProcess%umdr/3d1)

dHeat = dHeat / myThermodyn%density
dVol = dVol / myThermodyn%density

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN

 IF     (ADJUSTL(TRIM(myProcess%pTYPE)).eq."THROUGHPUT") THEN
  mySSE_covergence%Monitor(mySSE_covergence%iC) = dPressureDifference
 ELSEIF (ADJUSTL(TRIM(myProcess%pTYPE)).eq."PRESSUREDROP") THEN
  mySSE_covergence%Monitor(mySSE_covergence%iC) = dVolFlow2
 ELSE
  STOP
 END IF
 mySSE_covergence%iC = MOD(mySSE_covergence%iC,mySSE_covergence%nC) + 1

 IF (itns.gt.mySSE_covergence%nC) THEN
  mySSE_covergence%average = 0d0
  DO i=1,mySSE_covergence%nC
   mySSE_covergence%average = mySSE_covergence%average +  mySSE_covergence%Monitor(i)
  END DO
  mySSE_covergence%average = mySSE_covergence%average/DBLE(mySSE_covergence%nC)
  if (mySSE_covergence%average.lt.0d0) then
   mySSE_covergence%average = MIN(mySSE_covergence%average,(1d1*mySSE_covergence%dCharVisco*1d-5))
  else
   mySSE_covergence%average = MAX(mySSE_covergence%average,(1d1*mySSE_covergence%dCharVisco*1d-5))
  end if

  mySSE_covergence%std_dev = 0d0
  DO i=1,mySSE_covergence%nC
   mySSE_covergence%std_dev = mySSE_covergence%std_dev +  (mySSE_covergence%Monitor(i)-mySSE_covergence%average)**2d0
  END DO
  mySSE_covergence%std_dev = (mySSE_covergence%std_dev/DBLE(mySSE_covergence%nC))**0.5d0
  mySSE_covergence%std_dev = ABS(1d2*mySSE_covergence%std_dev/mySSE_covergence%average)
 END IF

 mySetup%bPressureConvergence = .false.

 if (itns.gt.mySSE_covergence%start) then
  if (mySSE_covergence%std_dev.lt.2d-1) then
   mySetup%bPressureConvergence = .true.
  end if
 end if

END IF

IF (myid.eq.showID) THEN
  WRITE(MTERM,5)
  WRITE(MFILE,5)
  write(mfile,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mterm,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mfile,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
  write(mterm,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
  IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE") THEN
    write(mfile,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
    write(mterm,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
  END IF
  IF (ADJUSTL(TRIM(mySigma%cType)).EQ."SSE") THEN
    write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
  END IF
  IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
    write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
  END IF

  write(mfile,'(A,7ES14.4)') "PressureDiffAtInflow[bar]_&_Area[mm2]: ",timens,(dIntPres3/dArea3)*1e-6,1d2*dArea3
  write(mterm,'(A,7ES14.4)') "PressureDiffAtInflow[bar]_&_Area[mm2]: ",timens,(dIntPres3/dArea3)*1e-6,1d2*dArea3

 IF (ADJUSTL(TRIM(mySigma%cType)).EQ."DIE") THEN
  DO iSeg=1,mySigma%NumberOfSeg
    dS = 0d0
    if (mySigma%mySegment(iSeg)%ObjectType.eq.'DIE') dS = -1d0
    if (mySigma%mySegment(iSeg)%ObjectType.eq.'OBSTACLE') dS = +1d0

    write(mfile,'(A,I0,A,4ES14.4)') "Force_acting_on_Segment",iSeg,"_[kN]:",timens,dS*1e-3*1e-1*1e-4*SegmentForce(:,iSeg)
    write(mterm,'(A,I0,A,4ES14.4)') "Force_acting_on_Segment",iSeg,"_[kN]:",timens,dS*1e-3*1e-1*1e-4*SegmentForce(:,iSeg)
  END DO
 END IF

END IF

5  FORMAT(100('-'))

END SUBROUTINE Calculate_Torque
!=========================================================================
!
!=========================================================================
SUBROUTINE DNA_GetTorques(mfile)
integer mfile
REAL*8 Torque1(3),daux

external E013

ilev = nlmax
call setlev(2)

call GetDNATorque(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                  LinSc%ValP(NLMAX)%x,BndrForce,& !How separate????
                  mg_mesh%level(ilev)%kvert,&
                  mg_mesh%level(ilev)%karea,&
                  mg_mesh%level(ilev)%kedge,&
                  mg_mesh%level(ilev)%dcorvg,&
                  Viscosity,Torque1, E013,1)

!daux = 1/(mu*rho*OMEGA*R_i^2*H)
daux = 1d0/(Properties%Viscosity(1)*Properties%Density(1)*(2*3.14d0*(1d0/60d0))*(0.2d0**2d0)*0.8d0)

IF (myid.eq.1) then
 write(mfile,'(A,4ES14.4)') "Torque acting on surface:",timens,daux*Torque1(1:3)
 write(mterm,'(A,4ES14.4)') "Torque acting on surface:",timens,daux*Torque1(1:3)
end if

END SUBROUTINE DNA_GetTorques
!=========================================================================
!
!=========================================================================
SUBROUTINE Get_DissipationIntegral(mfile)
integer mfile
REAL*8 Torque1(3),daux
real*8 :: area, G, L, U

external E013

ilev = nlmax
call setlev(2)

call CalculateDissipationIntegralNumerator(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                  mg_mesh%level(ilev)%kvert,&
                  mg_mesh%level(ilev)%karea,&
                  mg_mesh%level(ilev)%kedge,&
                  mg_mesh%level(ilev)%dcorvg,&
                  E013)


END SUBROUTINE Get_DissipationIntegral
!=========================================================================
!
!=========================================================================
SUBROUTINE DNA_GetSoosForce(mfile)
integer mfile
REAL*8 Torque1(3),daux
real*8 :: area, G, L, U

external E013

ilev = nlmax
call setlev(2)

call GetSoosForce(QuadSc%valU,QuadSc%valV,QuadSc%valW,&
                  LinSc%ValP(NLMAX)%x,BndrForce,& !How separate????
                  mg_mesh%level(ilev)%kvert,&
                  mg_mesh%level(ilev)%karea,&
                  mg_mesh%level(ilev)%kedge,&
                  mg_mesh%level(ilev)%dcorvg,&
                  Viscosity,Torque1, E013,1)

! Cube side length
L = 0.1

! Top surface area
area = L**2

! Shear velocity
U = 1.0

! Uniform shear rate
G = U / L

daux = 1./ (G * area)

IF (myid.eq.1) then
 write(mfile,'(A,4ES14.4)') "Shear Stress acting on top wall: ",timens,daux*Torque1(1:3)
 write(mterm,'(A,4ES14.4)') "Shear Stress acting on top wall: ",timens,daux*Torque1(1:3)
end if

END SUBROUTINE DNA_GetSoosForce
