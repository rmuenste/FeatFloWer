!=========================================================================
! QuadSc_force_torque_calc.f90
!
! Force and torque computation functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================

!=========================================================================
! FBM_GetForces - Compute forces using Fictitious Boundary Method
! Evaluates drag and lift forces on fictitious boundaries
!=========================================================================
subroutine FBM_GetForces()
  implicit none

  call EvaluateDragLift(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                        LinSc%valP(NLMAX)%x, Viscosity, FictKNPR)

end subroutine FBM_GetForces

!=========================================================================
! FAC_GetForcesParT - Compute forces for particle/object benchmarks
! Calculates dimensionless drag and lift coefficients for flow around particles
! Supports both Newtonian and non-Newtonian fluids with viscoelastic effects
!=========================================================================
subroutine FAC_GetForcesParT(mfile, iT)
  implicit none
  integer, intent(in) :: mfile, iT
  integer :: i, nn
  real(8) :: Force2(3), ForceV(3), ForceP(3), Force(3), Factor
  real(8) :: Scale
  real(8), parameter :: PI = 3.141592654d0
  real(8) :: U_mean = 1.0d0, H = 0.20d0, D = 1.0d0
  external E013

  ILEV = NLMAX
  call SETLEV(2)

  ! Note: Using non-ParT versions as ParT-specific functions don't exist
  ! The iT parameter is not used by these force calculation routines
  if (bNonNewtonian .and. myMatrixRenewal%S /= 0) then
    call EvaluateDragLift9_old(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                               LinSc%P_new, BndrForce, ForceV, ForceP)
  else
    call EvaluateDragLift_old(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                              LinSc%P_new, BndrForce, ForceV, ForceP)
  end if

  if (bViscoElastic) then
    Force = ForceV + ForceP + ViscoElasticForce
    if (b2DViscoBench) then
      Scale = 2d0 / (U_mean * U_mean * D * H)
      Force = Scale * Force
      ViscoElasticForce = (ViscoElasticForce) * Scale
    end if
    if (b3DViscoBench) then
      Scale = 6d0 * PI * postParams%Sc_Mu * postParams%Sc_U * postParams%Sc_a
      Force = (4d0 * Force) / Scale
      ViscoElasticForce = (4d0 * ViscoElasticForce) / Scale
    end if
  else
    Factor = 2d0 / (postParams%U_mean * postParams%U_mean * postParams%D * postParams%H)
    ForceP = Factor * ForceP
    ForceV = Factor * ForceV
  end if

  if (myid .eq. showID) then
    if (bViscoElastic) then
      write (MTERM, 5)
      write (MFILE, 5)
      if (b2DViscoBench) then
        write (mfile, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(1), (Force(1) - ViscoElasticForce(1)), Force(1)
        write (mterm, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(1), (Force(1) - ViscoElasticForce(1)), Force(1)
        write (666, '(10ES13.5)') timens, ViscoElasticForce, &
          (Force - ViscoElasticForce), Force
      end if
      if (b3DViscoBench) then
        write (mfile, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(3), (Force(3) - ViscoElasticForce(3)), Force(3)
        write (mterm, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(3), (Force(3) - ViscoElasticForce(3)), Force(3)
        write (666, '(10ES13.5)') timens, ViscoElasticForce, &
          (Force - ViscoElasticForce), Force
      end if
    else
      write (MTERM, 5)
      write (MFILE, 5)
      write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
  write (*, '(A7,7A15)') "Force: ", "Time", "C_D", "C_L", "ForceVx", "ForceVy", "ForcePx", "ForcePy"
    end if
    write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
  write (*, '(A7,7ES15.7E2)') "Force: ", timens, ForceV(1:2) + forceP(1:2), ForceV(1:2), forceP(1:2)
    write (mfile, '(A12,3ES15.7E2)') "BenchForce: ", timens, ForceV(1:2) + forceP(1:2)
    write (*, '(A12,3ES15.7E2)') "BenchForce: ", timens, ForceV(1:2) + forceP(1:2)
    write (MTERM, 5)
    write (MFILE, 5)
  end if

5 format(104('-'))

end subroutine FAC_GetForcesParT
!=========================================================================
!
!=========================================================================
subroutine FAC_GetSurfForces(mfile)

  integer mfile
  external E013

  ILEV = NLMAX
  call SETLEV(2)

  call GetForcesOnSubmeshX(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                           LinSc%P_new, &
                           mg_Mesh%level(ILEV)%kvert, &
                           mg_Mesh%level(ILEV)%karea, &
                           mg_Mesh%level(ILEV)%kedge, &
                           BndrForce, &
                           mg_Mesh%level(ILEV)%dcorvg, &
                           Properties%Viscosity(1), mfile, &
                           E013)

end subroutine FAC_GetSurfForces
!=========================================================================
!
!=========================================================================
subroutine FAC_GetForces(mfile)
  integer mfile
  real*8 :: Force2(3), ForceV(3), ForceP(3), Force(3), Factor, PI = 3.141592654d0
  real*8 :: Scale
  real*8 :: U_mean = 1.0d0, H = 0.20d0, D = 1d0
  integer i, nn
  external E013

  ILEV = NLMAX
  call SETLEV(2)

  if (bNonNewtonian .AND. myMatrixRenewal%S .NE. 0) then
    call EvaluateDragLift9_old(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                               LinSc%P_new, BndrForce, ForceV, ForceP)
  else
    call EvaluateDragLift_old(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                              LinSc%P_new, BndrForce, ForceV, ForceP)
  end if

  if (bViscoElastic) then
    Force = ForceV + ForceP + ViscoElasticForce
    if (b2DViscoBench) then
      Scale = 2d0 / (U_mean * U_mean * D * H)
      Force = Scale * Force
      ViscoElasticForce = (ViscoElasticForce) * Scale
    end if
    if (b3DViscoBench) then
      Scale = 6d0 * PI * postParams%Sc_Mu * postParams%Sc_U * postParams%Sc_a
      Force = (4d0 * Force) / Scale
      ViscoElasticForce = (4d0 * ViscoElasticForce) / Scale
    end if
  else
    Factor = 2d0 / (postParams%U_mean * postParams%U_mean * postParams%D * postParams%H)
    ForceP = Factor * ForceP
    ForceV = Factor * ForceV
  end if

  if (myid .eq. showID) then

    if (bViscoElastic) then
      write (MTERM, 5)
      write (MFILE, 5)
      if (b2DViscoBench) then
        write (mfile, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(1), (Force(1) - ViscoElasticForce(1)), Force(1)
        write (mterm, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(1), (Force(1) - ViscoElasticForce(1)), Force(1)
        write (666, '(10ES13.5)') timens, ViscoElasticForce, &
          (Force - ViscoElasticForce), Force
      end if
      if (b3DViscoBench) then
        write (mfile, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(3), (Force(3) - ViscoElasticForce(3)), Force(3)
        write (mterm, '(A,10ES13.5)') "TimevsForce(Visco,Hydro,Full):", &
          timens, ViscoElasticForce(3), (Force(3) - ViscoElasticForce(3)), Force(3)
        write (666, '(10ES13.5)') timens, ViscoElasticForce, &
          (Force - ViscoElasticForce), Force
      end if
    else
      write (MTERM, 5)
      write (MFILE, 5)
      write(mfile,'(A7,7A15)') "Force: ","Time","C_D","C_L","ForceVx","ForceVy","ForcePx","ForcePy"
  write (*, '(A7,7A15)') "Force: ", "Time", "C_D", "C_L", "ForceVx", "ForceVy", "ForcePx", "ForcePy"
    end if
    write(mfile,'(A7,7ES15.7E2)') "Force: ",timens,ForceV(1:2)+forceP(1:2),ForceV(1:2),forceP(1:2)
  write (*, '(A7,7ES15.7E2)') "Force: ", timens, ForceV(1:2) + forceP(1:2), ForceV(1:2), forceP(1:2)
    write (mfile, '(A12,3ES15.7E2)') "BenchForce: ", timens, ForceV(1:2) + forceP(1:2)
    write (*, '(A12,3ES15.7E2)') "BenchForce: ", timens, ForceV(1:2) + forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "FBMForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2)

!      write(mfile,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)
!      write(*,'(A12,3ES15.7E2)') "TotalForce: ",timens,myFBM%ParticleNew(1)%ResistanceForce(1:2) + &
!      ForceV(1:2)+forceP(1:2)

    write (666, '(10ES16.8)') Timens, ForceV + forceP, ForceV, forceP
    write (MTERM, 5)
    write (MFILE, 5)
  end if

5 format(104('-'))

end subroutine FAC_GetForces
!=========================================================================
!
!=========================================================================
subroutine Calculate_Torque(mfile)
  implicit none
  integer mfile, i, iSeg
  real * 8 Torque1(3), Torque2(3), dVolFlow1, dVolFlow2, myPI, daux
  real * 8 dHeat, Ml_i, Shear, Visco, dVol, dArea1, dArea2, dArea3
  real * 8 dIntPres1, dIntPres2, dIntPres3, dPressureDifference, zMin, zMax, dS
  real*8, allocatable :: SegmentForce(:, :)

  integer :: ilevel

  external E013

  ilevel = mg_mesh%nlmax

  if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "TSE") then
    call GetTorqueMixer(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                        LinSc%ValP(NLMAX)%x, MixerKNPR, & !How separate????
                        mg_mesh%level(ilevel)%kvert, &
                        mg_mesh%level(ilevel)%karea, &
                        mg_mesh%level(ilevel)%kedge, &
                        mg_mesh%level(ilevel)%dcorvg, &
                        Viscosity, Torque1, E013, 101)

    call GetTorqueMixer(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                        LinSc%ValP(NLMAX)%x, MixerKNPR, & !How separate????
                        mg_mesh%level(ilevel)%kvert, &
                        mg_mesh%level(ilevel)%karea, &
                        mg_mesh%level(ilevel)%kedge, &
                        mg_mesh%level(ilevel)%dcorvg, &
                        Viscosity, Torque2, E013, 102)
  end if

  if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "SSE") then
    call GetTorqueMixer(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                        LinSc%ValP(NLMAX)%x, MixerKNPR, & !How separate????
                        mg_mesh%level(ilevel)%kvert, &
                        mg_mesh%level(ilevel)%karea, &
                        mg_mesh%level(ilevel)%kedge, &
                        mg_mesh%level(ilevel)%dcorvg, &
                        Viscosity, Torque1, E013, 103)
  end if

  if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "DIE") then
    call GetTorqueMixer(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                        LinSc%ValP(NLMAX)%x, MixerKNPR, & !How separate????
                        mg_mesh%level(ilevel)%kvert, &
                        mg_mesh%level(ilevel)%karea, &
                        mg_mesh%level(ilevel)%kedge, &
                        mg_mesh%level(ilevel)%dcorvg, &
                        Viscosity, Torque1, E013, 103)

    if (.not. allocated(SegmentForce)) allocate (SegmentForce(3, mySigma%NumberOfSeg))
    do iSeg = 1, mySigma%NumberOfSeg
      call GetSegmentForce(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                           LinSc%ValP(NLMAX)%x, mySegmentIndicator, &
                           mg_mesh%level(ilevel)%kvert, &
                           mg_mesh%level(ilevel)%karea, &
                           mg_mesh%level(ilevel)%kedge, &
                           mg_mesh%level(ilevel)%dcorvg, &
                           Viscosity, SegmentForce(:, iSeg), E013, iSeg)

    end do

  end if

  if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "DIE") then
    if (myid .ne. 0) then
      call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg, &
                                  mg_mesh%level(ilevel)%karea, &
                                  mg_mesh%level(ilevel)%kvert, &
                                  mg_mesh%level(ilevel)%nel, &
                                  dVolFlow1, 0)

      call Integrate_DIE_Flowrate(mg_mesh%level(ilevel)%dcorvg, &
                                  mg_mesh%level(ilevel)%karea, &
                                  mg_mesh%level(ilevel)%kvert, &
                                  mg_mesh%level(ilevel)%nel, &
                                  dVolFlow2, 1)
    end if
  else
    if (myid .ne. 0) then
      call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg, &
                             mg_mesh%level(ilevel)%karea, &
                             mg_mesh%level(ilevel)%kvert, &
                             mg_mesh%level(ilevel)%nel, &
                             dVolFlow1, mySigma%L0)

      call IntegrateFlowrate(mg_mesh%level(ilevel)%dcorvg, &
                             mg_mesh%level(ilevel)%karea, &
                             mg_mesh%level(ilevel)%kvert, &
                             mg_mesh%level(ilevel)%nel, &
                             dVolFlow2, mySigma%L0 + mySigma%L)
    end if
  end if

  if (myid .ne. 0) then
    call IntegratePressureAtInflow(mg_mesh%level(ilevel)%dcorvg, &
                                   mg_mesh%level(ilevel)%karea, &
                                   mg_mesh%level(ilevel)%kvert, &
                                   mg_mesh%level(ilevel)%nel, &
                                   dIntPres3, dArea3)

    call IntegratePressure(mg_mesh%level(ilevel)%dcorvg, &
                           mg_mesh%level(ilevel)%karea, &
                           mg_mesh%level(ilevel)%kvert, &
                           mg_mesh%level(ilevel)%nel, &
                           dIntPres1, dArea1, mySigma%L0)

    call IntegratePressure(mg_mesh%level(ilevel)%dcorvg, &
                           mg_mesh%level(ilevel)%karea, &
                           mg_mesh%level(ilevel)%kvert, &
                           mg_mesh%level(ilevel)%nel, &
                           dIntPres2, dArea2, mySigma%L0 + mySigma%L)
  end if

  dHeat = 0d0
  dVol = 0d0

  if (myid .ne. 0) then
    do i = 1, QuadSc%ndof
      if (MixerKNPR(i) .eq. 0) then
        Shear = Shearrate(i)
        Visco = 0.1d0 * Viscosity(i)
        Ml_i = mg_MlRhoMat(NLMAX)%a(i) * 1e-6
        dHeat = dHeat + Ml_i * Shear * Shear * Visco
        dVol = dVol + mg_MlRhoMat(NLMAX)%a(i) * 1e-3
      end if
    end do
  end if

  call COMM_SUMM(dVolFlow1)
  call COMM_SUMM(dVolFlow2)
  call COMM_SUMM(dArea1)
  call COMM_SUMM(dArea2)
  call COMM_SUMM(dArea3)
  call COMM_SUMM(dIntPres1)
  call COMM_SUMM(dIntPres2)
  call COMM_SUMM(dIntPres3)
  call COMM_SUMM(dHeat)
  call COMM_SUMM(dVol)

  dPressureDifference = ((dIntPres2 / dArea2) - (dIntPres1 / dArea1)) * 1e-6 !!! [Bar]

  myPI = dATAN(1d0) * 4d0
  daux = 1D0 * 1e-7 * myPI * (myProcess%umdr / 3d1)

  dHeat = dHeat / myThermodyn%density
  dVol = dVol / myThermodyn%density

IF (ADJUSTL(TRIM(mySigma%cType)).EQ."TSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."SSE".or.ADJUSTL(TRIM(mySigma%cType)).EQ."XSE") THEN

    if (ADJUSTL(TRIM(myProcess%pTYPE)) .eq. "THROUGHPUT") then
      mySSE_covergence%Monitor(mySSE_covergence%iC) = dPressureDifference
    elseif (ADJUSTL(TRIM(myProcess%pTYPE)) .eq. "PRESSUREDROP") then
      mySSE_covergence%Monitor(mySSE_covergence%iC) = dVolFlow2
    else
      stop
    end if
    mySSE_covergence%iC = MOD(mySSE_covergence%iC, mySSE_covergence%nC) + 1

    if (itns .gt. mySSE_covergence%nC) then
      mySSE_covergence%average = 0d0
      do i = 1, mySSE_covergence%nC
        mySSE_covergence%average = mySSE_covergence%average + mySSE_covergence%Monitor(i)
      end do
      mySSE_covergence%average = mySSE_covergence%average / DBLE(mySSE_covergence%nC)
      if (mySSE_covergence%average .lt. 0d0) then
mySSE_covergence%average = MIN(mySSE_covergence%average, (1d1 * mySSE_covergence%dCharVisco * 1d-5))
      else
mySSE_covergence%average = MAX(mySSE_covergence%average, (1d1 * mySSE_covergence%dCharVisco * 1d-5))
      end if

      mySSE_covergence%std_dev = 0d0
      do i = 1, mySSE_covergence%nC
   mySSE_covergence%std_dev = mySSE_covergence%std_dev +  (mySSE_covergence%Monitor(i)-mySSE_covergence%average)**2d0
      end do
      mySSE_covergence%std_dev = (mySSE_covergence%std_dev / DBLE(mySSE_covergence%nC))**0.5d0
      mySSE_covergence%std_dev = ABS(1d2 * mySSE_covergence%std_dev / mySSE_covergence%average)
    end if

    mySetup%bPressureConvergence = .false.

    if (itns .gt. mySSE_covergence%start) then
      if (mySSE_covergence%std_dev .lt. 2d-1) then
        mySetup%bPressureConvergence = .true.
      end if
    end if

  end if

  if (myid .eq. showID) then
    write (MTERM, 5)
    write (MFILE, 5)
  write(mfile,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mterm,'(A,6ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_Volume_[l]_&_RT_[s]:",timens,dVolFlow1*3.6d0,dVolFlow1*3.6d0*myThermodyn%density,dVol,dVol/dVolFlow1*1000d0
  write(mfile,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
  write(mterm,'(A,7ES14.4)') "Throughput_[l/h]_&_[kg/h]_&_PressureDiff[bar]_stdERR[%]: ",timens,dVolFlow2*3.6d0,dVolFlow2*3.6d0*myThermodyn%density,dPressureDifference,mySSE_covergence%std_dev
    if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "TSE") then
    write(mfile,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
    write(mterm,'(A,3ES14.4,A,ES14.4)') "Power_acting_on_the_screws_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),1e-3*daux*Torque2(3),' & ',1e-3*dHeat
    end if
    if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "SSE") then
    write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    end if
    if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "DIE") then
    write(mfile,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    write(mterm,'(A,2ES14.4,A,ES14.4)') "Power_acting_on_the_screw_[kW]_&_heat_generation_rate_[kW]:",timens,1e-3*daux*Torque1(3),' & ',1e-3*dHeat
    end if

  write(mfile,'(A,7ES14.4)') "PressureDiffAtInflow[bar]_&_Area[mm2]: ",timens,(dIntPres3/dArea3)*1e-6,1d2*dArea3
  write(mterm,'(A,7ES14.4)') "PressureDiffAtInflow[bar]_&_Area[mm2]: ",timens,(dIntPres3/dArea3)*1e-6,1d2*dArea3

    if (ADJUSTL(TRIM(mySigma%cType)) .EQ. "DIE") then
      do iSeg = 1, mySigma%NumberOfSeg
        dS = 0d0
        if (mySigma%mySegment(iSeg)%ObjectType .eq. 'DIE') dS = -1d0
        if (mySigma%mySegment(iSeg)%ObjectType .eq. 'OBSTACLE') dS = +1d0

    write(mfile,'(A,I0,A,4ES14.4)') "Force_acting_on_Segment",iSeg,"_[kN]:",timens,dS*1e-3*1e-1*1e-4*SegmentForce(:,iSeg)
    write(mterm,'(A,I0,A,4ES14.4)') "Force_acting_on_Segment",iSeg,"_[kN]:",timens,dS*1e-3*1e-1*1e-4*SegmentForce(:,iSeg)
      end do
    end if

  end if

5 format(100('-'))

end subroutine Calculate_Torque
!=========================================================================
!
!=========================================================================
subroutine DNA_GetTorques(mfile)
  integer mfile
  real * 8 Torque1(3), daux

  external E013

  ilev = nlmax
  call setlev(2)

  call GetDNATorque(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                    LinSc%ValP(NLMAX)%x, BndrForce, & !How separate????
                    mg_mesh%level(ilev)%kvert, &
                    mg_mesh%level(ilev)%karea, &
                    mg_mesh%level(ilev)%kedge, &
                    mg_mesh%level(ilev)%dcorvg, &
                    Viscosity, Torque1, E013, 1)

!daux = 1/(mu*rho*OMEGA*R_i^2*H)
 daux = 1d0/(Properties%Viscosity(1)*Properties%Density(1)*(2*3.14d0*(1d0/60d0))*(0.2d0**2d0)*0.8d0)

  if (myid .eq. 1) then
    write (mfile, '(A,4ES14.4)') "Torque acting on surface:", timens, daux * Torque1(1:3)
    write (mterm, '(A,4ES14.4)') "Torque acting on surface:", timens, daux * Torque1(1:3)
  end if

end subroutine DNA_GetTorques
!=========================================================================
!
!=========================================================================
subroutine Get_DissipationIntegral(mfile)
  integer mfile
  real * 8 Torque1(3), daux
  real*8 :: area, G, L, U

  external E013

  ilev = nlmax
  call setlev(2)

  call CalculateDissipationIntegralNumerator(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                                             mg_mesh%level(ilev)%kvert, &
                                             mg_mesh%level(ilev)%karea, &
                                             mg_mesh%level(ilev)%kedge, &
                                             mg_mesh%level(ilev)%dcorvg, &
                                             E013)

end subroutine Get_DissipationIntegral
!=========================================================================
!
!=========================================================================
subroutine DNA_GetSoosForce(mfile)
  integer mfile
  real * 8 Torque1(3), daux
  real*8 :: area, G, L, U

  external E013

  ilev = nlmax
  call setlev(2)

  call GetSoosForce(QuadSc%valU, QuadSc%valV, QuadSc%valW, &
                    LinSc%ValP(NLMAX)%x, BndrForce, & !How separate????
                    mg_mesh%level(ilev)%kvert, &
                    mg_mesh%level(ilev)%karea, &
                    mg_mesh%level(ilev)%kedge, &
                    mg_mesh%level(ilev)%dcorvg, &
                    Viscosity, Torque1, E013, 1)

! Cube side length
  L = 0.1

! Top surface area
  area = L**2

! Shear velocity
  U = 1.0

! Uniform shear rate
  G = U / L

  daux = 1./(G * area)

  if (myid .eq. 1) then
    write (mfile, '(A,4ES14.4)') "Shear Stress acting on top wall: ", timens, daux * Torque1(1:3)
    write (mterm, '(A,4ES14.4)') "Shear Stress acting on top wall: ", timens, daux * Torque1(1:3)
  end if

end subroutine DNA_GetSoosForce
