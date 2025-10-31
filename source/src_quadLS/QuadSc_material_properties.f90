!=========================================================================
! QuadSc_material_properties.f90
!
! Material property and density update functions
! Extracted from QuadSc_main.f90 for better code organization
!=========================================================================
!
!=========================================================================
SUBROUTINE  GetNonNewtViscosity()
INTEGER i
REAL*8 ViscosityModel
integer :: ilevel

EXTERNAL E013

ilev   = mg_mesh%nlmax
ilevel = mg_mesh%nlmax

 IF (myid.ne.0) then
  QuadSc%defU = 0d0
  QuadSc%defV = 0d0
  QuadSc%defW = 0d0
  CALL L2ProjVisco(QuadSc%ValU,QuadSc%ValV,QuadSc%ValW,&
                   Temperature,MaterialDistribution(ilev)%x,&
                   QuadSc%defU,QuadSc%defV,QuadSc%defW,&
                   mg_mesh%level(ilevel)%kvert,&
                   mg_mesh%level(ilevel)%karea,&
                   mg_mesh%level(ilevel)%kedge,&
                   mg_mesh%level(ilevel)%dcorvg,E013)

  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
!   CALL E013Sum(QuadSc%defU)
!   CALL E013Sum(QuadSc%defV)
!   CALL E013Sum(QuadSc%defW)

  if(.not.allocated(Shearrate)) allocate(Shearrate(QuadSc%ndof))

  DO i=1,QuadSc%ndof
   Shearrate(i) = QuadSc%defV(i)/QuadSc%defW(i)
   Viscosity(i) = QuadSc%defU(i)/QuadSc%defW(i)
  END DO

 END IF

END SUBROUTINE  GetNonNewtViscosity
!=========================================================================
!
!=========================================================================
SUBROUTINE FilterColdElements(mfile)
integer, intent(in) :: mfile
integer :: i
real*8 :: dj

if (myid.ne.0) then
 dj = 0
 DO i=1,QuadSc%ndof
  if (MixerKNPR(i).eq.0.and.Temperature(i).lt.80d0) then
   MixerKNPR(i) = 105
   Shell(i) = -1d0
   dj = dj + 1
  end if
 END DO
end if

CALL COMM_SUMM(dj)

if (myid.eq.1) write(mfile,*) 'Number of solidified dofs: ',nint(dj)
if (myid.eq.1) write(mterm,*) 'Number of solidified dofs: ',nint(dj)

END SUBROUTINE FilterColdElements
!=========================================================================
!
!=========================================================================
SUBROUTINE  UpdateDensityDistribution_XSE(mfile)
 INTEGER i,iel,mfile
 REAL*8 daux,taux,dAlpha
 REAL*8 AlphaViscosityMatModel,WallSlip
 REAL*8 dMaxMat,dWSFactor
 integer ifld,iMat

 if (.not.bMasterTurnedOn) return

 if (myid.eq.1) WRITE(MTERM,*) "Update of the density distribution!"
 if (myid.eq.1) WRITE(MFILE,*) "Update of the density distribution!"

 DO ILEV=NLMIN,NLMAX

  DO iel=1,mg_mesh%level(ilev)%nel

   i = mg_mesh%level(ilev)%nvt + &
       mg_mesh%level(ilev)%net + &
       mg_mesh%level(ilev)%nat + &
       iel

   taux   = Temperature(i)

   IF (myMultiMat%nOfMaterials.gt.1) THEN

    iMat = myMultiMat%InitMaterial
    dMaxMat = 1d-5
    do iFld=2,GenLinScalar%nOfFields
     if (GenLinScalar%Fld(iFld)%val(i).gt.dMaxMat) then
      iMat = iFld-1
      dMaxMat = GenLinScalar%Fld(iFld)%val(i)
     end if
    end do

   ELSE

    iMat = 1

   END IF

   IF (ADJUSTL(TRIM(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)).eq."DENSITY") THEN
    mgDensity(ILEV)%x(iel) = myMultiMat%Mat(iMat)%Thermodyn%densityT0 - taux * myMultiMat%Mat(iMat)%Thermodyn%densitySteig
   END IF
   IF (ADJUSTL(TRIM(myMultiMat%Mat(iMat)%Thermodyn%DensityModel)).eq."SPECVOLUME") THEN
    mgDensity(ILEV)%x(iel) = 1d0/(myMultiMat%Mat(iMat)%Thermodyn%densityT0 + taux * myMultiMat%Mat(iMat)%Thermodyn%densitySteig)
   END IF

  END DO

 END DO

 ! send data to the master

ILEV = LinSc%prm%MGprmIn%MedLev
IF (LinSc%prm%MGprmIn%MedLev.ge.1.and.LinSc%prm%MGprmIn%CrsSolverType.le.4) THEN
 CALL E010GATHR_L1(mgDensity(1)%x,mg_mesh%level(1)%nel)
END IF

IF (LinSc%prm%MGprmIn%MedLev.ge.2.and.LinSc%prm%MGprmIn%CrsSolverType.le.4) THEN
 CALL E010GATHR_L2(mgDensity(2)%x,mg_mesh%level(2)%nel)
END IF

IF (LinSc%prm%MGprmIn%MedLev.ge.3.and.LinSc%prm%MGprmIn%CrsSolverType.le.4) THEN
 CALL E010GATHR_L3(mgDensity(3)%x,mg_mesh%level(3)%nel)
END IF

ILEV = NLMAX

END SUBROUTINE  UpdateDensityDistribution_XSE
!=========================================================================
!
!=========================================================================
SUBROUTINE  GetAlphaNonNewtViscosity_sse()
  INTEGER i
  REAL*8 daux,taux,dAlpha
  REAL*8 AlphaViscosityMatModel,WallSlip
  REAL*8 dMaxMat,dWSFactor
  integer ifld,iMat

  if (.not.bMasterTurnedOn) return

  ILEV = NLMAX
  CALL SETLEV(2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValU)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,1)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValV)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,2)

  CALL GetGradVelo_rhs(QuadSc,QuadSc%ValW)
  CALL E013Sum3(QuadSc%defU,QuadSc%defV,QuadSc%defW)
  CALL GetGradVelo_val(QuadSc,3)

  DO i=1,SIZE(QuadSc%ValU)

   daux = QuadSc%ValUx(i)**2d0 + QuadSc%ValVy(i)**2d0 + QuadSc%ValWz(i)**2d0 + &
          0.5d0*(QuadSc%ValUy(i)+QuadSc%ValVx(i))**2d0 + &
          0.5d0*(QuadSc%ValUz(i)+QuadSc%ValWx(i))**2d0 + &
          0.5d0*(QuadSc%ValVz(i)+QuadSc%ValWy(i))**2d0

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!! tempertaure is sampled from the Temperature@q2p1_see_temp output    !!!!!!!
   taux   = Temperature(i)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   IF (myMultiMat%nOfMaterials.gt.1) THEN
    iMat = myMultiMat%InitMaterial
    dMaxMat = 1d-5
    do iFld=2,GenLinScalar%nOfFields
     if (GenLinScalar%Fld(iFld)%val(i).gt.dMaxMat) then
      iMat = iFld-1
      dMaxMat = GenLinScalar%Fld(iFld)%val(i)
     end if
    end do
    Shearrate(i) = sqrt(2d0 * daux)
    Viscosity(i) = AlphaViscosityMatModel(daux,iMat,taux)
    if (myMultiMat%Mat(iMat)%Rheology%bWallSlip) then
     dWSFactor = WallSlip(shell(i),screw(i),iMat,Viscosity(i)*Shearrate(i))
     Viscosity(i) = dWSFactor*Viscosity(i)
    END IF
   else
    Shearrate(i) = sqrt(2d0 * daux)
    Viscosity(i) = AlphaViscosityMatModel(daux,1,taux)
    if (myMultiMat%Mat(1)%Rheology%bWallSlip) then
     dWSFactor = WallSlip(shell(i),screw(i),1,Viscosity(i)*Shearrate(i))
     Viscosity(i) = dWSFactor*Viscosity(i)
    END IF
   end if

  END DO

END SUBROUTINE  GetAlphaNonNewtViscosity_sse
!=========================================================================
!
!=========================================================================
SUBROUTINE UpdateMaterialProperties
integer ii,iMat,jel(8),ndof_nel,iMaxFrac
REAL*8 dMaxFrac
REAL*8, allocatable :: dFrac(:)


if (myid.eq.1) WRITE(*,*) 'updating Material Properties Distribution... '
IF (myid.ne.0) then

!     do iel = 1,mg_mesh%level(nlmax)%nel
!      ii = mg_mesh%level(nlmax)%nvt + mg_mesh%level(nlmax)%net + mg_mesh%level(nlmax)%nat + iel
!      if (mg_mesh%level(nlmax)%dcorvg(3,ii).gt.92d0) then
!       MaterialDistribution(nlmax)%x(iel) = 2
!      end if
!     end do

!  write(*,*) MaterialDistribution(NLMAX+0-1)%x(1:knel(NLMAX+0-1))
! pause

 IF (allocated(MaterialDistribution)) then

  if (istart.eq.2) then
    do iel = 1,mg_mesh%level(nlmax-1)%nel
      jel(1) = iel
      jel(2) = mg_mesh%level(ilev+1)%kadj(3,jel(1))
      jel(3) = mg_mesh%level(ilev+1)%kadj(3,jel(2))
      jel(4) = mg_mesh%level(ilev+1)%kadj(3,jel(3))
      jel(5) = mg_mesh%level(ilev+1)%kadj(6,jel(1))
      jel(6) = mg_mesh%level(ilev+1)%kadj(6,jel(2))
      jel(7) = mg_mesh%level(ilev+1)%kadj(6,jel(3))
      jel(8) = mg_mesh%level(ilev+1)%kadj(6,jel(4))

      iMat = MaterialDistribution(nlmax-1)%x(jel(1))
!       write(*,*) iMAt
      do ii=1,8
       MaterialDistribution(nlmax)%x(jel(1)) = iMat
      end do
     end do
  end if

  ndof_nel = (knvt(NLMAX) + knat(NLMAX) + knet(NLMAX))
  allocate(dFrac(myMultiMat%nOfMaterials))

  DO ilev=NLMAX-1,NLMIN,-1
   ndof_nel = (knvt(ilev+1) + knat(ilev+1) + knet(ilev+1))
   do iel = 1,mg_mesh%level(ilev)%nel
     jel(1) = iel
     jel(2) = mg_mesh%level(ilev+1)%kadj(3,jel(1))
     jel(3) = mg_mesh%level(ilev+1)%kadj(3,jel(2))
     jel(4) = mg_mesh%level(ilev+1)%kadj(3,jel(3))
     jel(5) = mg_mesh%level(ilev+1)%kadj(6,jel(1))
     jel(6) = mg_mesh%level(ilev+1)%kadj(6,jel(2))
     jel(7) = mg_mesh%level(ilev+1)%kadj(6,jel(3))
     jel(8) = mg_mesh%level(ilev+1)%kadj(6,jel(4))

!      write(*,*) jel
     dFrac = 0d0
     do ii=1,8
      iMat = MaterialDistribution(ilev+1)%x(jel(ii))
      dFrac(iMat) = dFrac(iMat) + 1d0 !mg_mesh%level(ilev+1)%dvol(jel(ii))
     end do
     dMaxFrac = 0d0
!      write(*,*) dFrac
     do ii=1,myMultiMat%nOfMaterials
      IF (dFrac(ii).gt.dMaxFrac) THEN
       dMaxFrac = dFrac(ii)
       iMaxFrac = ii
      END IF
     end do

     MaterialDistribution(ilev)%x(jel(1)) = iMaxFrac

    end do
  END DO

  deallocate(dFrac)
 end if

else
 DO ilev=NLMAX,NLMIN,-1
  MaterialDistribution(ilev)%x = myMultiMat%InitMaterial
 end do
! write(*,*) MaterialDistribution(nlmax)%x
end if

! write(*,*) MaterialDistribution(NLMAX)%x(1:knel(NLMAX))
! pause
END SUBROUTINE UpdateMaterialProperties
