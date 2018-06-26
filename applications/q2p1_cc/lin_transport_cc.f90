module transport_linear_cc
use Transport_Q1
use def_LinScalar

contains
!
! ----------------------------------------------
!
SUBROUTINE Init_LinScalar_cc

NLMAX = NLMAX + 1

IF (myid.ne.master) THEN

 ! Building up the matrix strucrures
 CALL Create_MatStruct()

 ! Building up the matrix strucrures
 CALL Create_AFCStruct()

 ! Building up linear operators
 ! Mass matrix
 CALL Create_MassMat()

! Mass matrix
 CALL Create_LMassMat()

! Convection matrix (only allocation)
 CALL Create_LKonvMat()

! Diffusion matrix 
 CALL Create_DiffMat(mgDiffCoeff(NLMAX)%x)

! Iteration matrix (only allocation)
 CALL Create_AMat()

END IF

! Initialize the scalar quantity
 CALL Initialize(Tracer)

! Set the types of boundary conditions (set up knpr)
 CALL Create_Knpr(LinSc_Knpr)

Tracer%cName = "Tracer"
Tracer%prm%SolvIter = 1
Tracer%prm%AFC = .TRUE.
IF (Tracer%prm%AFC) THEN
 Tracer%prm%NLmin = 2
ELSE
 Tracer%prm%NLmin = 2
END IF
Tracer%prm%NLmax   =20
Tracer%prm%defCrit =1d-6
Tracer%prm%epsCrit =1d-3
Tracer%prm%MinDef  =1d-16
Tracer%prm%SolvType=1

NLMAX = NLMAX - 1

END SUBROUTINE Init_LinScalar_cc
!
! ----------------------------------------------
!
SUBROUTINE InitCond_LinScalar_cc()
USE PP3D_MPI, ONLY : myid,master,showid,myMPI_Barrier

NLMAX = NLMAX + 1

ILEV=NLMAX
CALL SETLEV(2)

CALL LinSc_InitCond(mg_mesh%level(ilev)%dcorvg)

! Set boundary conditions
CALL Boundary_LinSc_Val()

NLMAX = NLMAX - 1

END SUBROUTINE InitCond_LinScalar_cc

end module transport_linear_cc
