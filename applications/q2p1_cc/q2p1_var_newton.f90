MODULE var_QuadScalar_newton

USE var_QuadScalar

IMPLICIT NONE

!Newton-Preconditioner
! REAL*8  , DIMENSION(:)  , POINTER :: barM11mat,barM22mat,barM33mat,barM12mat,barM13mat,barM23mat,barM21mat,barM31mat,barM32mat
REAL*8  , DIMENSION(:)  , POINTER :: AA11mat,AA22mat,AA33mat,AA12mat,AA13mat,AA23mat,AA21mat,AA31mat,AA32mat

! TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM11mat,mg_barM22mat,mg_barM33mat
! TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM12mat,mg_barM13mat,mg_barM23mat
! TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_barM21mat,mg_barM31mat,mg_barM32mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_AA11mat,mg_AA22mat,mg_AA33mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_AA12mat,mg_AA13mat,mg_AA23mat
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_AA21mat,mg_AA31mat,mg_AA32mat


! (possible) Force calculation
REAL*8  , DIMENSION(:)  , POINTER :: BXMat_new,BYMat_new,BZMat_new
TYPE (mg_Matrix), DIMENSION(:)  , ALLOCATABLE , TARGET :: mg_BXMat_new,mg_BYMat_new,mg_BZMat_new

TYPE tParamCC
 INTEGER NLmin,NLmax
 Real*8 :: Alpha,StoppingCriterion,ValAdap(2)
 INTEGER MinLev,MedLev,MinIterCycle,MaxIterCycle,nSmootherSteps
 REAL*8 Criterion,RLX 
 CHARACTER*1 CycleType
 integer :: vanka,BDF
END TYPE tParamCC

TYPE(tParamCC) ccParams

! needed for BDF(2) or BDF(3), better look than playing around with tstep
REAL*8 :: zeitstep
INTEGER :: tsm

END MODULE var_QuadScalar_newton
