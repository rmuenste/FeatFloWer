SUBROUTINE MeltFunction_MF(MF,T)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: MF
real*8, intent (in) :: T
!---------------------------------------------------
REAL*8 :: TM=403.15d0-273.15d0,TS=0.44d0,MS=0.15d0

 MF = (0.5d0*((tanh((T-TM)*TS)) + 1d0))**MS

END SUBROUTINE MeltFunction_MF
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Rho(Rho,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Rho
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: mV_s=0.0000005d0, nV_s=0.000887d0
REAL*8 :: mV_l=0.0000010d0, nV_l=0.000865d0
REAL*8 :: ES=1.1d0
REAL*8 :: V_s,V_l,V

 V_s = mV_s*(T - 273.15d0) + nV_s
 V_l = mV_l*(T - 273.15d0) + nV_l
 
 V = V_s*(1d0 - MF**(ES)) + V_s*(MF)**(ES)
 
 Rho = 1d0/V ! kg/m3 

END SUBROUTINE MeltFunction_Rho
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Cp(Cp,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Cp
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: Cp_s=2.5d3, nH_s=-738.5d3
REAL*8 :: Cp_l=3.4d3, nH_l=-912.0d3

 ! (Cp_s*T + n_s)*(1-MF) + (Cp_l*T + n_l)*MF =
 ! (1-MF)*Cp_s*T + (1-MF)*n_s + MF*Cp_l*T + MF*n_l = 
 ! [(1-MF)*Cp_s + MF*Cp_l]*T + [(1-MF)*n_s + MF*n_l] =
 ! Cp*T + H_ref
 
 Cp = (1d0-MF)*Cp_s + MF*Cp_l ! J/g/K

END SUBROUTINE MeltFunction_Cp
!
! ----------------------------------------------
!
SUBROUTINE MeltFunction_Lambda(Lambda,T,MF)
USE Transport_Q2P1, ONLY : Properties
USE Sigma_User, ONLY: myMultiMat,tRheology
real*8, intent(out) :: Lambda
real*8, intent (in) :: T,MF
!---------------------------------------------------
REAL*8 :: HS=0.43 !W/(m*K)
REAL*8 :: m=0.0 !W/(m*K*K)
REAL*8 :: n=0.258 !W/(m*K)
REAL*8 :: T_0=273.15 !K
REAL*8 :: WS=0.000022 !1/(K*K)
REAL*8 :: e=2.7182818284d0
REAL*8 :: Lambda_l,Lambda_s

 Lambda_s = HS*(e**(-WS*((T-273.15)-T_0)**2d0))
 Lambda_l = m*(T-273.15) + n
 
 Lambda = (1d0-MF)*Lambda_s + MF*Lambda_l ! W/(m*K)

END SUBROUTINE MeltFunction_Lambda
