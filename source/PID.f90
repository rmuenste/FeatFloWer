SUBROUTINE PID_controller(T,dt,myPID)
USE types, only : tPID

IMPLICIT NONE

TYPE(tPID) myPID
REAL*8 T,dt
!--------------------------------------------------
REAL*8 e

 e = myPID%T_set  - T  
 myPID%sumI = myPID%sumI + e*dt
 
 myPID%P = myPID%omega_P * e
 myPID%I = (myPID%omega_P/myPID%omega_I) * myPID%sumI
 myPID%D = (myPID%omega_P*myPID%omega_D) * (e-myPID%e_old)/dt
 
 myPID%e_old = e
 myPID%PID = myPID%P + myPID%I + myPID%D

END SUBROUTINE PID_controller
