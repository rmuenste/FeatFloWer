 REAL*8 function RotParabolicVelo2Dx(YR,ZR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,YR,ZR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((Y-YR)**2d0+(Z-ZR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dx = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dx = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 REAL*8 function RotParabolicVelo2Dy(XR,ZR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,ZR,XR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((X-XR)**2d0+(Z-ZR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dy = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dy = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 REAL*8 function RotParabolicVelo2Dz(XR,YR,DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux,YR,XR
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO)
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux
  DIST = SQRT((X-XR)**2d0+(Y-YR)**2d0)
 IF (DIST.LT.dR) THEN
  RotParabolicVelo2Dz = dScale*(DR - DIST)*(DR + DIST)
 ELSE
  RotParabolicVelo2Dz = 0d0
 END IF
 END 
!------------------------------------------------------------------------------
 function FlatVelo3D(DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR
 REAL*8  dNRM,dArea,dAvgVelo
 REAL*8, dimension(3) :: FlatVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 
  dArea = PI*DR**2d0
  dAvgVelo = dVolFlow/dArea

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)

  IF (DIST.LT.dR) THEN
   FlatVelo3D(:) = dNormal(:)*dAvgVelo
  ELSE
   FlatVelo3D(:) = 0d0
  END IF
  
  END 
!------------------------------------------------------------------------------
 function CurvedFlatVelo3D(DM,DRHO,dR2,dR3)
 REAL*8  dVolFlow,DM,DRHO,dR2,dR3
 REAL*8  dNRM,dAux1,dAux2,dAux3,dP1,dP2,dScale
 REAL*8, dimension(3) :: CurvedFlatVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)
  
  dP1 = dR3
  dP2 = (dR2 + dR3)/2d0
  
  dAux1 = -2d0*PI*( (dP1**4d0)/4d0 - ((dR2+dR3)*dP1**3d0)/3d0 + (dR2*dR3*dP1**2d0)/2d0 )
  dAux2 = -2d0*PI*( (dP2**4d0)/4d0 - ((dR2+dR3)*dP2**3d0)/3d0 + (dR2*dR3*dP2**2d0)/2d0 )
  dAux3 = 1d0*PI*( ((dR3-dR2)*(dR2+dR3)/4d0)**2d0)
  
  dScale = dVolFlow/(dAux1 - dAux2 + dAux3)

  IF (DIST.LT.(dR2+dR3)*0.5d0) THEN
   CurvedFlatVelo3D(:) =  dNormal(:)*dScale*0.25d0*(dR3-dR2)**2d0
  ELSE
   IF (DIST.LT.dR3) THEN
    CurvedFlatVelo3D(:) =  dNormal(:)*dScale*(DR3 - DIST)*(DIST - DR2)
   ELSE
    CurvedFlatVelo3D(:) = 0d0
   END IF
  END IF
  
  END 
!------------------------------------------------------------------------------
 function RotParabolicVelo3D(DM,DRHO,dR)
 REAL*8  dVolFlow,DM,DRHO,dR,daux
 REAL*8  dNRM
 REAL*8, dimension(3) :: RotParabolicVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 
  daux = (PI/2d0)*(DR**4d0)
  dScale = (dVolFlow/1d0)/daux

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)

  IF (DIST.LT.dR) THEN
   RotParabolicVelo3D(:) = dNormal(:)*dScale*(DR - DIST)*(DR + DIST)
  ELSE
   RotParabolicVelo3D(:) = 0d0
  END IF
  
  END 
!------------------------------------------------------------------------------
 function RotDoubleParabolicVelo3D(DM,DRHO,dR1,dR2)
 REAL*8  dVolFlow,DM,DRHO,dR1,dR2,daux
 REAL*8  dNRM
 REAL*8, dimension(3) :: RotDoubleParabolicVelo3D
 
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 
  daux = (PI/6d0)*(dR1+dR2)*((dR2-dR1)**3d0)
  dScale = (dVolFlow/1d0)/daux

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  dist = SQRT((X-(dCenter(1)))**2d0 + (Y-(dCenter(2)))**2d0 + (Z-(dCenter(3)))**2d0)

  IF (DIST.LT.dR2.and.DIST.GT.dR1) THEN
   RotDoubleParabolicVelo3D(:) = dNormal(:)*dScale*(DR2 - DIST)*(DIST - DR1)
  ELSE
   RotDoubleParabolicVelo3D(:) = 0d0
  END IF
  
  END 
!------------------------------------------------------------
 real function pulsativeProfile(simTime)
 implicit none
 real*8 :: simTime
 
  x_arr = (/0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)
  y_arr = (/0.8,3.0,3.6,3.3,2.2,2.0,1.5,1.3,1.1,0.9,0.8/)

  CC = (/0.6438, 1.0840, 1.1125, 0.8352, 0.6234, 0.5479, 0.4157, 0.3740, 0.3036, 0.2578, 0.0/)
  DD = (/7.9528, 0.8514, -0.2817, -5.2633, 1.0272, -2.5377, -0.1073, -0.7253, -0.6836, -0.2325, 0.0/)
  MM = (/0.0, -71.0134, -11.3310, -49.8165, 62.9046, -35.6481, 24.3032, -6.1801, 0.4173, 4.5111, 0.0/)

  normalizedTime = simTime - real(floor(simTime)) 

  j = 1
  do while (normalizedTime .gt. x_arr(j+1))
    j = j + 1
  end do

  h = x_arr(j+1) - x_arr(j)

  pulsativeProfile = CC(j) + DD(j) * (normalizedTime - 0.5 * (x_arr(j) + x_arr(j+1))) + &
                     1.0/(6.0 * h) * &
                     (MM(j+1) * (normalizedTime - x_arr(j))**3.0 - &
                     MM(j) * (normalizedTime - x_arr(j+1))**3.0) 
  
  end function pulsativeProfile
!------------------------------------------------------------
 function RectangleVelo3D(DM,DRHO,DA,DB)
 REAL*8  dVolFlow,DM,DRHO,dR2,dR3
 REAL*8  dNRM,dAux,dAux2,dAux3,dScale
 REAL*8  DA(3), DB(3), DP(3)
 REAL*8, dimension(3) :: RectangleVelo3D
  
  dP = (/X, Y, Z/)

  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  
  dAux2 = DOT_PRODUCT(DA-dCenter, DA-dCenter)**2-DOT_PRODUCT(dP-dCenter, DA-dCenter)**2
  dAux3 = DOT_PRODUCT(DB-dCenter, DB-dCenter)**2-DOT_PRODUCT(dP-dCenter, DB-dCenter)**2

  IF ( (dAux2.GE.0D0).and.(dAux3.GE.0D0) ) THEN
    dAux = dAux2*dAux3
    dScale = 9D0/16D0/ &
            DOT_PRODUCT(DA-dCenter, DA-dCenter)**2/ &
            DOT_PRODUCT(DB-dCenter, DB-dCenter)**2/ &
            NORM2(DA-dCenter)/NORM2(DB-dCenter)
    
    RectangleVelo3D(:) = dVolFlow*dScale*dAux*dNormal(:)
  ELSE
    RectangleVelo3D(:) = 0D0
  END IF


  end function RectangleVelo3D
!------------------------------------------------------------
 function CurvedRectangleVelo3D(DM,DRHO,DA,DB)
 REAL*8  dVolFlow,DM,DRHO
 REAL*8  dNRM,dScale
 REAL*8  DA(3), DB(3), DP(3)
 REAL*8  dVol, dAux, dR
 REAL*8, dimension(3) :: CurvedRectangleVelo3D
 REAL*8 dAdC, dBdC, dPdA, dPdB

  dP = (/X, Y, Z/)

  dScale = 8d0/3d0*NORM2(DB-dCenter)*NORM2(DA-dCenter)+ &
         PI/2*NORM2(DB-dCenter)**2
  
  dVolFlow = (1e3/3.6d3)*DM/(DRHO) 

  dNRM = SQRT(dNormal(1)*dNormal(1) +dNormal(2)*dNormal(2) + dNormal(3)*dNormal(3))
  dNormal(1) = dNormal(1)/dNRM
  dNormal(2) = dNormal(2)/dNRM
  dNormal(3) = dNormal(3)/dNRM
  
  dR = NORM2(DB-dCenter)

  ! some auxillary values
  dAdC = DOT_PRODUCT(dA-dCenter, DA-dCenter)
  dBdC = DOT_PRODUCT(dB-dCenter, DB-dCenter)

  dPdA = DOT_PRODUCT(dP-dCenter, DA-dCenter)
  dPdB = DOT_PRODUCT(dP-dCenter, DB-dCenter)

  IF ( (dPdA**2.LE.dAdC**2).and.(dPdB**2.LE.dBdC**2) ) THEN
      dAux = (dBdC**2-dPdB**2) / (dBdC**2)
  ELSEIF ( (dPdA.GE.dAdC).and.(NORM2(dP-DA).LE.dR) ) THEN
      dAux = (dR**2-NORM2(dP-DA)**2) / dR**2 
  ELSEIF ( (dPdA**2.GE.dAdC**2).and.(NORM2(dP-2*dCenter+DA).LE.dR) ) THEN
      dAux = (dR**2-NORM2(dP-2*dCenter+DA)**2) / dR**2
  ELSE
      dAux = 0D0
  END IF
  
  dScale = 1/dScale
  CurvedRectangleVelo3D(:) = dVolFlow*dScale*dAux*dNormal(:)


  end function CurvedRectangleVelo3D




