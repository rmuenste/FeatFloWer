      SUBROUTINE SlopeLimiterSub(dVal,dMid,dCorvg,kVert,NEL,NVT)
      REAL*8 dVal(4,*),dMid(3,*),dCorvg(3,*)
      INTEGER kVert(8,*),NVT,NEL,j
      REAL*8 U_min(NVT),U_max(NVT),dX(3),Alpha

      ! Collecting the maxima and minima
      U_min = +1d30
      U_max = -1d30
      DO IEL=1,NEL
       DO IVT=1,8
        iVert = kVert(ivt,iel)
        DD = dVal(1,iel)
        IF (DD.LT.U_min(iVert)) U_min(iVert) = DD
        IF (DD.GT.U_max(iVert)) U_max(iVert) = DD
       END DO
      END DO

      DO IEL=1,NEL
       Alpha = 1d0
       DO IVT=1,8
        iVert = kVert(ivt,iel)
!         IF (knpr(iVert).EQ.0) THEN
         dX(1) = dCorvg(1,iVert) - dMid(1,iel)
         dX(2) = dCorvg(2,iVert) - dMid(2,iel)
         dX(3) = dCorvg(3,iVert) - dMid(3,iel)
        ! dd = u_i - u_c
         DD = dVal(2,iel)*dX(1)
     *      + dVal(3,iel)*dX(2) + dVal(4,iel)*dX(3)
         IF (DD.GT.0d0) THEN
          Alpha = MIN(Alpha,(U_max(iVert)-dVal(1,iel))/DD)
         END IF
         IF (DD.LT.0d0) THEN
          Alpha = MIN(Alpha,(U_min(iVert)-dVal(1,iel))/DD)
         END IF
!         END IF
       END DO
       IF (Alpha.NE.1d0) THEN
        dVal(2,iel) = Alpha*dVal(2,iel)
        dVal(3,iel) = Alpha*dVal(3,iel)
        dVal(4,iel) = Alpha*dVal(4,iel)
!         write(*,'(8G12.4)') iel,alpha,(dVal(j,iel),j=1,4)
       END IF
      END DO

!       pause
      END
