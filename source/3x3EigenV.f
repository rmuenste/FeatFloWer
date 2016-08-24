* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVD3(A, Q, W)
* ----------------------------------------------------------------------------
* Calculates the eigenvalues and normalized eigenvectors of a symmetric
* 3x3 matrix A using Cuppen's Divide & Conquer algorithm.
* The function accesses only the diagonal and upper triangular parts of
* A. The access is read-only. 
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
* Dependencies:
*   DSYEV2(), SLVSEC3(), DSYTRD3()
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

*     .. Local Variables ..
      DOUBLE PRECISION R(3,3)
      DOUBLE PRECISION P(3,3)
      DOUBLE PRECISION D(3), E(2), Z(3)
      DOUBLE PRECISION C, S, T
      INTEGER          I, J, K

*     .. External Functions ..
      EXTERNAL         DSYEV2, SLVSEC3, DSYTRD3
      
*     Initialize Q
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 I = 1, N
        DO 11 J = 1, N
          Q(I, J) = 0.0D0
   11   CONTINUE
   10 CONTINUE
  
*     Transform A to real tridiagonal form by the Householder method
      CALL DSYTRD3(A, R, W, E)

      
*     "Divide"
*     -------------------------------
  
*     Detect matrices that factorize to avoid multiple eigenvalues in the Divide/Conquer algorithm
      DO 20 I = 1, N-1
        T = ABS(W(I)) + ABS(W(I+1))
        IF (ABS(E(I)) .LE. 8.0D0 * EPS * T) THEN
          IF (I .EQ. 1) THEN
            CALL DSYEV2(W(2), E(2), W(3), D(2), D(3), C, S)
            W(2) = D(2)
            W(3) = D(3)
*           --- The rest of this IF-branch can be omitted if only the eigenvalues are desired ---
            Q(1, 1) = 1.0D0
            DO 30 J = 2, N
              Q(J, 2) = S * R(J, 3) + C * R(J, 2)
              Q(J, 3) = C * R(J, 3) - S * R(J, 2)
   30       CONTINUE
          ELSE
            CALL DSYEV2(W(1), E(1), W(2), D(1), D(2), C, S)
            W(1)    = D(1)
            W(2)    = D(2)
*           --- The rest of this ELSE-branch can be omitted if only the eigenvalues are desired ---
            Q(1, 1) = C
            Q(1, 2) = -S
            Q(2, 1) = R(2, 2) * S
            Q(2, 2) = R(2, 2) * C
            Q(2, 3) = R(2, 3)
            Q(3, 1) = R(3, 2) * S
            Q(3, 2) = R(3, 2) * C
            Q(3, 3) = R(3, 3)
          END IF
          RETURN
        END IF
   20 CONTINUE

*     Calculate eigenvalues and eigenvectors of 2x2 block
      CALL DSYEV2(W(2) - E(1), E(2), W(3), D(2), D(3), C, S)
      D(1) = W(1) - E(1)

      
*     "Conquer"
*     -------------------------------

*     Determine coefficients of secular equation
      Z(1) = E(1)
      Z(2) = E(1) * C**2
      Z(3) = E(1) * S**2

*     Call SLVSEC3 with D sorted in ascending order. We make use of the
*     fact that DSYEV2 guarantees D[1] >= D[2].
      IF (D(1) .LT. D(3)) THEN
        CALL SLVSEC3(D, Z, W, P, 1, 3, 2)
      ELSE IF (D(1) .LT. D(2)) THEN
        CALL SLVSEC3(D, Z, W, P, 3, 1, 2)
      ELSE
        CALL SLVSEC3(D, Z, W, P, 3, 2, 1)
      END IF

*     --- The rest of this subroutine can be omitted if only the eigenvalues are desired ---

*     Calculate eigenvectors of matrix D + BETA * Z * Z^T and store them
*     in the columns of P
      Z(1) = SQRT(ABS(E(1)))
      Z(2) = C * Z(1)
      Z(3) = -S * Z(1)

*     Detect duplicate elements in D to avoid division by zero
      T = 8.0D0 * EPS * ( ABS(D(1)) + ABS(D(2)) + ABS(D(3)) )
      IF (ABS(D(2) - D(1)) .LE. T) THEN
        DO 40 J = 1, N
          IF (P(1, J) * P(2, J) .LE. 0.0D0) THEN
            P(1, J) = Z(2)
            P(2, J) = -Z(1)
            P(3, J) = 0.0D0
          ELSE
            DO 45 I = 1, N
              P(I, J) = Z(I) / P(I, J)
   45       CONTINUE
          END IF
   40   CONTINUE
      ELSE IF (ABS(D(3) - D(1)) .LE. T) THEN
        DO 50 J = 1, N
          IF (P(1, J) * P(3, J) .LE. 0.0D0) THEN
            P(1, J) = Z(3)
            P(2, J) = 0.0D0
            P(3, J) = -Z(1)
          ELSE
            DO 55 I = 1, N
              P(I, J) = Z(I) / P(I, J)
   55       CONTINUE
          END IF
   50   CONTINUE
      ELSE
        DO 60 J = 1, N
          DO 61 I = 1, N
            IF (P(I, J) .EQ. 0.0) THEN
              P(I, J) = 1.0
              P(1 + MOD(I, N), J)   = 0.0
              P(1 + MOD(I+1, N), J) = 0.0
              GO TO 60
            ELSE
              P(I, J) = Z(I) / P(I, J)
            END IF
   61     CONTINUE
   60   CONTINUE
      END IF

*     Normalize eigenvectors of D + BETA * Z * Z^T
      DO 70 J = 1, N
        T = P(1, J)**2 + P(2, J)**2 + P(3, J)**2
        T = 1.0D0 / SQRT(T)
        DO 75 I = 1, N
          P(I, J) = P(I, J) * T
   75   CONTINUE
   70 CONTINUE
  
*     Undo diagonalization of 2x2 block
      DO 80 J = 1, N
        T       = P(2, J)
        P(2, J) = C * T - S * P(3, J)
        P(3, J) = S * T + C * P(3, J)
   80 CONTINUE

*     Undo Householder transformation
      DO 90 J = 1, N
        DO 91 K = 1, N
          T       = P(K, J)
          DO 95 I = 1, N
            Q(I, J) = Q(I, J) + T * R(I, K)
   95     CONTINUE
   91   CONTINUE
   90 CONTINUE

      END SUBROUTINE
* End of subroutine DSYEVD3

C
C
C

* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYTRD3(A, Q, D, E)
* ----------------------------------------------------------------------------
* Reduces a symmetric 3x3 matrix to real tridiagonal form by applying
* (unitary) Householder transformations:
*            [ D[1]  E[1]       ]
*    A = Q . [ E[1]  D[2]  E[2] ] . Q^T
*            [       E[2]  D[3] ]
* The function accesses only the diagonal and upper triangular parts of
* A. The access is read-only.
* ---------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION D(3)
      DOUBLE PRECISION E(2)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )

*     .. Local Variables ..
      DOUBLE PRECISION U(N), P(N)
      DOUBLE PRECISION OMEGA, F
      DOUBLE PRECISION K, H, G
      INTEGER          I, J

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 I = 1, N
        Q(I,I) = 1.0D0
        DO 11, J = 1, I-1
          Q(I, J) = 0.0D0
          Q(J, I) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Bring first row and column to the desired form
      H = A(1,2)**2 + A(1,3)**2
      IF (A(1,2) .GT. 0.0D0) THEN
        G = -SQRT(H)
      ELSE
        G = SQRT(H)
      END IF
      E(1)  = G
      F     = G * A(1,2)
      U(2)  = A(1,2) - G
      U(3)  = A(1,3)

      OMEGA = H - F
      IF (OMEGA > 0.0D0) THEN
        OMEGA = 1.0D0 / OMEGA
        K     = 0.0D0
        DO 20 I = 2, N
          F    = A(2,I)*U(2) + A(I,3)*U(3)
          P(I) = OMEGA * F
          K    = K + U(I) * F
  20    CONTINUE
        K = 0.5D0 * K * OMEGA**2

        DO 30 I = 2, N
          P(I) = P(I) - K * U(I)
  30    CONTINUE

        D(1) = A(1,1)
        D(2) = A(2,2) - 2.0D0 * P(2) * U(2)
        D(3) = A(3,3) - 2.0D0 * P(3) * U(3)

*       Store inverse Householder transformation in Q
*       --- This loop can be omitted if only the eigenvalues are desired ---
        DO 40, J = 2, N
          F = OMEGA * U(J)
          DO 41 I = 2, N
            Q(I,J) = Q(I,J) - F * U(I)
   41     CONTINUE
   40   CONTINUE
            
*       Calculated updated A(2, 3) and store it in E(2)
        E(2) = A(2, 3) - P(2) * U(3) - U(2) * P(3)
      ELSE
        DO 50 I = 1, N
          D(I) = A(I, I)
  50    CONTINUE
        E(2) = A(2, 3)
      END IF
      
      END SUBROUTINE
* End of subroutine DSYTRD3

C
C
C

* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE DSYEV2(A, B, C, RT1, RT2, CS, SN)
* ----------------------------------------------------------------------------
* Calculates the eigensystem of a real symmetric 2x2 matrix
*    [ A  B ]
*    [ B  C ]
* in the form
*    [ A  B ]  =  [ CS  -SN ] [ RT1   0  ] [  CS  SN ]
*    [ B  C ]     [ SN   CS ] [  0   RT2 ] [ -SN  CS ]
* where RT1 >= RT2. Note that this convention is different from the one used
* in the LAPACK routine DLAEV2, where |RT1| >= |RT2|.
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A, B, C, RT1, RT2, CS, SN

*     .. Local Variables ..
      DOUBLE PRECISION SM, DF, RT, T

      SM = A + C
      DF = A - C
      RT = SQRT(DF**2 + 4.0D0 * B**2)

*     Calculate eigenvalues
      IF (SM .GT. 0.0D0) THEN
        RT1 = 0.5D0 * (SM + RT)
        T   = 1.0D0 / RT1
        RT2 = (A*T)*C - (B*T)*B
      ELSE IF (SM .LT. 0.0D0) THEN
        RT2 = 0.5D0 * (SM - RT)
        T   = 1.0D0 / RT2
        RT1 = (A*T)*C - (B*T)*B
      ELSE
*       This case needs to be treated separately to avoid DIV by 0
        RT1 = 0.5D0 * RT
        RT2 = -0.5D0 * RT
      END IF
      
*     Calculate eigenvectors
      IF (DF .GT. 0.0D0) THEN
        CS = DF + RT
      ELSE
        CS = DF - RT
      END IF

      IF (ABS(CS) .GT. 2.0D0 * ABS(B)) THEN
        T  = -2.0D0 * B / CS
        SN = 1.0D0 / SQRT(1.0D0 + T**2)
        CS = T * SN
      ELSE IF (ABS(B) .EQ. 0.0D0) THEN
        CS = 1.0D0
        SN = 0.0D0
      ELSE
        T  = -0.5D0 * CS / B
        CS = 1.0D0 / SQRT(1.0D0 + T**2)
        SN = T * CS
      END IF

      IF (DF .GT. 0.0D0) THEN
        T  = CS
        CS = -SN
        SN = T
      END IF
      
      END SUBROUTINE
* End of subroutine DSYEV2

C
C
C

* ----------------------------------------------------------------------------
* Numerical diagonalization of 3x3 matrcies
* Copyright (C) 2006  Joachim Kopp
* ----------------------------------------------------------------------------
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with this library; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
* ----------------------------------------------------------------------------


* ----------------------------------------------------------------------------
      SUBROUTINE SLVSEC3(D, Z, W, R, I1, I2, I3)
* ----------------------------------------------------------------------------
* Finds the three roots lambda_j of the secular equation
*   f(W_j) = 1 + Sum[ Z_i / (D_i - W_j) ]  ==  0.
* It is assumed that D_0 <= D_1 <= D_2, and that all Z_i have the same sign.
* The arrays R will contain the information required for the calculation
* of the eigenvectors:
*   R_ij = D_i - W_j.
* These differences can be obtained with better accuracy from intermediate
* results.
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION D(3), Z(3), W(3)
      DOUBLE PRECISION R(3,3)
      INTEGER          I1, I2, I3

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
      DOUBLE PRECISION SQRT3
      PARAMETER        ( SQRT3 = 1.73205080756887729352744634151D0 )
      DOUBLE PRECISION EPS
      PARAMETER        ( EPS = 2.2204460492503131D-16 )

*     .. Local Variables ..
      DOUBLE PRECISION A(4)
      DOUBLE PRECISION DELTA
      DOUBLE PRECISION DD(3)
      DOUBLE PRECISION XL, XH, X, X0(3)
      DOUBLE PRECISION DX, DXOLD
      DOUBLE PRECISION F, DF
      DOUBLE PRECISION ERROR
      DOUBLE PRECISION T(3)
      DOUBLE PRECISION ALPHA, BETA, GAMMA
      DOUBLE PRECISION P, SQRTP, Q, C, S, PHI
      INTEGER          I, J, NITER

*     Determine intervals which must contain the roots
      IF (Z(1) .GT. 0.0D0) THEN
        A(1) = D(I1)
        A(2) = D(I2)
        A(3) = D(I3)
        A(4) = ABS(D(1) + 3.0D0*Z(1))
     $           + ABS(D(2) + 3.0D0*Z(2))
     &           + ABS(D(3) + 3.0D0*Z(3))
      ELSE
        A(1) = -ABS(D(1) + 3.0D0*Z(1))
     $           - ABS(D(2) + 3.0D0*Z(2))
     &           - ABS(D(3) + 3.0D0*Z(3))
        A(2) = D(I1)
        A(3) = D(I2)
        A(4) = D(I3)
      END IF

*     Calculate roots of f(x) = 0 analytically (analogous to ZHEEVC3)
      T(1)  = D(2) * D(3)
      T(2)  = D(1) * D(3)
      T(3)  = D(1) * D(2)
      GAMMA = T(1) * D(1) + ( Z(1)*T(1) + Z(2)*T(2) + Z(3)*T(3) )
      BETA  = ( Z(1) * (D(2) + D(3) ) + Z(2) * (D(1) + D(3) )
     $            + Z(3) * (D(1) + D(2)) )
     $            + ( T(1) + T(2) + T(3) )
      ALPHA = ( Z(1) + Z(2) + Z(3) ) + ( D(1) + D(2) + D(3) )

      P     = ALPHA**2 - 3.0D0 * BETA
      Q     = ALPHA * (P - (3.0D0/2.0D0) * BETA) + (27.0D0/2.0D0)*GAMMA
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0D0 * ( 0.25D0 * BETA**2 * (P - BETA)
     $                       - GAMMA * (Q - (27.0D0/4.0D0) * GAMMA) )
      PHI   = (1.0D0/3.0D0) * ATAN2(SQRT(ABS(PHI)), Q)
      C     = SQRTP * COS(PHI)
      S     = (1.0D0/SQRT3) * SQRTP * ABS(SIN(PHI))

*     Make sure the roots are in ascending order
      X0(1) = (1.0D0/3.0D0) * (ALPHA - C)
      X0(2) = X0(1)
      X0(3) = X0(1)
      IF (C .GT. S) THEN
        X0(1) = X0(1) - S
        X0(2) = X0(2) + S
        X0(3) = X0(3) + C
      ELSE IF (C .LT. -S) THEN
        X0(1) = X0(1) + C
        X0(2) = X0(2) - S
        X0(3) = X0(3) + S
      ELSE
        X0(1) = X0(1) - S
        X0(2) = X0(2) + C
        X0(3) = X0(3) + S
      END IF

*     Refine roots with a combined Bisection/Newton-Raphson method
      DO 10 I = 1, N
        XL    = A(I)
        XH    = A(I+1)
        DX    = 0.5D0 * (XH - XL)
        DXOLD = DX

*       Make sure that XL != XH
        IF (DX .EQ. 0.0D0) THEN
          W(I) = XL
          DO 15 J = 1, N
            R(J, I) = D(J) - XL
   15     CONTINUE
          GO TO 10
        END IF

*       Shift the root close to zero to achieve better accuracy
        IF (X0(I) .GE. XH) THEN
          DELTA = XH
          X     = -DX
          DO 20 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J) - X
   20     CONTINUE
        ELSE IF (X0(I) .LE. XL) THEN
          DELTA = XL
          X     = DX
          DO 30 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J) - X
   30     CONTINUE
        ELSE
          DELTA = X0(I)
          X     = 0.0D0
          DO 40 J = 1, N
            DD(J)   = D(J) - DELTA
            R(J, I) = DD(J)
   40     CONTINUE
        END IF
        XL = XL - DELTA
        XH = XH - DELTA

*       Make sure that f(XL) < 0 and f(XH) > 0
        IF (Z(1) .LT. 0.0D0) THEN
          P  = XH
          XH = XL
          XL = P
        END IF

*       Main iteration loop
        DO 50 NITER = 1, 500
*         Evaluate f and f', and calculate an error estimate
          F     = 1.0D0
          DF    = 0.0D0
          ERROR = 1.0D0
          DO 60 J = 1, N
            T(1)  = 1.0D0 * (1.0 / R(J, I))
            T(2)  = Z(J) * T(1)
            T(3)  = T(2) * T(1)
            F     = F + T(2)
            ERROR = ERROR + ABS(T(2))
            DF    = DF + T(3)
   60     CONTINUE

*         Check for convergence and leave loop if applicable
          IF (ABS(F) .LE. EPS * (8.0D0 * ERROR + ABS(X * DF))) THEN
            GO TO 70
          END IF

*         Adjust interval boundaries
          IF (F .LT. 0.0D0) THEN
            XL   = X
          ELSE
            XH   = X
          END IF

*         Check, whether Newton-Raphson would converge fast enough.
*         If so, give it a try. If not, or if it would run out of
*         bounds, use bisection.
          IF ( ABS(2.0D0 * F) .LT. ABS(DXOLD * DF) ) THEN
            DXOLD = DX
            DX    = F * (1.0 / DF)
            X     = X - DX
            IF ( (X - XH) * (X - XL) .GE. 0.0D0 ) THEN
              DX = 0.5D0 * (XH - XL)
              X  = XL + DX
            END IF
          ELSE
            DX = 0.5D0 * (XH - XL)
            X  = XL + DX
          END IF
            
*         Prepare next iteration
          DO 80 J = 1, N
            R(J, I) = DD(J) - X
   80     CONTINUE
   50   CONTINUE

*       Un-shift result
   70   W(I) = X + DELTA
   10 CONTINUE

      END SUBROUTINE
* End of subroutine SLVSEC3

 

