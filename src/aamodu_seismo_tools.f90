!*******************************************************************************************************************************************
! PRELIMINARY SENSITIVITY ANALYSIS
! CAN BE PERFORMED BEFORE SOIL COLUMN INVERSION BY GENETIC ALGORITHM USING THOMSON-HASKELL MATRIX PROPAGATOR METHOD
!*******************************************************************************************************************************************
!
! COPYRIGHT   : - BRGM, FRANCE
!               - DPRI, JAPAN
!
! CONTRIBUTORS: - ARIANE DUCELLIER (a DOT ducellier AT brgm DOT fr)
!
! PURPOSE     : - THIS COMPUTER PROGRAM COMPUTES SOBOL COEFFICIENTS TO DETERMINE WHICH COLUMN PARAMETERS
!                 HAVE THE MOST INFLUENCE ON SPECTRAL RATIOS
!                 IT CAN BE USED BEFORE CARRYING OUT INVERSION BY GENETIC ALGORITHM OF 1-D SOIL COLUMN IN THE FREQUENCY DOMAIN
!                 USING THOMSON-HASKELL MATRIX PROPAGATOR METHOD
!
! VERSION     :  1.1
!
!*******************************************************************************************************************************************
!
MODULE AAMODU_SEISMO_TOOLS
!
INTERFACE SEISMO_TOOLS
   MODULE PROCEDURE SPEWIN
END INTERFACE SEISMO_TOOLS
!
CONTAINS
!
!*******************************************************************************************************************************************
! SMOOTH SPECTRA BY PARZEN'S SPECTRAL WINDOW
!*******************************************************************************************************************************************
   SUBROUTINE SPEWIN(X,DF,B)
!-------------------------------------------------------------------------------------------------------------------------------------------
! INPUT
!-------------------------------------------------------------------------------------------------------------------------------------------
!  X : WAVE IN FOURIER DOMAIN
!  B : SMOOTHING WINDOW
!  DF: FREQUENCY STEP DISCRETIZATION
!-------------------------------------------------------------------------------------------------------------------------------------------
! OUTPUT
!-------------------------------------------------------------------------------------------------------------------------------------------
!  X: SMOOTHED WAVE
!  B: SMOOTHING WINDOW
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   IMPLICIT NONE
!
   DOUBLE PRECISION, INTENT(INOUT), DIMENSION(:) :: X
   DOUBLE PRECISION, INTENT(INOUT)               :: B
   DOUBLE PRECISION, INTENT(IN)                  :: DF
   DOUBLE PRECISION, ALLOCATABLE  , DIMENSION(:) :: W,X1,X2 
   DOUBLE PRECISION                              :: BMAX,BMAX1,BMAX2,BMIN,DIF,S,UDF
!
   INTEGER                                       :: K,L,LE,LL,LMAX,LN,LT,N,ND
!
   N  = SIZE(X)
   ND = 2*(N-1)
   ALLOCATE(X1(ND+1),X2(ND+1))
!
   BMIN  = 560.0D0/151.0D0*DF
   BMAX1 = 140.0D0*DBLE(N-1)/151.0D0*DF
   BMAX2 = 14000.0D0/150.0D0*DF
   BMAX  = MIN(BMAX1,BMAX2)
!
   IF( (B.LT.BMIN) .OR. (B.GT.BMAX) ) THEN
      IF (ND.EQ.65536) THEN
         B = BMAX*1.0D0
      ELSEIF (ND.EQ.32768) THEN
         B = BMAX*0.5D0
      ELSEIF (ND.EQ.16384) THEN
         B = BMAX*0.25D0
      ELSEIF (ND.EQ.8192) THEN
         B = BMAX
      ELSE
         B = BMAX
      ENDIF
   ENDIF
!
   UDF = 1.854305D0/B*DF
!
   IF (UDF.GT.0.5D0) THEN
      WRITE(0,*) "ERROR IN SPEWIN: BANDWIDTH IS TOO NARROW"
      WRITE(0,'(E15.7,A,E15.7)') BMIN," < BAND WIDTH (HZ) < ",BMAX
   ENDIF
!
   LMAX = INT(2.0D0/UDF)+1
   IF (LMAX.GT.101) THEN
      WRITE(0,*) "ERROR IN SPEWIN: BANDWIDTH IS TOO WIDE"
      WRITE(0,'(E15.7,A,E15.7)') BMIN," < BAND WIDTH (HZ) < ",BMAX
   ENDIF
!
! SPECTRAL WINDOW
!
   ALLOCATE(W(LMAX))
   W(1) = 0.75D0*UDF
   DO L = 2,LMAX
      DIF  = 1.570796D0*DBLE(L-1)*UDF
      W(L) = W(1)*(SIN(DIF)/DIF)**4
   ENDDO
!
! SMOOTHING OF FOURIER SPECTRUM
!
   LL = LMAX*2-1
   LN =  LL-1    + N
   LT = (LL-1)*2 + N
   LE = LT-LMAX+1
!
   DO K = 1,LT
      X1(K) = 0.0D0
   ENDDO
   DO K = 1,N
      X1(LL-1+K) = X(K)
   ENDDO
!
   DO K = LMAX,LE
      S = W(1)*X1(K)
      DO L = 2,LMAX
         S = S + W(L)*(X1(K-L+1) + X1(K+L-1))
      ENDDO
      X2(K) = S
   ENDDO
!
   DO L = 2,LMAX
      X2(LL+L-1) = X2(LL+L-1) + X2(LL-L+1)
      X2(LN-L+1) = X2(LN-L+1) + X2(LN+L-1)
   ENDDO
!
   DO K = 1,N
      X(K) = X2(LL-1+K)
   ENDDO
!
   DEALLOCATE(X1,X2,W)
!
   RETURN
!*******************************************************************************************************************************************
   END SUBROUTINE SPEWIN
!*******************************************************************************************************************************************
!
END MODULE AAMODU_SEISMO_TOOLS
