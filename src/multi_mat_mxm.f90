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
SUBROUTINE MULTI_MAT_MXM(Z1,Z2,M)
!
!*******************************************************************************************************************************************
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: M
INTEGER             :: I,J,K
!
DOUBLE COMPLEX      :: Z1(M,M),Z2(M,M),ZWK(M,M)
!
!*******************************************************************************************************************************************
!
DO I = 1,M
   DO J = 1,M
      ZWK(J,I) = CMPLX(0.0D0,0.0D0)
      DO K = 1,M
         ZWK(J,I) = ZWK(J,I) + Z1(J,K)*Z2(K,I)
      ENDDO
   ENDDO
ENDDO
!
DO I = 1,M
   DO J = 1,M
      Z1(J,I) = ZWK(J,I)
   ENDDO
ENDDO
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE MULTI_MAT_MXM
