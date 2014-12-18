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
SUBROUTINE CHKOPN(F,IOS)
!
!*******************************************************************************************************************************************
!
IMPLICIT NONE
!
INTEGER      , INTENT(IN) :: IOS
!
CHARACTER(92), INTENT(IN) :: F
!
!*******************************************************************************************************************************************
!
IF (IOS.NE.0) THEN
   WRITE(*,'(3A)') "FILE ",TRIM(F)," NOT FOUND"
   STOP
ENDIF
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE CHKOPN
