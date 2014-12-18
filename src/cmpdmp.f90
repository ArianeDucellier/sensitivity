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
SUBROUTINE CMPDMP(OM,XA,XB,XG,VS,DAMP,INDD)
!
!*******************************************************************************************************************************************
!
IMPLICIT NONE
!
DOUBLE PRECISION, INTENT(INOUT) :: DAMP
DOUBLE PRECISION, INTENT(IN)    :: OM,VS,XA,XB,XG
DOUBLE PRECISION                :: F,PI
!
INTEGER         , INTENT(IN)    :: INDD
!
!*******************************************************************************************************************************************
!
PI = ACOS(-1.0D0)
F  = OM/(2.0D0*PI)
!
IF (INDD.EQ.0) THEN
   DAMP = XA/OM+XB
ELSEIF (INDD.EQ.1) THEN
   DAMP = XA/100.0D0
ELSEIF (INDD.EQ.2) THEN
   DAMP = XA/100.D0*F
ELSEIF (INDD.EQ.3) THEN
   DAMP = XA/100.0D0*F**(XB)
ELSEIF (INDD.EQ.5) THEN
   DAMP = XA/100.0D0*F**(XB)+XG/100.0D0
ELSEIF (INDD.EQ.6) THEN
   DAMP = XB/(2.0D0*VS)*F**(XA)
ELSE
   WRITE(*,'(A)') "ERROR: DAMPING TYPE NOT RECOGNIZED"
   STOP
ENDIF
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE CMPDMP
