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
SUBROUTINE READPF(C10)
!
!*******************************************************************************************************************************************
!
IMPLICIT NONE
!
INTEGER       :: IOS99
!
CHARACTER(10) :: C10
!
!*******************************************************************************************************************************************
!
OPEN(UNIT=99,FILE="control_file_genetic_algo",STATUS="OLD",IOSTAT=IOS99)
!
IF (IOS99.NE.0) THEN
   WRITE(*,'(A)') "FILE control_file_genetic_algo NOT FOUND"
   STOP
ELSE
   READ(UNIT=99,FMT=*) C10
ENDIF
!
CLOSE(99)
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE READPF
