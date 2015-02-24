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
SUBROUTINE INIVCO()
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGALPHA,DGBETA_,DGGAMMA,DGINANG                                                                                  &
                         ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVPS,DGINVTH,DGINVVP,DGINVVS                                                  &
                         ,DGLAHEI,DGPVELO,DGSVELO,DGWEIPS                                                                                  &
!
                         ,IGNLAYE
!
IMPLICIT NONE
!
INTEGER :: J
!
!*******************************************************************************************************************************************
!
!$OMP PARALLEL DO
DO J = 1,IGNLAYE
   DGINVVS(J) = DGSVELO(J)
   DGINVVP(J) = DGPVELO(J)
   DGINVAL(J) = DGALPHA(J)
   DGINVBE(J) = DGBETA_(J)
   DGINVGA(J) = DGGAMMA(J)
   DGINVIN(J) = DGINANG(J)
   DGINVPS(J) = DGWEIPS(J)
   DGINVTH(J) = DGLAHEI(J)
ENDDO
!$OMP END PARALLEL DO
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE INIVCO
