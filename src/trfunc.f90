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
SUBROUTINE TRFUNC(ILINDIV)
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGINVIN,DGINVNU,DGINVPS,DGINVVP,DGINVVS                                                                          &
                         ,DGRTHVB,DGRTHVT,DGRTPPX,DGRTPPZ,DGRTPSX,DGRTPSZ,DGRTSHY,DGRTSVX,DGRTSVZ                                          &
                         ,DGTRPZB,DGTRPZT,DGTRSHB,DGTRSHT                                                                                  &
!
                         ,IGNFOLD,IGNLAYE                                                                                                  &
!
                         ,CGIDPPX,CGIDPPZ,CGIDPSX,CGIDPSZ,CGIDSHY,CGIDSVX,CGIDSVZ,CGTYPIN
!
IMPLICIT NONE
!
DOUBLE PRECISION    :: WEIGHTPS
!
INTEGER, INTENT(IN) :: ILINDIV
INTEGER             :: I,ILACTIN
!
!*******************************************************************************************************************************************
! UPDATE OF POISSON RATIO
!*******************************************************************************************************************************************
!
DO I = 1,IGNLAYE
   DGINVNU(I) =(DGINVVP(I)**2-2.0D0*DGINVVS(I)**2)/(2.0D0*(DGINVVP(I)**2-DGINVVS(I)**2))
ENDDO
!
!*******************************************************************************************************************************************
! CALL THOMSON-HASKELL PROPAGATOR MATRIX
!*******************************************************************************************************************************************
!
ILACTIN = 0
IF (DGINVIN(IGNLAYE).EQ.0.0D0) THEN
   ILACTIN = 1
   DGINVIN(IGNLAYE) = 1.0D-5 !USED TO AVOID HAVING PROBLEMS WITH PPWAVE SUBROUTINE. THIS HAS NO INFLUENCE ON THE RESULTS BECAUSE THE ANGLE IS VERY SMALL
ENDIF
!
IF (CGIDSHY.EQ." SH") THEN
   CALL SHWAVE()
ENDIF
IF ( (CGIDSVX.EQ."SVX").OR.(CGIDSVZ.EQ."SVZ").OR.(CGIDPSX.EQ."PSX").OR.(CGIDPSZ.EQ."PSZ") ) THEN
   CALL SVWAVE()
ENDIF
IF ( (CGIDPPX .EQ." PX").OR.(CGIDPPZ .EQ." PZ").OR.(CGIDPSX.EQ."PSX").OR.(CGIDPSZ.EQ."PSZ") ) THEN
   CALL PPWAVE()
ENDIF
!
IF (ILACTIN.EQ.1) DGINVIN(IGNLAYE) = 0.0D0
!
!*******************************************************************************************************************************************
! MERGE SVX AND PX SPECTRAL RATIOS AND SVZ AND PZ SPECTRAL RATIOS TO CREATE PSX AND PSZ
!*******************************************************************************************************************************************
!
WEIGHTPS = DGINVPS(IGNLAYE)
!
IF ( (CGTYPIN.EQ."SPR") .AND. ((CGIDPSX.EQ."PSX") .OR. (CGIDPSZ.EQ."PSZ")) ) THEN
   IF (DGINVIN(IGNLAYE).GT.1.0D-4) THEN
      IF ( (WEIGHTPS.GT.0.0D0) .AND. (WEIGHTPS.LT.1.0D0) ) THEN
         IF (CGIDPSX.EQ."PSX") THEN
            DGRTPSX(1:IGNFOLD-1) = WEIGHTPS*DGRTSVX(1:IGNFOLD-1) + (1.0D0-WEIGHTPS)*DGRTPPX(1:IGNFOLD-1)
         ENDIF
         IF (CGIDPSZ.EQ."PSZ") THEN
            DGRTPSZ(1:IGNFOLD-1) = WEIGHTPS*DGRTSVZ(1:IGNFOLD-1) + (1.0D0-WEIGHTPS)*DGRTPPZ(1:IGNFOLD-1)
         ENDIF
      ELSEIF (WEIGHTPS.EQ.0.0D0) THEN
         IF (CGIDPSX.EQ."PSX") THEN
            DGRTPSX(1:IGNFOLD-1) = DGRTPPX(1:IGNFOLD-1)
         ENDIF
         IF (CGIDPSZ.EQ."PSZ") THEN
            DGRTPSZ(1:IGNFOLD-1) = DGRTPPZ(1:IGNFOLD-1)
         ENDIF
      ELSEIF (WEIGHTPS.EQ.1.0D0) THEN
         IF (CGIDPSX.EQ."PSX") THEN
            DGRTPSX(1:IGNFOLD-1) = DGRTSVX(1:IGNFOLD-1)
         ENDIF
         IF (CGIDPSZ.EQ."PSZ") THEN
            DGRTPSZ(1:IGNFOLD-1) = DGRTSVZ(1:IGNFOLD-1)
         ENDIF
      ENDIF
   ELSE
      IF (CGIDPSX.EQ."PSX") THEN
         DGRTPSX(1:IGNFOLD-1) = DGRTSVX(1:IGNFOLD-1)
      ENDIF
      IF (CGIDPSZ.EQ."PSZ") THEN
         DGRTPSZ(1:IGNFOLD-1) = DGRTPPZ(1:IGNFOLD-1)
      ENDIF
   ENDIF
ENDIF
!
!*******************************************************************************************************************************************
! INVERSION WITH HV RATIO (BOTTOM) OR COMBINATION OF SEVERAL RATIOS
!*******************************************************************************************************************************************
!
IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
   DGRTHVB(1:IGNFOLD-1) = SQRT(2*DGINVVP(IGNLAYE)/DGINVVS(IGNLAYE))*DGTRSHB(1:IGNFOLD-1)/DGTRPZB(1:IGNFOLD-1)
ENDIF
!
!*******************************************************************************************************************************************
! INVERSION WITH HV RATIO (TOP) OR COMBINATION OF SEVERAL RATIOS
!*******************************************************************************************************************************************
!
IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
   DGRTHVT(1:IGNFOLD-1) = SQRT(2*DGINVVP(IGNLAYE)/DGINVVS(IGNLAYE))*DGTRSHT(1:IGNFOLD-1)/DGTRPZT(1:IGNFOLD-1)
ENDIF
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE TRFUNC
