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
SUBROUTINE READCF()
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGBNDPR                                                                                                          &
                         ,DGDELTT,DGDFREQ,DGSPBAN                                                                                          &
!
                         ,IGLAYPR                                                                                                          &
                         ,IGNBPRM,IGNFOLD,IGNT___                                                                                          &
!
                         ,CGPRTYP                                                                                                          &
                         ,CGIDPPX,CGIDPPZ,CGIDPSX,CGIDPSZ,CGIDSHY,CGIDSVX,CGIDSVZ,CGTYPIN
!
IMPLICIT NONE
!
INTEGER       :: I,IOS99,J,K
!
CHARACTER(92) :: FTMP
!
!*******************************************************************************************************************************************
! INITIALIZATION
!*******************************************************************************************************************************************
!
FTMP = "control_file_genetic_algo"
OPEN(UNIT=99,FILE=FTMP,STATUS="OLD",IOSTAT=IOS99)
CALL CHKOPN(FTMP,IOS99)
!
!*******************************************************************************************************************************************
! READ TYPE OF TRANSFER FUNCTIONS
!*******************************************************************************************************************************************
!
READ(UNIT=99,FMT='(10X,I10,2F10.0,17X,A3)') IGNT___,DGDELTT,DGSPBAN,CGTYPIN
!
IF (.NOT.((CGTYPIN.EQ."HVB").OR.(CGTYPIN.EQ."HVT").OR.(CGTYPIN.EQ."SPR").OR.(CGTYPIN.EQ."SUM"))) THEN
   WRITE(*,'(A)') "ERROR IN control_file_genetic_algo: NEITHER SPR OR HV INVERSION IS CHOSEN"
   STOP
ENDIF
!
IGNFOLD = IGNT___/2+1
DGDFREQ = 1.0D0/(DBLE(IGNT___)*DGDELTT)
READ(UNIT=99,FMT='(7X,A3)') CGIDSHY
READ(UNIT=99,FMT='(7X,A3)') CGIDSVX
READ(UNIT=99,FMT='(7X,A3)') CGIDSVZ
READ(UNIT=99,FMT='(7X,A3)') CGIDPPX
READ(UNIT=99,FMT='(7X,A3)') CGIDPPZ
READ(UNIT=99,FMT='(7X,A3)') CGIDPSX
READ(UNIT=99,FMT='(7X,A3)') CGIDPSZ
!
IF (.NOT.((CGIDPPX.EQ." PX").OR.(CGIDPPZ.EQ." PZ").OR.(CGIDPSX.EQ."PSX").OR.(CGIDPSZ.EQ."PSZ").OR.                                         &
          (CGIDSHY.EQ." SH").OR.(CGIDSVX.EQ."SVX").OR.(CGIDSVZ.EQ."SVZ"))) THEN
   WRITE(*,'(A)') "ERROR IN control_file_genetic_algo: NO TRANSFER FUNCTION IS OPTIMIZED"
   STOP
ENDIF
!
IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
   IF (.NOT.((CGIDSHY.EQ." SH").AND.(CGIDSVX.EQ."   ").AND.(CGIDSVZ.EQ."   ").AND.(CGIDPPX.EQ."   ")                                       &
        .AND.(CGIDPPZ.EQ." PZ").AND.(CGIDPSX.EQ."   ").AND.(CGIDPSZ.EQ."   "))) THEN
      WRITE(*,'(A)') "ERROR IN READ_CONTROL_FILE: ONLY SH AND PZ TRANSFER FUNCTION ARE COMPUTED"
      STOP
   ENDIF
ENDIF
!
WRITE(10,'(A,I10)') "NUMBER OF TIME STEPS: ",IGNT___
WRITE(10,'(A,I10)') "NUMBER OF FREQUENCY STEPS: ",IGNFOLD
WRITE(10,'(A,F10.4)') "TIME STEP: ",DGDELTT
WRITE(10,'(A,F10.4)') "FREQUENCY STEP: ",DGDFREQ
WRITE(10,'(A,F10.4)') "SMOOTHING WINDOW: ",DGSPBAN
WRITE(10,'(A)') ""
!
IF (CGTYPIN.EQ."HVB") WRITE(10,'(A)') "SPECTRAL HV RATIO ON BOTTOM"
IF (CGTYPIN.EQ."HVT") WRITE(10,'(A)') "SPECTRAL HV RATIO ON TOP"
IF (CGTYPIN.EQ."SPR") WRITE(10,'(A)') "SPECTRAL RATIOS BETWEEN TOP AND BOTTOM"
IF (CGTYPIN.EQ."SUM") WRITE(10,'(A)') "COMBINATION OF SPECTRAL RATIOS"
WRITE(10,'(A)') ""
!
IF (CGIDPPX.EQ." PX") WRITE(10,'(A)') " PX WAVE TRANSFER FUNCTION"
IF (CGIDPPZ.EQ." PZ") WRITE(10,'(A)') " PZ WAVE TRANSFER FUNCTION"
IF (CGIDPSX.EQ."PSX") WRITE(10,'(A)') "PSX WAVE TRANSFER FUNCTION"
IF (CGIDPSZ.EQ."PSZ") WRITE(10,'(A)') "PSZ WAVE TRANSFER FUNCTION"
IF (CGIDSHY.EQ." SH") WRITE(10,'(A)') " SH WAVE TRANSFER FUNCTION"
IF (CGIDSVX.EQ."SVX") WRITE(10,'(A)') "SVX WAVE TRANSFER FUNCTION"
IF (CGIDSVZ.EQ."SVZ") WRITE(10,'(A)') "SVZ WAVE TRANSFER FUNCTION"
WRITE(10,'(A)') ""
!
!*******************************************************************************************************************************************
! READ NUMBER OF PARAMETERS TO INVERT
!*******************************************************************************************************************************************
!
READ(UNIT=99,FMT='(I10)') IGNBPRM
WRITE(10,'(A,I10)') "NUMBER OF INVERTED PARAMETERS: ",IGNBPRM
!
!*******************************************************************************************************************************************
! READ PARAMETERS TO INVERT: TYPE AND INTERVAL OF VARIATIONS
!*******************************************************************************************************************************************
!
ALLOCATE(CGPRTYP(IGNBPRM),IGLAYPR(IGNBPRM,2),DGBNDPR(IGNBPRM,2))
DO I = 1,IGNBPRM
   READ(UNIT=99,FMT='(A10,2I10,2F10.0)') CGPRTYP(I),(IGLAYPR(I,J),J=1,2),(DGBNDPR(I,K),K=1,2)
ENDDO
!
CLOSE(99)
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE READCF
