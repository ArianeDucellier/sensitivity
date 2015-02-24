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
SUBROUTINE FCTOVL(ILINDIV,IND,ILNBPRM)
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGVECT1,DGVECT2                                                                                                  &
                         ,DGALPHA,DGBETA_,DGDEPTH,DGGAMMA,DGINANG                                                                          &
                         ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVPS,DGINVTH,DGINVVP,DGINVVS                                                  &
                         ,DGLAHEI,DGPVELO,DGSVELO,DGWEIPS                                                                                  &
!
                         ,IGLAYPR                                                                                                          &
                         ,IGLAOBJ,IGLAREF,IGNBPRM,IGNLAYE,IGSOBOL                                                                          &
!
                         ,CGPRTYP
IMPLICIT NONE
!
DOUBLE PRECISION    :: D
!
INTEGER, INTENT(IN) :: ILINDIV,ILNBPRM,IND
INTEGER             :: I,IP,L
!
!*******************************************************************************************************************************************
! FIRST CASE: WE USE DATA OF FIRST SOBOL MATRIX
!*******************************************************************************************************************************************
!
IF (IND.EQ.1) THEN
!
   !$OMP PARALLEL DO PRIVATE(L)
   DO IP = 1,IGNBPRM
      SELECTCASE (TRIM(ADJUSTL(CGPRTYP(IP))))
      CASE("VS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVVS(L) = DGVECT1(ILINDIV,IP) * DGSVELO(L)
         ENDDO
      CASE("VP")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVVP(L) = DGVECT1(ILINDIV,IP) * DGPVELO(L)
         ENDDO
      CASE("XA")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVAL(L) = DGVECT1(ILINDIV,IP) * DGALPHA(L)
         ENDDO
      CASE("XB")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVBE(L) = DGVECT1(ILINDIV,IP) * DGBETA_(L)
         ENDDO
      CASE("XG")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVGA(L) = DGVECT1(ILINDIV,IP) * DGGAMMA(L)
         ENDDO
      CASE("IN")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVIN(L) = DGVECT1(ILINDIV,IP) * DGINANG(L)
         ENDDO
      CASE("PS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVPS(L) = DGVECT1(ILINDIV,IP) * DGWEIPS(L)
         ENDDO
      CASE("TH")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVTH(L) = DGVECT1(ILINDIV,IP) * DGLAHEI(L)
         ENDDO
      ENDSELECT
   ENDDO
   !$OMP END PARALLEL DO
   IF ((IGLAOBJ.NE.0).AND.(IGLAREF.NE.0)) THEN
      D = 0.0D0
      DO I = IGLAOBJ,IGLAREF
         D = D + DGINVTH(I)
      ENDDO
      DGINVTH(IGLAREF+1) = DGDEPTH(IGLAREF+1) - D
      IF (DGINVTH(IGLAREF+1).LE.0.0D0) THEN
         WRITE (*,'(A)') "ERROR: LAST LAYER HAS NEGATIVE THICKNESS"
         STOP
      ENDIF
   ENDIF
!
!*******************************************************************************************************************************************
! SECOND CASE: WE USE DATA OF SECOND SOBOL MATRIX
!*******************************************************************************************************************************************
!
ELSEIF (IND.EQ.2) THEN
!
   !$OMP PARALLEL DO PRIVATE(L)
   DO IP = 1,IGNBPRM
      SELECTCASE (TRIM(ADJUSTL(CGPRTYP(IP))))
      CASE("VS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVVS(L) = DGVECT2(ILINDIV,IP) * DGSVELO(L)
         ENDDO
      CASE("VP")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVVP(L) = DGVECT2(ILINDIV,IP) * DGPVELO(L)
         ENDDO
      CASE("XA")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVAL(L) = DGVECT2(ILINDIV,IP) * DGALPHA(L)
         ENDDO
      CASE("XB")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVBE(L) = DGVECT2(ILINDIV,IP) * DGBETA_(L)
         ENDDO
      CASE("XG")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVGA(L) = DGVECT2(ILINDIV,IP) * DGGAMMA(L)
         ENDDO
      CASE("IN")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVIN(L) = DGVECT2(ILINDIV,IP) * DGINANG(L)
         ENDDO
      CASE("PS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVPS(L) = DGVECT2(ILINDIV,IP) * DGWEIPS(L)
         ENDDO
      CASE("TH")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            DGINVTH(L) = DGVECT2(ILINDIV,IP) * DGLAHEI(L)
         ENDDO
      ENDSELECT
   ENDDO
   !$OMP END PARALLEL DO
   IF ((IGLAOBJ.NE.0).AND.(IGLAREF.NE.0)) THEN
      D = 0.0D0
      DO I = IGLAOBJ,IGLAREF
         D = D + DGINVTH(I)
      ENDDO
      DGINVTH(IGLAREF+1) = DGDEPTH(IGLAREF+1) - D
      IF (DGINVTH(IGLAREF+1).LE.0.0D0) THEN
         WRITE (*,'(A)') "ERROR: LAST LAYER HAS NEGATIVE THICKNESS"
         STOP
      ENDIF
   ENDIF
!
!*******************************************************************************************************************************************
! THIRD CASE: WE USE DATA OF FIRST AND SECOND SOBOL MATRICES
!*******************************************************************************************************************************************
!
ELSEIF (IND.EQ.3) THEN
!
   !$OMP PARALLEL DO PRIVATE(L)
   DO IP = 1,IGNBPRM
      SELECTCASE (TRIM(ADJUSTL(CGPRTYP(IP))))
      CASE("VS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVVS(L) = DGVECT2(ILINDIV,IP) * DGSVELO(L)
            ELSE
               DGINVVS(L) = DGVECT1(ILINDIV,IP) * DGSVELO(L)
            ENDIF
         ENDDO
      CASE("VP")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVVP(L) = DGVECT2(ILINDIV,IP) * DGPVELO(L)
            ELSE
               DGINVVP(L) = DGVECT1(ILINDIV,IP) * DGPVELO(L)
            ENDIF
         ENDDO
      CASE("XA")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVAL(L) = DGVECT2(ILINDIV,IP) * DGALPHA(L)
            ELSE
               DGINVAL(L) = DGVECT1(ILINDIV,IP) * DGALPHA(L)
            ENDIF
         ENDDO
      CASE("XB")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVBE(L) = DGVECT2(ILINDIV,IP) * DGBETA_(L)
            ELSE
               DGINVBE(L) = DGVECT1(ILINDIV,IP) * DGBETA_(L)
            ENDIF
         ENDDO
      CASE("XG")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVGA(L) = DGVECT2(ILINDIV,IP) * DGGAMMA(L)
            ELSE
               DGINVGA(L) = DGVECT1(ILINDIV,IP) * DGGAMMA(L)
            ENDIF
         ENDDO
      CASE("IN")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVIN(L) = DGVECT2(ILINDIV,IP) * DGINANG(L)
            ELSE
               DGINVIN(L) = DGVECT1(ILINDIV,IP) * DGINANG(L)
            ENDIF
         ENDDO
      CASE("PS")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVPS(L) = DGVECT2(ILINDIV,IP) * DGWEIPS(L)
            ELSE
               DGINVPS(L) = DGVECT1(ILINDIV,IP) * DGWEIPS(L)
            ENDIF
         ENDDO
      CASE("TH")
         DO L = IGLAYPR(IP,1),IGLAYPR(IP,2)
            IF (IP.EQ.ILNBPRM) THEN
               DGINVTH(L) = DGVECT2(ILINDIV,IP) * DGLAHEI(L)
            ELSE
               DGINVTH(L) = DGVECT1(ILINDIV,IP) * DGLAHEI(L)
            ENDIF
         ENDDO
      ENDSELECT
   ENDDO
   !$OMP END PARALLEL DO
   IF ((IGLAOBJ.NE.0).AND.(IGLAREF.NE.0)) THEN
      D = 0.0D0
      DO I = IGLAOBJ,IGLAREF
         D = D + DGINVTH(I)
      ENDDO
      DGINVTH(IGLAREF+1) = DGDEPTH(IGLAREF+1) - D
      IF (DGINVTH(IGLAREF+1).LE.0.0D0) THEN
         WRITE (*,'(A)') "ERROR: LAST LAYER HAS NEGATIVE THICKNESS"
         STOP
      ENDIF
   ENDIF
!
ENDIF
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE FCTOVL
