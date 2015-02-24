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
SUBROUTINE READSC()
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGBNDPR,DGPRINI                                                                                                 &
                         ,DGALPHA,DGBETA_,DGDENSI,DGDEPTH,DGGAMMA,DGINANG                                                                 &
                         ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVNU,DGINVPS,DGINVTH,DGINVVP,DGINVVS                                         &
                         ,DGLAHEI,DGPOISR,DGPVELO                                                                                         &
                         ,DGRTHVB,DGRTHVT,DGRTPPX,DGRTPPZ,DGRTPSX,DGRTPSZ,DGRTSHY,DGRTSVX,DGRTSVZ                                         &
                         ,DGSHEMO,DGSVELO                                                                                                 &
                         ,DGTRPZB,DGTRPZT,DGTRSHB,DGTRSHT                                                                                 &
                         ,DGUWEIG,DGWEIPS                                                                                                 &
!
                         ,IGLAYPR                                                                                                         &
                         ,IGIDDMP,IGLAOBJ,IGLAREF,IGNBPRM,IGNFOLD,IGNLAYE                                                                 &
!
                         ,CGPRTYP                                                                                                         &
                         ,CGFINPU
!
IMPLICIT NONE
!
DOUBLE PRECISION :: D
!
INTEGER          :: I,IOS99,J
!
!*******************************************************************************************************************************************
! READ INITIAL PARAMETERS OF THE SOIL COLUMN
!*******************************************************************************************************************************************
!
OPEN(UNIT=99,FILE=CGFINPU,STATUS="OLD",IOSTAT=IOS99)
CALL CHKOPN(CGFINPU,IOS99)
READ(UNIT=99,FMT='(3I10,10X,I10)') IGNLAYE,IGLAOBJ,IGLAREF,IGIDDMP
IF ((IGLAOBJ.NE.0).AND.(IGLAREF.NE.0).AND.(IGLAREF.GT.(IGNLAYE-2))) THEN
   WRITE(*,'(A)') "ERROR IN *.in: YOU CANNOT INVERT THE THICKNESS OF THE LAST LAYER"
   STOP
ENDIF
ALLOCATE(DGALPHA(IGNLAYE),DGBETA_(IGNLAYE),DGDENSI(IGNLAYE),DGGAMMA(IGNLAYE),DGINANG(IGNLAYE),DGLAHEI(IGNLAYE),DGPOISR(IGNLAYE)            &
        ,DGPVELO(IGNLAYE),DGSHEMO(IGNLAYE),DGSVELO(IGNLAYE),DGUWEIG(IGNLAYE),DGWEIPS(IGNLAYE)                                              &
        ,DGDEPTH(0:IGNLAYE),DGPRINI(IGNLAYE,IGNBPRM))
!
DGDEPTH(0) = 0.0D0
D          = 0.0D0
DO I = 1,IGNLAYE
   READ(UNIT=99,FMT='(9F10.0)') DGPOISR(I),DGUWEIG(I),DGALPHA(I),DGBETA_(I),DGLAHEI(I),DGSVELO(I),DGINANG(I),DGWEIPS(I),DGGAMMA(I)
   DGDENSI(I) = DGUWEIG(I)/9.80665D0*10**4
   DGSHEMO(I) = DGDENSI(I)*DGSVELO(I)**2
   DGPVELO(I) = DGSVELO(I)*SQRT((2.0D0-2.0D0*DGPOISR(I))/(1.0D0-2.0D0*DGPOISR(I)))
   D          = D + DGLAHEI(I)
   DGDEPTH(I) = D
ENDDO
CLOSE(99)
!
!*******************************************************************************************************************************************
! FILL MATRIX OF INITIAL SOIL COLUMN
!*******************************************************************************************************************************************
!
DO J = 1,IGNBPRM
   DO I = 1,IGNLAYE
      DGPRINI(I,J) = 0.0D0
   ENDDO
ENDDO
!
DO I = 1,IGNBPRM
   IF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."VS") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGSVELO(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."VP") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGPVELO(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."XA") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGALPHA(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."XB") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGBETA_(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."XG") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGGAMMA(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."IN") THEN
      IF (IGLAYPR(I,1).NE.IGNLAYE) IGLAYPR(I,1) = IGNLAYE
      IF (IGLAYPR(I,2).NE.IGNLAYE) IGLAYPR(I,2) = IGNLAYE
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGINANG(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."PS") THEN
      IF (IGLAYPR(I,1).NE.IGNLAYE) IGLAYPR(I,1) = IGNLAYE
      IF (IGLAYPR(I,2).NE.IGNLAYE) IGLAYPR(I,2) = IGNLAYE
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGWEIPS(J)
      ENDDO
   ELSEIF (TRIM(ADJUSTL(CGPRTYP(I))).EQ."TH") THEN
      DO J = IGLAYPR(I,1),IGLAYPR(I,2)
         DGPRINI(J,I)=DGLAHEI(J)
      ENDDO
   ENDIF
ENDDO
!
DO I = 1,IGNBPRM
   WRITE(10,'(2A)') "TYPE OF PARAMETER: ",TRIM(ADJUSTL(CGPRTYP(I)))
   DO J = IGLAYPR(I,1),IGLAYPR(I,2)
      WRITE(10,'(A,I10,A,F10.4,A,F10.4)') "   LAYER ",J," - VARIATIONS BETWEEN "                                                           &
      ,DGPRINI(J,I)*DGBNDPR(I,1)," AND ",DGPRINI(J,I)*DGBNDPR(I,2)
   ENDDO
ENDDO
!
!*******************************************************************************************************************************************
! MEMORY ALLOCATION FOR INVERSION PARAMETERS AND THEORETICAL TRANSFER FUNCTIONS
!*******************************************************************************************************************************************
!
ALLOCATE(DGRTHVB(IGNFOLD),DGRTHVT(IGNFOLD),DGRTPPX(IGNFOLD),DGRTPPZ(IGNFOLD),DGRTPSX(IGNFOLD),DGRTPSZ(IGNFOLD)                             &
        ,DGRTSHY(IGNFOLD),DGRTSVX(IGNFOLD),DGRTSVZ(IGNFOLD))
!
DGRTHVB(:) = 0.0D0
DGRTHVT(:) = 0.0D0
DGRTPPX(:) = 0.0D0
DGRTPPZ(:) = 0.0D0
DGRTPSX(:) = 0.0D0
DGRTPSZ(:) = 0.0D0
DGRTSHY(:) = 0.0D0
DGRTSVX(:) = 0.0D0
DGRTSVZ(:) = 0.0D0
!
ALLOCATE(DGTRPZB(IGNFOLD),DGTRPZT(IGNFOLD),DGTRSHB(IGNFOLD),DGTRSHT(IGNFOLD))
!
DGTRPZB(:) = 0.0D0
DGTRPZT(:) = 0.0D0
DGTRSHB(:) = 0.0D0
DGTRSHT(:) = 0.0D0
!
ALLOCATE(DGINVAL(IGNLAYE),DGINVBE(IGNLAYE),DGINVGA(IGNLAYE),DGINVIN(IGNLAYE),DGINVNU(IGNLAYE),DGINVPS(IGNLAYE)                             &
        ,DGINVVP(IGNLAYE),DGINVVS(IGNLAYE),DGINVTH(IGNLAYE))
!
CLOSE(99)
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE READSC
