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
SUBROUTINE SHWAVE()
!
!*******************************************************************************************************************************************
!
USE AAMODU_SEISMO_TOOLS, ONLY : SPEWIN
!
USE AAMODU_GLOBVA, ONLY       : DGDENSI,DGDEPTH                                                                                            &
                               ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVTH,DGINVVS                                                            &
                               ,DGRTSHY                                                                                                    &
                               ,DGTRSHB,DGTRSHT                                                                                            &
                               ,DGDFREQ,DGSPBAN                                                                                            &
!
                               ,IGIDDMP,IGNFOLD,IGNLAYE                                                                                    &
!
                               ,CGTYPIN
!
IMPLICIT NONE
!
DOUBLE COMPLEX  , DIMENSION(:), ALLOCATABLE :: ZSHSPR,ZV_HS,ZV_OBJ,ZW_HS_D,ZW_HS_U,ZW_OBJ_U
DOUBLE COMPLEX                              :: Z_INI(2,2),ZF_HS(2,2),ZF_OBJ(2,2),ZPROPAG(2,2)
DOUBLE COMPLEX                              :: Z_ETA,ZE,ZG,ZG_HS,ZI,ZP,ZVS,ZVS_HS
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AMZV_HS,AMZV_OBJ
DOUBLE PRECISION                            :: D,DAMP,DW,OPT_INC,PI,TI0,W
!
INTEGER                                     :: I,IFR,IL
!
!*******************************************************************************************************************************************
! INITIALIZATION
!*******************************************************************************************************************************************
!
ZI = CMPLX(0.0D0,1.0D0)
PI = ACOS(-1.0D0)
!
TI0     = 0.0D0
DW      = 2.0D0*PI*DGDFREQ
OPT_INC = DGINVIN(IGNLAYE)*PI/180.0D0
!
ALLOCATE(ZSHSPR(IGNFOLD),ZW_HS_U(IGNFOLD),ZW_HS_D(IGNFOLD),ZW_OBJ_U(IGNFOLD),ZV_HS(IGNFOLD),ZV_OBJ(IGNFOLD),AMZV_OBJ(IGNFOLD),AMZV_HS(IGNFOLD))
!
ZV_OBJ(:)= CMPLX(0.0D0,0.0D0)
ZV_HS(:) = CMPLX(0.0D0,0.0D0)
ZSHSPR(:)= CMPLX(1.0D0,0.0D0)
!
!*******************************************************************************************************************************************
! START LOOP ON PULSATION (FREQUENCY = 0.0 GIVES NaN)
!*******************************************************************************************************************************************
!
!$OMP PARALLEL DO PRIVATE(W,DAMP,ZG_HS,ZVS_HS,ZP,ZG,ZVS,Z_ETA,D,ZF_HS,Z_INI,IL,ZPROPAG,ZF_OBJ,ZE)
DO IFR = 1,(IGNFOLD-1)
!
   W = IFR*DW
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! INITIALIZATION FOR HALF-SPACE (REFERENCE)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   DAMP     = 0.0D0
   CALL       CMPDMP(W,DGINVAL(IGNLAYE),DGINVBE(IGNLAYE),DGINVGA(IGNLAYE),DGINVVS(IGNLAYE),DAMP,IGIDDMP)
   ZG_HS    = DGDENSI(IGNLAYE)*DGINVVS(IGNLAYE)**2*CMPLX(1.0D0,2.0D0*DAMP)
   ZVS_HS   = SQRT(ZG_HS/CMPLX(DGDENSI(IGNLAYE),0.0D0))
   ZP       = SIN(CMPLX(OPT_INC,0.0D0))/ZVS_HS
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF MATRIX Fn(Zn)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   DAMP     = 0.0D0
   CALL       CMPDMP(W,DGINVAL(IGNLAYE-1),DGINVBE(IGNLAYE-1),0.0D0,DGINVVS(IGNLAYE-1),DAMP,IGIDDMP)
   ZG       = DGDENSI(IGNLAYE-1)*DGINVVS(IGNLAYE-1)**2*CMPLX(1.0D0,2.0D0*DAMP) 
   ZVS      = SQRT(ZG/CMPLX(DGDENSI(IGNLAYE-1),0.0D0))
!
   Z_ETA    = SQRT(CMPLX(1.0D0,0.0D0)/ZVS**2-ZP**2)
   IF (AIMAG(Z_ETA).LT.0.0D0) Z_ETA = -Z_ETA
!
   D = 0.0D0
   DO I = 1,(IGNLAYE-1)
      D = D + DGINVTH(I)
   ENDDO
!
   ZF_HS(1,1) =                EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(2,1) = +ZI*W*ZG*Z_ETA*EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(1,2) =                EXP(-ZI*W*Z_ETA*(D-TI0))
   ZF_HS(2,2) = -ZI*W*ZG*Z_ETA*EXP(-ZI*W*Z_ETA*(D-TI0))
!
   Z_INI(1,1)    =            1.0D0/2.0D0*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(2,1)    =            1.0D0/2.0D0*EXP(+ZI*W*Z_ETA*(D-TI0))
   Z_INI(1,2)    = -ZI/(2.0D0*W*ZG*Z_ETA)*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(2,2)    = +ZI/(2.0D0*W*ZG*Z_ETA)*EXP(+ZI*W*Z_ETA*(D-TI0))
!
   DO IL = IGNLAYE-1,1,-1
!
!-------------------------------------------------------------------------------------------------------------------------------------------
!                                       -1
! CALCULUS OF THE PROPAGATOR MATRIX : Fn(Zn) x P(Zn,Z0)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
      DAMP     = 0.0D0
      CALL       CMPDMP(W,DGINVAL(IL),DGINVBE(IL),0.0D0,DGINVVS(IL),DAMP,IGIDDMP)
      ZG       = DGDENSI(IL)*DGINVVS(IL)**2*CMPLX(1.0D0,2.0D0*DAMP) 
      ZVS      = SQRT(ZG/CMPLX(DGDENSI(IL),0.0D0))
!
      Z_ETA    = SQRT(CMPLX(1.0D0,0.0D0)/ZVS**2-ZP**2)
      IF (AIMAG(Z_ETA).LT.0.0D0) Z_ETA = -Z_ETA
!
      ZPROPAG(1,1) =                      COS(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(2,1) =         (W*ZG*Z_ETA)*SIN(W*Z_ETA*DGINVTH(IL)) 
      ZPROPAG(1,2) =  -1.0D0/(W*ZG*Z_ETA)*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(2,2) =                      COS(W*Z_ETA*DGINVTH(IL))
!
      CALL MULTI_MAT_MXM(Z_INI,ZPROPAG,2) 
!
   ENDDO
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF P(Zn,Z0) x F1(Z0)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZF_OBJ(1,1) =                EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(1,2) =                EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(2,1) = +ZI*W*ZG*Z_ETA*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(2,2) = -ZI*W*ZG*Z_ETA*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
!
   CALL MULTI_MAT_MXM(Z_INI,ZF_OBJ,2)
!
   ZE = EXP(-2.0D0*ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! INPUT MOTION
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZW_HS_U(IFR) = DCMPLX(1.0D0,0.0D0)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! SOLUTION OF THE SYSTEM 2 UNKNOWNS - 2 EQUATIONS
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZW_OBJ_U(IFR) = ZW_HS_U(IFR)/(Z_INI(2,1)*ZE+Z_INI(2,2))
   ZW_HS_D(IFR)  = ZW_OBJ_U(IFR)*(Z_INI(1,1)*ZE+Z_INI(1,2))
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF DISPLACEMENT
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZV_OBJ(IFR)  = (ZW_OBJ_U(IFR)*ZE*ZF_OBJ(1,1) + ZW_OBJ_U(IFR)*ZF_OBJ(1,2))
   ZV_HS(IFR)   = (ZW_HS_D(IFR)*ZF_HS(1,1)      + ZW_HS_U(IFR)*ZF_HS(1,2))   
   ZSHSPR(IFR)  = ZV_OBJ(IFR)/ZV_HS(IFR)
ENDDO
!$OMP END PARALLEL DO
!
!*******************************************************************************************************************************************
! PARZEN'S SPECTRAL WINDOW
!*******************************************************************************************************************************************
!
AMZV_OBJ(:) = ABS(ZV_OBJ(:))
AMZV_HS(:)  = ABS(ZV_HS(:))
IF (DGSPBAN.GT.0.0D0) THEN
   CALL SPEWIN(AMZV_OBJ,DGDFREQ,DGSPBAN)
   CALL SPEWIN(AMZV_HS ,DGDFREQ,DGSPBAN)
ENDIF
!
! INVERSION WITH SPECTRAL RATIO
! SH TRANSFER FUNCTION
!
IF ( (CGTYPIN.EQ."SPR") .OR. (CGTYPIN.EQ."SUM") ) THEN
   !$OMP PARALLEL DO
   DO I = 1,(IGNFOLD-1)
      DGRTSHY(I) = AMZV_OBJ(I)/AMZV_HS(I)
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
!
! INVERSION WITH HV RATIO (BOTTOM) OR COMBINATION OF SEVERAL RATIOS
!
IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
   !$OMP PARALLEL DO
   DO I = 1,(IGNFOLD-1)
      DGTRSHB(I) = AMZV_HS(I)
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
!
! INVERSION WITH HV RATIO (TOP) OR COMBINATION OF SEVERAL RATIOS
!
IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
   !$OMP PARALLEL DO
   DO I = 1,(IGNFOLD-1)
      DGTRSHT(I) = AMZV_OBJ(I)
   ENDDO
   !$OMP END PARALLEL DO
ENDIF
!
DEALLOCATE(ZW_HS_U,ZW_HS_D,ZW_OBJ_U,ZV_HS,ZV_OBJ,ZSHSPR,AMZV_OBJ,AMZV_HS)
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE SHWAVE
