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
SUBROUTINE PPWAVE()
!
!*******************************************************************************************************************************************
!
USE AAMODU_SEISMO_TOOLS, ONLY : SPEWIN
!
USE AAMODU_GLOBVA, ONLY       : DGDENSI,DGDEPTH                                                                                            &
                               ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVNU,DGINVTH,DGINVVS                                                    &
                               ,DGRTPPX,DGRTPPZ                                                                                            &
                               ,DGTRPZB,DGTRPZT                                                                                            &
                               ,DGDFREQ,DGSPBAN                                                                                            &
!
                               ,IGIDDMP,IGNFOLD,IGNLAYE                                                                                    &
!
                               ,CGIDPPX,CGIDPPZ,CGIDPSX,CGIDPSZ,CGTYPIN
!
IMPLICIT NONE
!
DOUBLE COMPLEX  , DIMENSION(:), ALLOCATABLE :: ZPDOWN_HS,ZPDOWN_OBJ,ZPUP_HS,ZPUP_OBJ,ZSDOWN_HS,ZSDOWN_OBJ,ZSUP_HS,ZSUP_OBJ                 &
                                              ,ZGA,ZGB,ZGC,ZGD,ZGE,ZGF,ZGG,ZGH,ZGI,ZSVSPRX,ZSVSPRZ,ZUX_HS,ZUX_OBJ,ZUZ_HS,ZUZ_OBJ
DOUBLE COMPLEX                              :: Z_INI(4,4),ZF_HS(4,4),ZF_OBJ(4,4),ZPROPAG(4,4)
DOUBLE COMPLEX                              :: Z_ETA,Z_XI,ZG,ZG_HS,ZI,ZP,ZVP,ZVP_HS,ZVS,ZVS_HS
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AMZUX_HS,AMZUX_OBJ,AMZUZ_HS,AMZUZ_OBJ 
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
ALLOCATE(ZSVSPRX(IGNFOLD),ZSVSPRZ(IGNFOLD),ZSUP_HS(IGNFOLD),ZPUP_HS(IGNFOLD),ZSDOWN_HS(IGNFOLD),ZPDOWN_HS(IGNFOLD)                          &
        ,ZSUP_OBJ(IGNFOLD),ZPUP_OBJ(IGNFOLD),ZSDOWN_OBJ(IGNFOLD),ZPDOWN_OBJ(IGNFOLD),ZUX_OBJ(IGNFOLD)                                       &
        ,ZUZ_OBJ(IGNFOLD),ZUX_HS(IGNFOLD),ZUZ_HS(IGNFOLD),AMZUX_HS(IGNFOLD),AMZUX_OBJ(IGNFOLD),AMZUZ_HS(IGNFOLD),AMZUZ_OBJ(IGNFOLD))
ALLOCATE(ZGA(IGNFOLD),ZGB(IGNFOLD),ZGC(IGNFOLD),ZGD(IGNFOLD),ZGE(IGNFOLD),ZGF(IGNFOLD),ZGG(IGNFOLD),ZGH(IGNFOLD),ZGI(IGNFOLD))
!
ZUX_OBJ(:) = CMPLX(0.0D0,0.0D0)
ZUX_HS(:)  = CMPLX(0.0D0,0.0D0)
!
!*******************************************************************************************************************************************
! START LOOP ON PULSATION (FREQUENCY = 0.0 GIVES NaN)
!*******************************************************************************************************************************************
!
!$OMP PARALLEL DO PRIVATE(W,DAMP,ZG_HS,ZVS_HS,ZVP_HS,ZP,ZG,ZVS,ZVP,Z_ETA,Z_XI,D,ZF_HS,Z_INI,IL,ZPROPAG,ZF_OBJ)
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
   ZVP_HS   = ZVS_HS*SQRT((2.0D0-2.0D0*DGINVNU(IGNLAYE))/(1.0D0-2.0D0*DGINVNU(IGNLAYE)))
   ZP       = SIN(CMPLX(OPT_INC,0.0D0))/ZVP_HS
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF MATRIX Fn(Zn)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   DAMP     = 0.0D0
   CALL       CMPDMP(W,DGINVAL(IGNLAYE-1),DGINVBE(IGNLAYE-1),0.0D0,DGINVVS(IGNLAYE-1),DAMP,IGIDDMP)
   ZG       = DGDENSI(IGNLAYE-1)*DGINVVS(IGNLAYE-1)**2*CMPLX(1.0D0,2.0D0*DAMP)
   ZVS      = SQRT(ZG/CMPLX(DGDENSI(IGNLAYE-1),0.0D0))
   ZVP      = ZVS*SQRT((2.0D0-2.0D0*DGINVNU(IGNLAYE-1))/(1.0D0-2.0D0*DGINVNU(IGNLAYE-1)))
!
   Z_ETA    = SQRT(CMPLX(1.0D0,0.0D0)/ZVS**2-ZP**2)
   IF (AIMAG(Z_ETA).LT.0.0D0) Z_ETA = -Z_ETA
!
   Z_XI     = SQRT(CMPLX(1.0D0,0.0D0)/ZVP**2-ZP**2)
   IF (AIMAG(Z_XI).LT.0.0D0) Z_XI = -Z_XI
!
   D = 0.0D0
   DO I = 1,(IGNLAYE-1)
      D = D + DGINVTH(I)
   ENDDO
!
   ZF_HS(1,1) =                                                 +ZVP*ZP*EXP(+ZI*W*Z_XI *(D-TI0))
   ZF_HS(2,1) =                                               +ZVP*Z_XI*EXP(+ZI*W*Z_XI *(D-TI0))
   ZF_HS(3,1) =       +2.0D0*ZI*W*DGDENSI(IGNLAYE-1)*ZVP*ZVS**2*ZP*Z_XI*EXP(+ZI*W*Z_XI *(D-TI0))
   ZF_HS(4,1) = +ZI*W*DGDENSI(IGNLAYE-1)*ZVP*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(+ZI*W*Z_XI *(D-TI0))
   ZF_HS(1,2) =                                              +ZVS*Z_ETA*EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(2,2) =                                                 -ZVS*ZP*EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(3,2) = +ZI*W*DGDENSI(IGNLAYE-1)*ZVS*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(4,2) =          -2.0D0*ZI*W*DGDENSI(IGNLAYE-1)*ZVS**3*ZP*Z_ETA*EXP(+ZI*W*Z_ETA*(D-TI0))
   ZF_HS(1,3) =                                                 +ZVP*ZP*EXP(-ZI*W*Z_XI *(D-TI0))
   ZF_HS(2,3) =                                               -ZVP*Z_XI*EXP(-ZI*W*Z_XI *(D-TI0))
   ZF_HS(3,3) =       -2.0D0*ZI*W*DGDENSI(IGNLAYE-1)*ZVP*ZVS**2*ZP*Z_XI*EXP(-ZI*W*Z_XI *(D-TI0))
   ZF_HS(4,3) = +ZI*W*DGDENSI(IGNLAYE-1)*ZVP*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_XI *(D-TI0))
   ZF_HS(1,4) =                                              +ZVS*Z_ETA*EXP(-ZI*W*Z_ETA*(D-TI0))
   ZF_HS(2,4) =                                                 +ZVS*ZP*EXP(-ZI*W*Z_ETA*(D-TI0))
   ZF_HS(3,4) = -ZI*W*DGDENSI(IGNLAYE-1)*ZVS*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_ETA*(D-TI0))
   ZF_HS(4,4) =          -2.0D0*ZI*W*DGDENSI(IGNLAYE-1)*ZVS**3*ZP*Z_ETA*EXP(-ZI*W*Z_ETA*(D-TI0)) !VERIF1 31/03/08 OK
!
   Z_INI(1,1) =                                    (ZVS**2*ZP)/(ZVP)*EXP(-ZI*W*Z_XI *(D-TI0))
   Z_INI(2,1) =         (1.0D0-2.0D0*ZVS**2*ZP**2)/(2.0D0*ZVS*Z_ETA)*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(3,1) =                                    (ZVS**2*ZP)/(ZVP)*EXP(+ZI*W*Z_XI *(D-TI0))
   Z_INI(4,1) =         (1.0D0-2.0D0*ZVS**2*ZP**2)/(2.0D0*ZVS*Z_ETA)*EXP(+ZI*W*Z_ETA*(D-TI0))
   Z_INI(1,2) =          (1.0D0-2.0D0*ZVS**2*ZP**2)/(2.0D0*ZVP*Z_XI)*EXP(-ZI*W*Z_XI *(D-TI0))
   Z_INI(2,2) =                                              -ZVS*ZP*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(3,2) = (-1.0D0*(1.0D0-2.0D0*ZVS**2*ZP**2))/(2.0D0*ZVP*Z_XI)*EXP(+ZI*W*Z_XI *(D-TI0))
   Z_INI(4,2) =                                               ZVS*ZP*EXP(+ZI*W*Z_ETA*(D-TI0))
   Z_INI(1,3) =       (-ZI*ZP)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVP*Z_XI)*EXP(-ZI*W*Z_XI *(D-TI0))
   Z_INI(2,3) =               (-ZI)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVS)*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(3,3) =        (ZI*ZP)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVP*Z_XI)*EXP(+ZI*W*Z_XI *(D-TI0))
   Z_INI(4,3) =                (ZI)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVS)*EXP(+ZI*W*Z_ETA*(D-TI0))
   Z_INI(1,4) =               (-ZI)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVP)*EXP(-ZI*W*Z_XI *(D-TI0))
   Z_INI(2,4) =       (ZI*ZP)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVS*Z_ETA)*EXP(-ZI*W*Z_ETA*(D-TI0))
   Z_INI(3,4) =               (-ZI)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVP)*EXP(+ZI*W*Z_XI *(D-TI0))
   Z_INI(4,4) =       (ZI*ZP)/(2.0D0*W*DGDENSI(IGNLAYE-1)*ZVS*Z_ETA)*EXP(+ZI*W*Z_ETA*(D-TI0)) !VERIF1 31/03/08 OK
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
      ZVP      = ZVS*SQRT((2.0D0-2.0D0*DGINVNU(IL))/(1.0D0-2.0D0*DGINVNU(IL)))
!
      Z_ETA    = SQRT(CMPLX(1.0D0,0.0D0)/ZVS**2-ZP**2)
      IF (AIMAG(Z_ETA).LT.0.0D0) Z_ETA = -Z_ETA
!
      Z_XI     = SQRT(CMPLX(1.0D0,0.0D0)/ZVP**2-ZP**2)
      IF (AIMAG(Z_XI).LT.0.0D0) Z_XI = -Z_XI
!
      ZPROPAG(1,1) =                                           2.0D0*ZVS**2*ZP**2*COS(W*Z_XI*DGINVTH(IL)) +                            (1.0D0-2.0D0*ZVS**2*ZP**2)*COS(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(2,1) =                                      2.0D0*ZI*ZVS**2*ZP*Z_XI*SIN(W*Z_XI*DGINVTH(IL)) -            (ZI*ZP)/(Z_ETA)*(1.0D0-2.0D0*ZVS**2*ZP**2)*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(3,1) =                       -4.0D0*W*DGDENSI(IL)*ZVS**4*ZP**2*Z_XI*SIN(W*Z_XI*DGINVTH(IL)) - (W*DGDENSI(IL))/(Z_ETA)*(1.0D0-2.0D0*ZVS**2*ZP**2)**2*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(4,1) = 2.0D0*ZI*W*DGDENSI(IL)*ZVS**2*ZP*(1.0D0-2.0D0*ZVS**2*ZP**2)*(COS(W*Z_XI*DGINVTH(IL)) -                                                       COS(W*Z_ETA*DGINVTH(IL)))
      ZPROPAG(1,2) =                    (ZI*ZP)/(Z_XI)*(1.0D0-2.0D0*ZVS**2*ZP**2)*SIN(W*Z_XI*DGINVTH(IL)) -                              2.0D0*ZI*ZVS**2*ZP*Z_ETA*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(2,2) =                                   (1.0D0-2.0D0*ZVS**2*ZP**2)*COS(W*Z_XI*DGINVTH(IL)) +                                    2.0D0*ZVS**2*ZP**2*COS(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(3,2) = 2.0D0*ZI*W*DGDENSI(IL)*ZVS**2*ZP*(1.0D0-2.0D0*ZVS**2*ZP**2)*(COS(W*Z_XI*DGINVTH(IL)) -                                                       COS(W*Z_ETA*DGINVTH(IL)))
      ZPROPAG(4,2) =        -(W*DGDENSI(IL))/(Z_XI)*(1.0D0-2.0D0*ZVS**2*ZP**2)**2*SIN(W*Z_XI*DGINVTH(IL)) -                4.0D0*W*DGDENSI(IL)*ZVS**4*ZP**2*Z_ETA*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(1,3) =                                 (ZP**2)/(W*DGDENSI(IL)*Z_XI)*SIN(W*Z_XI*DGINVTH(IL)) +                               (Z_ETA)/(W*DGDENSI(IL))*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(2,3) =                                    -(ZI*ZP)/(W*DGDENSI(IL))*(COS(W*Z_XI*DGINVTH(IL)) -                                                       COS(W*Z_ETA*DGINVTH(IL)))
      ZPROPAG(3,3) =                                           2.0D0*ZVS**2*ZP**2*COS(W*Z_XI*DGINVTH(IL)) +                            (1.0D0-2.0D0*ZVS**2*ZP**2)*COS(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(4,3) =                    (ZI*ZP)/(Z_XI)*(1.0D0-2.0D0*ZVS**2*ZP**2)*SIN(W*Z_XI*DGINVTH(IL)) -                              2.0D0*ZI*ZVS**2*ZP*Z_ETA*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(1,4) =                                    -(ZI*ZP)/(W*DGDENSI(IL))*(COS(W*Z_XI*DGINVTH(IL)) -                                                       COS(W*Z_ETA*DGINVTH(IL)))
      ZPROPAG(2,4) =                                       (Z_XI)/(W*DGDENSI(IL))*SIN(W*Z_XI*DGINVTH(IL)) +                         (ZP**2)/(W*DGDENSI(IL)*Z_ETA)*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(3,4) =                                      2.0D0*ZI*ZVS**2*ZP*Z_XI*SIN(W*Z_XI*DGINVTH(IL)) -            (ZI*ZP)/(Z_ETA)*(1.0D0-2.0D0*ZVS**2*ZP**2)*SIN(W*Z_ETA*DGINVTH(IL))
      ZPROPAG(4,4) =                                   (1.0D0-2.0D0*ZVS**2*ZP**2)*COS(W*Z_XI*DGINVTH(IL)) +                                    2.0D0*ZVS**2*ZP**2*COS(W*Z_ETA*DGINVTH(IL)) !VERIF1 OK 31/03/08
!
      CALL MULTI_MAT_MXM(Z_INI,ZPROPAG,4)
!
   ENDDO
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF P(Zn,Z0) x F1(Z0)
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZF_OBJ(1,1) =                                         +ZVP*ZP*EXP(+ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(2,1) =                                       +ZVP*Z_XI*EXP(+ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(3,1) =       +2.0D0*ZI*W*DGDENSI(1)*ZVP*ZVS**2*ZP*Z_XI*EXP(+ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(4,1) = +ZI*W*DGDENSI(1)*ZVP*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(+ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(1,2) =                                      +ZVS*Z_ETA*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(2,2) =                                         -ZVS*ZP*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(3,2) = +ZI*W*DGDENSI(1)*ZVS*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(4,2) =          -2.0D0*ZI*W*DGDENSI(1)*ZVS**3*ZP*Z_ETA*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(1,3) =                                         +ZVP*ZP*EXP(-ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(2,3) =                                       -ZVP*Z_XI*EXP(-ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(3,3) =       -2.0D0*ZI*W*DGDENSI(1)*ZVP*ZVS**2*ZP*Z_XI*EXP(-ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(4,3) = +ZI*W*DGDENSI(1)*ZVP*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_XI *(DGDEPTH(0)-TI0))
   ZF_OBJ(1,4) =                                      +ZVS*Z_ETA*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(2,4) =                                         +ZVS*ZP*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(3,4) = -ZI*W*DGDENSI(1)*ZVS*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
   ZF_OBJ(4,4) =          -2.0D0*ZI*W*DGDENSI(1)*ZVS**3*ZP*Z_ETA*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))
!
   CALL MULTI_MAT_MXM(Z_INI,ZF_OBJ,4)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! INPUT MOTION
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZSUP_HS(IFR) = CMPLX(0.0D0,0.0D0)
   ZPUP_HS(IFR) = CMPLX(1.0D0,0.0D0)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! COMPUTATION OF SOME INTERMEDIATES TO SIMPLIFY FORMULAE
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   IF (OPT_INC.NE.0.0D0)  THEN
!
      ZGA(IFR) = (-(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0))) &
                       /(2.0D0*ZVP*ZVS*ZP*Z_XI*EXP(+ZI*W*Z_XI*(DGDEPTH(0)-TI0)))
!
      ZGC(IFR) = (+(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))) &
                       /(2.0D0*ZVP*ZVS*ZP*Z_XI*EXP(+ZI*W*Z_XI*(DGDEPTH(0)-TI0)))
   ELSE
!
      ZGA(IFR) = 0.0D0
      ZGC(IFR) = 0.0D0
!
   ENDIF
!
   ZGB(IFR) = EXP(-ZI*W*Z_XI*(DGDEPTH(0)-TI0))/EXP(+ZI*W*Z_XI*(DGDEPTH(0)-TI0))
!
   ZGD(IFR) =    (4.0D0*ZVP*ZVS*ZP*Z_XI*(1.0D0-2.0D0*ZVS**2*ZP**2)*EXP(-ZI*W*Z_XI *(DGDEPTH(0)-TI0))) &
   /(((1.0D0-2.0D0*ZVS**2*ZP**2)**2+4.0D0*ZVS**4*ZP**2*Z_XI*Z_ETA)*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0)))
!
   ZGE(IFR) = (((1.0D0-2.0D0*ZVS**2*ZP**2)**2-4.0D0*ZVS**4*ZP**2*Z_XI*Z_ETA)*EXP(-ZI*W*Z_ETA*(DGDEPTH(0)-TI0))) &
             /(((1.0D0-2.0D0*ZVS**2*ZP**2)**2+4.0D0*ZVS**4*ZP**2*Z_XI*Z_ETA)*EXP(+ZI*W*Z_ETA*(DGDEPTH(0)-TI0)))
!
   ZGF(IFR) = Z_INI(3,1)*(ZGA(IFR)*ZGE(IFR)+ZGC(IFR)) + Z_INI(3,2)*ZGE(IFR) + Z_INI(3,4)
!
   ZGG(IFR) = Z_INI(3,1)*(ZGA(IFR)*ZGD(IFR)+ZGB(IFR)) + Z_INI(3,2)*ZGD(IFR) + Z_INI(3,3)
!
   ZGH(IFR) = Z_INI(4,1)*(ZGA(IFR)*ZGE(IFR)+ZGC(IFR)) + Z_INI(4,2)*ZGE(IFR) + Z_INI(4,4)
!
   ZGI(IFR) = Z_INI(4,1)*(ZGA(IFR)*ZGD(IFR)+ZGB(IFR)) + Z_INI(4,2)*ZGD(IFR) + Z_INI(4,3)
!
   IF ( ZGG(IFR)                            .EQ.CMPLX(0.0D0,0.0D0)) STOP "ERROR IN PPWAVE: ZGG=0"
   IF ((ZGG(IFR)*ZGH(IFR)-ZGF(IFR)*ZGI(IFR)).EQ.CMPLX(0.0D0,0.0D0)) STOP "ERROR IN PPWAVE: ZGG*ZGH-ZGF*ZGI=0"
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF PUP AND SUP FOR FREE SURFACE
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   IF (OPT_INC.NE.0.0D0) THEN
      ZSUP_OBJ(IFR) = (ZSUP_HS(IFR)-ZGI(IFR)/ZGG(IFR)*ZPUP_HS(IFR))*(ZGG(IFR)/(ZGG(IFR)*ZGH(IFR)-ZGF(IFR)*ZGI(IFR)))
   ELSE
      ZSUP_OBJ(IFR) = CMPLX(0.0D0,0.0D0)
   ENDIF
   ZPUP_OBJ(IFR) = (ZPUP_HS(IFR)-ZSUP_OBJ(IFR)*ZGF(IFR))/ZGG(IFR)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF PDOWN AND SDOWN FOR FREE SURFACE
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   IF (OPT_INC.NE.0.0D0) THEN
      ZSDOWN_OBJ(IFR) = ZGD(IFR)*ZPUP_OBJ(IFR) + ZGE(IFR)*ZSUP_OBJ(IFR)
   ELSE
      ZSDOWN_OBJ(IFR) = CMPLX(0.0D0,0.0D0)
   ENDIF
   ZPDOWN_OBJ(IFR) = ZGA(IFR)*ZSDOWN_OBJ(IFR) + ZGB(IFR)*ZPUP_OBJ(IFR) + ZGC(IFR)*ZSUP_OBJ(IFR)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF PDOWN AND SDOWN FOR HALF SPACE
!-------------------------------------------------------------------------------------------------------------------------------------------
!
   ZPDOWN_HS(IFR) = Z_INI(1,1)*ZPDOWN_OBJ(IFR)+Z_INI(1,2)*ZSDOWN_OBJ(IFR)+Z_INI(1,3)*ZPUP_OBJ(IFR)+Z_INI(1,4)*ZSUP_OBJ(IFR)
   ZSDOWN_HS(IFR) = Z_INI(2,1)*ZPDOWN_OBJ(IFR)+Z_INI(2,2)*ZSDOWN_OBJ(IFR)+Z_INI(2,3)*ZPUP_OBJ(IFR)+Z_INI(2,4)*ZSUP_OBJ(IFR)
!
!-------------------------------------------------------------------------------------------------------------------------------------------
! CALCULUS OF DISPLACEMENT
!-------------------------------------------------------------------------------------------------------------------------------------------
!
! FOR FREE SURFACE
!
   ZUX_OBJ(IFR) = ZF_OBJ(1,1)*ZPDOWN_OBJ(IFR)+ZF_OBJ(1,2)*ZSDOWN_OBJ(IFR)+ZF_OBJ(1,3)*ZPUP_OBJ(IFR)+ZF_OBJ(1,4)*ZSUP_OBJ(IFR)
   ZUZ_OBJ(IFR) = ZF_OBJ(2,1)*ZPDOWN_OBJ(IFR)+ZF_OBJ(2,2)*ZSDOWN_OBJ(IFR)+ZF_OBJ(2,3)*ZPUP_OBJ(IFR)+ZF_OBJ(2,4)*ZSUP_OBJ(IFR)
!
! FOR HALF SPACE
!
   ZUX_HS(IFR)  = ZF_HS(1,1)*ZPDOWN_HS(IFR)+ZF_HS(1,2)*ZSDOWN_HS(IFR)+ZF_HS(1,3)*ZPUP_HS(IFR)+ZF_HS(1,4)*ZSUP_HS(IFR)
   ZUZ_HS(IFR)  = ZF_HS(2,1)*ZPDOWN_HS(IFR)+ZF_HS(2,2)*ZSDOWN_HS(IFR)+ZF_HS(2,3)*ZPUP_HS(IFR)+ZF_HS(2,4)*ZSUP_HS(IFR)
!
   ZSVSPRX(IFR) = ZUX_OBJ(IFR)/ZUX_HS(IFR)
   ZSVSPRZ(IFR) = ZUZ_OBJ(IFR)/ZUZ_HS(IFR)
ENDDO
!$OMP END PARALLEL DO
!
!*******************************************************************************************************************************************
! PARZEN'S SPECTRAL WINDOW
!*******************************************************************************************************************************************
!
! INVERSION WITH SPECTRAL RATIO
! PX OR PSX TRANSFER FUNCTION
!
IF ( (CGTYPIN.EQ."SPR") .AND. (OPT_INC.GT.1.0D-4) ) THEN
   IF ( (CGIDPPX.NE." PX") .AND. (CGIDPSX.NE."PSX") ) THEN
      DGRTPPX(:) = 0.0D0
   ELSE
      AMZUX_OBJ(:) = ABS(ZUX_OBJ(:))
      AMZUX_HS(:)  = ABS(ZUX_HS(:))
      IF (DGSPBAN.GT.0.0D0) THEN
         CALL SPEWIN(AMZUX_OBJ,DGDFREQ,DGSPBAN)
         CALL SPEWIN(AMZUX_HS ,DGDFREQ,DGSPBAN)
      ENDIF
      !$OMP PARALLEL DO
      DO I = 1,(IGNFOLD-1)
         DGRTPPX(I) = AMZUX_OBJ(I)/AMZUX_HS(I)
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
ELSE
   DGRTPPX(:) = 0.0D0
ENDIF
!
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPPX.NE." PZ") .AND. (CGIDPSX.NE."PSZ") ) THEN
   DGRTPPZ(:) = 0.0D0
ELSE
   AMZUZ_OBJ(:) = ABS(ZUZ_OBJ(:))
   AMZUZ_HS(:)  = ABS(ZUZ_HS(:))
   IF (DGSPBAN.GT.0.0D0) THEN
      CALL SPEWIN(AMZUZ_OBJ,DGDFREQ,DGSPBAN)
      CALL SPEWIN(AMZUZ_HS ,DGDFREQ,DGSPBAN)
   ENDIF
!
! INVERSION WITH SPECTRAL RATIO
! PZ OR PSZ TRANSFER FUNCTION
!
   IF ( (CGTYPIN.EQ."SPR") .OR. (CGTYPIN.EQ."SUM") ) THEN
      !$OMP PARALLEL DO
      DO I = 1,(IGNFOLD-1)
         DGRTPPZ(I) = AMZUZ_OBJ(I)/AMZUZ_HS(I)
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
!
! INVERSION WITH HV RATIO (BOTTOM) OR COMBINATION OF SEVERAL RATIOS
!
   IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
      !$OMP PARALLEL DO
      DO I = 1,(IGNFOLD-1)
         DGTRPZB(I) = AMZUZ_HS(I)
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
!
! INVERSION WITH HV RATIO (TOP) OR COMBINATION OF SEVERAL RATIOS
!
   IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
      !$OMP PARALLEL DO
      DO I = 1,(IGNFOLD-1)
         DGTRPZT(I) = AMZUZ_OBJ(I)
      ENDDO
      !$OMP END PARALLEL DO
   ENDIF
!
ENDIF
!
DEALLOCATE(ZSUP_HS,ZPUP_HS,ZSDOWN_HS,ZPDOWN_HS,ZSUP_OBJ,ZPUP_OBJ,ZSDOWN_OBJ,ZPDOWN_OBJ,ZUX_OBJ,ZUZ_OBJ,ZUX_HS,ZUZ_HS,ZSVSPRX,ZSVSPRZ       &
          ,AMZUX_HS,AMZUX_OBJ,AMZUZ_OBJ,AMZUZ_HS)
DEALLOCATE(ZGA,ZGB,ZGC,ZGD,ZGE,ZGF,ZGG,ZGH,ZGI)
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE PPWAVE
