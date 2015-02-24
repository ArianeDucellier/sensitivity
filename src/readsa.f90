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
SUBROUTINE READSA()
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGBNDPR                                                                                                          &
                         ,DGPVHVB,DGPVHVT,DGPVPPX,DGPVPPZ,DGPVPSX,DGPVPSZ,DGPVSHY,DGPVSVX,DGPVSVZ                                          &
                         ,DGSFHVB,DGSFHVT,DGSFPPX,DGSFPPZ,DGSFPSX,DGSFPSZ,DGSFSHY,DGSFSVX,DGSFSVZ                                          &
                         ,DGSUHVB,DGSUHVT,DGSUPPX,DGSUPPZ,DGSUPSX,DGSUPSZ,DGSUSHY,DGSUSVX,DGSUSVZ                                          &
                         ,DGUVHVB,DGUVHVT,DGUVPPX,DGUVPPZ,DGUVPSX,DGUVPSZ,DGUVSHY,DGUVSVX,DGUVSVZ                                          &
                         ,DGVECT1,DGVECT2                                                                                                  &
                         ,DGMNHVB,DGMNHVT,DGMNPPX,DGMNPPZ,DGMNPSX,DGMNPSZ,DGMNSHY,DGMNSVX,DGMNSVZ                                          &
                         ,DGSCHVB,DGSCHVT,DGSCPPX,DGSCPPZ,DGSCPSX,DGSCPSZ,DGSCSHY,DGSCSVX,DGSCSVZ                                          &
                         ,DGSIHVB,DGSIHVT,DGSIPPX,DGSIPPZ,DGSIPSX,DGSIPSZ,DGSISHY,DGSISVX,DGSISVZ                                          &
                         ,DGVTHVB,DGVTHVT,DGVTPPX,DGVTPPZ,DGVTPSX,DGVTPSZ,DGVTSHY,DGVTSVX,DGVTSVZ                                          &
!
                         ,IGNBPRM,IGNFOLD,IGSOBOL                                                                                          &
!
                         ,CGFINSA
!
IMPLICIT NONE
!
DOUBLE PRECISION                   :: X1,X2
!
INTEGER, DIMENSION(:), ALLOCATABLE :: SEED
INTEGER                            :: I,IOS99,ISEED,J
!
!*******************************************************************************************************************************************
! READ PARAMETERS FOR SENSITIVITY ANALYSIS
!*******************************************************************************************************************************************
!
OPEN(UNIT=99,FILE=CGFINSA,STATUS="OLD",IOSTAT=IOS99)
CALL CHKOPN(CGFINSA,IOS99)
!
READ(UNIT=99,FMT='(2I10)') IGSOBOL,ISEED
WRITE(10,'(A,I10)') "NUMBER OF SETS OF PARAMETERS: ",IGSOBOL
CALL RANDOM_SEED(SIZE=ISEED)
ALLOCATE(SEED(ISEED))
DO I = 1,ISEED
   READ(UNIT=99,FMT='(I10)') SEED(I)
ENDDO
CALL RANDOM_SEED(PUT=SEED)
!
ALLOCATE(DGVECT1(IGSOBOL,IGNBPRM),DGVECT2(IGSOBOL,IGNBPRM))
!
CLOSE(99)
!
!*******************************************************************************************************************************************
! FILL MATRICES RANDOMLY
!*******************************************************************************************************************************************
!
DO J = 1,IGNBPRM
   DO I = 1,IGSOBOL
      CALL RANDOM_NUMBER(X1)
      DGVECT1(I,J) = X1*(DGBNDPR(J,2)-DGBNDPR(J,1)) + DGBNDPR(J,1)
      CALL RANDOM_NUMBER(X2)
      DGVECT2(I,J) = X2*(DGBNDPR(J,2)-DGBNDPR(J,1)) + DGBNDPR(J,1)
   ENDDO
ENDDO
!
!*******************************************************************************************************************************************
! INITIALIZE VECTORS FOR MEAN
!*******************************************************************************************************************************************
!
ALLOCATE(DGMNHVB(IGNFOLD),DGMNHVT(IGNFOLD),DGMNPPX(IGNFOLD),DGMNPPZ(IGNFOLD),DGMNPSX(IGNFOLD),DGMNPSZ(IGNFOLD)                             &
        ,DGMNSHY(IGNFOLD),DGMNSVX(IGNFOLD),DGMNSVZ(IGNFOLD))
!
DGMNHVB(:) = 0.0D0
DGMNHVT(:) = 0.0D0
DGMNPPX(:) = 0.0D0
DGMNPPZ(:) = 0.0D0
DGMNPSX(:) = 0.0D0
DGMNPSZ(:) = 0.0D0
DGMNSHY(:) = 0.0D0
DGMNSVX(:) = 0.0D0
DGMNSVZ(:) = 0.0D0
!
!*******************************************************************************************************************************************
! INITIALIZE VECTORS FOR TOTAL VARIANCE
!*******************************************************************************************************************************************
!
ALLOCATE(DGVTHVB(IGNFOLD),DGVTHVT(IGNFOLD),DGVTPPX(IGNFOLD),DGVTPPZ(IGNFOLD),DGVTPSX(IGNFOLD),DGVTPSZ(IGNFOLD)                             &
        ,DGVTSHY(IGNFOLD),DGVTSVX(IGNFOLD),DGVTSVZ(IGNFOLD))
!
DGVTHVB(:) = 0.0D0
DGVTHVT(:) = 0.0D0
DGVTPPX(:) = 0.0D0
DGVTPPZ(:) = 0.0D0
DGVTPSX(:) = 0.0D0
DGVTPSZ(:) = 0.0D0
DGVTSHY(:) = 0.0D0
DGVTSVX(:) = 0.0D0
DGVTSVZ(:) = 0.0D0
!
!*******************************************************************************************************************************************
! INITIALIZE MATRICES FOR PARTIAL VARIANCES
!*******************************************************************************************************************************************
!
ALLOCATE(DGPVHVB(IGNBPRM,IGNFOLD),DGPVHVT(IGNBPRM,IGNFOLD),DGPVPPX(IGNBPRM,IGNFOLD),DGPVPPZ(IGNBPRM,IGNFOLD)                               &
        ,DGPVPSX(IGNBPRM,IGNFOLD),DGPVPSZ(IGNBPRM,IGNFOLD),DGPVSHY(IGNBPRM,IGNFOLD),DGPVSVX(IGNBPRM,IGNFOLD),DGPVSVZ(IGNBPRM,IGNFOLD))
ALLOCATE(DGUVHVB(IGNBPRM,IGNFOLD),DGUVHVT(IGNBPRM,IGNFOLD),DGUVPPX(IGNBPRM,IGNFOLD),DGUVPPZ(IGNBPRM,IGNFOLD)                               &
        ,DGUVPSX(IGNBPRM,IGNFOLD),DGUVPSZ(IGNBPRM,IGNFOLD),DGUVSHY(IGNBPRM,IGNFOLD),DGUVSVX(IGNBPRM,IGNFOLD),DGUVSVZ(IGNBPRM,IGNFOLD))
!
DGPVHVB(:,:) = 0.0D0
DGPVHVT(:,:) = 0.0D0
DGPVPPX(:,:) = 0.0D0
DGPVPPZ(:,:) = 0.0D0
DGPVPSX(:,:) = 0.0D0
DGPVPSZ(:,:) = 0.0D0
DGPVSHY(:,:) = 0.0D0
DGPVSVX(:,:) = 0.0D0
DGPVSVZ(:,:) = 0.0D0
!
DGUVHVB(:,:) = 0.0D0
DGUVHVT(:,:) = 0.0D0
DGUVPPX(:,:) = 0.0D0
DGUVPPZ(:,:) = 0.0D0
DGUVPSX(:,:) = 0.0D0
DGUVPSZ(:,:) = 0.0D0
DGUVSHY(:,:) = 0.0D0
DGUVSVX(:,:) = 0.0D0
DGUVSVZ(:,:) = 0.0D0
!
!*******************************************************************************************************************************************
! INITIALIZE MATRICES FOR SOBOL COEFFICIENTS
!*******************************************************************************************************************************************
!
ALLOCATE(DGSFHVB(IGNBPRM,IGNFOLD),DGSFHVT(IGNBPRM,IGNFOLD),DGSFPPX(IGNBPRM,IGNFOLD),DGSFPPZ(IGNBPRM,IGNFOLD)                               &
        ,DGSFPSX(IGNBPRM,IGNFOLD),DGSFPSZ(IGNBPRM,IGNFOLD),DGSFSHY(IGNBPRM,IGNFOLD),DGSFSVX(IGNBPRM,IGNFOLD),DGSFSVZ(IGNBPRM,IGNFOLD))
ALLOCATE(DGSUHVB(IGNBPRM,IGNFOLD),DGSUHVT(IGNBPRM,IGNFOLD),DGSUPPX(IGNBPRM,IGNFOLD),DGSUPPZ(IGNBPRM,IGNFOLD)                               &
        ,DGSUPSX(IGNBPRM,IGNFOLD),DGSUPSZ(IGNBPRM,IGNFOLD),DGSUSHY(IGNBPRM,IGNFOLD),DGSUSVX(IGNBPRM,IGNFOLD),DGSUSVZ(IGNBPRM,IGNFOLD))
!
DGSFHVB(:,:) = 0.0D0
DGSFHVT(:,:) = 0.0D0
DGSFPPX(:,:) = 0.0D0
DGSFPPZ(:,:) = 0.0D0
DGSFPSX(:,:) = 0.0D0
DGSFPSZ(:,:) = 0.0D0
DGSFSHY(:,:) = 0.0D0
DGSFSVX(:,:) = 0.0D0
DGSFSVZ(:,:) = 0.0D0
!
DGSUHVB(:,:) = 0.0D0
DGSUHVT(:,:) = 0.0D0
DGSUPPX(:,:) = 0.0D0
DGSUPPZ(:,:) = 0.0D0
DGSUPSX(:,:) = 0.0D0
DGSUPSZ(:,:) = 0.0D0
DGSUSHY(:,:) = 0.0D0
DGSUSVX(:,:) = 0.0D0
DGSUSVZ(:,:) = 0.0D0
!
ALLOCATE(DGSCHVB(IGNBPRM),DGSCHVT(IGNBPRM),DGSCPPX(IGNBPRM),DGSCPPZ(IGNBPRM),DGSCPSX(IGNBPRM),DGSCPSZ(IGNBPRM)                             &
        ,DGSCSHY(IGNBPRM),DGSCSVX(IGNBPRM),DGSCSVZ(IGNBPRM))
ALLOCATE(DGSIHVB(IGNBPRM),DGSIHVT(IGNBPRM),DGSIPPX(IGNBPRM),DGSIPPZ(IGNBPRM),DGSIPSX(IGNBPRM),DGSIPSZ(IGNBPRM)                             &
        ,DGSISHY(IGNBPRM),DGSISVX(IGNBPRM),DGSISVZ(IGNBPRM))
!
DGSCHVB(:) = 0.0D0
DGSCHVT(:) = 0.0D0
DGSCPPX(:) = 0.0D0
DGSCPPZ(:) = 0.0D0
DGSCPSX(:) = 0.0D0
DGSCPSZ(:) = 0.0D0
DGSCSHY(:) = 0.0D0
DGSCSVX(:) = 0.0D0
DGSCSVZ(:) = 0.0D0
!
DGSIHVB(:) = 0.0D0
DGSIHVT(:) = 0.0D0
DGSIPPX(:) = 0.0D0
DGSIPPZ(:) = 0.0D0
DGSIPSX(:) = 0.0D0
DGSIPSZ(:) = 0.0D0
DGSISHY(:) = 0.0D0
DGSISVX(:) = 0.0D0
DGSISVZ(:) = 0.0D0
!
!*******************************************************************************************************************************************
!
RETURN
END SUBROUTINE READSA
