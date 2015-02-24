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
MODULE AAMODU_GLOBVA
!
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: DGBNDPR,DGPRINI                                                                           &
                                                ,DGPVHVB,DGPVHVT,DGPVPPX,DGPVPPZ,DGPVPSX,DGPVPSZ,DGPVSHY,DGPVSVX,DGPVSVZ                   &
                                                ,DGSFHVB,DGSFHVT,DGSFPPX,DGSFPPZ,DGSFPSX,DGSFPSZ,DGSFSHY,DGSFSVX,DGSFSVZ                   &
                                                ,DGSUHVB,DGSUHVT,DGSUPPX,DGSUPPZ,DGSUPSX,DGSUPSZ,DGSUSHY,DGSUSVX,DGSUSVZ                   &
                                                ,DGUVHVB,DGUVHVT,DGUVPPX,DGUVPPZ,DGUVPSX,DGUVPSZ,DGUVSHY,DGUVSVX,DGUVSVZ                   &
                                                ,DGVECT1,DGVECT2
!
DOUBLE PRECISION, DIMENSION(:)  , ALLOCATABLE :: DGALPHA,DGBETA_,DGDENSI,DGDEPTH,DGGAMMA,DGINANG                                           &
                                                ,DGINVAL,DGINVBE,DGINVGA,DGINVIN,DGINVNU,DGINVPS,DGINVTH,DGINVVP,DGINVVS                   &
                                                ,DGLAHEI                                                                                   &
                                                ,DGMNHVB,DGMNHVT,DGMNPPX,DGMNPPZ,DGMNPSX,DGMNPSZ,DGMNSHY,DGMNSVX,DGMNSVZ                   &
                                                ,DGPOISR,DGPVELO                                                                           &
                                                ,DGRTHVB,DGRTHVT,DGRTPPX,DGRTPPZ,DGRTPSX,DGRTPSZ,DGRTSHY,DGRTSVX,DGRTSVZ                   &
                                                ,DGSCHVB,DGSCHVT,DGSCPPX,DGSCPPZ,DGSCPSX,DGSCPSZ,DGSCSHY,DGSCSVX,DGSCSVZ                   &
                                                ,DGSHEMO                                                                                   &
                                                ,DGSIHVB,DGSIHVT,DGSIPPX,DGSIPPZ,DGSIPSX,DGSIPSZ,DGSISHY,DGSISVX,DGSISVZ                   &
                                                ,DGSVELO                                                                                   &
                                                ,DGTRPZB,DGTRPZT,DGTRSHB,DGTRSHT                                                           &
                                                ,DGVTHVB,DGVTHVT,DGVTPPX,DGVTPPZ,DGVTPSX,DGVTPSZ,DGVTSHY,DGVTSVX,DGVTSVZ                   &
                                                ,DGUWEIG,DGWEIPS
!
DOUBLE PRECISION                              :: DGDELTT,DGDFREQ                                                                           &
                                                ,DGMMHVB,DGMMHVT,DGMMPPX,DGMMPPZ,DGMMPSX,DGMMPSZ,DGMMSHY,DGMMSVX,DGMMSVZ                   &
                                                ,DGSPBAN                                                                                   &
                                                ,DGVMHVB,DGVMHVT,DGVMPPX,DGVMPPZ,DGVMPSX,DGVMPSZ,DGVMSHY,DGVMSVX,DGVMSVZ
!
INTEGER         , DIMENSION(:,:), ALLOCATABLE :: IGLAYPR
!
INTEGER                                       :: IGIDDMP,IGINDPR,IGLAOBJ,IGLAREF,IGNBPRM,IGNFOLD,IGNLAYE,IGNT___,IGNTHRD                   &
                                                ,IGSOBOL
!
CHARACTER(10)   , DIMENSION(:)  , ALLOCATABLE :: CGPRTYP
!
CHARACTER(3)                                  :: CGIDPPX,CGIDPPZ,CGIDPSX,CGIDPSZ,CGIDSHY,CGIDSVX,CGIDSVZ,CGTYPIN
!
CHARACTER(10)                                 :: CGPREFI
!
CHARACTER(92)                                 :: CGFINPU,CGFINSA,CGFLLOG                                                                   &
                                                ,CGSFHVB,CGSFHVT,CGSFPPX,CGSFPPZ,CGSFPSX,CGSFPSZ,CGSFSHY,CGSFSVX,CGSFSVZ                   &
                                                ,CGSOBCO
!
END MODULE AAMODU_GLOBVA
