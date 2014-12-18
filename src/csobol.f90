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
PROGRAM CSOBOL
!
!*******************************************************************************************************************************************
!
USE AAMODU_GLOBVA, ONLY : DGPVHVB,DGPVHVT,DGPVPPX,DGPVPPZ,DGPVPSX,DGPVPSZ,DGPVSHY,DGPVSVX,DGPVSVZ                                          &
                         ,DGSFHVB,DGSFHVT,DGSFPPX,DGSFPPZ,DGSFPSX,DGSFPSZ,DGSFSHY,DGSFSVX,DGSFSVZ                                          &
                         ,DGSUHVB,DGSUHVT,DGSUPPX,DGSUPPZ,DGSUPSX,DGSUPSZ,DGSUSHY,DGSUSVX,DGSUSVZ                                          &
                         ,DGUVHVB,DGUVHVT,DGUVPPX,DGUVPPZ,DGUVPSX,DGUVPSZ,DGUVSHY,DGUVSVX,DGUVSVZ                                          &
                         ,DGMNHVB,DGMNHVT,DGMNPPX,DGMNPPZ,DGMNPSX,DGMNPSZ,DGMNSHY,DGMNSVX,DGMNSVZ                                          &
                         ,DGRTHVB,DGRTHVT,DGRTPPX,DGRTPPZ,DGRTPSX,DGRTPSZ,DGRTSHY,DGRTSVX,DGRTSVZ                                          &
                         ,DGSCHVB,DGSCHVT,DGSCPPX,DGSCPPZ,DGSCPSX,DGSCPSZ,DGSCSHY,DGSCSVX,DGSCSVZ                                          &
                         ,DGSIHVB,DGSIHVT,DGSIPPX,DGSIPPZ,DGSIPSX,DGSIPSZ,DGSISHY,DGSISVX,DGSISVZ                                          &
                         ,DGVTHVB,DGVTHVT,DGVTPPX,DGVTPPZ,DGVTPSX,DGVTPSZ,DGVTSHY,DGVTSVX,DGVTSVZ                                          &
                         ,DGDFREQ                                                                                                          &
                         ,DGMMHVB,DGMMHVT,DGMMPPX,DGMMPPZ,DGMMPSX,DGMMPSZ,DGMMSHY,DGMMSVX,DGMMSVZ                                          &
                         ,DGVMHVB,DGVMHVT,DGVMPPX,DGVMPPZ,DGVMPSX,DGVMPSZ,DGVMSHY,DGVMSVX,DGVMSVZ                                          &
!
                         ,IGLAYPR                                                                                                          &
                         ,IGINDPR,IGNBPRM,IGNFOLD,IGSOBOL                                                                                  &
!
                         ,CGPRTYP                                                                                                          &
                         ,CGIDPPX,CGIDPPZ,CGIDPSX,CGIDPSZ,CGIDSHY,CGIDSVX,CGIDSVZ,CGTYPIN                                                  &
                         ,CGPREFI                                                                                                          &
                         ,CGFINPU,CGFINSA,CGFLLOG                                                                                          &
                         ,CGSFHVB,CGSFHVT,CGSFPPX,CGSFPPZ,CGSFPSX,CGSFPSZ,CGSFSHY,CGSFSVX,CGSFSVZ                                          &
                         ,CGSOBCO
!
IMPLICIT NONE
!
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DLR1HVB,DLR1HVT,DLR1PPX,DLR1PPZ,DLR1PSX,DLR1PSZ,DLR1SHY,DLR1SVX,DLR1SVZ
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: DLR2HVB,DLR2HVT,DLR2PPX,DLR2PPZ,DLR2PSX,DLR2PSZ,DLR2SHY,DLR2SVX,DLR2SVZ
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TIME4
DOUBLE PRECISION                            :: TIME_TEMP1,TIME_TEMP2,TIME_TEMP3,TIME_TEMP4,TIME_TEMP5,TIME1,TIME2,TIME3,TIME5,TIME6,TIME7
!
INTEGER                                     :: I,ILSOBOL,IW,J
!
!*******************************************************************************************************************************************
! READ PREFIX - DEFINE FILE NAMES
!*******************************************************************************************************************************************
!
CALL CPU_TIME(TIME1)
!
CALL READPF(CGPREFI)
CGFINPU = TRIM(ADJUSTL(CGPREFI))//".in"
CGFINSA = TRIM(ADJUSTL(CGPREFI))//".in.sobol"
CGFLLOG = TRIM(ADJUSTL(CGPREFI))//".log"
CGSFHVB = TRIM(ADJUSTL(CGPREFI))//".sobol.hvb"
CGSFHVT = TRIM(ADJUSTL(CGPREFI))//".sobol.hvt"
CGSFPPX = TRIM(ADJUSTL(CGPREFI))//".sobol.ppx"
CGSFPPZ = TRIM(ADJUSTL(CGPREFI))//".sobol.ppz"
CGSFPSX = TRIM(ADJUSTL(CGPREFI))//".sobol.psx"
CGSFPSZ = TRIM(ADJUSTL(CGPREFI))//".sobol.psz"
CGSFSHY = TRIM(ADJUSTL(CGPREFI))//".sobol.shy"
CGSFSVX = TRIM(ADJUSTL(CGPREFI))//".sobol.svx"
CGSFSVZ = TRIM(ADJUSTL(CGPREFI))//".sobol.svz"
CGSOBCO = TRIM(ADJUSTL(CGPREFI))//".sobol.coeff"
OPEN(UNIT=10,FILE=CGFLLOG)
!
!*******************************************************************************************************************************************
! READ CONTROL FILE
!*******************************************************************************************************************************************
!
WRITE(10,'(A)') "****************************************"
WRITE(10,'(A)') " READING FILE control_file_genetic_algo"
WRITE(10,'(A)') "****************************************"
CALL READCF()
!
!*******************************************************************************************************************************************
! READ INITIAL PARAMETERS OF THE SOIL COLUMN
!*******************************************************************************************************************************************
!
WRITE(10,'(A)') "**************************"
WRITE(10,'(A)') " READING SOIL COLUMN FILE"
WRITE(10,'(A)') "**************************"
CALL READSC()
!
!*******************************************************************************************************************************************
! READ PARAMETERS FOR SENSITIVITY ANALYSIS
!*******************************************************************************************************************************************
!
WRITE(10,'(A)') "*******************************"
WRITE(10,'(A)') " READING SOBOL PARAMETERS FILE"
WRITE(10,'(A)') "*******************************"
CALL READSA()
!
!*******************************************************************************************************************************************
! INITIALIZATION OF VECTORS PARAMETERS
!*******************************************************************************************************************************************
!
CALL INIVCO()
!
ALLOCATE(DLR1HVB(IGNFOLD),DLR1HVT(IGNFOLD),DLR1PPX(IGNFOLD),DLR1PPZ(IGNFOLD),DLR1PSX(IGNFOLD),DLR1PSZ(IGNFOLD)                             &
        ,DLR1SHY(IGNFOLD),DLR1SVX(IGNFOLD),DLR1SVZ(IGNFOLD))
ALLOCATE(DLR2HVB(IGNFOLD),DLR2HVT(IGNFOLD),DLR2PPX(IGNFOLD),DLR2PPZ(IGNFOLD),DLR2PSX(IGNFOLD),DLR2PSZ(IGNFOLD)                             &
        ,DLR2SHY(IGNFOLD),DLR2SVX(IGNFOLD),DLR2SVZ(IGNFOLD))
!
DLR1HVB(:) = 0.0D0
DLR1HVT(:) = 0.0D0
DLR1PPX(:) = 0.0D0
DLR1PPZ(:) = 0.0D0
DLR1PSX(:) = 0.0D0
DLR1PSZ(:) = 0.0D0
DLR1SHY(:) = 0.0D0
DLR1SVX(:) = 0.0D0
DLR1SVZ(:) = 0.0D0
!
DLR2HVB(:) = 0.0D0
DLR2HVT(:) = 0.0D0
DLR2PPX(:) = 0.0D0
DLR2PPZ(:) = 0.0D0
DLR2PSX(:) = 0.0D0
DLR2PSZ(:) = 0.0D0
DLR2SHY(:) = 0.0D0
DLR2SVX(:) = 0.0D0
DLR2SVZ(:) = 0.0D0
!
!*******************************************************************************************************************************************
! BEGINNING OF LOOP
!*******************************************************************************************************************************************
!
CALL CPU_TIME(TIME2)
!
TIME3 = 0.0D0
ALLOCATE(TIME4(IGNBPRM))
TIME4(:) = 0.0D0
TIME5 = 0.0D0
!
DO ILSOBOL = 1,IGSOBOL
!
   CALL CPU_TIME(TIME_TEMP1)
!
   CALL FCTOVL(ILSOBOL,1,0)
   CALL TRFUNC(ILSOBOL)
!
   DLR1HVB(:) = DGRTHVB(:)
   DLR1HVT(:) = DGRTHVT(:)
   DLR1PPX(:) = DGRTPPX(:)
   DLR1PPZ(:) = DGRTPPZ(:)
   DLR1PSX(:) = DGRTPSX(:)
   DLR1PSZ(:) = DGRTPSZ(:)
   DLR1SHY(:) = DGRTSHY(:)
   DLR1SVX(:) = DGRTSVX(:)
   DLR1SVZ(:) = DGRTSVZ(:)
!
!*******************************************************************************************************************************************
! COMPUTATION OF MEAN
!*******************************************************************************************************************************************
!
   DO IW = 1,(IGNFOLD-1)
      DGMNHVB(IW) = DGMNHVB(IW) + DLR1HVB(IW)/IGSOBOL
      DGMNHVT(IW) = DGMNHVT(IW) + DLR1HVT(IW)/IGSOBOL
      DGMNPPX(IW) = DGMNPPX(IW) + DLR1PPX(IW)/IGSOBOL
      DGMNPPZ(IW) = DGMNPPZ(IW) + DLR1PPZ(IW)/IGSOBOL
      DGMNPSX(IW) = DGMNPSX(IW) + DLR1PSX(IW)/IGSOBOL
      DGMNPSZ(IW) = DGMNPSZ(IW) + DLR1PSZ(IW)/IGSOBOL
      DGMNSHY(IW) = DGMNSHY(IW) + DLR1SHY(IW)/IGSOBOL
      DGMNSVX(IW) = DGMNSVX(IW) + DLR1SVX(IW)/IGSOBOL
      DGMNSVZ(IW) = DGMNSVZ(IW) + DLR1SVZ(IW)/IGSOBOL
   ENDDO
!
!*******************************************************************************************************************************************
! COMPUTATION OF TOTAL VARIANCE
!*******************************************************************************************************************************************
!
   DO IW = 1,(IGNFOLD-1)
      DGVTHVB(IW) = DGVTHVB(IW) + DLR1HVB(IW)*DLR1HVB(IW)/(IGSOBOL-1)
      DGVTHVT(IW) = DGVTHVT(IW) + DLR1HVT(IW)*DLR1HVT(IW)/(IGSOBOL-1)
      DGVTPPX(IW) = DGVTPPX(IW) + DLR1PPX(IW)*DLR1PPX(IW)/(IGSOBOL-1)
      DGVTPPZ(IW) = DGVTPPZ(IW) + DLR1PPZ(IW)*DLR1PPZ(IW)/(IGSOBOL-1)
      DGVTPSX(IW) = DGVTPSX(IW) + DLR1PSX(IW)*DLR1PSX(IW)/(IGSOBOL-1)
      DGVTPSZ(IW) = DGVTPSZ(IW) + DLR1PSZ(IW)*DLR1PSZ(IW)/(IGSOBOL-1)
      DGVTSHY(IW) = DGVTSHY(IW) + DLR1SHY(IW)*DLR1SHY(IW)/(IGSOBOL-1)
      DGVTSVX(IW) = DGVTSVX(IW) + DLR1SVX(IW)*DLR1SVX(IW)/(IGSOBOL-1)
      DGVTSVZ(IW) = DGVTSVZ(IW) + DLR1SVZ(IW)*DLR1SVZ(IW)/(IGSOBOL-1)
   ENDDO
!
   CALL CPU_TIME(TIME_TEMP2)
!
   CALL FCTOVL(ILSOBOL,2,0)
   CALL TRFUNC(ILSOBOL)
   DLR2HVB(:) = DGRTHVB(:)
   DLR2HVT(:) = DGRTHVT(:)
   DLR2PPX(:) = DGRTPPX(:)
   DLR2PPZ(:) = DGRTPPZ(:)
   DLR2PSX(:) = DGRTPSX(:)
   DLR2PSZ(:) = DGRTPSZ(:)
   DLR2SHY(:) = DGRTSHY(:)
   DLR2SVX(:) = DGRTSVX(:)
   DLR2SVZ(:) = DGRTSVZ(:)
!
   DO IGINDPR = 1,IGNBPRM
!
      CALL CPU_TIME(TIME_TEMP3)
!
      CALL FCTOVL(ILSOBOL,3,IGINDPR)
      CALL TRFUNC(ILSOBOL)
!
!*******************************************************************************************************************************************
! COMPUTATION OF FIRST ORDER PARTIAL VARIANCES
!*******************************************************************************************************************************************
!
      DO IW = 1,(IGNFOLD-1)
         DGUVHVB(IGINDPR,IW) = DGUVHVB(IGINDPR,IW) + DLR2HVB(IW)*(DGRTHVB(IW)-DLR1HVB(IW))/(IGSOBOL-1)
         DGUVHVT(IGINDPR,IW) = DGUVHVT(IGINDPR,IW) + DLR2HVT(IW)*(DGRTHVT(IW)-DLR1HVT(IW))/(IGSOBOL-1)
         DGUVPPX(IGINDPR,IW) = DGUVPPX(IGINDPR,IW) + DLR2PPX(IW)*(DGRTPPX(IW)-DLR1PPX(IW))/(IGSOBOL-1)
         DGUVPPZ(IGINDPR,IW) = DGUVPPZ(IGINDPR,IW) + DLR2PPZ(IW)*(DGRTPPZ(IW)-DLR1PPZ(IW))/(IGSOBOL-1)
         DGUVPSX(IGINDPR,IW) = DGUVPSX(IGINDPR,IW) + DLR2PSX(IW)*(DGRTPSX(IW)-DLR1PSX(IW))/(IGSOBOL-1)
         DGUVPSZ(IGINDPR,IW) = DGUVPSZ(IGINDPR,IW) + DLR2PSZ(IW)*(DGRTPSZ(IW)-DLR1PSZ(IW))/(IGSOBOL-1)
         DGUVSHY(IGINDPR,IW) = DGUVSHY(IGINDPR,IW) + DLR2SHY(IW)*(DGRTSHY(IW)-DLR1SHY(IW))/(IGSOBOL-1)
         DGUVSVX(IGINDPR,IW) = DGUVSVX(IGINDPR,IW) + DLR2SVX(IW)*(DGRTSVX(IW)-DLR1SVX(IW))/(IGSOBOL-1)
         DGUVSVZ(IGINDPR,IW) = DGUVSVZ(IGINDPR,IW) + DLR2SVZ(IW)*(DGRTSVZ(IW)-DLR1SVZ(IW))/(IGSOBOL-1)
      ENDDO
!
!*******************************************************************************************************************************************
! COMPUTATION OF TOTAL PARTIAL VARIANCES
!*******************************************************************************************************************************************
!
      DO IW = 1,(IGNFOLD-1)
         DGPVHVB(IGINDPR,IW) = DGPVHVB(IGINDPR,IW) + ((DGRTHVB(IW)-DLR1HVB(IW))**2)/(2*(IGSOBOL-1))
         DGPVHVT(IGINDPR,IW) = DGPVHVT(IGINDPR,IW) + ((DGRTHVT(IW)-DLR1HVT(IW))**2)/(2*(IGSOBOL-1))
         DGPVPPX(IGINDPR,IW) = DGPVPPX(IGINDPR,IW) + ((DGRTPPX(IW)-DLR1PPX(IW))**2)/(2*(IGSOBOL-1))
         DGPVPPZ(IGINDPR,IW) = DGPVPPZ(IGINDPR,IW) + ((DGRTPPZ(IW)-DLR1PPZ(IW))**2)/(2*(IGSOBOL-1))
         DGPVPSX(IGINDPR,IW) = DGPVPSX(IGINDPR,IW) + ((DGRTPSX(IW)-DLR1PSX(IW))**2)/(2*(IGSOBOL-1))
         DGPVPSZ(IGINDPR,IW) = DGPVPSZ(IGINDPR,IW) + ((DGRTPSZ(IW)-DLR1PSZ(IW))**2)/(2*(IGSOBOL-1))
         DGPVSHY(IGINDPR,IW) = DGPVSHY(IGINDPR,IW) + ((DGRTSHY(IW)-DLR1SHY(IW))**2)/(2*(IGSOBOL-1))
         DGPVSVX(IGINDPR,IW) = DGPVSVX(IGINDPR,IW) + ((DGRTSVX(IW)-DLR1SVX(IW))**2)/(2*(IGSOBOL-1))
         DGPVSVZ(IGINDPR,IW) = DGPVSVZ(IGINDPR,IW) + ((DGRTSVZ(IW)-DLR1SVZ(IW))**2)/(2*(IGSOBOL-1))
      ENDDO
!
      CALL CPU_TIME(TIME_TEMP4)
      TIME4(IGINDPR) = TIME4(IGINDPR) + (TIME_TEMP4-TIME_TEMP3)
!
   ENDDO
!
   CALL CPU_TIME(TIME_TEMP5)
   TIME3 = TIME3 + (TIME_TEMP2-TIME_TEMP1)
   TIME5 = TIME5 + (TIME_TEMP5-TIME_TEMP2)
!
ENDDO
!
DO IW = 1,(IGNFOLD-1)
   DGVTHVB(IW) = DGVTHVB(IW) - DGMNHVB(IW)*DGMNHVB(IW)*(1+1/(IGSOBOL-1))
   DGVTHVT(IW) = DGVTHVT(IW) - DGMNHVT(IW)*DGMNHVT(IW)*(1+1/(IGSOBOL-1))
   DGVTPPX(IW) = DGVTPPX(IW) - DGMNPPX(IW)*DGMNPPX(IW)*(1+1/(IGSOBOL-1))
   DGVTPPZ(IW) = DGVTPPZ(IW) - DGMNPPZ(IW)*DGMNPPZ(IW)*(1+1/(IGSOBOL-1))
   DGVTPSX(IW) = DGVTPSX(IW) - DGMNPSX(IW)*DGMNPSX(IW)*(1+1/(IGSOBOL-1))
   DGVTPSZ(IW) = DGVTPSZ(IW) - DGMNPSZ(IW)*DGMNPSZ(IW)*(1+1/(IGSOBOL-1))
   DGVTSHY(IW) = DGVTSHY(IW) - DGMNSHY(IW)*DGMNSHY(IW)*(1+1/(IGSOBOL-1))
   DGVTSVX(IW) = DGVTSVX(IW) - DGMNSVX(IW)*DGMNSVX(IW)*(1+1/(IGSOBOL-1))
   DGVTSVZ(IW) = DGVTSVZ(IW) - DGMNSVZ(IW)*DGMNSVZ(IW)*(1+1/(IGSOBOL-1))
ENDDO
!
!*******************************************************************************************************************************************
! COMPUTE AND WRITE SOBOL INDICES
!*******************************************************************************************************************************************
!
CALL CPU_TIME(TIME6)
!
WRITE(10,'(A)') "***********************************"
WRITE(10,'(A)') " COMPUTATION OF SOBOL COEFFICIENTS"
WRITE(10,'(A)') "***********************************"
!
IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
   OPEN(UNIT=21,FILE=CGSFHVB)
ENDIF
IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
   OPEN(UNIT=22,FILE=CGSFHVT)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPPX.EQ." PX") ) THEN
   OPEN(UNIT=23,FILE=CGSFPPX)
ENDIF
IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDPPZ.EQ." PZ")) .OR. (CGTYPIN.EQ."SUM") ) THEN
   OPEN(UNIT=24,FILE=CGSFPPZ)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSX.EQ."PSX") ) THEN
   OPEN(UNIT=25,FILE=CGSFPSX)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSZ.EQ."PSZ") ) THEN
   OPEN(UNIT=26,FILE=CGSFPSZ)
ENDIF
IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDSHY.EQ." SH")) .OR. (CGTYPIN.EQ."SUM") ) THEN
   OPEN(UNIT=27,FILE=CGSFSHY)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVX.EQ."SVX") ) THEN
   OPEN(UNIT=28,FILE=CGSFSVX)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVZ.EQ."SVZ") ) THEN
   OPEN(UNIT=29,FILE=CGSFSVZ)
ENDIF
!
DO IW = 1,(IGNFOLD-1)
   DO IGINDPR = 1,IGNBPRM
      DGSUHVB(IGINDPR,IW) = DGUVHVB(IGINDPR,IW)/DGVTHVB(IW)
      DGSUHVT(IGINDPR,IW) = DGUVHVT(IGINDPR,IW)/DGVTHVT(IW)
      DGSUPPX(IGINDPR,IW) = DGUVPPX(IGINDPR,IW)/DGVTPPX(IW)
      DGSUPPZ(IGINDPR,IW) = DGUVPPZ(IGINDPR,IW)/DGVTPPZ(IW)
      DGSUPSX(IGINDPR,IW) = DGUVPSX(IGINDPR,IW)/DGVTPSX(IW)
      DGSUPSZ(IGINDPR,IW) = DGUVPSZ(IGINDPR,IW)/DGVTPSZ(IW)
      DGSUSHY(IGINDPR,IW) = DGUVSHY(IGINDPR,IW)/DGVTSHY(IW)
      DGSUSVX(IGINDPR,IW) = DGUVSVX(IGINDPR,IW)/DGVTSVX(IW)
      DGSUSVZ(IGINDPR,IW) = DGUVSVZ(IGINDPR,IW)/DGVTSVZ(IW)
!
      DGSFHVB(IGINDPR,IW) = DGPVHVB(IGINDPR,IW)/DGVTHVB(IW)
      DGSFHVT(IGINDPR,IW) = DGPVHVT(IGINDPR,IW)/DGVTHVT(IW)
      DGSFPPX(IGINDPR,IW) = DGPVPPX(IGINDPR,IW)/DGVTPPX(IW)
      DGSFPPZ(IGINDPR,IW) = DGPVPPZ(IGINDPR,IW)/DGVTPPZ(IW)
      DGSFPSX(IGINDPR,IW) = DGPVPSX(IGINDPR,IW)/DGVTPSX(IW)
      DGSFPSZ(IGINDPR,IW) = DGPVPSZ(IGINDPR,IW)/DGVTPSZ(IW)
      DGSFSHY(IGINDPR,IW) = DGPVSHY(IGINDPR,IW)/DGVTSHY(IW)
      DGSFSVX(IGINDPR,IW) = DGPVSVX(IGINDPR,IW)/DGVTSVX(IW)
      DGSFSVZ(IGINDPR,IW) = DGPVSVZ(IGINDPR,IW)/DGVTSVZ(IW)
   ENDDO
!
   IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
      WRITE(UNIT=21,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNHVB(IW),DGVTHVB(IW),(DGSUHVB(I,IW),I=1,IGNBPRM),(DGSFHVB(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
      WRITE(UNIT=22,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNHVT(IW),DGVTHVT(IW),(DGSUHVT(I,IW),I=1,IGNBPRM),(DGSFHVT(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPPX.EQ." PX") ) THEN
      WRITE(UNIT=23,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNPPX(IW),DGVTPPX(IW),(DGSUPPX(I,IW),I=1,IGNBPRM),(DGSFPPX(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDPPZ.EQ." PZ")) .OR. (CGTYPIN.EQ."SUM") ) THEN
      WRITE(UNIT=24,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNPPZ(IW),DGVTPPZ(IW),(DGSUPPZ(I,IW),I=1,IGNBPRM),(DGSFPPZ(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSX.EQ."PSX") ) THEN
      WRITE(UNIT=25,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNPSX(IW),DGVTPSX(IW),(DGSUPSX(I,IW),I=1,IGNBPRM),(DGSFPSX(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSZ.EQ."PSZ") ) THEN
      WRITE(UNIT=26,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNPSZ(IW),DGVTPSZ(IW),(DGSUPSZ(I,IW),I=1,IGNBPRM),(DGSFPSZ(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDSHY.EQ." SH")) .OR. (CGTYPIN.EQ."SUM") ) THEN
      WRITE(UNIT=27,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNSHY(IW),DGVTSHY(IW),(DGSUSHY(I,IW),I=1,IGNBPRM),(DGSFSHY(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVX.EQ."SVX") ) THEN
      WRITE(UNIT=28,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNSVX(IW),DGVTSVX(IW),(DGSUSVX(I,IW),I=1,IGNBPRM),(DGSFSVX(J,IW),J=1,IGNBPRM)
   ENDIF
   IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVZ.EQ."SVZ") ) THEN
      WRITE(UNIT=29,FMT='(F10.4,2(E15.7E3),1000(1X,F10.4))') IW*DGDFREQ,DGMNSVZ(IW),DGVTSVZ(IW),(DGSUSVZ(I,IW),I=1,IGNBPRM),(DGSFSVZ(J,IW),J=1,IGNBPRM)
   ENDIF
ENDDO
!
IF ( (CGTYPIN.EQ."HVB") .OR. (CGTYPIN.EQ."SUM") ) THEN
   CLOSE(21)
ENDIF
IF ( (CGTYPIN.EQ."HVT") .OR. (CGTYPIN.EQ."SUM") ) THEN
   CLOSE(22)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPPX.EQ." PX") ) THEN
   CLOSE(23)
ENDIF
IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDPPZ.EQ." PZ")) .OR. (CGTYPIN.EQ."SUM") ) THEN
   CLOSE(24)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSX.EQ."PSX") ) THEN
   CLOSE(25)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDPSZ.EQ."PSZ") ) THEN
   CLOSE(26)
ENDIF
IF ( ((CGTYPIN.EQ."SPR") .AND. (CGIDSHY.EQ." SH")) .OR. (CGTYPIN.EQ."SUM") ) THEN
   CLOSE(27)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVX.EQ."SVX") ) THEN
   CLOSE(28)
ENDIF
IF ( (CGTYPIN.EQ."SPR") .AND. (CGIDSVZ.EQ."SVZ") ) THEN
   CLOSE(29)
ENDIF
!
OPEN(UNIT=20,FILE=CGSOBCO)
!
DGMMHVB = SUM(DGMNHVB(:))/DBLE(IGNFOLD-1)
DGMMHVT = SUM(DGMNHVT(:))/DBLE(IGNFOLD-1)
DGMMPPX = SUM(DGMNPPX(:))/DBLE(IGNFOLD-1)
DGMMPPZ = SUM(DGMNPPZ(:))/DBLE(IGNFOLD-1)
DGMMPSX = SUM(DGMNPSX(:))/DBLE(IGNFOLD-1)
DGMMPSZ = SUM(DGMNPSZ(:))/DBLE(IGNFOLD-1)
DGMMSHY = SUM(DGMNSHY(:))/DBLE(IGNFOLD-1)
DGMMSVX = SUM(DGMNSVX(:))/DBLE(IGNFOLD-1)
DGMMSVZ = SUM(DGMNSVZ(:))/DBLE(IGNFOLD-1)
!
DGVMHVB = SUM(DGVTHVB(:))/DBLE(IGNFOLD-1)
DGVMHVT = SUM(DGVTHVT(:))/DBLE(IGNFOLD-1)
DGVMPPX = SUM(DGVTPPX(:))/DBLE(IGNFOLD-1)
DGVMPPZ = SUM(DGVTPPZ(:))/DBLE(IGNFOLD-1)
DGVMPSX = SUM(DGVTPSX(:))/DBLE(IGNFOLD-1)
DGVMPSZ = SUM(DGVTPSZ(:))/DBLE(IGNFOLD-1)
DGVMSHY = SUM(DGVTSHY(:))/DBLE(IGNFOLD-1)
DGVMSVX = SUM(DGVTSVX(:))/DBLE(IGNFOLD-1)
DGVMSVZ = SUM(DGVTSVZ(:))/DBLE(IGNFOLD-1)
!
WRITE(UNIT=20,FMT='(A,9(1X,E15.7E3))') 'MEAN:',DGMMHVB,DGMMHVT,DGMMPPX,DGMMPPZ,DGMMPSX,DGMMPSZ,DGMMSHY,DGMMSVX,DGMMSVZ
WRITE(UNIT=20,FMT='(A,9(1X,E15.7E3))') 'VARIANCE:',DGVMHVB,DGVMHVT,DGVMPPX,DGVMPPZ,DGVMPSX,DGVMPSZ,DGVMSHY,DGVMSVX,DGVMSVZ
!
DO IGINDPR = 1,IGNBPRM
   DGSIHVB(IGINDPR) = SUM(DGSUHVB(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSIHVT(IGINDPR) = SUM(DGSUHVT(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSIPPX(IGINDPR) = SUM(DGSUPPX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSIPPZ(IGINDPR) = SUM(DGSUPPZ(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSIPSX(IGINDPR) = SUM(DGSUPSX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSIPSZ(IGINDPR) = SUM(DGSUPSZ(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSISHY(IGINDPR) = SUM(DGSUSHY(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSISVX(IGINDPR) = SUM(DGSUSVX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSISVZ(IGINDPR) = SUM(DGSUSVZ(IGINDPR,:))/DBLE(IGNFOLD-1)
!
   DGSCHVB(IGINDPR) = SUM(DGSFHVB(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCHVT(IGINDPR) = SUM(DGSFHVT(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCPPX(IGINDPR) = SUM(DGSFPPX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCPPZ(IGINDPR) = SUM(DGSFPPZ(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCPSX(IGINDPR) = SUM(DGSFPSX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCPSZ(IGINDPR) = SUM(DGSFPSZ(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCSHY(IGINDPR) = SUM(DGSFSHY(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCSVX(IGINDPR) = SUM(DGSFSVX(IGINDPR,:))/DBLE(IGNFOLD-1)
   DGSCSVZ(IGINDPR) = SUM(DGSFSVZ(IGINDPR,:))/DBLE(IGNFOLD-1)
!
   WRITE(UNIT=20,FMT='(A,1X,I5,1X,I5,18(1X,F10.4))') TRIM(ADJUSTL(CGPRTYP(IGINDPR))),IGLAYPR(IGINDPR,1),IGLAYPR(IGINDPR,2)                 &
   ,DGSIHVB(IGINDPR),DGSIHVT(IGINDPR),DGSIPPX(IGINDPR),DGSIPPZ(IGINDPR),DGSIPSX(IGINDPR),DGSIPSZ(IGINDPR)                                  &
   ,DGSISHY(IGINDPR),DGSISVX(IGINDPR),DGSISVZ(IGINDPR)                                                                                     &
   ,DGSCHVB(IGINDPR),DGSCHVT(IGINDPR),DGSCPPX(IGINDPR),DGSCPPZ(IGINDPR),DGSCPSX(IGINDPR),DGSCPSZ(IGINDPR)                                  &
   ,DGSCSHY(IGINDPR),DGSCSVX(IGINDPR),DGSCSVZ(IGINDPR)
ENDDO
!
CLOSE(20)
!
WRITE(10,'(A)') "****************"
WRITE(10,'(A)') " END OF PROGRAM"
WRITE(10,'(A)') "****************"
!
CLOSE(10)
!
CALL CPU_TIME(TIME7)
!
WRITE(*,'(A,E15.7)') "CPU TIME FOR INITIALIZATION: ",TIME2-TIME1
WRITE(*,'(A,E15.7)') "CPU TIME FOR LOOP: ",TIME6-TIME2
WRITE(*,'(A,E15.7)') "CPU TIME FOR COMPUTATION OF MEAN AND TOTAL VARIANCE: ",TIME3
WRITE(*,'(A,E15.7)') "CPU TIME FOR COMPUTATION OF PARTIAL VARIANCES: ",TIME5
DO IGINDPR = 1,IGNBPRM
   WRITE(*,'(A,I10,A,E15.7)') "CPU TIME FOR PARAMETER ",IGINDPR,": ",TIME4(IGINDPR)
ENDDO
WRITE(*,'(A,E15.7)') "CPU TIME FOR COMPUTATION OF SOBOL COEFFICIENTS: ",TIME7-TIME6
WRITE(*,'(A,E15.7)') "TOTAL CPU TIME: ",TIME7-TIME1
!
!*******************************************************************************************************************************************
!
END PROGRAM CSOBOL
