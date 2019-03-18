SUBROUTINE UPDATE_SWE(DT)
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Brian Henn, as part of FUSE snow model implementation, 6/2013
! Based on subroutines QSATEXCESS and UPDATSTATE, by Martyn Clark
! Modified by Nans Addor to enable distributed modeling, 9/2016
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the snow accumulation and melt from forcing data
! Then updates the SWE band states based on the fluxes
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multibands -- SWE bands stored in MODULE multibands
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE model_defn                                        ! model definition structure
USE model_defnames                                    ! integer model definitions
USE multiparam                                        ! model parameters
USE multibands                                        ! model basin band structure
USE multiforce                                        ! model forcing
USE multistate                                        ! model states
USE multi_flux                                        ! model fluxes
IMPLICIT NONE
! input
REAL(SP), INTENT(IN)                   :: DT             ! length of the time step
! internal variables
LOGICAL(LGT)                           :: LEAP           ! leap year flag
REAL(SP)                               :: JDAY           ! Julian day of year
REAL(SP), DIMENSION(12)                :: CUMD           ! cumulative days of year
REAL(SP)                               :: DZ             ! vert. distance from forcing
REAL(SP)                               :: TEMP_Z         ! band temperature at timestep
REAL(SP)                               :: PRECIP_Z       ! band precipitation at timestep
REAL(SP)                               :: MF             ! melt factor (mm/deg.C-6hr)
INTEGER(I4B)                           :: ISNW           ! loop through snow model bands
! Snow redistribution variables (TODO some of these should be parameters?)
REAL(SP)                               :: Z_REDIST_UP = 4000._sp   ! snow is redistributed from levels above this hight
REAL(SP)                               :: Z_REDIST_LOW = 1900._sp  ! snow is not redistributed below this height
REAL(SP)                               :: SWE_MAX = 500._sp
REAL(SP)                               :: SWE_REDIST_BULK ! SWE redistributed (bulk over whole cell) 
REAL(SP)                               :: SWE_REDIST      ! SWE redistributed (over fraction)
REAL(SP)                               :: SWE_TOT0  ! Diagnostic: total swe before redistribution
REAL(SP)                               :: SWE_TOT1  ! Diagnostic: total swe after redistribution
INTEGER(I4B),DIMENSION(:),ALLOCATABLE  :: REDIST_TARGET  ! 'logical' value used for determining bands to redistribute from/to
REAL(SP),DIMENSION(:),ALLOCATABLE      :: REDIST_AMOUNTS ! SWE removed from each band in redistribution
INTEGER(I4B)                           :: IERR        ! error code

! allocate data stuctures (TODO) maybe allocate in get_mbands to speed up?
ALLOCATE(REDIST_TARGET(N_BANDS),STAT=IERR)
ALLOCATE(REDIST_AMOUNTS(N_BANDS),STAT=IERR)
! Initialise redist amounts
REDIST_AMOUNTS = 0._sp

! ---------------------------------------------------------------------------------------
! snow accumulation and melt calculations for each band
! also calculates effective precipitation
! ---------------------------------------------------------------------------------------
! first calculate day of year for melt factor calculation
CUMD = real((/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /),sp)
IF (MOD(timDat%IY,4).EQ.0) THEN
 LEAP = .TRUE.
 CUMD(3:12) = CUMD(3:12) + 1._sp
ELSE
 LEAP = .FALSE.
ENDIF

JDAY = CUMD(timDat%IM) + timDat%ID

IF (LEAP) THEN   ! calculate melt factor from MFMAX, MFMIN and day of year
 MF = ((0.5_sp*SIN(((JDAY-81._sp)*2._sp*PI)/366._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ELSE
 MF = ((0.5_sp*SIN(((JDAY-80._sp)*2._sp*PI)/365._sp))+0.5_sp)*(MPARAM%MFMAX - MPARAM%MFMIN) + MPARAM%MFMIN
ENDIF

! loop through model bands
DO ISNW=1,N_BANDS

 ! calculate forcing data for each band
 DZ = MBANDS(ISNW)%Z_MID - Z_FORCING
 TEMP_Z = MFORCE%TEMP + DZ*MPARAM%LAPSE/1000._sp ! adjust for elevation using lapse rate
 IF (DZ.GT.0._sp) THEN ! adjust for elevation using OPG
  PRECIP_Z = MFORCE%PPT * (1._sp + DZ*MPARAM%OPG/1000._sp)
 ELSE
  PRECIP_Z = MFORCE%PPT / (1._sp - DZ*MPARAM%OPG/1000._sp)
 ENDIF
 IF ((MBANDS(ISNW)%SWE.GT.0._sp).AND.(TEMP_Z.GT.MPARAM%MBASE)) THEN
  ! calculate the initial snowmelt rate from the melt factor and the temperature
  MBANDS(ISNW)%SNOWMELT = MF*(TEMP_Z - MPARAM%MBASE) ! MBANDS%SNOWMELT has units of mm day-1
 ELSE
  MBANDS(ISNW)%SNOWMELT = 0.0_sp
 ENDIF

 ! calculate the accumulation rate from the forcing data
 IF (TEMP_Z.LT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
    MBANDS(ISNW)%SNOWACCMLTN = MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
    MBANDS(ISNW)%SNOWACCMLTN = PRECIP_Z * MPARAM%RFERR_MLT
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  MBANDS(ISNW)%SNOWACCMLTN = 0.0_sp
 ENDIF

 ! update SWE, and check to ensure non-negative values
 MBANDS(ISNW)%DSWE_DT = MBANDS(ISNW)%SNOWACCMLTN - MBANDS(ISNW)%SNOWMELT
 IF ((MBANDS(ISNW)%SWE + MBANDS(ISNW)%DSWE_DT*DT).GE.0._sp) THEN
  MBANDS(ISNW)%SWE = MBANDS(ISNW)%SWE + MBANDS(ISNW)%DSWE_DT*DT
 ELSE ! reduce melt rate in case of negative SWE
  MBANDS(ISNW)%SNOWMELT = MBANDS(ISNW)%SWE/DT + MBANDS(ISNW)%SNOWACCMLTN
  MBANDS(ISNW)%SWE = 0.0_sp
 ENDIF

 ! calculate rainfall plus snowmelt
 IF (TEMP_Z.GT.MPARAM%PXTEMP) THEN
  SELECT CASE(SMODL%iRFERR)
   CASE(iopt_additive_e) ! additive rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * &
   (MAX(0.0_sp, PRECIP_Z + MPARAM%RFERR_ADD) + MBANDS(ISNW)%SNOWMELT)
   CASE(iopt_multiplc_e) ! multiplicative rainfall error
   M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * &
   (PRECIP_Z * MPARAM%RFERR_MLT +  MBANDS(ISNW)%SNOWMELT)
   CASE DEFAULT       ! check for errors
    print *, "SMODL%iRFERR must be either iopt_additive_e or iopt_multiplc_e"
    STOP
  END SELECT
 ELSE
  M_FLUX%EFF_PPT = M_FLUX%EFF_PPT + MBANDS(ISNW)%AF * MBANDS(ISNW)%SNOWMELT
 ENDIF
END DO

! Calculate amount of snow to redistribute and remove from each band
REDIST_TARGET = (MBANDS%Z_MID .GT. Z_REDIST_UP .AND. MBANDS%SWE .GT. SWE_MAX) *-1
!print*,'DEBUG, redist_target: removal',REDIST_TARGET
IF (SUM(REDIST_TARGET).GT.0) THEN
 SWE_TOT0 = SUM(MBANDS%SWE * MBANDS%AF)
ENDIF
WHERE(REDIST_TARGET.EQ.1)
 REDIST_AMOUNTS = (MBANDS%SWE - SWE_MAX) * MBANDS%AF
 MBANDS%SWE = SWE_MAX
END WHERE
! SWE mass to redistribute
SWE_REDIST_BULK = SUM(REDIST_AMOUNTS)

! Transfer redistributed snow to lower levels
REDIST_TARGET = (MBANDS%Z_MID .GT. Z_REDIST_LOW .AND. MBANDS%Z_MID .LT. Z_REDIST_UP) *-1
!print*,'DEBUG, redist_target: deposit',REDIST_TARGET
IF (SWE_REDIST > 0._sp) THEN
 IF (SUM(REDIST_TARGET).EQ.0) THEN
  ! All levels are above Z_REDIST_LOW: put snow into bottom level
  print *,'WARNING, All bands are above Z_REDIST_LOW, distributing snow to bottom band',MBANDS(1)%Z_MID,SWE_REDIST
  REDIST_TARGET(1)=1
 END IF
 ! For each band increase SWE by same height, reduce by fraction of cell
 SWE_REDIST = SWE_REDIST_BULK / SUM(REDIST_TARGET*MBANDS%AF) 
 ! Add 
 WHERE(REDIST_TARGET.EQ.1) MBANDS%SWE = MBANDS%SWE + SWE_REDIST
 !print *,'DEBUG: tot_swe bfore/after redist',SWE_TOT0,SUM(MBANDS%SWE * MBANDS%AF),SWE_REDIST
 SWE_TOT1 = SUM(MBANDS%SWE * MBANDS%AF)
 IF( SWE_TOT1-SWE_TOT0 .GT. 1E-6) print*,'WARNING snow redistribution not conserving: swe before/after,distributed',SWE_TOT0,SWE_TOT1,SWE_REDIST_BULK

END IF ! SWE_REDIST

! DEALLOCATE vars
DEALLOCATE(REDIST_TARGET, STAT=IERR)
IF (IERR.NE.0) THEN
 print *,'ERROR! update_swe/problem deallocating REDIST_TARGET'
END IF
DEALLOCATE(REDIST_AMOUNTS, STAT=IERR)
IF (IERR.NE.0) THEN
 print *,'ERROR! update_swe/problem deallocating REDIST_AMOUNTS'
END IF

END SUBROUTINE UPDATE_SWE
