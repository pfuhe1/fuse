SUBROUTINE MEAN_TIPOW()
! ---------------------------------------------------------------------------------------
! Creator:
! --------
! Martyn Clark, 2007
! ---------------------------------------------------------------------------------------
! Purpose:
! --------
! Computes the mean of the power-transformed topographic index
! ---------------------------------------------------------------------------------------
! Modules Modified:
! -----------------
! MODULE multiparam -- mean topographic index stored in MODULE multiparam
! ---------------------------------------------------------------------------------------
USE nrtype                                            ! variable types, etc.
USE nr, ONLY : gammp                                  ! interface for the incomplete gamma function
USE multiparam                                        ! model parameters
IMPLICIT NONE
! internal variables
INTEGER(I4B)                           :: IBIN        ! loop through bins
INTEGER(I4B), PARAMETER                :: NBINS=2000  ! number of bins in PDF of topo index
REAL(SP), PARAMETER                    :: TI_MAX=50._SP  ! maximum possible log-transformed index
REAL(SP)                               :: TI_OFF      ! offset in the Gamma distribution
REAL(SP)                               :: TI_SHP      ! shape of the Gamma distribution
REAL(SP)                               :: TI_CHI      ! CHI, see Sivapalan et al., 1987
REAL(SP)                               :: LOWERV      ! lower value of frequency bin
REAL(SP)                               :: UPPERV      ! upper value of frequency bin
REAL(SP)                               :: LOWERP      ! cumulative probability of the lower value
REAL(SP)                               :: UPPERP      ! cumulative probability of the upper value
REAL(SP)                               :: GMARG2      ! 2nd argument to the incomplete Gamma function
REAL(SP)                               :: PROBIN      ! probability of the current bin
REAL(SP)                               :: LOGVAL      ! log-transformed index for the current bin
REAL(SP)                               :: POWVAL      ! power-transformed index for the current bin
REAL(SP)                               :: AVEPOW      ! average power-transformed index
! ---------------------------------------------------------------------------------------
! preliminaries -- get parameters of the Gamma distribution (save typing)
TI_OFF = 3._SP              ! offset in the Gamma distribution (the "3rd" parameter)
TI_SHP = MPARAM%TISHAPE  ! shape of the Gamma distribution (the "2nd" parameter)
TI_CHI = (MPARAM%LOGLAMB - TI_OFF) / MPARAM%TISHAPE ! Chi -- loglamb is the first parameter (mean)
! loop through the frequency distribution
LOWERV = 0._SP
LOWERP = 0._SP
AVEPOW = 0._SP
DO IBIN=1,NBINS
 ! get probability for the current bin
 UPPERV = (REAL(IBIN)/REAL(NBINS)) * TI_MAX          ! upper value in frequency bin
 GMARG2 = MAX(0._SP, UPPERV - TI_OFF) / TI_CHI          ! 2nd argument to the Gamma function
 UPPERP = GAMMP(TI_SHP, GMARG2)                      ! GAMMP is the incomplete Gamma function
 PROBIN = UPPERP-LOWERP                              ! probability of the current bin
 ! get the scaled topographic index value
 LOGVAL = 0.5_SP*(LOWERV+UPPERV)                        ! log-transformed index for the current bin
 POWVAL = (EXP(LOGVAL))**(1._SP/MPARAM%QB_POWR)         ! power-transformed index for the current bin
 AVEPOW = AVEPOW + POWVAL*PROBIN                     ! average power-transformed index
 ! save the lower value and probability
 LOWERV = UPPERV                                     ! lower value for the next bin
 LOWERP = UPPERP                                     ! cumulative probability for the next bin
END DO  ! (looping through bins)
DPARAM%MAXPOW  = POWVAL
DPARAM%POWLAMB = AVEPOW
! ---------------------------------------------------------------------------------------
END SUBROUTINE MEAN_TIPOW
