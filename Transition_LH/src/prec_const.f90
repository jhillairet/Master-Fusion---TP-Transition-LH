!---------------------------------------------
! file : Prec_Const.f90
! date : 8/02/2002
! Constant and parameter initialisation
!---------------------------------------------

module prec_const
  implicit none
!*** Precision for real
  integer, parameter :: float = SELECTED_REAL_KIND(13,99)
  
!*** Some useful constants
  REAL(float),PARAMETER :: ZE     = 0.0_float
  REAL(float),PARAMETER :: HF     = 0.5_float
  REAL(float),PARAMETER :: ON     = 1.0_float
  REAL(float),PARAMETER :: TW     = 2.0_float
  REAL(float),PARAMETER :: TH     = 3.0_float
  REAL(float),PARAMETER :: FO     = 4.0_float
  REAL(float),PARAMETER :: PI     = 3.141592653589793238462643383279502884197_float
   
  COMPLEX(float),PARAMETER :: CZ  = (0._float, 0._float)
  COMPLEX(float),PARAMETER :: CO  = (1._float, 0._float)
  COMPLEX(float),PARAMETER :: IC  = (0._float, 1._float)

  COMPLEX(float),PARAMETER :: cst  = (1.E-0_float, 0._float)
end module prec_const

