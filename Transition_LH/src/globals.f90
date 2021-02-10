!-------------------------------------------------------------------
! file : globals.f90
! date : 09/03/2010
! global variable definition
!-------------------------------------------------------------------
module globals
  use prec_const
  
  implicit none
  
! Geometry data - Time
  integer    , save  :: Nx     = 2**7
  integer    , save  :: Xscale = 2**0	       ! data saved = Nx / Xscale
  real(float), save  :: Lx     = 10._float
  real(float), save  :: dt     = 5.E-2_float
  real(float), save  :: dtDIAG = 5.E+1_float
  integer,     save  :: nbDIAG = 500
  real(float) :: T, dx
  integer     :: ITER, step
  
! Physical parameters
  integer    , save :: versionVEprime = 2
  integer    , save :: versionSource  = 3
  integer    , save :: hyperdiffusion = 0
  real(float), save :: Nu             = 1.E-3_float
  real(float), save :: alphan         = 5.E-2_float
  real(float), save :: alphap         = 5.E-2_float
  real(float), save :: D0             = 2.E-1_float
  real(float), save :: D1             = 1.E-0_float
  real(float), save :: Chi0           = 2.E-1_float
  real(float), save :: Chi1           = 1.E-0_float
  real(float), save :: coef           = 5.E+1_float
  real(float), save :: LD             = 2.E-1_float
  real(float), save :: Sn0            = 3.E-1_float
  real(float), save :: LSn            = 1.E-1_float
  real(float), save :: rhoSn          = 7.E-1_float
  real(float), save :: Sp0            = 2.5E-1_float
  real(float), save :: LSp            = 1.E-1_float
  real(float), save :: rhoSp          = 7.E-1_float

! Knobs
  logical    , save :: RST         = .false.	! ReSTart
  logical    , save :: singleField = .true.	! 1 dynamical field only

! Loops
  integer*8   :: ix, it

! Boundary Conditions
  real(float) :: DensBC_right, PresBC_right

! Input and output unit
  integer, parameter :: urst        = 50
  integer, parameter :: uprm        = 100
  integer, parameter :: uDensdat    = 200
  integer, parameter :: udxDensdat  = 210
  integer, parameter :: uPresdat    = 300
  integer, parameter :: udxPresdat  = 310
  integer, parameter :: ud2xPresdat = 320
  integer, parameter :: ud3xPresdat = 330
  integer, parameter :: uDntotdat   = 400
  integer, parameter :: uChiptotdat = 410

! Main Arrays
  real(float), dimension(:), pointer :: Dens, Pres
  real(float), dimension(:), pointer :: DDX0, DDX1
  real(float), dimension(:), pointer :: CCX0, CCX1
  real(float), dimension(:), pointer :: Sn, Sp

  real(float), dimension(:,:), pointer :: diagDens, diagdxDens
  real(float), dimension(:,:), pointer :: diagPres, diagdxPres, diagd2xPres
  real(float), dimension(:,:), pointer :: diagd3xPres
  real(float), dimension(:,:), pointer :: diagDntot, diagChiptot

end module globals
