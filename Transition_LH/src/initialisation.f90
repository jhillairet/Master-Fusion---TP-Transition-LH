!----------------------------------------
! file : initialisation.f90
! date : 09/03/2010
! used for initialisation
!----------------------------------------
!-----------------------------------------------------------------
!  Initialize Boundary Conditions
!-----------------------------------------------------------------
subroutine init_BC
  use prec_const
  use globals

  DensBC_right = ON
  PresBC_right = ON

end subroutine init_BC


!-----------------------------------------------------------------
!  Initialize Density and Pressure Profiles
!-----------------------------------------------------------------
subroutine init_prof
  use prec_const
  use globals

  real(float), dimension(Nx) :: dxDens, dxPres
  real(float) :: X, rho

  if (D1.ne.ZE) then
    do ix = 1,Nx
!      dxDens(ix) = - sum(Sn(1:ix))/DDX1(ix) * dx
!      Dens(ix)   = sum(dxDens(1:ix))*dx
      Dens(ix)   = DensBC_right
    end do
    Dens = Dens - Dens(Nx) + DensBC_right
  else
    do ix = 1,Nx
      rho      = dfloat(ix-1)/dfloat(Nx-1) 
      Dens(ix) = DensBC_right+ON-rho**2
    end do
  end if
  if (Chi1.ne.ZE) then
    do ix = 1,Nx
!      dxPres(ix) = - sum(Sp(1:ix))/CCX1(ix) * dx
!      Pres(ix)   = sum(dxPres(1:ix))*dx
      Pres(ix)   = PresBC_right
    end do
    Pres = Pres - Pres(Nx) + PresBC_right
  else
    do ix = 1,Nx
      rho      = dfloat(ix-1)/dfloat(Nx-1) 
      Pres(ix) = PresBC_right+ON-rho**2
    end do
  end if
  
end subroutine init_prof

!-----------------------------------------------------------------
!  Initialize Diffusion coefficients
!-----------------------------------------------------------------
subroutine init_diffus
  use prec_const
  use globals

  implicit none
  !real(float) :: X, varTANH, coef, LD
  real(float) :: X, varTANH

  !coef = 50._float
  !LD   = Lx/30._float

  do ix = 1,Nx
     X  = dfloat(ix-1)*dx
     ! varTANH   = TW - tanh(X/LD) - tanh((Lx-X)/LD)
     varTANH   = ON - tanh((Lx-X)/LD)
     DDX0(ix)  = D0*(ON+coef*varTANH)
     DDX1(ix)  = D1
     CCX0(ix)  = Chi0*(ON+coef*varTANH)
     CCX1(ix)  = Chi1
  end do
  
end subroutine init_diffus

!-----------------------------------------------------------------
!  Initialize Source terms
!-----------------------------------------------------------------
subroutine init_Sterms
  use prec_const
  use globals

  implicit none
  real(float) :: rho, rho1, rho2, rho3, X, Xmax, frac
  
  Xmax = dfloat(Nx-1)*dx
  frac = 1/6._float
  
  if(versionSource .eq. 1) then 
    do ix = 1,Nx
       X      = dfloat(ix-1)*dx
       rho    = X/Xmax
       Sn(ix) = Sn0 * exp(-((rho-rhoSn)/LSn)**2)
       Sp(ix) = Sp0 * exp(-((rho-rhoSp)/LSp)**2)
    end do
  else if  (versionSource .eq. 2) then
    do ix = 1,Nx
      X      = dfloat(ix-1)*dx
      rho    = X/Xmax
      Sn(ix) = Sn0 *tanh(-((rho-rhoSn)/LSn))
      Sp(ix) = Sp0 *tanh(-((rho-rhoSp)/LSp))
    end do
  else if  (versionSource .eq. 3) then
    do ix = 1,int(frac*Nx)
      Sn(ix) = ZE
      Sp(ix) = ZE
    end do
    rho1  = dfloat(int(frac*Nx)-1)*dx/Xmax
    do ix = int(frac*Nx),int(TW*frac*Nx)
      X      = dfloat(ix-1)*dx
      rho    = X/Xmax
      Sn(ix) = Sn0*(rho-rho1)/LSn  
      Sp(ix) = Sp0*(rho-rho1)/LSp
    end do
    rho2  = dfloat(int(TW*frac*Nx)-1)*dx/Xmax
    do ix = int(TW*frac*Nx),int(FO*frac*Nx)
      Sn(ix) = Sn0*(rho2-rho1)/LSn
      Sp(ix) = Sp0*(rho2-rho1)/LSp
    end do
    rho3  = dfloat(int(FO*frac*Nx)-1)*dx/Xmax
    do ix = int(FO*frac*Nx),int((FO+ON)*frac*Nx)
      X      = dfloat(ix-1)*dx
      rho    = X/Xmax
      Sn(ix) = Sn0*((rho2-rho1)/LSn-(rho-rho3)/LSn)
      Sp(ix) = Sp0*((rho2-rho1)/LSp-(rho-rho3)/LSp)
    end do
    do ix = int((FO+ON)*frac*Nx),Nx
      Sn(ix) = ZE
      Sp(ix) = ZE
    end do
  end if

end subroutine init_Sterms

!-----------------------------------------------------------------
! Restart
!-----------------------------------------------------------------
subroutine restart
  use prec_const
  use globals

  open(urst, file = 'transitionLH.rst', status = 'OLD', &
       form = 'UNFORMATTED')
     read(urst) Dens
     read(urst) Pres
     read(urst) DensBC_right
     read(urst) PresBC_right
  close(urst)
end subroutine restart

!-----------------------------------------------------------------
! Main initialisation
!-----------------------------------------------------------------
subroutine initialise
  use prec_const
  use globals, only : RST

  call init_BC
  call init_diffus
  call init_Sterms
  IF (RST) THEN
     print*,' RESTART '
     call restart
  ELSE
     call init_prof
  ENDIF
end subroutine initialise

