!----------------------------------------
! file : evolution.f90
! date : 09/03/2010
! used for time evolution
!----------------------------------------
subroutine evolve
   use prec_const
   use globals
   use tools_module

   implicit none
   real(float), dimension(Nx) :: dxDens, dxPres, d2xPres, &
                                 coefDn_tot, coefCp_tot
   real(float) :: varn, varp, VEprime
 
!----------------------------------------
! Computation of the total transport coefficients:
! Chi_tot = Chi_0 + Chi1/(1+alpha*VEprime^2)
! with VEprime = -d[(dP/dx)/n]/dx
!----------------------------------------

   call DERIV1(Dens,dxDens,Nx,dx,0)
   call DERIV1(Pres,dxPres,Nx,dx,0)
   
   if (versionVEprime.eq.1) then
   !-------------------------------------
   ! VEprime = -[P''-n'P'/n]/n
   !-------------------------------------
     call DERIV2(Pres,d2xPres,Nx,dx,0)
     do ix=1,Nx
       VEprime        = -(d2xPres(ix)-dxDens(ix)*dxPres(ix)/Dens(ix))/Dens(ix)
       varn           = ON+alphan*VEprime*VEprime
       varp           = ON+alphap*VEprime*VEprime
       coefDn_tot(ix) = DDX0(ix) + DDX1(ix)/varn
       coefCp_tot(ix) = CCX0(ix) + CCX1(ix)/varp
     end do
   elseif (versionVEprime.eq.2) then
   !-------------------------------------
   ! VEprime = (n'/n)^2		for density
   ! VEprime =  p'^2		for pressure
   !-------------------------------------
     do ix=1,Nx
       VEprime        = dxDens(ix)*dxDens(ix)/(Dens(ix)*Dens(ix))
       varn           = ON+alphan*VEprime*VEprime
       VEprime        = dxPres(ix)*dxPres(ix)
       varp           = ON+alphap*VEprime*VEprime
       coefDn_tot(ix) = DDX0(ix) + DDX1(ix)/varn
       coefCp_tot(ix) = CCX0(ix) + CCX1(ix)/varp
     end do
   end if

!----------------------------------------
! Non linear diffusion operator
!----------------------------------------
   if (hyperdiffusion.eq.0) then 
     call DiffusionNL(coefDn_tot,coefCp_tot)
   else if (hyperdiffusion.eq.1) then 
     call DiffusionNL_hyperdiffusion(coefDn_tot,coefCp_tot)
   end if

end subroutine evolve
