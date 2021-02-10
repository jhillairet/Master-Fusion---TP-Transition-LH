!----------------------------------------
! file : sauvegarde.f90
! date : 09/03/2010
! used for savings
!----------------------------------------
subroutine saveDATA(itDIAG)
   use prec_const
   use globals
   use tools_module

   implicit none

   real(float), dimension(Nx) :: dxDens, dxPres, d2xPres, d3xPres
   real(float)                :: varn, varp, VEprime
   integer                    :: itDIAG

   call DERIV1(Dens,dxDens,Nx,dx,0)
   call DERIV1(Pres,dxPres,Nx,dx,0)
   call DERIV2(Pres,d2xPres,Nx,dx,0)
   call DERIV1(d2xPres,d3xPres,Nx,dx,0)
  
   do ix=1,Nx
     diagDens(itDIAG,ix)    = Dens(ix)
     diagPres(itDIAG,ix)    = Pres(ix)
     diagdxDens(itDIAG,ix)  = dxDens(ix)
     diagdxPres(itDIAG,ix)  = dxPres(ix)
     diagd2xPres(itDIAG,ix) = d2xPres(ix)
     diagd3xPres(itDIAG,ix) = d3xPres(ix)
   end do
   if (versionVEprime.eq.1) then
   !-------------------------------------
   ! VEprime = -[P''-n'P'/n]/n
   !-------------------------------------
     do ix=1,Nx
       VEprime = -(d2xPres(ix)-dxDens(ix)*dxPres(ix)/Dens(ix))/Dens(ix)
       varn    = ON+alphan*VEprime*VEprime
       varp    = ON+alphap*VEprime*VEprime
       diagDntot(itDIAG,ix)   = DDX0(ix) + DDX1(ix)/varn
       diagChiptot(itDIAG,ix) = CCX0(ix) + CCX1(ix)/varp
     end do
   elseif (versionVEprime.eq.2) then
   !-------------------------------------
   ! VEprime = (n'/n)^2		for density
   ! VEprime =  p'^2		for pressure
   !-------------------------------------
     do ix=1,Nx
       VEprime = dxDens(ix)*dxDens(ix)/(Dens(ix)*Dens(ix))
       varn    = ON+alphan*VEprime*VEprime
       VEprime = dxPres(ix)*dxPres(ix)
       varp    = ON+alphap*VEprime*VEprime
       diagDntot(itDIAG,ix)   = DDX0(ix) + DDX1(ix)/varn
       diagChiptot(itDIAG,ix) = CCX0(ix) + CCX1(ix)/varp
     end do
   end if

end subroutine saveDATA
