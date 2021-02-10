!-----------------------------------------------------------------------
! Code derived from Hysteresis_2Fields/*.f90
!
!
! Resolution of the following system:
!
!    d/dt N - d/dx{[D0(x) + Dturb] d/dx N}     = S(x)
!    d/dt P - d/dx{[Chi0(x) + Chiturb] d/dx P} = H(x)
!    Dturb   = D1(x)   / (1 + alpha VEprime^2)
!    Chiturb = Chi1(x) / (1 + alpha VEprime^2)
!    VEprime = d/dx (dP/dx / N)
!
!   where:  N  ->  Density
!	    P  ->  Pressure
!
!   Numerical scheme:  Trapezoidal leap-frog in time
!                      Crank-Nicholson for diffusive terms
!-----------------------------------------------------------------------

PROGRAM LH_transition

use prec_const
use globals

implicit none

!-----------------------------------------------------------------------------
!	Definition of variables
!-----------------------------------------------------------------------------

integer  itDIAG

itDIAG = 1

call read_input

print*,' Nx     = ',Nx
print*,' nbDIAG = ',nbDIAG
print*,' '

!===============================================================
!call read_input
call initialise
call saveDATA(itDIAG)

DO it=1,ITER
   call evolve
   IF (mod(it,step) .EQ. 0) THEN
      itDIAG = itDIAG+1
      print*,' itDIAG = ',itDIAG
      call saveDATA(itDIAG)
   ENDIF
END DO
      
!  Create restart file
open(urst, file = 'transitionLH.rst', form = 'UNFORMATTED', &
	 status='UNKNOWN')
   write(urst) Dens
   write(urst) Pres
   write(urst) DensBC_right
   write(urst) PresBC_right
close(urst)
! ==============================================================
!---     STORAGE

open(uDensdat,file='uDens.dat',status='unknown')
open(uPresdat,file='uPres.dat',status='unknown')
open(udxDensdat,file='udxDens.dat',status='unknown')
open(udxPresdat,file='udxPres.dat',status='unknown')
open(ud2xPresdat,file='ud2xPres.dat',status='unknown')
open(ud3xPresdat,file='ud3xPres.dat',status='unknown')
open(uDntotdat,file='uDntot.dat',status='unknown')
open(uChiptotdat,file='uChiptot.dat',status='unknown')

do it=1,nbDIAG
   write(uDensdat,2000)    (diagDens(it,ix),   ix=1,Nx,Xscale)
   write(uPresdat,2000)    (diagPres(it,ix),   ix=1,Nx,Xscale)
   write(udxDensdat,2000)  (diagdxDens(it,ix), ix=1,Nx,Xscale)
   write(udxPresdat,2000)  (diagdxPres(it,ix), ix=1,Nx,Xscale)
   write(ud2xPresdat,2000) (diagd2xPres(it,ix),ix=1,Nx,Xscale)
   write(ud3xPresdat,2000) (diagd3xPres(it,ix),ix=1,Nx,Xscale)
   write(uDntotdat,2000)   (diagDntot(it,ix),  ix=1,Nx,Xscale)
   write(uChiptotdat,2000) (diagChiptot(it,ix),ix=1,Nx,Xscale)
end do

close(uDensdat)
close(uPresdat)
close(udxDensdat)
close(udxPresdat)
close(ud2xPresdat)
close(ud3xPresdat)
close(uDntotdat)
close(uChiptotdat)

2000  format(1pE15.6)
2001  format(1pE15.6, 1pE15.6)

stop
end

! ==============================================================
