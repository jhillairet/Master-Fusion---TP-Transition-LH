!----------------------------------------
! file : diffusionNL.f90
! date : 09/03/2010
! compute non linear diffusion terms
!----------------------------------------
subroutine DiffusionNL(Dn,Cp)
   use prec_const
   use globals
   use tools_module

   implicit none
   real(float), intent(in), dimension(Nx)   :: Dn, Cp
   real(float),             dimension(2,Nx) :: RR, AA, BB, CC
   real(float) :: alpha, coefN, coefP, &
     		      coefNa, coefNb, coefNc, &
     		      coefPa, coefPb, coefPc

   alpha = dt/(FO*dx**2)

!----------------------------------------
!---  Diffusion Matrix
   do ix=2,Nx-1
      coefNa = (Dn(ix-1)+Dn(ix))*alpha
      coefNc = (Dn(ix)+Dn(ix+1))*alpha
      coefNb = coefNa + coefNc
      coefPa = (Cp(ix-1)+Cp(ix))*alpha
      coefPc = (Cp(ix)+Cp(ix+1))*alpha
      coefPb = coefPa + coefPc

      RR(1,ix) = coefNa*Dens(ix-1)+(ON-coefNb)*Dens(ix)+coefNc*Dens(ix+1)+&
                 dt*Sn(ix)
      AA(1,ix) = -coefNa
      BB(1,ix) = ON+coefNb
      CC(1,ix) = -coefNc

      RR(2,ix) = coefPa*Pres(ix-1)+(ON-coefPb)*Pres(ix)+coefPc*Pres(ix+1)+&
                 dt*Sp(ix)
      AA(2,ix) = -coefPa
      BB(2,ix) = ON+coefPb
      CC(2,ix) = -coefPc
   end do

!----------------------------------------
!---  BOUNDARY CONDITIONS

! Density    
! Vanishing gradient at xmin
! (assuming vanishing Laplacian of diffus. coeff. at xmin)
   coefN    = FO*Dn(1)*alpha
   RR(1,1)  = (ON-coefN)*Dens(1)+coefN*Dens(2)+dt*Sn(1)
   BB(1,1)  = ON+coefN
   CC(1,1)  = -coefN

! Prescribed value "DensBCright" at xmax
!   RR(1,Nx) = DensBC_right
!   BB(1,Nx) = ON
!   AA(1,Nx) = ZE
   coefN    = FO*Dn(Nx)*alpha
   RR(1,Nx) = TW*alpha*(Dn(Nx-1)-Dn(Nx))*Dens(Nx-1) + (ON-coefN)*Dens(Nx) + &
              FO*alpha*(TH*Dn(Nx)-Dn(Nx-1))*DensBC_right
   BB(1,Nx) = ON+coefN
   AA(1,Nx) = -TW*alpha*(Dn(Nx-1)-Dn(Nx))

! Pressure   
! Vanishing gradient at xmin
   coefP    = FO*Cp(1)*alpha
   RR(2,1)  = (ON-coefP)*Pres(1)+coefP*Pres(2)+dt*Sp(1)
   BB(2,1)  = ON+coefP
   CC(2,1)  = -coefP

! Prescribed value "PresBCright" at xmax
!   RR(2,Nx) = PresBC_right
!   BB(2,Nx) = ON
!   AA(2,Nx) = ZE
   coefP    = FO*Cp(Nx)*alpha
   RR(2,Nx) = TW*alpha*(Cp(Nx-1)-Cp(Nx))*Pres(Nx-1) + (ON-coefP)*Pres(Nx) + &
              FO*alpha*(TH*Cp(Nx)-Cp(Nx-1))*PresBC_right
   BB(2,Nx) = ON+coefP
   AA(2,Nx) = -TW*alpha*(Cp(Nx-1)-Cp(Nx))

!----------------------------------------
!---  INVERSION TRIDIAGONAL MATRIX
   call TRIDIAG(AA,BB,CC,RR,Dens,Pres,Nx)

end subroutine DiffusionNL
!----------------------------------------


!----------------------------------------
! file : diffusionNL.f90
! date : 13/04/2010
! compute non linear diffusion terms and add of hyperdiffusion
!----------------------------------------
subroutine DiffusionNL_hyperdiffusion(Dn,Cp)
   use prec_const
   use globals
   use tools_module

   implicit none
   real(float), intent(in), dimension(Nx)   :: Dn,  Cp 
   real(FLOAT),             dimension(Nx)   :: RRn, RRp 
   real(float),             dimension(Nx,5) :: an,  ap
   real(float),             dimension(Nx,2) :: alP, alN
   
   real(float)  ::  alpha, coefNa, coefNb, coefNc,coefNd  ! coefNa=coefNe 
   real(float)  ::  coefPa, coefPb, coefPc,coefPd
   real(float)  ::  coefNj1, coefNjJ, coefPj1, coefPjJ,d,coefNj2,coefPj2
                     
   integer      ::  m1,m2, indxP(Nx), indxN(Nx)
		      
   ap    = 0.
   an    = 0.
   alP   = 0.
   alN   = 0.
   RRn   = 0.
   RRp   = 0.
   indxP = 0
   indxN = 0

   alpha = dt/(FO*dx**4)
   m1=2
   m2=2
!----------------------------------------
!---  Diffusion Matrix with term of hyperdiffusion

   do ix=3,Nx-2
      coefNa = TW * alpha * Nu 
      coefNb = - TW * alpha *((Dn(ix)+Dn(ix+1))*(dx**2)/TW + FO * Nu )
      coefNc =   ON + alpha * (( TW *Dn(ix) + Dn(ix+1) + Dn(ix-1))*dx**2 + &
           TW*TW*TH*Nu)
      coefNd = - TW * alpha *((Dn(ix)+Dn(ix-1))*(dx**2)/TW + FO * Nu )
      
      coefPa = TW * alpha * Nu 
      coefPb = - TW * alpha *((Cp(ix)+Cp(ix+1))*(dx**2)/TW + FO * Nu )
      coefPc =   ON + alpha * (( TW *Cp(ix) + Cp(ix+1) + Cp(ix-1))*dx**2 + &
           TW*TW*TH *Nu)
      coefPd = - TW * alpha *((Cp(ix)+Cp(ix-1))*(dx**2)/TW + FO * Nu )
      
      RRn(ix)= - coefNa*Dens(ix+2) - coefNb*Dens(ix+1) + &
           (TW-coefNc)*Dens(ix) - coefNd*Dens(ix-1) - &
           coefNa*Dens(ix-2) + dt*Sn(ix)
      an(ix,m1+1) = coefNc ! Diagonal
      an(ix, m1)  = coefNd ! subdiagonal
      an(ix,m1-1) = coefNa ! sub_subdiagonal
      an(ix,m1+2) = coefNb ! superdiagonal
      an(ix,m1+3) = coefNa ! super_superdiagonal
     
      RRp(ix) = - coefPb*Pres(ix+1) - coefPa*Pres(ix+2) + &
           (TW-coefPc)*Pres(ix) - coefPd*Pres(ix-1) - &
           coefPa*Pres(ix-2) + dt*Sp(ix)
      ap(ix,m1+1) = coefPc ! Diagonal
      ap(ix, m1)  = coefPd ! subdiagonal
      ap(ix,m1-1) = coefPa ! sub_subdiagonal
      ap(ix,m1+2) = coefPb ! superdiagonal
      ap(ix,m1+3) = coefPa ! super_superdiagonal
     
   end do

!----------------------------------------------------
!---  BOUNDARY CONDITIONS with term of hyperdiffusion

! Density
! Vanishing gradient at xmin
!-------initial-------------------------------------------
!(assuming vanishing gradient of diffus. coeff. at xmin=1)
     coefNj1    =  FO*alpha * Dn(1)*(dx**2)
   
     an(1,m1+1) = (ON + coefNj1 + FO*TH*alpha*Nu)!Diagonal
     an(1,m1+2) =  -coefNj1 - FO*FO*alpha*Nu    ! superdiagonal
     an(1,m1+3) = FO*alpha*Nu
     RRn(1)     = - FO*alpha*Nu*Dens(3) + &
          (coefNj1 + FO*FO*alpha*Nu ) * Dens(2) + &
          (ON - coefNj1 - FO*TH*alpha*Nu) * Dens(1)+ dt* Sn(1)

!(assuming vanishing gradient of diffus. coeff. at xmin=2)
     coefNa = TW * alpha * Nu 
     coefNb = - alpha *((Dn(2)+Dn(3))*(dx**2) + TW * FO * Nu )
     coefNc =   ON + &
          alpha * (( TW *Dn(2) + Dn(3) + Dn(1))*(dx**2) + (FO*TH+TW)*Nu)
     coefNd = - alpha *((Dn(2)+Dn(1))*(dx**2) + TW * FO * Nu )
   
     RRn(2)     = - coefNb*Dens(3) - coefNa*Dens(4) +&
          (TW-coefNc)*Dens(2) - coefNd*Dens(1) + &
          dt * Sn(2)
     an(2,m1+1) = coefNc ! Diagonal
     an(2, m1)  = coefNd ! subdiagonal
     an(2,m1+2) = coefNb ! superdiagonal
     an(2,m1+3) = coefNa ! super_superdiagonal

!--------final----------------------------------------------
!Prescribed value "Dens_BCright" at xmax=Nx

     coefNjJ    = FO*alpha * (Dn(Nx)*(dx**2)+TH*Nu)
   
     an(Nx,m1+1)= (ON + coefNjJ)                         ! Diagonal
     an(Nx, m1) = TW  * alpha * ((dx**2)*(Dn(Nx)-Dn(Nx-1)) - FO*Nu) ! subdiagonal 
     an(Nx,m1-1)= FO*alpha*Nu 
     RRn(Nx)    = (ON - coefNjJ)*Dens(Nx) - &
          TW*alpha*((dx**2)*(Dn(Nx)-Dn(Nx-1)) - FO*Nu) * Dens(Nx-1) - &
          FO*alpha*Nu*Dens(Nx-2) + &
          FO*alpha* ((dx**2)*(TH*Dn(Nx)-Dn(Nx-1)) + FO*Nu)* DensBC_right + &
          dt * Sn(Nx)
   
     ! Prescribed value "DensBCright" at xmax=Nx-1
     coefNa = TW * alpha * Nu 
     coefNb = - alpha *((Dn(Nx-1)+Dn(Nx))*(dx**2) + TW * FO * Nu )
     coefNc =   ON + &
          alpha * (( TW *Dn(Nx-1) + Dn(Nx) + Dn(Nx-2))*(dx**2) + (TW*FO+TW)*Nu)
     coefNd = - alpha *((Dn(Nx-1)+Dn(Nx-2))*(dx**2) + TW * FO * Nu )
      
     RRn(Nx-1)     = - coefNb*Dens(Nx)  + (TW-coefNc)*Dens(Nx-1) - &
          coefNd*Dens(Nx-2) - coefNa*Dens(Nx-3) + &
          dt*Sn(Nx-1) - TW*FO*alpha*Nu*DensBC_right

     an(Nx-1,m1+1) = coefNc ! Diagonal
     an(Nx-1, m1)  = coefNd ! subdiagonal
     an(Nx-1,m1-1) = coefNa ! sub_subdiagonal
     an(Nx-1,m1+2) = coefNb ! superdiagonal
     
! Pressure   
! Vanishing gradient at xmin
!-------initial-------------------------------------------
!(assuming vanishing gradient of diffus. coeff. at xmin=1)
     coefPj1    =  FO*alpha * Cp(1)*(dx**2)
   
     ap(1,m1+1) = (ON + coefPj1 + FO*TH*alpha*Nu)!Diagonal
     ap(1,m1+2) =  -coefPj1 - FO*FO*alpha*Nu    ! superdiagonal
     ap(1,m1+3) = FO*alpha*Nu
     RRp(1)     = - FO*alpha*Nu*Pres(3) + &
          (coefPj1 + FO*FO*alpha*Nu ) * Pres(2) + &
          (ON - coefPj1 - FO*TH*alpha*Nu) * Pres(1)+ dt* Sp(1)

!(assuming vanishing gradient of diffus. coeff. at xmin=2)
     coefPa = TW * alpha * Nu 
     coefPb = - alpha *((Cp(2)+Cp(3))*(dx**2) + TW * FO * Nu )
     coefPc =   ON + &
          alpha * (( TW *Cp(2) + Cp(3) + Cp(1))*(dx**2) + (FO*TH+TW)*Nu)
     coefPd = - alpha *((Cp(2)+Cp(1))*(dx**2) + TW * FO * Nu )
   
     RRp(2)     = - coefPb* Pres(3) - coefPa*Pres(4) +&
          (TW-coefPc)*Pres(2) - coefPd*Pres(1) + &
          dt * Sp(2)
     ap(2,m1+1) = coefPc ! Diagonal
     ap(2, m1)  = coefPd ! subdiagonal
     ap(2,m1+2) = coefPb ! superdiagonal
     ap(2,m1+3) = coefPa ! super_superdiagonal

!--------final----------------------------------------------
!Prescribed value "PresBCright" at xmax=Nx

    coefPjJ    = FO*alpha * (Cp(Nx)*(dx**2)+TH*Nu)
   
    ap(Nx,m1+1)= (ON + coefPjJ)                         ! Diagonal
    ap(Nx, m1) = TW  * alpha * ((dx**2)*(Cp(Nx)-Cp(Nx-1)) - FO*Nu) ! subdiagonal 
    ap(Nx,m1-1)= FO*alpha*Nu 
    RRp(Nx)    = (ON - coefPjJ)*Pres(Nx) - &
         TW*alpha*((dx**2)*(Cp(Nx)-Cp(Nx-1)) - FO*Nu) * Pres(Nx-1) - &
         FO*alpha*Nu*Pres(Nx-2) + &
         FO*alpha* ((dx**2)*(TH*Cp(Nx)-Cp(Nx-1)) + FO*Nu)* PresBC_right + &
         dt * Sp(Nx)
   
! Prescribed value "PresBCright" at xmax=Nx-1
    coefPa = TW * alpha * Nu 
    coefPb = - alpha *((Cp(Nx-1)+Cp(Nx))*(dx**2) + TW * FO * Nu )
    coefPc =   ON + &
         alpha * (( TW *Cp(Nx-1) + Cp(Nx) + Cp(Nx-2))*(dx**2) + (TW*FO+TW)*Nu)
    coefPd = - alpha *((Cp(Nx-1)+Cp(Nx-2))*(dx**2) + TW * FO * Nu )
      
    RRp(Nx-1)     = - coefPb*Pres(Nx)  + (TW-coefPc)*Pres(Nx-1) - &
         coefPd*Pres(Nx-2) - coefPa*Pres(Nx-3) + &
         dt*Sp(Nx-1) - TW*FO*alpha*Nu*PresBC_right

    ap(Nx-1,m1+1) = coefPc ! Diagonal
    ap(Nx-1, m1)  = coefPd ! subdiagonal
    ap(Nx-1,m1-1) = coefPa ! sub_subdiagonal
    ap(Nx-1,m1+2) = coefPb ! superdiagonal

!------------------------------------------------
!---  INVERSION  PentDIAGONAL MATRIX
  
 !Density
    call bandec(an,Nx,m1,m2,Nx,m1+m2+1,alN,m1,indxN,d) 
    call banbks(an,Nx,m1,m2,Nx,m1+m2+1,alN,m1,indxN,RRn)
    Dens = RRn
 !Pression
    call bandec(ap,Nx,m1,m2,Nx,m1+m2+1,alP,m1,indxP,d)
    call banbks(ap,Nx,m1,m2,Nx,m1+m2+1,alP,m1,indxP,RRp)
    Pres = RRp

end subroutine DiffusionNL_hyperdiffusion
!----------------------------------------
