!-----------------------------------------------------------------------------
! file : tools.f90
! date : 8/02/2002
! utility subroutines
!-----------------------------------------------------------------------------
module tools_module

implicit none
private
public :: TRIDIAG, BANDEC, BANBKS, DERIV1, DERIV2

interface TRIDIAG
   module procedure TRIDIAG
end interface

interface DERIV1
   module procedure DERIV1_r, DERIV1_c
end interface

interface DERIV2
   module procedure DERIV2_r, DERIV2_c
end interface

contains
!-----------------------------------------------------------------------------
! Solves for a vector U of length Nx the tridiagonal linear set given
! by the equation M*U = R , where A,B,C and R are input vectors and
! are not modified
subroutine TRIDIAG(A,B,C,R,U,V,Nx)
   use prec_const

   implicit none
   integer                           :: ix, Nx
   real(float), dimension(2,Nx+1) :: gam
   real(float), dimension(2,Nx)   :: A, B, C, R
   real(float), dimension(Nx)     :: U, V
   real(float)                    :: bet1, bet2

   if ((B(1,1).eq.0.).OR.(B(2,1).eq.0.)) print*,'error in tridiag'
   bet1 = B(1,1)
   bet2 = B(2,1)
   U(1) = R(1,1) / bet1
   V(1) = R(2,1) / bet2 

   do ix=2,Nx
      gam(1,ix) = C(1,ix-1)/bet1
      bet1  	= B(1,ix) - A(1,ix) * gam(1,ix)
      gam(2,ix) = C(2,ix-1)/bet2
      bet2  	= B(2,ix) - A(2,ix) * gam(2,ix)
      if ((bet1.eq.0.).OR.(bet2.eq.0.)) print*,'error in tridiag'
      U(ix)   = (R(1,ix) - A(1,ix) * U(ix-1)) / bet1
      V(ix)   = (R(2,ix) - A(2,ix) * V(ix-1)) / bet2		
   end do
   
   do ix=Nx-1,1,-1
      U(ix) = U(ix) - gam(1,ix+1) * U(ix+1)
      V(ix) = V(ix) - gam(2,ix+1) * V(ix+1)
   end do
end subroutine TRIDIAG

!-----------------------------------------------------------------------------------------
! Numerical Recipies routine to make LU decompostion of multiband matrix
subroutine bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
  use prec_const
  implicit none
  integer     :: m1,m2,mp,mpl,n,np,indx(n)
  real(float) :: d,a(np,mp) ,al(np,mpl) ,TINY
  parameter (TINY=1.e-20)
!Given an n x n band diagonal matrix A with m1 subdiagonal rows and m2 superdiagonal
!rows, compactly stored in the array a(1:n, 1:m1+m2-I-1) as described in the comment for
!routine banmul, this routine constructs an LU decomposition of a rcwwise permutation
!of A. The upper triangular matrix replaces a, while the lower triangular matrix is returned
!in al(1:n,1:m1). indx(1:n) is an output vector which records the row permutation
!effected by the partial pivoting; d is output as :l:1 depending on whether the number of
!row interchanges was even or odd, respectively. This routine is used in combination with
!banbks to solve band-diagonal sets of equations.
  integer     :: i,j,k,l,mm
  real(float) :: dum

  mm=m1+m2+1
  if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) print*,'bad args in bandec'
  l=m1
  do i=1,m1 ! Rearrange the storage a bit.
     do j=m1+2-i,mm
        a(i,j-l)=a(i,j)
     enddo
     l=l-1
     do j=mm-l,mm
        a(i,j)=0.
     end do
  enddo
  d=1.
  l=m1
  do k=1,n !For each row...
     dum=a(k,1)
     i=k
     if (l.lt.n) l=l+1
     do j=k+1,l ! Find the pivot element.
        if (abs(a(j,1)) .gt. abs(dum) )then
           dum=a(j,1)
           i=j
        endif
     enddo
     indx(k)=i
     if(dum.eq.0.) a(k,1)=TINY
!Matrix is algorithmically singular, but proceed anyway with TINY pivot (desirable in some
!applications).
     if (i.ne.k) then !Interchange rows.
        d=-d
        do j=1,mm
           dum=a(k,j)
           a(k,j)=a(i,j)
           a(i,j)=dum
        enddo
     endif
     do i=k+1,l ! Do the elimination.
        dum=a(i,1)/a(k,1)
        al(k,i-k)=dum
        do j=2,mm
           a(i,j-1)=a(i,j)-dum*a(k,j)
        enddo
        a(i,mm)=0.
     enddo
  enddo
  return
end subroutine bandec

!-----------------------------------------------------------------------------------------
! Numerical Recipies routine to solve linear multiband sytems
subroutine banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
  use prec_const
  implicit none
  integer     :: m1,m2,mp,mpl,n,np,indx(n)
  real(float) :: a(np,mp), al(np,mpl), b(n)
!-----------------------------------------------------------------------------------------
!Given the arrays a, al, and indx as returned from bandec, and given a right-hand side
!vector b(1:n), solves the band diagonal linear equations A.x = b. The solution vector x
!overwrites b(1:n). The other input arrays are not modified, and can be left in place for
!successive calls with different right-hand sides.
!-----------------------------------------------------------------------------------------
  integer     ::i,k,l,mm
  real(float) :: dum

  mm=m1+m2+1
  if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) print*,'bad args in banbks'
  l=m1
  do k=1,n     !Forward substitution, unscrambling the permuted rows as we go
     i=indx(k) 
     if (i .ne. k) then
        dum=b(k)
        b(k)=b(i)
        b(i)=dum
     endif
     if (l .lt. n) l=l+1
     do i=k+1,l
        b(i)=b(i)-al(k,i-k)*b(k)
     enddo
  enddo
  l=1
  do i=n,1,-1 ! Backsubstitution.
     dum=b(i)
     do k=2,l
        dum=dum-a(i,k)*b(k+i-1)
     enddo
     b(i)=dum/a(i,1)
     if (l .lt. mm) l=l+1
  enddo
  return
end subroutine banbks


!-----------------------------------------------------------------------------
! date : 30/10/2003
! derivee premiere calculee en differences finies a l'ordre 4
! 
!   input:	- U:  fonction a deriver (REELLE)
! 			- dx: pas du maillage
! 			- nx: nb de points de grille
! 			- BC: "Boundary Conditions": 0=non-periodique, 1=periodique
! 
!   output: - dxU: derivee premiere de U selon x (REELLE)
! 
subroutine DERIV1_r(U,dxU,nx,dx,BC)
  use prec_const
  
  implicit none 
  integer     :: nx, i, BC
  real(float) :: dx, cdx, c1a, c1b, c2a, c2b
  real(float), dimension(nx) :: U, dxU

  cdx = TW/(TH*dx)
  c1a = 25._float/12._float
  c1b = FO/TH
  c2a = 5._float/6._float
  c2b = TH/TW

  do i=3,nx-2
     dxU(i) = cdx*(U(i+1)-U(i-1) + (U(i-2)-U(i+2))/8._float)
  end do
  IF (BC.eq.0) THEN
     dxU(1) = (-c1a*U(1)+FO*U(2)-TH*U(3)+c1b*U(4)-U(5)/FO)/dx
     dxU(2) = (-U(1)/FO-c2a*U(2)+c2b*U(3)-HF*U(4)+U(5)/12._float)/dx
     dxU(nx)   = (c1a*U(nx)-FO*U(nx-1)+TH*U(nx-2)-c1b*U(nx-3)+U(nx-4)/FO)/dx
     dxU(nx-1) = (U(nx)/FO+c2a*U(nx-1)-c2b*U(nx-2)+HF*U(nx-3)-U(nx-4)/12._float)/dx
  ELSEIF (BC.eq.1) THEN
     dxU(1) = cdx*( U(2) - U(nx) + (U(nx-1)- U(3))/8._float )
     dxU(2) = cdx*( U(3) - U(1)  + (U(nx)  - U(4))/8._float )
     dxU(nx)   = cdx*( U(1)  - U(nx-1) + (U(nx-2) - U(2))/8._float )
     dxU(nx-1) = cdx*( U(nx) - U(nx-2) + (U(nx-3) - U(1))/8._float )
  ELSE
     print*,'Erreur dans le choix de BC'
     stop
  END IF

end subroutine DERIV1_r

!-----------------------------------------------------------------------------
! date : 30/10/2003
! derivee premiere calculee en differences finies a l'ordre 4
! 
!   input:	- U:  fonction a deriver (COMPLEXE)
! 			- dx: pas du maillage
! 			- nx: nb de points de grille
! 			- BC: "Boundary Conditions": 0=non-periodique, 1=periodique
! 
!   output: - dxU: derivee premiere de U selon x (COMPLEXE)
! 
subroutine DERIV1_c(U,dxU,nx,dx,BC)
  use prec_const
  
  implicit none 
  integer     :: nx, i, BC
  real(float) :: dx, cdx, c1a, c1b, c2a, c2b
  complex(float), dimension(nx) :: U, dxU

  cdx = TW/(TH*dx)
  c1a = 25._float/12._float
  c1b = FO/TH
  c2a = 5._float/6._float
  c2b = TH/TW

  do i=3,nx-2
     dxU(i) = cdx*(U(i+1)-U(i-1) + (U(i-2)-U(i+2))/8._float)
  end do
  IF (BC.eq.0) THEN
     dxU(1) = (-c1a*U(1)+FO*U(2)-TH*U(3)+c1b*U(4)-U(5)/FO)/dx
     dxU(2) = (-U(1)/FO-c2a*U(2)+c2b*U(3)-HF*U(4)+U(5)/12._float)/dx
     dxU(nx)   = (c1a*U(nx)-FO*U(nx-1)+TH*U(nx-2)-c1b*U(nx-3)+U(nx-4)/FO)/dx
     dxU(nx-1) = (U(nx)/FO+c2a*U(nx-1)-c2b*U(nx-2)+HF*U(nx-3)-U(nx-4)/12._float)/dx
  ELSEIF (BC.eq.1) THEN
     dxU(1) = cdx*( U(2) - U(nx) + (U(nx-1)- U(3))/8._float )
     dxU(2) = cdx*( U(3) - U(1)  + (U(nx)  - U(4))/8._float )
     dxU(nx)   = cdx*( U(1)  - U(nx-1) + (U(nx-2) - U(2))/8._float )
     dxU(nx-1) = cdx*( U(nx) - U(nx-2) + (U(nx-3) - U(1))/8._float )
  ELSE
     print*,'Erreur dans le choix de BC'
     stop
  END IF

end subroutine DERIV1_c

!-----------------------------------------------------------------------------
! date : 3/03/2004
! derivee seconde calculee en differences finies a l'ordre 4
! 
!   input:	- U:  fonction a deriver (REELLE)
! 			- dx: pas du maillage
! 			- nx: nb de points de grille
! 			- BC: "Boundary Conditions": 0=non-periodique, 1=periodique
! 
!   output: - dxxU: derivee seconde de U selon x (REELLE)
! 
subroutine DERIV2_r(U,dxxU,nx,dx,BC)
  use prec_const
  
  implicit none 
  integer     :: nx, i, BC
  real(float) :: dx, dx2, cdx2, c1a, c1b, c1c, c1d, c1e, c2b
  real(float), dimension(nx) :: U, dxxU

  dx2  = dx*dx
  cdx2 = ON/(12._float*dx2)
  c1a   = 35._float/12._float
  c1b   = 26._float/TH
  c1c   = 19._float/TW
  c1d   = 14._float/TH
  c1e   = 11._float/12._float
  c2b   = 5._float/TH

  do i=3,nx-2
     dxxU(i) = cdx2*(-30._float*U(i)+16._float*(U(i+1)+U(i-1))-U(i+2)-U(i-2))
  end do
  IF (BC.eq.0) THEN
     dxxU(1) = (c1a*U(1)-c1b*U(2)+c1c*U(3)-c1d*U(4)+c1e*U(5))/dx2
     dxxU(2) = (c1e*U(1)-c2b*U(2)+HF*U(3)+U(4)/TH-U(5)/12._float)/dx2
     dxxU(nx)   = (c1a*U(nx)-c1b*U(nx-1)+c1c*U(nx-2)-c1d*U(nx-3)+c1e*U(nx-4))/dx2
     dxxU(nx-1) = (c1e*U(nx)-c2b*U(nx-1)+HF*U(nx-2)+U(nx-3)/TH-U(nx-4)/12._float)/dx2
  ELSEIF (BC.eq.1) THEN
     dxxU(1) = cdx2*(-30._float*U(1)+16._float*(U(2)+U(nx))-U(3)-U(nx-1))
     dxxU(2) = cdx2*(-30._float*U(2)+16._float*(U(3)+U(1))-U(4)-U(nx))
     dxxU(nx)   = cdx2*(-30._float*U(nx)+16._float*(U(1)+U(nx-1))-U(2)-U(nx-2))
     dxxU(nx-1) = cdx2*(-30._float*U(nx-1)+16._float*(U(nx)+U(nx-2))-U(1)-U(nx-3))
  ELSE
     print*,'Erreur dans le choix de BC'
     stop
  END IF

end subroutine DERIV2_r

!-----------------------------------------------------------------------------
! date : 3/03/2004
! derivee seconde calculee en differences finies a l'ordre 4
! 
!   input:	- U:  fonction a deriver (COMPLEXE)
! 			- dx: pas du maillage
! 			- nx: nb de points de grille
! 			- BC: "Boundary Conditions": 0=non-periodique, 1=periodique
! 
!   output: - dxxU: derivee seconde de U selon x (COMPLEXE)
! 
subroutine DERIV2_c(U,dxxU,nx,dx,BC)
  use prec_const
  
  implicit none 
  integer     :: nx, i, BC
  real(float) :: dx, dx2, cdx2, c1a, c1b, c1c, c1d, c1e, c2b
  complex(float), dimension(nx) :: U, dxxU

  dx2  = dx*dx
  cdx2 = ON/(12._float*dx2)
  c1a   = 35._float/12._float
  c1b   = 26._float/TH
  c1c   = 19._float/TW
  c1d   = 14._float/TH
  c1e   = 11._float/12._float
  c2b   = 5._float/TH

  do i=3,nx-2
     dxxU(i) = cdx2*(-30._float*U(i)+16._float*(U(i+1)+U(i-1))-U(i+2)-U(i-2))
  end do
  IF (BC.eq.0) THEN
     dxxU(1) = (c1a*U(1)-c1b*U(2)+c1c*U(3)-c1d*U(4)+c1e*U(5))/dx2
     dxxU(2) = (c1e*U(1)-c2b*U(2)+HF*U(3)+U(4)/TH-U(5)/12._float)/dx2
     dxxU(nx)   = (c1a*U(nx)-c1b*U(nx-1)+c1c*U(nx-2)-c1d*U(nx-3)+c1e*U(nx-4))/dx2
     dxxU(nx-1) = (c1e*U(nx)-c2b*U(nx-1)+HF*U(nx-2)+U(nx-3)/TH-U(nx-4)/12._float)/dx2
  ELSEIF (BC.eq.1) THEN
     dxxU(1) = cdx2*(-30._float*U(1)+16._float*(U(2)+U(nx))-U(3)-U(nx-1))
     dxxU(2) = cdx2*(-30._float*U(2)+16._float*(U(3)+U(1))-U(4)-U(nx))
     dxxU(nx)   = cdx2*(-30._float*U(nx)+16._float*(U(1)+U(nx-1))-U(2)-U(nx-2))
     dxxU(nx-1) = cdx2*(-30._float*U(nx-1)+16._float*(U(nx)+U(nx-2))-U(1)-U(nx-3))
  ELSE
     print*,'Erreur dans le choix de BC'
     stop
  END IF

end subroutine DERIV2_c
!-----------------------------------------------------------------------------

end module tools_module
