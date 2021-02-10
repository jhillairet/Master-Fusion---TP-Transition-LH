!----------------------------------------
! file : read_data.f90
! date : 09/03/2010
! used for the reading of :
! - mesh data
! - the knobs data
!----------------------------------------

!------------------------------------
! mesh data reading
!------------------------------------
subroutine mesh_data
  use prec_const
  use globals
  implicit none
  
  namelist /MESH/ Nx, Xscale, Lx, dt, dtDIAG, nbDIAG

  read(*,MESH)
  write(*,MESH)
  
  dx     = Lx/dfloat(Nx-1)
  T      = dfloat((nbDIAG-1))*dtDIAG
  ITER   = ifix(T/dt)
  step   = ifix(dtDIAG/dt)

  print*,'    '
  print*,'Tmax      = ',T
  print*,'dt        = ',dt
  print*,'nb iterations    = ',ITER
  print*,'nb diagnostics   = ',nbDIAG
  print*,'    '
  
end subroutine mesh_data

!----------------------------------------
! Physical parameters reading
!----------------------------------------
subroutine physparam_data
  use prec_const
  use globals
  implicit none
  real(float) :: CFL
  
  namelist /PhysParam/ versionVEprime, versionSource, hyperdiffusion, &
                       Nu, alphan, alphap, &
                       D0, D1, Chi0, Chi1, coef, LD, &
                       Sn0, LSn, rhoSn, Sp0, LSp, rhoSp
  
  read(*,PhysParam)
  write(*,PhysParam)

  !--------------
  ! CFL Condition
  !--------------
  CFL = dt*Nu/(dx**4)
  if (CFL .gt. 1) then
     print*,'------------------------------------'
     print*,'WARNING, CFL CRITERION NOT SATISFIED'
     print*,'CFL VALUE :',CFL
     print*,'------------------------------------'
  end if
end subroutine physparam_data

!----------------------------------------
! knobs data reading
!----------------------------------------
subroutine knobs_data
  use prec_const
  use globals
  implicit none
  
  namelist /KNOBS/ RST, singleField
  
  read(*,KNOBS)
  write(*,KNOBS)

  print*,'    '
  if (RST) then
     print*,'This is a RESTART simulation'
     print*,'    '
  endif
end subroutine knobs_data

!-----------------------------------------------------------------
!  Allocate dimensions of arrays
!-----------------------------------------------------------------
subroutine create_array
  use prec_const
  use globals

  allocate(Dens(Nx))
  allocate(Pres(Nx))

  allocate(DDX0(Nx))
  allocate(DDX1(Nx))
  allocate(CCX0(Nx))
  allocate(CCX1(Nx))

  allocate(Sn(Nx))
  allocate(Sp(Nx))

  allocate(diagDens(nbDIAG,Nx))
  allocate(diagPres(nbDIAG,Nx))
  allocate(diagdxDens(nbDIAG,Nx))
  allocate(diagdxPres(nbDIAG,Nx))
  allocate(diagd2xPres(nbDIAG,Nx))
  allocate(diagd3xPres(nbDIAG,Nx))
  allocate(diagDntot(nbDIAG,Nx))
  allocate(diagChiptot(nbDIAG,Nx))

end subroutine create_array

!-----------------------------------------------------------------
!  Save parameters
!-----------------------------------------------------------------
subroutine save_param
  use prec_const
  use globals
  
  implicit none
  integer :: iRST, isingleField
    
  iRST=0
  if (RST) then
    iRST=1
  end if

  isingleField=0
  if (singleField) then
    isingleField=1
  end if

  open(uprm,file='uprm.dat',status='unknown')
     write(uprm,5000) Nx, Xscale, nbDIAG, versionVEprime, versionSource, hyperdiffusion
     write(uprm,5001) dx, Lx, dt, dtDIAG, coef, LD
     write(uprm,5001) D0, D1, Chi0, Chi1, Nu, ZE
     write(uprm,5001) Sn0, LSn, rhoSn, alphan, ZE, ZE
     write(uprm,5001) Sp0, LSp, rhoSp, alphap, ZE, ZE
     write(uprm,5000) iRST, isingleField, isingleField, isingleField, 0, 0
  close(uprm)

5000  format(I6,' ',I6,' ',I6,' ',I6,' ',I6,' ',I6)
5001  format(G10.4e2,' ',G10.4e2,' ',G10.4e2,' ',G10.4e2,' ',G10.4e2,' ',G10.4e2)

end subroutine save_param

!----------------------------------------
!  reading of the input file
!----------------------------------------
subroutine read_input
  implicit none
  
  call mesh_data
  call physparam_data
  call knobs_data
  call create_array
  call save_param
end subroutine read_input

!-----------------------------------------------------------------
!  De-allocate dimensions of arrays
!-----------------------------------------------------------------
subroutine delete_array
  use prec_const
  use globals

  deallocate(Dens)
  deallocate(Pres)

  deallocate(DDX0)
  deallocate(DDX1)
  deallocate(CCX0)
  deallocate(CCX1)

  deallocate(Sn)
  deallocate(Sp)

  deallocate(diagDens)
  deallocate(diagPres)
  deallocate(diagdxDens)
  deallocate(diagdxPres)
  deallocate(diagd2xPres)
  deallocate(diagd3xPres)
  deallocate(diagDntot)
  deallocate(diagChiptot)
  
end subroutine delete_array
