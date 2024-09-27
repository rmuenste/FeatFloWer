! interface for the fbm_update function that handles
! the forces and dynamics of the fictitious boundary
! objects
interface 
  subroutine update_fbm_handler(rho, dt, simTime, g, mfile, pid)

    implicit none
    
    real*8, intent(in) :: rho      ! fluid density
    real*8, intent(in) :: dt       ! time step
    real*8, intent(in) :: simTime  ! simulation time
    real*8, dimension(3), intent(in) :: g  ! simulation time

    integer, intent(in) :: mfile   ! prot file unit
    integer, intent(in) :: pid    ! process id
  end subroutine update_fbm_handler
end interface 

interface 
subroutine fbm_force_wrapper(u,v,w,p)
  use PP3D_MPI, only: myMPI_Barrier
  use cinterface
  
  ! U/V/W velocity components
  REAL*8, dimension(:), intent(inout) :: u,v,w
  
  ! Pressure, viscosity
  Real*8, dimension(:), intent(inout) :: p


  end subroutine fbm_force_wrapper
end interface 
