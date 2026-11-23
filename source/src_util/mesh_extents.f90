!=======================================================================
! mesh_extents: global bounding box of the CFD mesh.
!
! global_mesh_extents returns the extents of the finest mesh level over
! all worker ranks (collective on MPI_COMM_WORLD; the master rank 0 holds
! no mesh and contributes neutral elements). Used to hand the CFD column
! to the rigid-body library, whose parallel-mode setups decompose exactly
! this box (pe/interface/setup/dkt.h). Same reduction as the
! get_global_domain_extents of the EL applications.
!=======================================================================
module mesh_extents

  implicit none
  private
  public :: global_mesh_extents

contains

  subroutine global_mesh_extents(xmin, xmax, ymin, ymax, zmin, zmax)
    use def_FEAT
    use PP3D_MPI, only: myid, showid
    use var_QuadScalar, only: mg_mesh
    implicit none
    real*8, intent(out) :: xmin, xmax, ymin, ymax, zmin, zmax
    real*8 :: xmin_local, xmax_local, ymin_local, ymax_local
    real*8 :: zmin_local, zmax_local
    integer :: ierr
    include 'mpif.h'

    if (myid .ne. 0) then
      xmin_local = minval(mg_mesh%level(NLMAX)%dcorvg(1,:))
      xmax_local = maxval(mg_mesh%level(NLMAX)%dcorvg(1,:))
      ymin_local = minval(mg_mesh%level(NLMAX)%dcorvg(2,:))
      ymax_local = maxval(mg_mesh%level(NLMAX)%dcorvg(2,:))
      zmin_local = minval(mg_mesh%level(NLMAX)%dcorvg(3,:))
      zmax_local = maxval(mg_mesh%level(NLMAX)%dcorvg(3,:))
    else
      xmin_local =  huge(1.0d0)
      ymin_local =  huge(1.0d0)
      zmin_local =  huge(1.0d0)
      xmax_local = -huge(1.0d0)
      ymax_local = -huge(1.0d0)
      zmax_local = -huge(1.0d0)
    end if

    call MPI_Allreduce(xmin_local, xmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(xmax_local, xmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(ymin_local, ymin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(ymax_local, ymax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(zmin_local, zmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(zmax_local, zmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    if (myid .eq. showid) then
      write(mterm, '(A,6ES14.6)') 'Global CFD domain extents: ', xmin, xmax, ymin, ymax, zmin, zmax
      write(mfile1, '(A,6ES14.6)') 'Global CFD domain extents: ', xmin, xmax, ymin, ymax, zmin, zmax
    end if

  end subroutine global_mesh_extents

end module mesh_extents
