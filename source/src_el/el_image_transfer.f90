!=========================================================================
! EL_IMAGE_TRANSFER - the coupling operations of one particle record on one rank
!
! EL_PARTICLE_MESH_PASS (EL_TRANSFER) couples a particle record to the elements of a rank in
! three operations: the kernel sample of the fluid fields (EL_SAMPLE_RECORD), the additive
! pieces of the kernel-averaged Laplacian for the viscous force of a momentum model
! (EL_SAMPLE_RECORD_VISCOUS) and the deposit of solid volume, reaction force and drag
! coefficient (EL_DEPOSIT_RECORD). This module is the record level of the three: each is the
! single-kernel routine of EL_QUADRATURE / EL_VISCOUS applied to every periodic image of the
! record that reaches the rank (EL_PERIODIC: EL_IMAGE_SHIFTS), the results added.
!
! Why adding is right. All three single-kernel results are additive in the kernel: the sample
! holds kernel-weighted integrals with the plain kernel integral in slot
! EL_SAMPLE_NORMALIZATION (the owner divides after the cross-rank sum, EL_HALO), the viscous
! record is additive by construction (EL_VISCOUS), and the deposit adds
! kernel-weighted integrals divided by the normalisation the owner sends back. The kernel of
! a particle on a periodic lattice is the sum of the kernels of its images, hence the sums
! below; the normalisation then is the kernel integral over the real plus image support, and
! the weights of a particle sum to one exactly as in the interior.
!
! box is the bounding box of the elements of this rank (EL_HALO: EL_HALO_LOCAL_BOX in the
! application, the mesh box in a unit test); it only selects the images, the element loop of
! the single-kernel routine does the exact culling.
!
! With the lattice inactive (SimPar@ELPeriodicKernel = off, the default) each routine is one
! call of the single-kernel routine with the record as it is: no image list, no copy, no sum.
! An image differs from its record in the position only (velocity, radius, identity are those
! of the particle).
!=========================================================================
MODULE EL_IMAGE_TRANSFER

  USE TYPES, ONLY: tMultiMesh
  USE EL_CONFIG, ONLY: EL_KERNEL_SUPPORT_RADIUS
  USE EL_HALO, ONLY: tElParticleRecord
  USE EL_FIELDS, ONLY: tElFieldStorage
  USE EL_PERIODIC, ONLY: EL_PERIODIC_ACTIVE, EL_IMAGE_SHIFTS, EL_MAX_IMAGES
  USE EL_QUADRATURE, ONLY: EL_SAMPLE_SIZE, EL_INTEGRATE_PARTICLE, EL_DEPOSIT_PARTICLE
  USE EL_VISCOUS, ONLY: EL_INTEGRATE_PARTICLE_VISCOUS, EL_VIS_SIZE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: EL_SAMPLE_RECORD, EL_SAMPLE_RECORD_VISCOUS, EL_DEPOSIT_RECORD, EL_RECORD_IMAGES

CONTAINS

  ! The images of a record that reach the box: copies of the record at position + shift.
  SUBROUTINE EL_RECORD_IMAGES(record, box, images, n)

    TYPE(tElParticleRecord), INTENT(IN) :: record
    REAL*8, INTENT(IN) :: box(6)
    TYPE(tElParticleRecord), INTENT(OUT) :: images(EL_MAX_IMAGES)
    INTEGER, INTENT(OUT) :: n

    REAL*8 :: shifts(3,EL_MAX_IMAGES)
    INTEGER :: k

    CALL EL_IMAGE_SHIFTS(record%position, EL_KERNEL_SUPPORT_RADIUS(record%radius), box, shifts, n)
    DO k = 1, n
      images(k) = record
      images(k)%position = record%position + shifts(:,k)
    END DO

  END SUBROUTINE EL_RECORD_IMAGES

  ! Kernel sample of the fluid fields for one record (layout: EL_QUADRATURE, EL_SAMPLE_*).
  SUBROUTINE EL_SAMPLE_RECORD(mesh, ilev, velocity_u, velocity_v, velocity_w, pressure, &
                              epsilon_f, record, box, sample)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: velocity_u(:), velocity_v(:), velocity_w(:)
    REAL*8, INTENT(IN) :: pressure(:), epsilon_f(:)
    TYPE(tElParticleRecord), INTENT(IN) :: record
    REAL*8, INTENT(IN) :: box(6)
    REAL*8, INTENT(OUT) :: sample(EL_SAMPLE_SIZE)

    TYPE(tElParticleRecord) :: images(EL_MAX_IMAGES)
    REAL*8 :: part(EL_SAMPLE_SIZE)
    INTEGER :: n, k

    IF (.NOT.EL_PERIODIC_ACTIVE()) THEN
      CALL EL_INTEGRATE_PARTICLE(mesh, ilev, velocity_u, velocity_v, velocity_w, pressure, &
                                 epsilon_f, record, sample)
      RETURN
    END IF

    CALL EL_RECORD_IMAGES(record, box, images, n)
    sample = 0.0d0
    DO k = 1, n
      CALL EL_INTEGRATE_PARTICLE(mesh, ilev, velocity_u, velocity_v, velocity_w, pressure, &
                                 epsilon_f, images(k), part)
      sample = sample + part
    END DO

  END SUBROUTINE EL_SAMPLE_RECORD

  ! Additive pieces of the kernel-averaged Laplacian for one record (layout: EL_VISCOUS).
  ! boundary_face_dof marks the faces that cut the kernel; a periodic plane must not be one
  ! of them (EL_PARTICLE_MESH_PASS verifies the mask of the fluid adapter in every pass).
  SUBROUTINE EL_SAMPLE_RECORD_VISCOUS(mesh, ilev, coor, velocity_u, velocity_v, velocity_w, &
                                      boundary_face_dof, record, box, vis)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    REAL*8, INTENT(IN) :: coor(:,:), velocity_u(:), velocity_v(:), velocity_w(:)
    LOGICAL, INTENT(IN) :: boundary_face_dof(:)
    TYPE(tElParticleRecord), INTENT(IN) :: record
    REAL*8, INTENT(IN) :: box(6)
    REAL*8, INTENT(OUT) :: vis(EL_VIS_SIZE)

    TYPE(tElParticleRecord) :: images(EL_MAX_IMAGES)
    REAL*8 :: part(EL_VIS_SIZE)
    INTEGER :: n, k

    IF (.NOT.EL_PERIODIC_ACTIVE()) THEN
      CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, ilev, coor, velocity_u, velocity_v, velocity_w, &
        boundary_face_dof, record%position, EL_KERNEL_SUPPORT_RADIUS(record%radius), vis)
      RETURN
    END IF

    CALL EL_RECORD_IMAGES(record, box, images, n)
    vis = 0.0d0
    DO k = 1, n
      CALL EL_INTEGRATE_PARTICLE_VISCOUS(mesh, ilev, coor, velocity_u, velocity_v, velocity_w, &
        boundary_face_dof, images(k)%position, EL_KERNEL_SUPPORT_RADIUS(record%radius), part)
      vis = vis + part
    END DO

  END SUBROUTINE EL_SAMPLE_RECORD_VISCOUS

  ! Deposit of one record: solid volume (alpha_p and its kernel gradient), reaction force and,
  ! when present, the drag coefficient (with the HRZ weights when hrz is .TRUE., see
  ! EL_QUADRATURE:EL_DEPOSIT_PARTICLE). normalization is the kernel integral of the particle
  ! over all ranks and images, as the owner summed it.
  SUBROUTINE EL_DEPOSIT_RECORD(mesh, ilev, record, box, normalization, force, fields, drag_b, &
                               hrz)

    TYPE(tMultiMesh), INTENT(IN) :: mesh
    INTEGER, INTENT(IN) :: ilev
    TYPE(tElParticleRecord), INTENT(IN) :: record
    REAL*8, INTENT(IN) :: box(6), normalization, force(3)
    TYPE(tElFieldStorage), INTENT(INOUT) :: fields
    REAL*8, INTENT(IN), OPTIONAL :: drag_b
    LOGICAL, INTENT(IN), OPTIONAL :: hrz

    TYPE(tElParticleRecord) :: images(EL_MAX_IMAGES)
    INTEGER :: n, k

    IF (.NOT.EL_PERIODIC_ACTIVE()) THEN
      CALL EL_DEPOSIT_PARTICLE(mesh, ilev, record, normalization, force, fields, drag_b, hrz)
      RETURN
    END IF

    CALL EL_RECORD_IMAGES(record, box, images, n)
    DO k = 1, n
      CALL EL_DEPOSIT_PARTICLE(mesh, ilev, images(k), normalization, force, fields, drag_b, &
                               hrz)
    END DO

  END SUBROUTINE EL_DEPOSIT_RECORD

END MODULE EL_IMAGE_TRANSFER
