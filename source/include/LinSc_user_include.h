
interface 
  subroutine LinSc_IC(dcorvg)

    implicit none
   
    ! Coordinates of the query point 
    REAL*8, dimension(:,:), pointer :: dcorvg

  end subroutine

  subroutine LinSc_BC()

   
  end subroutine

  subroutine LinSc_SRC()


  end subroutine
end interface
