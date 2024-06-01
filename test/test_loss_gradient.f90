module test_loss_gradient_mod
   implicit none
   integer, parameter :: low_dimension = 2, number_points = 3

contains

! dataset - store answer (LJ13, 100 points)
! benchmarking: test_quick, test-full
! memory usage (ram- benchmarking- instrument)
! timing profiling (timing tool)

   subroutine test_loss_gradient()
      implicit none

      real(sp) :: exaggeration, cost, expected_cost
      real(sp), dimension(low_dimension, number_points) :: gradient_matrix, expected_gradient_matrix
      real(sp), dimension(low_dimension, number_points) :: low_dimension_position
      real(sp), dimension(number_points, number_points) :: pij
      integer, dimension(number_points) :: point_count
      real(sp), dimension(number_points) :: point_radius
      real(sp) :: core_strength

      exaggeration = 1.0_sp
      low_dimension_position = [0.0_sp, 0.0_sp, 1.0_sp, 1.0_sp, 2.0_sp, 2.0_sp]
      pij = [0.0_sp, 0.1_sp, 0.2_sp, 0.1_sp, 0.0_sp, 0.3_sp]
      point_count = [1, 1, 1, 1, 1, 1]
      point_radius = [0.5_sp, 0.5_sp, 0.5_sp]
      core_strength = 1.0_sp

      ! Expected results (to be calculated or known from correct implementation)
      expected_cost = -2.0_sp  ! Example value, you need to calculate or know this
      expected_gradient_matrix = reshape([1.0_sp, 2.0_sp, -1.0_sp, -2.0_sp, 0.0_sp, 0.0_sp], shape(gradient_matrix))

      ! Call the subroutine
      call loss_gradient(gradient_matrix, cost, exaggeration)

      ! Check the results
      if (abs(cost - expected_cost) < 1.0e-6_sp) then
         print *, 'Cost test passed.'
      else
         print *, 'Cost test failed. Expected:', expected_cost, 'Got:', cost
      end if

      if (all(abs(gradient_matrix - expected_gradient_matrix) < 1.0e-6_sp)) then
         print *, 'Gradient matrix test passed.'
      else
         print *, 'Gradient matrix test failed.'
         print *, 'Expected:', expected_gradient_matrix
         print *, 'Got:', gradient_matrix
      end if
   end subroutine test_loss_gradient

end module test_loss_gradient_mod

program test_loss_gradient_program
   use test_loss_gradient_mod
   implicit none

   call test_loss_gradient()

end program test_loss_gradient_program
