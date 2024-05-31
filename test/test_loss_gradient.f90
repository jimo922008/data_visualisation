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

      real(dp) :: exageration, cost, expected_cost
      real(dp), dimension(low_dimension, number_points) :: gradient_matrix, expected_gradient_matrix
      real(dp), dimension(low_dimension, number_points) :: low_dimension_position
      real(dp), dimension(number_points, number_points) :: pij
      integer, dimension(number_points) :: point_count
      real(dp), dimension(number_points) :: point_radius
      real(dp) :: core_strength

      exageration = 1.0_dp
      low_dimension_position = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]
      pij = [0.0_dp, 0.1_dp, 0.2_dp, 0.1_dp, 0.0_dp, 0.3_dp]
      point_count = [1, 1, 1, 1, 1, 1]
      point_radius = [0.5_dp, 0.5_dp, 0.5_dp]
      core_strength = 1.0_dp

      ! Expected results (to be calculated or known from correct implementation)
      expected_cost = -2.0_dp  ! Example value, you need to calculate or know this
      expected_gradient_matrix = reshape([1.0_dp, 2.0_dp, -1.0_dp, -2.0_dp, 0.0_dp, 0.0_dp], shape(gradient_matrix))

      ! Call the subroutine
      call loss_gradient(gradient_matrix, cost, exageration)

      ! Check the results
      if (abs(cost - expected_cost) < 1.0e-6_dp) then
         print *, 'Cost test passed.'
      else
         print *, 'Cost test failed. Expected:', expected_cost, 'Got:', cost
      end if

      if (all(abs(gradient_matrix - expected_gradient_matrix) < 1.0e-6_dp)) then
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
