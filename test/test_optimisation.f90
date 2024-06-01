program test_optimisation

   USE parameters
   USE data_reader
   USE high_dimension
   USE low_dimension_probability
   USE optimisation
   implicit none

   call run_tests

contains

   subroutine run_tests
      call test_read_file(trim('test/LJ13.vec'))
      call test_high_dimension_distance()
      call normalisation(data_vec, number_points, number_features)
      call test_remove_duplicates(data_vec, high_dist_matrix, number_points, similar_threshold, energy_threshold)
      call test_find_sigma(perplexity)
      call test_high_dimension_distribution()

   end subroutine run_tests

   subroutine test_read_file(file_name)
      character(len=*) :: file_name

      call read_file(file_name)

      ! Check the values read from the file
      if (data_vec(1, 1) /= 0.0 .or. data_vec(1, 2) /= 0.0 .or. data_vec(1, 3) /= 0.0) then
         print *, 'Test failed: data_vec wrong.'
      else if (number_points /= 10000) then
         print *, 'Test failed: number_points /= 10000'
      else
         print *, '- read_file test PASSED.'
      end if

   end subroutine test_read_file

   subroutine test_remove_duplicates(data_vec, high_dist_matrix, number_points, similar_threshold, energy_threshold)
      REAL(kind=sp), intent(in) :: data_vec(:, :)
      REAL(kind=sp), intent(inout) :: similar_threshold, energy_threshold
      REAL(kind=sp), intent(in) :: high_dist_matrix(:, :)
      INTEGER, intent(inout)    :: number_points

      call remove_duplicates(data_vec, high_dist_matrix, number_points, similar_threshold, energy_threshold)
      if (reduced_number_points /= 793) then
         print *, 'Test failed: reduced_number_points /= 793'
      else
         print *, '- remove_duplicates test PASSED.'
      end if
   end subroutine test_remove_duplicates

   subroutine test_high_dimension_distance()
      call high_dimension_distance()

      if (high_dist_matrix(1, 1) /= 0.00000000 .or. &
          high_dist_matrix(1, 2) /= 40.3325195 .or. &
          high_dist_matrix(1, 3) /= 8.93066406) then
         print *, 'Test failed: high_dist_matrix incorrect'
      else
         print *, '- high_dimension_distance test PASSED.'
      end if

   end subroutine test_high_dimension_distance

   subroutine test_find_sigma(perplexity)
      real(kind=sp), intent(in) :: perplexity

      call find_sigma(perplexity)
      if (sigma(1) /= 0.631218076 .or. sigma(2) /= 1.10291600 .or. sigma(3) /= 0.636737883) then
         print *, 'Test failed: sigma incorrect'
      else
         print *, '- find_sigma test PASSED.'
      end if
   end subroutine test_find_sigma

   subroutine test_high_dimension_distribution()
      integer :: i, j
      real(kind=sp) :: pij(reduced_number_points, reduced_number_points)

      do i = 1, reduced_number_points
         do j = 1, reduced_number_points
            if (i == j) cycle
            pij(j, i) = exp(-high_dist_matrix_clean(j, i)/(sigma(i)*sigma(i)*2.0_sp))
         end do
         pij(:, i) = pij(:, i)/sum(pij(:, i))
      end do

      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
            pij(i, j) = (pij(i, j) + pij(j, i))/real(2.0_sp*reduced_number_points)
            pij(j, i) = pij(i, j)
         end do
      end do

      if ((pij(1, 1)) /= 0.00000000 .or. pij(1, 2) /= 4.08979650E-09 .or. pij(1, 3) /= 1.85994082E-08) then
         print *, 'Test failed: pij incorrect'
      else
         print *, '- high_dimension_distribution test PASSED.'
      end if
   end subroutine test_high_dimension_distribution

   subroutine test_initialize_variables
      ! Test the subroutine initialize_variables
      integer :: i, j
      real(kind=sp) :: gradient_norm, running_gradient_norm

      ! Initialize necessary variables
      reduced_number_points = 10  ! Example value for testing
      allocate (pij(reduced_number_points, reduced_number_points))
      pij = 0.0_sp  ! Dummy values for testing

      ! Call the subroutine to initialize variables
      call initialize_variables(i, j, gradient_norm, running_gradient_norm)

      ! Check the initialized values
      if (i /= 0) print *, 'Test failed: i /= 0'
      if (j /= 0) print *, 'Test failed: j /= 0'
      if (gradient_norm /= huge(1.0_sp)) print *, 'Test failed: gradient_norm /= huge(1.0_sp)'
      if (running_gradient_norm /= 0.0_sp) print *, 'Test failed: running_gradient_norm /= 0.0_sp'
      if (cost_zero /= 0.0_sp) print *, 'Test failed: cost_zero /= 0.0_sp'
      if (z /= real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)) print *, 'Test failed: z is incorrect'
      if (inv_z /= 1.0_sp/z) print *, 'Test failed: inv_z is incorrect'

      print *, 'All tests passed.'

   end subroutine test_initialize_variables

end program test_optimisation
