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
      call test_hard_sphere_initialisation()

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
         stop
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
         print *, '-FAILED: pij incorrect'
      else
         print *, '-PASSED: high_dimension_distribution test.'
      end if
   end subroutine test_high_dimension_distribution

   subroutine test_hard_sphere_initialisation()

      call hard_sphere_initialisation()

      if (sum(point_radius) /= 21.0950317 .or. &
          point_radius(1) /= 2.99999993E-02 .or. &
          point_radius(2) /= 0.194935888) then
         print *, '-FAILED: hard_sphere_initialisation test.'
      else
         print *, '-PASSED: hard_sphere_initialisation test.'
      end if

   end subroutine test_hard_sphere_initialisation

end program test_optimisation
