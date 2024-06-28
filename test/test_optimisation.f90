program test_optimisation

   USE parameters
   USE data_reader
   USE high_dimension
   USE low_dimension_probability
   USE optimisation
   implicit none

   TYPE(high_dim_results)     :: results
   TYPE(low_dim_results)      :: low_results
   TYPE(optimisation_parameters) :: optimisation_params
   TYPE(high_dim_parameters)  :: high_dim_params
   TYPE(low_dim_parameters)   :: low_dim_params
   TYPE(file_data)            :: data

   real(kind=sp), allocatable :: data_vec(:, :)
   real(kind=sp), allocatable :: high_dist_matrix(:, :)
   real(kind=sp), allocatable :: high_dist_matrix_clean(:, :)
   real(kind=sp), allocatable :: sigma(:)
   real(kind=sp)              :: similar_threshold, energy_threshold
   integer                    :: number_points, number_features, reduced_number_points
   real(kind=sp)              :: perplexity
   real(kind=sp), allocatable :: point_radius(:)

   similar_threshold = high_dim_params%similar_threshold
   energy_threshold = high_dim_params%energy_threshold
   perplexity = high_dim_params%perplexity

   call run_tests

contains

   subroutine run_tests
      call test_read_file(data, trim('test/LJ13.vec'))
      call test_high_dimension_distance(data, results)
      call normalisation(data)
      call test_remove_duplicates(data, results, high_dim_params)
      call test_find_sigma(results, high_dim_params)
      call test_high_dimension_distribution(data, results, high_dim_params)
      call test_hard_sphere_initialisation(results, low_results, low_dim_params)

   end subroutine run_tests

   subroutine test_read_file(data, file_name)

      type(file_data), intent(inout) :: data
      character(len=*), intent(in)  :: file_name

      call read_file(file_name, data)

      ! Check the values read from the file
      if (data%data_vec(1, 1) /= 0.0 .or. data%data_vec(1, 2) /= 0.0 .or. data%data_vec(1, 3) /= 0.0) then
         print *, 'Test failed: data_vec wrong.'
         stop
      else if (data%number_points /= 10000) then
         print *, 'Test failed: number_points /= 10000'
         stop
      else
         print *, '- read_file test PASSED.'
      end if

   end subroutine test_read_file

   subroutine test_remove_duplicates(data, results, high_dim_params)

      TYPE(high_dim_results), intent(inout) :: results
      TYPE(file_data), intent(inout) :: data
      TYPE(high_dim_parameters), intent(in) :: high_dim_params

      call remove_duplicates(data, results, high_dim_params)
      if (results%reduced_number_points /= 793) then
         print *, 'Test failed: reduced_number_points /= 793'
         stop
      else
         print *, '- remove_duplicates test PASSED.'
      end if
   end subroutine test_remove_duplicates

   subroutine test_high_dimension_distance(data, results)

      type(high_dim_results), intent(inout) :: results
      type(file_data), intent(in) :: data

      call high_dimension_distance(data, results)

      print *, results%high_dist_matrix(1, 1), results%high_dist_matrix(1, 2), results%high_dist_matrix(1, 3)

      if ((results%high_dist_matrix(1, 1) /= 0.00000000) .or. &
          (results%high_dist_matrix(1, 2) - 40.332) > 0.001 .or. &
          (results%high_dist_matrix(1, 3) - 8.930) > 0.001) then
         print *, 'Test failed: high_dist_matrix incorrect'
         stop
      else
         print *, '- high_dimension_distance test PASSED.'
      end if

   end subroutine test_high_dimension_distance

   subroutine test_find_sigma(results, high_dim_params)

      type(high_dim_results), intent(inout) :: results
      type(high_dim_parameters), intent(in) :: high_dim_params

      call find_sigma(results, high_dim_params)
      if ((results%sigma(1) - 0.631) > 0.001 .or. (results%sigma(2) - 1.102) > 0.001 .or. (results%sigma(3) - 0.636) > 0.001) then
         print *, 'Test failed: sigma incorrect'
         stop
      else
         print *, '- find_sigma test PASSED.'
      end if
   end subroutine test_find_sigma

   subroutine test_high_dimension_distribution(data, results, high_dim_params)

      type(file_data), intent(inout) :: data
      type(high_dim_results), intent(inout) :: results
      type(high_dim_parameters), intent(in) :: high_dim_params

      integer :: i, j

      allocate (results%pij(results%reduced_number_points, results%reduced_number_points))

      do i = 1, results%reduced_number_points
         do j = 1, results%reduced_number_points
            if (i == j) cycle
            results%pij(j, i) = exp(-results%high_dist_matrix_clean(j, i)/(results%sigma(i)*results%sigma(i)*2.0_sp))
         end do
         results%pij(:, i) = results%pij(:, i)/sum(results%pij(:, i))
      end do

      do i = 1, results%reduced_number_points
         do j = i + 1, results%reduced_number_points
            results%pij(i, j) = (results%pij(i, j) + results%pij(j, i))/real(2.0_sp*results%reduced_number_points)
            results%pij(j, i) = results%pij(i, j)
         end do
      end do

  if ((results%pij(1, 1)) /= 0.00000000 .or. (results%pij(1, 2)-4.087E-09>0.01E-9) .or. (results%pij(1, 3) - 1.86E-08>0.01E-8)) then
         print *, '-FAILED: pij incorrect'
         stop
      else
         print *, '-PASSED: high_dimension_distribution test.'
      end if
   end subroutine test_high_dimension_distribution

   subroutine test_hard_sphere_initialisation(results, low_results, low_dim_params)

      type(high_dim_results), intent(in) :: results
      type(low_dim_results), intent(inout) :: low_results
      type(low_dim_parameters), intent(in) :: low_dim_params

      call hard_sphere_initialisation(results, low_results, low_dim_params)

      if ((sum(low_results%point_radius) - 21.09 > 0.01) .or. &
          (low_results%point_radius(1) - 2.99E-02 > 0.01E-2) .or. &
          (low_results%point_radius(2) - 0.19 > 0.01)) then
         print *, '-FAILED: hard_sphere_initialisation test.'
         stop
      else
         print *, '-PASSED: hard_sphere_initialisation test.'
      end if

   end subroutine test_hard_sphere_initialisation

end program test_optimisation
