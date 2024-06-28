!> @file high_dimension.f90
!> @brief This file contains the module for the high dimension distance calculation and the perplexity calculation.
!> @details The high dimension distance calculation is done using the matrix multiplication and the perplexity calculation is done using the bisection method.
!! The high dimension distance matrix is calculated using the formula -2*data_vec*data_vec' + sum(data_vec*data_vec, dim=1) + sum(data_vec*data_vec, dim=1)' + 1e-10.
!! The sigma calculation is done using the bisection method. The pij calculation is done using the high dimension distance matrix and the sigma value. sgmm is used for the matrix multiplication.
!> @param(in) data_vec The data vector.
!> @param(in) number_points The number of points.
!> @param(in) number_features The number of features.
!> @param(in) perplexity The perplexity value.
!> @param(out) reduced_number_points The reduced number of points.
!> @param(out) high_dist_matrix The high dimension distance matrix.
!> @param(out) sigma The sigma value.
!> @param(out) pij The pij value.

MODULE high_dimension

   USE omp_lib
   USE parameters
   USE data_reader
   USE initialisation

   IMPLICIT NONE

   PUBLIC :: high_dimension_distance
   PUBLIC :: find_sigma
   PUBLIC :: calculating_perplexity
   PUBLIC :: calculating_pji

contains

   subroutine high_dimension_distribution(data, results, high_dim_params)
      !> @brief This subroutine calculates the pij matrix.
      !> @details The high dimension distance matrix is calculated using the formula -2*data_vec*data_vec' + sum(data_vec*data_vec, dim=1) + sum(data_vec*data_vec, dim=1)' + 1e-10.
      !! The sigma calculation is done using the bisection method. The pij calculation is done using the high dimension distance matrix and the sigma value.
      !! sgmm is used for the matrix multiplication.
      !! The pij matrix is calculated using the formula exp(-high_dist_matrix_clean(j, i)/(sigma(i)*sigma(i)*2.0_sp)). The pij matrix is then normalised.
      !> @param[in] data The file_data structure containing the data vector and other information.
      !> @param[inout] results The high_dim_results structure containing the high dimension distance matrix, sigma, distribution and other information.
      !> @param[in] high_dim_params The high_dim_params structure containing high dimension parameters.

      IMPLICIT NONE

      ! Input arguments
      TYPE(file_data), intent(inout)               :: data
      TYPE(high_dim_results), intent(inout)        :: results
      TYPE(high_dim_parameters), intent(in)        :: high_dim_params

      ! Local variables
      integer :: i, j

      write (*, *) 'calculating high dimension distance......'

      call high_dimension_distance(data, results)

      call remove_duplicates(data, results, high_dim_params)

      write (*, *) 'finding sigma......'

      call find_sigma(results, high_dim_params)

      write (*, *) 'calculating pij.......'

      allocate (results%pij(results%reduced_number_points, results%reduced_number_points))

      write (*, *) 'allocated pij'

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

   end subroutine high_dimension_distribution

   subroutine high_dimension_distance(data, results)

      !> @brief This subroutine calculates the high dimension distance matrix.
      !> @param[in] data The file_data structure containing the data vector and other information.
      !> @param[inout] results The high_dim_results structure containing the high dimension distance matrix.

      implicit none

      ! Input arguments
      TYPE(file_data), intent(in)           :: data
      TYPE(high_dim_results), intent(inout) :: results

      ! Local variables
      integer                                :: i
      real(sp), allocatable                  :: length_2(:)

      allocate (results%high_dist_matrix(data%number_points, data%number_points))
      allocate (length_2(data%number_points))

      length_2 = sum(data%data_vec*data%data_vec, dim=1)
      results%high_dist_matrix = 0.0_sp

      call sgemm('T', 'N', data%number_points, data%number_points, data%number_features, -2.0_sp, data%data_vec(:, :), &
              data%number_features, data%data_vec(:, :), data%number_features, 1.0_sp, results%high_dist_matrix, data%number_points)

      !$omp parallel do
      do i = 1, data%number_points
         results%high_dist_matrix(i, :) = results%high_dist_matrix(i, :) + length_2 + length_2(i) + 1e-10_sp
         results%high_dist_matrix(i, i) = 0.0_sp
      end do
      !$end parallel do

      deallocate (length_2)

   end subroutine high_dimension_distance

   subroutine find_sigma(results, high_dim_params)

      !> @brief This subroutine calculates the sigma value.
      !> @details The sigma calculation is done using the bisection method.
      !> @param[in] data The file_data structure containing the data vector and other information.
      !> @param[inout] results The high_dim_results structure containing the high dimension distance matrix.
      !> @param[in] high_dim_params The high_dim_params structure containing high dimension parameters.

      implicit none

      ! Input arguments
      type(high_dim_parameters), intent(in)     :: high_dim_params
      type(high_dim_results), intent(inout)     :: results

      ! Local variables
      real(sp), allocatable :: sigma(:)
      real(sp), allocatable :: high_dist_matrix(:, :)
      real(sp) :: low_sigma, high_sigma
      real(sp) :: low_perplexity, high_perplexity, perplexity, gr
      real(sp) :: tolerance
      real(sp) :: factor
      integer  :: i, index

      allocate (results%sigma(results%reduced_number_points))
      allocate (sigma(results%reduced_number_points))
      allocate (high_dist_matrix(results%reduced_number_points, results%reduced_number_points))

      tolerance = high_dim_params%find_sigma_tolerance
      index = results%reduced_number_points
      perplexity = high_dim_params%perplexity
      gr = high_dim_params%gr
      high_dist_matrix = results%high_dist_matrix_clean

      !$omp parallel do private(i, high_sigma, high_perplexity, low_sigma, low_perplexity, factor) schedule(dynamic)
      do i = 1, index

         high_sigma = 1.0_sp
         high_perplexity = calculating_perplexity(high_sigma, i, high_dist_matrix)

         factor = merge(1.0_sp/gr, gr, high_perplexity > perplexity)

         low_sigma = high_sigma*factor
         low_perplexity = calculating_perplexity(low_sigma, i, high_dist_matrix)

         do while ((low_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_sp)
            high_sigma = low_sigma
            high_perplexity = low_perplexity
            low_sigma = low_sigma*factor
            low_perplexity = calculating_perplexity(low_sigma, i, high_dist_matrix)
         end do

         call bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, high_dist_matrix, sigma(i))
         print *, 'sigma', i, sigma(i)
      end do
      !$omp end parallel do

      results%sigma = sigma

   end subroutine find_sigma

   function calculating_perplexity(sigma, i, high_dist_matrix) result(perplexity)
      implicit none

      ! Input arguments
      real(sp), intent(in)               :: sigma
      integer, intent(in)                :: i
      real(sp), intent(in)               :: high_dist_matrix(:, :)

      ! Local variables
      integer                            :: j, index
      real(sp), allocatable              :: pji_conditional(:)
      real(sp)                           :: entropy, perplexity

      allocate (pji_conditional(size(high_dist_matrix, 1)))

      entropy = 0.0_sp
      perplexity = 0.0_sp
      pji_conditional = calculating_pji(sigma, i, high_dist_matrix)
      index = size(high_dist_matrix, 1)

      !omp parallel do private(j) reduction(+:entropy)
      do j = 1, index
         if ((i == j) .or. (pji_conditional(j) < tiny(1.0_sp))) cycle
         entropy = entropy - pji_conditional(j)*log(pji_conditional(j))/log(2.0_sp)
      end do
      !omp end parallel do

      perplexity = 2.0_sp**(entropy)

      deallocate (pji_conditional)

   end function calculating_perplexity

   function calculating_pji(sigma, i, high_dist_matrix) result(pji_conditional)
      implicit none

      ! Input arguments
      integer, intent(in)                :: i
      real(kind=sp), intent(in)          :: sigma
      real(kind=sp), intent(in)          :: high_dist_matrix(:, :)

      ! Output arguments
      real(kind=sp), allocatable         :: pji_conditional(:)

      ! Local variables
      integer                            :: j, index

      index = size(high_dist_matrix, 1)

      allocate (pji_conditional(size(high_dist_matrix, 1)))

      !$omp parallel do private(j)
      do j = 1, index
         if (i == j) cycle
         pji_conditional(j) = exp(-high_dist_matrix(j, i)/(sigma*sigma*2.0_sp))
      end do
      !$omp end parallel do

      pji_conditional(:) = pji_conditional(:)/sum(pji_conditional)

   end function calculating_pji

   subroutine bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, high_dist_matrix, sigma)
      implicit none

      ! Input arguments
      integer, intent(in) :: i
      real(sp), intent(in) :: perplexity, tolerance
      real(sp), intent(inout) :: sigma
      real(sp), intent(inout) :: low_sigma, high_sigma
      real(sp), intent(in) :: high_dist_matrix(:, :)

      ! Local variables
      real(sp):: low_perplexity, mid_perplexity, high_perplexity, mid_sigma
      logical :: condition

      low_perplexity = calculating_perplexity(low_sigma, i, high_dist_matrix)
      high_perplexity = calculating_perplexity(high_sigma, i, high_dist_matrix)
      mid_perplexity = huge(1.0_sp)

      do while ((abs(high_sigma - low_sigma) > tolerance) .and. (abs(mid_perplexity - perplexity) > tolerance))
         mid_sigma = (low_sigma + high_sigma)/2.0_sp
         if ((abs(mid_sigma - low_sigma) < tolerance) .or. (abs(mid_sigma - high_sigma) < tolerance)) exit

         mid_perplexity = calculating_perplexity(mid_sigma, i, high_dist_matrix)

         condition = (mid_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_sp

         high_sigma = merge(mid_sigma, high_sigma, condition)
         high_perplexity = merge(mid_perplexity, high_perplexity, condition)
         low_sigma = merge(low_sigma, mid_sigma, condition)
         low_perplexity = merge(low_perplexity, mid_perplexity, condition)
      end do

      sigma = mid_sigma

   end subroutine bisection_method

END MODULE high_dimension
