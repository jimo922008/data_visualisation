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

   real(dp), allocatable, dimension(:, :)   :: high_dist_matrix(:, :)
   real(dp), allocatable, dimension(:)      :: sigma(:)
   real(dp), allocatable, dimension(:, :)   :: pij(:, :)

contains

   subroutine high_dimension_distribution()
      implicit none
      integer :: i, j

      write (*, *) 'calculating high dimension distance......'

      call high_dimension_distance()

      call remove_duplicates(data_vec, high_dist_matrix, number_points, 0.5_dp, .10000E+09_dp)

      write (*, *) 'finding sigma......'

      call find_sigma(perplexity)

      write (*, *) 'calculating pij.......'

      allocate (pij(reduced_number_points, reduced_number_points))

      write (*, *) 'allocated pij'

      do i = 1, reduced_number_points
         do j = 1, reduced_number_points
            if (i == j) cycle
            pij(j, i) = exp(-high_dist_matrix_clean(j, i)/(sigma(i)*sigma(i)*2.0_dp))
         end do
         pij(:, i) = pij(:, i)/sum(pij(:, i))
      end do

      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
            pij(i, j) = (pij(i, j) + pij(j, i))/real(2.0_dp*reduced_number_points)
            pij(j, i) = pij(i, j)
         end do
      end do

   end subroutine high_dimension_distribution

   subroutine high_dimension_distance()
      implicit none
      integer :: i
      real(dp), allocatable, dimension(:)   :: length_2(:)

      allocate (length_2(number_points), &
                high_dist_matrix(number_points, number_points))

      length_2 = sum(data_vec*data_vec, dim=1)

      high_dist_matrix = spread(length_2, 2, number_points) + spread(length_2, 1, number_points) + 1e-10_dp

      call dgemm('T', 'N', number_points, number_points, number_features, -2.0_dp, data_vec(:, :), &
                 number_features, data_vec(:, :), number_features, 1.0_dp, high_dist_matrix, number_points)

      !$omp parallel do
      do i = 1, number_points
         high_dist_matrix(i, i) = 0.0_dp
      end do
      !$omp end parallel do

      deallocate (length_2)

   end subroutine high_dimension_distance

   subroutine find_sigma(perplexity)

      real(dp), intent(in) :: perplexity
      real(dp) :: low_sigma, high_sigma
      real(dp) :: low_perplexity, high_perplexity
      real(dp) :: tolerance, factor
      integer  :: i

      allocate (sigma(reduced_number_points))

      tolerance = 1.0e-13_dp

      !$omp parallel do private(i, high_sigma, high_perplexity, low_sigma, low_perplexity, factor) schedule(dynamic)
      do i = 1, reduced_number_points

         high_sigma = 1.0_dp
         high_perplexity = calculating_perplexity(high_sigma, i)

         factor = merge(1.0_dp/gr, gr, high_perplexity > perplexity)

         low_sigma = high_sigma*factor
         low_perplexity = calculating_perplexity(low_sigma, i)

         do while ((low_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_dp)
            high_sigma = low_sigma
            high_perplexity = low_perplexity
            low_sigma = low_sigma*factor
            low_perplexity = calculating_perplexity(low_sigma, i)
         end do

         call bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma(i))

      end do
      !$omp end parallel do

   end subroutine find_sigma

   function calculating_perplexity(sigma, i) result(perplexity)
      implicit none
      real(dp), intent(in) :: sigma
      integer, intent(in) :: i
      integer              :: j
      real(dp)             :: pji_conditional(reduced_number_points)
      real(dp)             :: entropy, perplexity

      entropy = 0.0_dp
      perplexity = 0.0_dp
      pji_conditional = 0.0_dp

      pji_conditional = calculating_pji(sigma, i)

      !$omp parallel do private(j) reduction(+:entropy)
      do j = 1, reduced_number_points
         if ((i == j) .or. (pji_conditional(j) < tiny(1.0_dp))) cycle
         entropy = entropy - pji_conditional(j)*log(pji_conditional(j))/log(2.0_dp)
      end do
      !$omp end parallel do

      perplexity = 2.0_dp**(entropy)

   end function calculating_perplexity

   function calculating_pji(sigma, i) result(pji_conditional)
      implicit none
      real(kind=dp)             :: pji_conditional(reduced_number_points)
      real(kind=dp), intent(in) :: sigma
      integer, intent(in) :: i
      integer              :: j

      !$omp parallel do private(j)
      do j = 1, reduced_number_points
         if (i == j) cycle
         pji_conditional(j) = exp(-high_dist_matrix_clean(j, i)/(sigma*sigma*2.0_dp))
      end do
      !$omp end parallel do

      pji_conditional(:) = pji_conditional(:)/sum(pji_conditional)

   end function calculating_pji

   subroutine bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma)

      integer, intent(in) :: i
      real(dp), intent(in) :: perplexity, tolerance
      real(dp), intent(inout) :: low_sigma, high_sigma
      real(dp), intent(out) :: sigma
      real(dp) :: mid_sigma, low_perplexity, mid_perplexity, high_perplexity
      logical :: condition

      low_perplexity = calculating_perplexity(low_sigma, i)
      high_perplexity = calculating_perplexity(high_sigma, i)
      mid_perplexity = huge(1.0_dp)

      do while ((abs(high_sigma - low_sigma) > tolerance) .and. (abs(mid_perplexity - perplexity) > tolerance))
         mid_sigma = (low_sigma + high_sigma)/2.0_dp
         if ((abs(mid_sigma - low_sigma) < tolerance) .or. (abs(mid_sigma - high_sigma) < tolerance)) exit

         mid_perplexity = calculating_perplexity(mid_sigma, i)

         condition = (mid_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_dp

         high_sigma = merge(mid_sigma, high_sigma, condition)
         high_perplexity = merge(mid_perplexity, high_perplexity, condition)
         low_sigma = merge(low_sigma, mid_sigma, condition)
         low_perplexity = merge(low_perplexity, mid_perplexity, condition)
      end do

      sigma = mid_sigma

   end subroutine bisection_method

END MODULE high_dimension
