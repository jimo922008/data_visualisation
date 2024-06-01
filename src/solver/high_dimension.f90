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

   real(sp), allocatable, dimension(:, :)   :: high_dist_matrix(:, :)
   real(sp), allocatable, dimension(:)      :: sigma(:)
   real(sp), allocatable, dimension(:, :)   :: pij(:, :)

contains

   subroutine high_dimension_distribution()
      implicit none
      integer :: i, j

      write (*, *) 'calculating high dimension distance......'

      call high_dimension_distance()

      call remove_duplicates(data_vec, high_dist_matrix, number_points, 0.5_sp, .10000E+09_sp)

      write (*, *) 'finding sigma......'

      call find_sigma(perplexity)

      write (*, *) 'calculating pij.......'

      allocate (pij(reduced_number_points, reduced_number_points))

      write (*, *) 'allocated pij'

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

   end subroutine high_dimension_distribution

   subroutine high_dimension_distance()
      implicit none
      integer :: i
      real(sp), allocatable, dimension(:)   :: length_2(:)

      allocate (length_2(number_points))
      allocate (high_dist_matrix(number_points, number_points))

      length_2 = sum(data_vec*data_vec, dim=1)

      high_dist_matrix = 0.0_sp
      call sgemm('T', 'N', number_points, number_points, number_features, -2.0_sp, data_vec(:, :), &
                 number_features, data_vec(:, :), number_features, 1.0_sp, high_dist_matrix, number_points)

      !$omp parallel do
      do i = 1, number_points
         high_dist_matrix(i, :) = high_dist_matrix(i, :) + length_2 + length_2(i) + 1e-10_sp
         high_dist_matrix(i, i) = 0.0_sp
      end do
      !$end parallel do

      deallocate (length_2)

   end subroutine high_dimension_distance

   subroutine find_sigma(perplexity)

      real(sp), intent(in) :: perplexity
      real(sp) :: low_sigma, high_sigma
      real(sp) :: low_perplexity, high_perplexity
      real(sp) :: tolerance, factor
      integer  :: i

      allocate (sigma(reduced_number_points))

      tolerance = 1.0e-13_sp

      !$omp parallel do private(i, high_sigma, high_perplexity, low_sigma, low_perplexity, factor) schedule(dynamic)
      do i = 1, reduced_number_points

         high_sigma = 1.0_sp
         high_perplexity = calculating_perplexity(high_sigma, i)

         factor = merge(1.0_sp/gr, gr, high_perplexity > perplexity)

         low_sigma = high_sigma*factor
         low_perplexity = calculating_perplexity(low_sigma, i)

         do while ((low_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_sp)
            high_sigma = low_sigma
            high_perplexity = low_perplexity
            low_sigma = low_sigma*factor
            low_perplexity = calculating_perplexity(low_sigma, i)
         end do

         !call bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma(i))
         call brent_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma(i))
      end do
      !$omp end parallel do

   end subroutine find_sigma

   function calculating_perplexity(sigma, i) result(perplexity)
      implicit none
      real(sp), intent(in) :: sigma
      integer, intent(in) :: i
      integer              :: j
      real(sp)             :: pji_conditional(reduced_number_points)
      real(sp)             :: entropy, perplexity

      entropy = 0.0_sp
      perplexity = 0.0_sp
      pji_conditional = 0.0_sp

      pji_conditional = calculating_pji(sigma, i)

      !$omp parallel do private(j) reduction(+:entropy)
      do j = 1, reduced_number_points
         if ((i == j) .or. (pji_conditional(j) < tiny(1.0_sp))) cycle
         entropy = entropy - pji_conditional(j)*log(pji_conditional(j))/log(2.0_sp)
      end do
      !$omp end parallel do

      perplexity = 2.0_sp**(entropy)

   end function calculating_perplexity

   function calculating_pji(sigma, i) result(pji_conditional)
      implicit none
      real(kind=sp)             :: pji_conditional(reduced_number_points)
      real(kind=sp), intent(in) :: sigma
      integer, intent(in) :: i
      integer              :: j

      !$omp parallel do private(j)
      do j = 1, reduced_number_points
         if (i == j) cycle
         pji_conditional(j) = exp(-high_dist_matrix_clean(j, i)/(sigma*sigma*2.0_sp))
      end do
      !$omp end parallel do

      pji_conditional(:) = pji_conditional(:)/sum(pji_conditional)

   end function calculating_pji

   subroutine bisection_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma)

      integer, intent(in) :: i
      real(sp), intent(in) :: perplexity, tolerance
      real(sp), intent(inout) :: low_sigma, high_sigma
      real(sp), intent(out) :: sigma
      real(sp) :: mid_sigma, low_perplexity, mid_perplexity, high_perplexity
      logical :: condition

      low_perplexity = calculating_perplexity(low_sigma, i)
      high_perplexity = calculating_perplexity(high_sigma, i)
      mid_perplexity = huge(1.0_sp)

      do while ((abs(high_sigma - low_sigma) > tolerance) .and. (abs(mid_perplexity - perplexity) > tolerance))
         mid_sigma = (low_sigma + high_sigma)/2.0_sp
         if ((abs(mid_sigma - low_sigma) < tolerance) .or. (abs(mid_sigma - high_sigma) < tolerance)) exit

         mid_perplexity = calculating_perplexity(mid_sigma, i)

         condition = (mid_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_sp

         high_sigma = merge(mid_sigma, high_sigma, condition)
         high_perplexity = merge(mid_perplexity, high_perplexity, condition)
         low_sigma = merge(low_sigma, mid_sigma, condition)
         low_perplexity = merge(low_perplexity, mid_perplexity, condition)
      end do

      sigma = mid_sigma

   end subroutine bisection_method

   subroutine brent_method(i, perplexity, low_sigma, high_sigma, tolerance, sigma)
      implicit none
      integer, intent(in) :: i
      real(sp), intent(in) :: perplexity, tolerance
      real(sp), intent(inout) :: low_sigma, high_sigma
      real(sp), intent(out) :: sigma

      real(sp) :: mid_sigma, d, e, low_perp, high_perp, mid_perp, tol1, p, low_ratio, high_ratio, ratio

      low_perp = calculating_perplexity(low_sigma, i) - perplexity
      high_perp = calculating_perplexity(high_sigma, i) - perplexity

      mid_sigma = low_sigma
      d = high_sigma - low_sigma
      e = d
      mid_perp = low_perp

      do while (abs(high_sigma - low_sigma) > tolerance)
         if (high_perp*mid_perp > 0.0_sp) then
            mid_sigma = low_sigma
            mid_perp = low_perp
            d = high_sigma - low_sigma
            e = d
         end if

         if (abs(mid_perp) < abs(high_perp)) then
            low_sigma = high_sigma
            high_sigma = mid_sigma
            mid_sigma = low_sigma
            low_perp = high_perp
            high_perp = mid_perp
            mid_perp = low_perp
         end if

         tol1 = 2.0_sp*1.0e-10_sp*abs(high_sigma) + 0.5_sp*tolerance

         if (abs(0.5_sp*(mid_sigma - high_sigma)) <= tol1 .or. high_perp == 0.0_sp) exit

         if (abs(e) >= tol1 .and. abs(low_perp) > abs(high_perp)) then
            ratio = high_perp/low_perp
            if (low_sigma == mid_sigma) then
               p = 2.0_sp*(0.5_sp*(mid_sigma - high_sigma))*ratio
               low_ratio = 1.0_sp - ratio
            else
               low_ratio = low_perp/mid_perp
               high_ratio = high_perp/mid_perp
                p = ratio * (2.0_sp * (0.5_sp * (mid_sigma - high_sigma)) * low_ratio * (low_ratio - high_ratio) - (high_sigma - low_sigma) * (high_ratio - 1.0_sp))
               low_ratio = (low_ratio - 1.0_sp)*(high_ratio - 1.0_sp)*(ratio - 1.0_sp)
            end if
            if (p > 0.0_sp) low_ratio = -low_ratio
            p = abs(p)
            if (2.0_sp*p < min(3.0_sp*(0.5_sp*(mid_sigma - high_sigma))*low_ratio - abs(tol1*low_ratio), abs(e*low_ratio))) then
               e = d
               d = p/low_ratio
            else
               d = (0.5_sp*(mid_sigma - high_sigma))
               e = d
            end if
         else
            d = (0.5_sp*(mid_sigma - high_sigma))
            e = d
         end if

         low_sigma = high_sigma
         low_perp = high_perp

         if (abs(d) > tol1) then
            high_sigma = high_sigma + d
         else
            high_sigma = high_sigma + sign(tol1, abs(0.5_sp*(mid_sigma - high_sigma)))
         end if

         high_perp = calculating_perplexity(high_sigma, i) - perplexity
      end do

      sigma = high_sigma
   end subroutine brent_method

END MODULE high_dimension
