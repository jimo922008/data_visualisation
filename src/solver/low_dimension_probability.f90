MODULE low_dimension_probability

   USE parameters
   USE data_reader
   USE initialisation
   USE high_dimension

   IMPLICIT NONE
   PUBLIC :: low_dimension_distribution
   PUBLIC :: calculating_low_dimension_distance
   PUBLIC :: calculating_qij
   PUBLIC :: random_add_noise

contains

   subroutine low_dimension_distribution(data, high_dim_params, low_dim_params, results, low_results)

      !> @brief Initialising the distribution of the low dimension position
      !> @param[in] data structure containing the source data and other information.
      !> @param[in] high_dim_params structure containing the high dimension parameters.
      !> @param[in] low_dim_params structure containing the low dimension parameters.
      !> @param[inout] results structure containing the results of the high dimension.
      !> @param[inout] low_results structure containing the results of the low dimension.

      IMPLICIT NONE

      ! Input variables
      TYPE(file_data), intent(in)             :: data
      TYPE(high_dim_parameters), intent(in)   :: high_dim_params
      TYPE(low_dim_parameters), intent(in)    :: low_dim_params
      TYPE(high_dim_results), intent(inout)   :: results
      TYPE(low_dim_results), intent(inout)    :: low_results

      call low_dimension_init(data, low_dim_params, results, low_results)
      write (*, *) 'low dimension init done.'
      call low_dimension_normalisation(results, low_results, low_dim_params)
      write (*, *) 'low dimension normalisation done.'
      call hard_sphere_initialisation(results, low_results, low_dim_params)
      write (*, *) 'hard sphere initialisation done.'

   end subroutine low_dimension_distribution

   subroutine low_dimension_init(data, low_dim_params, results, low_results)

      !> @brief Initialisation of the low dimension position
      !> @param[in] data structure containing the source data and other information.
      !> @param[in] high_dim_params structure containing the high dimension parameters.
      !> @param[in] low_dim_params structure containing the low dimension parameters.
      !> @param[inout] results structure containing the results of the high dimension.
      !> @param[inout] low_results structure containing the results of the low dimension.

      IMPLICIT NONE

      ! Input variables
      TYPE(file_data), intent(in)             :: data
      TYPE(low_dim_parameters), intent(in)    :: low_dim_params
      TYPE(high_dim_results), intent(inout)   :: results
      TYPE(low_dim_results), intent(inout)    :: low_results

      ! Local variables
      REAL(kind=sp), allocatable              :: random_matrix(:, :)
      REAL(kind=sp)                           :: coeff, sigma_average
      INTEGER                                 :: i, j

      allocate (random_matrix(low_dim_params%low_dimension, data%number_features))

      random_matrix = 0.0_sp

      do i = 1, low_dim_params%low_dimension
         call random_add_noise(random_matrix(i, :), 1.0_sp)
         random_matrix(i, :) = random_matrix(i, :)/norm2(random_matrix(i, :))
         if (i > 1) then
            do j = 1, i - 1
               coeff = dot_product(random_matrix(i - j, :), random_matrix(i, :))
               random_matrix(i, :) = random_matrix(i, :) - coeff*random_matrix(i - j, :)
            end do
         end if
         random_matrix(i, :) = random_matrix(i, :)/norm2(random_matrix(i, :))
      end do

      sigma_average = sum(results%sigma)/real(results%reduced_number_points, sp)

  low_results%low_dimension_position = matmul(random_matrix, results%data_clean)/sigma_average/low_dim_params%projection_compression

      write (*, *) 'low dimension position done.'

      deallocate (random_matrix)

   end subroutine low_dimension_init

   subroutine low_dimension_init_pca(low_dim_params, results, low_results)

      !> @brief Initialisation of the low dimension position using PCA
      !> @param[in] data structure containing the source data and other information.
      !> @param[in] high_dim_params structure containing the high dimension parameters.
      !> @param[in] low_dim_params structure containing the low dimension parameters.
      !> @param[inout] results structure containing the results of the high dimension.
      !> @param[inout] low_results structure containing the results of the low dimension.

      IMPLICIT NONE

      ! Input variables
      TYPE(low_dim_parameters), intent(in)    :: low_dim_params
      TYPE(high_dim_results), intent(inout)   :: results
      TYPE(low_dim_results), intent(inout)    :: low_results

      ! Local variables
      integer :: e_num, info, high_dimension, lwork
      integer, dimension(:), allocatable :: iwork, ifail
      real(kind=dp)                              :: vl, vu, sigma_average
      real(kind=dp), allocatable, dimension(:)   :: work, e_val, work_0
      real(kind=dp), dimension(:, :), allocatable :: covariance_matrix, e_vec

      high_dimension = size(results%data_clean, 1)
      allocate (covariance_matrix(high_dimension, high_dimension))
      allocate (e_vec(high_dimension, high_dimension))
      allocate (e_val(high_dimension))
      allocate (iwork(5*high_dimension))
      allocate (ifail(high_dimension))

      sigma_average = sum(results%sigma)/real(results%reduced_number_points, sp)

      covariance_matrix = matmul(results%data_clean, transpose(results%data_clean))

      ifail = 0
      allocate (work_0(1))
      call dsyevx('V','I','U',high_dimension,covariance_matrix,high_dimension,vl,vu,high_dimension-low_dim_params%low_dimension+1,high_dimension,epsilon(1.0),&
                  e_num, e_val, e_vec, high_dimension, work_0, -1, iwork, ifail, info)

      lwork = int(work_0(1))
      allocate (work(lwork))

      write (*, *) 'lwork = ', lwork

      call dsyevx('V','I','U',high_dimension,covariance_matrix,high_dimension,vl,vu,high_dimension-low_dim_params%low_dimension+1,high_dimension,epsilon(1.0),&
                  e_num, e_val, e_vec, high_dimension, work, lwork, iwork, ifail, info)

      write (*, *) 'eigenvalues done.'

      low_results%low_dimension_position = matmul(transpose(e_vec(:, low_dim_params%low_dimension:1:-1)), results%data_clean)/low_dim_params%projection_compression/sigma_average

      write (*, *) 'low dimension position done.'

   end subroutine low_dimension_init_pca

   subroutine low_dimension_normalisation(results, low_results, low_dim_params)

      !> @brief Normalisation of the low dimension position
      !> @param[inout] low_results structure containing the results of the low dimension.

      IMPLICIT NONE

      ! Input variables
      TYPE(low_dim_parameters), intent(in)  :: low_dim_params
      TYPE(high_dim_results), intent(in)    :: results
      TYPE(low_dim_results), intent(inout)  :: low_results

      allocate (low_results%low_dimension_mean(low_dim_params%low_dimension))

      low_results%low_dimension_mean = sum(low_results%low_dimension_position, dim=2)/real(results%reduced_number_points, sp)

      low_results%low_dimension_position = low_results%low_dimension_position - spread(low_results%low_dimension_mean, 2, results%reduced_number_points)

   end subroutine low_dimension_normalisation

   subroutine calculating_low_dimension_distance(results, low_results, low_dim_params)

      !> @brief Calculation of the low dimension distance
      !> @param[in] results structure containing the results of the high dimension.
      !> @param[inout] low_results structure containing the results of the low dimension.

      IMPLICIT NONE
      ! Input variables
      TYPE(high_dim_results), intent(in)    :: results
      TYPE(low_dim_results), intent(inout)  :: low_results
      TYPE(low_dim_parameters), intent(in)  :: low_dim_params

      ! Local variables
      integer :: i, j
      real(kind=sp), allocatable            :: vec(:)

      allocate (vec(low_dim_params%low_dimension))
      allocate (low_results%low_dimension_distance(results%reduced_number_points, results%reduced_number_points))

      low_results%low_dimension_distance = 0.0_sp

      !$omp parallel do schedule(static)
      do i = 1, results%reduced_number_points
         do j = i + 1, results%reduced_number_points
            vec(:) = low_results%low_dimension_position(:, i) - low_results%low_dimension_position(:, j)
            low_results%low_dimension_distance(j, i) = dot_product(vec, vec)
            low_results%low_dimension_distance(i, j) = low_results%low_dimension_distance(j, i)
         end do
      end do
      !$omp end parallel do

   end subroutine calculating_low_dimension_distance

   function calculating_qij(i, j, results, low_results, low_dim_params) result(qij)
      !> @brief Calculation of the qij value
      !> @param[in] i index of the first point.
      !> @param[in] j index of the second point.
      !> @param[in] results structure containing the results of the high dimension.
      !> @param[in] low_results structure containing the results of the low dimension.
      !> @return qij value.

      IMPLICIT NONE

      ! Input variables
      integer, intent(in) :: i, j
      TYPE(high_dim_results), intent(in)   :: results
      TYPE(low_dim_results), intent(in)    :: low_results
      TYPE(low_dim_parameters), intent(in) :: low_dim_params

      ! Local variables
      real(kind=sp)                        :: qij
      real(kind=sp)                        :: rij2, z
      real(kind=sp), allocatable           :: pos(:)

      allocate (pos(low_dim_params%low_dimension))

      z = real(results%reduced_number_points, sp)*(real(results%reduced_number_points, sp) - 1.0_sp)
      pos(:) = low_results%low_dimension_position(:, i) - low_results%low_dimension_position(:, j)
      rij2 = dot_product(pos, pos)
      qij = 1.0_sp/(1.0_sp + rij2)/z

      deallocate (pos)

   end function calculating_qij

   subroutine random_add_noise(vec, variance)

      !> @brief Add noise to the vector
      !> @param[inout] vec vector to which the noise is added.
      !> @param[in] variance variance of the noise.

      IMPLICIT NONE

      ! Input variables
      real(kind=sp), intent(in)                  :: variance
      real(kind=sp), dimension(:), intent(inout) :: vec
      real(kind=sp), dimension(size(vec))        :: u, v

      call random_number(u)
      call random_number(v)

      vec = vec + sqrt(-2.0_sp*log(u))*cos(2.0_sp*acos(-1.0_sp)*v)*variance

   end subroutine random_add_noise

   subroutine hard_sphere_initialisation(results, low_results, low_dim_params)

      !> @brief Initialisation of the hard sphere radius
      !> @param[in] data structure containing the source data and other information.
      !> @param[in] results structure containing the results of the high dimension.
      !> @param[inout] low_results structure containing the results of the low dimension.
      !> @param[in] low_dim_params structure containing the low dimension parameters.

      IMPLICIT NONE

      ! Input variables
      TYPE(high_dim_results), intent(in)      :: results
      TYPE(low_dim_results), intent(inout)    :: low_results
      TYPE(low_dim_parameters), intent(in)    :: low_dim_params

      integer :: i

      allocate (low_results%point_radius(results%reduced_number_points))

      low_results%point_radius = 0.0_sp

      !$omp parallel do schedule(static)
      do i = 1, results%reduced_number_points
         low_results%point_radius(i) = low_dim_params%sphere_radius*real(results%point_count_clean(i), sp)**(1.0_sp/real(low_dim_params%low_dimension, sp))
      end do
      !$omp end parallel do

   end subroutine hard_sphere_initialisation

end module low_dimension_probability
