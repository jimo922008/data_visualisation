MODULE data_writer

   USE parameters
   USE data_reader
   USE high_dimension
   USE low_dimension_probability
   USE optimisation

   IMPLICIT NONE

contains
   subroutine write_file(filename, data, results, low_results, low_dim_params)
      !> @brief Write the data to a file
      !> @param filename The name of the file to write to
      !> @param data structure containing the original data to write
      !> @param low_results structure containing the results of the low dimension optimisation

      IMPLICIT NONE

      !Input variables
      character(len=*)                       :: filename
      type(file_data), intent(in)            :: data
      type(high_dim_results), intent(in)     :: results
      type(low_dim_results), intent(inout)   :: low_results
      type(low_dim_parameters), intent(in)   :: low_dim_params

      !Local variables
      character(len=20):: fmt
      integer :: unit_number
      integer :: i, j

      unit_number = 20
      open (unit=unit_number, file=trim(filename), status='replace', action='write')

      write (unit_number, '(i0)') results%reduced_number_points
      write (unit_number, '(a,i0,a,f16.10)') 'dim ', low_dim_params%low_dimension

      print *, 'Writing data to file ', trim(filename)

      call data_pca(data, results, low_results, low_dim_params)

      do i = 1, results%reduced_number_points
         write (fmt, *) low_dim_params%low_dimension
            write(unit_number,'(a2,'//trim(adjustl(fmt))//'g25.13,a50,i6,a15,a10,g25.13,g25.13,i6,g25.13)') 'H',(low_results%low_dimension_position(j,i),j=1,low_dim_params%low_dimension),&
            trim(results%ion_label_clean(i)), results%ion_num_clean(i), trim(results%ion_formula_clean(i)), trim(results%ion_symmetry_clean(i)), results%ion_volume_clean(i), results%ion_energy_clean(i), &
            results%point_count_clean(i), max(low_dim_params%sphere_radius, low_results%point_radius(i))
      end do

      flush (unit_number)
      close (unit_number)

   end subroutine write_file

   subroutine data_pca(data, results, low_results, low_dim_params)
      IMPLICIT NONE

      !Input variables
      type(file_data), intent(in)            :: data
      type(high_dim_results), intent(in)     :: results
      type(low_dim_results), intent(inout)   :: low_results
      type(low_dim_parameters), intent(in)   :: low_dim_params

      ! Local variables
      integer :: e_num, info, dimension, lwork
      integer, dimension(:), allocatable :: iwork, ifail
      real(kind=dp)                              :: vl, vu, sigma_average
      real(kind=dp), allocatable, dimension(:)   :: work, e_val, work_0
      real(kind=dp), dimension(:, :), allocatable :: covariance_matrix, e_vec

      dimension = low_dim_params%low_dimension
      allocate (covariance_matrix(dimension, dimension))
      allocate (e_vec(dimension, dimension))
      allocate (e_val(dimension))
      allocate (iwork(5*dimension))
      allocate (ifail(dimension))

      covariance_matrix = matmul(low_results%low_dimension_position, transpose(low_results%low_dimension_position))

      ifail = 0
      allocate (work_0(1))
      call dsyevx('V', 'I', 'U', dimension, covariance_matrix, dimension, vl, vu, 1, dimension, epsilon(1.0), &
                  e_num, e_val, e_vec, dimension, work_0, -1, iwork, ifail, info)

      lwork = int(work_0(1))
      allocate (work(lwork))

      write (*, *) 'lwork = ', lwork

      call dsyevx('V', 'I', 'U', dimension, covariance_matrix, dimension, vl, vu, 1, dimension, epsilon(1.0), &
                  e_num, e_val, e_vec, dimension, work, lwork, iwork, ifail, info)

      write (*, *) 'eigenvalues done.'

      low_results%low_dimension_position = matmul(transpose(e_vec(:, low_dim_params%low_dimension:1:-1)), low_results%low_dimension_position)

   end subroutine data_pca

END MODULE data_writer
