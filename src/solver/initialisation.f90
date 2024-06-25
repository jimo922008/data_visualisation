MODULE initialisation

   USE parameters
   USE data_reader

   IMPLICIT NONE

   PUBLIC :: normalisation
   PUBLIC :: remove_duplicates

contains
   subroutine normalisation(data)

      !> @brief Normalise the source data.
      !> @param[inout] data The file_data structure containing the source data.
      !> @details The source data is normalised by subtracting the mean and dividing by the standard deviation.

      ! Input arguments
      TYPE(file_data), intent(inout)  :: data

      allocate (data%data_std(data%number_features))
      allocate (data%data_mean(data%number_features))

      WRITE (*, *) 'Centering source data.'

      data%data_mean = sum(data%data_vec, dim=2)/real(data%number_points, sp)

      data%data_vec = data%data_vec - spread(data%data_mean, 2, data%number_points)

      data%data_std = sqrt(sum(data%data_vec**2, dim=2)/real(data%number_points - 1, sp))

      write (*, *) 'Source data normalised.'

   end subroutine normalisation

   subroutine remove_duplicates(data, results, high_dim_params)

      !> @brief Remove duplicate points from the source data.
      !> @param[in] data The source data.
      !> @param[in] high_dist_matrix The high dimensional distance matrix.

      ! Input arguments
      TYPE(file_data), intent(inout)            :: data
      TYPE(high_dim_results), intent(inout)     :: results
      TYPE(high_dim_parameters), intent(in)     :: high_dim_params

      ! Local arguments
      INTEGER, allocatable       :: valid_points(:)
      REAL(kind=sp)              :: s_threshold, e_threshold
      INTEGER                    :: i, j, l, m

      s_threshold = (high_dim_params%similar_threshold*sum(data%data_std)/100_sp)**2
      e_threshold = high_dim_params%energy_threshold

      allocate (valid_points(count(data%point_count /= 0)))

      valid_points = pack([(i, i=1, data%number_points)], data%point_count /= 0)

      do l = 1, size(valid_points)
         i = valid_points(l)
         do m = l + 1, size(valid_points)
            j = valid_points(m)
         if ((results%high_dist_matrix(j, i) < s_threshold) .and. (abs(data%ion_energy(i) - data%ion_energy(j)) < e_threshold)) then

               data%point_count(i) = data%point_count(i) + merge(data%point_count(j), 0, data%ion_energy(i) <= data%ion_energy(j))

               data%point_count(i) = data%point_count(i)*merge(1, 0, data%ion_energy(i) <= data%ion_energy(j))

               data%point_count(j) = data%point_count(j) + merge(data%point_count(i), 0, data%ion_energy(i) > data%ion_energy(j))

               data%point_count(j) = data%point_count(j)*merge(1, 0, data%ion_energy(i) > data%ion_energy(j))

               if (data%ion_energy(i) > data%ion_energy(j)) exit
            end if
         end do
      end do

      results%reduced_number_points = count(data%point_count(1:data%number_points) > 0)

      write (*, *) 'number of points after removing duplicates: ', results%reduced_number_points

      allocate (results%data_clean(results%reduced_number_points, size(data%data_vec, 2)))
      allocate (results%high_dist_matrix_clean(results%reduced_number_points, results%reduced_number_points))
      allocate (results%point_count_clean(results%reduced_number_points))

      allocate (results%ion_num_clean(results%reduced_number_points))
      allocate (results%ion_energy_clean(results%reduced_number_points))
      allocate (results%ion_volume_clean(results%reduced_number_points))
      allocate (results%ion_label_clean(results%reduced_number_points))
      allocate (results%ion_symmetry_clean(results%reduced_number_points))
      allocate (results%ion_formula_clean(results%reduced_number_points))

      results%data_clean = data%data_vec(:, pack([(i, i=1, data%number_points)], data%point_count > 0))

      results%high_dist_matrix_clean = results%high_dist_matrix(pack([(i, i=1, data%number_points)], data%point_count > 0), pack([(i, i=1, data%number_points)], data%point_count > 0))

      results%point_count_clean = data%point_count(pack([(i, i=1, data%number_points)], data%point_count > 0))

      results%ion_num_clean = data%ion_num(pack([(i, i=1, data%number_points)], data%point_count > 0))
      results%ion_energy_clean = data%ion_energy(pack([(i, i=1, data%number_points)], data%point_count > 0))
      results%ion_volume_clean = data%ion_volume(pack([(i, i=1, data%number_points)], data%point_count > 0))
      results%ion_label_clean = data%ion_label(pack([(i, i=1, data%number_points)], data%point_count > 0))
      results%ion_symmetry_clean = data%ion_symmetry(pack([(i, i=1, data%number_points)], data%point_count > 0))
      results%ion_formula_clean = data%ion_formula(pack([(i, i=1, data%number_points)], data%point_count > 0))

      deallocate (valid_points)

   end subroutine remove_duplicates

end module initialisation

