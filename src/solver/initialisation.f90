MODULE initialisation 

    USE constants
    USE data_reader

    IMPLICIT NONE

    PUBLIC :: normalisation
    PUBLIC :: remove_duplicates

    private
    REAL(kind=dp), dimension(:), allocatable  :: mass_centre, std_dev
    real(dp), allocatable, dimension(:,:), PUBLIC  :: data_clean(:,:)
    real(dp), allocatable, dimension(:,:), PUBLIC  :: high_dist_matrix_clean(:,:)
    real(dp), allocatable, dimension(:), PUBLIC    :: point_count_clean(:)

    INTEGER, PUBLIC :: reduced_number_points

contains
    subroutine normalisation(data_vec, number_points, number_features)
        
        REAL(kind=dp), intent(inout) :: data_vec(:,:)
        INTEGER, intent(in) :: number_points, number_features

        WRITE(stderr, *) 'centering source data'

        allocate(mass_centre(number_features))

        allocate(std_dev(number_features))

        mass_centre = sum(data_vec, dim=2)/real(number_points, dp)

        data_vec = data_vec - spread(mass_centre, 2, number_points)

        std_dev = sqrt(sum(data_vec**2, dim=2)/real(number_points-1,dp))

        write(stderr, *) 'normalising source data'

    end subroutine normalisation

    subroutine remove_duplicates(data_vec, high_dist_matrix, number_points, similar_threshold_, energy_threshold)
        
        REAL(kind=dp), intent(in) :: data_vec(:,:)
        REAL(kind=dp), intent(in) :: similar_threshold_, energy_threshold
        REAL(kind=dp), intent(in) :: high_dist_matrix(:,:)
        REAL(kind=dp)             :: similar_threshold
        INTEGER, dimension(:), allocatable     :: valid_points
        INTEGER, intent(inout)    :: number_points
        integer :: i, j, l, m

        similar_threshold=similar_threshold_*sum(std_dev)/100_dp

        allocate (valid_points(count(point_count /= 0)))
        
        valid_points = pack([(i, i=1, number_points)], point_count /= 0)

        !$omp parallel do shared(ion_energy, high_dist_matrix) private(i, j)
        do l = 1, size(valid_points)
            i = valid_points(l)
            do m= l+1, size(valid_points)               
                j = valid_points(m)                 
                if ((high_dist_matrix(j,i) < similar_threshold**2) .and. (abs(ion_energy(i)-ion_energy(j)) < energy_threshold)) then
                    
                    point_count(i) = point_count(i) + merge(point_count(j), 0, ion_energy(i) <= ion_energy(j))
                    point_count(i) = point_count(i) * merge(1, 0, ion_energy(i) <= ion_energy(j))
                    point_count(j) = point_count(j) + merge(point_count(i), 0, ion_energy(i) > ion_energy(j))
                    point_count(j) = point_count(j) * merge(1, 0, ion_energy(i) > ion_energy(j))
                        
                    if (ion_energy(i) > ion_energy(j)) exit
                end if
            end do
        end do
        !$omp end parallel do

        reduced_number_points = count(point_count(1:number_points) > 0)

        write (stderr, *) 'number of points after removing duplicates: ', reduced_number_points

        allocate(data_clean(reduced_number_points, size(data_vec,2)))
        allocate(high_dist_matrix_clean(reduced_number_points, reduced_number_points))
        allocate(point_count_clean(reduced_number_points))

        data_clean = data_vec(:, pack([(i, i=1, number_points)], point_count > 0))

        high_dist_matrix_clean = high_dist_matrix(pack([(i, i=1, number_points)], point_count > 0), pack([(i, i=1, number_points)], point_count > 0))

        point_count_clean = point_count(pack([(i, i=1, number_points)], point_count>0))

        deallocate(valid_points)

    end subroutine remove_duplicates

end module initialisation

    