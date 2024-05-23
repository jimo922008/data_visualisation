MODULE initialisation 

    USE constants
    USE data_reader

    IMPLICIT NONE

    PUBLIC :: normalisation
    PUBLIC :: remove_duplicates

    private
    REAL(kind=dp), dimension(:), allocatable :: mass_centre, std_dev

    INTEGER, PUBLIC :: reduced_number_points

contains
    subroutine normalisation(data_vec, number_points, number_features)
        implicit none
        REAL(kind=dp), intent(inout) :: data_vec(:,:)
        INTEGER, intent(in) :: number_points, number_features
        integer ::  i, j

        WRITE(stderr, *) 'centering source data'

        ALLOCATE(mass_centre(number_features))
        ALLOCATE(std_dev(number_features))

        mass_centre = 0.0_dp
        std_dev = 0.0_dp

        !$omp parallel do reduction(+:mass_centre) schedule(dynamic)
        do i=1, number_points
            mass_centre = mass_centre + data_vec(1:number_features, i)
        end do
        !$omp end parallel do 

        mass_centre = mass_centre/real(number_points, dp)

        write (stderr, *) 'normalising data'

        do i=1, number_points
            data_vec(1:number_features, i) = data_vec(1:number_features, i) - mass_centre(:)
            !std_dev = std_dev + data_vec(1:number_features, col)**2
        end do 

        !std_dev = sqrt(std_dev/real(number_points-1,dp))

        !do row=1, number_features
        !    if(std_dev(row)== 0.0_dp) cycle
        !    data_vec(row,:) = data_vec(row, :)/std_dev(row)
        !end do

    end subroutine normalisation

    subroutine remove_duplicates(data_vec, high_dist_matrix, number_points, similar_threshold_, energy_threshold)
        implicit none
        REAL(kind=dp), intent(in) :: data_vec(:,:)
        REAL(kind=dp), intent(in) :: similar_threshold_, energy_threshold
        REAL(kind=dp), intent(in) :: high_dist_matrix(:,:)
        REAL(kind=dp)             :: similar_threshold
        INTEGER, intent(inout)    :: number_points
        integer :: i, j

        similar_threshold=similar_threshold_*sum(std_dev)/100_dp

        do i = 1, number_points
            if (point_count(i) == 0) cycle
            do j= i+1, number_points                
                if (point_count(j) ==0) cycle                  
                if ((high_dist_matrix(j,i) < similar_threshold**2) .and. (abs(ion_energy(i)-ion_energy(j)) < energy_threshold)) then
                    if (ion_energy(i) <= ion_energy(j)) then
                        point_count(i) = point_count(i) + point_count(j)
                        point_count(j) = 0
                    else
                        point_count(j) = point_count(j) + point_count(i)
                        point_count(i) = 0
                        exit
                    end if
                end if
            end do
        end do

        reduced_number_points = count(point_count(1:number_points) > 0)
        write (stderr, *) 'number of points after removing duplicates: ', reduced_number_points 

    end subroutine remove_duplicates

end module initialisation

    