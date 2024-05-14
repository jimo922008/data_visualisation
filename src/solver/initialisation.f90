MODULE initialisation 

    USE constants
    USE data_reader

    IMPLICIT NONE

    PUBLIC :: normalisation

contains
    subroutine normalisation(data_vec, number_points, number_features)
        implicit none
        REAL(kind=dp), intent(inout) :: data_vec(:,:)
        REAL(kind=dp), dimension(:), allocatable :: mass_centre, std_dev
        INTEGER, intent(in) :: number_points, number_features
        integer ::  row, col

        WRITE(stderr, *) 'initialising'
        WRITE(stderr, *) 'centering source data'

        ALLOCATE(mass_centre(number_features))
        ALLOCATE(std_dev(number_features))

        mass_centre = 0.0_dp
        std_dev = 0.0_dp

        !$omp parallel do reduction(+:number_points) schedule(dynamic)
        do col=1, number_points
            mass_centre = mass_centre + data_vec(:, col)
        end do
        !$omp end parallel do 

        mass_centre = mass_centre/real(number_points, dp)

        write (stderr, *) 'normalising data'

        do col=1, number_points
            data_vec(:, col) = data_vec(:, col) - mass_centre(:)
            std_dev = std_dev + data_vec(:, col)* data_vec(:, col)
        end do 

        std_dev = sqrt(std_dev/real(number_points,dp))

        do row=1, number_features
            if(std_dev(row).eq.0.0_dp) cycle
            data_vec(row,:) = data_vec(row, :)/std_dev(row)
        end do

    end subroutine normalisation

end module initialisation

    