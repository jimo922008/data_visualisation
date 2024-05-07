MODULE initialisation 

    USE constants
    USE data_reader

    IMPLICIT NONE

    REAL(kind=dp), allocatable :: mass_centre(:), std_dev(:)

    PUBLIC :: normalisation

contains
    subroutine normalisation(data_vec, number_points, number_features)
        implicit none
        REAL(kind=dp), intent(inout) :: data_vec(:,:)
        INTEGER, intent(in) :: number_points, number_features
        integer ::  row, col

        WRITE(stderr, *) 'initialising'
        WRITE(stderr, *) 'centring source data'

        ALLOCATE(mass_centre(number_features), std_dev(number_features))

        mass_centre = 0.0_dp
        std_dev = 0.0_dp

        !$omp parallel do reduction(+:number_points) schedule(dynamic)
        do col=1, number_points
            mass_centre = mass_centre + data_vec(1:number_features, col)
        end do
        !$omp end parallel do 

        write (stderr, *) 'normalising data'

        do col=1, number_points
            data_vec(1:number_features, col) = data_vec(1:number_features, col) - mass_centre(:)
            std_dev = std_dev + data_vec(1:number_features, col)* data_vec(1:number_features, col)
        end do 

        std_dev = sqrt(std_dev/real(number_points-1,dp))

        do row=1, number_features
            data_vec(row,:) = data_vec(row, :)/std_dev(row)
        end do
    end subroutine normalisation
end module initialisation

    