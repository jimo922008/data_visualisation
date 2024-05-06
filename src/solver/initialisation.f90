MODULE initialisation 

    USE constants
    USE data_io

    IMPLICIT NONE

    REAL(), allocatable :: mass_centre, std_dev

    PUBLIC :: centre_of_mass

    
contains
    subroutine normalisation()

        integer ::  row, col
        real(kind=dp) :: stddev


        WRITE(stderr, *) 'initialising'

        WRITE(stderr, *) 'centring source data'

        ALLOCATE(
            mass_centre (col_count)
            std_dev (col_count)
        )

        mass_centre = 0.0_dp

        !$omp parallel do reduction(+:row_count) schedule(dynamic)
        do row=1, row_count
            mass_centre = mass_centre + data_vec(1:col_count, row)
        end do
        !$omp end parallel do 

        write (stderr, *) 'normalising data'

        stddev = 0.0_dp
        do row=1, row_count
            data_vec(1:row_count, ni) = data_vec(1:row_count, row) - centre_of_mass(:)
            stddev = stddev + data_vec(1:row_count, row)* data_vec(1:row_count, row)
        end do 

        stddev = sqrt(stddev/ral(npoints-1,dp))

        do ni=1, hdim
            data_vec(ni,:) = data_vec(ni, :)/stddev(ni)
        end do


    