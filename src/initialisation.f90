module initialisation 

contains
    subroutine centre_of_mass()
        integer ::  ni,nj

        write (stderr, *) 'initialising'

        write (stderr, *) 'centring source data'

        hcom = 0.0_dp
        !$omp parallel do reduction(+:hcom) schedule(dynamic)
        do ni=1, npoints
            hcom = hcom + data_vec(1:hdim, ni)
        end do
        !$omp end parallel do 

    subroutine normalisation()

        real(kind=dp) :: stddev

        write (stderr, *) 'normalising data'

        stddev = 0.0_dp
        do ni=1, npoints
            data_vec(1:hdim, ni) = data_vec(1:hdim, ni) - hcom(:)
            stddev = stddev + data_vec(1:hdim, ni)* data_vec(1:hdim, ni)
        end do 

        stddev = sqrt(stddev/ral(npoints-1,dp))

        do ni=1, hdim
            data_vec(ni,:) = data_vec(ni, :)/stddev(ni)
        end do


    subroutine high_dimension_distance()

        write (stderr, *) 'calculating Pi|j'
        do ni=1, npoints
            pij_conditional= dot_product(data_vec(1:hdim, ni), data_vec(1:hdim, ni))
        end do 

        !$omp parallel do private(nj) schedule(dynamic)
        do ni=1, npoints
            do nj=1, npoints, 3
                pij(ni,nj) = pij_conditional(ni) + pij_conditional(nj)
                pij(ni,nj+1) = pij_conditional(ni) + pij_conditional(nj+1)
                pij(ni,nj+2) = pij_conditional(ni) + pij_conditional(nj+2)
            end do
        end do

    