program main

    USE omp_lib
    USE data_reader
    USE data_writer
    USE initialisation
    USE high_dimension
    USE low_dimension_probability
    USE optimisation
        
    integer :: iostat
    character(len=256) :: filename

    ! Initialize the filename

    call omp_set_num_threads(4)

    call get_command_argument(1, filename)
    if (trim(filename) == '') then
        print *, 'No filename provided.'
        stop
    end if

    print *, 'Reading data from file:', trim(filename)

    ! Read the data from the file
    call read_file(trim(filename))

    if (allocated(data_vec)) then
        print *, 'Data read from file:'
     
        if (size(data_vec, 1) > 0 ) then
            print *, number_points, 'points read.'
        else
            print *, 'Data array is empty or dimensions are zero.'
        end if

    else
        print *, 'No data read or allocation failed.'
    end if

    call normalisation(data_vec, number_points, number_features)
    
    print *, 'Data normalised.'
    print *, 'Insert the perplexity value: '
    read *, perplexity

    call high_dimension_distribution()
    write (*,*) sum(data_vec(:,1)), sum(data_vec(:,2)), sum(data_vec(:,3)), sum(data_vec(:,4)), sum(data_vec(:,5))
    write (*,*) high_dist_matrix(1,1), high_dist_matrix(1,2), high_dist_matrix(1,3), high_dist_matrix(1,4), high_dist_matrix(1,5)
    write (*,*) sigma(1), sigma(2), sigma(3), sigma(4), sigma(5)
    write (*,*) pij(1,1), pij(1,2), pij(1,3), pij(1,4), pij(1,5)
    call low_dimension_distribution()

    write (*,*) 'Optimising the low dimension distribution'

    call tpsd (1e-8_dp, 10000, 1.0_dp)


end program main
