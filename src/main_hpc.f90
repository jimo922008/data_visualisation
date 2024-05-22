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

    call omp_set_num_threads(8)

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

    call low_dimension_distribution()

    write (*,*) 'high_dimension_distribution', pij(1,1), pij(1,2), pij(1,3), pij(1,4), pij(1,5)
    write (*,*) 'low_dimension_position', low_dimension_position(1,1), low_dimension_position(1,2), low_dimension_position(1,3), low_dimension_position(1,4), low_dimension_position(1,5)
    write (*,*) 'hard_sphere_initialisation', point_radius(1), point_radius(2), point_radius(3), point_radius(4), point_radius(5)

    write (*,*) 'Optimising the low dimension distribution'

    call tpsd(1e-8_dp, 10000)
    call write_file(trim('lj13.xyz'))


end program main
