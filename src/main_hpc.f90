program main

    USE data_reader
    USE initialisation
    USE high_dimension
    USE low_dimension_probability
        
    integer :: iostat
    character(len=256) :: filename

    ! Initialize the filename
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

    print *, 'Program finished.'


end program main
