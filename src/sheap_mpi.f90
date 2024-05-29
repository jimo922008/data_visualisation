program sheap_mpi

include 'mpif.h'

    USE omp_lib
    USE data_reader
    USE data_writer
    USE initialisation
    USE high_dimension
    USE low_dimension_probability
    USE optimisation
    USE mpi_model 
        
    character(len=256) :: filename
    integer ierr, tag

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

    ! Initialize the filename

    call omp_set_num_threads(8)

    if (rank == 0) then
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
        write (*,*) 'Optimising the low dimension distribution'

        

    else 
        call tpsd_mpi(threshold, maxsteps, nranks, ierr)

    end if

    call MPI_Finalize(ierr)

    write(*,*) 'final cost = ', cost

    if (rank == 0) then
        call write_file(trim('LJ13-sheap.xyz'))
    end if

end program sheap_mpi