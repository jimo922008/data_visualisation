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
    integer ierr, tag, rank, nranks
    integer, parameter  :: master = 0 

    call omp_set_num_threads(8)

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

    if (rank == master) then 

        call get_command_argument(1, filename)
        if (trim(filename) == '') then
            print *, 'No filename provided.'
            stop
        end if

        call read_file(trim(filename))
        call initialise()
        call high_dimension_probability()
        call low_dimension_probability()

        call normalisation(data_vec, number_points, number_features)
        
        print *, 'Data normalised.'
        print *, 'Insert the perplexity value: '
        read *, perplexity

        call high_dimension_distribution()
        call low_dimension_distribution()
        write (*,*) 'Optimising the low dimension distribution'

    end if

    if (rank == master) then 
        call master_process(size, threshold, maxsteps)
    else
        call worker_process(rank, threshold, maxsteps)
    end if

    call MPI_Finalize(ierr)

    if (rank == master) then
        write(*,*) 'final cost = ', cost
        call write_file(trim('LJ13-sheap.xyz'))
    end if

end program sheap_mpi