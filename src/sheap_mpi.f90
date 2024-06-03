program sheap_mpi

   !USE mpi
   USE omp_lib
   USE data_reader
   USE data_writer
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   USE optimisation
   !USE mpi_model

   character(len=256) :: filename
   integer             :: ierr, rank, nranks

   !call omp_set_num_threads(8)

   !call MPI_Init(ierr)
   !call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   !call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

   if (rank == 0) then

      call get_command_argument(1, filename)
      if (trim(filename) == '') then
         print *, 'No filename provided.'
         stop
      end if

      call read_file(trim(filename))

      call normalisation(data_vec, number_points, number_features)

      print *, 'Data normalised.'

      call high_dimension_distribution()
      call low_dimension_distribution()
      write (*, *) 'Optimising the low dimension distribution'

   end if

   !call tpsd_mpi(1e-8_sp, 10000, rank, nranks)

   if (rank == 0) then
      call write_file(trim('LJ13-sheap.xyz'))
   end if

   !call MPI_Finalize(ierr)

end program sheap_mpi
