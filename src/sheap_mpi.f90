program sheap_mpi

   USE mpi
   USE parameters
   USE omp_lib
   USE data_reader
   USE data_writer
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   USE optimisation
   USE mpi_model

   character(len=256)  :: filename
   integer             :: ierr, rank, nranks
   real(kind=sp), dimension(:), allocatable ::pij_1d, low_pos_vec, point_radius_1d

   call omp_set_num_threads(8)
   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)

   if (rank == 0) then

      call get_command_argument(1, filename)
      if (trim(filename) == '') then
         print *, 'No filename provided.'
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

      call read_file(trim(filename))

      call normalisation(data_vec, number_points, number_features)

      print *, 'Data normalised.'

      call high_dimension_distribution()
      call low_dimension_distribution()
      write (*, *) 'Optimising the low dimension distribution'

      allocate (pij_1d(reduced_number_points*reduced_number_points))
      allocate (low_pos_vec(low_dimension*reduced_number_points))

      pij_1d = reshape(pij, shape(pij_1d))
      low_pos_vec = reshape(low_dimension_position, (/low_dimension*reduced_number_points/))
      point_radius_1d = point_radius
   end if

   call MPI_Bcast(reduced_number_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (rank /= 0) then
      allocate (pij_1d(reduced_number_points*reduced_number_points))
      allocate (low_pos_vec(low_dimension*reduced_number_points))
      allocate (point_radius_1d(reduced_number_points))
   end if

   call MPI_Bcast(pij_1d, reduced_number_points*reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(point_radius_1d, reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   call tpsd_mpi(pij_1d, point_radius_1d, low_pos_vec, tolerance, max_iteration, rank, nranks)

   if (rank == 0) then
      call write_file(trim('LJ13-sheap.xyz'))
   end if

   call MPI_Finalize(ierr)

end program sheap_mpi
