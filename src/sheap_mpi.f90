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

   ! Declare variables
   character(len=256)            :: filename
   type(file_data)               :: data
   type(high_dim_parameters)     :: high_dim_params
   type(low_dim_parameters)      :: low_dim_params
   type(high_dim_results)        :: results
   type(low_dim_results)         :: low_results
   type(optimisation_parameters) :: optimisation_params

   ! MPI variables
   integer                       :: ierr, rank, nranks
   real(kind=sp), allocatable    :: pij_1d(:), low_pos_vec(:), point_radius_1d(:)

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

      call read_file(trim(filename), data)

      call normalisation(data)

      print *, 'Data normalised.'

      call high_dimension_distribution(data, results, high_dim_params)
      call low_dimension_distribution(data, high_dim_params, low_dim_params, results, low_results)
      write (*, *) 'Optimising the low dimension distribution'

      allocate (pij_1d((results%reduced_number_points)**2))
      allocate (low_pos_vec((low_dim_params%low_dimension)*(results%reduced_number_points)))

      pij_1d = reshape(results%pij, shape(pij_1d))
      low_pos_vec = reshape(low_results%low_dimension_position, (/(low_dim_params%low_dimension)*(results%reduced_number_points)/))
      point_radius_1d = low_results%point_radius
   end if

   call MPI_Bcast(results%reduced_number_points, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

   if (rank /= 0) then
      allocate (pij_1d(results%reduced_number_points*results%reduced_number_points))
      allocate (low_pos_vec(low_dim_params%low_dimension*results%reduced_number_points))
      allocate (point_radius_1d(results%reduced_number_points))
   end if

   call MPI_Bcast(pij_1d, results%reduced_number_points*results%reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(low_pos_vec, low_dim_params%low_dimension*results%reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Bcast(point_radius_1d, results%reduced_number_points, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   call tpsd_mpi(pij_1d, point_radius_1d, low_pos_vec, results, low_results, optimisation_params, low_dim_params, rank, nranks)

   if (rank == 0) then
      low_results%point_radius = point_radius_1d
      call write_file(trim('LJ13-sheap.xyz'), data, results, low_results, low_dim_params)
   end if

   call MPI_Finalize(ierr)

end program sheap_mpi
