program main

   USE omp_lib
   USE parameters
   USE data_reader
   USE data_writer
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   USE optimisation

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

      if (size(data_vec, 1) > 0) then
         print *, number_points, 'points read.'
      else
         print *, 'Data array is empty or dimensions are zero.'
      end if

   else
      print *, 'No data read or allocation failed.'
   end if

   call normalisation(data_vec, number_points, number_features)

   print *, 'Data normalised.'

   call high_dimension_distribution()

   call low_dimension_distribution()

   !call tpsd(pij, point_radius, low_dimension_position, exaggeration_init, tolerance, max_iteration)
   call adam(pij, point_radius, low_dimension_position, tolerance, max_iteration)
   call write_file(trim('LJ13-sheap.xyz'))

end program main
