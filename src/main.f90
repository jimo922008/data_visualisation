!> @brief Main program for the TPSD code
!> @details Stochastic Hyperspace Embedding And Projection (SHEAP) is a dimensionality reduction method designed for visualising potential energy surfaces.
!! Computational structure prediction can assist the discovery of new materials. One searches for the most stable configurations of a given set of atomic building blocks, 
!! which correspond to the deepest regions of an energy landscape—the system's energy as a function of the relative positions of its atoms. 
!! To explore these landscapes efficiently, it is important to understand their topologies. However, they exist in spaces with very large numbers of dimensions, 
!! making them difficult to visualise. SHEAP uses dimensionality reduction through manifold learning to effectively visualise the distribution of stable structures across a high-dimensional energy landscape.
!! This program reads in a file containing the pairwise distance of a set of atoms and their corresponding properties.
!! It then uses SHEAP to reduce the dimensionality of the data and writes the results to a file.



program main

   USE omp_lib
   USE parameters
   USE data_reader
   USE data_writer
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   USE optimisation

   implicit none

   ! Declare variables
   character(len=256)            :: filename
   type(file_data)               :: data
   type(high_dim_parameters)     :: high_dim_params
   type(low_dim_parameters)      :: low_dim_params
   type(high_dim_results)        :: results
   type(low_dim_results)         :: low_results
   type(optimisation_parameters) :: optimisation_params

   ! Initialize the filename

   call omp_set_num_threads(8)

   call get_command_argument(1, filename)
   if (trim(filename) == '') then
      print *, 'No filename provided.'
      stop
   end if

   print *, 'Reading data from file:', trim(filename)

   ! Read the data from the file
   call read_file(trim(filename), data)

   call normalisation(data)

   print *, 'Data normalised.'

   call high_dimension_distribution(data, results, high_dim_params)

   call low_dimension_distribution(data, high_dim_params, low_dim_params, results, low_results)

   call tpsd(data, results, low_results, low_dim_params, optimisation_params)
   !call adam(results, low_results, low_dim_params, optimisation_params)
   call write_file(trim('LJ13-sheap.xyz'), data, results, low_results, low_dim_params)

end program main
