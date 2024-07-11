!> \mainpage SHEAP
!> @brief SHEAP MPI implementation
!> @details Stochastic Hyperspace Embedding And Projection (SHEAP) is a dimensionality reduction method designed for visualising potential energy surfaces.
!! Computational structure prediction can assist the discovery of new materials. One searches for the most stable configurations of a given set of atomic building blocks, 
!! which correspond to the deepest regions of an energy landscape—the system's energy as a function of the relative positions of its atoms. 
!! To explore these landscapes efficiently, it is important to understand their topologies. However, they exist in spaces with very large numbers of dimensions, 
!! making them difficult to visualise. SHEAP uses dimensionality reduction through manifold learning to effectively visualise the distribution of stable structures across a high-dimensional energy landscape.
!! This program reads in a file containing the pairwise distance of a set of atoms and their corresponding properties.
!! It then uses SHEAP to reduce the dimensionality of the data and writes the results to a file.
!! The program is parallelised using OpenMP and MPI. The MPI implementation is used to distribute the data across multiple nodes, while OpenMP is used to parallelise the calculations on each node.
!! The program is designed to be run on a high-performance computing cluster, as well as on a PC.
!! \verbatim
!! '                            8888              '
!! '               8888888  8888    88            '
!! '              88      88          888         '
!! '             8                      8         '
!! '             88                     88        '
!! '              88 888         88   888         '
!! '         .8888888---888888888--88888888.      '
!! '       8888  888-----------------88   888     '
!! '      88     888--8888-----8888--888    88    '
!! '     88     888---88 8-----88 8--888      8   '
!! '    88      8 8---8888-----8888---88       8  '
!! '    8       8 8-------------------88       8  '
!! '    88      8 8-------------------8 88    88  '
!! '     88    88 8------------------8  88   88   '
!! '       8888   88----88-----88----8    8888    '
!! '               88-----8---8----88             '
!! '                88------------88              '
!! '                  88--------88                '
!! '                    "888888"                  '
!! '                                              '
!! '           888                                '
!! '           888                                '
!! '           888                                '
!! '  .d8888b  88888b.   .d88b.   8888b.  88888b. '
!! '  888      888 "88b d88  88b     "88b 888 "88b'
!! '  "888888. 888  888 88888888 .d888888 888  888'
!! '       888 888  888 88b.     888  888 888 d88P'
!! '   888888" 888  888  "88888  "8888888 88888P" '
!! '                                      888     '
!! '                                      888     '
!! '                                      888     '
!! '                                              '
!! '   Authors: Ben Shires (bs511@cam.ac.uk)      '
!! '            Chris Pickard (cjp20@cam.ac.uk)   '
!! '            Mo Ji (mj425@cam.ac.uk)           '
!! '                 2019-2024 (c)                '
!! '                                              '
!! \endverbatim

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
