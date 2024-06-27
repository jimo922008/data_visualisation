module parameters

   implicit none

   integer, parameter                 :: dp = 8
   integer, parameter                 :: sp = 4

   type :: high_dim_parameters
      real(kind=sp)                   :: similar_threshold = 0.5_sp
      real(kind=sp)                   :: energy_threshold = 0.10000E+09_sp
      real(kind=sp)                   :: perplexity = 30.0_sp
      real(kind=sp)                   :: gr = (sqrt(5.0_sp) + 1.0_sp)/2.0_sp
      real(kind=sp)                   :: find_sigma_tolerance = 1.0e-13_sp
   end type high_dim_parameters

   type :: low_dim_parameters
      integer                         :: low_dimension = 3
      real(kind=sp)                   :: projection_compression = 10.0_sp
      real(kind=sp)                   :: sphere_radius = 0.01_sp
   end type low_dim_parameters

   type :: optimisation_parameters
      integer                         :: maxsteps = 10000
      real(kind=sp)                   :: threshold = 1e-8_sp
      integer                         :: exag_cutoff = 200
      integer                         :: growth_steps = 500
      real(kind=sp)                   :: exaggeration_init = 5.0_sp
      real(kind=sp)                   :: core_strength = 1.0_sp
      real(kind=sp)                   :: growth_coeff = 2.0_sp
      INTEGER                         :: block_size = 16
      real(kind=sp)                   :: momentum_coeff = 0.9_sp
   end type optimisation_parameters

   type :: file_data
      integer                         :: number_points
      integer                         :: number_features
      integer, allocatable            :: point_count(:)
      real(kind=sp), allocatable      :: data_vec(:, :)
      real(kind=sp), allocatable      :: data_mean(:)
      real(kind=sp), allocatable      :: data_std(:)

      integer, allocatable            :: ion_num(:)
      real(kind=sp), allocatable      :: ion_energy(:)
      real(kind=sp), allocatable      :: ion_volume(:)
      character(len=256), allocatable :: ion_label(:)
      character(len=256), allocatable :: ion_symmetry(:)
      character(len=256), allocatable :: ion_formula(:)

   end type file_data

   type :: high_dim_results
      integer                         :: reduced_number_points
      REAL(kind=sp), allocatable      :: high_dist_matrix(:, :)
      REAL(kind=sp), allocatable      :: high_dist_matrix_clean(:, :)
      INTEGER, allocatable            :: point_count_clean(:)
      REAL(kind=sp), allocatable      :: data_clean(:, :)
      REAL(kind=sp), allocatable      :: sigma(:)
      REAL(kind=sp), allocatable      :: pij(:, :)

      integer, allocatable            :: ion_num_clean(:)
      real(kind=sp), allocatable      :: ion_energy_clean(:)
      real(kind=sp), allocatable      :: ion_volume_clean(:)
      character(len=256), allocatable :: ion_label_clean(:)
      character(len=256), allocatable :: ion_symmetry_clean(:)
      character(len=256), allocatable :: ion_formula_clean(:)
   end type high_dim_results

   type :: low_dim_results
      real(kind=sp), allocatable      :: low_dimension_position(:, :)
      real(kind=sp), allocatable      :: low_dimension_distance(:, :)
      real(kind=sp), allocatable      :: low_dimension_mean(:)
      real(kind=sp), allocatable      :: low_dimension_std(:)
      real(kind=sp), allocatable      :: qij(:, :)
      real(kind=sp), allocatable      :: point_radius(:)
   end type low_dim_results

end module parameters
