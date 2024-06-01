module parameters

   implicit none

   integer, parameter           :: dp = 8
   integer, parameter           :: sp = 4

   ! high dimension parameters
   real(kind=sp), parameter     :: perplexity = 30.0_sp
   real(kind=sp), parameter     :: gr = (sqrt(5.0_sp) + 1.0_sp)/2.0_sp

   ! Low dimension parameters
   integer, parameter           :: low_dimension = 2
   real(kind=sp), parameter     :: projection_compression = 10.0_sp
   real(kind=sp), parameter     :: sphere_radius = 0.01_sp

   ! Optimisation parameters
   integer, parameter           :: max_iteration = 10000
   real(kind=sp), parameter     :: tolerance = 1e-8_sp
   integer, parameter           :: exag_cutoff = 200
   integer, parameter           :: growth_steps = 500
   real(kind=sp), parameter     :: exaggeration_init = 5.0_sp
   real(kind=sp), parameter     :: core_strength = 1.0_sp
   real(kind=sp), parameter     :: growth_coeff = 2.0_sp

end module parameters
