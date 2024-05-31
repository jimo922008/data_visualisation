module parameters

   implicit none

   integer, parameter           :: dp = 8

   ! high dimension parameters
   real(kind=dp), parameter     :: perplexity = 30.0_dp
   real(kind=dp), parameter     :: gr = (sqrt(5.0_dp) + 1.0_dp)/2.0_dp

   ! Low dimension parameters
   integer, parameter           :: low_dimension = 2
   real(kind=dp), parameter     :: projection_compression = 10.0_dp
   real(kind=dp), parameter     :: sphere_radius = 0.01_dp

   ! Optimisation parameters
   integer, parameter           :: max_iteration = 10000
   real(kind=dp), parameter     :: tolerance = 1e-8_dp
   integer, parameter           :: exag_cutoff = 200
   integer, parameter           :: growth_steps = 500
   real(kind=dp), parameter     :: exageration_init = 5.0_dp
   real(kind=dp), parameter     :: core_strength = 1.0_dp
   real(kind=dp), parameter     :: growth_coeff = 2.0_dp

end module parameters
