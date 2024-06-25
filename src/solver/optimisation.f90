MODULE optimisation

   USE parameters
   USE data_reader
   USE high_dimension
   USE low_dimension_probability
   USE timing_module

   IMPLICIT NONE

   PUBLIC :: adam
   PUBLIC :: tpsd
   PUBLIC :: calculate_stepsize
   PUBLIC :: initialize_variables
   PUBLIC :: loss_gradient_position
   PUBLIC :: loss_gradient_core
   PUBLIC :: gradient_vec_addnoise

   real(kind=sp) :: z, inv_z

CONTAINS

   subroutine adam(results, low_results, low_dim_params, optimisation_params)

      !> @brief Optimises the low dimensional positions using the Adam optimisation algorithm
      !> @param results structure containing the results of the high dimensional data and other parameters
      !> @param low_results structure containing the results of the low dimensional data and other parameters
      !> @param optimisation_params structure containing the parameters of the optimisation algorithm

      IMPLICIT NONE

      ! Input variables
      TYPE(high_dim_results), intent(in)             :: results
      TYPE(low_dim_results), intent(inout)           :: low_results
      TYPE(low_dim_parameters), intent(in)           :: low_dim_params
      TYPE(optimisation_parameters), intent(in)      :: optimisation_params

      ! Internal variables
      integer                      :: i, j
      real(kind=sp)                :: exaggeration, cost_criteria
      real(kind=sp)                :: gradient_norm, running_gradient_norm
      real(kind=sp), save          :: previous_gradient_norm
      real(kind=sp), allocatable   :: low_pos_vec(:)
      real(kind=sp), allocatable   :: gradient_vec(:)
      real(kind=sp), allocatable   :: step_vec(:)
      logical                      :: growth_step_limit = .true.

      ! Adam variables
      real(kind=sp)                :: beta1, beta2, beta1_t, beta2_t, epsilon, learning_rate, decay
      real(kind=sp), parameter     :: lr_factor = 0.8_sp
      real(kind=sp), allocatable   :: m(:), v(:), m_hat(:), v_hat(:)

      allocate (low_pos_vec((low_dim_params%low_dimension)*(results%reduced_number_points)))
      allocate (gradient_vec(size(low_pos_vec)))
      allocate (step_vec(size(low_pos_vec)))
      allocate (m(size(low_pos_vec)))
      allocate (v(size(low_pos_vec)))
      allocate (m_hat(size(low_pos_vec)))
      allocate (v_hat(size(low_pos_vec)))

      call start_timer()

      call initialize_variables(i, j, gradient_norm, running_gradient_norm, z, inv_z, results)

      low_pos_vec = reshape(low_results%low_dimension_position, (/(low_dim_params%low_dimension)*(results%reduced_number_points)/))

      beta1 = 0.9_sp
      beta2 = 0.999_sp
      beta1_t = 1.0_sp
      beta2_t = 1.0_sp

      epsilon = 1e-10_sp
      learning_rate = 0.01_sp
      decay = 0.00001_sp

      m = 0.0_sp
      v = 0.0_sp

      cost_criteria = (optimisation_params%threshold)*(optimisation_params%growth_coeff)

      do while ((((running_gradient_norm > log10(cost_criteria) .or. (i < 100 + (optimisation_params%exag_cutoff)))) .and. (i < (optimisation_params%maxsteps))))

         i = i + 1

         beta1_t = beta1_t*beta1

         beta2_t = beta2_t*beta2

         learning_rate = learning_rate*(1.0_sp - decay*real(i, sp))

         exaggeration = merge(1.0_sp, (optimisation_params%exaggeration_init), i > (optimisation_params%exag_cutoff))

         call loss_gradient_position(low_pos_vec, gradient_vec, exaggeration, results, low_dim_params)

         m = beta1*m + (1.0_sp - beta1)*gradient_vec

         v = beta2*v + (1.0_sp - beta2)*gradient_vec**2

         m_hat = m/(1.0_sp - beta1_t)

         v_hat = v/(1.0_sp - beta2_t)

         step_vec = learning_rate*m_hat/(sqrt(v_hat) + epsilon)

         low_pos_vec = low_pos_vec - step_vec

         gradient_norm = dot_product(step_vec, gradient_vec)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         write (*, *) ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm, learning_rate

      end do

      write (*, *) 'Growth phase'

      beta1_t = 1.0_sp
      beta2_t = 1.0_sp
      m = 0.0_sp
      v = 0.0_sp

      learning_rate = 0.01_sp

      do while (((running_gradient_norm > log10(optimisation_params%threshold)) .or. (growth_step_limit)) .and. (i + j <optimisation_params%maxsteps))

         j = j + 1

         beta1_t = beta1_t*beta1

         beta2_t = beta2_t*beta2

         learning_rate = learning_rate*(1.0_sp - decay*real(i, sp))

         call loss_gradient_core(low_pos_vec, gradient_vec, results, low_dim_params, low_results, optimisation_params)

         m = beta1*m + (1.0_sp - beta1)*gradient_vec

         v = beta2*v + (1.0_sp - beta2)*gradient_vec**2

         m_hat = m/(1.0_sp - beta1_t)

         v_hat = v/(1.0_sp - beta2_t)

         step_vec = learning_rate*m_hat/(sqrt(v_hat) + epsilon)

         low_pos_vec = low_pos_vec - step_vec

         gradient_norm = dot_product(step_vec, gradient_vec)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/100

         call handle_growth_phase(j, low_results%point_radius, optimisation_params, growth_step_limit)

         write (*, *) 'radius', sum(low_results%point_radius), ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm

      end do

      low_results%low_dimension_position = reshape(low_pos_vec, (/low_dim_params%low_dimension, results%reduced_number_points/))

      call stop_timer()
      write (*, *) 'Time: ', elapsed_time()

   end subroutine adam

   subroutine tpsd(data, results, low_results, low_dim_params, optimisation_params)

      !> @brief Optimises the low dimensional positions using the TPSD optimisation algorithm
      !> @param data structure containing the data
      !> @param results structure containing the results of the high dimensional data
      !> @param low_results structure containing the results of the low dimensional data
      !> @param low_dim_params structure containing the parameters of the low dimensional data
      !> @param optimisation_params structure containing the parameters of the optimisation algorithm

      IMPLICIT NONE

      ! Input variables
      TYPE(file_data), intent(in)                    :: data
      TYPE(high_dim_results), intent(in)             :: results
      TYPE(low_dim_results), intent(inout)           :: low_results
      TYPE(low_dim_parameters), intent(in)           :: low_dim_params
      TYPE(optimisation_parameters), intent(in)      :: optimisation_params

      ! Internal variables
      integer                                        :: i, j
      real(kind=sp)                                  :: exaggeration, cost_criteria
      logical                                        :: growth_step_limit = .true.
      real(kind=sp)                                  :: step_size, gradient_norm, running_gradient_norm
      real(kind=sp), allocatable                     :: low_pos_vec(:), perturbed_pos_vec(:)
      real(kind=sp), allocatable                     :: gradient_vec(:), gradient_vec_noise(:)
      real(kind=sp), allocatable                     :: gradient_vec_current(:), gradient_vec_perturbed(:), refined_gradient_vec(:)
      character(len=128)                             :: filename
      real(kind=sp), parameter                       :: delta = 1e-2_sp
      real(kind=sp), parameter                       :: delta_g = 1.0_sp

      allocate (low_pos_vec((low_dim_params%low_dimension)*(results%reduced_number_points)))
      allocate (gradient_vec(size(low_pos_vec)))
      allocate (gradient_vec_noise(size(low_pos_vec)))
      allocate (perturbed_pos_vec(size(low_pos_vec)))
      allocate (gradient_vec_current(size(low_pos_vec)))
      allocate (refined_gradient_vec(size(low_pos_vec)))
      allocate (gradient_vec_perturbed(size(low_pos_vec)))

      call start_timer()

      call initialize_variables(i, j, gradient_norm, running_gradient_norm, z, inv_z, results)

      low_pos_vec = reshape(low_results%low_dimension_position, (/(low_dim_params%low_dimension)*(results%reduced_number_points)/))

      cost_criteria = (optimisation_params%threshold)*(optimisation_params%growth_coeff)

      do while ((((running_gradient_norm > log10(cost_criteria) .or. (i < 100 + (optimisation_params%exag_cutoff)))) .and. (i < (optimisation_params%maxsteps))))

         i = i + 1

         exaggeration = merge(1.0_sp, (optimisation_params%exaggeration_init), i > (optimisation_params%exag_cutoff))

         perturbed_pos_vec = low_pos_vec + delta*gradient_vec_current

         call loss_gradient_position(perturbed_pos_vec, gradient_vec_perturbed, exaggeration, results, low_dim_params)

         refined_gradient_vec = (gradient_vec_perturbed - gradient_vec_current)/delta

         call gradient_vec_addnoise(refined_gradient_vec, gradient_vec_noise, 1e-2_sp)

  call calculate_stepsize(low_pos_vec, gradient_vec_noise, step_size, init=((i == 1) .or. (i == (optimisation_params%exag_cutoff))))

         low_pos_vec = low_pos_vec - step_size*refined_gradient_vec

         gradient_norm = dot_product(step_size*refined_gradient_vec, gradient_vec_noise)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         write (*, *) ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm, ' Step size: ', step_size

      end do

      write (*, *) 'Growth phase'

      do while (((running_gradient_norm > log10(optimisation_params%threshold)) .or. (growth_step_limit)) .and. (i + j < optimisation_params%maxsteps))
         j = j + 1

         perturbed_pos_vec = low_pos_vec + delta_g*gradient_vec_current

       call loss_gradient_core(perturbed_pos_vec, gradient_vec_perturbed, results, low_dim_params, low_results, optimisation_params)

         refined_gradient_vec = (gradient_vec_perturbed - gradient_vec_current)/delta_g

         call gradient_vec_addnoise(refined_gradient_vec, gradient_vec_noise, 1e-2_sp)

         call calculate_stepsize(low_pos_vec, gradient_vec_noise, step_size, init=(j == 1))

         low_pos_vec = low_pos_vec - step_size*refined_gradient_vec

         gradient_norm = dot_product(step_size*refined_gradient_vec, refined_gradient_vec)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/100

         call handle_growth_phase(j, low_results%point_radius, optimisation_params, growth_step_limit)

         write (*, *) ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm, ' Step size: ', step_size
         write (*, *) 'point radius', sum(low_results%point_radius)

      end do

      low_results%low_dimension_position = reshape(low_pos_vec, (/low_dim_params%low_dimension, results%reduced_number_points/))

      call stop_timer()
      write (*, *) 'Time: ', elapsed_time()
   end subroutine tpsd

   subroutine initialize_variables(i, j, gradient_norm, running_gradient_norm, z, inv_z, results)
      !> @brief Initialises variables for the optimisation process
      !> @param i The current iteration
      !> @param j The current growth phase iteration
      !> @param gradient_norm The norm of the gradient
      !> @param running_gradient_norm The running average of the gradient norm
      !> @param results structure containing the results of the high dimensional data

      IMPLICIT NONE

      ! Input variables
      TYPE(high_dim_results), intent(in) :: results

      ! Internal variables
      integer, intent(out) :: i, j
      real(kind=sp), intent(out) :: gradient_norm, running_gradient_norm
      real(kind=sp), intent(inout) :: z, inv_z

      i = 0
      j = 0
      gradient_norm = huge(1.0_sp)
      running_gradient_norm = 0.0_sp
      z = real(results%reduced_number_points, sp)*(real(results%reduced_number_points, sp) - 1.0_sp)
      inv_z = 1.0_sp/z
   end subroutine initialize_variables

   subroutine loss_gradient_position(low_pos_vec, gradient_vec, exaggeration, results, low_dim_params)

      !> @brief Calculates the gradient of the loss function with respect to the low dimensional positions
      !> @param low_pos_vec The low dimensional positions
      !> @param gradient_vec The gradient of the loss function
      !> @param exaggeration The exaggeration factor
      !> @param results structure containing the results of the high dimensional data
      !> @param low_dim_params structure containing the parameters of the low dimensional data

      IMPLICIT NONE

      ! Input variables
      real(kind=sp), intent(in)                 :: exaggeration
      TYPE(high_dim_results), intent(in)        :: results
      TYPE(low_dim_parameters), intent(in)      :: low_dim_params

      ! Ouput variables
      real(kind=sp), intent(inout)              :: gradient_vec(:)
      real(kind=sp), intent(inout)              :: low_pos_vec(:)

      ! Internal variables
      real(kind=sp), allocatable                :: vec(:), pos(:)
      real(kind=sp)                             :: qij, rij2
      integer                                   :: i, j, index_i, index_ii, index_j, index_jj, index

      allocate (pos(low_dim_params%low_dimension))
      allocate (vec(low_dim_params%low_dimension))

      gradient_vec = 0.0_sp
      index = results%reduced_number_points

      !$omp parallel do private(pos, rij2, qij, vec, index_i, index_ii, index_j, index_jj) reduction(+:gradient_vec) schedule(dynamic)
      do i = 1, index
         index_i = (i - 1)*low_dim_params%low_dimension + 1
         index_ii = i*low_dim_params%low_dimension
         do j = i + 1, index
            index_j = (j - 1)*low_dim_params%low_dimension + 1
            index_jj = j*low_dim_params%low_dimension
            pos(:) = low_pos_vec(index_i:index_ii) - low_pos_vec(index_j:index_jj)
            rij2 = dot_product(pos, pos)
            qij = 1.0_sp/(1.0_sp + rij2)*inv_z
            vec(:) = 4.0_sp*z*(exaggeration*(results%pij(j, i)) - (1 - (results%pij(j, i)))/(1 - qij)*qij)*qij*pos(:)
            gradient_vec(index_i:index_ii) = gradient_vec(index_i:index_ii) + vec(:)
            gradient_vec(index_j:index_jj) = gradient_vec(index_j:index_jj) - vec(:)
         end do
      end do
      !$omp end parallel do

      deallocate (pos)
      deallocate (vec)

   end subroutine loss_gradient_position

   subroutine loss_gradient_core(low_pos_vec, gradient_vec, results, low_dim_params, low_results, optimisation_params)
      !> @brief Calculates the gradient of the loss function with respect to the low dimensional positions
      !> @param low_pos_vec The low dimensional positions
      !> @param gradient_vec The gradient of the loss function

      IMPLICIT NONE

      ! Input variables
      TYPE(high_dim_results), intent(in)        :: results
      TYPE(low_dim_parameters), intent(in)      :: low_dim_params
      TYPE(low_dim_results), intent(in)         :: low_results
      TYPE(optimisation_parameters), intent(in) :: optimisation_params

      ! Output variables
      real(kind=sp), intent(inout)  :: low_pos_vec(:)
      real(kind=sp), intent(inout)  :: gradient_vec(:)

      ! Internal variables
      real(kind=sp), allocatable    :: vec(:), pos(:)
      real(kind=sp)                 :: rij2, qij, dist
      integer                       :: i, j, index_i, index_ii, index_j, index_jj, index

      allocate (pos(low_dim_params%low_dimension))
      allocate (vec(low_dim_params%low_dimension))

      gradient_vec = 0.0_sp
      index = results%reduced_number_points

      !$omp parallel do private(pos, rij2, qij, dist, vec, index_i, index_ii, index_j, index_jj) reduction(+:gradient_vec) schedule(dynamic)
      do i = 1, index
         index_i = (i - 1)*low_dim_params%low_dimension + 1
         index_ii = i*low_dim_params%low_dimension
         do j = i + 1, index
            index_j = (j - 1)*low_dim_params%low_dimension + 1
            index_jj = j*low_dim_params%low_dimension
            pos(:) = low_pos_vec(index_i:index_ii) - low_pos_vec(index_j:index_jj)
            rij2 = dot_product(pos, pos)
            qij = 1.0_sp/(1.0_sp + rij2)*inv_z
            vec(:) = 4.0_sp*z*((results%pij(j, i)) - (1 - (results%pij(j, i)))/(1 - qij)*qij)*qij*pos(:)
            gradient_vec(index_i:index_ii) = gradient_vec(index_i:index_ii) + vec(:)
            gradient_vec(index_j:index_jj) = gradient_vec(index_j:index_jj) - vec(:)

            dist = sqrt(rij2)
            if (dist < low_results%point_radius(i) + low_results%point_radius(j)) then
               vec(:) = -pos/dist
               dist = (low_results%point_radius(i) + low_results%point_radius(j) - dist)/2.0_sp
            gradient_vec(index_i:index_ii) = gradient_vec(index_i:index_ii) + vec(:)*dist*(optimisation_params%core_strength)/2.0_sp
            gradient_vec(index_j:index_jj) = gradient_vec(index_j:index_jj) - vec(:)*dist*(optimisation_params%core_strength)/2.0_sp
            end if
         end do
      end do
      !$omp end parallel do

   end subroutine loss_gradient_core

   subroutine calculate_stepsize(current_position, current_gradient, step_size, init)

      ! Input variables
      logical, optional, intent(in)                                    :: init
      real(kind=sp), dimension(:), intent(in)                          :: current_position, current_gradient
      real(kind=sp), intent(out)                                       :: step_size

      ! Internal variables
      real(kind=sp), save, allocatable, dimension(:)                   :: previous_position, previous_gradient
      real(kind=sp), dimension(size(current_position))                 :: position_change, gradient_change
      real(kind=sp)                                                    :: gradient_change_magnitude, position_gradient_dot_product
      real(kind=sp), parameter                                         :: default_step_size = 1e-1_sp

      if (init) then
         if (allocated(previous_position)) deallocate (previous_position)
         if (allocated(previous_gradient)) deallocate (previous_gradient)
         allocate (previous_position(size(current_position)))
         allocate (previous_gradient(size(current_gradient)))
         step_size = default_step_size

      else

         position_change = current_position - previous_position
         gradient_change = current_gradient - previous_gradient

         gradient_change_magnitude = dot_product(gradient_change, gradient_change)
         position_gradient_dot_product = dot_product(position_change, gradient_change)

         step_size = abs(position_gradient_dot_product/gradient_change_magnitude)

         step_size = max(step_size, default_step_size)

      end if

      previous_position = current_position
      previous_gradient = current_gradient

   end subroutine calculate_stepsize

   subroutine gradient_vec_addnoise(gradient_vec, gradient_vec_noise, r)

      !> @brief Adds noise to the gradient vector
      !> @param gradient_vec The gradient vector
      !> @param gradient_vec_noise The gradient vector with noise
      !> @param r The noise factor, i.e. the magnitude of the noise
      !> @param results structure containing the results of the high dimensional data
      !> @param low_dim_params structure containing the parameters of the low dimensional data

      IMPLICIT NONE

      ! Input variables
      real(kind=sp), intent(in)                  :: r
      real(kind=sp), intent(in)                  :: gradient_vec(:)

      ! Output variables
      real(kind=sp), allocatable, intent(out)    :: gradient_vec_noise(:)

      ! Internal variables
      real(kind=sp), allocatable                 :: noise_vector(:)
      real(kind=sp)                              :: d

      allocate (gradient_vec_noise(size(gradient_vec)))
      allocate (noise_vector(size(gradient_vec)))
      gradient_vec_noise = 0.0_sp

      noise_vector = 0.0_sp

      call random_number(d)

      d = r*d**(1/real(size(gradient_vec), sp))

      call random_add_noise(noise_vector, 1.0_sp)

      noise_vector = d*noise_vector/norm2(noise_vector)

      gradient_vec_noise = gradient_vec*(1 + noise_vector)

   end subroutine gradient_vec_addnoise

   subroutine handle_growth_phase(j, point_radius, optimisation_params, growth_step_limit)
      implicit none

      ! Input variables
      TYPE(optimisation_parameters), intent(in) :: optimisation_params
      integer, intent(in)                       :: j
      real(kind=sp), intent(inout)              :: point_radius(:)
      logical, intent(inout)                    :: growth_step_limit

      if (j < optimisation_params%growth_steps) then
         if (j < 2) then
            point_radius = point_radius*(real(j))/real(optimisation_params%growth_steps)
         else
            point_radius = point_radius*real(j)/((real(j) - 1.0_sp))
         end if
      end if

      if (j > (optimisation_params%growth_steps + 100)) growth_step_limit = .false.

   end subroutine handle_growth_phase

end MODULE optimisation

