MODULE optimisation

   USE parameters
   USE data_reader
   USE high_dimension
   USE low_dimension_probability
   USE timing_module

   IMPLICIT NONE

   PUBLIC :: tpsd
   PUBLIC :: loss_gradient
   PUBLIC :: calculate_stepsize
   PUBLIC :: initialize_variables
   PUBLIC :: loss_gradient_position
   PUBLIC :: loss_gradient_core

   real(kind=sp) :: cost_zero, final_cost, z, inv_z

contains

   subroutine tpsd(pij, point_radius, low_dimension_position, exaggeration_init, threshold, maxsteps)

      implicit none

      ! Input variables
      integer, intent(in)          :: maxsteps
      real(kind=sp), intent(in)    :: threshold
      real(kind=sp), intent(in)    :: exaggeration_init
      real(kind=sp), intent(in)    :: pij(reduced_number_points, reduced_number_points)
      real(kind=sp), intent(inout) :: point_radius(reduced_number_points)
      real(kind=sp), intent(inout) :: low_dimension_position(low_dimension, reduced_number_points)

      ! Internal variables
      integer                      :: i, j
      real(kind=sp)                :: exaggeration
      logical                      :: growth_step_limit = .true.
      real(kind=sp)                :: step_size, gradient_norm, running_gradient_norm
      real(kind=sp)                :: low_pos_vec(low_dimension*reduced_number_points)
      real(kind=sp)                :: gradient_vec(low_dimension*reduced_number_points)

      call start_timer()

      call initialize_variables(i, j, gradient_norm, running_gradient_norm)

      low_pos_vec = reshape(low_dimension_position, (/low_dimension*reduced_number_points/))

      do while ((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100 + exag_cutoff))) .and. (i < maxsteps)))

         i = i + 1

         exaggeration = merge(1.0_sp, exaggeration_init, i > exag_cutoff)

         call loss_gradient_position(low_pos_vec, gradient_vec, exaggeration)

         call calculate_stepsize(low_pos_vec, gradient_vec, step_size, init=((i == 1) .or. (i == exag_cutoff)))

         print *, 'Step size: ', step_size

         low_pos_vec = low_pos_vec - step_size*gradient_vec

         gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         write (*, *) ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm, ' Step size: ', step_size

      end do

      write (*, *) 'Growth phase'

      do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i + j < maxsteps))
         j = j + 1

         call loss_gradient_core(low_pos_vec, gradient_vec)

         call calculate_stepsize(low_pos_vec, gradient_vec, step_size, init=(j == 1))

         low_pos_vec = low_pos_vec - step_size*gradient_vec

         gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/100

         call handle_growth_phase(j, point_radius, growth_step_limit)

         write (*, *) ' Gradient norm: ', gradient_norm, 'running gradient norm', running_gradient_norm, ' Step size: ', step_size
         write (*, *) 'point radius', sum(point_radius)

      end do

      low_dimension_position = reshape(low_pos_vec, (/low_dimension, reduced_number_points/))

      call stop_timer()
      write (*, *) 'Time: ', elapsed_time()
   end subroutine tpsd

   subroutine initialize_variables(i, j, gradient_norm, running_gradient_norm)
      implicit none

      integer, intent(out) :: i, j
      real(kind=sp), intent(out) :: gradient_norm, running_gradient_norm

      i = 0
      j = 0
      gradient_norm = huge(1.0_sp)
      running_gradient_norm = 0.0_sp
      cost_zero = calculating_cost_zero(pij)
      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)
      inv_z = 1.0_sp/z
   end subroutine initialize_variables

   subroutine loss_gradient_position(low_pos_vec, gradient_vec, exaggeration)
      implicit none

      real(kind=sp), intent(in)  :: exaggeration
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)      :: gradient_vec, low_pos_vec
      real(kind=sp), dimension(low_dimension) :: vec, pos
      real(kind=sp) :: qij, rij2
      integer :: i, j

      gradient_vec = 0.0_sp

      !$omp parallel do private(pos, rij2, qij, vec) reduction(+:gradient_vec) schedule(dynamic)
      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
        pos(:) = low_pos_vec(((i - 1)*low_dimension + 1):i*low_dimension) - low_pos_vec(((j - 1)*low_dimension + 1):j*low_dimension)
            qij = 1.0_sp/(1.0_sp + dot_product(pos, pos))*inv_z
            vec(:) = 4.0_sp*z*(exaggeration*pij(j, i) - (1 - pij(j, i))/(1 - qij)*qij)*qij*pos(:)
      gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) + vec(:)
      gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) - vec(:)
         end do
      end do
      !$omp end parallel do

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

   end subroutine loss_gradient_position

   subroutine loss_gradient_core(low_pos_vec, gradient_vec)
      implicit none
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)  :: low_pos_vec, gradient_vec
      real(kind=sp), dimension(low_dimension)               :: vec, pos
      real(kind=sp)                                         :: rij2, qij, dist
      integer                                               :: i, j

      gradient_vec = 0.0_sp

      !$omp parallel do private(pos, rij2, qij, dist, vec) reduction(+:gradient_vec) schedule(dynamic)
      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
        pos(:) = low_pos_vec(((i - 1)*low_dimension + 1):i*low_dimension) - low_pos_vec(((j - 1)*low_dimension + 1):j*low_dimension)
            rij2 = dot_product(pos, pos)
            qij = 1.0_sp/(1.0_sp + rij2)*inv_z
            vec(:) = 4.0_sp*z*(pij(j, i) - (1 - pij(j, i))/(1 - qij)*qij)*qij*pos(:)
      gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) + vec(:)
      gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) - vec(:)

            dist = sqrt(rij2)
            if (dist < point_radius(i) + point_radius(j)) then
               vec(:) = -pos/dist
               dist = (point_radius(i) + point_radius(j) - dist)/2.0_sp
               gradient_vec(((i-1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i-1)*low_dimension + 1):i*low_dimension) + vec(:)*dist*core_strength/2.0_sp
               gradient_vec(((j-1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j-1)*low_dimension + 1):j*low_dimension) - vec(:)*dist*core_strength/2.0_sp
            end if
         end do
      end do
      !$omp end parallel do

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

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

      end if

      previous_position = current_position
      previous_gradient = current_gradient

   end subroutine calculate_stepsize

   function calculating_cost_zero(pij) result(cost_zero)

      implicit none
      real(kind=sp), dimension(reduced_number_points, reduced_number_points), intent(in) :: pij
      real(kind=sp) :: cost_zero
      integer:: i, j

      cost_zero = 0.0_sp

      !$omp parallel do reduction(+:cost_zero) collapse(2)
      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
            if (pij(j, i) > 0 .and. pij(j, i) < 1) then
               cost_zero = cost_zero + pij(j, i)*log(pij(j, i))*2.0_sp + (1 - pij(j, i))*log(1 - pij(j, i))*2.0_sp
            else if (pij(i, j) > 1) then
               cost_zero = cost_zero + pij(j, i)*log(pij(j, i))*2.0_sp
            else if (pij(i, j) < 0) then
               cost_zero = cost_zero + (1 - pij(j, i))*log(1 - pij(j, i))*2.0_sp
            end if
         end do
      end do
      !$omp end parallel do

   end function calculating_cost_zero

   subroutine gradient_matrix_addnoise(gradient_matrix, r)

      real(kind=sp)                              :: d
      real(kind=sp), intent(in)    :: r
      real(kind=sp), intent(inout)  :: gradient_matrix(low_dimension, reduced_number_points)
      real(kind=sp)                             :: noise_vector(reduced_number_points)

      noise_vector = 0.0_sp

      call random_number(d)

      d = r*d**(1/real(reduced_number_points, sp))

      call random_add_noise(noise_vector, 1.0_sp)

      noise_vector = d*noise_vector/norm2(noise_vector)

      gradient_matrix = gradient_matrix*(1.0_sp + spread(noise_vector, 1, low_dimension))

   end subroutine gradient_matrix_addnoise

   subroutine handle_growth_phase(j, point_radius, growth_step_limit)
      implicit none

      integer, intent(in)             :: j
      real(kind=sp), intent(inout)    :: point_radius(reduced_number_points)
      logical, intent(inout)          :: growth_step_limit

      if (j < growth_steps) then
         if (j < 2) then
            point_radius = point_radius*(real(j))/real(growth_steps)
         else
            point_radius = point_radius*(real(j)*1.0_sp/real(growth_steps))/((real(j) - 1.0_sp)*1.0_sp/real(growth_steps))
         end if
      end if

      if (j > (growth_steps + 100)) growth_step_limit = .false.

   end subroutine handle_growth_phase

   subroutine loss_gradient(gradient_matrix, cost, exaggeration, growth_switch)
      implicit none
      logical, intent(in)  :: growth_switch
      real(kind=sp), intent(in)  :: exaggeration
      real(kind=sp), intent(out) :: cost
      real(kind=sp), dimension(low_dimension, reduced_number_points), intent(out) :: gradient_matrix

      real(kind=sp), dimension(low_dimension)               :: vec, pos
      real(kind=sp)                                         :: z, rij2, qij, dist
      integer                                               :: i, j

      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)

      cost = cost_zero
      gradient_matrix = 0.0_sp

      !omp parallel do private(pos, rij2, qij, vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
            pos(:) = low_dimension_position(:, i) - low_dimension_position(:, j)
            qij = calculating_qij(i, j)
            vec(:) = 4.0_sp*z*(exaggeration*pij(j, i) - (1 - pij(j, i))/(1 - qij)*qij)*qij*pos(:)
            gradient_matrix(:, i) = gradient_matrix(:, i) + vec(:)
            gradient_matrix(:, j) = gradient_matrix(:, j) - vec(:)
            cost = cost - pij(j, i)*log(qij)*2.0_sp - (1 - pij(j, i))*log(1 - qij)*2.0_sp
         end do
      end do
      !omp end parallel do

      if (growth_switch) then
         !$omp parallel do private(pos, rij2, dist, vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
         do i = 1, reduced_number_points
            do j = i + 1, reduced_number_points
               pos(:) = low_dimension_position(:, i) - low_dimension_position(:, j)
               rij2 = dot_product(pos, pos)
               dist = sqrt(rij2)
               if (dist < point_radius(i) + point_radius(j)) then
                  vec(:) = -pos/dist
                  dist = (point_radius(i) + point_radius(j) - dist)/2.0_sp
                  gradient_matrix(:, i) = gradient_matrix(:, i) + vec(:)*dist*core_strength/2.0_sp
                  gradient_matrix(:, j) = gradient_matrix(:, j) - vec(:)*dist*core_strength/2.0_sp
                  cost = cost + dist**2/2.0_sp*core_strength
               end if
            end do
         end do
         !$omp end parallel do
         write (*, *) 'Growth phase, cost: ', cost
      end if

   end subroutine loss_gradient

   subroutine loss_gradient_vectorisation(gradient_matrix, cost)

      implicit none
      real(kind=sp), intent(inout)                                          :: gradient_matrix(low_dimension, reduced_number_points)
      real(kind=sp), dimension(low_dimension, reduced_number_points - 1)    :: vec_matrix, pos_matrix
      real(kind=sp), dimension(reduced_number_points - 1)                   :: rij2_vector, pij_vector, qij_vector, factor
      real(kind=sp), intent(out)                                            :: cost
      real(kind=sp), dimension(:), allocatable                              :: point_radius_packed, dist_packed
      real(kind=sp), dimension(:, :), allocatable                           :: vec_matrix_packed
      logical, dimension(:), allocatable                                    :: overlap_mask
      integer:: i, j

      cost = cost_zero
      gradient_matrix = 0.0_sp

      !omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor, dist_packed, overlap_mask, point_radius_packed, vec_matrix_packed, j) reduction(+:gradient_matrix,cost) schedule(dynamic)

      do i = 1, reduced_number_points
         j = reduced_number_points - i

         pos_matrix(:, :j) = spread(low_dimension_position(:, i), 2, j) - low_dimension_position(:, i + 1:reduced_number_points)
         rij2_vector(:j) = sum(pos_matrix(:, :j)*pos_matrix(:, :j), dim=1)
         qij_vector(:j) = inv_z/(1.0_sp + rij2_vector(:j))
         pij_vector(:j) = pij(i + 1:reduced_number_points, i)
         factor(:j) = 4.0_sp*z*(pij_vector(:j) - (1 - pij_vector(:j))/(1 - qij_vector(:j))*qij_vector(:j))*qij_vector(:j)
         vec_matrix(:, :j) = spread(factor(:j), 1, low_dimension)*pos_matrix(:, :j)
         gradient_matrix(:, i) = gradient_matrix(:, i) + sum(vec_matrix(:, :j), dim=2)
         gradient_matrix(:, i + 1:reduced_number_points) = gradient_matrix(:, i + 1:reduced_number_points) - vec_matrix(:, :j)
         cost = cost - sum(pij_vector(:j)*log(qij_vector(:j))*2.0_sp - (1 - pij_vector(:j))*log(1 - qij_vector(:j))*2.0_sp)

         overlap_mask = sqrt(rij2_vector(:j)) < (point_radius(i) + point_radius(i + 1:reduced_number_points))

         if (any(overlap_mask)) then
            allocate (dist_packed(count(overlap_mask)))
            allocate (point_radius_packed(count(overlap_mask)))
            allocate (vec_matrix_packed(low_dimension, count(overlap_mask)))

            dist_packed = pack(sqrt(rij2_vector(:j)), overlap_mask)

            point_radius_packed = pack(point_radius(i + 1:reduced_number_points), overlap_mask)

 vec_matrix_packed = -pos_matrix(:, pack([(j, j=1, reduced_number_points - i)], overlap_mask))/spread(dist_packed, 1, low_dimension)

            dist_packed = (point_radius(i) + point_radius_packed - dist_packed)/2.0_sp

            gradient_matrix(:,i)= gradient_matrix(:,i)+ sum(vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength/2.0_sp, dim=2)
            gradient_matrix(:,pack([(i+j, j=1, reduced_number_points-i)],overlap_mask))= gradient_matrix(:,pack([(i+j, j=1, reduced_number_points-i)],overlap_mask))-vec_matrix_packed* spread(dist_packed, 1, low_dimension)*core_strength/2.0_sp

            cost = cost + sum(dist_packed*dist_packed/2.0_sp*core_strength)

            deallocate (vec_matrix_packed)
            deallocate (dist_packed)
            deallocate (point_radius_packed)
         end if
      end do

      !omp end parallel do

      call gradient_matrix_addnoise(gradient_matrix, 1e-2_sp)

   end subroutine loss_gradient_vectorisation

end MODULE optimisation

