MODULE mpi_model

   USE mpi
   USE data_reader
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   use optimisation
   use parameters

   IMPLICIT NONE

contains

   subroutine tpsd_mpi(threshold, maxsteps, rank, nranks)

      integer, intent(in)       :: maxsteps, rank, nranks
      real(kind=sp), intent(in) :: threshold
      real(kind=sp)             :: exaggeration
      real(kind=sp)             :: cost, step_size, gradient_norm, running_gradient_norm
      real(kind=sp)             :: gradient_matrix(low_dimension, reduced_number_points)
      integer                   :: i, j
      logical                   :: growth_step_limit = .true.

      if (rank == 0) then
         call start_timer()
      end if

      call initialize_variables(i, j, gradient_norm, running_gradient_norm)

      ! Main optimization loop
      do while ((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100 + exag_cutoff))) .and. (i < maxsteps)))

         i = i + 1

         exaggeration = merge(1.0_sp, exaggeration_init, i > exag_cutoff)

         call loss_gradient_position_mpi(gradient_matrix, cost, exaggeration, rank, nranks)

         step_size = calculate_stepsize(low_dimension_position, gradient_matrix, init=((i == 1) .or. (i == exag_cutoff)))

         low_dimension_position = low_dimension_position - step_size*gradient_matrix

         gradient_norm = abs(sum(step_size*gradient_matrix*gradient_matrix))

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

      end do

      if (rank == 0) then
         write (*, *) 'Growth phase'
      end if

      do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i + j < maxsteps))
         j = j + 1

         call loss_gradient_core_mpi(gradient_matrix, cost, rank, nranks)

         step_size = calculate_stepsize(low_dimension_position, gradient_matrix, init=((i == 1) .or. (i == exag_cutoff)))

         low_dimension_position = low_dimension_position - step_size*gradient_matrix

         gradient_norm = abs(sum(step_size*gradient_matrix*gradient_matrix))

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         call handle_growth_phase(j, growth_step_limit)
      end do

      if (rank == 0) then
         final_cost = cost
         write (*, *) 'Final cost: ', final_cost
         call stop_timer()
         write (*, *) 'Time: ', elapsed_time()
      end if
   end subroutine tpsd_mpi

   subroutine loss_gradient_position_mpi(gradient_matrix, cost, exaggeration, rank, nranks)
      implicit none
      real(kind=sp), intent(in)    :: exaggeration
      integer, intent(in)          :: rank, nranks
      real(kind=sp), intent(inout) :: cost
      real(kind=sp), intent(inout) :: gradient_matrix(low_dimension, reduced_number_points)
      real(kind=sp)                :: local_cost
      real(kind=sp)                :: local_gradient_matrix(low_dimension, reduced_number_points)

      real(kind=sp), dimension(:, :), allocatable :: vec_matrix, pos_matrix
      real(kind=sp), dimension(:), allocatable :: rij2_vector, pij_vector, qij_vector, factor
      integer :: i, j, k, ierr, status(MPI_STATUS_SIZE), source, tag, task, reply, finished, start, end
      logical :: flag

      cost = 0.0_sp
      gradient_matrix = 0.0_sp
      finished = 0

      allocate (pos_matrix(low_dimension, reduced_number_points))
      allocate (vec_matrix(low_dimension, reduced_number_points))
      allocate (rij2_vector(reduced_number_points))
      allocate (pij_vector(reduced_number_points))
      allocate (qij_vector(reduced_number_points))
      allocate (factor(reduced_number_points))

      if (rank == 0) then
         task = 0
         do k = 1, nranks - 1
            call MPI_Send(task, 1, MPI_INTEGER, k, k, MPI_COMM_WORLD, ierr)
            task = task + 1
         end do

         do while (finished < nranks - 1)
            call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, status, ierr)
            if (flag) then
               source = status(MPI_SOURCE)
               tag = status(MPI_TAG)

     call MPI_Recv(local_gradient_matrix, low_dimension*reduced_number_points, MPI_REAL4, source, tag, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_cost, 1, MPI_REAL4, source, tag, MPI_COMM_WORLD, status, ierr)

               if (task < nranks - 1) then
                  call MPI_Send(task, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  task = task + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  finished = finished + 1
               end if
            end if
         end do

      else
         do
            call MPI_Recv(task, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, status, ierr)
            if (task == -1) exit

            call work_distribution(rank, nranks, reduced_number_points, start, end)

            local_cost = 0.0_sp
            !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor) reduction(+:local_gradient_matrix, local_cost) schedule(dynamic)
            do i = start, end

               pos_matrix = spread(low_dimension_position(:, i), 2, reduced_number_points - i) - low_dimension_position(:, i + 1:reduced_number_points)
               rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
               qij_vector = inv_z/(1.0_sp + rij2_vector)
               pij_vector = pij(i + 1:reduced_number_points, i)
               factor = 4.0_sp*z*(exaggeration*pij_vector - (1 - pij_vector)/(1 - qij_vector)*qij_vector)*qij_vector
               vec_matrix = spread(factor, 1, low_dimension)*pos_matrix
               local_gradient_matrix(:, i) = local_gradient_matrix(:, i) + sum(vec_matrix, dim=2)
          local_gradient_matrix(:, i + 1:reduced_number_points) = local_gradient_matrix(:, i + 1:reduced_number_points) - vec_matrix
               local_cost = local_cost - sum(pij(i + 1:reduced_number_points, i) * log(qij_vector) * 2.0_sp - (1 - pij(i + 1:reduced_number_points, i)) * log(1 - qij_vector) * 2.0_sp)

            end do
            !$omp end parallel do

            call MPI_Send(local_gradient_matrix, low_dimension*reduced_number_points, MPI_REAL4, 0, tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_cost, 1, MPI_REAL4, 0, tag, MPI_COMM_WORLD, ierr)
         end do
      end if

      deallocate (pos_matrix)
      deallocate (vec_matrix)
      deallocate (rij2_vector)
      deallocate (pij_vector)
      deallocate (qij_vector)
      deallocate (factor)

      ! Reduce results across all processes
    call MPI_Allreduce(MPI_IN_PLACE, gradient_matrix, low_dimension*reduced_number_points, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cost, 1, MPI_REAL4, MPI_SUM, MPI_COMM_WORLD, ierr)

      call gradient_matrix_addnoise(gradient_matrix, 1e-2_sp)
   end subroutine loss_gradient_position_mpi

   subroutine loss_gradient_core_mpi(gradient_matrix, cost, rank, nranks)

      integer, intent(in)                         :: rank, nranks
      real(kind=sp), intent(inout)                :: cost
      real(kind=sp), intent(inout)                :: gradient_matrix(low_dimension, reduced_number_points)
      real(kind=sp), dimension(:, :), allocatable :: vec_matrix, pos_matrix
      real(kind=sp), dimension(:), allocatable    :: rij2_vector, pij_vector, qij_vector, factor
      real(kind=sp), dimension(:), allocatable    :: point_radius_packed, dist_packed
      real(kind=sp), dimension(:, :), allocatable :: vec_matrix_packed
      logical, dimension(:), allocatable          :: overlap_mask
      real(kind=sp)                               :: local_cost
      real(kind=sp)                               :: local_gradient_matrix(low_dimension, reduced_number_points)
      integer :: i, j, k, ierr, status(MPI_STATUS_SIZE), source, tag, task, reply, finished, start, end
      logical :: flag

      cost = 0.0_sp
      gradient_matrix = 0.0_sp
      finished = 0

      allocate (pos_matrix(low_dimension, reduced_number_points))
      allocate (vec_matrix(low_dimension, reduced_number_points))
      allocate (rij2_vector(reduced_number_points))
      allocate (pij_vector(reduced_number_points))
      allocate (qij_vector(reduced_number_points))
      allocate (factor(reduced_number_points))

      if (rank == 0) then
         task = 0
         do k = 1, nranks - 1
            call MPI_Send(task, 1, MPI_INTEGER, k, k, MPI_COMM_WORLD, ierr)
            task = task + 1
         end do

         do while (finished < nranks - 1)
            call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, flag, status, ierr)
            if (flag) then
               source = status(MPI_SOURCE)
               tag = status(MPI_TAG)

     call MPI_Recv(local_gradient_matrix, low_dimension*reduced_number_points, MPI_REAL4, source, tag, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_cost, 1, MPI_REAL4, source, tag, MPI_COMM_WORLD, status, ierr)

               if (task < nranks - 1) then
                  call MPI_Send(task, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  task = task + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  finished = finished + 1
               end if
            end if
         end do

      else
         do
            call MPI_Recv(task, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, status, ierr)
            if (task == -1) exit

            call work_distribution(rank, nranks, reduced_number_points, start, end)

            local_cost = 0.0_sp
            !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor, dist_packed, overlap_mask, point_radius_packed, vec_matrix_packed) reduction(+:gradient_matrix, local_cost) schedule(dynamic)
            do i = start, end

                    pos_matrix = spread(low_dimension_position(:, i), 2, reduced_number_points - i) - low_dimension_position(:, i + 1:reduced_number_points)
               rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
               qij_vector = inv_z/(1.0_sp + rij2_vector)
               pij_vector = pij(i + 1:reduced_number_points, i)
               factor = 4.0_sp*z*(pij_vector - (1 - pij_vector)/(1 - qij_vector)*qij_vector)*qij_vector
               vec_matrix = spread(factor, 1, low_dimension)*pos_matrix
               local_gradient_matrix(:, i) = local_gradient_matrix(:, i) + sum(vec_matrix, dim=2)
          local_gradient_matrix(:, i + 1:reduced_number_points) = local_gradient_matrix(:, i + 1:reduced_number_points) - vec_matrix
                    local_cost = local_cost - sum(pij(i + 1:reduced_number_points, i) * log(qij_vector) * 2.0_sp - (1 - pij(i + 1:reduced_number_points, i)) * log(1 - qij_vector) * 2.0_sp)

               overlap_mask = sqrt(rij2_vector) < (point_radius(i) + point_radius(i + 1:reduced_number_points))

               if (any(overlap_mask)) then
                  allocate (dist_packed(count(overlap_mask)))
                  allocate (point_radius_packed(count(overlap_mask)))
                  allocate (vec_matrix_packed(low_dimension, count(overlap_mask)))

                  dist_packed = pack(sqrt(rij2_vector), overlap_mask)
                  point_radius_packed = pack(point_radius(i + 1:reduced_number_points), overlap_mask)

 vec_matrix_packed = -pos_matrix(:, pack([(j, j=1, reduced_number_points - i)], overlap_mask))/spread(dist_packed, 1, low_dimension)
                  dist_packed = (point_radius(i) + point_radius_packed - dist_packed)/2.0_sp

                        local_gradient_matrix(:, i) = local_gradient_matrix(:, i) + sum(vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength / 2.0_sp, dim = 2)
                        local_gradient_matrix(:, pack([(i + j, j = 1, reduced_number_points - i)], overlap_mask)) = local_gradient_matrix(:, pack([(i + j, j = 1, reduced_number_points - i)], overlap_mask)) - vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength / 2.0_sp

                  local_cost = local_cost + sum(dist_packed*dist_packed/2.0_sp*core_strength)

                  deallocate (vec_matrix_packed)
                  deallocate (dist_packed)
                  deallocate (point_radius_packed)
               end if
            end do
            !$omp end parallel do
         end do

         call MPI_Send(local_gradient_matrix, low_dimension*reduced_number_points, MPI_REAL4, 0, tag, MPI_COMM_WORLD, ierr)
         call MPI_Send(local_cost, 1, MPI_REAL4, 0, tag, MPI_COMM_WORLD, ierr)
      end if

    call MPI_Allreduce(MPI_IN_PLACE, gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cost, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      call gradient_matrix_addnoise(gradient_matrix, 1e-2_sp)
   end subroutine loss_gradient_core_mpi

   subroutine work_distribution(rank, nranks, reduced_number_points, start, end)
      implicit none
      integer, intent(in)  :: rank, nranks, reduced_number_points
      integer, intent(out) :: start, end
      integer              :: i, current_load, total_load, load_per_rank

      total_load = reduced_number_points*(reduced_number_points - 1)
      load_per_rank = total_load/nranks
      current_load = 0
      start = 1
      end = 0

      do i = 1, reduced_number_points - 1
         current_load = current_load + (reduced_number_points - i)
         if (current_load >= load_per_rank*rank .and. end == 0) then
            start = i
         end if
         if (current_load >= load_per_rank*(rank + 1)) then
            end = i
            exit
         end if
      end do

      if (end == 0) end = reduced_number_points - 1
   end subroutine work_distribution

end module mpi_model
