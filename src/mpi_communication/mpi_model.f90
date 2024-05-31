MODULE mpi_model

   include 'mpif.h'

   USE data_reader
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   use optimisation
   use parameters

   IMPLICIT NONE

contains

   subroutine tpsd_mpi(threshold, maxsteps, rank, nranks)

      implicit none
      integer, intent(in) :: maxsteps, rank, nranks
      real(kind=dp), intent(in) :: threshold
      real(kind=dp) :: cost, step_size, gradient_norm, running_gradient_norm
      real(kind=dp) :: gradient_matrix(low_dimension, reduced_number_points)
      integer :: i, j, ierr, work, start, end, flag, source, dest
      integer :: status(MPI_STATUS_SIZE)
      logical :: growth_step_limit = .true.

      if (rank == 0) then
         call start_timer()
      end if

      call initialize_variables(i, j, gradient_norm, running_gradient_norm)

      ! Main optimization loop
      do while ((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100 + exag_cutoff))) .and. (i < maxsteps)))

         i = i + 1
         exageration = merge(1.0_dp, exageration, i > exag_cutoff)

         call loss_gradient_position_mpi(gradient_matrix, cost, exageration, rank, nranks, i, maxsteps)

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

         call loss_gradient_core_mpi(gradient_matrix, cost, rank, nranks, i, maxsteps)

         step_size = calculate_stepsize(low_dimension_position, gradient_matrix, init=((i == 1) .or. (i == exag_cutoff)))

         low_dimension_position = low_dimension_position - step_size*gradient_matrix

         gradient_norm = abs(sum(step_size*gradient_matrix*gradient_matrix))

         running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         call handle_growth_phase(i, j, growth_step_limit)
      end do

      if (rank == 0) then
         final_cost = cost
         write (*, *) 'Final cost: ', final_cost
         call stop_timer()
         write (*, *) 'Time: ', elapsed_time()
      end if
   end subroutine tpsd_mpi

   subroutine loss_gradient_position_mpi(gradient_matrix, cost, exageration, rank, nranks, iter, maxsteps)
      implicit none
      real(kind=dp), intent(in) :: exageration
      real(kind=dp), intent(out) :: cost
      real(kind=dp), dimension(low_dimension, reduced_number_points), intent(inout) :: gradient_matrix
      real(kind=dp), dimension(:, :), allocatable :: vec_matrix, pos_matrix
      real(kind=dp), dimension(:), allocatable :: rij2_vector, pij_vector, qij_vector, factor
      integer :: i, ierr, work, start, end, status(MPI_STATUS_SIZE), source, dest
      integer, parameter :: tag = 1
      real(kind=dp) :: local_cost

      cost = 0.0_dp
      gradient_matrix = 0.0_dp

      if (rank == 0) then
         ! Master process
         work = 0
         do dest = 1, nranks - 1
            call MPI_Send(work, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
            work = work + 1
         end do

         ! Receive results from worker processes
         do while (work < reduced_number_points .or. flag > 0)
            call MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, flag, status, ierr)
            if (flag > 0) then
               source = status(MPI_SOURCE)
           call MPI_Recv(gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, source, tag, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_cost, 1, MPI_REAL8, source, tag, MPI_COMM_WORLD, status, ierr)
               cost = cost + local_cost

               if (work < reduced_number_points) then
                  call MPI_Send(work, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  work = work + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
               end if
            end if
         end do

      else
         ! Worker process
         do
            call MPI_Recv(work, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, status, ierr)
            if (work == -1) exit

            start = work*(reduced_number_points/nranks) + 1
            end = min((work + 1)*(reduced_number_points/nranks), reduced_number_points)

            local_cost = 0.0_dp
            !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor) reduction(+:gradient_matrix, local_cost) schedule(dynamic)
            do i = start, end
               allocate (pos_matrix(low_dimension, reduced_number_points - i))
               allocate (vec_matrix(low_dimension, reduced_number_points - i))
               allocate (rij2_vector(reduced_number_points - i))
               allocate (pij_vector(reduced_number_points - i))
               allocate (qij_vector(reduced_number_points - i))
               allocate (factor(reduced_number_points - i))

                    pos_matrix = spread(low_dimension_position(:, i), 2, reduced_number_points - i) - low_dimension_position(:, i + 1:reduced_number_points)
               rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
               qij_vector = inv_z/(1.0_dp + rij2_vector)
               pij_vector = pij(i + 1:reduced_number_points, i)
               factor = 4.0_dp*z*(exageration*pij_vector - (1 - pij_vector)/(1 - qij_vector)*qij_vector)*qij_vector
               vec_matrix = spread(factor, 1, low_dimension)*pos_matrix
               gradient_matrix(:, i) = gradient_matrix(:, i) + sum(vec_matrix, dim=2)
               gradient_matrix(:, i + 1:reduced_number_points) = gradient_matrix(:, i + 1:reduced_number_points) - vec_matrix
                    local_cost = local_cost - sum(pij(i + 1:reduced_number_points, i) * log(qij_vector) * 2.0_dp - (1 - pij(i + 1:reduced_number_points, i)) * log(1 - qij_vector) * 2.0_dp)

               deallocate (pos_matrix)
               deallocate (vec_matrix)
               deallocate (rij2_vector)
               deallocate (pij_vector)
               deallocate (qij_vector)
               deallocate (factor)
            end do
            !$omp end parallel do

            call MPI_Send(gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, 0, tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_cost, 1, MPI_REAL8, 0, tag, MPI_COMM_WORLD, ierr)
         end do
      end if

      ! Reduce results across all processes
    call MPI_Allreduce(MPI_IN_PLACE, gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cost, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      call gradient_matrix_addnoise(gradient_matrix, 1e-2_dp)
   end subroutine loss_gradient_position_mpi

   subroutine loss_gradient_core_mpi(gradient_matrix, cost, rank, size, iter, maxsteps)
      implicit none
      real(kind=dp), intent(inout) :: gradient_matrix(low_dimension, reduced_number_points)
      real(kind=dp), dimension(:, :), allocatable :: vec_matrix, pos_matrix
      real(kind=dp), dimension(:), allocatable :: rij2_vector, pij_vector, qij_vector, factor
      real(kind=dp), intent(out) :: cost
      real(kind=dp), dimension(:), allocatable :: point_radius_packed, dist_packed
      real(kind=dp), dimension(:, :), allocatable :: vec_matrix_packed
      logical, dimension(:), allocatable :: overlap_mask
      integer :: i, j, ierr, work, start, end, status(MPI_STATUS_SIZE), source, dest
      integer, parameter :: tag = 2
      real(kind=dp) :: local_cost

      cost = 0.0_dp
      gradient_matrix = 0.0_dp

      if (rank == 0) then
         ! Master process
         work = 0
         do dest = 1, size - 1
            call MPI_Send(work, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
            work = work + 1
         end do

         ! Receive results from worker processes
         do while (work < reduced_number_points .or. flag > 0)
            call MPI_Iprobe(MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, flag, status, ierr)
            if (flag > 0) then
               source = status(MPI_SOURCE)
           call MPI_Recv(gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, source, tag, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_cost, 1, MPI_REAL8, source, tag, MPI_COMM_WORLD, status, ierr)
               cost = cost + local_cost

               if (work < reduced_number_points) then
                  call MPI_Send(work, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
                  work = work + 1
               else
                  call MPI_Send(-1, 1, MPI_INTEGER, source, tag, MPI_COMM_WORLD, ierr)
               end if
            end if
         end do

      else
         ! Worker process
         do
            call MPI_Recv(work, 1, MPI_INTEGER, 0, tag, MPI_COMM_WORLD, status, ierr)
            if (work == -1) exit

            start = work*(reduced_number_points/size) + 1
            end = min((work + 1)*(reduced_number_points/size), reduced_number_points)

            local_cost = 0.0_dp
            !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor, dist_packed, overlap_mask, point_radius_packed, vec_matrix_packed) reduction(+:gradient_matrix, local_cost) schedule(dynamic)
            do i = start, end
               allocate (pos_matrix(low_dimension, reduced_number_points - i))
               allocate (vec_matrix(low_dimension, reduced_number_points - i))
               allocate (rij2_vector(reduced_number_points - i))
               allocate (pij_vector(reduced_number_points - i))
               allocate (qij_vector(reduced_number_points - i))
               allocate (factor(reduced_number_points - i))

                    pos_matrix = spread(low_dimension_position(:, i), 2, reduced_number_points - i) - low_dimension_position(:, i + 1:reduced_number_points)
               rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
               qij_vector = inv_z/(1.0_dp + rij2_vector)
               pij_vector = pij(i + 1:reduced_number_points, i)
               factor = 4.0_dp*z*(exageration*pij_vector - (1 - pij_vector)/(1 - qij_vector)*qij_vector)*qij_vector
               vec_matrix = spread(factor, 1, low_dimension)*pos_matrix
               gradient_matrix(:, i) = gradient_matrix(:, i) + sum(vec_matrix, dim=2)
               gradient_matrix(:, i + 1:reduced_number_points) = gradient_matrix(:, i + 1:reduced_number_points) - vec_matrix
                    local_cost = local_cost - sum(pij(i + 1:reduced_number_points, i) * log(qij_vector) * 2.0_dp - (1 - pij(i + 1:reduced_number_points, i)) * log(1 - qij_vector) * 2.0_dp)

               overlap_mask = sqrt(rij2_vector) < (point_radius(i) + point_radius(i + 1:reduced_number_points))

               if (any(overlap_mask)) then
                  allocate (dist_packed(count(overlap_mask)))
                  allocate (point_radius_packed(count(overlap_mask)))
                  allocate (vec_matrix_packed(low_dimension, count(overlap_mask)))

                  dist_packed = pack(sqrt(rij2_vector), overlap_mask)
                  point_radius_packed = pack(point_radius(i + 1:reduced_number_points), overlap_mask)

 vec_matrix_packed = -pos_matrix(:, pack([(j, j=1, reduced_number_points - i)], overlap_mask))/spread(dist_packed, 1, low_dimension)
                  dist_packed = (point_radius(i) + point_radius_packed - dist_packed)/2.0_dp

                        gradient_matrix(:, i) = gradient_matrix(:, i) + sum(vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength / 2.0_dp, dim = 2)
                        gradient_matrix(:, pack([(i + j, j = 1, reduced_number_points - i)], overlap_mask)) = gradient_matrix(:, pack([(i + j, j = 1, reduced_number_points - i)], overlap_mask)) - vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength / 2.0_dp

                  local_cost = local_cost + sum(dist_packed*dist_packed/2.0_dp*core_strength)

                  deallocate (vec_matrix_packed)
                  deallocate (dist_packed)
                  deallocate (point_radius_packed)
               end if

               deallocate (overlap_mask)
               deallocate (pos_matrix)
               deallocate (vec_matrix)
               deallocate (rij2_vector)
               deallocate (pij_vector)
               deallocate (qij_vector)
               deallocate (factor)
            end do
            !$omp end parallel do

            call MPI_Send(gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, 0, tag, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_cost, 1, MPI_REAL8, 0, tag, MPI_COMM_WORLD, ierr)
         end do
      end if

      ! Reduce results across all processes
    call MPI_Allreduce(MPI_IN_PLACE, gradient_matrix, low_dimension*reduced_number_points, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(MPI_IN_PLACE, cost, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)

      call gradient_matrix_addnoise(gradient_matrix, 1e-2_dp)
   end subroutine loss_gradient_core
