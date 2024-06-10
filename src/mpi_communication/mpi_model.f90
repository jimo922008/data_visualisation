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
      real(kind=sp), dimension(low_dimension*reduced_number_points) :: gradient_vec, low_pos_vec, local_gradient_vec, local_low_pos_vec, velocity
      integer                   :: i, j, k, start, end, ierr, task, status(MPI_STATUS_SIZE)
      logical                   :: growth_step_limit = .true.

      if (rank == 0) then
         call start_timer()
         task = 1
      end if

      call initialize_variables(i, j, gradient_norm, running_gradient_norm)
      call work_distribution(rank, nranks, reduced_number_points, start, end)
      low_pos_vec = reshape(low_dimension_position, (/low_dimension*reduced_number_points/))

      if (rank == 0) then

         do while ((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100 + exag_cutoff))) .and. (i < maxsteps)))

            i = i + 1

            exaggeration = merge(1.0_sp, exaggeration_init, i > exag_cutoff)

            call MPI_Bcast(exaggeration, 1, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

            do k = 1, nranks - 1
               call MPI_Send(task, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, ierr)
            end do

            call loss_gradient_position_mpi(low_pos_vec, gradient_vec, exaggeration, start, end)

            do k = 1, nranks - 1
               call MPI_Recv(local_low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 1, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 1, MPI_COMM_WORLD, status, ierr)
               gradient_vec = gradient_vec + local_gradient_vec
               low_pos_vec = low_pos_vec + local_low_pos_vec
            end do

            call MPI_Bcast(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

            step_size = calculate_stepsize(low_pos_vec, gradient_vec, init=((i == 1) .or. (i == exag_cutoff)))

            velocity = momentum_coeff*velocity - step_size*gradient_vec

            low_pos_vec = low_pos_vec + velocity

            gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

            running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

         end do

         write (*, *) 'Growth phase...'

         do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i + j < maxsteps))

            j = j + 1

            do k = 1, nranks - 1
               call MPI_Send(task, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, ierr)
            end do

            call loss_gradient_core_mpi(low_pos_vec, gradient_vec, start, end)

            do k = 1, nranks - 1
               call MPI_Recv(local_low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 1, MPI_COMM_WORLD, status, ierr)
               call MPI_Recv(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 1, MPI_COMM_WORLD, status, ierr)
               gradient_vec = gradient_vec + local_gradient_vec
               low_pos_vec = low_pos_vec + local_low_pos_vec
            end do

            call MPI_Bcast(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, MPI_COMM_WORLD, ierr)

            step_size = calculate_stepsize(low_pos_vec, gradient_vec, init=(j == 1))

            velocity = momentum_coeff*velocity - step_size*gradient_vec

            low_pos_vec = low_pos_vec + velocity

            gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

            running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/100

            call handle_growth_phase(j, growth_step_limit)

         end do

         task = -1
         do k = 1, nranks - 1
            call MPI_Send(task, 1, MPI_INTEGER, k, 1, MPI_COMM_WORLD, ierr)
         end do

      else

         do
            call MPI_Recv(task, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierr)
            if (task == -1) exit
            call loss_gradient_position_mpi(local_low_pos_vec, local_gradient_vec, exaggeration, start, end)
            call MPI_Send(local_low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, ierr)
         end do

         do
            call MPI_Recv(task, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, status, ierr)
            if (task == -1) exit
            call loss_gradient_core_mpi(local_low_pos_vec, local_gradient_vec, start, end)
            call MPI_Send(local_low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, ierr)
            call MPI_Send(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, ierr)
         end do

         final_cost = cost
         write (*, *) 'Final cost: ', final_cost
         call stop_timer()
         write (*, *) 'Time: ', elapsed_time()
      end if
   end subroutine tpsd_mpi

   subroutine loss_gradient_position_mpi(low_pos_vec, gradient_vec, exaggeration, start, end)
      implicit none
      real(kind=sp), intent(in)                                                           :: exaggeration
      integer, intent(in)                                                                 :: start, end
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)        :: gradient_vec, low_pos_vec
      real(kind=sp), dimension(low_dimension)                                             :: vec, pos
      real(kind=sp) :: qij
      integer :: i, j, ierr

      gradient_vec = 0.0_sp

      !omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor) reduction(+:local_gradient_matrix, local_cost) schedule(dynamic)
      do i = start, end
         do j = i + 1, reduced_number_points
        pos(:) = low_pos_vec(((i - 1)*low_dimension + 1):i*low_dimension) - low_pos_vec(((j - 1)*low_dimension + 1):j*low_dimension)
            qij = 1.0_sp/(1.0_sp + dot_product(pos, pos))*inv_z
            vec(:) = 4.0_sp*z*(exaggeration*pij(j, i) - (1 - pij(j, i))/(1 - qij)*qij)*qij*pos(:)
      gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) + vec(:)
      gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) - vec(:)
         end do

      end do
      !omp end parallel do

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

   end subroutine loss_gradient_position_mpi

   subroutine loss_gradient_core_mpi(low_pos_vec, gradient_vec, start, end)

      integer, intent(in)                                                                 :: start, end
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)       :: gradient_vec, low_pos_vec
      real(kind=sp), dimension(low_dimension)                                             :: vec, pos
      real(kind=sp)                                                                       :: qij, rij2, dist
      integer :: i, j, ierr

      gradient_vec = 0.0_sp

      !omp parallel private(pos, rij2, qij, dist, vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
      do i = start, end
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
      !omp end do
      !omp end parallel

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

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
