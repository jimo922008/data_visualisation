MODULE mpi_model

   USE mpi
   USE data_reader
   USE initialisation
   USE high_dimension
   USE low_dimension_probability
   use optimisation
   use parameters
   use timing_module

   IMPLICIT NONE

contains

   subroutine tpsd_mpi(pij_1d, low_pos_vec, threshold, maxsteps, rank, nranks)

      ! Arguments
      integer, intent(in)                                                               :: maxsteps, rank, nranks
      real(kind=sp), intent(in)                                                         :: threshold
      real(kind=sp), dimension(reduced_number_points*reduced_number_points), intent(in) :: pij_1d
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)      :: low_pos_vec

      ! Local variables
      real(kind=sp)                                                                     :: exaggeration, step_size, gradient_norm, running_gradient_norm, z, inv_z
      real(kind=sp), dimension(:), allocatable                                          :: gradient_vec, velocity
      integer                                                                           :: i, j
      logical                                                                           :: growth_step_limit = .true.

      ! MPI variables
      real(kind=sp), dimension(:), allocatable                                          :: local_gradient_vec
      real(kind=sp), dimension(:, :), allocatable                                        :: recv_gradient_vec
      integer, dimension(nranks - 1)                                                      :: requests
      integer                                                                           :: request
      integer                                                                           :: k, start, end, ierr, task, status(MPI_STATUS_SIZE), recv_status(MPI_STATUS_SIZE)

      task = 1
      i = 0
      j = 0
      gradient_norm = huge(1.0_sp)
      running_gradient_norm = 0.0_sp
      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)
      inv_z = 1.0_sp/z

      call work_distribution(rank, nranks, reduced_number_points, start, end)

      if (rank == 0) then

         allocate (local_gradient_vec(low_dimension*reduced_number_points))
         allocate (velocity(low_dimension*reduced_number_points))
         allocate (gradient_vec(low_dimension*reduced_number_points))
         allocate (recv_gradient_vec(low_dimension*reduced_number_points, (nranks - 1)))

         velocity = 0.0_sp
         gradient_vec = 0.0_sp
         local_gradient_vec = 0.0_sp

         do while ((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100 + exag_cutoff))) .and. (i < maxsteps)))

            i = i + 1

            exaggeration = merge(1.0_sp, exaggeration_init, i > exag_cutoff)

            do k = 1, nranks - 1
               call MPI_Isend(exaggeration, 1, MPI_REAL4, k, 0, MPI_COMM_WORLD, requests(k), ierr)
            end do

            call loss_gradient_position_mpi(pij_1d, low_pos_vec, local_gradient_vec, exaggeration, start, end)

            gradient_vec = local_gradient_vec

            do k = 1, nranks - 1
  call MPI_Irecv(recv_gradient_vec(:, k), (low_dimension*reduced_number_points), MPI_REAL4, k, 1, MPI_COMM_WORLD, requests(k), ierr)
            end do

            call MPI_Waitall(nranks - 1, requests, MPI_STATUSES_IGNORE, ierr)

            gradient_vec = gradient_vec + sum(recv_gradient_vec, dim=2)

            step_size = calculate_stepsize(low_pos_vec, gradient_vec, init=((i == 1) .or. (i == exag_cutoff)))

            velocity = momentum_coeff*velocity - step_size*gradient_vec

            low_pos_vec = low_pos_vec + velocity

            do k = 1, nranks - 1
               call MPI_Isend(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 2, MPI_COMM_WORLD, requests(k), ierr)
            end do

            call MPI_Waitall(nranks - 1, requests, MPI_STATUSES_IGNORE, ierr)

            gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

            running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/min(i, 100)

            write (*, *) 'Iteration: ', i, ' Gradient norm: ', running_gradient_norm, 'step size: ', step_size

         end do

         write (*, *) 'Growth phase...'

         exaggeration = -1

         do k = 1, nranks - 1
            call MPI_Isend(exaggeration, 1, MPI_INTEGER, k, 0, MPI_COMM_WORLD, requests(k), ierr)
         end do

         call MPI_Waitall(nranks - 1, requests, MPI_STATUSES_IGNORE, ierr)

         do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i + j < maxsteps))

            j = j + 1

            do k = 1, nranks - 1
               call MPI_Isend(task, 1, MPI_REAL4, k, 0, MPI_COMM_WORLD, requests(k), ierr)
            end do

            print *, 'Rank: ', rank, 'Task sent'

            call loss_gradient_core_mpi(pij_1d, low_pos_vec, local_gradient_vec, start, end)

            print *, 'Rank: ', rank, 'Gradient calculated'

            gradient_vec = local_gradient_vec

            do k = 1, nranks - 1
  call MPI_Irecv(recv_gradient_vec(:, k), (low_dimension*reduced_number_points), MPI_REAL4, k, 1, MPI_COMM_WORLD, requests(k), ierr)
            end do

            call MPI_Waitall(nranks - 1, requests, MPI_STATUSES_IGNORE, ierr)

            gradient_vec = gradient_vec + sum(recv_gradient_vec, dim=2)

            step_size = calculate_stepsize(low_pos_vec, gradient_vec, init=(j == 1))

            velocity = momentum_coeff*velocity - step_size*gradient_vec

            low_pos_vec = low_pos_vec + velocity

            do k = 1, nranks - 1
               call MPI_Isend(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, k, 2, MPI_COMM_WORLD, requests(k), ierr)
            end do

            call MPI_Waitall(nranks - 1, requests, MPI_STATUSES_IGNORE, ierr)

            gradient_norm = dot_product(step_size*gradient_vec, gradient_vec)

            running_gradient_norm = running_gradient_norm + (log10(gradient_norm) - running_gradient_norm)/100

            call handle_growth_phase(j, growth_step_limit)

            ! write (*, *) 'Iteration: ', i+j, ' Gradient norm: ', running_gradient_norm, 'step size: ', step_size, 'point radius: ', sum(point_radius)

         end do

         do k = 1, nranks
            call MPI_Isend(-1, 1, MPI_INTEGER, k, 0, MPI_COMM_WORLD, requests(k), ierr)
         end do

         call MPI_Waitall(nranks, requests, MPI_STATUSES_IGNORE, ierr)

         low_dimension_position = reshape(low_pos_vec, (/low_dimension, reduced_number_points/))

      else

         allocate (local_gradient_vec(low_dimension*reduced_number_points))
         local_gradient_vec = 0.0_sp

         do
            call MPI_Irecv(exaggeration, 1, MPI_REAL4, 0, 0, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)
            if (exaggeration == -1) exit

            call loss_gradient_position_mpi(pij_1d, low_pos_vec, local_gradient_vec, exaggeration, start, end)

            call MPI_Isend(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)

            call MPI_Irecv(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 2, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)

         end do

         do

            print *, 'Rank: ', rank, 'Waiting for task'
            call MPI_Irecv(task, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)
            if (task == -1) exit

            print *, 'Rank: ', rank, 'Task received'

            call loss_gradient_core_mpi(pij_1d, low_pos_vec, local_gradient_vec, start, end)

            print *, 'Rank: ', rank, 'Gradient calculated'

            call MPI_Isend(local_gradient_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 1, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)

            call MPI_Irecv(low_pos_vec, low_dimension*reduced_number_points, MPI_REAL4, 0, 2, MPI_COMM_WORLD, request, ierr)
            call MPI_Wait(request, MPI_STATUSES_IGNORE, ierr)

         end do
      end if

   end subroutine tpsd_mpi

   subroutine loss_gradient_position_mpi(pij_1d, low_pos_vec, gradient_vec, exaggeration, start, end)
      implicit none
      real(kind=sp), intent(in)                                                           :: exaggeration
      integer, intent(in)                                                                 :: start, end
      real(kind=sp), dimension(reduced_number_points*reduced_number_points), intent(in)   :: pij_1d
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)        :: low_pos_vec, gradient_vec
      real(kind=sp), dimension(low_dimension)                                             :: vec, pos
      real(kind=sp) :: qij
      integer :: i, j, ierr

      gradient_vec = 0.0_sp
      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)
      inv_z = 1.0_sp/z

      !$omp parallel do private(pos, qij, vec) reduction(+:gradient_vec) schedule(dynamic)
      do i = start, end
         do j = i + 1, reduced_number_points
        pos(:) = low_pos_vec(((i - 1)*low_dimension + 1):i*low_dimension) - low_pos_vec(((j - 1)*low_dimension + 1):j*low_dimension)
            qij = 1.0_sp/(1.0_sp + dot_product(pos, pos))*inv_z
            vec(:) = 4.0_sp*z*(exaggeration*pij_1d((j - 1) * reduced_number_points + i) - (1 - pij_1d((j - 1) * reduced_number_points + i))/(1 - qij)*qij)*qij*pos(:)
      gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) + vec(:)
      gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) - vec(:)
         end do
      end do
      !$omp end parallel do

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

   end subroutine loss_gradient_position_mpi

   subroutine loss_gradient_core_mpi(pij_1d, low_pos_vec, gradient_vec, start, end)
      implicit none
      integer, intent(in)                                                                 :: start, end
      real(kind=sp), dimension(reduced_number_points*reduced_number_points), intent(in)   :: pij_1d
      real(kind=sp), dimension(low_dimension*reduced_number_points), intent(inout)        :: gradient_vec, low_pos_vec
      real(kind=sp), dimension(low_dimension)                                             :: vec, pos
      real(kind=sp)                                                                       :: qij, rij2, dist
      integer :: i, j, ierr

      gradient_vec = 0.0_sp
      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)
      inv_z = 1.0_sp/z

      !$omp parallel do private(pos, rij2, qij, dist, vec) reduction(+:gradient_vec) schedule(dynamic)
      do i = start, end
         do j = i + 1, reduced_number_points
        pos(:) = low_pos_vec(((i - 1)*low_dimension + 1):i*low_dimension) - low_pos_vec(((j - 1)*low_dimension + 1):j*low_dimension)
            rij2 = dot_product(pos, pos)
            qij = 1.0_sp/(1.0_sp + rij2)*inv_z
            vec(:) = 4.0_sp*z*(pij_1d((j - 1) * reduced_number_points + i) - (1 - pij_1d((j - 1) * reduced_number_points + i))/(1 - qij)*qij)*qij*pos(:)
      gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i - 1)*low_dimension + 1):i*low_dimension) + vec(:)
      gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j - 1)*low_dimension + 1):j*low_dimension) - vec(:)

            print *, 'gradient_vec: ', dot_product(gradient_vec, gradient_vec)

            dist = sqrt(rij2)
            if (dist < point_radius(i) + point_radius(j)) then
               vec(:) = -pos/dist
               dist = (point_radius(i) + point_radius(j) - dist)/2.0_sp
               gradient_vec(((i-1)*low_dimension + 1):i*low_dimension) = gradient_vec(((i-1)*low_dimension + 1):i*low_dimension) + vec(:)*dist*core_strength/2.0_sp
               gradient_vec(((j-1)*low_dimension + 1):j*low_dimension) = gradient_vec(((j-1)*low_dimension + 1):j*low_dimension) - vec(:)*dist*core_strength/2.0_sp
            end if

            print *, 'gradient_vec: ', dot_product(gradient_vec, gradient_vec)
         end do
      end do
      !$omp end parallel do

      call gradient_matrix_addnoise(gradient_vec, 1e-2_sp)

   end subroutine loss_gradient_core_mpi

   subroutine work_distribution(rank, nranks, reduced_number_points, start, end)
      implicit none
      integer, intent(in)  :: rank, nranks, reduced_number_points
      integer, intent(out) :: start, end
      integer              :: i, total_load, load_per_rank, cumulative_load, remainder

      total_load = reduced_number_points*(reduced_number_points - 1)/2
      load_per_rank = total_load/nranks
      remainder = total_load - load_per_rank*nranks

      start = 1
      end = 0
      cumulative_load = 0

      do i = 1, reduced_number_points

         cumulative_load = cumulative_load + (reduced_number_points - i)

         if (cumulative_load >= load_per_rank*(rank) .and. start == 1 .and. rank /= 0) then
            start = i
         end if

         if (cumulative_load >= load_per_rank*(rank + 1)) then
            end = i - 1
         end if

         if (end /= 0) exit

      end do

      if (rank == nranks - 1) end = reduced_number_points

   end subroutine work_distribution

end module mpi_model
