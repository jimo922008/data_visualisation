MODULE mpi_model 
    
    USE data_reader
    USE initialisation
    USE high_dimension
    USE low_dimension_probability
    use optimisation
    use constants

    IMPLICIT NONE

contains 

    SUBROUTINE master_process(size, threshold, maxsteps)

        USE MPI

        IMPLICIT NONE

        INTEGER, INTENT(IN)       :: size, maxsteps
        REAL(kind=dp), INTENT(IN) :: threshold
        INTEGER :: ierr, task, source, tag, completed_tasks
        INTEGER :: i, j
        INTEGER :: status(MPI_STATUS_SIZE)
        REAL(kind=dp) :: cost, step_size, gradient_norm, running_gradient_norm
        REAL(kind=dp) :: gradient_vector(low_dimension*reduced_number_points), position_vector(low_dimension*reduced_number_points)

        task = 0
        completed_tasks = 0 

        call initialise_variables(cost, gradient_norm, running_gradient_norm, i, position_vector)

        do while((running_gradient_norm > log10(threshold) .or. (i < 100+nexag)) .and. (i < maxsteps))
                
            i = i + 1

            exageration = merge(1.0_dp, exageration, i > nexag)

            call distribute_tasks(position_vector, gradient_vector, nranks)

            call gather_results(gradient_vector, cost, local_gradient_matrix, local_cost)

            step_size = calculate_stepsize(position_vector, gradient_vector, init = ((i == 1) .or. (i == nexag)))

            position_vector=position_vector - step_size*gradient_vector

            low_dimension_position=reshape(position_vector,(/low_dimension,reduced_number_points/))

            gradient_norm=abs(dot_product(step_size*gradient_vector,gradient_vector))

            running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

        end do

        growth_switch = .true.
        write (*,*) 'Growth phase'
        growth_start_step = i
    
        do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i<maxsteps))
            i = i + 1

            call distribute_tasks(position_vector, gradient_vector, nranks)

            call gather_results(gradient_vector, cost, local_gradient_matrix, local_cost)

            step_size = calculate_stepsize(position_vector, gradient_vector)

            position_vector=position_vector - step_size*gradient_vector

            low_dimension_position=reshape(position_vector,(/low_dimension,reduced_number_points/))

            gradient_norm=abs(dot_product(step_size*gradient_vector,gradient_vector))

            running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

            call handle_growth_phase(i, growth_start_step, growth_step_limit)

        end do

        write (*,*) 'Final cost: ', final_cost

    end subroutine master_process
    
    
    subroutine worker_process(rank, threshold, maxsteps)
        implicit none

        integer, intent(in) :: rank, maxsteps
        real(kind=dp), intent(in) :: threshold

        integer :: ierr, task, source, tag
        integer :: i, j
        integer :: status(MPI_STATUS_SIZE)
        real(kind=dp) :: cost, step_size, gradient_norm, running_gradient_norm
        real(kind=dp) :: gradient_vector(low_dimension*reduced_number_points), position_vector(low_dimension*reduced_number_points)
        logical :: growth_switch = .false.
        logical :: growth_step_limit = .true.

        call initialise_variables(cost, gradient_norm, running_gradient_norm, i, position_vector)

        call MPI_Recv(task, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)

        call loss_gradient_vectorisation(position_vector, gradient_vector, cost)

        call MPI_Send(gradient_vector, low_dimension*reduced_number_points, MPI_REAL, 0, 0, MPI_COMM_WORLD, ierr)
    
    end subroutine worker_process

    subroutine distribute_work(position_vector, current_iteration, nranks)
        
        implicit none
        real(kind=dp), intent(in) :: position_vector(:)
        
        integer, intent(in) :: current_iteration, nranks
        integer :: ierr, tag, task, j

        if (task <= reduced_number_points) then
            call MPI_Send(task, 1, MPI_INTEGER, j, 0, MPI_COMM_WORLD, ierr)
            task = task + 1
        else 
            call MPI_Send(-1, 1, MPI_INTEGER, j, 0, MPI_COMM_WORLD, ierr)
        end if
    end subroutine distribute_work
    
    subroutine gather_results(gradient_vector, cost, local_gradient_matrix, local_cost)
        
        implicit none
        real(kind=dp), intent(inout) :: gradient_vector(:)
        real(kind=dp), intent(inout) :: cost
        real(kind=dp), intent(in) :: local_gradient_matrix(:,:), local_cost(:)
        
        integer :: ierr, tag, source
        integer :: status(MPI_STATUS_SIZE)
        integer :: i, j
        
        call MPI_Recv(gradient_vector, low_dimension*reduced_number_points, MPI_REAL, source, 0, MPI_COMM_WORLD, status, ierr)
        call MPI_Recv(cost, 1, MPI_REAL, source, 0, MPI_COMM_WORLD, status, ierr)
        
        completed_tasks = completed_tasks + 1
        
    end subroutine gather_results
END MODULE mpi_model