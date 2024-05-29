MODULE mpi_model 
    
    USE data_reader
    USE initialisation
    USE high_dimension
    USE low_dimension_probability
    use optimisation
    use constants


    IMPLICIT NONE

    integer                                    :: rank, nranks
    integer                                    :: i, j
    integer                                    :: pair_index, start_index, end_index, number_pairs, pairs_per_rank
    integer, dimension(:,:), allocatable       :: pairs
    real(kind=dp)                              :: cost_local, cost_global
    real(kind=dp), dimension(:,:), allocatable :: gradient_matrix_local, gradient_matrix_global

CONTAINS

    SUBROUTINE tpsd_mpi(threshold, maxsteps, nranks, ierr)
        implicit none

        integer,       intent(in) :: maxsteps, rank, nranks
        real(kind=dp), intent(in) :: threshold
        integer,       intent(out):: ierr
        real(kind=dp)             :: cost, step_size, gradient_norm, running_gradient_norm
        real(kind=dp)             :: gradient_vector(low_dimension*reduced_number_points)
        real(kind=dp)             :: ptb_vec(low_dimension*reduced_number_points), local_gptb(low_dimension*reduced_number_points)
        real(kind=dp)             :: position_vector(low_dimension*reduced_number_points)
        integer                   :: i, start_growth
        logical                   :: growth_switch = .false.
        logical                   :: growing = .true.

        ! * Initialise MPI parameters
        integer :: task, completed_tasks, worker_rank, status(MPI_STATUS_SIZE)
        integer, dimension(:), allocatable :: tasks 
        real(kind=dp), dimension(:), allocatable ::local_gradient_matrix
        real(kind=dp) :: local_cost
    
        !initialise values
        i=0
        gradient_norm=huge(1.0_dp)
        running_gradient_norm=0.0_dp

        cost_zero=calculating_cost_zero(pij)
    
        position_vector=reshape(low_dimension_position,(/low_dimension*reduced_number_points/))

        if (rank == 0) then 
            call start_timer()  

            task = 1

            completed_tasks = 0 

            do while((((running_gradient_norm > log10(threshold) .or. (i < 100+nexag))) .and. (i < maxsteps)) .or. (growing))
                
                i = i + 1

                if (i > nexag) exageration = 1.0_dp

                do j = 1, nranks-1
                    if (task <= reduced_number_points) then 
                        call MPI_Send(task, 1, MPI_INTEGER, j, 0, MPI_COMM_WORLD, ierr)
                        task = task + 1
                    else 
                        call MPI_Send(-1, 1, MPI_INTEGER, j, 0, MPI_COMM_WORLD, ierr)
                    end if 
                end do 

                do j = 1, nranks-1
                    call MPI_Recv(local_gradient_matrix, low_dimension*reduced_number_points, MPI_DOUBLE_PRECISION, j, 0, MPI_COMM_WORLD, status, ierr)
                    call MPI_Recv(local_cost, 1, MPI_DOUBLE_PRECISION, j, 0, MPI_COMM_WORLD, status, ierr)
                    gradient_matrix = gradient_matrix + local_gradient_matrix
                    cost_local = cost_local + local_cost
                end do
                
                

                call random_point_in_hypersphere(ptb_vec, 1e-2_dp)

                gptb = gradient_vector*(1.0_dp+ ptb_vec)

                step_size = abs(calculate_stepsize(position_vector, gptb, init = ((i == 1) .or. (i == nexag))))

                position_vector=position_vector - step_size*gptb

                low_dimension_position=reshape(position_vector,(/low_dimension,reduced_number_points/))

                gradient_norm=abs(dot_product(step_size*gptb,gradient_vector))

                running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

                if ((running_gradient_norm < log10(threshold* growth_tol_coeff)) .and. (i> nexag+100) .and. (.not.growth_switch)) then
                    growth_switch = .true.
                    write (*,*) 'Growth phase'
                    start_growth = i
                end if 

                if ((growth_switch) .and. (i > start_growth) .and. (i < start_growth + growth_steps)) then
                    if (i < start_growth + 2) then
                        point_radius = point_radius * ((real(i)-real(start_growth))/real(growth_steps))
                    else
                        point_radius = point_radius * ((real(i)-real(start_growth))*1.0_dp/real(growth_steps))&
                        /((real(i)-real(start_growth)-1.0_dp)*1.0_dp/real(growth_steps))
                    end if

                end if

                if(growth_switch .and. i > (start_growth+growth_steps+100)) growing=.false.

                write (*,*) 'start growth: ', start_growth
                write (*,*) 'Step: ', step_size, ' cost: ', cost, ' Running gradient norm: ', running_gradient_norm
                write (*,*) 'point_radius: ', sum(point_radius)
            end do

            call stop_timer()
            ! * Store final map cost
            final_cost=cost

            write (*,*) 'Final cost: ', final_cost

            write (*,*) 'Total time: ', elapsed_time()

        end if 

    end subroutine tpsd_mpi


END MODULE mpi_model