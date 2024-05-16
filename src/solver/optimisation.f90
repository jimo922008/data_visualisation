MODULE optimisation

    USE constants
    USE data_reader
    USE high_dimension
    USE low_dimension_probability

    IMPLICIT NONE

    PUBLIC :: tpsd
    PUBLIC :: loss_gradient
    PUBLIC :: calculate_stepsize

    real(kind=dp) :: cost_zero

contains

    subroutine tpsd(threshold,maxsteps)

        implicit none
       
        integer,       intent(in) :: maxsteps 
        real(kind=dp), intent(in) :: threshold
        real(kind=dp)             :: cost, final_cost, step_size, gradient_norm, running_gradient_norm
        real(kind=dp)             :: gradient_matrix(low_dimension, number_points), gradient_vector(low_dimension*number_points)
        real(kind=dp)             :: position_vector(low_dimension*number_points), centre_of_mass(low_dimension)
        integer                   :: i, log_interval, steps

        steps=0
        gradient_norm=huge(1.0_dp)
        running_gradient_norm=0.0_dp

        cost_zero=calculating_cost_zero(pij)
    
        position_vector=reshape(low_dimension_position,(/low_dimension*number_points/))

        do while((running_gradient_norm > threshold .or. i < 100) .and. i < maxsteps)
            i = i + 1

            call loss_gradient(gradient_matrix, cost)

            gradient_vector=reshape(gradient_matrix,(/low_dimension*number_points/))

            step_size = calculate_stepsize(position_vector, gradient_vector)

            write (*,*) 'Step size: ', step_size

            position_vector=position_vector - step_size*gradient_vector

            low_dimension_position=reshape(position_vector,(/low_dimension,number_points/))

            gradient_norm=abs(dot_product(step_size*gradient_vector,gradient_vector))

            running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

            write (*,*) 'Step: ', i, ' Gradient norm: ', gradient_norm, ' Running gradient norm: ', running_gradient_norm
        end do

        ! * Store final map cost
        final_cost=cost
        
    end subroutine tpsd

    function calculating_cost_zero (pij) result(cost_zero)

        implicit none
        real(kind=dp), dimension(number_points,number_points), intent(in) :: pij
        real(kind=dp) :: cost_zero
        integer:: i, j

        cost_zero=0.0_dp

        do i=1,number_points
            if (point_count(i)==0) cycle
            do j=i+1,number_points
                if (point_count(j)==0) cycle
                if (pij(i,j)>0 .and. pij(i,j)<1) then
                    cost_zero = cost_zero + pij(i, j) * log(pij(i, j)) * 2.0_dp + (1 - pij(i, j)) * log(1 - pij(i, j)) * 2.0_dp
                end if
            end do
        end do

    end function calculating_cost_zero

    subroutine loss_gradient(gradient_matrix, cost)
        implicit none
        real(kind=dp),                          intent(out) :: cost
        real(kind=dp), dimension(low_dimension, number_points), intent(out):: gradient_matrix
        real(kind=dp), dimension(low_dimension) :: vec 
        integer:: i, j
        
        gradient_matrix=0.0_dp
        cost = cost_zero

        call calculating_low_dimension_distance()
        call calculating_qij()

        !$omp parallel do private(vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
        do i=1,number_points
            if (point_count(i)==0) cycle
            do j=i+1,number_points
                if (point_count(j)==0) cycle
                vec(:) = 4.0_dp * number_points * (number_points-1) * ((pij(i,j) - (1-pij(i,j)) / (1-qij(i,j)) * qij(i,j)) * qij(i,j)) * (low_dimension_position(:,i) - low_dimension_position(:,j))
                gradient_matrix(:,i)=gradient_matrix(:,i)+vec(:)
                gradient_matrix(:,j)=gradient_matrix(:,j)-vec(:)
                cost=cost-pij(i,j)*log(qij(i,j))*2.0_dp-(1-pij(i,j))*log(1-qij(i,j))*2.0_dp
            end do
        end do
        !$omp end parallel do

        write (*,*) 'Cost: ', cost
        
    end subroutine loss_gradient


    function calculate_stepsize(current_position, current_gradient) result(step_size)

        real(kind=dp), dimension(:), intent(in)                       :: current_position, current_gradient
        real(kind=dp), save, allocatable, dimension(:)                :: previous_position, previous_gradient
        real(kind=dp), dimension(size(current_position))              :: position_change, gradient_change
        real(kind=dp)                                                 :: gradient_change_magnitude, position_gradient_dot_product
        real(kind=dp)                                                 :: step_size
        real(kind=dp), parameter                                      :: default_step_size=1e-1_dp

        if(.not.allocated(previous_position)) then
            allocate(previous_position(size(current_position)))
            allocate(previous_gradient(size(current_gradient)))
            previous_position=current_position
            previous_gradient=current_gradient
            step_size=default_step_size
        end if

        position_change=current_position-previous_position
        gradient_change=current_gradient-previous_gradient

        gradient_change_magnitude=sqrt(dot_product(gradient_change,gradient_change))
        position_gradient_dot_product=dot_product(position_change,gradient_change)
        
        if (gradient_change_magnitude==0.0_dp) then
            step_size=default_step_size
        else
            step_size=position_gradient_dot_product/gradient_change_magnitude
        end if 

        previous_position=current_position
        previous_gradient=current_gradient

    end function calculate_stepsize


end MODULE optimisation
