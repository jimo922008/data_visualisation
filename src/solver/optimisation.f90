MODULE optimisation

    USE constants
    USE data_reader
    USE high_dimension
    USE low_dimension_probability
    USE timing_module

    IMPLICIT NONE

    PUBLIC :: tpsd
    PUBLIC :: loss_gradient
    PUBLIC :: calculate_stepsize
    PUBLIC :: initialize_variables

    real(kind=dp) :: cost_zero, final_cost, z, inv_z

    integer                   :: exag_cutoff = 200
    real(kind=dp)             :: exageration = 5.0_dp 
    integer                   :: growth_steps= 500
    real(kind=dp)             :: core_strength=1.0_dp
    real(kind=dp)             :: growth_coeff=2.0_dp
    
contains

    subroutine tpsd(threshold, maxsteps)

        implicit none

        integer,       intent(in) :: maxsteps 
        real(kind=dp), intent(in) :: threshold
        real(kind=dp)             :: cost, step_size, gradient_norm, running_gradient_norm
        real(kind=dp)             :: gradient_vector(low_dimension*reduced_number_points)
        real(kind=dp)             :: position_vector(low_dimension*reduced_number_points)
        integer                   :: i, growth_start_step
        logical                   :: growth_switch = .false.
        logical                   :: growth_step_limit = .true.

  
        call initialize_variables(cost_zero, gradient_norm, running_gradient_norm, i, position_vector)

        do while((((running_gradient_norm > log10(threshold*growth_coeff) .or. (i < 100+exag_cutoff))) .and. (i < maxsteps)))
            
            i = i + 1

            exageration = merge (1.0_dp, exageration, i > exag_cutoff)

            call loss_gradient_position(gradient_vector, cost, exageration)

            step_size = calculate_stepsize(position_vector, gradient_vector, init = ((i == 1) .or. (i == exag_cutoff)))

            position_vector=position_vector - step_size*gradient_vector

            low_dimension_position=reshape(position_vector,(/low_dimension,reduced_number_points/))

            gradient_norm=abs(dot_product(step_size*gradient_vector,gradient_vector))

            running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

            write (*,*) 'Step: ', step_size, 'Cost: ', cost, 'Gradient norm: ', running_gradient_norm   

        end do 

        write (*,*) 'Growth phase'
        growth_start_step = i
    
        do while (((running_gradient_norm > log10(threshold)) .or. (growth_step_limit)) .and. (i<maxsteps))
            i = i + 1

            call loss_gradient_core(gradient_vector, cost)

            step_size = calculate_stepsize(position_vector, gradient_vector)

            position_vector=position_vector - step_size*gradient_vector

            low_dimension_position=reshape(position_vector,(/low_dimension,reduced_number_points/))

            gradient_norm=abs(dot_product(step_size*gradient_vector,gradient_vector))

            running_gradient_norm=running_gradient_norm + (log10(gradient_norm)-running_gradient_norm)/min(i,100)

            call handle_growth_phase(i, growth_start_step, growth_step_limit)

            write (*,*) 'Step: ', i, 'Cost: ', cost, 'Gradient norm: ', running_gradient_norm   
            write (*,*) 'Point radius: ', sum(point_radius)

        end do

        final_cost=cost

        write (*,*) 'Final cost: ', final_cost
        
    end subroutine tpsd

    subroutine initialize_variables(cost_zero, gradient_norm, running_gradient_norm, i, position_vector)
        implicit none
        real(kind=dp), intent(out) :: cost_zero, gradient_norm, running_gradient_norm
        integer, intent(out)       :: i
        real(kind=dp), dimension(low_dimension*reduced_number_points), intent(out) :: position_vector

        i=0
        gradient_norm=huge(1.0_dp)
        running_gradient_norm=0.0_dp
        cost_zero = calculating_cost_zero(pij)
        z = real(reduced_number_points, dp) * (real(reduced_number_points, dp)-1.0_dp)
        inv_z = 1.0_dp/z
        position_vector = reshape(low_dimension_position,(/low_dimension*reduced_number_points/))
    end subroutine initialize_variables

    subroutine loss_gradient_position(gradient_vector, cost, exageration)
        implicit none
        
        real(kind=dp),  intent(in)  :: exageration
        real(kind=dp),  intent(out) :: cost
        real(kind=dp),  intent(out) :: gradient_vector (low_dimension* reduced_number_points)

        real(kind=dp), dimension(low_dimension, reduced_number_points)      :: gradient_matrix
        real(kind=dp), dimension(:, :), allocatable                         :: vec_matrix, pos_matrix
        real(kind=dp), dimension(:), allocatable                            :: rij2_vector, pij_vector, qij_vector, factor
        integer                                                             :: i
        
        cost = cost_zero
        gradient_matrix=0.0_dp

        !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor) reduction(+:gradient_matrix,cost) schedule(dynamic)
        do i=1,reduced_number_points
            allocate(pos_matrix(low_dimension, reduced_number_points-i))
            allocate(vec_matrix(low_dimension, reduced_number_points-i))
            allocate(rij2_vector(reduced_number_points-i))
            allocate(pij_vector(reduced_number_points-i))
            allocate(qij_vector(reduced_number_points-i))
            allocate(factor(reduced_number_points-i))

            pos_matrix = spread(low_dimension_position(:,i), 2, reduced_number_points-i) -low_dimension_position(:,i+1:reduced_number_points)
            rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
            qij_vector = inv_z/(1.0_dp+rij2_vector)
            pij_vector = pij(i+1:reduced_number_points,i)
            factor = 4.0_dp * z * (exageration*pij_vector - (1-pij_vector) / (1-qij_vector) * qij_vector) * qij_vector             
            vec_matrix = spread(factor, 1, low_dimension) * pos_matrix
            gradient_matrix(:,i)=gradient_matrix(:,i)+sum(vec_matrix, dim=2)
            gradient_matrix(:,i+1:reduced_number_points)=gradient_matrix(:,i+1:reduced_number_points)-vec_matrix
            cost=cost-sum(pij(i+1:reduced_number_points,i)*log(qij_vector)*2.0_dp-(1-pij(i+1:reduced_number_points,i))*log(1-qij_vector)*2.0_dp)

            deallocate(pos_matrix)
            deallocate(vec_matrix)
            deallocate(rij2_vector)
            deallocate(pij_vector)
            deallocate(qij_vector)
            deallocate(factor)
    
        end do
        !$omp end parallel do

        gradient_vector = reshape(gradient_matrix,(/low_dimension*reduced_number_points/))

        call gradient_vector_addnoise(gradient_vector, 1e-2_dp)

    end subroutine loss_gradient_position

    subroutine loss_gradient_core (gradient_vector, cost)

        implicit none
        real(kind=dp)                                                         :: gradient_matrix(low_dimension, reduced_number_points)
        real(kind=dp), dimension(:, :), allocatable                           :: vec_matrix, pos_matrix
        real(kind=dp), dimension(:), allocatable                              :: rij2_vector, pij_vector, qij_vector, factor
        real(kind=dp), intent(out)                                            :: cost
        real(kind=dp), dimension(:), allocatable                              :: point_radius_packed, dist_packed
        real(kind=dp), dimension(:,:), allocatable                            :: vec_matrix_packed
        logical, dimension(:), allocatable                                    :: overlap_mask
        real(kind=dp), intent(out)                                            :: gradient_vector(low_dimension*reduced_number_points)
        integer:: i, j

        cost = cost_zero
        gradient_matrix=0.0_dp

        !$omp parallel do private(pos_matrix, rij2_vector, qij_vector, vec_matrix, pij_vector, factor, dist_packed, overlap_mask, point_radius_packed, vec_matrix_packed) reduction(+:gradient_matrix,cost) schedule(dynamic)
        do i=1,reduced_number_points
            allocate(pos_matrix(low_dimension, reduced_number_points-i))
            allocate(vec_matrix(low_dimension, reduced_number_points-i))
            allocate(rij2_vector(reduced_number_points-i))
            allocate(pij_vector(reduced_number_points-i))
            allocate(qij_vector(reduced_number_points-i))
            allocate(factor(reduced_number_points-i))

            pos_matrix = spread(low_dimension_position(:,i), 2, reduced_number_points-i) -low_dimension_position(:,i+1:reduced_number_points)
            rij2_vector = sum(pos_matrix*pos_matrix, dim=1)
            qij_vector = inv_z/(1.0_dp+rij2_vector)
            pij_vector = pij(i+1:reduced_number_points,i)
            factor = 4.0_dp * z * (exageration*pij_vector - (1-pij_vector) / (1-qij_vector) * qij_vector) * qij_vector             
            vec_matrix = spread(factor, 1, low_dimension) * pos_matrix
            gradient_matrix(:,i)=gradient_matrix(:,i)+sum(vec_matrix, dim=2)
            gradient_matrix(:,i+1:reduced_number_points)=gradient_matrix(:,i+1:reduced_number_points)-vec_matrix
            cost=cost-sum(pij(i+1:reduced_number_points,i)*log(qij_vector)*2.0_dp-(1-pij(i+1:reduced_number_points,i))*log(1-qij_vector)*2.0_dp)

            overlap_mask = sqrt(rij2_vector) < (point_radius(i)+point_radius(i+1:reduced_number_points))
        
            if (any(overlap_mask)) then
                allocate(dist_packed(count(overlap_mask)))
                allocate(point_radius_packed(count(overlap_mask)))
                allocate(vec_matrix_packed(low_dimension, count(overlap_mask)))

                dist_packed = pack(sqrt(rij2_vector), overlap_mask)

                point_radius_packed = pack(point_radius(i+1:reduced_number_points), overlap_mask)

                vec_matrix_packed= -pos_matrix(:, pack([(j, j=1, reduced_number_points-i)],overlap_mask))/spread(dist_packed, 1, low_dimension)
                
                dist_packed=(point_radius(i)+point_radius_packed-dist_packed)/2.0_dp
                
                gradient_matrix(:,i)= gradient_matrix(:,i)+ sum(vec_matrix_packed * spread(dist_packed, 1, low_dimension) * core_strength/2.0_dp, dim=2)
                gradient_matrix(:,pack([(i+j, j=1, reduced_number_points-i)],overlap_mask))= gradient_matrix(:,pack([(i+j, j=1, reduced_number_points-i)],overlap_mask))-vec_matrix_packed* spread(dist_packed, 1, low_dimension)*core_strength/2.0_dp

                cost=cost+sum(dist_packed*dist_packed/2.0_dp*core_strength)

                deallocate(vec_matrix_packed)
                deallocate(dist_packed)
                deallocate(point_radius_packed)
            end if 

            deallocate(overlap_mask)
            deallocate(pos_matrix)
            deallocate(vec_matrix)
            deallocate(rij2_vector)
            deallocate(pij_vector)
            deallocate(qij_vector)
            deallocate(factor)
        end do
        !$omp end parallel do

        gradient_vector = reshape(gradient_matrix,(/low_dimension*reduced_number_points/))

        call gradient_vector_addnoise(gradient_vector, 1e-2_dp)
        
        
    end subroutine loss_gradient_core


    function calculate_stepsize(current_position, current_gradient, init) result(step_size)

        logical, optional, intent(in)                                 :: init
        real(kind=dp), dimension(:), intent(in)                       :: current_position, current_gradient
        real(kind=dp), save, allocatable, dimension(:)                :: previous_position, previous_gradient
        real(kind=dp), dimension(size(current_position))              :: position_change, gradient_change
        real(kind=dp)                                                 :: gradient_change_magnitude, position_gradient_dot_product
        real(kind=dp)                                                 :: step_size
        real(kind=dp), parameter                                      :: default_step_size=1e-1_dp

        if(init) then
            if (allocated(previous_position)) deallocate(previous_position)
            if (allocated(previous_gradient)) deallocate(previous_gradient)
            allocate(previous_position(size(current_position)))
            allocate(previous_gradient(size(current_gradient)))
            step_size=default_step_size
        
        else
        
            position_change=current_position- previous_position
            gradient_change=current_gradient- previous_gradient

            gradient_change_magnitude=dot_product(gradient_change,gradient_change)
            position_gradient_dot_product=dot_product(position_change,gradient_change)
        
            step_size=abs(position_gradient_dot_product/gradient_change_magnitude)
        
        end if
        
        previous_position=current_position
        previous_gradient=current_gradient

    end function calculate_stepsize

    function calculating_cost_zero (pij) result(cost_zero)

        implicit none
        real(kind=dp), dimension(reduced_number_points,reduced_number_points), intent(in) :: pij
        real(kind=dp) :: cost_zero
        integer:: i, j

        cost_zero=0.0_dp

        !$omp parallel do reduction(+:cost_zero) collapse(2)
        do i=1,reduced_number_points
            do j=i+1,reduced_number_points
                if (pij(j,i)>0 .and. pij(j,i)<1) then
                    cost_zero = cost_zero + pij(j,i) * log(pij(j,i)) * 2.0_dp + (1 - pij(j,i)) * log(1 - pij(j,i)) * 2.0_dp
                else if (pij(i,j)>1) then
                    cost_zero = cost_zero + pij(j,i)*log(pij(j,i))*2.0_dp
                else if (pij(i,j)<0) then
                    cost_zero = cost_zero + (1-pij(j,i))*log(1-pij(j,i))*2.0_dp
                end if
            end do
        end do
        !$omp end parallel do

    end function calculating_cost_zero

    subroutine gradient_vector_addnoise(gradient_vector, r) 

        real(kind=dp)                              :: d
        real(kind=dp),               intent(in)    :: r
        real(kind=dp), dimension(:), intent(in) :: gradient_vector
        real(kind=dp), dimension(size(gradient_vector)) :: noise_vector
        real(kind=dp), dimension(size(gradient_vector)) :: gradient_vector_noise

        noise_vector=0.0_dp
        gradient_vector_noise=0.0_dp

        call random_number(d)
        
        d= r*d**(1/real(size(gradient_vector),dp))

        call random_add_noise(noise_vector,1.0_dp)

        noise_vector = d*noise_vector/norm2(noise_vector)

        gradient_vector_noise = gradient_vector*(1.0_dp+noise_vector)

    end subroutine gradient_vector_addnoise

    subroutine handle_growth_phase(i, growth_start_step, growth_step_limit)
        implicit none

        integer, intent(in)       :: i
        integer, intent(in)       :: growth_start_step
        logical, intent(inout)    :: growth_step_limit

        if ((i > growth_start_step) .and. (i < growth_start_step + growth_steps)) then
            if (i < growth_start_step + 2) then
                point_radius = point_radius * ((real(i)-real(growth_start_step))/real(growth_steps))
            else
                point_radius = point_radius * ((real(i)-real(growth_start_step))*1.0_dp/real(growth_steps))&
                /((real(i)-real(growth_start_step)-1.0_dp)*1.0_dp/real(growth_steps))
            end if
        end if

        if(i > (growth_start_step+growth_steps+100)) growth_step_limit=.false.

    end subroutine handle_growth_phase
    
    subroutine loss_gradient(gradient_matrix, cost, exageration, growth_switch)
        implicit none
        logical,                                                intent(in)  :: growth_switch 
        real(kind=dp),                                          intent(in)  :: exageration
        real(kind=dp),                                          intent(out) :: cost
        real(kind=dp), dimension(low_dimension, reduced_number_points), intent(out) :: gradient_matrix

        real(kind=dp), dimension(low_dimension)               :: vec, pos
        real(kind=dp)                                         :: z, rij2, qij, dist
        integer                                               :: i, j
        
        z = real(reduced_number_points, dp) * (real(reduced_number_points, dp)-1.0_dp)

        cost = cost_zero
        gradient_matrix=0.0_dp

        !$omp parallel do private(pos, rij2, qij, vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
        do i=1,reduced_number_points
            do j=i+1,reduced_number_points
                pos(:)=low_dimension_position(:,i)-low_dimension_position(:,j)
                qij= calculating_qij(i,j)                
                vec(:) = 4.0_dp * z * (exageration*pij(j,i) - (1-pij(j,i)) / (1-qij) * qij) * qij * pos(:)
                gradient_matrix(:,i)=gradient_matrix(:,i)+vec(:)
                gradient_matrix(:,j)=gradient_matrix(:,j)-vec(:)
                cost=cost-pij(j,i)*log(qij)*2.0_dp-(1-pij(j,i))*log(1-qij)*2.0_dp    
            end do
        end do
        !$omp end parallel do
                
        if (growth_switch) then
            !$omp parallel do private(pos, rij2, dist, vec) reduction(+:gradient_matrix,cost) schedule(dynamic)
            do i=1,reduced_number_points
                do j=i+1,reduced_number_points
                    pos(:)=low_dimension_position(:,i)-low_dimension_position(:,j)
                    rij2=dot_product(pos,pos)
                    dist=sqrt(rij2)
                    if(dist < point_radius(i)+point_radius(j)) then
                        vec(:)=-pos/dist
                        dist=(point_radius(i)+point_radius(j)-dist)/2.0_dp
                        gradient_matrix(:,i)=gradient_matrix(:,i)+vec(:)*dist*core_strength/2.0_dp
                        gradient_matrix(:,j)=gradient_matrix(:,j)-vec(:)*dist*core_strength/2.0_dp
                        cost=cost+dist**2/2.0_dp*core_strength
                    end if
                end do
            end do
            !$omp end parallel do
            write (*,*) 'Growth phase, cost: ', cost
        end if


    end subroutine loss_gradient



end MODULE optimisation
