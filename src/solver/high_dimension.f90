MODULE high_dimension

    USE omp_lib
    USE constants
    USE data_reader
    USE initialisation

    IMPLICIT NONE

    PUBLIC :: high_dimension_distance
    PUBLIC :: find_sigma
    PUBLIC :: calculating_perplexity
    PUBLIC :: calculating_pji
    
    real(dp):: perplexity, perplexity_target
    real(dp), allocatable, dimension(:,:)   :: high_dist_matrix(:,:)
    real(dp), allocatable, dimension(:)     :: sigma(:)
    real(dp), allocatable, dimension(:,:)   :: pij(:,:)
    

contains


    subroutine high_dimension_distance()
        implicit none
        integer :: i, j
        real(dp), allocatable, dimension(:)   :: length_2(:)

        allocate(length_2(number_points), &
        high_dist_matrix(number_points, number_points))

        do i = 1, number_points
            length_2(i) = dot_product(data_vec(:, i), data_vec(:, i))
        end do

        !$omp parallel do private(j) schedule(dynamic)
        do i = 1, number_points
            do j = 1, number_points
                high_dist_matrix(j,i) = length_2(i) + length_2(j) + 1e-10_dp
            end do
            high_dist_matrix (i,i) = 0.0_dp  
        end do
        !$omp end parallel do

        call dgemm('T','N',number_points,number_points,number_features,-2.0_dp,data_vec(1:number_features,1:number_points),&
               number_features,data_vec(1:number_features,1:number_points),number_features,1.0_dp,high_dist_matrix,number_points)

        ! _ single precision version??? 
        ! 

        do i =1, number_points
            high_dist_matrix(i,i) = 0.0_dp
        end do
        
        deallocate(length_2)  

    end subroutine high_dimension_distance

    subroutine find_sigma(perplexity)
        implicit none
    
        real(dp), intent(in) :: perplexity
        real(dp) :: low_sigma, high_sigma, mid_sigma
        real(dp) :: low_perplexity, high_perplexity, mid_perplexity
        real(dp) :: tolerance, factor
        integer  :: i

        allocate(sigma(number_points))

        tolerance = 1.0e-13_dp

        !omp parallel do private(i, high_sigma, high_perplexity, low_sigma, low_perplexity, factor, mid_sigma, mid_perplexity) schedule(dynamic)
        do i = 1, number_points
            
            if (point_count(i) .eq. 0) cycle
            
            factor = gr

            high_sigma = 1.0_dp
            high_perplexity = calculating_perplexity(high_sigma, i)

            if (high_perplexity > perplexity) then
                factor = 1.0_dp/factor
            end if
            
            low_sigma = high_sigma * factor
            low_perplexity = calculating_perplexity(low_sigma, i)

            do while ((low_perplexity - perplexity)*(high_perplexity - perplexity) > 0.0_dp)
                high_sigma = low_sigma 
                high_perplexity = low_perplexity
                low_sigma = low_sigma * factor
                low_perplexity = calculating_perplexity(low_sigma, i)      
            end do
           
            mid_perplexity = huge(1.0_dp)
            do while ((abs(high_sigma - low_sigma) > tolerance) .and. (abs(mid_perplexity - perplexity) > tolerance))

                mid_sigma = (low_sigma + high_sigma) / 2.0_dp

                if ((mid_sigma == low_sigma) .or. (mid_sigma == high_sigma)) exit

                mid_perplexity = calculating_perplexity(mid_sigma, i)

                if ((mid_perplexity- perplexity)*(high_perplexity-perplexity)>0) then
                    high_sigma = mid_sigma
                    high_perplexity = mid_perplexity
                else
                    low_sigma = mid_sigma
                    low_perplexity = mid_perplexity
                end if

            end do
            sigma(i) = mid_sigma

        end do
        !omp end parallel do

    end subroutine find_sigma
    
    function calculating_perplexity(sigma, i) result(perplexity)
        implicit none
        real(dp), intent(in) :: sigma
        integer, intent (in) :: i
        integer              :: j
        real(dp)             :: pji_conditional(number_points)
        real(dp)             :: entropy, perplexity

        entropy = 0.0_dp
        perplexity = 0.0_dp
        pji_conditional = 0.0_dp

        pji_conditional = calculating_pji(sigma,i)

        !omp parallel do private(j) reduction(+:entropy)
        do j = 1, number_points
            if ((i == j) .or. (point_count(j) == 0.0_dp) .or. (pji_conditional(j)<tiny(1.0_dp))) cycle
            entropy = entropy - pji_conditional(j) * log(pji_conditional(j))/log(2.0_dp)
        end do
        !omp end parallel do

        perplexity = 2.0_dp**(entropy)

    end function calculating_perplexity

    function calculating_pji(sigma,i) result(pji_conditional)
        implicit none
        real(kind=dp)             :: pji_conditional(number_points)
        real(kind=dp), intent(in) :: sigma
        integer, intent (in) :: i
        integer              :: j

        !omp parallel do private(j)
        do j = 1, number_points
            if ((i == j).or.(point_count(j) == 0)) cycle
            pji_conditional(j) = exp(-high_dist_matrix(j,i)/(sigma*sigma*2.0_dp))
        end do
        !omp end parallel do

        pji_conditional(:) = pji_conditional(:) / sum(pji_conditional)

    end function calculating_pji

    subroutine high_dimension_distribution()
        implicit none
        integer :: i, j
        real(kind=dp) :: pi_j(number_points, number_points)
        
        write (stderr, *) 'calculating high dimension distance......'

        call high_dimension_distance()

        call remove_duplicates(data_vec, high_dist_matrix, number_points, 0.5_dp, .10000E+09_dp)

        write (stderr, *) 'finding sigma......'

        call find_sigma(perplexity)

        write (stderr, *) 'calculating pij.......'

        allocate(pij(number_points, number_points))
        
        pi_j = 0.0_dp
        !$omp parallel do schedule(dynamic)
        do i = 1, number_points
            if (point_count(i) ==0) cycle
            do j = 1, number_points
                if ((point_count(j) ==0) .or. (i == j)) cycle
                pi_j(j,i) = exp(-high_dist_matrix(j,i)/(sigma(i)**2*2.0_dp))
            end do
            pi_j(:,i) = pi_j(:,i) / sum(pi_j(:,i))
        end do   
        !$omp end parallel do

        !$omp parallel do schedule(dynamic)
        do i = 1, number_points
            if (point_count(i) ==0) cycle
            do j = i+1, number_points
                if (point_count(j) ==0) cycle
                pij(i,j) = (pi_j(i,j)+pi_j(j,i))/real(2.0_dp*reduced_number_points)
                pij(j,i) = pij(i,j)
            end do
        end do
        !$omp end parallel do
       

    end subroutine high_dimension_distribution

END MODULE high_dimension