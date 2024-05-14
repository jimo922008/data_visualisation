MODULE high_dimension

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

        write (stderr, *) 'calculating distance'
        do i = 1, number_points
            do j = i + 1, number_points
                high_dist_matrix(i,j) = sqrt(length_2(i) + length_2(j) - 2.0_dp * dot_product(data_vec(:, i), data_vec(:, j)))
                high_dist_matrix(j,i) = high_dist_matrix(i,j)  
            end do
            high_dist_matrix(i,i) = 0.0_dp  
        end do
        
        deallocate(length_2)  

    end subroutine high_dimension_distance

    subroutine find_sigma(perplexity)
        implicit none
    
        real(dp), intent(in) :: perplexity
        real(dp) :: low_sigma, high_sigma, mid_sigma
        real(dp) :: low_perplexity, high_perplexity, mid_perplexity
        real(dp) :: tolerance
        integer  :: i, iter, max_iter

        allocate(sigma(number_points))

        tolerance = 1.0e-5_dp
        max_iter = 1000

        do i = 1, number_points

            low_sigma = 0.0001_dp
            high_sigma = 1.0_dp

            iter = 0

            low_perplexity = calculating_perplexity(low_sigma, i)
            high_perplexity = calculating_perplexity(high_sigma, i)

            do 
                if (iter >= max_iter .or. abs(high_sigma - low_sigma) <= tolerance) exit

                mid_sigma = (low_sigma + high_sigma) / 2.0_dp
                mid_perplexity = calculating_perplexity(mid_sigma, i)

                if (abs(mid_perplexity - perplexity) < tolerance) then
                    sigma(i) = mid_sigma
                    exit
                end if

                if (mid_perplexity > perplexity) then
                    high_sigma = mid_sigma
                else
                    low_sigma = mid_sigma
                end if

                iter = iter + 1 
            end do

            sigma(i) = mid_sigma
        end do

    end subroutine find_sigma
    
    function calculating_perplexity(sigma, i) result(perplexity)
        implicit none
        real(dp), intent(in) :: sigma
        integer, intent (in) :: i
        integer              :: j
        real(dp), allocatable, dimension(:) :: pji_conditional
        real(dp)             :: entropy, perplexity

        entropy = 0.0_dp
        perplexity = 0.0_dp

        pji_conditional = calculating_pji(sigma,i)

        do j = 1, number_points
            if (i == j) then
                cycle
            end if
            if (pji_conditional(j) == 0.0_dp) then
                cycle
            end if
            entropy = entropy + pji_conditional(j) * log(pji_conditional(j))
        end do

        perplexity = 2.0_dp ** (-entropy / log(2.0_dp))

    end function calculating_perplexity

    function calculating_pji(sigma,i) result(pji_conditional)
        implicit none
        real(kind=dp), allocatable, dimension(:) :: pji_conditional
        real(kind=dp), intent(in) :: sigma
        integer, intent (in) :: i
        integer              :: j

        allocate(pji_conditional(number_points))

        do j = 1, number_points
            pji_conditional(j) = exp(-high_dist_matrix(i,j)/(sigma*sigma*2.0_dp))
            if (i == j) then
                pji_conditional(j) = 0.0_dp
            end if
        end do

        pji_conditional = pji_conditional / sum(pji_conditional)

    end function calculating_pji

    function calculating_pij(sigma,j) result(pij_conditional)
        implicit none
        real(dp), allocatable, dimension(:) :: pij_conditional
        real(kind=dp), intent(in) :: sigma
        integer, intent (in) :: j
        integer              :: i

        allocate(pij_conditional(number_points))

        do i = 1, number_points
            pij_conditional(i) = exp(-high_dist_matrix(i,j)/(sigma*sigma*2.0_dp))
            if (i == j) then
                pij_conditional(i) = 0.0_dp
            end if
        end do
    end function calculating_pij

    subroutine high_dimension_distribution()
        implicit none
        integer :: i, j
        real(kind=dp) :: pij_1
        real(kind=dp) :: pji_1
        
        write (stderr, *) 'calculating high dimension distance......'

        call high_dimension_distance()

        write (stderr, *) 'finding sigma......'

        call find_sigma(perplexity)

        write (stderr, *) 'calculating pij.......'

        allocate(pij(number_points, number_points))

        do i = 1, number_points
            do j = i+1, number_points
                pij_1 = exp(-high_dist_matrix(i,j)/(sigma(j)*sigma(j)*2.0_dp))
                pji_1 = exp(-high_dist_matrix(i,j)/(sigma(i)*sigma(i)*2.0_dp))
                pij(i,j) = (pij_1+pji_1)/(2.0_dp*number_points)
                pij(j,i) = pij(i,j)
            end do
            pij(i,i) = 0.0_dp
        end do   
    end subroutine high_dimension_distribution

END MODULE high_dimension