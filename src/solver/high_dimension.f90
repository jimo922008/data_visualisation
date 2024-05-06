MODULE high_dimension_distribution

    USE constants
    USE data_io
    USE initialisation

    IMPLICIT NONE

    PRIVATE

    PUBLIC :: 

    subroutine high_dimension_distance()

        integer :: i, j
        real(dp), allocatable :: length_2(:)
        real(dp), allocatable :: high_dimension_distance(:,:)

        allocate(
            length_2(row_count),
            high_dimension_distance(row_count, row_count))

        do i = 1, row_count
            length_2(i) = dot_product(data_vec(1:col_count, i), data_vec(1:col_count, i))
        end do

        write (stderr, *) 'calculating distance'
        do i = 1, row_count
            do j = i + 1, row_count
                distance(i,j) = sqrt(length2(i) + length2(j) - 2.0_dp * dot_product(data_vec(1:col_count, i), data_vec(1:col_count, j)))
                distance(j,i) = distance(i,j)  
            end do
            distance(i,i) = 0.0_dp  
        end do

    end subroutine high_dimension_distance

    subroutine find_sigma(i, sigma, perplexity)
        integer, intent(in) :: i
        real(dp), intent(in) :: perplexity
        real(dp), intent(out) :: sigma

        real(dp) :: low_sigma, high_sigma, mid_sigma
        real(dp) :: low_perplexity, high_perplexity, mid_perplexity
        real(dp) :: tolerance
        integer :: iter, max_iter

        low_sigma = 1.0_dp
        high_sigma = 10.0_dp
        tolerance = 1.0e-5_dp
        max_iter = 1000

        low_perplexity = calculating_perplexity(low_sigma, i)
        high_perplexity = calculating_perplexity(high_sigma, i)

        iter = 0
        do while (iter < max_iter) .and. (abs(high_sigma - low_sigma) > tolerance)
            mid_sigma = (low_sigma + high_sigma) / 2.0_dp
            mid_perplexity = calculating_perplexity(mid_sigma, i)

            if (abs(mid_perplexity - perplexity) < tolerance) then
                sigma = mid_sigma
                return
            end if

            if (mid_perplexity > perplexity) then
                high_sigma = mid_sigma
            else
                low_sigma = mid_sigma
            end if

            iter = iter + 1
        end do
        if (abs(perplexity - low_perplexity) < tolerance) then
            sigma = low_sigma
            return
        end if
    
    real(dp) function calculating_perplexity(i, sigma) 
        implicit none
        integer, intent :: i
        real(dp), intent(in) :: sigma
        real(dp) :: perplexity

        perplexity = 2.0_dp ** (sum(-calculating_pji(i, sigma) * log(calculating_pji(i, sigma))) / log(2.0_dp))
    end function


    real(dp) function calculating_pji(i, sigma)
            
        integer, intent (in) :: i
        integer :: j
        real(dp), intent (in) :: sigma
        real(dp), allocatable :: pji(:)

        allocate(pji(row_count))

        write (stderr, *) 'calculating pji'

        do j = 1, row_count
            pji(j) = exp(-distance(i,j)/(sigma^2*2.0_dp))
        end do

        pji(:) = pji(:) / sum(pji(:))

    end function calculating_pji

END MODULE model