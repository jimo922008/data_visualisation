MODULE mpi_optimisation
    
    use optimisation
    use constants


    IMPLICIT NONE

    integer                                    :: i, j
    integer                                    :: pair_index, start_index, end_index, number_pairs
    integer, dimension(:,:), allocatable       :: pairs, pairs_per_rank
    real(kind=dp)                              :: cost_local, cost_global
    real(kind=dp), dimension(:,:), allocatable :: gradient_matrix_local, gradient_matrix_global

CONTAINS

    SUBROUTINE loss_gradient_mpi()

        implicit none

        integer :: number_points, i, j

        do pair_index = start_index, end_index
            i = pairs(1, pair_index)
            j = pairs(2, pair_index)
            
        end do


    END  SUBROUTINE loss_gradient_mpi

    SUBROUTINE counting_pairs (number_points)

        implicit none

        integer :: number_points, i, j
        
        number_pairs = (number_points * (number_points - 1)) / 2
        
        allocate(pairs(2, number_pairs))

        pair_index = 0

        do i = 1, number_points
            do j = i + 1, number_points
                pair_index = pair_index + 1
                pairs(1, pair_index) = i
                pairs(2, pair_index) = j
            end do
        end do

        pairs_per_rank = number_pairs / nranks
        start_index = rank * pairs_per_rank + 1

        if (rank == nranks-1) then
            end_index = number_pairs
        else
            end_index = start_index + pairs_per_rank - 1
        end if

    END SUBROUTINE counting_pairs

END MODULE mpi_optimisation