MODULE low_dimension_probability

    USE constants
    USE data_reader
    USE initialisation
    USE high_dimension

    IMPLICIT NONE
    PUBLIC :: low_dimension_distribution
    PUBLIC :: calculating_low_dimension_distance
    PUBLIC :: calculating_qij
    
    integer :: low_dimension  
    real(kind=dp), allocatable, dimension(:,:)  :: low_dimension_position(:,:)
    real(kind=dp), allocatable, dimension(:,:)  :: low_dimension_distance(:,:)
    real(kind=dp), allocatable, dimension(:,:)  :: qij(:,:)

contains

    subroutine low_dimension_distribution()

        call low_dimension_init()
        call calculating_low_dimension_distance()
        call calculating_qij()

    end subroutine low_dimension_distribution

    subroutine low_dimension_init()

        IMPLICIT NONE

        real(kind=dp), allocatable :: random_matrix(:,:)
        real(kind=dp)              :: coeff, sigma_average
        integer :: i, j

        write (*,*) 'Please insert the dimension of the low dimensional data.'
        read (*,*) low_dimension

        allocate (low_dimension_position(low_dimension, number_points))
        allocate (low_dimension_distance(number_points, number_points))
        allocate (qij(number_points, number_points))
        allocate (random_matrix(low_dimension, number_features))
        
        do i = 1, low_dimension
            call random_add_noise(random_matrix(i,:), 1.0_dp)
            random_matrix(i,:) = random_matrix(i,:) / sqrt(sum(random_matrix(i,:)**2))
            if (i>1) then
                do j = 1, i-1
                    coeff = dot_product(random_matrix(i-j,:), random_matrix(i,:))
                    random_matrix(i,:) = random_matrix(i,:) - coeff * random_matrix(i-j,:)
                end do
            end if
        end do 
        write (*,*) 'random matrix done.'
        
        sigma_average = sum(sigma) / number_points

        low_dimension_position = matmul(random_matrix, data_vec) / sigma_average

        write (*,*) 'low dimension position done.'

        deallocate(random_matrix)
    
    end subroutine low_dimension_init

    subroutine calculating_low_dimension_distance()

        integer :: i, j

        low_dimension_distance = 0.0_dp

        do i = 1, number_points
            do j = i+1, number_points
                low_dimension_distance(i,j) = sum((low_dimension_position(:,i) - low_dimension_position(:,j))**2)
                low_dimension_distance(j,i) = low_dimension_distance(i,j)
            end do
        end do  

        write (*,*) 'low dimension distance done.'

    end subroutine calculating_low_dimension_distance

    subroutine calculating_qij() 
        implicit none
        integer :: i, j

        qij = 0.0_dp
        do i = 1, number_points
            do j = i+1, number_points
                qij(i,j)= 1/((1.0_dp+ low_dimension_distance(i,j))*number_points*(number_points-1))
                qij(j,i) = qij(i,j)
            end do   
        end do
    
    end subroutine calculating_qij

    subroutine random_add_noise(vec,sigma)

        integer                                    :: n

        real(kind=dp)                              :: U,V
        real(kind=dp), intent(in)                  :: sigma
        real(kind=dp), dimension(:), intent(inout) :: vec

        call random_seed()
        do n=1,size(vec)
            
            call random_number(U)
            call random_number(V)

            vec(n)=vec(n)+sqrt(-2*log(U))*cos(tpi*V)*sigma

        end do

    end subroutine random_add_noise

end module low_dimension_probability