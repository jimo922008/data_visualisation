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
    real(kind=dp), allocatable, dimension(:)    :: point_radius(:)

    real(kind=dp), public     :: projection_compression =10.0_dp
    real(kind=dp), public     :: packing_fraction=0.5_dp

    integer, public           :: trace_step=1
    integer, public           :: niter=10000
    integer, public           :: nexag=200
    integer, public           :: nneg=1

    real(kind=dp), public     :: sexag=5.0_dp
    real(kind=dp), public     :: alpha=-1.0_dp
    real(kind=dp), public     :: tolerance=1e-8_dp

    ! * Hard sphere growth    
    integer, public           :: growth_steps=500
    real(kind=dp), public     :: sphere_radius=0.01_dp
    real(kind=dp), public     :: core_strength=1.0_dp
    real(kind=dp), public     :: growth_tol_coeff=2.0_dp

contains

    subroutine low_dimension_distribution()

        call low_dimension_init()
        call low_dimension_normalisation()
        call calculating_low_dimension_distance()
        call calculating_qij()
        call hard_sphere_initialisation()

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
            random_matrix(i,:) = random_matrix(i,:) /norm2(random_matrix(i,:))
            if (i>1) then
                do j = 1, i-1
                    coeff = dot_product(random_matrix(i-j,:), random_matrix(i,:))
                    random_matrix(i,:) = random_matrix(i,:) - coeff * random_matrix(i-j,:)
                end do
            end if
            random_matrix(i,:) = random_matrix(i,:) / norm2(random_matrix(i,:))
        end do 
        
        sigma_average = sum(sigma) / real(reduced_number_points,dp)

        low_dimension_position = matmul(random_matrix, data_vec) / sigma_average / projection_compression 

        write (*,*) 'low dimension position done.'

        deallocate(random_matrix)
    
    end subroutine low_dimension_init

    subroutine low_dimension_normalisation()

        implicit none

        integer :: i
        real(kind=dp), allocatable, dimension(:) :: mass_centre

        allocate (mass_centre(low_dimension))

        do i = 1, number_points
            if (point_count(i) == 0) cycle
            mass_centre = mass_centre + low_dimension_position(:,i)
        end do

        mass_centre = mass_centre / number_points

        do i = 1, number_points
            if (point_count(i) == 0) cycle
            low_dimension_position(:,i) = low_dimension_position(:,i) - mass_centre
        end do
    
    end subroutine low_dimension_normalisation

    subroutine calculating_low_dimension_distance()

        integer :: i, j
        real(kind=dp), dimension(low_dimension) :: vec

        low_dimension_distance = 0.0_dp

        !$omp parallel do collapse(2)
        do i = 1, number_points
            do j = i+1, number_points
                if ((point_count(i) /= 0).and.(point_count(j) /= 0)) then
                    vec = low_dimension_position(:,i) - low_dimension_position(:,j)
                    low_dimension_distance(j,i) = dot_product(vec,vec)
                    low_dimension_distance(i,j) = low_dimension_distance(j,i)
                end if
            end do
        end do  

    end subroutine calculating_low_dimension_distance

    subroutine calculating_qij() 
        implicit none
        integer :: i, j

        qij = 0.0_dp

        !$omp parallel do collapse(2)
        do i = 1, number_points
            do j = i+1, number_points
                if ((point_count(i) /= 0).and.(point_count(j) /= 0)) then
                    qij(j,i)= 1/((1.0_dp+ low_dimension_distance(j,i))*number_points*(number_points-1))
                    qij(i,j) = qij(j,i)
                end if
            end do   
        end do
        !$omp end parallel do
    
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

    subroutine hard_sphere_initialisation()
        implicit none

        integer :: i

        real(kind=dp) :: radius

        allocate (point_radius(number_points))

        point_radius(i) = 0.0_dp

        !$omp parallel do schedule(static)
        do i=1,number_points
            if(point_count(i)==0) cycle
            point_radius(i)=sphere_radius*real(point_count(i),dp)**(1.0_dp/real(low_dimension,dp))
        end do
        !$omp end parallel do
    
    end subroutine hard_sphere_initialisation


end module low_dimension_probability