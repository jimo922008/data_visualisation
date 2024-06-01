MODULE low_dimension_probability

   USE parameters
   USE data_reader
   USE initialisation
   USE high_dimension

   IMPLICIT NONE
   PUBLIC :: low_dimension_distribution
   PUBLIC :: calculating_low_dimension_distance
   PUBLIC :: calculating_qij
   PUBLIC :: random_add_noise

   real(kind=sp), allocatable, dimension(:, :)  :: low_dimension_position(:, :)
   real(kind=sp), allocatable, dimension(:, :)  :: low_dimension_distance(:, :)
   real(kind=sp), allocatable, dimension(:, :)  :: qij(:, :)
   real(kind=sp), allocatable, dimension(:)     :: point_radius(:)

contains

   subroutine low_dimension_distribution()

      call low_dimension_init()
      write (*, *) 'low dimension init done.'
      call low_dimension_normalisation()
      write (*, *) 'low dimension normalisation done.'
      call hard_sphere_initialisation()
      write (*, *) 'hard sphere initialisation done.'

   end subroutine low_dimension_distribution

   subroutine low_dimension_init()

      IMPLICIT NONE

      real(kind=sp), allocatable :: random_matrix(:, :)
      real(kind=sp)              :: coeff, sigma_average
      integer :: i, j

      allocate (low_dimension_position(low_dimension, reduced_number_points))
      allocate (low_dimension_distance(reduced_number_points, reduced_number_points))
      allocate (qij(reduced_number_points, reduced_number_points))
      allocate (random_matrix(low_dimension, number_features))

      random_matrix = 0.0_sp

      do i = 1, low_dimension
         call random_add_noise(random_matrix(i, :), 1.0_sp)
         random_matrix(i, :) = random_matrix(i, :)/norm2(random_matrix(i, :))
         if (i > 1) then
            do j = 1, i - 1
               coeff = dot_product(random_matrix(i - j, :), random_matrix(i, :))
               random_matrix(i, :) = random_matrix(i, :) - coeff*random_matrix(i - j, :)
            end do
         end if
         random_matrix(i, :) = random_matrix(i, :)/norm2(random_matrix(i, :))
      end do

      sigma_average = sum(sigma)/real(reduced_number_points, sp)

      low_dimension_position = matmul(random_matrix, data_clean)/sigma_average/projection_compression

      write (*, *) 'low dimension position done.'

      deallocate (random_matrix)

   end subroutine low_dimension_init

   subroutine low_dimension_normalisation()

      real(kind=sp), allocatable, dimension(:) :: mass_centre

      allocate (mass_centre(low_dimension))

      mass_centre = sum(low_dimension_position, dim=2)/real(reduced_number_points, sp)

      low_dimension_position = low_dimension_position - spread(mass_centre, 2, reduced_number_points)

   end subroutine low_dimension_normalisation

   subroutine calculating_low_dimension_distance()

      integer :: i, j
      real(kind=sp), dimension(low_dimension) :: vec

      low_dimension_distance = 0.0_sp

      !$omp parallel do schedule(static)
      do i = 1, reduced_number_points
         do j = i + 1, reduced_number_points
            vec(:) = low_dimension_position(:, i) - low_dimension_position(:, j)
            low_dimension_distance(j, i) = dot_product(vec, vec)
            low_dimension_distance(i, j) = low_dimension_distance(j, i)
         end do
      end do
      !$omp end parallel do

   end subroutine calculating_low_dimension_distance

   function calculating_qij(i, j) result(qij)

      implicit none
      integer, intent(in) :: i, j
      real(kind=sp)                      :: qij
      real(kind=sp)                      :: rij2, z
      real(kind=sp), dimension(low_dimension) :: pos

      z = real(reduced_number_points, sp)*(real(reduced_number_points, sp) - 1.0_sp)

      pos(:) = low_dimension_position(:, i) - low_dimension_position(:, j)
      rij2 = dot_product(pos, pos)
      qij = 1.0_sp/(1.0_sp + rij2)/z

   end function calculating_qij

   subroutine random_add_noise(vec, variance)

      real(kind=sp), intent(in)                  :: variance
      real(kind=sp), dimension(:), intent(inout) :: vec
      real(kind=sp), dimension(size(vec))        :: u, v

      call random_number(u)
      call random_number(v)

      vec = vec + sqrt(-2.0_sp*log(u))*cos(2.0_sp*acos(-1.0_sp)*v)*variance

   end subroutine random_add_noise

   subroutine hard_sphere_initialisation()
      implicit none

      integer :: i

      allocate (point_radius(reduced_number_points))

      point_radius = 0.0_sp

      !$omp parallel do schedule(static)
      do i = 1, reduced_number_points
         point_radius(i) = sphere_radius*real(point_count_clean(i), sp)**(1.0_sp/real(low_dimension, sp))
      end do
      !$omp end parallel do

      point_radius = point_radius/100.0_sp

   end subroutine hard_sphere_initialisation

end module low_dimension_probability
