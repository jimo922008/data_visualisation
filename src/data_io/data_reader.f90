!> @brief Read data from a file with .vec extension.
!> @param[in] filename The name of the file to read.
!> @param[out] data_vec The data vector.
!> @param[out] ion_energy The ion energy.
!> @param[out] ion_volume The ion volume.
!> @param[out] ion_label The ion label.
!> @param[out] ion_symmetry The ion symmetry.
!> @param[out] ion_formula The ion formula.
!> @param[out] ion_num The ion number.
!> @param[out] point_count The point count.
!> @param[out] number_points The number of points.
!> @param[out] number_features The number of features.
!> @details This module reads data from a file. The file format is as follows:
!! The first line contains the number of features. The next lines contain the data vector, the ion label, the ion number, the ion formula, the ion symmetry, the ion volume, the ion energy, and the point count.
!! The number of features is the number of elements in the data vector. The number of points is the number of rows in the file divided by 3.
!! The data vector is a 2D array with the number of features as the first dimension and the number of points as the second dimension.
!! The ion label, ion formula, ion symmetry, ion volume, ion energy, and point count are 1D arrays with the number of points as the first dimension.
!! The ion volume and ion energy are divided by the ion number. The number of points and the number of features are public variables.
!! The data vector is initialized to zero. The data is read from the file and stored in the arrays.
!! The file is closed after reading the data. The data is printed to the console. If the data vector is allocated, the number of points is printed. Otherwise, an error message is printed.

MODULE data_reader

   USE parameters

   IMPLICIT NONE

   public :: read_file

contains

   subroutine read_file(filename, data)

      !> @brief Read data from a file with .vec extension.
      !> @param[in] filename The name of the file to read.
      !> @param[out] data The file_data structure, please refer to parameters module.
      !> @details This module reads data from a file. The file format is as follows:
      !! The first line contains the number of features. The next line contains the data vector, and the third line contains the ion label, the ion number, the ion formula,
      !! the ion symmetry, the ion volume, the ion energy, and the point count.

      ! * Declare variables
      character(len=*), intent(in)                :: filename
      type(file_data), intent(out)                :: data

      ! * Declare local variables
      integer                       :: unit, stat, col, number_of_rows
      character(len=256)            :: iomsg
      character(len=:), allocatable :: line

      ! * Open the file
      open (newunit=unit, file=filename, status='old', action='read', iostat=stat, iomsg=iomsg)

      print *, 'Reading file: ', filename

      ! * Count the number of points.

      data%number_features = 0

      read (unit, *, iostat=stat) data%number_features
      if (stat /= 0) then
         print *, 'Error reading the number of features'
         print *, 'IOSTAT value:', stat
         print *, 'IOMSG value:', iomsg
         stop
      end if

      number_of_rows = 1
      data%number_points = 0

      allocate (character(len=data%number_features) :: line)

      do
         read (unit, *, iostat=stat, iomsg=iomsg) line
         if (stat /= 0) then
            exit
         else
            number_of_rows = number_of_rows + 1
         end if
      end do

      data%number_points = number_of_rows/3

      rewind (unit)

      allocate (data%data_vec(data%number_features, data%number_points))
      allocate (data%ion_energy(data%number_points))
      allocate (data%ion_volume(data%number_points))
      allocate (data%ion_label(data%number_points))
      allocate (data%ion_symmetry(data%number_points))
      allocate (data%ion_formula(data%number_points))
      allocate (data%ion_num(data%number_points))
      allocate (data%point_count(data%number_points))

      ! * Read the data into arrays.
      do col = 1, data%number_points

         ! * Read the column length for each row
         read (unit, *) data%number_features

         ! * Read the data vector
         data%data_vec(:, col) = 0.0_sp
         read (unit, *) data%data_vec(1:data%number_features, col)

         ! * Read meta-data vector
         read (unit, *) data%ion_label(col), data%ion_num(col), data%ion_formula(col), &
            data%ion_symmetry(col), data%ion_volume(col), &
            data%ion_energy(col), data%point_count(col)

         data%ion_volume(col) = data%ion_volume(col)/real(data%ion_num(col), sp)
         data%ion_energy(col) = data%ion_energy(col)/real(data%ion_num(col), sp)

      end do

      close (unit)

      if (allocated(data%data_vec)) then
         print *, 'Data read from file:'

         if (size(data%data_vec, 1) > 0) then
            print *, data%number_points, 'points read.'
         else
            print *, 'Data array is empty or dimensions are zero. Programme stopping now.'
            stop
         end if

      else
         print *, 'No data read or allocation failed. Programme stopping now.'
         stop
      end if

   end subroutine read_file

end module data_reader
