module test_read_file
   use pFUnit
   use data_io
   implicit none

contains

   @test
   subroutine test_read_file()
      integer :: iostat
      character(len=256) :: filename
      type(data_vector), allocatable :: expected_data(:)

      ! Setup
      filename = 'LJ13.vec'  ! Assume this test data file is pre-populated with known values

      ! Execute
      call read_file(filename)

      @assertEqual(number_points, 10000, "Number of points mismatch")

   end subroutine test_read_file

end module test_read_file
