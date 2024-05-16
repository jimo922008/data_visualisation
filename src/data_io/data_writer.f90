MODULE data_writer
    
    USE constants

    IMPLICIT NONE 

contains
    subroutine write_file (low_dimension_position)
        ! Declare variables
        real(dp), dimension(:,:), intent(in) :: low_dimension_position
        integer :: unit_number
        integer :: i,j

        unit_number = 20
        open(unit=unit_number, file='output.txt', status='replace', action='write')

      ! Write array to file
        do i = 1, size(low_dimension_position, 1)
            ! Write each column of the current row to file
            do j = 1, size(low_dimension_position, 2)
                write(unit_number, '(F8.2, 1X)', advance='no') low_dimension_position(i, j)
            end do
            write(unit_number, *) ! New line after each row
        end do

        ! Close the file
        close(unit_number)
    end subroutine write_file

END MODULE data_writer