MODULE data_reader

    USE constants
    
    IMPLICIT NONE

    public :: read_file

    real(kind=dp), dimension(:,:), allocatable  :: data_vec
    real(kind=dp), dimension(:), allocatable    :: ion_energy, ion_volume
    integer, dimension (:), allocatable         :: ion_num, point_count
    character(len= 256), dimension(:), allocatable :: ion_label, ion_symmetry, ion_formula
        
    integer, public                             :: number_points, number_features

contains 

    subroutine read_file(filename)

        character(len=*), intent(in) :: filename
        integer :: i, unit, stat, row, col, number_of_rows

        ! * Open the file.
        open(newunit=unit, file=filename, status='old', action='read')

        ! * Count the number of points.

        number_features=0

        read(unit, *, iostat= stat) number_features
        if (stat /= 0) then
            print *, 'Error reading the number of features'
            stop
        end if

        number_of_rows  = 1
        number_points = 0

        do
            read(unit, *, iostat= stat) i
            if (stat /= 0) exit 
            number_of_rows = number_of_rows +1
        end do 

        number_points = number_of_rows/3

        rewind(unit) 
        
        allocate(data_vec(number_features, number_points), &
        ion_energy(number_points), ion_volume(number_points), &
        ion_label(number_points), ion_symmetry(number_points), &
        ion_formula(number_points), ion_num(number_points), &
        point_count(number_points))

        ! * Read the data into arrays. 
        do col=1, number_points

            ! * Read the column length for each row
            read(unit,*) number_features

            ! * Read the data vector 
            data_vec(:,col)= 0.0_dp
            read(unit,*) data_vec(1:number_features,col)

            ! * Read meta-data vector
            read(unit,*) ion_label(col),ion_num(col),ion_formula(col), &
            ion_symmetry(col),ion_volume(col), &
            ion_energy(col),point_count(col)

            ion_volume(col)=ion_volume(col)/real(ion_num(col),dp)
            ion_energy(col)=ion_energy(col)/real(ion_num(col),dp)

        end do

        close(unit)

    end subroutine read_file

end module data_reader