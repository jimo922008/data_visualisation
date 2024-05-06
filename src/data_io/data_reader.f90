MODULE data_io

    USE constants
    
    IMPLICIT NONE

    public :: read_vec, write_xyz

    type data_vector
        real(kind=dp), dimension(:,:), allocatable  :: data_vec
        real(kind=dp), dimension(:), allocatable    :: ion_energy, ion_volume
        integer, dimension (:), allocatable         :: ion_num, point_count
        character(len=*), dimension(:), allocatable :: ion_label, ion_symmetry, ion_formula
        type(data_vector), pointer                  :: next => null()
    end type data_vector

    type(data_vector) :: vectors
        
    integer, public, allocatable, dimension(:)      :: number_points, number_features

contains 

    subroutine read_vec(filename)

        character(len=*) intent(in) :: filename
        integer :: i, unit, stat,row, col, number_of_rows

        ! * Open the file.
        open(newunit=unit, file=filename, status='old', action='read')

        ! * Count the number of rows.
        number_points=0
        number_features=0

        do
            read(unit, *, iostat= stat) i
            if (stat /= 0) exit 
            number_of_rows = number_of_rows +1
        end do 

        number_points = number_of_rows/3

        rewind(unit) 
        
        ALLOCATE (
            data_vec(number_points, number_features),
            ion_energy(number_points),
            ion_volume(number_points),
            ion_label(number_points),
            ion_symmetry(number_points),
            ion_formula(number_points),
            ion_num(number_points),
            point_count(number_points)
        )

        ! * Read the data into arrays. 
        do row=1, number_points

            ! * Read the column length for each row
            read(stdin,*) number_features

            ! * Read the data vector 
            data_vec(:,row)= 0.0_dp
            read(stdin,*) data_vec(1:col_count,row)
            row = row + 1

            ! * Read meta-data vector
            read(stdin,*) ion_label(row),ion_num(row),ion_formula(row),ion_symmetry(row),ion_volume(row),ion_energy(row),point_count(row)

            ion_volume(row)=ion_volume(row)/real(ion_num(row),dp)
            ion_energy(row)=ion_energy(row)/real(ion_num(row),dp)

        end do

        allocate(
            point_pos(col_count,row_count),
            point_radius(row_count),
            dij2(row_count,row_count)
        )

        total_point_count=sum(point_count(1:row_count))

    end subroutine read_vec
    
end module data_reader