MODULE data_io

    USE constants
    
    IMPLICIT NONE

    public :: read_vec
    public :: write_xyz
    integer, public, allocatable, dimension(:)      :: row_count, col_count

contains 

    subroutine read_vec()

        integer :: i, unit, stat
        integer :: row, col

        type data_vector
            real(kind=dp), dimension(:,:), allocatable  :: data_vec
            real(kind=dp), dimension(:), allocatable    :: ion_energy, ion_volume
            integer, dimension (:), allocatable         :: ion_num, point_count
            character(len=*), dimension(:), allocatable :: ion_label, ion_symmetry, ion_formula
            type(data_vector), pointer                  :: next => null()
        end type data_vector
        
        ! * Count the number of rows.
        row_count=0
        col_count=0

        unit = 10
        do
            read(unit, *, iostat= stat) i
            if (stat /= 0) exit 
            number_of_rows = number_of_rows +1
        end do 

        row_count = number_of_rows/3

        rewind(unit) 
        
        ALLOCATE (
            data_vec(col_count, row_count),
            ion_energy(row_count),
            ion_volume(row_count),
            ion_label(row_count),
            ion_symmetry(row_count),
            ion_formula(row_count),
            ion_num(row_count),
            point_count(row_count)
        )

        ! * Read the data into arrays. 
        do row=1, row_count

            ! * Read the column length for each row
            read(stdin,*) col_count

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

    subroutine write_xyz (cost)
        integer                   :: ni,i
        character(len=20)         :: fmt
        real(kind=dp), intent(in) :: cost

        write(stdout,'(i0)') reduced_point_count
        write(stdout,'(a,i0,a,f16.10)') 'dim ',ldim,' : cost', cost

        do ni=1,npoints
            if(point_count(ni).eq.0) cycle
                write(fmt,*) ldim
                write(stdout,'(a2,'//trim(adjustl(fmt))//'g25.13,a50,i6,a15,a10,g25.13,g25.13,i6,g25.13)') 'H',(point_pos(i,ni),i=1,ldim),&
                trim(ion_label(ni)),ion_num(ni),trim(ion_formula(ni)),trim(ion_symmetry(ni)),ion_volume(ni),ion_energy(ni),&
                point_count(ni),max(sphere_radius,point_radius(ni))
        end do

        flush(stdout)

    end subroutine model_write
    
end module data_io