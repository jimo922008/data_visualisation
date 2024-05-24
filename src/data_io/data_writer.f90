MODULE data_writer
    
    USE constants
    USE data_reader
    USE high_dimension
    USE low_dimension_probability
    USE optimisation

    IMPLICIT NONE 

contains
    subroutine write_file (filename)
        character(len=*) :: filename
        character(len=20):: fmt
        integer :: unit_number
        integer :: i,j

        unit_number = 20
        open(unit=unit_number, file=trim(filename), status='replace', action='write')

        write (unit_number, '(i0)') reduced_number_points
        write (unit_number,'(a,i0,a,f16.10)') 'dim ',low_dimension,' : cost', final_cost
        do i=1,reduced_number_points
            write(fmt,*) low_dimension
            write(unit_number,'(a2,'//trim(adjustl(fmt))//'g25.13,a50,i6,a15,a10,g25.13,g25.13,i6,g25.13)') 'H',(low_dimension_position(j,i),j=1,low_dimension),&
            trim(ion_label(i)),ion_num(i),trim(ion_formula(i)),trim(ion_symmetry(i)),ion_volume(i),ion_energy(i),&
            point_count(i),max(sphere_radius,point_radius(i))
        end do

        flush(unit_number)
        close(unit_number)


    end subroutine write_file

END MODULE data_writer