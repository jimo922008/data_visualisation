MODULE data_io

    USE constants
    
    IMPLICIT NONE   
    
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