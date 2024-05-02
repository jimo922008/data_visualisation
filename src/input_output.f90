module input_output

  ! * Prerequisite modules
  use parameters

  ! * All variables must be declared
  implicit none

  ! * All variables private unless declared public
  private

  ! * Public subroutines
  public :: read_vec
  public :: write_xyz

  ! * Public variables
  integer, public, allocatable, dimension(:)      :: n_points, n_points_max

  contains 

    subroutine read_vec()

        integer :: ni,length,min_length
        integer                       :: indx,stat,met2,met7
        character(len=50)             :: met1,met3,met4
        real(kind=dp)                 :: met5,met6

        min_length=huge(1)

        npoints=0
        do ni=1,npoints_max

        ! * Read length of vector, update minimum length if necessary

        read(stdin,*,end=1001) length
        if(length.lt.min_length) min_length=length
        if(ni.eq.1) then

        allocate(
            data_vec(min_length,npoints_max),
            ion_energy(npoints_max),
            ion_volume(npoints_max),
            ion_label(npoints_max),&
            ion_symmetry(npoints_max),
            ion_formula(npoints_max),
            ion_num(npoints_max),
            point_count(npoints_max))
        end if

        ! * Read descriptor vector
        data_vec(:,ni)=0.0_dp
        read(stdin,*) data_vec(1:min_length,ni)
        npoints=npoints+1

        ! * Read meta-data from .vec file

        read(stdin,*) ion_label(ni),ion_num(ni),ion_formula(ni),ion_symmetry(ni),ion_volume(ni),ion_energy(ni),point_count(ni)

        end if
        ion_volume(ni)=ion_volume(ni)/real(ion_num(ni),dp)
        ion_energy(ni)=ion_energy(ni)/real(ion_num(ni),dp)
        end do
        1001 continue

        if(npoints.eq.npoints_max) then
        write(stderr,*) 'ERROR: number of points in .vec file exceeds maximum number allocated. Increase this with -np.'
        end if

        hdim=min_length
        allocate(point_pos(ldim,npoints),point_radius(npoints),dij2(npoints,npoints))

        total_point_count=sum(point_count(1:npoints))

    end subroutine read_vec

    subroutine write_xyz (cost)
        integer                   :: ni,i
        character(len=20)         :: fmt
        real(kind=dp), intent(in) :: cost

        !! * For paper: computes #contacts with distances descriptor
!       call compute_ion_nn

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
    
    subroutine compute_ion_nn
  
        integer :: ni,nj
    
        if(.not.allocated(ion_nn)) allocate(ion_nn(npoints))
    
        ion_nn_dist=2.3
    
        ion_nn=0
        do ni=1,npoints
            if(point_count(ni).eq.0) cycle
            do nj=1,hdim
                if(data_vec(nj,ni).lt.ion_nn_dist) then
                    ion_nn(ni)=ion_nn(ni)+1
                else
                    exit
                end if
            end do
        ion_nn(ni)=ion_nn(ni)/real(ion_num(ni),dp)
        end do
  
    end subroutine compute_ion_nn

end module input_output