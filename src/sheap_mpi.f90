program sheap_mpi

include 'mpif.h'

integer rank, nranks, ierr, tag

call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)


    


call MPI_Finalize(ierr)

end program sheap_mpi