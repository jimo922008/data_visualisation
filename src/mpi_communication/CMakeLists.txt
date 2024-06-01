find_package(MPI REQUIRED)

add_library(datasync_lib 
            STATIC
            mpi_model.f90
            )

target_link_libraries(datasync_lib
            PRIVATE
            datain_lib
            parameters_lib
            solver_lib
            timing_lib
            PUBLIC MPI::MPI_Fortran
            PUBLIC OpenMP::OpenMP_Fortran
            ${BLAS_LIBRARIES}
)

target_include_directories(datasync_lib
            PUBLIC
            ${CMAKE_BINARY_DIR}/src/parameters
            ${CMAKE_BINARY_DIR}/src/data_io
            ${CMAKE_BINARY_DIR}/src/solver
            ${CMAKE_BINARY_DIR}/src/timing
)