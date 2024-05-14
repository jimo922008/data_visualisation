add_subdirectory(constants)
add_subdirectory(data_io)
add_subdirectory(solver)
add_subdirectory(timing)
add_subdirectory(mpi_communication)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

add_executable(sheap main_hpc.f90)

target_link_libraries(
    sheap
    PRIVATE 
    constants_lib
    dataio_lib
    solver_lib
    timing_lib
    PUBLIC MPI::MPI_Fortran
    PUBLIC OpenMP::OpenMP_Fortran
)


target_include_directories(sheap
                 PUBLIC
                 ${NETCDF_INCLUDE_DIRS}
                 ${CMAKE_BINARY_DIR}/src/constants
                ${CMAKE_BINARY_DIR}/src/data_io
                ${CMAKE_BINARY_DIR}/src/solver
)