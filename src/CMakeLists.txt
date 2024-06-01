add_subdirectory(parameters)
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
    parameters_lib
    datain_lib
    dataout_lib
    solver_lib
    timing_lib
    
    PUBLIC MPI::MPI_Fortran
    PUBLIC OpenMP::OpenMP_Fortran
    ${BLAS_LIBRARIES}
)


target_include_directories(sheap
                 PUBLIC
                 ${CMAKE_BINARY_DIR}/src/parameters
                ${CMAKE_BINARY_DIR}/src/data_io
                ${CMAKE_BINARY_DIR}/src/solver
                ${CMAKE_BINARY_DIR}/src/timing
                
)

add_executable(sheap_fast sheap_mpi.f90)

target_link_libraries(
    sheap_fast
    PRIVATE 
    parameters_lib
    datain_lib
    dataout_lib
    solver_lib
    timing_lib
    datasync_lib
    
    PUBLIC MPI::MPI_Fortran
    PUBLIC OpenMP::OpenMP_Fortran
    ${BLAS_LIBRARIES}
)


target_include_directories(sheap_fast
                 PUBLIC
                ${CMAKE_BINARY_DIR}/src/parameters
                ${CMAKE_BINARY_DIR}/src/data_io
                ${CMAKE_BINARY_DIR}/src/solver
                ${CMAKE_BINARY_DIR}/src/timing
                ${CMAKE_BINARY_DIR}/src/mpi_communication
                
)