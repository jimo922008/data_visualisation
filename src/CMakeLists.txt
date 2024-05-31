add_subdirectory(parameters)
add_subdirectory(data_io)
add_subdirectory(solver)
add_subdirectory(timing)


find_package(MPI REQUIRED)
#find_package(OpenMP REQUIRED)

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
    #PUBLIC OpenMP::OpenMP_Fortran
    ${BLAS_LIBRARIES}
)


target_include_directories(sheap
                 PUBLIC
                 ${CMAKE_BINARY_DIR}/src/parameters
                ${CMAKE_BINARY_DIR}/src/data_io
                ${CMAKE_BINARY_DIR}/src/solver
                ${CMAKE_BINARY_DIR}/src/timing
                
)