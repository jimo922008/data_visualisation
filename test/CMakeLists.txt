include(CTest)
enable_testing()

add_executable(test_optimisation test_optimisation.f90)

target_link_libraries(
    test_optimisation
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


target_include_directories(test_optimisation
                 PUBLIC
                ${CMAKE_BINARY_DIR}/src/parameters
                ${CMAKE_BINARY_DIR}/src/data_io
                ${CMAKE_BINARY_DIR}/src/solver
                ${CMAKE_BINARY_DIR}/src/timing
                
)

add_test(NAME test_optimisation  COMMAND test_optimisation)