add_library(solver_lib 
            STATIC
            high_dimension.f90
            low_dimension_probability.f90
            initialisation.f90
            optimisation.f90
            )

target_link_libraries(solver_lib
            PRIVATE
            parameters_lib
            datain_lib
            timing_lib
            )

target_include_directories(solver_lib
                 PUBLIC
                 ${CMAKE_BINARY_DIR}/src/parameters
                 ${CMAKE_BINARY_DIR}/src/data_io
                 ${CMAKE_BINARY_DIR}/src/timing
)