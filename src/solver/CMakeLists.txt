add_library(solver_lib 
            STATIC
            high_dimension.f90
            low_dimension_distribution.f90
            initialisation.f90
            optimisation.f90
            )

target_link_libraries(solver_lib
            PRIVATE
            constants_lib
            dataio_lib
            ${NETCDF_LIBRARIES})

target_include_directories(solver_lib
                 PUBLIC
                 ${NETCDF_INCLUDE_DIRS}
                 ${CMAKE_BINARY_DIR}/src/constants
                ${CMAKE_BINARY_DIR}/src/data_io
)