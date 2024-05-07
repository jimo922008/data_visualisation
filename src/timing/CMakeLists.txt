add_library(timing_lib 
            STATIC
            timing.f90
            )

target_link_libraries(solver_lib
            PRIVATE
            constants_lib
            ${NETCDF_LIBRARIES})

target_include_directories(solver_lib
                 PUBLIC
                 ${NETCDF_INCLUDE_DIRS}
                 ${CMAKE_BINARY_DIR}/src/constants
)