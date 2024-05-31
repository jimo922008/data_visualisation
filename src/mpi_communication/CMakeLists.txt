add_library(datasync_lib 
            STATIC
            mpi_model.f90
            )

target_link_libraries(datasync_lib
            PRIVATE
            dataio_lib
            parameters_lib
            solver_lib
            ${NETCDF_LIBRARIES})

target_include_directories(datasync_lib
            PUBLIC
            ${NETCDF_INCLUDE_DIRS}
            ${CMAKE_BINARY_DIR}/src/parameters
            ${CMAKE_BINARY_DIR}/src/data_io
            ${CMAKE_BINARY_DIR}/src/solver
)