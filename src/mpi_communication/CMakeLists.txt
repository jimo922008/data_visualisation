add_library(datasync_lib 
            STATIC
            mpi_model.f90
            )

target_link_libraries(dataio_lib
            PRIVATE
            constants_lib
            ${NETCDF_LIBRARIES})

target_include_directories(dataio_lib
                 PUBLIC
                 ${NETCDF_INCLUDE_DIRS}
                 ${CMAKE_BINARY_DIR}/src/constants
)