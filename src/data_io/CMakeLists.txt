add_library(dataio_lib 
            STATIC
            data_reader.f90
            data_writer.f90)

target_link_libraries(dataio_lib
                      PRIVATE
                      constants_lib
                      solver_lib 
                      ${NETCDF_LIBRARIES})

target_include_directories(dataio_lib
                        PUBLIC
                        ${NETCDF_INCLUDE_DIRS}
                        ${CMAKE_BINARY_DIR}/src/constants
                        ${CMAKE_BINARY_DIR}/src/solver
)