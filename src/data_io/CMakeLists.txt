add_library(datain_lib 
            STATIC
            data_reader.f90)

target_link_libraries(datain_lib
                      PRIVATE
                      parameters_lib    
                      ${NETCDF_LIBRARIES})

target_include_directories(datain_lib
                        PUBLIC
                        ${NETCDF_INCLUDE_DIRS}
                        ${CMAKE_BINARY_DIR}/src/parameters
)

add_library(dataout_lib 
            STATIC
            data_writer.f90)

target_link_libraries(dataout_lib
                      PRIVATE
                      parameters_lib
                      solver_lib
                      ${NETCDF_LIBRARIES})

target_include_directories(dataout_lib
                        PUBLIC
                        ${NETCDF_INCLUDE_DIRS}
                        ${CMAKE_BINARY_DIR}/src/parameters
                        ${CMAKE_BINARY_DIR}/src/solver
                       
)