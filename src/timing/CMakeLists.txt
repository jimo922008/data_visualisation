add_library(timing_lib 
            STATIC
            timing.f90
            )

target_link_libraries(timing_lib
            PRIVATE
            parameters_lib
            )

target_include_directories(timing_lib
                 PUBLIC
                 ${CMAKE_BINARY_DIR}/src/parameters
)