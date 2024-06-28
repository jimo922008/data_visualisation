message("* Adding build type Custom...")

SET(CMAKE_Fortran_FLAGS_CUSTOM
    "-O3 -g -march=native -fimplicit-none -Wall -Wextra -fopenmp "
    CACHE STRING "Flags used by the Fortran compiler during Custom builds."
    FORCE )

# Be careful if useing two custom build types this will overwrite as it does not append.
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "Custom")

set(CMAKE_Fortran_FLAGS_CUSTOM "${CMAKE_Fortran_FLAGS_CUSTOM}" CACHE STRING "Flags used by the Fortran compiler during custom build type." FORCE)