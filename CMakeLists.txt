cmake_minimum_required(VERSION 3.15)

project(Sheap
        VERSION 1.0 
        DESCRIPTION "Sheap project using Fortran" 
        LANGUAGES Fortran)

# Compiler options for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fexternal-blas -march=native -ffast-math -funroll-loops")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

# Setup for the executable and linking the libraries
add_executable(sheap
  src/data_parsing.f90
  src/initialisation.f90
  src/parameters.f90
)

# Link libraries
target_link_libraries(
    sheap 
    PUBLIC 
    MPI::MPI_Fortran
    OpenMP::OpenMP_Fortran
)

install(TARGETS sheap DESTINATION bin)



