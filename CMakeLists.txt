cmake_minimum_required(VERSION 3.15)

project(sheap
        VERSION 1.0 
        DESCRIPTION "Sheap project using Fortran" 
        LANGUAGES Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
enable_testing()

# Compiler options for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fexternal-blas -march=native -ffast-math -funroll-loops")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(src)
add_subdirectory(test)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

install(TARGETS sheap DESTINATION bin)



