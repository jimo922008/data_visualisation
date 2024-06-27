cmake_minimum_required(VERSION 3.15)

project(sheap
        VERSION 1.0 
        DESCRIPTION "Sheap project using Fortran" 
        LANGUAGES Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -fexternal-blas -march=native -ffast-math -funroll-loops")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

# Compiler options for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fexternal-blas -march=native -ffast-math -funroll-loops ")

# set(CMAKE_BUILD_TYPE "DEBUG")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

include(CTest)
enable_testing()

add_subdirectory(src)
# add_subdirectory(test)
