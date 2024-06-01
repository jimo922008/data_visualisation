cmake_minimum_required(VERSION 3.15)

project(sheap
        VERSION 1.0 
        DESCRIPTION "Sheap project using Fortran" 
        LANGUAGES Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -g -O0")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")

# Compiler options for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fexternal-blas -march=native -ffast-math -funroll-loops")

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
find_package(BLAS REQUIRED)

if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_Fortran_LINK_FLAGS}")
    set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_Fortran_LIBRARIES})
endif()

include(CTest)
enable_testing()

add_subdirectory(src)
add_subdirectory(test)


