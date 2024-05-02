cmake_minimum_required(VERSION 3.15)

project(Sheap
        VERSION 1.0 
        DESCRIPTION "Sheap project using Fortran" 
        LANGUAGES Fortran)

# Compiler options for different build types
set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -fimplicit-none -Wall -Wextra")
set(CMAKE_Fortran_FLAGS_RELEASE "-Ofast -fexternal-blas -march=native -ffast-math -funroll-loops")

find_package(OpenMP REQUIRED)
find_package(MKL REQUIRED)
find_package(MPI REQUIRED)

if(MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
endif()

# Setup for the executable and linking the libraries
add_executable(sheap
  src/sheap.f90
  src/optim.f90
  src/loss.f90
  src/model.f90
  src/control.f90
  src/random.f90
  src/rng.f90
  src/constants.f90
)

# Link libraries
if(OpenMP_CXX_FOUND)
  target_link_libraries(]
    sheap 
    PUBLIC 
    OpenMP::OpenMP_CXX)
endif()

if(MKL_FOUND)
  target_link_libraries(
    sheap 
    PUBLIC 
    ${MKL_LIBRARIES})

  target_link_options(
    sheap 
    PUBLIC 
    ${MKL_LINK_FLAGS})
endif()

# Additional libraries
find_library(
    SYMS_LIB 
    NAMES symspg 
    PATHS ${PREFIX}/lib)

target_link_libraries(
    sheap 
    ${SYMS_LIB})

# Installation rules
install(
    TARGETS sheap 
    DESTINATION bin)

# Custom targets mimicking Makefile behaviour
add_custom_target(
    neat COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/cleanup.cmake)

add_custom_command(
  OUTPUT ${CMAKE_BINARY_DIR}/cleanup.cmake
  COMMAND echo "file(REMOVE_RECURSE \${OBJ} *.mod)" > ${CMAKE_BINARY_DIR}/cleanup.cmake
)
