include(FetchContent)

FetchContent_Declare(
    PFUNIT
    GIT_REPOSITORY https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git
    GIT_TAG        main
)

FetchContent_MakeAvailable(PFUNIT)

enable_testing()

find_package(MPI REQUIRED)

set(test_sources 
    test_data_reader.pf)

add_pfunit_ctest(
    sheap_test
    TEST_SOURCES ${test_sources}
    LINK_LIBRARIES 
    dataio_lib
    EXTRA_USE mpi_pfunit_mod
)