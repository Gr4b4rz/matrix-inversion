cmake_minimum_required (VERSION 2.8.9)
project (matrix-inversion-tests)

enable_testing()
set (CMAKE_CXX_STANDARD 17)

find_program(MPIC /usr/lib64/openmpi/bin/mpicc)
if (NOT MPIC)
    message(FATAL_ERROR "mpicc not found! Forgot to install mpi package?")
endif()
set(MPI_C_COMPILER /usr/lib64/openmpi/bin/mpicc)

find_program(MPIC++ /usr/lib64/openmpi/bin/mpic++)
if (NOT MPIC++)
    message(FATAL_ERROR "mpi++ not found! Forgot to install mpi package?")
endif()
set(MPI_CXX_COMPILER /usr/lib64/openmpi/bin/mpic++)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(SOURCES ../main.cpp ../Matrix.h)
set(CMAKE_CXX_FLAGS "-lpthread -floop-parallelize-all -ftree-parallelize-loops=4")
find_package(Boost COMPONENTS unit_test_framework system  REQUIRED)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
add_executable(tests main_test.cpp)
target_link_libraries(tests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_SYSTEM_LIBRARY})
target_link_libraries(tests MPI::MPI_CXX)

ADD_TEST (tests tests)
