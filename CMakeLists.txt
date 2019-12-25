cmake_minimum_required (VERSION 3.14)
project (matrix-inversion)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)

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
set(CMAKE_CXX_FLAGS "-lpthread -floop-parallelize-all -ftree-parallelize-loops=4")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(SOURCES main.cpp Matrix.h)
add_executable(matrix-inversion main.cpp)
target_link_libraries(matrix-inversion MPI::MPI_CXX)

