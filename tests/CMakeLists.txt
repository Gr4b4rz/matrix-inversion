cmake_minimum_required (VERSION 2.8.9)
project (matrix-inversion-tests)

enable_testing()
set (CMAKE_CXX_STANDARD 17)
set(SOURCES ../main.cpp ../Matrix.h)
set(CMAKE_CXX_FLAGS "-lpthread -floop-parallelize-all -ftree-parallelize-loops=4")
find_package(Boost COMPONENTS unit_test_framework system  REQUIRED)
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK)
add_executable(tests main_test.cpp)
target_link_libraries(tests ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_SYSTEM_LIBRARY})

ADD_TEST (tests tests)
