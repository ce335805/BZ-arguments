#specefiy required version - this is the newest one
cmake_minimum_required(VERSION 3.11)
set (CMAKE_CXX_EXTENSIONS OFF)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_COMPILER "mpiicpc")
set (CMAKE_C_COMPILER "mpiicc")

#set as inlcude dir the current directory
#one could devide the project into more directorys later
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(project_sources
  victory.cpp
  vertex.cpp
	Index.cpp
	output.inc
	FF.cpp
	util.cpp
        )

#include gtest
include_directories("/home/chris/gtest/googletest/googletest/include")
#include MKL
include_directories("/opt/intel/mkl/include")

#add compiler flags inlcuding cblas
#last 4 options for hdf5 files
set(C_flags "-O3 -Wall -mkl -Wextra -pedantic -lsz -lz -ldl -lm")
set(CMAKE_CXX_FLAGS	"${CMAKE_CXX_FLAGS} ${C_flags}")

#set project name
project (victory)
#add an executable - using variables
add_executable(${PROJECT_NAME} ${project_sources})

#link to gtest
target_link_libraries(${PROJECT_NAME} "/home/chris/gtest/build/lib/libgtest.a")
target_link_libraries(${PROJECT_NAME} "/home/chris/gtest/build/lib/libgtest_main.a")
#place where libfabric is installed - needed for IntelMPI-2019
target_link_libraries(${PROJECT_NAME} "/usr/local/lib/libfabric.so.1")
#link to pthreads
target_link_libraries(${PROJECT_NAME} "pthread")
#link to hdf5_libraries
target_link_libraries(${PROJECT_NAME} "/usr/lib/libhdf5_hl_cpp.a")
target_link_libraries(${PROJECT_NAME} "/usr/lib/libhdf5_cpp.a")
target_link_libraries(${PROJECT_NAME} "/usr/lib/libhdf5_hl.a")
target_link_libraries(${PROJECT_NAME} "/usr/lib/libhdf5.a")
