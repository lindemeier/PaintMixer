cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
include(ExternalProject)

project(3rdParty)

# ceres variables
set(BUILD_TESTING OFF)
set(BUILD_EXAMPLES OFF)
set(BUILD_SHARED_LIBS OFF)

if(WITH_TESTS)
  add_subdirectory(Catch2)  
endif(WITH_TESTS)

add_subdirectory(glog)
add_subdirectory(ceres-solver)
