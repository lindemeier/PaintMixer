cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
include(ExternalProject)

project(PaintMixer)

option(WITH_TESTS "Whether tests should be built and run" ON)

# 3rdParty dependencies
find_package(OpenMP REQUIRED)
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
add_subdirectory(3rdParty)


# build the main target
add_executable(${PROJECT_NAME} 
  src/main.cpp

  src/PaintCoeff.hxx
  src/PaintCoeff.cxx

  src/Palette.hxx
  src/Palette.cxx

  src/PaintMixer.hxx
  src/PaintMixer.cxx
)

target_include_directories(${PROJECT_NAME} 
  PRIVATE ${Eigen3_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
  ceres
  glog
  OpenMP::OpenMP_CXX
)

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD_REQUIRED ON)