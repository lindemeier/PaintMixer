cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
include(ExternalProject)

project(PaintMixer)

option(WITH_TESTS "Whether tests should be built and run" ON)

# 3rdParty dependencies
find_package(OpenCV 3.4 REQUIRED core imgproc imgcodecs highgui)
find_package(OpenMP REQUIRED)
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
add_subdirectory(3rdParty)



# build the main target
add_library(${PROJECT_NAME}Lib STATIC
  src/PaintCoeff.hxx
  src/PaintCoeff.cxx

  src/Palette.hxx
  src/Palette.cxx

  src/PaintMixer.hxx
  src/PaintMixer.cxx
)

target_include_directories(${PROJECT_NAME}Lib 
  PUBLIC ${Eigen3_INCLUDE_DIRS}
  PRIVATE ${OpenCV_INCLUDE_DIRS}
)

target_link_libraries(${PROJECT_NAME}Lib PRIVATE
  ceres
  glog
  ${OpenCV_LIBRARIES}
  OpenMP::OpenMP_CXX
)

set_target_properties(${PROJECT_NAME}Lib PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME}Lib PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}Lib PROPERTIES CXX_STANDARD_REQUIRED ON)

# build the main target
add_executable(${PROJECT_NAME} 
  src/main.cpp
)

target_link_libraries(${PROJECT_NAME}  PRIVATE
  ${PROJECT_NAME}Lib 
)

set_target_properties(${PROJECT_NAME}  PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME}  PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}  PROPERTIES CXX_STANDARD_REQUIRED ON)