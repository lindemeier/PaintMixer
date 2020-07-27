cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
include(ExternalProject)

project(PaintMixer)

option(WITH_TESTS "Whether tests should be built and run" ON)

# 3rdParty dependencies
find_package(OpenCV REQUIRED core imgproc imgcodecs highgui)
find_package(OpenMP REQUIRED)
find_package(PkgConfig)
pkg_search_module(Eigen3 REQUIRED eigen3)
add_subdirectory(3rdParty)
set (ColorConverter_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/3rdParty/ColorConverter/include)
set (cereal_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/3rdParty/cereal/include)
set (catch2_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/3rdParty/Catch2/single_include)
set (cxxopts_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/3rdParty/cxxopts/include)
# ceres solver library
find_package(Ceres REQUIRED)

# build the main target
add_library(${PROJECT_NAME}Lib STATIC
  src/include/PaintMixer/PaintCoeff.hxx
  src/PaintCoeff.cxx

  src/include/PaintMixer/KubelkaMunk.hxx

  src/include/PaintMixer/Palette.hxx
  src/Palette.cxx

  src/include/PaintMixer/Serialization.hxx
  src/Serialization.cxx

  src/include/PaintMixer/PaintMixer.hxx
  src/PaintMixer.cxx

  src/ColorExtraction.hxx
  src/ColorExtraction.cxx
)

target_include_directories(${PROJECT_NAME}Lib
  PUBLIC src/include
  PUBLIC ${Eigen3_INCLUDE_DIRS}
  PUBLIC ${OpenCV_INCLUDE_DIRS}
  PRIVATE ${ColorConverter_INCLUDE_DIRS}
  PUBLIC ${cereal_INCLUDE_DIRS}
  PUBLIC $<BUILD_INTERFACE:${CERES_INCLUDE_DIRS}>
)

target_link_libraries(${PROJECT_NAME}Lib PRIVATE
  ${OpenCV_LIBRARIES}
  OpenMP::OpenMP_CXX
  ${CERES_LIBRARIES}
)

set_target_properties(${PROJECT_NAME}Lib PROPERTIES COMPILE_FLAGS "-Wall -Weverything    -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME}Lib PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}Lib PROPERTIES CXX_STANDARD_REQUIRED ON)

# build the main target
add_executable(${PROJECT_NAME}
  src/main.cpp
)

target_include_directories(${PROJECT_NAME}
  PRIVATE $<BUILD_INTERFACE:${cxxopts_INCLUDE_DIRS}>
)

target_link_libraries(${PROJECT_NAME}  PRIVATE
  ${PROJECT_NAME}Lib
)

set_target_properties(${PROJECT_NAME}  PROPERTIES COMPILE_FLAGS "-Wall -Weverything -Wno-c++98-compat -Wno-padded -Wno-documentation")
set_target_properties(${PROJECT_NAME}  PROPERTIES CXX_STANDARD 17)
set_target_properties(${PROJECT_NAME}  PROPERTIES CXX_STANDARD_REQUIRED ON)

# tests
if (WITH_TESTS)
  add_executable(${PROJECT_NAME}_UnitTests test/main.cxx)

  target_include_directories(${PROJECT_NAME}_UnitTests
    PUBLIC ${catch2_INCLUDE_DIRS}
    )

  target_link_libraries(${PROJECT_NAME}_UnitTests
    PRIVATE ${PROJECT_NAME}Lib)

endif(WITH_TESTS)
