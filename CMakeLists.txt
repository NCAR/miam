cmake_minimum_required(VERSION 3.24)

project(miam VERSION 0.0.0 
        LANGUAGES CXX
        DESCRIPTION "Model-Independent Aerosol Module"
        HOMEPAGE_URL "https://github.com/NCAR/miam")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add custom CMake modules to the list where CMake searches for modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

# Force Debug for coverage
if(MIAM_ENABLE_COVERAGE AND NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "" FORCE)
endif()

include(GNUInstallDirs)

################################################################################
# Project wide setup options and variables

option(MIAM_ENABLE_TESTS "Build the test suite" ON)
option(MIAM_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(MIAM_ENABLE_COVERAGE "Enable code coverage output" OFF)
option(MIAM_BUILD_DOCS "Build the documentation" OFF)

set(MIAM_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(MIAM_LIB_DIR ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

################################################################################
# Set options for the third party libraries that MIAM needs

set(MICM_DEFAULT_VECTOR_SIZE "4" CACHE STRING "Set MICM vector-ordered matrix dimension")

include(dependencies)  # finds cmake/dependencies.cmake

################################################################################
# MIAM library

add_library(miam)
add_library(musica::miam ALIAS miam)

add_subdirectory(src)

target_compile_features(miam PUBLIC cxx_std_20)

target_include_directories(miam
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${MIAM_INSTALL_INCLUDE_DIR}>
)

target_link_libraries(miam
  PUBLIC
    musica::micm
)

set_target_properties(miam PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

if(MIAM_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)

  if(MIAM_ENABLE_COVERAGE)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(miam PRIVATE --coverage)
    target_link_options(miam PRIVATE --coverage)

    include(CodeCoverage)
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE ctest
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY ${PROJECT_SOURCE_DIR}/src
    )
    else()
      message(FATAL_ERROR "Coverage is only supported with GCC or Clang")
    endif()
  endif()
endif()

if(MIAM_BUILD_DOCS)
  add_subdirectory(docs)
endif()

