################################################################################
# Preamble
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
# Dependencies

include(dependencies)

################################################################################
# MIAM library and docs

add_subdirectory(src)

if(MIAM_BUILD_DOCS)
  add_subdirectory(docs)
endif()


################################################################################
# Tests

if(PROJECT_IS_TOP_LEVEL AND MIAM_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)

  if(MIAM_ENABLE_COVERAGE)
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE ctest
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY ${PROJECT_SOURCE_DIR}/src
        LCOV_ARGS --ignore-errors mismatch --ignore-errors unused
    )
  endif()
endif()

################################################################################
# Packaging

# only include packaging if we are the top level project being built
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()