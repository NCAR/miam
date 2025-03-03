cmake_minimum_required(VERSION 3.21)

project(miam VERSION 0.0.0 LANGUAGES CXX)

# Add custom CMake modules to the list where CMake searches for modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)

################################################################################
# Projet wide setup options and variables

option(MIAM_ENABLE_TESTS "Build the tests" ON)
option(MIAM_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(MIAM_ENABLE_COVERAGE "Enable code coverage output" OFF)
option(MIAM_BUILD_DOCS "Build the documentation" OFF)

set(MIAM_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(MIAM_LIB_DIR ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

include(dependencies)

################################################################################
add_subdirectory(src)

if(MIAM_ENABLE_TESTS)
  if(MIAM_ENABLE_COVERAGE)
    include(CodeCoverage)
    append_coverage_compiler_flags()
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY ${PROJECT_SOURCE_DIR}/src
    )
  endif()

  enable_testing()
  add_subdirectory(test)
endif()

if(MIAM_BUILD_DOCS)
  add_subdirectory(docs)
endif()

# Packaging

if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()