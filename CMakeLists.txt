cmake_minimum_required(VERSION 3.21)

project(miam VERSION 0.0.0 LANGUAGES CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

include(GNUInstallDirs)

################################################################################
# Projet wide setup options

option(MIAM_ENABLE_TESTS "Build the tests" ON)
option(MIAM_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)

# Project wide setup variables
set(MIAM_INSTALL_INCLUDE_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(MIAM_INCLUDE_DIR ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
set(MIAM_LIB_DIR ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

add_subdirectory(src)

if(MIAM_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()