cmake_minimum_required(VERSION 3.21)

project(miam VERSION 0.0.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

message(STATUS "CMake build configuration for MIAM (${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}")

include(GNUInstallDirs)

################################################################################
# Projet wide setup options

option(MIAM_ENABLE_TESTS "Build the tests" ON)

add_subdirectory(src)

if(MIAM_ENABLE_TESTS)
  add_subdirectory(test)
endif()