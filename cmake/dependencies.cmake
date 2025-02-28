include(FetchContent)

################################################################################
# Google test

if(MIAM_ENABLE_TESTS)
  FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.16.0
  )

  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  # Don't run clang-tidy on google test
  set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")
  set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
endif()

################################################################################
# Docs

if(MIAM_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

################################################################################
# Memory check

if(MIAM_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
endif()