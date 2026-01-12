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
# MICM

FetchContent_Declare(micm
    GIT_REPOSITORY https://github.com/NCAR/micm.git
    GIT_TAG rework_system_params
    GIT_PROGRESS NOT ${FETCHCONTENT_QUIET}
    FIND_PACKAGE_ARGS NAMES micm
)

set(MICM_ENABLE_TESTS OFF)
set(MICM_ENABLE_EXAMPLES OFF)

FetchContent_MakeAvailable(micm)

################################################################################
# Docs

if(MIAM_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()