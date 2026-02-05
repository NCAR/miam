################################################################################
# Test utilities

include(test_util)

################################################################################
# Disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

add_subdirectory(integration)
add_subdirectory(unit)