# Disable clang tidy for test files
set(CMAKE_CXX_CLANG_TIDY "")

if(MIAM_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
endif()

add_subdirectory(util)