################################################################################
# Find Valgrind

if(MIAM_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1" "--trace-children=yes" "--leak-check=full" "--gen-suppressions=all")
endif()

################################################################################
# Create a standard test that links to the MIAM library

function(create_standard_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES LIBRARIES)
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} " " "${singleValues}" "${multiValues}" ${ARGN})

  add_executable(test_${TEST_NAME} ${TEST_SOURCES})
  target_link_libraries(test_${TEST_NAME} PUBLIC miam GTest::gtest_main)

  # Link additional libraries
  foreach(library ${TEST_LIBRARIES})
    target_link_libraries(test_${TEST_NAME} PUBLIC ${library})
  endforeach()

  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  endif()

  add_test(NAME ${TEST_NAME}
          COMMAND test_${TEST_NAME}
          WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY})

  if(MIAM_ENABLE_MEMCHECK AND MEMORYCHECK_COMMAND)
    add_test(NAME memcheck_${TEST_NAME}
            COMMAND ${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS} ${TEST_WORKING_DIRECTORY}/test_${TEST_NAME}
            WORKING_DIRECTORY ${TEST_WORKING_DIRECTORY})
  endif()

endfunction(create_standard_test)
