################################################################################
# Kokkos-backed integration tests
#
# Gated on MIAM_ENABLE_KOKKOS. When enabled, MICM is fetched with Kokkos
# support (see cmake/dependencies.cmake), which makes the Kokkos::kokkos
# target available for linking.
#
# Kokkos tests supply their own main() that calls Kokkos::initialize /
# ::finalize around RUN_ALL_TESTS, so they link against gtest (not
# gtest_main) and don't go through create_standard_test.

function(create_kokkos_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES LIBRARIES)
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} "" "${singleValues}" "${multiValues}" ${ARGN})

  add_executable(test_${TEST_NAME} ${TEST_SOURCES})
  target_link_libraries(test_${TEST_NAME} PUBLIC miam GTest::gtest Kokkos::kokkos)
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
endfunction()

create_kokkos_test(
  NAME kokkos_cam_cloud_chemistry
  SOURCES test_kokkos_cam_cloud_chemistry.cpp
)
