################################################################################
# MIAM version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/miam/version.hpp @ONLY)

################################################################################
# MIAM library
add_library(miam)
add_library(musica::miam ALIAS miam)

target_compile_features(miam PUBLIC cxx_std_20)

# Check for Intel compiler
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  # Set compiler flags for Intel C++
  if (WIN32)
    # For Windows
    target_compile_options(miam PUBLIC "/fp:precise")
  else ()
    # For other operating systems (e.g., Linux, macOS)
    target_compile_options(miam PUBLIC "-ffp-model=precise")
  endif ()
endif ()

target_sources(miam
  PRIVATE
    ${CMAKE_BINARY_DIR}/version.cpp
)

target_include_directories(miam
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${MIAM_INSTALL_INCLUDE_DIR}>
)

set_target_properties(miam PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

add_subdirectory(util)