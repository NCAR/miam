################################################################################
# miam version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/miam/version.hpp @ONLY)

################################################################################
# miam library

add_library(miam INTERFACE)
add_library(musica::miam ALIAS miam)

set_property(TARGET miam PROPERTY INTERFACE_CXX_STANDARD 20)
set_property(TARGET miam PROPERTY INTERFACE_CXX_STANDARD_REQUIRED YES)

target_include_directories(miam
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${MIAM_INSTALL_INCLUDE_DIR}>
)

target_link_libraries(miam
  INTERFACE
    musica::micm
)

set_target_properties(miam PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)