################################################################################
# MIAM library
add_library(miam)
add_library(musica::miam ALIAS miam)

target_compile_features(miam PUBLIC cxx_std_20)

target_include_directories(miam
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${MIAM_INSTALL_INCLUDE_DIR}>
)

target_link_libraries(miam
  PUBLIC
    musica::micm
)

set_target_properties(miam PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${MIAM_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

add_subdirectory(model)
add_subdirectory(util)