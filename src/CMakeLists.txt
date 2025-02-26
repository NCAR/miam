################################################################################
# miam version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/miam/version.hpp @ONLY)

################################################################################
# miam solver
add_library(miam INTERFACE)
add_library(musica::miam ALIAS miam)

target_compile_features(miam INTERFACE cxx_std_20)

target_include_directories(miam
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Check for Intel compiler
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" OR CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  # Set compiler flags for Intel C++
  if (WIN32)
    # For Windows
    target_compile_options(miam INTERFACE "/fp:precise")
  else ()
    # For other operating systems (e.g., Linux, macOS)
    target_compile_options(miam INTERFACE "-ffp-model=precise")
  endif ()
endif ()