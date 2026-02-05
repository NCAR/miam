################################################################################
# miam version
configure_file(version.hpp.in ${PROJECT_SOURCE_DIR}/include/miam/version.hpp @ONLY)

################################################################################
# miam library

target_sources(miam
PRIVATE
    stub.cpp
)