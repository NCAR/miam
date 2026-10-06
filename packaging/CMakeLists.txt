include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

# MICM is either fetched and built here (then MIAM installs it with its own targets) or found
# as an installed package (then miamConfig.cmake finds it again with find_dependency).
set(MIAM_INSTALL_TARGETS miam)
set(MIAM_BUNDLES_MICM OFF)
if(TARGET micm)
  get_target_property(_micm_is_imported micm IMPORTED)
  if(NOT _micm_is_imported)
    set(MIAM_BUNDLES_MICM ON)
    list(APPEND MIAM_INSTALL_TARGETS micm)
    if(MIAM_ENABLE_KOKKOS)
      list(APPEND MIAM_INSTALL_TARGETS micm_kokkos)
    endif()
  endif()
endif()

install(
  TARGETS
    ${MIAM_INSTALL_TARGETS}
  EXPORT
    miam_Exports
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/include/miam
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}
)

if(MIAM_BUNDLES_MICM)
  install(
    DIRECTORY
      ${micm_SOURCE_DIR}/include/micm
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}
  )
endif()

# install the cmake config files
set(cmake_config_install_location ${CMAKE_INSTALL_LIBDIR}/cmake/miam)

install(
  EXPORT
    miam_Exports
  DESTINATION
    ${cmake_config_install_location}
  NAMESPACE musica::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/miamConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/miamConfig.cmake"
  INSTALL_DESTINATION
    ${cmake_config_install_location}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/miamConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/miamConfig.cmake
    ${PROJECT_BINARY_DIR}/miamConfigVersion.cmake
  DESTINATION
    ${cmake_config_install_location}
)

######################################################################
# uninstall target

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
