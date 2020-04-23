# - Use CMake's module to help generating relocatable config files
include(CMakePackageConfigHelpers)

# - Versioning
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/installtree/GarfieldConfigVersion.cmake
  VERSION ${Garfield_VERSION}
  COMPATIBILITY SameMajorVersion)

# - Install time config and target files
configure_package_config_file(${CMAKE_CURRENT_LIST_DIR}/GarfieldConfig.cmake.in
  "${PROJECT_BINARY_DIR}/installtree/GarfieldConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Garfield"
  PATH_VARS
    CMAKE_INSTALL_BINDIR
    CMAKE_INSTALL_INCLUDEDIR
    CMAKE_INSTALL_LIBDIR
  )

# - install and export
install(FILES
  "${PROJECT_BINARY_DIR}/installtree/GarfieldConfigVersion.cmake"
  "${PROJECT_BINARY_DIR}/installtree/GarfieldConfig.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Garfield"
  )
install(EXPORT GarfieldTargets
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Garfield"
  )

