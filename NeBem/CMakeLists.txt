set(libname NeBem)
set(lib_apiver 1)

# Creating the target as a shared library
message(STATUS "Creating a new library named ${libname} with API version ${lib_apiver}")
set(nebem_sources
    ComputeProperties.c
    Isles.c
    luc.c
    neBEM.c
    neBEMInterface.c
    nrutil.c
    ReTriM.c
    svdcmp.c
    Vector.c)
set_source_files_properties(${nebem_sources} PROPERTIES LANGUAGE CXX)

add_library(${libname} SHARED "")
add_library(Garfield::${libname} ALIAS ${libname})
target_sources(${libname} PRIVATE ${nebem_sources})

garfield_set_all_default_properties(${libname})
# Sets the API version of the library
set_property(TARGET ${libname} PROPERTY SOVERSION ${lib_apiver})
set_property(TARGET ${libname} PROPERTY CXX_STANDARD 11)
set_property(TARGET ${libname} PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)
set_property(TARGET ${libname} PROPERTY CXX_EXTENSIONS OFF)

target_link_libraries(${libname} PUBLIC GSL::gsl)

target_include_directories(${libname}
                           PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(
  TARGETS ${libname}
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
