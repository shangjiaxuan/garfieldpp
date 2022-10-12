set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)
set(libname "magboltz")
set(lib_fullver 11.9)
set(lib_apiver 11)

# Creating the target as a shared library
message(STATUS "Creating a new library named ${libname} with API version ${lib_apiver}")
add_library(${libname} SHARED "")
add_library(Garfield::${libname} ALIAS ${libname})
target_sources(${libname} PRIVATE magboltz.f)
target_compile_definitions(${libname} PRIVATE BUILDING_MEGABOLTZ=1)

# windows dll building need either __declspec(dllexport) or a .def file to specify symbols
if(CMAKE_HOST_WIN32)
  target_link_options(${libname} PRIVATE /DEF:${CMAKE_CURRENT_SOURCE_DIR}/magboltz.def)
  target_compile_options(${libname} PRIVATE "/nofpp")
endif()

# assume only intel compilers do not use "-std=legacy"
if(NOT (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM"))
  target_compile_options(${libname} PRIVATE -std=legacy)
endif()

# Sets the API version of the library
set_property(TARGET ${libname} PROPERTY SOVERSION ${lib_apiver})

install(
  TARGETS ${libname}
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
