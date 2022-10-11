set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)
set(libname "magboltz")
set(lib_fullver 11.9)
set(lib_apiver 11)

# Creating the target as a shared library
message(STATUS "Creating a new library named ${libname} with API version ${lib_apiver}")
add_library(${libname} SHARED "")
add_library(Garfield::${libname} ALIAS ${libname})
target_sources(${libname} PRIVATE magboltz.f)
target_link_options(${libname} PRIVATE /DEF:${CMAKE_CURRENT_SOURCE_DIR}/magboltz.def)
target_compile_definitions(${libname} PRIVATE BUILDING_MEGABOLTZ=1)
# Sets the API version of the library

install(
  TARGETS ${libname}
  EXPORT "${PROJECT_NAME}Targets"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
