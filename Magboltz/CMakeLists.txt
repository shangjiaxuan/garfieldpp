# required if building on msvc+gfortran, subproject uses mingw build system,
# but when using addsubdir logic from cmake_add_fortran_subdirectory, library
# will need extra work with new sub-project here
if(GARFIELD_USE_SUBPROJ_MINGW)
  # need this in external project
  cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
  message(STATUS "building Megaboltz shared library as subproject")
  project(Magboltz_mingw LANGUAGES Fortran)
else()
  # added with direct add_subdirectory, variables inherited
  if(MSVC)
    # ifort automatically exported symbols will be in the form:
    # MEGABOLTZ, instead of megaboltz_, and will not link from cpp code
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS OFF)
    # need 3.13 to use target_link_options, or symbols will not be same even with CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS
    cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
  endif()
  message(STATUS "building Megaboltz library as a library target in parent Garfield project")
endif()
set(libname "magboltz")
set(lib_fullver 11.9)
set(lib_apiver 11)

# Creating the target as a shared library
message(STATUS "Creating a new library named ${libname} with API version ${lib_apiver}")
add_library(${libname} SHARED "")
target_sources(${libname} PRIVATE magboltz.f)
target_compile_definitions(${libname} PRIVATE BUILDING_MEGABOLTZ=1)

# Windows dll building need either __declspec(dllexport), a .def file, or CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS to specify symbols.
# Automatically generated symbols does not use gnu style xxx_, but is in XXX
if(MSVC)
  target_link_options(${libname} PRIVATE /DEF:${CMAKE_CURRENT_SOURCE_DIR}/magboltz.def)
endif()

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" OR CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
# ifort will not compile with hard-coded cmake /fpp, thus need this
  target_compile_options(${libname} PRIVATE "-nofpp")
else()
# assume only intel compilers do not use "-std=legacy"
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
