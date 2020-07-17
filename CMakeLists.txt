###############################################################################
# build Garfield++ library and examples
# depend on: ROOT
# @author Klaus Zenker, Pere Mato
cmake_minimum_required(VERSION 3.3)
################################################################################
project(Garfield LANGUAGES CXX Fortran VERSION 0.3.0)

#---- Use GNU-style hierarchy for installing build products
include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)  # overwrite this one

#--- Dependencies---------------------------------------------------------------
if(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
elseif(DEFINED ROOTSYS)
    list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
elseif(DEFINED ROOT_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${ROOT_CMAKE_DIR} )
else()
    message(FATAL_ERROR "Please consider to set ROOTSYS or use -DROOTSYS=..."
                        "If still there is a problem, point to the directory which includes FindROOT.cmake using -DROOT_CMAKE_DIR=...")
endif(DEFINED ENV{ROOTSYS})

find_package(ROOT 6.0 REQUIRED COMPONENTS Geom Gdml)
#--- Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

find_package(OpenMP)

#--- If you cannot create symbolic links turn FALSE the library version
set(BUILD_WITH_VERSION_NUMBER TRUE)
if(BUILD_WITH_VERSION_NUMBER)
    set( Garfield_SOVERSION ${Garfield_VERSION} )
endif()

#--- Add Garfield---------------------------------------------------------------
aux_source_directory(Source garfield_sources)

#--- Add Magboltz
set(magboltz_sources Source/magboltz-11.9.f )
set_source_files_properties(${magboltz_sources} PROPERTIES LANGUAGE Fortran)

#--- Add Heed-------------------------------------------------------------------
set( heed_dirs 
    Heed
    Heed/heed++/code 
    Heed/wcpplib/clhep_units 
    Heed/wcpplib/geometry 
    Heed/wcpplib/ioniz 
    Heed/wcpplib/math 
    Heed/wcpplib/matrix 
    Heed/wcpplib/matter 
    Heed/wcpplib/particle 
    Heed/wcpplib/random 
    Heed/wcpplib/safetl 
    Heed/wcpplib/stream 
    Heed/wcpplib/util )
foreach(heed_dir ${heed_dirs})
  aux_source_directory(${heed_dir} heed_sources)
endforeach()
set_source_files_properties(${heed_sources} PROPERTIES LANGUAGE CXX)

#--- Add neBEM------------------------------------------------------------------
aux_source_directory(NeBem nebem_sources)
set_source_files_properties(${nebem_sources} PROPERTIES LANGUAGE CXX)

#--- Flags to pass to the compiler----------------------------------------------
if(APPLE)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -Wno-long-long -Woverloaded-virtual -fno-common" )
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -ansi -Wno-long-long -Woverloaded-virtual -fno-common" )
endif()

#--- Set the C++ standard (depending on the ROOT version)-----------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
  if(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+([0-9][0-9])")
    set(CMAKE_CXX_STANDARD ${CMAKE_MATCH_1})
  elseif(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+1z")
    set(CMAKE_CXX_STANDARD 17)
  endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#--- Set default installation prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED ENV{GARFIELD_HOME})
    set(CMAKE_INSTALL_PREFIX $ENV{GARFIELD_HOME}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  endif()
endif()

#--- Add switches used in HEED--------------------------------------------------
# -DNOT_INCLUDE_GASLIB_IN_PACSLIB is used in Heed/heed++/code/PhotoAbsCS.c
set_source_files_properties(${heed_sources} PROPERTIES COMPILE_FLAGS "-DNOT_INCLUDE_GASLIB_IN_PACSLIB -DGARFIELD_HEED_INTERFACE -DFUNNAMESTACK")

#--- create ROOT dictionary-----------------------------------------------------
if(ROOT_VERSION VERSION_LESS 6.17)
  file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
  ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                           LINKDEF Include/Garfield/LinkDef.h
                           OPTIONS -I${PROJECT_SOURCE_DIR}/Include -I${PROJECT_SOURCE_DIR}/Heed)
else()
  file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
  include_directories(Include Heed)
  ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                           LINKDEF Include/Garfield/LinkDef.h)
endif()
list(APPEND garfield_sources GarfieldDict.cxx)

#--- Print information on considered files and the root dictionary include directories
message(STATUS "-----------------------------------------------------------" )
message(STATUS "Garfield source files: ${garfield_sources}" )
message(STATUS "-----------------------------------------------------------" )
message(STATUS "Heed source files: ${heed_sources}" )
message(STATUS "-----------------------------------------------------------" )
message(STATUS "Magboltz source files: ${magboltz_sources}" )
message(STATUS "-----------------------------------------------------------" )
message(STATUS "neBEM source files: ${nebem_sources}" )
message(STATUS "-----------------------------------------------------------" )
message(STATUS "Directory used to find ROOT dict files: ${ROOT_DICT_INCLUDE_DIRS}" )
message(STATUS "-----------------------------------------------------------" )

#--- build the Garfield library
add_library(Garfield SHARED ${garfield_sources} ${magboltz_sources} ${heed_sources} ${nebem_sources})
target_link_libraries(Garfield PUBLIC ROOT::Geom ROOT::Gdml ROOT::Graf3d ${GSL_LIBRARIES})
if(OpenMP_CXX_FOUND)
  target_link_libraries(Garfield PRIVATE OpenMP::OpenMP_CXX)
endif()
target_compile_features(Garfield INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
target_include_directories(Garfield PUBLIC $<INSTALL_INTERFACE:include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Heed>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/NeBem>)
if(BUILD_WITH_VERSION_NUMBER)
  set_target_properties(Garfield PROPERTIES
                        VERSION ${Garfield_VERSION}
                        SOVERSION ${Garfield_SOVERSION} )
endif()

#--- target for Doxygen documentation ------------------------------------------
option(WITH_DOC "Whether or not to create doxygen doc target." OFF)
if(WITH_DOC)
  include(cmake/GarfieldDoxygen.cmake)
endif()

#--- add CMake infrastructure --------------------------------------------------
export(PACKAGE Garfield)
include(CMake/GarfieldCreateConfig.cmake)

configure_file(${CMAKE_SOURCE_DIR}/CMake/setup.sh.in
               ${CMAKE_BINARY_DIR}/setup.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/CMake/GarfieldConfig.dummy.in
               ${CMAKE_BINARY_DIR}/GarfieldConfig.cmake @ONLY)

#--- Build the examples --------------------------------------------------------
option(WITH_EXAMPLES "Build Garfield++ examples" ON)
if(WITH_EXAMPLES)
  message(STATUS "Garfield++ examples will be built. In order disable this option use -DWITH_EXAMPLES=OFF")
  add_subdirectory(Examples)
else()
  message(STATUS "Building of examples is switched off.")
endif()

#--- Installation of all artifacts----------------------------------------------
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ...")
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

#--- Where the library will end up and to write info in GarfieldTargets.cmake--
install(TARGETS Garfield DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT GarfieldTargets)

#--- Install pcm file in case of root 6
if(ROOT_VERSION VERSION_GREATER 6.06)
    install( FILES ${PROJECT_BINARY_DIR}/libGarfieldDict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_LIBDIR}
             COMPONENT dict)
endif()

#--- install the Heed database--------------------------------------------------
install(DIRECTORY Heed/heed++/database DESTINATION Heed/heed++)

#--- install headers------------------------------------------------------------
install(DIRECTORY Include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} REGEX LinkDef EXCLUDE)
