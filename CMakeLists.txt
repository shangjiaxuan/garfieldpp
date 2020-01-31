########################################################
# build Garfield++ library and examples
# depend on: ROOT
# @author Klaus Zenker, Pere Mato
cmake_minimum_required(VERSION 3.3)
########################################################

project(Garfield LANGUAGES CXX Fortran VERSION 0.3.0)

#---- Use GNU-style hierarchy for installing build products
include(GNUInstallDirs)
set(CMAKE_INSTALL_LIBDIR lib)  # overwrite this one

#--- Dependencies
if(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
elseif(DEFINED ROOTSYS)
    list(APPEND CMAKE_PREFIX_PATH ${ROOTSYS})
elseif(DEFINED ROOT_CMAKE_DIR)
    list(APPEND CMAKE_PREFIX_PATH ${ROOT_CMAKE_DIR} )
else()
    message(FATAL_ERROR "Please consider to set ROOTSYS or use -DROOTSYS=..."
                        "If still there is a proplem, point to the directory which includes FindROOT.cmake using -DROOT_CMAKE_DIR=...")
endif(DEFINED ENV{ROOTSYS})


find_package(ROOT 5.0 REQUIRED COMPONENTS Geom Gdml)
#--- Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#--- Internal ROOT dictionary building works starting with version 6.02
if(ROOT_VERSION VERSION_LESS 6.02)
    set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake )
    include( MacroCheckPackageLibs )
    include( MacroRootDict )
    message(STATUS "Using Garfield's FindRoot.cmake")
endif()

#--- If you cannot create symbolic links turn FALSE the library version
set(BUILD_WITH_VERSION_NUMBER TRUE)
if(BUILD_WITH_VERSION_NUMBER)
    set( Garfield_SOVERSION ${Garfield_VERSION} )
endif()

#--- Add ROOT libraries
if(DEFINED ROOT_INCLUDE_DIRS)
    # ROOT 6
    include_directories ( ${ROOT_INCLUDE_DIRS} )
    set(ROOT_COMPONENT_LIBRARIES ROOT::Geom ROOT::Gdml ROOT::Graf3d) 
else()
    # ROOT 5
    include_directories ( ${ROOT_INCLUDE_DIR} )
    # COMPONENT mechanism in the find_package is not supported 
    set(ROOT_COMPONENT_LIBRARIES "-lGeom -lGdml")
endif()

#--- Add Garfield
aux_source_directory(Source garfield_sources)

#--- Add Magboltz
set(magboltz_sources Source/magboltz-11.9.f )
set_source_files_properties(${magboltz_sources} PROPERTIES LANGUAGE Fortran )

#--- Add Heed
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
set_source_files_properties( ${heed_sources} PROPERTIES LANGUAGE CXX)

#--- Flags to pass to the compiler
if(APPLE)
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -ansi -Wabi -Wno-long-long -Woverloaded-virtual -fno-common" )
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic -ansi -Wabi -Wno-long-long -Woverloaded-virtual -fno-common -pg" )
endif()

#--- Set the C++ standard (depending on the ROOT version)
if(NOT DEFINED CMAKE_CXX_STANDARD)
  if(ROOT_VERSION VERSION_GREATER 6.06)
    set(CMAKE_CXX_STANDARD 14)
  elseif(ROOT_VERSION VERSION_GREATER 5)
    set(CMAKE_CXX_STANDARD 98)
  endif()
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#--- Add switches used in HEED
#    -DNOT_INCLUDE_GASLIB_IN_PACSLIB is used in Heed/heed++/code/PhotoAbsCS.c
set_source_files_properties(${heed_sources} PROPERTIES COMPILE_FLAGS "-DNOT_INCLUDE_GASLIB_IN_PACSLIB -DGARFIELD_HEED_INTERFACE -DUSE_SRANLUX -DEXCLUDE_FUNCTIONS_WITH_HISTDEF -DINS_CRETURN -DFUNNAMESTACK")

#--- create ROOT dictionary
if(ROOT_VERSION VERSION_LESS 6.02)
    list(APPEND ROOT_DICT_INCLUDE_DIRS Include/Garfield Heed)
    set(ROOT_DICT_CINT_DEFINITIONS "-DHANDLE_DICT_EXCEPTIONS=IGNORED_FOR_CINT" )
    PREPARE_ROOT_DICT_HEADERS( Include/Garfield )
    GEN_ROOT_DICT_SOURCES( GarfieldDict.cxx )
    list(APPEND garfield_sources ${ROOT_DICT_OUTPUT_SOURCES})
elseif(ROOT_VERSION VERSION_LESS 6.17)
    file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
    ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                             LINKDEF Include/Garfield/LinkDef.h
                             OPTIONS -I${PROJECT_SOURCE_DIR}/Include -I${PROJECT_SOURCE_DIR}/Heed)
    list(APPEND garfield_sources GarfieldDict.cxx)
else()
    file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
    include_directories(Include Heed)
    ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                             LINKDEF Include/Garfield/LinkDef.h )
    list(APPEND garfield_sources GarfieldDict.cxx)
endif()

#--- Print information on considered files and the root dictionary include directories
message(VERBOSE "-----------------------------------------------------------" )
message(VERBOSE "Considered Garfield files: ${garfield_sources}" )
message(VERBOSE "-----------------------------------------------------------" )
message(VERBOSE "Considered Heed files: ${heed_sources}" )
message(VERBOSE "-----------------------------------------------------------" )
message(VERBOSE "Considered Magboltz files: ${magboltz_sources}" )
message(VERBOSE "-----------------------------------------------------------" )
message(VERBOSE "Directory used to find root dict files: ${ROOT_DICT_INCLUDE_DIRS}" )
message(VERBOSE "-----------------------------------------------------------" )

#--- Set output path for libraries and executables
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

#--- build the Garfield library
add_library(Garfield SHARED ${garfield_sources} ${magboltz_sources} ${heed_sources})
target_link_libraries(Garfield ${ROOT_COMPONENT_LIBRARIES})
target_compile_features(Garfield INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
target_include_directories(Garfield PUBLIC $<INSTALL_INTERFACE:include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Heed>)
if(BUILD_WITH_VERSION_NUMBER)
    set_target_properties(Garfield PROPERTIES
                          VERSION ${Garfield_VERSION}
                          SOVERSION ${Garfield_SOVERSION} )
endif()

#--- Build the examples --------------------------------------------------------
option(WITH_EXAMPLES "Build Garfield++ examples" ON)
if(WITH_EXAMPLES)
    message(STATUS "Garfield++ examples will be built. In order disable this option use -DWITH_EXAMPLES=OFF")
    add_subdirectory(Examples)
else()
    message(STATUS "Building of examples is switched off.")
endif()
configure_file(${CMAKE_SOURCE_DIR}/CMake/setup.sh.in
               ${CMAKE_BINARY_DIR}/setup.sh @ONLY)


#--- target for Doxygen documentation ------------------------------------------
option(WITH_DOC "Whether or not to create doxygen doc target." OFF)
if(WITH_DOC)
  include(cmake/GarfieldDoxygen.cmake)
endif()

#--- add CMake infrastructure --------------------------------------------------
export(PACKAGE Garfield)
include(CMake/GarfieldCreateConfig.cmake)

#--- Installation of all artifacts
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ...")
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

#--- Where the library will end up and to write info in  GarfieldTargets.cmake
install(TARGETS Garfield DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT GarfieldTargets)

#--- Install pcm file in case of root 6
if(ROOT_VERSION VERSION_GREATER 6.06)
    install( FILES ${PROJECT_BINARY_DIR}/libGarfieldDict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_LIBDIR}
             COMPONENT dict)
elseif(ROOT_VERSION VERSION_GREATER 5)
    install( FILES ${PROJECT_BINARY_DIR}/GarfieldDict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_LIBDIR}
             COMPONENT dict)
endif()

#--- install the Heed database
install(DIRECTORY Heed/heed++/database DESTINATION  Heed/heed++)
#--- install headers
install(DIRECTORY Include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} REGEX LinkDef EXCLUDE)
