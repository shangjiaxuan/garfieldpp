###############################################################################
# build Garfield++ library and examples
# depend on: ROOT
# @author Klaus Zenker, Pere Mato, Stefano Caiazza
# Version 3.9 is required by the ROOT Use file
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
################################################################################

# Setting a default build type to be release with debug information.
# Otherwise the default mode is debug.
set(default_build_type "Release" )
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  	set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      	STRING "Choose the type of build.")
  	# Set the possible values of build type for cmake-gui
  	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
		"Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

project(Garfield LANGUAGES CXX Fortran VERSION 0.3.0)

#--- Set default installation prefix
# TODO Possibly change the default install prefix simply to GARFIELD_HOME
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(DEFINED ENV{GARFIELD_HOME})
    set(CMAKE_INSTALL_PREFIX $ENV{GARFIELD_HOME}/install CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
  endif()
else()
  message(STATUS "The installation folder is set to its default value ${CMAKE_INSTALL_PREFIX}")
endif()

#Adding the CMake folder to the search path for modules
list( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

#---- Use GNU-style hierarchy for installing build products
include(GNUInstallDirs)
include(BuildUtils)
force_color_output()

#--- Add neBEM------------------------------------------------------------------
add_subdirectory(NeBem)

#--- Add Magboltz
add_subdirectory(Magboltz)

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
if(${ROOT_VERSION} VERSION_LESS "6.20.0")
  #--- Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
  # it only makes sense for ROOT < 6.20. Afterwards it is loaded automatically by the find_package
  get_filename_component(Root_Use_File_Dir ${ROOT_USE_FILE} DIRECTORY)
  include("${Root_Use_File_Dir}/modules/RootNewMacros.cmake")
endif()

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

find_package(OpenMP)

# Creating the main Garfield library
add_library(Garfield SHARED "")
target_sources(Garfield PRIVATE Source/AvalancheMC.cc
                                Source/AvalancheMicroscopic.cc
                                Source/Component.cc
                                Source/ComponentAnalyticField.cc
                                Source/ComponentAnsys121.cc
                                Source/ComponentAnsys123.cc
                                Source/ComponentCST.cc
                                Source/ComponentComsol.cc
                                Source/ComponentConstant.cc
                                Source/ComponentElmer.cc
                                Source/ComponentElmer2D.cc
                                Source/ComponentFieldMap.cc
                                Source/ComponentGrid.cc
                                Source/ComponentNeBem2d.cc
                                Source/ComponentNeBem3d.cc
                                Source/ComponentNeBem3dMap.cc
                                Source/ComponentTcad2d.cc
                                Source/ComponentTcad3d.cc
                                Source/ComponentUser.cc
                                Source/ComponentVoxel.cc
                                Source/DriftLineRKF.cc
                                Source/GeometryRoot.cc
                                Source/GeometrySimple.cc
                                Source/KDTree.cc
                                Source/Medium.cc
                                Source/MediumCdTe.cc
                                Source/MediumGaAs.cc
                                Source/MediumGaN.cc
                                Source/MediumGas.cc
                                Source/MediumMagboltz.cc
                                Source/MediumSilicon.cc
                                Source/NeBemInterface.cpp
                                Source/Numerics.cc
                                Source/OpticalData.cc
                                Source/PlottingEngineRoot.cc
                                Source/Polygon.cc
                                Source/Random.cc
                                Source/RandomEngineRoot.cc
                                Source/Sensor.cc
                                Source/Shaper.cc
                                Source/Solid.cc
                                Source/SolidBox.cc
                                Source/SolidExtrusion.cc
                                Source/SolidHole.cc
                                Source/SolidRidge.cc
                                Source/SolidSphere.cc
                                Source/SolidTube.cc
                                Source/SolidWire.cc
                                Source/TetrahedralTree.cc
                                Source/Track.cc
                                Source/TrackBichsel.cc
                                Source/TrackElectron.cc
                                Source/TrackHeed.cc
                                Source/TrackPAI.cc
                                Source/TrackSimple.cc
                                Source/TrackSrim.cc
                                Source/ViewBase.cc
                                Source/ViewCell.cc
                                Source/ViewDrift.cc
                                Source/ViewFEMesh.cc
                                Source/ViewField.cc
                                Source/ViewGeometry.cc
                                Source/ViewIsochrons.cc
                                Source/ViewMedium.cc
                                Source/ViewSignal.cc
)

#--- create ROOT dictionary-----------------------------------------------------
if(ROOT_VERSION VERSION_LESS 6.17)
  file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
  ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                           LINKDEF Include/Garfield/LinkDef.h
                           OPTIONS -I${CMAKE_CURRENT_SOURCE_DIR}/Include -I${PROJECT_SOURCE_DIR}/Heed)
else()
  file (GLOB garfield_header RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/Include Include/Garfield/*.hh)
  include_directories(Include Heed)
  ROOT_GENERATE_DICTIONARY(GarfieldDict ${garfield_header}
                           LINKDEF Include/Garfield/LinkDef.h)
endif()
target_sources(Garfield PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/GarfieldDict.cxx)

#--- Install pcm file in case of root 6
if(ROOT_VERSION VERSION_GREATER 6.06)
    install( FILES ${CMAKE_CURRENT_BINARY_DIR}/libGarfieldDict_rdict.pcm 
             DESTINATION ${CMAKE_INSTALL_LIBDIR}
             COMPONENT dict)
endif()

garfield_set_all_default_properties(Garfield)

#--- Set the C++ standard (depending on the ROOT version)-----------------------
if(NOT DEFINED CMAKE_CXX_STANDARD)
  if(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+([0-9][0-9])")
    set_property( TARGET Garfield PROPERTY CXX_STANDARD ${CMAKE_MATCH_1})
  elseif(${ROOT_CXX_FLAGS} MATCHES ".*-std=c\\+\\+1z")
    set_property( TARGET Garfield PROPERTY CXX_STANDARD 17)
  endif()
endif()

set_property( TARGET Garfield PROPERTY CMAKE_CXX_STANDARD_REQUIRED ON)
set_property( TARGET Garfield PROPERTY CXX_EXTENSIONS OFF)

#--- If you cannot create symbolic links turn FALSE the library version
option(BUILD_WITH_VERSION_NUMBER "If you cannot create symbolic links turn FALSE the library version" TRUE)
if(BUILD_WITH_VERSION_NUMBER)
  set_target_properties(Garfield PROPERTIES SOVERSION ${Garfield_VERSION} )
endif()
mark_as_advanced(BUILD_WITH_VERSION_NUMBER)

#--- Add Heed sources-------------------------------------------------------------------
add_subdirectory(Heed)

target_link_libraries(Garfield PUBLIC ROOT::Geom ROOT::Gdml ROOT::Graf3d ${GSL_LIBRARIES} NeBem magboltz)
if(OpenMP_CXX_FOUND)
  target_link_libraries(Garfield PRIVATE OpenMP::OpenMP_CXX)
endif()
target_include_directories(Garfield PUBLIC $<INSTALL_INTERFACE:include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Include>
                                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Heed>)

#--- target for Doxygen documentation ------------------------------------------
option(WITH_DOC "Whether or not to create doxygen doc target." OFF)
if(WITH_DOC)
  add_subdirectory(Doc)
endif()

#--- Build the examples --------------------------------------------------------
option(WITH_EXAMPLES "Build Garfield++ examples" ON)
if(WITH_EXAMPLES)
  message(STATUS "Garfield++ examples will be built. In order disable this option use -DWITH_EXAMPLES=OFF")
  add_subdirectory(Examples)
else()
  message(STATUS "Building of examples is switched off.")
endif()

#--- add CMake infrastructure --------------------------------------------------
export(PACKAGE Garfield)
include(CMake/GarfieldCreateConfig.cmake)

configure_file(${CMAKE_SOURCE_DIR}/CMake/setup.sh.in
               ${CMAKE_BINARY_DIR}/setup.sh @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/CMake/GarfieldConfig.dummy.in
               ${CMAKE_BINARY_DIR}/GarfieldConfig.cmake @ONLY)



#--- Installation of all artifacts----------------------------------------------
message(STATUS "Install path: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "If you want to change this path call cmake -DCMAKE_INSTALL_PREFIX=my_install_path ...")
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

#--- Where the library will end up and to write info in GarfieldTargets.cmake--
install(TARGETS Garfield DESTINATION ${CMAKE_INSTALL_LIBDIR} EXPORT GarfieldTargets)

#--- install the Heed database--------------------------------------------------
install(DIRECTORY Heed/heed++/database DESTINATION Heed/heed++)

#--- install headers------------------------------------------------------------
install(DIRECTORY Include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} REGEX LinkDef EXCLUDE)
