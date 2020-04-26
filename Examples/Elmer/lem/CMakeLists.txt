#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(LEM)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(avalanche avalanche.cc)
target_link_libraries(avalanche Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file gemcell/mesh.header gemcell/mesh.elements gemcell/mesh.nodes gemcell/dielectrics.dat gemcell/gemcell.result)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

