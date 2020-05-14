#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(WIRE2D)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(wire2d_example wire2d.cc)
target_link_libraries(wire2d_example Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file wire2d/mesh.header wire2d/mesh.elements wire2d/mesh.nodes wire2d/dielectrics.dat wire2d/wire2d.result)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

