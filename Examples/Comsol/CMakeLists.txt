#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(GemComsol)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(gem gem.C)
target_link_libraries(gem Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file dielectrics.dat field.txt mesh.mphtxt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

