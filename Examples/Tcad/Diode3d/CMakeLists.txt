#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Diode3d)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(diode diode.cpp)
target_link_libraries(diode Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file diode_des.dat diode_msh.grd)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

