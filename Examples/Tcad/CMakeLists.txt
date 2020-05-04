#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(TCAD)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(signal signal.C)
target_link_libraries(signal Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file pixel_des.dat pixel_des.grd)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
