#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(RPC)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(Weight Weight.C)
target_link_libraries(Weight Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file c2h2f4_93-5_iso_5_sf6_1-5_bis.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
