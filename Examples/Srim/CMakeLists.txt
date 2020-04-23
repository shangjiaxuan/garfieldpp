#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Srim)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(srim srim.C)
target_link_libraries(srim Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file alpha_ArCO2_70_30.txt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

