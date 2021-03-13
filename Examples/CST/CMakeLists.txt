#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(CST)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(cst cst.C)
target_link_libraries(cst Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file drift_std_200_320_2500.efm)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

