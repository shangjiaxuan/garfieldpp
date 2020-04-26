#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(DriftTube)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(signal signal.C)
target_link_libraries(signal Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file ar_93_co2_7_3bar.gas mdt_elx_delta.txt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()

