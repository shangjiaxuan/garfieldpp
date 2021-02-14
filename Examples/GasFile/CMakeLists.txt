#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(GasFile)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Define executables----------------------------------------------------------
add_executable(generate generate.C)
target_link_libraries(generate Garfield)

add_executable(read read.C)
target_link_libraries(read Garfield)

add_executable(merge merge.C)
target_link_libraries(merge Garfield)

add_executable(penning penning.C)
target_link_libraries(penning Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file ar_80_co2_20_0T.gas ar_80_co2_20_2T.gas ar_93_co2_7.gas)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
