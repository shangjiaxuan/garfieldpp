#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Grid)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Define executables----------------------------------------------------------
add_executable(plotfield plotfield.C)
target_link_libraries(plotfield Garfield)

add_executable(savefield savefield.C)
target_link_libraries(savefield Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file Efield.txt)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
