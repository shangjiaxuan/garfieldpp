#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Import)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(Test Simulatipn.cpp)
target_link_libraries(Test Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file Wfield.txt Efield.txt Attachment.txt)
  configure_file(${_file} ${_file} ${_file} COPYONLY)
endforeach()
