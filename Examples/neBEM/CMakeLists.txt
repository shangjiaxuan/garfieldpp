#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(neBEM)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executables-----------------------------------------------------------
add_executable(test test.C)
target_link_libraries(test Garfield)

