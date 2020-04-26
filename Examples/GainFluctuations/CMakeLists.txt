#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(GainFluctuations)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(arco2 arco2.C)
target_link_libraries(arco2 Garfield)

