#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Paschen)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executables-----------------------------------------------------------
add_executable(gastable gastable.C)
target_link_libraries(gastable Garfield)

add_executable(paschenariso paschenariso.C)
target_link_libraries(paschenariso Garfield)
