#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(neBEM2d)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executables-----------------------------------------------------------
add_executable(triangle triangle.C)
target_link_libraries(triangle Garfield)

add_executable(dielectric dielectric.C)
target_link_libraries(dielectric Garfield)

add_executable(wire wire.C)
target_link_libraries(wire Garfield)
