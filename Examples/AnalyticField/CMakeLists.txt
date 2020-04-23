#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(AnalyticField)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Define executables----------------------------------------------------------
add_executable(gallery gallery.C)
target_link_libraries(gallery Garfield)

add_executable(integrate integrate.C)
target_link_libraries(integrate Garfield)

add_executable(dipole dipole.C)
target_link_libraries(dipole Garfield)
