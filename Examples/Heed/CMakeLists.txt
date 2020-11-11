#---Setup the example project---------------------------------------------------
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(Heed)

#---Find Garfield package-------------------------------------------------------
find_package(Garfield REQUIRED)

#---Build executable------------------------------------------------------------
add_executable(edep edep.C)
target_link_libraries(edep Garfield)

add_executable(qdepSi qdepSi.C)
target_link_libraries(qdepSi Garfield)

add_executable(rangeSi rangeSi.C)
target_link_libraries(rangeSi Garfield)

add_executable(fe55 fe55.C)
target_link_libraries(fe55 Garfield)

add_executable(plotdedx plotdedx.C)
target_link_libraries(plotdedx Garfield)
