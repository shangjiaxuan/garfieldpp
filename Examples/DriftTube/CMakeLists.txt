# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(DriftTube)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

#---Build executable------------------------------------------------------------
add_executable(mdt mdt.C)
target_link_libraries(mdt Garfield::Garfield)

add_executable(mdt_mt mdt_mt.C)
target_link_libraries(mdt_mt Garfield::Garfield)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(mdt_mt OpenMP::OpenMP_CXX)
endif()

# ---Copy all files locally to the build directory-------------------------------
foreach(_file ar_93_co2_7_3bar.gas mdt_elx_delta.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
