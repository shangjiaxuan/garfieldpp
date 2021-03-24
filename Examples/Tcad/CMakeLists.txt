# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(TCAD)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Build executable------------------------------------------------------------
add_executable(signal signal.C)
target_link_libraries(signal Garfield::Garfield)

# ---Copy all files locally to the build directory-------------------------------
foreach(_file pixel_des.dat pixel_des.grd)
  configure_file(${_file} ${_file} COPYONLY)
endforeach()
