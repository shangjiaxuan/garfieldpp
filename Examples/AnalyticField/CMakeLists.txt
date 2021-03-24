# ########################################################################################
# Use the code in this section for standalone projects
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(AnalyticField)
if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()
# ########################################################################################

# ---Define executables----------------------------------------------------------
add_executable(gallery gallery.C)
target_link_libraries(gallery Garfield::Garfield)

add_executable(integrate integrate.C)
target_link_libraries(integrate Garfield::Garfield)

add_executable(dipole dipole.C)
target_link_libraries(dipole Garfield::Garfield)

add_executable(isochrons isochrons.C)
target_link_libraries(isochrons Garfield::Garfield)

add_executable(polar polar.C)
target_link_libraries(polar Garfield::Garfield)

add_executable(wedge wedge.C)
target_link_libraries(wedge Garfield::Garfield)

add_executable(strip strip.C)
target_link_libraries(strip Garfield::Garfield)

foreach(_file ar_50_c2h6_50_B_angle.gas)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
