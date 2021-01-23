#---Define executables----------------------------------------------------------
add_executable(gallery gallery.C)
target_link_libraries(gallery Garfield)

add_executable(integrate integrate.C)
target_link_libraries(integrate Garfield)

add_executable(dipole dipole.C)
target_link_libraries(dipole Garfield)

add_executable(isochrons isochrons.C)
target_link_libraries(isochrons Garfield)

add_executable(polar polar.C)
target_link_libraries(polar Garfield)

add_executable(wedge wedge.C)
target_link_libraries(wedge Garfield)

foreach(_file ar_50_c2h6_50_B_angle.gas)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
