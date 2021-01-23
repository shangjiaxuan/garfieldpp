#---Define executables----------------------------------------------------------
add_executable(plotfield plotfield.C)
target_link_libraries(plotfield Garfield)

add_executable(savefield savefield.C)
target_link_libraries(savefield Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file Efield.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
