#---Build executable------------------------------------------------------------
add_executable(srim srim.C)
target_link_libraries(srim Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file alpha_ArCO2_70_30.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

