#---Define executables----------------------------------------------------------
add_executable(generate generate.C)
target_link_libraries(generate Garfield)

add_executable(read read.C)
target_link_libraries(read Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file ar_80_co2_20_0T.gas ar_80_co2_20_2T.gas)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()
