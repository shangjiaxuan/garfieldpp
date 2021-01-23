#---Build executable------------------------------------------------------------
add_executable(signal signal.C)
target_link_libraries(signal Garfield)

add_executable(signal_mt signal_mt.C)
target_link_libraries(signal_mt Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file ar_93_co2_7_3bar.gas mdt_elx_delta.txt)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

