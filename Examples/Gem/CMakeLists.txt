#---Build executable------------------------------------------------------------
add_executable(gem gem.C)
target_link_libraries(gem Garfield)

#---Copy all files locally to the build directory-------------------------------
foreach(_file ELIST.lis MPLIST.lis NLIST.lis PRNSOL.lis)
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

