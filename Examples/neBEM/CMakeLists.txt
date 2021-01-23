#---Build executables-----------------------------------------------------------
add_executable(RPCbyneBEM RPCbyneBEM.C)
target_link_libraries(RPCbyneBEM Garfield)

add_executable(MWPCbyneBEM MWPCbyneBEM.C)
target_link_libraries(MWPCbyneBEM Garfield)

add_executable(BiconicalGemByneBEM BiconicalGemByneBEM.C)
target_link_libraries(BiconicalGemByneBEM Garfield)

add_executable(parallelPlates parallelPlates.C)
target_link_libraries(parallelPlates Garfield)

add_executable(crossingWires crossingWires.C)
target_link_libraries(crossingWires Garfield)

add_executable(parallelPlatesWeightingField parallelPlatesWeightingField.C)
target_link_libraries(parallelPlatesWeightingField Garfield)

