########################################################
# @author Klaus Zenker
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
########################################################
project(GarfieldExamples)

#--- GEM
add_subdirectory(Gem)

#--- Drift tube 
add_subdirectory(DriftTube)

#--- Gas files
add_subdirectory(GasFile)

#--- Analytic fields
add_subdirectory(AnalyticField)

#--- Regular grid
add_subdirectory(Grid)

#--- Heed
add_subdirectory(Heed)

#--- SRIM
add_subdirectory(Srim)

#--- Geant4-Garfield interface
add_subdirectory(Geant4GarfieldInterface)

#--- neBEM
add_subdirectory(neBEM)
