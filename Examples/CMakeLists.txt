########################################################
# @author Klaus Zenker
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
########################################################
project(GarfieldExamples)

#### AnalyticField         #############################
add_subdirectory(AnalyticField)

#### Voxel                 #############################
add_subdirectory(Voxel)

#### GEM                   #############################
add_subdirectory(Gem)

#### GasFile               #############################
add_subdirectory(GasFile)

#### Heed                   ############################
add_subdirectory(Heed)

#### SRIM                   ############################
add_subdirectory(Srim)

#### Geant4GarfieldInterface ###########################
add_subdirectory(Geant4GarfieldInterface)

