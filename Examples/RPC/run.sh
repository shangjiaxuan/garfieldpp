#!/bin/bash

echo "Running Weight script"

cwd=$(pwd)

echo "Current path:"

echo $cwd

export GARFIELD_HOME=/Users/djunesjanssens/Documents/Ph.D./Garfield

cd $GARFIELD_HOME

bash run.sh

cd $cwd

rm -rf build

mkdir build

cd build

cmake ..

make

./Weight
