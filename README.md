# Garfield++

Garfield++ is a toolkit for the detailed simulation of detectors which use gases or semi-conductors as sensitive medium. The main area of application is currently in micropattern gaseous detectors.

Garfield++ shares functionality with [Garfield](http://cern.ch/garfield). The main differences are the more up-to-date treatment of electron transport in gases and the user interface, which is derived from ROOT.

[More...](http://garfieldpp.web.cern.ch/garfieldpp/about)

## Building the project
For simplicity, we define an environment variable `$GARFIELD_HOME` 
pointing to the directory to which we cloned the repository. 
Assuming that `$GARFIELD_HOME` is your current working directory,
you can build and install Garfield++ as follows: 

    mkdir build
    cd build
    cmake [-DCMAKE_INSTALL_PREFIX=<installdir>] [-DWITH_DOCS=ON] [-DWITH_EXAMPLES=ON] <path to sources>
    make -j<number of cores on your machine>
    make install

If `CMAKE_INSTALL_PREFIX` is not provided in the configuration command, `$GARFIELD_HOME/install` will be used as installation prefix. The `WITH_DOCS` variable is optional, and should be passed if you wish to
build the Doxygen based API documentation. Please note that this requires an existing
installation of [Doxygen](http://www.doxygen.org/index.html). If CMake cannot locate
Doxygen, its install location should be added into `CMAKE_PREFIX_PATH`.
For further details please have a look at [the CMake tutorial](http://www.cmake.org/cmake-tutorial/).

## Building the documentation

The documentation of the project is based on doxygen. To build the documentation,
the project must have been configured with `Garfield_BUILD_DOCS` enabled, as
described earlier. It can then be built and installed:

    make doc
    make install

By default, this installs the documentation into `<installdir>/share/doc/HSFTEMPLATE/share/doc`.

## Running the examples

You can run the examples from the build directory (if WITH_EXAMPLES has been turned on) but you need to setup a running environment defining some variables, in particular for the HEED database.

In the following lines we use the `Gem` example:
```
source setupGarfield.sh
cd Examples/Gem
./gem
```
## Building and running examples using an installed version of Garfield

Make sure that all required environment variables are set by sourcing the script `setupGarfield.sh`:
```
source $GARFIELD_HOME/install/share/Garfield/setupGarfield.sh
```

To get started, it can be useful to copy one of the examples to 
a local directory, modify it, build it against an installed version of Garfield and run it. 
```
cp -r $GARFIELD_HOME/Examples/Gem .
mkdir Gem/build; cd Gem/build
cmake ..
make
./gem
```  

## Inclusion into other projects

If you want to build your own project against Garfield, CMake may be the best option for you. Just add its location to _CMAKE_PREFIX_PATH_ and call _find_package(Garfield)_ within your CMakeLists.txt.

```
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(test)

find_package(Garfield REQUIRED)

add_executable(test test.C)
target_link_libraries(test Garfield)
```

## About Windows Build:
Currently only tested with cmake+ifort. gfortran should use original one. Also, cmake by default adds /fpp to CMAKE_Fortran_FLAGS, 
which will break build. One should remove that option from cmake cache or specify correct options.
Also requires to manually add the following definitions:
`/D__STDC__ /D_USE_MATH_DEFINES /D_CRT_SECURE_NO_DEPRECATE`
