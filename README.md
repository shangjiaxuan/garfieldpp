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

In the following lines we use the `Gem` example (assuming that the build folder is our current working directory):
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
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(test)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(test test.C)
target_link_libraries(test Garfield::Garfield)
```

## About Windows Build:
Currently only tested 64 bit with cmake visual studio msvc + ifort 
(from [`Intel oneAPI HPC Toolkit`](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html), 
need registering and EULA to download) and msvc + gfortran (mingw binaries), using root from public release binary, 
and self-compiled geant4 with Qt 5.15 and GDML (xerces-c installed with vcpkg).

The compiled Garfield++ library have root generated header paths that loads on start, thus source directory needs to be
constant, and cannot make portable install binary. Dependency binaries can be found in the releases section (excluding Garfield.dll).

Since the library requires large memory resources, 64 bit compilation is required. Running GEM example uses 2.1GB on typical 
machine:
![GEM example](GEM-win64.PNG)

The following environment variables need to be set to run (running without the shell scripts):
```
GARFIELD_INSTALL (Where to look for gas datasets)
optional:
    HEED_DATABASE (By default a location relative to GARFIELD_INSTALL)
    G4XXXXDATA (The Geant4 data variables for datasets)
```

### Required dependencies:
#### Root
vs2022 binary available at https://root.cern/install/all_releases/, installer version sets environment variables. 
Tested with 6.26/06, should work with all versions built with vs2022.
#### Geant4
vs2019 binary available at https://geant4.web.cern.ch/support/download, binary compatible with vs2022.
Tested with self-compiled 6.11.03 with Qt and GDML(xerces-c) using vs2022. 
Qt5.15 binary in vs2019 is compatible with vs2022 compiler.
#### gsl
`vckpg install gsl`, and just use cmake with vcpkg. Currently you need to specify the gsl and gslblas library files manually. 
Located at `<vcpkg-dir>\installed\<compiler-setting(x64-windows)>`
