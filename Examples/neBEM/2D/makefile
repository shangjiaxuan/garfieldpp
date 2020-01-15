INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

triangle: triangle.C 
	$(CXX) $(CFLAGS) triangle.C
	$(CXX) -o triangle triangle.o $(LDFLAGS)
	rm triangle.o

dielectric: dielectric.C 
	$(CXX) $(CFLAGS) dielectric.C
	$(CXX) -o dielectric dielectric.o $(LDFLAGS)
	rm dielectric.o

