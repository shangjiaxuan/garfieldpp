INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

plotfield: plotfield.C 
	$(CXX) $(CFLAGS) plotfield.C
	$(CXX) -o plotfield plotfield.o $(LDFLAGS)
	rm plotfield.o

savefield: savefield.C 
	$(CXX) $(CFLAGS) savefield.C
	$(CXX) -o savefield savefield.o $(LDFLAGS)
	rm savefield.o

