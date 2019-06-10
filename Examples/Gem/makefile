INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

gem: gem.C 
	$(CXX) $(CFLAGS) gem.C
	$(CXX) `root-config --cflags` -o gem gem.o $(LDFLAGS)
	rm gem.o
