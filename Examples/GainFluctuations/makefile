INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

arco2: arco2.C 
	$(CXX) $(CFLAGS) arco2.C
	$(CXX) `root-config --cflags` -o arco2 arco2.o $(LDFLAGS)
	rm arco2.o
