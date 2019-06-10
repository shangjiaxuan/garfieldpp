INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield

avalanche:	avalanche.cc
	$(CXX) $(CFLAGS) avalanche.cc
	$(CXX) -o avalanche avalanche.o $(LDFLAGS)
	rm avalanche.o
