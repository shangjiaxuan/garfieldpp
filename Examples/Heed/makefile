INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = -L$(LIBDIR) -lGarfield
LDFLAGS += `root-config --glibs` -lGeom -lgfortran -lm

edep: edep.C 
	$(CXX) $(CFLAGS) edep.C
	$(CXX) -o edep edep.o $(LDFLAGS)
	rm edep.o

qdepSi: qdepSi.C 
	$(CXX) $(CFLAGS) qdepSi.C
	$(CXX) -o qdepSi qdepSi.o $(LDFLAGS)
	rm qdepSi.o

rangeSi: rangeSi.C 
	$(CXX) $(CFLAGS) rangeSi.C
	$(CXX) -o rangeSi rangeSi.o $(LDFLAGS)
	rm rangeSi.o

fe55: fe55.C 
	$(CXX) $(CFLAGS) fe55.C
	$(CXX) -o fe55 fe55.o $(LDFLAGS)
	rm fe55.o
