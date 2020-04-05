INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library

# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
	`root-config --cflags` \
	-O3 -fno-common -c \
	-I$(INCDIR)

LDFLAGS = `root-config --glibs` -lGeom -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield

full_tpc_sim:	full_tpc_sim.C
	$(CXX) $(CFLAGS) full_tpc_sim.C
	$(CXX) -o full_tpc_sim full_tpc_sim.o $(LDFLAGS)
	rm full_tpc_sim.o

ion_signal:	ion_signal.C
	$(CXX) $(CFLAGS) ion_signal.C
	$(CXX) -o ion_signal ion_signal.o $(LDFLAGS)
	rm ion_signal.o

plot_gas:	plot_gas.C
	$(CXX) $(CFLAGS) plot_gas.C
	$(CXX) -o plot_gas plot_gas.o $(LDFLAGS)
	rm plot_gas.o
