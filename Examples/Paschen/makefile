gastable: gastable.C
	g++ -c -O3 `root-config --cflags` -I$(GARFIELD_HOME)/Include gastable.C
	g++ `root-config --glibs` -lgfortran -lm -o gastable gastable.o $(GARFIELD_HOME)/Library/libGarfield.a
	rm gastable.o
paschenariso: paschenariso.C
	g++ -c -O3 `root-config --cflags` -I$(GARFIELD_HOME)/Include paschenariso.C
	g++ `root-config --glibs` -lgfortran -lm -o paschenariso paschenariso.o $(GARFIELD_HOME)/Library/libGarfield.a
	rm paschenariso.o
