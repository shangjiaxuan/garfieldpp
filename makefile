OBJDIR = $(GARFIELD_HOME)/Object
SRCDIR = $(GARFIELD_HOME)/Source
INCDIR = $(GARFIELD_HOME)/Include
LIBDIR = $(GARFIELD_HOME)/Library
HEEDDIR = $(GARFIELD_HOME)/Heed

HEADERS = $(wildcard $(INCDIR)/Garfield/*.hh)

SOURCES = $(wildcard $(SRCDIR)/*.cc)

OBJECTS = $(subst $(SRCDIR),$(OBJDIR),$(SOURCES:.cc=.o))
OBJECTS += $(OBJDIR)/magboltz.o 

TARGETS = $(OBJECTS)
TARGETS += heed

# Fortran compiler
FC = gfortran

# Compilation flags
CFLAGS = -std=c++11 -Wall -Wextra -pedantic -ansi -Wabi -Wno-long-long -Woverloaded-virtual -Wshadow \
	 `root-config --cflags` \
        -fpic -fno-common -c \
	-I$(INCDIR) -I$(HEEDDIR) -DINS_CRETURN 

FFLAGS = -fpic -c

# Optimization flags
# CFLAGS += -Os
# FFLAGS += -Os
CFLAGS += -O2
FFLAGS += -O2

# Debug flags
# CFLAGS += -g
# FFLAGS += -g
# Profiling flag
# CFLAGS += -pg

# Linking flags
LDFLAGS = `root-config --glibs` `root-config --ldflags`-lGeom \
	-lgfortran -lm

all:	$(TARGETS)
	@echo Creating library libGarfield...
	@ar rcs $(LIBDIR)/libGarfield.a $(OBJECTS) \
	$(wildcard $(OBJDIR)/Heed/*.o)
	@ranlib $(LIBDIR)/libGarfield.a
	@touch $(OBJDIR)/last_updated_on
	@echo Finished.

.PHONY:	heed

installdirs : 
	@if [ ! -d $(GARFIELD_HOME)/Library/$(BFARCH) ] ; then mkdir -p $(GARFIELD_HOME)/Library/$(BFARCH); \
	    echo "   >>>> Create $(GARFIELD_HOME)/Library/$(BFARCH)"; fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); \
	    echo "   >>>> Create $(OBJDIR)"; \
        else echo " $(OBJDIR) already exists"; fi
	@if [ ! -d $(OBJDIR)/Heed/ ]; then mkdir -p $(OBJDIR)/Heed; \
	    echo "   >>>> Create $(OBJDIR)/Heed"; \
        else echo " $(OBJDIR)/Heed already exists"; fi
        
heed:	
	@echo Compiling Heed...
	@cd $(HEEDDIR); make; cd $(GARFIELD_HOME)
	@touch $(OBJDIR)/last_updated_on

clean:
	@echo Removing object files...
	@$(RM) $(OBJDIR)/*.o
	@echo Removing libraries...
	@$(RM) $(LIBDIR)/*.a
	@cd $(HEEDDIR); make clean
	@echo Removing dictionary...
	@$(RM) $(SRCDIR)/GarfieldDict.C 

$(OBJDIR)/AvalancheMicroscopic.o: \
	$(SRCDIR)/AvalancheMicroscopic.cc \
	$(INCDIR)/Garfield/AvalancheMicroscopic.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh $(INCDIR)/Garfield/GarfieldConstants.hh \
	$(INCDIR)/Garfield/Random.hh \
	$(INCDIR)/Garfield/Sensor.hh $(INCDIR)/Garfield/Medium.hh $(INCDIR)/Garfield/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/AvalancheMC.o: \
	$(SRCDIR)/AvalancheMC.cc $(INCDIR)/Garfield/AvalancheMC.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh $(INCDIR)/Garfield/GarfieldConstants.hh \
	$(INCDIR)/Garfield/Random.hh \
	$(INCDIR)/Garfield/Sensor.hh $(INCDIR)/Garfield/Medium.hh $(INCDIR)/Garfield/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@      
$(OBJDIR)/DriftLineRKF.o: \
	$(SRCDIR)/DriftLineRKF.cc $(INCDIR)/Garfield/DriftLineRKF.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh \
	$(INCDIR)/Garfield/Sensor.hh $(INCDIR)/Garfield/Medium.hh $(INCDIR)/Garfield/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
 
$(OBJDIR)/Track.o: \
	$(SRCDIR)/Track.cc $(INCDIR)/Garfield/Track.hh \
	$(INCDIR)/Garfield/ViewDrift.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/TrackBichsel.o: \
	$(SRCDIR)/TrackBichsel.cc $(INCDIR)/Garfield/TrackBichsel.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@       
$(OBJDIR)/TrackPAI.o: \
	$(SRCDIR)/TrackPAI.cc $(INCDIR)/Garfield/TrackPAI.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/TrackSimple.o: \
	$(SRCDIR)/TrackSimple.cc $(INCDIR)/Garfield/TrackSimple.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/TrackHeed.o: \
	$(SRCDIR)/TrackHeed.cc $(INCDIR)/Garfield/TrackHeed.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc \
	$(HEEDDIR)/HeedChamber.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/TrackElectron.o: \
	$(SRCDIR)/TrackElectron.cc $(INCDIR)/Garfield/TrackElectron.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/TrackSrim.o: \
	$(SRCDIR)/TrackSrim.cc $(INCDIR)/Garfield/TrackSrim.hh \
	$(INCDIR)/Garfield/Track.hh $(SRCDIR)/Track.cc 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/ComponentBase.o: \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh \
	$(INCDIR)/Garfield/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentConstant.o: \
	$(SRCDIR)/ComponentConstant.cc $(INCDIR)/Garfield/ComponentConstant.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentUser.o: \
	$(SRCDIR)/ComponentUser.cc $(INCDIR)/Garfield/ComponentUser.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@       
$(OBJDIR)/ComponentAnalyticField.o: \
	$(SRCDIR)/ComponentAnalyticField.cc \
	$(INCDIR)/Garfield/ComponentAnalyticField.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentNeBem2d.o: \
	$(SRCDIR)/ComponentNeBem2d.cc $(INCDIR)/Garfield/ComponentNeBem2d.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/ComponentNeBem3d.o: \
	$(SRCDIR)/ComponentNeBem3d.cc $(INCDIR)/Garfield/ComponentNeBem3d.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/ComponentNeBem3dMap.o: \
	$(SRCDIR)/ComponentNeBem3dMap.cc $(INCDIR)/Garfield/ComponentNeBem3dMap.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        
$(OBJDIR)/ComponentFieldMap.o: \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/Garfield/ComponentFieldMap.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentAnsys121.o: \
	$(SRCDIR)/ComponentAnsys121.cc $(INCDIR)/Garfield/ComponentAnsys121.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/Garfield/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentAnsys123.o: \
	$(SRCDIR)/ComponentAnsys123.cc $(INCDIR)/Garfield/ComponentAnsys123.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/Garfield/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentCST.o: \
	$(SRCDIR)/ComponentCST.cc $(INCDIR)/Garfield/ComponentCST.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/Garfield/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentElmer.o: \
	$(SRCDIR)/ComponentElmer.cc $(INCDIR)/Garfield/ComponentElmer.hh \
	$(SRCDIR)/ComponentFieldMap.cc $(INCDIR)/Garfield/ComponentFieldMap.hh 
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentComsol.o: \
	$(SRCDIR)/ComponentComsol.cc $(INCDIR)/Garfield/ComponentComsol.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ComponentTcad2d.o: \
	$(SRCDIR)/ComponentTcad2d.cc $(INCDIR)/Garfield/ComponentTcad2d.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentTcad3d.o: \
	$(SRCDIR)/ComponentTcad3d.cc $(INCDIR)/Garfield/ComponentTcad3d.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentVoxel.o: \
	$(SRCDIR)/ComponentVoxel.cc $(INCDIR)/Garfield/ComponentVoxel.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/ComponentUserMapBase.o: \
	$(SRCDIR)/ComponentUserMapBase.cc $(INCDIR)/Garfield/ComponentUserMapBase.hh \
	$(SRCDIR)/ComponentBase.cc $(INCDIR)/Garfield/ComponentBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    
	
$(OBJDIR)/GeometrySimple.o: \
	$(SRCDIR)/GeometrySimple.cc $(INCDIR)/Garfield/GeometrySimple.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@   
$(OBJDIR)/GeometryRoot.o: \
	$(SRCDIR)/GeometryRoot.cc $(INCDIR)/Garfield/GeometryRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@    

$(OBJDIR)/ViewBase.o: \
	$(SRCDIR)/ViewBase.cc $(INCDIR)/Garfield/ViewBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewFEMesh.o: \
	$(SRCDIR)/ViewFEMesh.cc $(INCDIR)/Garfield/ViewFEMesh.hh \
	$(SRCDIR)/ViewDrift.cc $(INCDIR)/Garfield/ViewDrift.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/ComponentFieldMap.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewField.o: \
	$(SRCDIR)/ViewField.cc $(INCDIR)/Garfield/ViewField.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewDrift.o: \
	$(SRCDIR)/ViewDrift.cc $(INCDIR)/Garfield/ViewDrift.hh \
	$(INCDIR)/Garfield/ViewBase.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewMedium.o: \
	$(SRCDIR)/ViewMedium.cc $(INCDIR)/Garfield/ViewMedium.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewSignal.o: \
	$(SRCDIR)/ViewSignal.cc $(INCDIR)/Garfield/ViewSignal.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/Sensor.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewCell.o: \
	$(SRCDIR)/ViewCell.cc $(INCDIR)/Garfield/ViewCell.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/ComponentAnalyticField.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/ViewGeometry.o: \
	$(SRCDIR)/ViewGeometry.cc $(INCDIR)/Garfield/ViewGeometry.hh \
	$(INCDIR)/Garfield/ViewBase.hh \
	$(INCDIR)/Garfield/GeometrySimple.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Medium.o: \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh \
	$(INCDIR)/Garfield/Numerics.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumGas.o: \
	$(SRCDIR)/MediumGas.cc $(INCDIR)/Garfield/MediumGas.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumMagboltz.o: \
	$(SRCDIR)/MediumMagboltz.cc $(INCDIR)/Garfield/MediumMagboltz.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh $(SRCDIR)/OpticalData.cc \
	$(SRCDIR)/MediumGas.cc $(INCDIR)/Garfield/MediumGas.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh $(INCDIR)/Garfield/Random.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/magboltz.o: \
	$(SRCDIR)/magboltz-11.9.f
	@echo $@
	@$(FC) $(FFLAGS) $< -o $@
$(OBJDIR)/MediumSilicon.o: \
	$(SRCDIR)/MediumSilicon.cc $(INCDIR)/Garfield/MediumSilicon.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumCdTe.o: \
	$(SRCDIR)/MediumCdTe.cc $(INCDIR)/Garfield/MediumCdTe.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/MediumGaAs.o: \
	$(SRCDIR)/MediumGaAs.cc $(INCDIR)/Garfield/MediumGaAs.hh \
	$(SRCDIR)/Medium.cc $(INCDIR)/Garfield/Medium.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/OpticalData.o: \
	$(SRCDIR)/OpticalData.cc $(INCDIR)/Garfield/OpticalData.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Solid.o: \
	$(SRCDIR)/Solid.cc $(INCDIR)/Garfield/Solid.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidBox.o: \
	$(SRCDIR)/SolidBox.cc $(INCDIR)/Garfield/SolidBox.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidTube.o: \
	$(SRCDIR)/SolidTube.cc $(INCDIR)/Garfield/SolidTube.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidSphere.o: \
	$(SRCDIR)/SolidSphere.cc $(INCDIR)/Garfield/SolidSphere.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidHole.o: \
	$(SRCDIR)/SolidHole.cc $(INCDIR)/Garfield/SolidHole.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/SolidRidge.o: \
	$(SRCDIR)/SolidRidge.cc $(INCDIR)/Garfield/SolidRidge.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Random.o: \
	$(SRCDIR)/Random.cc $(INCDIR)/Garfield/Random.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  
$(OBJDIR)/RandomEngineGSL.o: \
	$(SRCDIR)/RandomEngineGSL.cc $(INCDIR)/Garfield/RandomEngineGSL.hh
	$(CXX) $(CFLAGS) $< -o $@
$(OBJDIR)/RandomEngineRoot.o: \
	$(SRCDIR)/RandomEngineRoot.cc $(INCDIR)/Garfield/RandomEngineRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/PlottingEngineRoot.o: \
	$(SRCDIR)/PlottingEngineRoot.cc $(INCDIR)/Garfield/PlottingEngineRoot.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@        

$(OBJDIR)/Numerics.o: \
	$(SRCDIR)/Numerics.cc $(INCDIR)/Garfield/Numerics.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@  

$(OBJDIR)/Sensor.o: \
	$(SRCDIR)/Sensor.cc $(INCDIR)/Garfield/Sensor.hh \
	$(INCDIR)/Garfield/ComponentBase.hh $(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/Shaper.o: \
	$(SRCDIR)/Shaper.cc $(INCDIR)/Garfield/Shaper.hh \
	$(INCDIR)/Garfield/FundamentalConstants.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/TetrahedralTree.o: \
	$(SRCDIR)/TetrahedralTree.cc $(INCDIR)/Garfield/TetrahedralTree.hh
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@

$(OBJDIR)/GarfieldDict.o: \
	$(SRCDIR)/GarfieldDict.C
	@echo $@
	@$(CXX) $(CFLAGS) -DDICT_SKIP_HEED $< -o $@

$(SRCDIR)/GarfieldDict.C: $(HEADERS) $(INCDIR)/Garfield/LinkDef.h
	@echo Creating dictionary...
	@rootcint -f $@ -c $(CFLAGS) -p $^ 
