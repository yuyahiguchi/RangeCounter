CXX = g++ 
ROOTGLIBS = $(shell root-config --cflags --glibs)
CXXLIBS    = $(ROOTGLIBS)

INCLUDE = include/functions.h include/ParticleID.h include/DetectorClass.h include/VirtualDetectorClass.h include/NewParticleClass.h include/BeamLossClass.h
CFLAGS = -Wall
DEST = /usr/local/bin
LDFLAGS = -L/usr/local/lib

all : DegraderLengthAna TrueAna CompressRootFile

DegraderLengthAna : src/DegraderLengthAna.C $(INCLUDE)
	$(CXX) $(CXXLIBS) -o $@ $<

TrueAna : src/TrueAna.C $(INCLUDE)
	$(CXX) $(CXXLIBS) -o $@ $<

CompressRootFile : src/CompressRootFile.C $(INCLUDE)
	$(CXX) $(CXXLIBS) -o $@ $<	

clean:
	rm -f DegraderLengthAna TrueAna CompressRootFile
