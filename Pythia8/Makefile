
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

CXXFLAGS     = -std=c++11
CXXFLAGS     += `root-config --cflags`
CXXFLAGS     += `fastjet-config --cxxflags`
CXXFLAGS     += -I$(PYTHIA8)/include

LDFLAGS      = 
LDFLAGS     += `root-config --libs`
LDFLAGS     += `fastjet-config --libs --plugins`
LDFLAGS     += -L$(PYTHIA8)/lib -lpythia8

#all : PythiaHepMC

all: runBeamShape.exe runBeamShapeHepMC.exe

HEPMC_DIR = $(EICDIRECTORY)

#PythiaHepMC : PythiaHepMC.o
#	$(CXX) PythiaHepMC.o -o PythiaHepMC $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HE#PMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaHepMC.o : PythiaHepMC.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaHepMC.cxx -o PythiaHepMC.o -I$(HEPMC_DIR)/include -I.

#runPhoto.exe: PythiaPhotoProduction.o
#	$(CXX) PythiaPhotoProduction.o -o runPhoto.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaPhotoProduction.o : PythiaPhotoProduction.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaPhotoProduction.cxx -o PythiaPhotoProduction.o -I$(HEPMC_DIR)/include -I.

#runDIS.exe: PythiaDIS.o
#	$(CXX) PythiaDIS.o -o runDIS.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

#PythiaDIS.o : PythiaDIS.cxx CustomPythia8ToHepMC3.h
#	$(CXX) $(CXXFLAGS) -c PythiaDIS.cxx -o PythiaDIS.o -I$(HEPMC_DIR)/include -I.

runBeamShape.exe: PythiaBeamShape.o eicBeamShape.o
	$(CXX) PythiaBeamShape.o eicBeamShape.o -o runBeamShape.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaBeamShape.o : PythiaBeamShape.cxx eicBeamShape.h
	$(CXX) $(CXXFLAGS) -c PythiaBeamShape.cxx -o PythiaBeamShape.o -I$(HEPMC_DIR)/include -I.

runBeamShapeHepMC.exe: PythiaBeamShapeHepMC.o eicBeamShape.o
	$(CXX) PythiaBeamShapeHepMC.o eicBeamShape.o -o runBeamShapeHepMC.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaBeamShapeHepMC.o : PythiaBeamShapeHepMC.cxx eicBeamShape.o
	$(CXX) $(CXXFLAGS) -c PythiaBeamShapeHepMC.cxx -o PythiaBeamShapeHepMC.o -I$(HEPMC_DIR)/include -I.

eicBeamShape.o: eicBeamShape.cxx
	$(CXX) $(CXXFLAGS) -c eicBeamShape.cxx -o eicBeamShape.o

clean :
	rm -vf *.o *.exe *~
