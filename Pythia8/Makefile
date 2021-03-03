
ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = c++
endif

CXXFLAGS     = -std=c++11
CXXFLAGS     += `root-config --cflags`
CXXFLAGS     += -I$(PYTHIA8)/include

LDFLAGS      = 
LDFLAGS     += `root-config --libs`
LDFLAGS     += -L$(PYTHIA8)/lib -lpythia8

#all : PythiaHepMC

all: runPhoto.exe runDIS.exe runCrossing.exe

HEPMC_DIR = $(EICDIRECTORY)

PythiaHepMC : PythiaHepMC.o
	$(CXX) PythiaHepMC.o -o PythiaHepMC $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaHepMC.o : PythiaHepMC.cxx CustomPythia8ToHepMC3.h
	$(CXX) $(CXXFLAGS) -c PythiaHepMC.cxx -o PythiaHepMC.o -I$(HEPMC_DIR)/include -I.

runPhoto.exe: PythiaPhotoProduction.o
	$(CXX) PythiaPhotoProduction.o -o runPhoto.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaPhotoProduction.o : PythiaPhotoProduction.cxx CustomPythia8ToHepMC3.h
	$(CXX) $(CXXFLAGS) -c PythiaPhotoProduction.cxx -o PythiaPhotoProduction.o -I$(HEPMC_DIR)/include -I.

runDIS.exe: PythiaDIS.o
	$(CXX) PythiaDIS.o -o runDIS.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaDIS.o : PythiaDIS.cxx CustomPythia8ToHepMC3.h
	$(CXX) $(CXXFLAGS) -c PythiaDIS.cxx -o PythiaDIS.o -I$(HEPMC_DIR)/include -I.

runCrossing.exe: PythiaXAngleDIS.o
	$(CXX) PythiaXAngleDIS.o -o runCrossing.exe $(LDFLAGS) -L$(HEPMC_DIR)/lib -L$(HEPMC_DIR)/lib64 $(LDFLAGS) -lHepMC3

PythiaXAngleDIS.o : PythiaXAngleDIS.cxx CustomPythia8ToHepMC3.h
	$(CXX) $(CXXFLAGS) -c PythiaXAngleDIS.cxx -o PythiaXAngleDIS.o -I$(HEPMC_DIR)/include -I.

clean :
	rm -vf PythiaHepMC *.o *.exe *~
