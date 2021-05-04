This directory conains methods for implementing various beam effects which 
will be present at the EIC (crossing angle, divergence, crabbing kicks, etc) 
in Pythia8 using the native BeamShape functionality. 

eicBeamShape.h/cxx: Definition and implementation of beam effects

PythiaBeamShape.cxx: Main code to run Pythia8 with beam effects implemented and generate various diagnostic plots in ROOT

PythiaBeamShapeHepMC.cxx: Code to run Pythia8 and pass output as HepMC file. No analysis is performed.

dis.steer: Dummy steering file to set Pythia8 configuration when running tests

steerFiles/: Directory containing steering files for the beam energy combinations and configurations (high divergence or high acceptance) listed in the CDR


To run PythiaBeamShape.cxx:

make

./runBeamShape.cxx steerFile configFlag hadEnergy lepEnergy x-Angle out.hist.root