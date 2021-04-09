This directory conains methods for implementing various beam effects which 
will be present at the EIC (crossing angle, divergence, crabbing kicks, etc) 
in Pythia8 using the native BeamShape functionality. 

eicBeamShape.h/cxx: Definition and implementation of beam effects

PythiaBeamShape.cxx: Main code to run Pythia8 with beam effects implemented and generate various diagnostic plots in ROOT

PythiaBeamShapeHepMC.cxx: Code to run Pythia8 and pass output as HepMC file. No analysis is performed.

dis.steer: Steering file to set Pythia8 configuration