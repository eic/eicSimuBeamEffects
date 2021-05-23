This directory conains methods for implementing various beam effects which 
will be present at the EIC (crossing angle, divergence, crabbing kicks, etc) 
in Pythia8 using the native BeamShape functionality. 

eicBeamShape.h/cxx: Definition and implementation of beam effects

PythiaBeamShape.cxx: Main code to run Pythia8 with beam effects implemented and generate various diagnostic plots in ROOT

PythiaBeamShapeHepMC.cxx: Code to run Pythia8 and pass output as HepMC file. No analysis is performed.

dis.steer: Dummy steering file to set Pythia8 configuration when running tests

steerFiles/: Directory containing steering files for the beam energy combinations and configurations (high divergence or high acceptance) listed in the CDR. Have also added files with the title *_JinHeadOn_* in which divergence and energy spread are turned off and the beam energies are scaled by Cos(xAngle/2). These files can be used to test Jin's afterburner code. Files *_noDivergence_* and *_noEnergySpread_* zero out either divergence or beam energy spread.


To run PythiaBeamShape.cxx:

make

./runBeamShape.cxx steerFile configFlag crabFlag hadEnergy lepEnergy x-Angle out.hist.root out.hepmc

steerFile is one of the files in steeringFiles/

configFlag takes three values: 1 = hiDiv steering file, 2 = hiAcc steering files, 0 = special flag to disable steering file - energy setting consistency (use with JinHeadOn, noDivergence, and noEnergySpread steering files

crabFlag turns off (0) or on (1) the momentum kick from crabbing

hadEnergy, lepEnergy, x-Angle, out.hist.root, and out.hepmc are the hadron and lepton beam energies, crossing angle, output root file an output hepmc file

Output:

Currently, output root and hepmc files live at /gpfs02/eic/bpage/home/eicBeamSimu/Pythia8/headonTestJin . The file names should be self documenting with the second section describing what effects are included. For example, test_crossDivNrgCrab_25mRad_18x275.hist.root was run with a 25 mRad crossing angle and 18x275 beam energy, divergence effects, beam energy smearing, and momentum kicks from crabbing.