// Modification of Pythia Class BeamShape for EIC Beam Conditions

#ifndef EICBEAMSHAPE
#define EICBEAMSHAPE

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"

using namespace Pythia8;

// Class to Define Beam Properties
class eicBeamShape : public BeamShape {

 public:

  // Constructor
  eicBeamShape(double ion, double lepton, double xAngle, double hadCrab, double lepCrab);

  // Set Beam Properties to Pass to Pythia
  virtual void pick();

 protected:

  double mIonBeamEnergy;
  double mLeptonBeamEnergy;
  double mXAngle;
  double mHadronCrabSize;
  double mLeptonCrabSize;

};

#endif
