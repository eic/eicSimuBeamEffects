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
  eicBeamShape(int config, double ion, double lepton, double xAngle);

  // Set Beam Properties to Pass to Pythia
  virtual void pick();

  // Rotation Utility
  void RotY(double theta, double xin, double yin, double zin, double *xout, double *yout, double *zout);

  void RotXY(double theta, double phi, double xin, double yin, double zin, double *xout, double *yout, double *zout);

  // Getters for configuration parameters
  int getDivAcc() const { return mDivAcc; }
  double getIonBeamEnergy() const { return mIonBeamEnergy; }
  double getLeptonBeamEnergy() const { return mLeptonBeamEnergy; }
  double getXAngle() const { return mXAngle; }
  double getXAngleY() const { return mXAngleY; }

 protected:

  int mDivAcc;
  double mIonBeamEnergy;
  double mLeptonBeamEnergy;
  double mXAngle;
  double mXAngleY;
  int mKill;

};

#endif
