// Modification of Pythia Class BeamShape for EIC Beam Conditions

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "eicBeamShape.h"

#include "TMath.h"

using namespace Pythia8;

eicBeamShape::eicBeamShape(double ion, double lepton, double xAngle, double hadCrab, double lepCrab) 
  {
    mIonBeamEnergy = ion;
    mLeptonBeamEnergy = lepton;
    mXAngle = xAngle;
    mHadronCrabSize = hadCrab;
    mLeptonCrabSize = lepCrab;
  }

void eicBeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Change in Beam Pz
  double tmpPzA, tmpPzB;
  tmpPzA = tmpPzB = 0.;

  // Set beam vertex location by a two-dimensional Gaussian.
  if(allowVertexSpread) 
    {
      if(sigmaVertexZ > 0.) 
	{
	  vertexZ = sigmaVertexZ * rndmPtr->gauss();
	}
    }

  // Smear Beam Energy to Set Scale of Other Effects
  if(allowMomentumSpread)
    {
      double gaussZA, gaussZB;
      gaussZA = gaussZB = 0.;
      if(sigmaPzA > 0.)
	{
	  gaussZA = rndmPtr->gauss();
	  tmpPzA = mIonBeamEnergy * (sigmaPzA * gaussZA);
	}
      if(sigmaPzB > 0.)
	{
	  gaussZB = rndmPtr->gauss();
	  tmpPzB = mLeptonBeamEnergy * (sigmaPzB * gaussZB);
	}

      // Modify Px Due to Crossing Angle
      deltaPxA += (mIonBeamEnergy + tmpPzA)*TMath::Sin(mXAngle);

      // Modify Pz Due to Crossing Angle
      deltaPzA += (mIonBeamEnergy + tmpPzA)*TMath::Cos(mXAngle) - mIonBeamEnergy;
      //deltaPzA += tmpPzA;

      // Electron Beam Has 0 Crossing Angle
      deltaPzB += tmpPzB;
    }

  if(allowMomentumSpread) // allowMomentumSpread
    {
      double gaussXA, gaussYA, gaussXB, gaussYB;

      // Implement Z-Dependent Momentum Kicks from Crabbing
      deltaPxA += mHadronCrabSize*vertexZ;
      deltaPxB += mLeptonCrabSize*vertexZ;

      if(sigmaPxA > 0.)
	{
	  // Ion Beam X Divergence
	  gaussXA = rndmPtr->gauss();
	  double div = sigmaPxA * gaussXA;
	  double pxLocal = (mIonBeamEnergy + tmpPzA)*TMath::Sin(div); // Dispersion in Beam Frame
	  deltaPxA += pxLocal*TMath::Cos(mXAngle); // Compensate for Crossing Angle to Lab Frame
	}
      if(sigmaPyA > 0.)
	{
	  // Ion Beam Y Divergence
	  gaussYA = rndmPtr->gauss();
	  double div = sigmaPyA * gaussYA;
	  deltaPyA += (mIonBeamEnergy + tmpPzA)*TMath::Sin(div); // No Crossing Angle in Y direction
	}
      if(sigmaPxB > 0.)
	{
	  // Lepton Beam X Divergence
	  gaussXB = rndmPtr->gauss();
	  double div = sigmaPxB * gaussXB;
	  deltaPxB += (mLeptonBeamEnergy + tmpPzB)*TMath::Sin(div); // Electron Beam along Z Axis
	}

      if(sigmaPyB > 0.)
	{
	  // Lepton Beam Y Divergence
	  gaussYB = rndmPtr->gauss();
	  double div = sigmaPyB * gaussYB;
	  deltaPyB += (mLeptonBeamEnergy + tmpPzB)*TMath::Sin(div);
	}
    }
}
