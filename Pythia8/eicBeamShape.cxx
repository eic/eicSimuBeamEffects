// Modification of Pythia Class BeamShape for EIC Beam Conditions

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "eicBeamShape.h"

#include "TMath.h"

using namespace Pythia8;

eicBeamShape::eicBeamShape(double ion, double lepton, double xAngle) 
  {
    mIonBeamEnergy = ion;
    mLeptonBeamEnergy = lepton;
    mXAngle = xAngle;
  }

void eicBeamShape::pick() {

  //=============================================
  // (Re)set Values
  //=============================================

  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;


  //=============================================
  // Set Collision Vertex and Time
  //=============================================

  // Collision vertex and time are determined by randomly sampling the z-positions of the interacting lepton and hadron within their respective bunches. Assume z=0 is middle of bunch, hadron bunch moves left to right and lepton bunch moves right to left, positive z within hadron bunch is in direction of travel and positive z within lepton bunch is opposite direction of travel

  //double xIPOffset = 0.;
  double hadronPartPos = 0.;
  double leptonPartPos = 0.;

  if(allowVertexSpread)
    {
      // RMS Bunch Length [mm]
      double hadronBL = 0.;
      double leptonBL = 0.;
      
      // Set different values depending on energy
      if(mIonBeamEnergy == 275.0) hadronBL = 60.0;
      if(mLeptonBeamEnergy == 18.0) leptonBL = 9.0;
      
      // Set particle positions
      hadronPartPos = hadronBL*rndmPtr->gauss();
      leptonPartPos = leptonBL*rndmPtr->gauss();

      // Find Collision time, z position of collision, z position of center of hadron bunch at collision time and x position of collision
      // Quantities are found via a system of parametric equations
      // Z position of colliding hadron as a function of time: Z_h = Cos(0.5*theta_c)*t + Z_h_offset (Z_h_offset = distance from center of bunch)
      // Z position of colliding lepton as a function of time: Z_l = -Cos(0.5*theta_c)*t + Z_l_offset
      // Collision time when Z_h = Z_l -> t_int = (Z_l_offset - Z_h_offset)/(2*Cos(0.5*theta_c))
      // Z position of collision: Z_int = (Z_l_offset + Z_h_offset)/2
      // Z position of center of hadron bunch a collision time - this gives x position of collision via relation x = z*Tan(0.5*theta_c)
      // Z_bunch_int = Cos(0.5*theta_c)*t_int
      // X_int = Z_bunch_int * Tan(0.5*theta_c)

      double c_c = TMath::Cos(mXAngle/2.0);
      double s_c = TMath::Sin(mXAngle/2.0);
      double t_c = TMath::Tan(mXAngle/2.0);

      double t_int = (leptonPartPos - hadronPartPos)/(2.0*c_c);
      double z_int = (leptonPartPos + hadronPartPos)/2.0;
      double z_bunch_int = c_c*t_int;
      double x_int = z_bunch_int*t_c;

      // x_int is the x position of collision assuming the colliding particles are in the center of the bunch. Sample random x position according to x-width of bunches. x_int is then an offset to this. Get y position as well.

      double y_int = 0.;
      if(sigmaVertexX > 0.)
	{
	  x_int += sigmaVertexX * rndmPtr->gauss();
	}

      if(sigmaVertexY > 0.)
	{
	  y_int += sigmaVertexY * rndmPtr->gauss();
	}

      // We now have the x-y-z position of the collision in the accelerator frame, but we want it in the detector frame. Rotate by 0.5*theta_c to get to accelerator frame

      double tmpVtxX, tmpVtxY, tmpVtxZ;
      tmpVtxX = tmpVtxY = tmpVtxZ = 0.;

      RotY(mXAngle/2.0,x_int,y_int,z_int,&tmpVtxX,&tmpVtxY,&tmpVtxZ);

      vertexT = t_int;
      vertexX = tmpVtxX;
      vertexY = tmpVtxY;
      vertexZ = tmpVtxZ;
    }


  //=============================================
  // Smear Beam Energy
  //=============================================

  // Change in Beam Pz
  double tmpPzA, tmpPzB;
  tmpPzA = tmpPzB = 0.;

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
    }


  //=============================================
  // Crossing Angle Effects
  //=============================================

  if(allowMomentumSpread)
    {
      // Modify Px Due to Crossing Angle
      deltaPxA += (mIonBeamEnergy + tmpPzA)*TMath::Sin(mXAngle);

      // Modify Pz Due to Crossing Angle
      deltaPzA += (mIonBeamEnergy + tmpPzA)*TMath::Cos(mXAngle) - mIonBeamEnergy;
      //deltaPzA += tmpPzA;

      // Electron Beam Has 0 Crossing Angle
      deltaPzB += tmpPzB;
    }


  //=============================================
  // Divergence Effects
  //=============================================

  if(allowMomentumSpread) // allowMomentumSpread
    {
      double gaussXA, gaussYA, gaussXB, gaussYB;

      if(sigmaPxA > 0.)
	{
	  // Ion Beam X Divergence
	  gaussXA = rndmPtr->gauss();
	  double div = sigmaPxA * gaussXA;
	  double pxLocal = (mIonBeamEnergy + tmpPzA)*TMath::Sin(div); // Dispersion in Beam Frame

	  // Rotate into Detector Frame
	  double divPxA, divPyA, divPzA;
	  divPxA = divPyA = divPzA = 0.;

	  RotY(mXAngle,pxLocal,0.,0.,&divPxA,&divPyA,&divPzA);

	  deltaPxA += divPxA;
	  deltaPzA += divPzA;
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


  //=============================================
  // Crabbing Momentum Kicks
  //=============================================

  if(allowMomentumSpread)
    {
      double betaCrabHad = 1300000.0; // Beta Crab in mm for 275 hadron beam
      double betaStarHad = 800.0; // Beta Star in mm for 275 hadron beam

      double betaCrabLep = 150000.0; // Beta Crab in mm for 18 lepton beam
      double betaStarLep = 590.0; // Beta Star in mm for 18 lepton beam

      // Calculate angular deflection
      double crabAngHad = ((mXAngle/2.0)*hadronPartPos)/(TMath::Sqrt(betaCrabHad*betaStarHad));
      double crabAngLep = ((mXAngle/2.0)*leptonPartPos)/(TMath::Sqrt(betaCrabLep*betaStarLep));

      // Calculate Magnitude of Px kick
      double crabKickHad = (mIonBeamEnergy + tmpPzA)*TMath::Sin(crabAngHad);
      double crabKickLep = (-1.0)*(mLeptonBeamEnergy + tmpPzB)*TMath::Sin(crabAngLep); // 

      // Rotate Momentum Kick into Detector Frame
      double tmpVertPxA, tmpVertPyA, tmpVertPzA;
      double tmpVertPxB, tmpVertPyB, tmpVertPzB;
      tmpVertPxA = tmpVertPyA = tmpVertPzA = tmpVertPxB = tmpVertPyB = tmpVertPzB = 0.;

      RotY(mXAngle/2.0,crabKickHad,0.,0.,&tmpVertPxA,&tmpVertPyA,&tmpVertPzA);
      RotY(mXAngle/2.0,crabKickLep,0.,0.,&tmpVertPxB,&tmpVertPyB,&tmpVertPzB);

      deltaPxA += tmpVertPxA;
      deltaPzA += tmpVertPzA;
      deltaPxB += tmpVertPxB;
      deltaPzB += tmpVertPzB;
    }
}


void eicBeamShape::RotY(double theta, double xin, double yin, double zin, double *xout, double *yout, double *zout) {

  *xout = xin*TMath::Cos(theta) + zin*TMath::Sin(theta);
  *yout = yin;
  *zout = zin*TMath::Cos(theta) - xin*TMath::Sin(theta);
}
