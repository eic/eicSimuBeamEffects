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

  if(allowVertexSpread)
    {
      // RMS Bunch Length [mm]
      double hadronBL = 0.;
      double leptonBL = 0.;
      
      // Set different values depending on energy
      if(mIonBeamEnergy == 275.0) hadronBL = 60.0;
      if(mLeptonBeamEnergy == 18.0) leptonBL = 9.0;
      
      // Set particle positions
      double hadronPartPos = hadronBL*rndmPtr->gauss();
      double leptonPartPos = leptonBL*rndmPtr->gauss();

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

      vertexT = t_int;
      vertexX = x_int*c_c + z_int*s_c; // rotate to detector frame
      //vertexX = x_int*c_c - z_int*s_c; // rotate to hadron beam frame
      //vertexX = x_int;
      vertexY = y_int;
      vertexZ = x_int*(-1.0)*s_c + z_int*c_c; // rotate to detector frame
      //vertexZ = x_int*s_c + z_int*c_c; // rotate to hadron beam frame
      //vertexZ = z_int;


      /*
      // Calculate interaction point
      // Distance between lepton and hadron at z=0
      double dist = (hadronPartPos - leptonPartPos)/2.0;

      // Interaction Point
      double intPtH = hadronPartPos - dist;
      double intPtL = leptonPartPos + dist;
      //if(intPtH != intPtL) 
      if(TMath::Abs(intPtH - intPtL) > 0.0001)
	{
	  cout << "SOMETHING VERY WRONG IN VERTEX CALC" << endl;
	  cout << hadronPartPos << " " << leptonPartPos << endl;
	  cout << dist << " " << intPtH << " " << intPtL << endl;
	  cout << endl;
	} 

      // Interaction Time
      //double c = 299.792458; // mm/ns
      double c = 1.0; // c=1 in pythia coordinates
      double intT = (-1.0*intPtH)/c;

      // Set Vertex Z and t
      vertexZ = intPtH; // in mm
      vertexT = intT; // in mm

      // Once z-vertex has been establized, determine the collision point offset in x due to crabbing
      // Offset ~ theta_c/2k * sin(kz) ~ theta_c/2 * z
      xIPOffset = (mXAngle/2.0)*intPtH; // in mm
      */

      /*
      // Finally, find x and y vertex
      if(sigmaVertexX > 0.)
	{
	  vertexX += sigmaVertexX * rndmPtr->gauss();
	}
      vertexX += xIPOffset;

      if(sigmaVertexY > 0.)
	{
	  vertexY += sigmaVertexY * rndmPtr->gauss();
	}
      */
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
	  deltaPxA += pxLocal*TMath::Cos(mXAngle); // Compensate for Crossing Angle to Lab Frame

	  deltaPzA += pxLocal*TMath::Sin(mXAngle); // Projection of the x component of divergence onto the z-axis due to the crossing angle
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

  /*
  // Set beam vertex location by a two-dimensional Gaussian.
  if(allowVertexSpread) 
    {
      if(sigmaVertexZ > 0.) 
	{
	  //vertexZ = sigmaVertexZ * rndmPtr->gauss();
	}
    }
  */

  /*
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
  */

  /*
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

	  deltaPzA += pxLocal*TMath::Sin(mXAngle); // Projection of the x component of divergence onto the z-axis due to the crossing angle
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
  */

  /*
    


   */

  // _
  //(
  //-
}
