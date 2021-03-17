// Create pythia eP collisions and save to HepMC
// based on example/main36.cc

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"


#include<iostream>
using std::cout;
using std::cerr;
using std::endl;

#include<string>
using std::string;

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"

using namespace Pythia8;

// A derived class to set beam momentum and interaction vertex spread. See main23.cc
class MyBeamShape : public BeamShape {

public:

  // Constructor.
  MyBeamShape(double ion, double lepton, double xAngle, double hadCrab, double lepCrab) 
  {
    mIonBeamEnergy = ion;
    mLeptonBeamEnergy = lepton;
    mXAngle = xAngle;
    mHadronCrabSize = hadCrab;
    mLeptonCrabSize = lepCrab;
  }

  // Initialize beam parameters.
  // In this particular example we will reuse the existing settings names
  // but with modified meaning, so init() in the base class can be kept.
  //virtual void init( Settings& settings, Rndm* rndmPtrIn);

  // Set the two beam momentum deviations and the beam vertex.
  virtual void pick();

protected:

  double mIonBeamEnergy;
  double mLeptonBeamEnergy;
  double mXAngle;
  double mHadronCrabSize;
  double mLeptonCrabSize;

};


// Set the two beam momentum deviations and the beam vertex.
// Note that momenta are in units of GeV and vertices in mm,
// always with c = 1, so that e.g. time is in mm/c.

void MyBeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Change in Beam Pz
  double tmpPzA, tmpPzB;
  tmpPzA = tmpPzB = 0.;

  //cout << mMyTest << endl;

  // Set beam vertex location by a two-dimensional Gaussian.
  if(allowVertexSpread) 
    {
      // Set beam B longitudinal momentum as a triangular shape.
      // This corresponds to two step-function beams colliding.
      // Reuse sigmaVertexZ to represent maximum deviation in this case.
      if(sigmaVertexZ > 0.) 
	{
	  vertexZ = sigmaVertexZ * rndmPtr->gauss();
	  //vertexZ     = sigmaVertexZ * ( 1. - sqrt(rndmPtr->flat()) );
	  //if (rndmPtr->flat() < 0.5) vertexZ = -vertexZ;
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

      //cout << tmpPzA << endl;

      // Modify Px Due to Crossing Angle
      deltaPxA += (mIonBeamEnergy + tmpPzA)*TMath::Sin(mXAngle);

      // Modify Pz Due to Crossing Angle
      deltaPzA += (mIonBeamEnergy + tmpPzA)*TMath::Cos(mXAngle) - mIonBeamEnergy;
      //deltaPzA += tmpPzA;

      // Electron Beam Has 0 Crossing Angle
      deltaPzB += tmpPzB;

      //cout << setprecision(10) << mIonBeamEnergy << " " << TMath::Cos(mXAngle) << " " << deltaPzA << endl;
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
	
  /*
  // Set beam A transverse momentum deviation by a two-dimensional Gaussian.
  if (allowMomentumSpread) {
    //cout << sigmaPxA << " " << sigmaPyA << " " << sigmaPzA << endl;
    //cout << maxDevA << endl;

    cout << "Enter pick()" << endl;

    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaPxA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxA  = sigmaPxA * gauss;
        totalDev += gauss * gauss;
	cout << "X gauss = " << gauss << endl;
	cout << "X totalDev = " << totalDev << endl;
      }
      if (sigmaPyA > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyA  = sigmaPyA * gauss;
        totalDev += gauss * gauss;
	cout << "Y gauss = " << gauss << endl;
	cout << "Y totalDev = " << totalDev << endl;
      }
    } while (totalDev > maxDevA * maxDevA);

    cout << endl;

    //cout << deltaPxA << " " << deltaPyA << endl;

    // Set beam A longitudinal momentum as a triangular shape.
    // Reuse sigmaPzA to represent maximum deviation in this case.
    if (sigmaPzA > 0.) {
      deltaPzA    = sigmaPzA * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) deltaPzA = -deltaPzA;
    }

    // Set beam B transverse momentum deviation by a two-dimensional Gaussian.
    do {
      totalDev = 0.;
      if (sigmaPxB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPxB  = sigmaPxB * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaPyB > 0.) {
        gauss     = rndmPtr->gauss();
        deltaPyB  = sigmaPyB * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevB * maxDevB);

    // Set beam B longitudinal momentum as a triangular shape.
    // Reuse sigmaPzB to represent maximum deviation in this case.
    if (sigmaPzB > 0.) {
      deltaPzB = sigmaPzB * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) deltaPzB = -deltaPzB;
    }
  }

  // Set beam vertex location by a two-dimensional Gaussian.
  if (allowVertexSpread) {
    double totalDev, gauss;
    do {
      totalDev = 0.;
      if (sigmaVertexX > 0.) {
        gauss     = rndmPtr->gauss();
        vertexX   = sigmaVertexX * gauss;
        totalDev += gauss * gauss;
      }
      if (sigmaVertexY > 0.) {
        gauss     = rndmPtr->gauss();
        vertexY   = sigmaVertexY * gauss;
        totalDev += gauss * gauss;
      }
    } while (totalDev > maxDevVertex * maxDevVertex);

    // Set beam B longitudinal momentum as a triangular shape.
    // This corresponds to two step-function beams colliding.
    // Reuse sigmaVertexZ to represent maximum deviation in this case.
    if (sigmaVertexZ > 0.) {
      vertexZ     = sigmaVertexZ * ( 1. - sqrt(rndmPtr->flat()) );
      if (rndmPtr->flat() < 0.5) vertexZ = -vertexZ;

      // Set beam collision time flat between +-(sigmaVertexZ - |vertexZ|).
      // This corresponds to two step-function beams colliding (with v = c).
      vertexT = (2. * rndmPtr->flat() - 1.) * (sigmaVertexZ - abs(vertexZ));
    }

    cout << offsetX << " " << offsetY << " " << offsetZ << " " << offsetT << endl;

    // Add offset to beam vertex.
    vertexX      += offsetX;
    vertexY      += offsetY;
    vertexZ      += offsetZ;
    vertexT      += offsetT;
  }
  */

}


int main(int argc, char* argv[])
{
  //const char* steerFile = 

  // Open Root File
  TFile *ofile = TFile::Open("test.hist.root","recreate");

  // Histos
  TH1D *q2Hist = new TH1D("q2Hist","",100,-10.,0.);
  TH1D *yHist = new TH1D("yHist","",1000,-5.,0.);
  TH2D *phaseSpaceHist = new TH2D("phaseSpaceHist","",100,-4.,1.,100,0.,5.);

  // Beam Shape
  TH1D *eCM = new TH1D("eCM","",10000,-0.05,0.05);
  TH2D *pXY1 = new TH2D("pXY1","",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd1 = new TH2D("pXZProd1","",100,-40.,40.,10000,-10.0,10.0);
  TH2D *pYZProd1 = new TH2D("pYZProd1","",100,-40.,40.,10000,-10.0,10.0);
  TH2D *pXY2 = new TH2D("pXY2","",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd2 = new TH2D("pXZProd2","",100,-40.,40.,10000,-10.0,10.0);
  TH2D *pYZProd2 = new TH2D("pYZProd2","",100,-40.,40.,10000,-10.0,10.0);
  TH1D *pZ1 = new TH1D("pZ1","",10000,270.0,280.0);
  TH1D *pZ2 = new TH1D("pZ2","",10000,-18.05,-17.95);

  TH1D *vtxX = new TH1D("vtxX","",100,-1.0,1.0);
  TH1D *vtxY = new TH1D("vtxY","",100,-1.0,1.0);
  TH1D *vtxZ = new TH1D("vtxZ","",100,-40.0,40.0);
  TH1D *vtxT = new TH1D("vtxT","",100,-100.0,100.0);

  // Particle Quantities
  TH1D *partPtHist = new TH1D("partPt","",500,0.,50.);
  TH1D *partEtaHist = new TH1D("partEta","",100,-5.,5.);
  TH1D *partPhiHist = new TH1D("partPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsEtaHist = new TH2D("partPhiVsEta","",100,-5.,5.,100,-1.0*TMath::Pi(),TMath::Pi());
  TH2D *partPhiVsQ2Hist = new TH2D("partPhiVsQ2","",25,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  // Hi Pt
  TH1D *partPtHiHist = new TH1D("partPtHi","",500,0.,50.);
  TH1D *partEtaHiHist = new TH1D("partEtaHi","",100,-5.,5.);
  TH1D *partPhiHiHist = new TH1D("partPhiHi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsEtaHiHist = new TH2D("partPhiVsEtaHi","",100,-5.,5.,100,-1.0*TMath::Pi(),TMath::Pi());
  TH2D *partPhiVsQ2HiHist = new TH2D("partPhiVsQ2Hi","",25,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  // "Jet" Quantities
  TH2D *jetEtaVsP = new TH2D("jetEtaVsP","",200,0.,200.,100,-5.,5.);

  TH2D *jetQ2VsP = new TH2D("jetQ2VsP","",200,0.,200.,100,0.,5.);

  TH1D *jetEtaHist = new TH1D("jetEta","",100,-5.,5.);
  TH1D *jetPhiHist = new TH1D("jetPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  // "Jet" - Electron Quantities
  TH1D *jetElecPhiHist = new TH1D("jetElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *jetKtElecPhiHist = new TH1D("jetKtElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *jetKtElecPhiAbsHist = new TH1D("jetKtElecPhiAbs","",100,-1.0*TMath::Pi(),TMath::Pi());


  // Set Up Pythia Event
  Pythia8::Pythia p8;
  Pythia8::Event &event = p8.event;

  // A class to generate beam parameters according to own parametrization.
  BeamShapePtr myBeamShape = make_shared<MyBeamShape>(275.0,18.0,0.025,0.0030,-0.0015);

  // Hand pointer to Pythia.
  // If you comment this out you get internal Gaussian-style implementation.
  p8.setBeamShapePtr(myBeamShape);

  // Read in Steering File & Make Other Settings
  p8.readFile(argv[1]);
  p8.readString("Main:timesAllowErrors = 10000"); // allow more errors, eP is brittle

  // Initialize Pythia
  p8.init();


  // Record Nominal COM Energy
  double eCMnom = p8.info.eCM();


  // Run
  int nevents = p8.mode("Main:numberOfEvents");
  for(int ev=0; ev<nevents; ev++)
    {
      if(!p8.next()) continue;
      // p8.event.list();

      // Modified COM Energy
      double eCMnow = p8.info.eCM();

      //cout << setprecision(10) << ev << " " << eCMnom << " " << eCMnow << endl;
      //cout << setprecision(10) << ev << " " << p8.process[1].px() << " " << p8.process[1].py() << " " << p8.process[1].pz() << " " << p8.process[1].e() << endl;
      //cout << setprecision(10) << ev << " " << p8.process[2].px() << " " << p8.process[2].py() << " " << p8.process[2].pz() << " " << p8.process[2].e() << endl;
      //cout << endl;

      // Beam Shape
      eCM->Fill(eCMnow - eCMnom);
      pXY1->Fill(p8.process[1].px(),p8.process[1].py());
      pXZProd1->Fill(p8.process[1].zProd(),p8.process[1].px());
      pYZProd1->Fill(p8.process[1].zProd(),p8.process[1].py());
      pXY2->Fill(p8.process[2].px(),p8.process[2].py());
      pXZProd2->Fill(p8.process[2].zProd(),p8.process[2].px());
      pYZProd2->Fill(p8.process[2].zProd(),p8.process[2].py());
      pZ1->Fill(p8.process[1].pz());
      pZ2->Fill(p8.process[2].pz());
      
      vtxX->Fill(p8.process[0].xProd());
      vtxY->Fill(p8.process[0].yProd());
      vtxZ->Fill(p8.process[0].zProd());
      vtxT->Fill(p8.process[0].tProd());

      // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
      Pythia8::Vec4 pProton = event[1].p();
      Pythia8::Vec4 peIn    = event[4].p();
      Pythia8::Vec4 peOut   = event[6].p();
      Pythia8::Vec4 pPhoton = peIn - peOut;

      // Q2, W2, Bjorken x, y.
      double Q2    = - pPhoton.m2Calc();
      double W2    = (pProton + pPhoton).m2Calc();
      double x     = Q2 / (2. * pProton * pPhoton);
      double y     = (pProton * pPhoton) / (pProton * peIn);
      double jetPt = event[5].pT();
      double jetEta = event[5].eta();
      double jetPhi = event[5].phi();
      double jetP = event[5].pAbs();
      double elecPhi = event[6].phi();
      double jetKtPhi = event[8].phi();

      // Event Level
      q2Hist->Fill(std::log10(Q2));
      yHist->Fill(std::log10(y));
      phaseSpaceHist->Fill(std::log10(x),std::log10(Q2));
      
      // Loop Through Particles
      for(int i=0; i < p8.event.size(); i++)
	{
	  bool partFin = p8.event[i].isFinal();
	  double partPt = p8.event[i].pT();
	  double partEta = p8.event[i].eta();
	  double partPhi = p8.event[i].phi();

	  if(partFin && partEta>-3.5 && partEta<3.5 && y<0.95 && y>0.01 && i > 7)
	    {
	      partPtHist->Fill(partPt);
	      partPhiHist->Fill(partPhi);
	      partEtaHist->Fill(partEta);

	      partPhiVsEtaHist->Fill(partEta,partPhi);
	      partPhiVsQ2Hist->Fill(std::log10(Q2),partPhi);

	      if(partPt > 1.0)
		{
		  partPtHiHist->Fill(partPt);
		  partPhiHiHist->Fill(partPhi);
		  partEtaHiHist->Fill(partEta);
		  
		  partPhiVsEtaHiHist->Fill(partEta,partPhi);
		  partPhiVsQ2HiHist->Fill(std::log10(Q2),partPhi);
		}
	    }
	}
      //cout << endl;

      if(y < 0.95 && y > 0.01 && jetP > 10.0)
	{
	  jetEtaVsP->Fill(jetP,jetEta);

	  jetQ2VsP->Fill(jetP,std::log10(Q2));

	  jetEtaHist->Fill(jetEta);
	  jetPhiHist->Fill(jetPhi);

	  jetElecPhiHist->Fill(jetPhi - elecPhi - TMath::Pi());
	  jetKtElecPhiHist->Fill(jetKtPhi - elecPhi - TMath::Pi());
	  jetKtElecPhiAbsHist->Fill(TMath::Abs(jetKtPhi - elecPhi - TMath::Pi()));
	}
    }

  // List Statistics  
  p8.stat();
  
  // Write and Close Root File
  ofile->Write();
  ofile->Close();
  // Done.
  return 0;

					  
}
