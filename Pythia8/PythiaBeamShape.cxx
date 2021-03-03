// Create pythia eP collisions and save to HepMC
// based on example/main36.cc

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
//#include "Pythia8Plugins/HepMC2.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "CustomPythia8ToHepMC3.h"
// #include "Pythia8ToHepMC3.h"


#include<iostream>
//#include<fstream>
//#include<cstdio>
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

// Test Event Record
void translate(Pythia8::Event&, Pythia8::Event&);


// A derived class to set beam momentum and interaction vertex spread. See main23.cc
class MyBeamShape : public BeamShape {

public:

  // Constructor.
  MyBeamShape() {}

  // Initialize beam parameters.
  // In this particular example we will reuse the existing settings names
  // but with modified meaning, so init() in the base class can be kept.
  //virtual void init( Settings& settings, Rndm* rndmPtrIn);

  // Set the two beam momentum deviations and the beam vertex.
  virtual void pick();

};


// Set the two beam momentum deviations and the beam vertex.
// Note that momenta are in units of GeV and vertices in mm,
// always with c = 1, so that e.g. time is in mm/c.

void MyBeamShape::pick() {

  // Reset all values.
  deltaPxA = deltaPyA = deltaPzA = deltaPxB = deltaPyB = deltaPzB
    = vertexX = vertexY = vertexZ = vertexT = 0.;

  // Set beam vertex location by a two-dimensional Gaussian.
  if(allowVertexSpread) 
    {
      // Set beam B longitudinal momentum as a triangular shape.
      // This corresponds to two step-function beams colliding.
      // Reuse sigmaVertexZ to represent maximum deviation in this case.
      if(sigmaVertexZ > 0.) 
	{
	  vertexZ     = sigmaVertexZ * ( 1. - sqrt(rndmPtr->flat()) );
	  if (rndmPtr->flat() < 0.5) vertexZ = -vertexZ;
	}
    }

  if(allowMomentumSpread)
    {
      double gaussX, gaussY;

      deltaPxA = -1.0*vertexZ;
      deltaPyA = 0.5*vertexZ;

      if(sigmaPxA > 0.)
	{
	  gaussX = rndmPtr->gauss();
	  deltaPxA += sigmaPxA * gaussX;
	}
      if(sigmaPyA > 0.)
	{
	  gaussY = rndmPtr->gauss();
	  deltaPyA += sigmaPyA * gaussY;
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
  TH1D *numCharged = new TH1D("numCharged","",50,0.,50.);

  TH1D *q2Hist = new TH1D("q2Hist","",100,-10.,0.);

  TH1D *yHist = new TH1D("yHist","",1000,-5.,0.);

  TH2D *phaseSpaceHist = new TH2D("phaseSpaceHist","",100,-4.,1.,100,0.,5.);

  TH2D *phaseSpaceJetEtaBin0 = new TH2D("phaseSpaceJetEtaBin0","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin1 = new TH2D("phaseSpaceJetEtaBin1","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin2 = new TH2D("phaseSpaceJetEtaBin2","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin3 = new TH2D("phaseSpaceJetEtaBin3","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin4 = new TH2D("phaseSpaceJetEtaBin4","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin5 = new TH2D("phaseSpaceJetEtaBin5","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin6 = new TH2D("phaseSpaceJetEtaBin6","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin7 = new TH2D("phaseSpaceJetEtaBin7","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin8 = new TH2D("phaseSpaceJetEtaBin8","",100,-4.,1.,100,0.,5.);
  TH2D *phaseSpaceJetEtaBin9 = new TH2D("phaseSpaceJetEtaBin9","",100,-4.,1.,100,0.,5.);

  // Beam Shape
  //TH1D *eCM = new TH1D("eCM","",100,-20.,20.);
  TH2D *pXY1 = new TH2D("pXY1","",100,-20.,20.,100,-20.,20.);
  TH2D *pXZProd1 = new TH2D("pXZProd1","",100,-20.,20.,100,-20.,20.);
  TH2D *pYZProd1 = new TH2D("pYZProd1","",100,-20.,20.,100,-20.,20.);
  TH2D *pXY2 = new TH2D("pXY2","",100,-20.,20.,100,-20.,20.);
  TH2D *pXZProd2 = new TH2D("pXZProd2","",100,-20.,20.,100,-20.,20.);
  TH2D *pYZProd2 = new TH2D("pYZProd2","",100,-20.,20.,100,-20.,20.);
  TH1D *pZ = new TH1D("pZ","",100,-20.,20.);

  TH1D *vtxX = new TH1D("vtxX","",100,-1.0,1.0);
  TH1D *vtxY = new TH1D("vtxY","",100,-1.0,1.0);
  TH1D *vtxZ = new TH1D("vtxZ","",100,-20.0,20.0);
  TH1D *vtxT = new TH1D("vtxT","",100,-100.0,100.0);

  TH2D *jetEtaVsP = new TH2D("jetEtaVsP","",200,0.,200.,100,-5.,5.);

  TH2D *jetQ2VsP = new TH2D("jetQ2VsP","",200,0.,200.,100,0.,5.);

  TH1D *jetEtaHist = new TH1D("jetEta","",100,-5.,5.);
  TH1D *jetPhiHist = new TH1D("jetPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH1D *partPtHist = new TH1D("partPt","",500,0.,50.);
  TH1D *partEtaHist = new TH1D("partEta","",100,-5.,5.);
  TH1D *partPhiHist = new TH1D("partPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsQ2Hist = new TH2D("partPhiVsQ2","",25,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH1D *jetElecPhiHist = new TH1D("jetElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *jetKtElecPhiHist = new TH1D("jetKtElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *jetKtElecPhiAbsHist = new TH1D("jetKtElecPhiAbs","",100,-1.0*TMath::Pi(),TMath::Pi());

  //TH2D *eVspzHist = new TH2D("eVspzHist","",5000,-250.,250.,2500,0.,250.);

  //TH2D *sHatVspTHat = new TH2D("sHatVspTHat","",1000,0.,100.,1000,0.,100.);

  //TH1D *resolvedEta = new TH1D("resolvedEta","",50,-5.,5.);
  //TH1D *directEta = new TH1D("directEta","",50,-5.,5.);

  Pythia8::Pythia p8;
  Pythia8::Event &event = p8.event;

  // A class to generate beam parameters according to own parametrization.
  BeamShapePtr myBeamShape = make_shared<MyBeamShape>();

  // Hand pointer to Pythia.
  // If you comment this out you get internal Gaussian-style implementation.
  p8.setBeamShapePtr(myBeamShape);


  //string basename = "nc";
  string basename = "dis";

  //p8.readFile("./photoproduction.steer");
  p8.readFile(argv[1]);
  
  p8.readString("Main:timesAllowErrors = 10000"); // allow more errors, eP is brittle


  // PDF Selection
  //p8.readString("PDF:pset = 2");
  //p8.settings.parm("PDF:pset", 2);


  // Uncomment to allow charged current.
  // p8.readString("WeakBosonExchange:ff2ff(t:W) = on");
  if ( basename=="nc" ){
    // Neutral current (with gamma/Z interference).
    //p8.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  }
  if ( basename=="dis" ){
    //p8.readString("PDF:lepton2gamma = on"); // DIS
    // // requires additional libraries
    // PDFPtr photonFlux = make_shared<Lepton2gamma2>(-11);
    // p8.setPhotonFluxPtr(photonFlux, 0);
  }
  // Phase-space cut: minimal Q2 of process.
  // changed configuration style for for 8.3?
  //p8.settings.parm("PhaseSpace:Q2Min", 25.0);
  // p8.readString("PhaseSpace:Q2Min = 25");


  // Set dipole recoil on. Necessary for DIS + shower.
  //p8.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit,
  // since rate known to match well to matrix elements everywhere.
  //p8.readString("SpaceShower:pTmaxMatch = 2");


  // QED radiation off lepton not handled yet by the new procedure.
  //p8.readString("PDF:lepton = off");
  //p8.readString("TimeShower:QEDshowerByL = off");
    
  // More Process details
  // ISR, hadronization, etc.
  /*
  p8.readString("HadronLevel:Decay = off");
  p8.readString("HadronLevel:all = on");
  p8.readString("PartonLevel:ISR = off");
  p8.readString("PartonLevel:MPI = off");
  p8.readString("PartonLevel:FSR = off");
  p8.readString("PromptPhoton:all=on");
  
  // Display
  p8.readString("Init:showProcesses = on");
  p8.readString("Init:showChangedSettings = on");
  p8.readString("Init:showMultipartonInteractions = on");
  p8.readString("Init:showChangedParticleData = on");

  p8.readString("Next:numberShowInfo = 1"); 
  p8.readString("Next:numberShowProcess = 1"); 
  p8.readString("Next:numberShowEvent = 5"); 

  // Random seed
  p8.readString("Random:setSeed = on");
  p8.readString("Random:seed = 42");
  */
  
  // initialize
  p8.init();

  Pythia8::Event hard;
  hard.init("(Pythia 6 conventions)", &p8.particleData);

//   // HepMC 2 
//   HepMC::Pythia8ToHepMC ToHepMC2;
//   HepMC::IO_GenEvent ascii_io( basename + "_hepmc2.dat", std::ios::out);

    // HepMC 3
  HepMC3::Pythia8ToHepMC3 ToHepMC3;

  std::shared_ptr<HepMC3::GenRunInfo> run = std::make_shared<HepMC3::GenRunInfo>();
  struct HepMC3::GenRunInfo::ToolInfo generator={
    std::string("Pythia8"),
    std::to_string(PYTHIA_VERSION).substr(0,5),
    std::string("Used generator")
  };
  run->tools().push_back(generator);
  // Can be used to save the name of the control card, if
  // pythia.readFile( ...) is used
  struct HepMC3::GenRunInfo::ToolInfo config={std::string(argv[1]),"1.0",std::string("Control cards")};
  run->tools().push_back(config);
  std::vector<std::string> names;
  for (int iWeight=0; iWeight < p8.info.nWeights(); ++iWeight) {
    std::string s=p8.info.weightLabel(iWeight);
    if (!s.length()) s=std::to_string((long long int)iWeight);
    names.push_back(s);
  }
  if (!names.size()) names.push_back("default");
  run->set_weight_names(names);
  // ############# Uncomment to write hepmc file
  //HepMC3::WriterAscii file(basename + "_hepmc3.dat",run);

  // Run
  //int nevents = 100;
  int numCharm = 0;
  int nevents = p8.mode("Main:numberOfEvents");
  for(int ev=0; ev<nevents; ev++)
    {
      if(!p8.next()) continue;
      // p8.event.list();

      // Beam Shape
      //eCM->Fill(eCMnow - eCMnom);
      pXY1->Fill(p8.process[1].px(),p8.process[1].py());
      pXZProd1->Fill(p8.process[1].zProd(),p8.process[1].px());
      pYZProd1->Fill(p8.process[1].zProd(),p8.process[1].py());
      pXY2->Fill(p8.process[2].px(),p8.process[2].py());
      pXZProd2->Fill(p8.process[2].zProd(),p8.process[2].px());
      pYZProd2->Fill(p8.process[2].zProd(),p8.process[2].py());
      pZ->Fill(p8.process[0].pz());
      
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
      
      int nCharged = 0;
      int indexeOut = -1;
      Vec4 totFinal;
      for(int i=0; i < p8.event.size(); i++)
	{
	  //cout << setprecision(15) << i << " " << p8.event[i].id() << " " << p8.event[i].status() << " " << p8.event[i].p().px() << " " << p8.event[i].p().py() << " " << p8.event[i].p().pz() << " " << p8.event[i].p().e() << " " << p8.event[i].eta() << endl;

	  if(p8.event[i].id() == 11 && p8.event[i].mother1() == 2 && p8.event[i].mother2() == 0) indexeOut = i;

	  if(p8.event[i].isFinal() && p8.event[i].isCharged()) 
	    {
	      nCharged++;
	    }

	  if(p8.event[i].isFinal() && p8.event[i].eta() > -3.5 && p8.event[i].eta() < 3.5 && i > 7)
	    {
	      //totFinal += p8.event[i].p();
	      partPtHist->Fill(p8.event[i].pT());
	      partPhiHist->Fill(p8.event[i].phi());
	      partEtaHist->Fill(p8.event[i].eta());

	      partPhiVsQ2Hist->Fill(std::log10(Q2),p8.event[i].phi());
	    }

	  if(p8.event[i].status() == -23)
	    {
	      totFinal += p8.event[i].p();
	      //if(p8.info.code() < 200) resolvedEta->Fill(p8.event[i].eta());
	      //if(p8.info.code() > 200) directEta->Fill(p8.event[i].eta());
	    }
	}
      //cout << endl;
      
      /*
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
      */

      if(TMath::Abs(event[5].id()) == 4)
	{
	  numCharm++;
	  //cout << p8.info.name() << " " << p8.info.code() << endl;
	  for(int i=0; i < p8.event.size(); i++)
	    {
	      //cout << i << " " << p8.event[i].id() << " " << p8.event[i].status() << " " << p8.event[i].mother1() << " " << p8.event[i].mother2() << endl;
	    }
	}

      //cout << setprecision(15) << "TEST Q2 x y " << Q2 << " " << x << " " << y << endl;
      //cout << ev << " " << totFinal.e() << " " << totFinal.pz() << " " << totFinal.pT() << endl;
      //cout << "TEST Q2 Q2Fac Q2Ren pTHat sHat rapidity " << Q2 << " " << p8.info.Q2Fac() << " " << p8.info.Q2Ren() << " " << p8.info.pTHat() << " " << p8.info.sHat() << " " << p8.info.y() << endl;
      //cout << "TEST id1 id2 x1 x2 " << p8.info.id1() << " " << p8.info.id2() << " " << p8.info.x1() << " " << p8.info.x2() << endl;
      //cout << setprecision(10) << "TEST y yJB " << y << " " << (totFinal.e() - totFinal.pz())/(2.0*peIn.e()) << endl;

      q2Hist->Fill(std::log10(Q2));
      yHist->Fill(std::log10(y));
      phaseSpaceHist->Fill(std::log10(x),std::log10(Q2));

      if(y < 0.95 && y > 0.01 && jetP > 10.0)
	{
	  if(jetEta < -1.0) phaseSpaceJetEtaBin0->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > -1.0 && jetEta < -0.5) phaseSpaceJetEtaBin1->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > -0.5 && jetEta < 0.0) phaseSpaceJetEtaBin2->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 0.0 && jetEta < 0.5) phaseSpaceJetEtaBin3->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 0.5 && jetEta < 1.0) phaseSpaceJetEtaBin4->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 1.0 && jetEta < 1.5) phaseSpaceJetEtaBin5->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 1.5 && jetEta < 2.0) phaseSpaceJetEtaBin6->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 2.0 && jetEta < 2.5) phaseSpaceJetEtaBin7->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 2.5 && jetEta < 3.0) phaseSpaceJetEtaBin8->Fill(std::log10(x),std::log10(Q2));
	  if(jetEta > 3.0 && jetEta < 3.5) phaseSpaceJetEtaBin9->Fill(std::log10(x),std::log10(Q2));

	  jetEtaVsP->Fill(jetP,jetEta);

	  jetQ2VsP->Fill(jetP,std::log10(Q2));

	  jetEtaHist->Fill(jetEta);
	  jetPhiHist->Fill(jetPhi);

	  jetElecPhiHist->Fill(jetPhi - elecPhi - TMath::Pi());
	  jetKtElecPhiHist->Fill(jetKtPhi - elecPhi - TMath::Pi());
	  jetKtElecPhiAbsHist->Fill(TMath::Abs(jetKtPhi - elecPhi - TMath::Pi()));
	}

      //eVspzHist->Fill(totFinal.pz(),totFinal.e());
      //sHatVspTHat->Fill(p8.info.pTHat(),std::sqrt(p8.info.sHat()));

      /*
      if(y>0.2 && y<0.8 && Q2 > 0.00001)
	{
	  for(int i=0; i < p8.event.size(); i++)
	    {
	      if(p8.event[i].status() == -23)
		{
		  if(p8.info.code() < 200 && p8.event[i].pT() > 10.0) resolvedEta->Fill(p8.event[i].eta());
		  if(p8.info.code() > 200 && p8.event[i].pT() > 10.0) directEta->Fill(p8.event[i].eta());
		}
	    }
	}
      */
    
      /*
      numCharged->Fill(nCharged);
      cout << "Event " << ev << " was type " << p8.info.name() << " (" << p8.info.code() << ")" << endl;
      cout << "Event " << ev << " photon mode " << p8.info.photonMode() << endl;
      cout << "Event " << ev << " is resolved " << p8.info.isResolved() << endl;
      cout << "Event " << ev << " Q2 : x : y : w2 " << Q2 << " : " << x << " : " << y << " : " << W2 << endl;
      cout << "Event " << ev << " has " << nCharged << " charged particles" << endl;
      cout << endl;
      */


      //     // Construct new empty HepMC2 event and fill it.
      //     // Units will be as chosen for HepMC build; but can be changed
      //     // by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
      //     HepMC::GenEvent* hepmc2evt = new HepMC::GenEvent();
      //     ToHepMC2.fill_next_event( p8, hepmc2evt );
      //     // write it out
      //     ascii_io << hepmc2evt;
      //     delete hepmc2evt;    
      
      // Construct new empty HepMC3 event and fill it.
      HepMC3::GenEvent hepmc3evt( HepMC3::Units::GEV, HepMC3::Units::MM );
      ToHepMC3.fill_next_event(p8.event, &hepmc3evt, -1, &p8.info);
      // write it out
      // ################# Uncomment to Write hepmc file
      //file.write_event(hepmc3evt);
      
      if(ev==0) 
	{
	  std::cout << "First event: " << std::endl;
	  HepMC3::Print::listing(hepmc3evt);
	} 
    }

  // Test Translated record
  //translate(event,hard);
  //hard.list();

  // #################### hepmc
  //file.close();
  
  p8.stat();

  cout << "Num Charm = " << numCharm << endl;
  
  // Write and Close Root File
  ofile->Write();
  ofile->Close();
  // Done.
  return 0;

					  
}

// Routine to identify the hard process of the event after
// accounting for the affects of ISR and resonance decay.

void translate(Pythia8::Event& event, Pythia8::Event& hard) {

    // Reset translated event record.
    hard.reset();

    // Identify primordial partons with kT smearing.
    int iL1 = event[1].daughter1();
    int iL2 = event[2].daughter1();
    int iLast = min(iL1,iL2);

    // Identify initial hard partons.
    int iP1 = 3;
    int iP2 = 4;
    int iEv=-1;

    // Counters to identify ISR radiator, recoil, and radiation.
    int i41 = 0, i42 = 0, i43 = 0;
    while ( iEv < iLast ) {
      iEv++;

      // Identify hard process by status in the 20s.
      int iStatus = event[iEv].status();
      if ( abs(iStatus) > 20 && abs(iStatus) < 30 ) {
        hard.append(event[iEv]);
        if ( iStatus==-21 ) {
          hard[hard.size()-1].mothers(0,0);
        } else {
          hard[hard.size()-1].mothers(1,2);
        }
        continue;
      }

      // Follow the flow of the hard event.
      int iDa1 = event[iEv].daughter2();
      // Test to determine origin of this shower.
      bool hardTest = (iDa1 == iP1) || (iDa1 == iP2);
      // Status -41 is the ISR radiator
      // Status -42 is the ISR recoil
      // Status -61 is the primordial parton with kT smearing
      // Status -53 is a recoil in initial-final dipole radiation
      if ( iStatus == -41 && hardTest ) {
        i41 = iEv;
      } else if ( iStatus == -42 && hardTest ) {
        i42 = iEv;
      } else if ( iStatus == -61 && hardTest ) {
        if ( i41 == 0 ) {
          i41 = iEv;
        } else {
          i42 = iEv;
        }
      } else if ( iStatus == -53 ) {
        if ( iDa1 == iP1 ) {
          iP1 = iEv;
        } else if ( iDa1 == iP2 ) {
          iP2 = iEv;
        }
      }

      if ( !(i41 > 0 && i42 > 0) ) continue;
      int ik2 = event[i41].daughter2();
      int ik1 = event[i42].daughter2();

      if ( event[ik2].pz() > 0 ) {
        int iTemp=ik1;
        ik1=ik2;
        ik2=iTemp;
      }

      // Boost to CM frame of hard partons.
      RotBstMatrix toCMS;
      toCMS.toCMframe(event[ik1].p(),event[ik2].p());

      // Momentum of off-shell incoming parton.
      i43 = event[i41].daughter1();
      Vec4 pBoost = event[i41].p()-event[i43].p();
      RotBstMatrix toHard;
      if ( event[i41].pz() > 0 ) {
        toHard.fromCMframe(pBoost,event[i42].p());
      } else {
        toHard.fromCMframe(event[i42].p(),pBoost);
      }

      // Boost to CM frame of old initiators,
      // then boost from frame of new initiators.
      for( int i = 0; i< hard.size(); ++i ) {
        hard[i].rotbst(toCMS);
        hard[i].rotbst(toHard);
      }

      // Update counter to location of new parton
      iEv = i43;
      // Update event history
      iDa1 = event[i41].daughter2();
      if ( iDa1 == iP1 ) {
        iP1 = i41;
        iP2 = i42;
      } else {
        iP2 = i41;
        iP1 = i42;
      }
      i41=0; i42=0;

    }

    // Handle kT smearing of initial partons.
    RotBstMatrix ref0, ref1;
    if ( event[iP1].pz() > 0 ) {
      ref0.toCMframe(event[iP1].p(),event[iP2].p());
    } else {
      ref0.toCMframe(event[iP2].p(),event[iP1].p());
    }
    if ( event[iL1].pz() > 0 ) {
      ref1.fromCMframe(event[iL1].p(),event[iL2].p());
    } else {
      ref1.fromCMframe(event[iL2].p(),event[iL1].p());
    }
    for( int i=0; i< hard.size(); ++i ) {
      hard[i].rotbst(ref0);
      hard[i].rotbst(ref1);
    }

    hard[1].daughters(3,hard.size()-1);
    hard[2].daughters(3,hard.size()-1);
    int iMax = hard.size();

    // Add resonance decays; start here.
    for( int i = 3; i < iMax && i < event.size(); ++i ) {
      int ilast = -1;
      int ida1, ida2 = -1;
      if ( hard[i].status()!=-22 ) {
        hard[i].statusPos();
        continue;
      }
      if ( hard[i].mother1()==3 ) {
        ilast = hard[i].daughter1();
        ida1 = event[ilast].daughter1();
        ida2 = event[ilast].daughter2();
      } else {
        ilast = i;
        ida1 = hard[ilast].daughter1();
        ida2 = hard[ilast].daughter2();
      }

      // Resonance decays occur when there are multiple daughters and
      // it is NOT FSR.
      while( ilast > 0 && ilast < event.size() ) {
        if ( ida1 != ida2 && event[ida1].status() != -51 ) break;
        ilast = ida1;
        ida1 = event[ilast].daughter1();
        ida2 = event[ilast].daughter2();
      }

      // Add daughters to the event record, boosting to the frame
      // of the mother.
      if ( ilast > 0 ) {
        Vec4 pall = Vec4(0,0,0,0);
        for(int ida = ida1; ida <= ida2; ++ida) {
          pall += event[ida].p();
        }
        RotBstMatrix toResonance;
        toResonance.bst( pall, hard[i].p() );
        int nDau = 0;
        for(int ida = ida1; ida <= ida2; ++ida) {
          Particle tmp = Particle(event[ida]);
          tmp.rotbst(toResonance);
          hard.append( tmp );
          hard[hard.size()-1].mothers(i,i);
          iMax++;
          nDau++;
        }
        int ip1 = hard.size() - nDau;
        hard[i].daughters(ip1,ip1+nDau-1);
      }
    }
}
