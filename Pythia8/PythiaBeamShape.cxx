// Generate Pythia events and apply EIC Beam Effects

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "eicBeamShape.h"

#include "fastjet/ClusterSequence.hh"

#include<iostream>
using std::cout;
using std::cerr;
using std::endl;

#include<string>
using std::string;

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TVector2.h"

using namespace Pythia8;
using namespace fastjet;


int main(int argc, char* argv[])
{
  //const char* steerFile = 
  
  if(argc != 9)
    {
      cerr << "Wrong number of arguments" << endl;
      cerr << "program.exe steer configuration Allow Crab hadronE leptonE xangle out.hist.root out.hepmc" << endl;
      exit(EXIT_FAILURE);
    }

  const int config     = atoi(argv[2]);
  const int crab       = atoi(argv[3]);
  const double hadE    = atof(argv[4]);
  const double lepE    = atof(argv[5]);
  const double xing    = atof(argv[6]);
  const char* rootOut  = argv[7];
  const char* hepmcOut = argv[8];

  cout << "Steering File = " << argv[1] << endl;
  if(config == 0) cout << "Checking of Steering File Settings Suspended - Components Zeroed out to Turn off Specific effects" << endl;
  if(config == 1) cout << "Configuration = High Divergence" << endl;
  if(config == 2) cout << "Configuration = High Acceptance" << endl; 
  if(crab == 0) cout << "Disabling Momentum Kick from Crab Crossing" << endl;
  else cout << "Including Momentum Kick from Crab Crossing" << endl;
  cout << "Hadron Energy = " << hadE << endl;
  cout << "Lepton Energy = " << lepE << endl;
  cout << "Beam Crossing Angle (Radians) = " << xing << endl;
  cout << "Root Output = " << rootOut << endl;
  cout << "HepMC Output = " << hepmcOut << endl;

  // Open Root File
  TFile *ofile = TFile::Open(rootOut,"recreate");

  // Histos
  TH1D *q2Hist = new TH1D("q2Hist","Log Q2",100,-10.,0.);
  TH1D *yHist = new TH1D("yHist","Inelasticity",1000,-5.,0.);
  TH2D *phaseSpaceHist = new TH2D("phaseSpaceHist","Log Q2 Vs Log x",100,-4.,1.,100,0.,5.);

  // Beam Shape
  TH1D *eCM = new TH1D("eCM","Modified - Nominal CM Energy",10000,-1.0,1.0);
  TH2D *pXY1 = new TH2D("pXY1","Hadron Beam Py Vs Px",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd1 = new TH2D("pXZProd1","Hadron Beam Px Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pYZProd1 = new TH2D("pYZProd1","Hadron Beam Py Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pXY2 = new TH2D("pXY2","Lepton Beam Py Vs Px",10000,-10.0,10.0,10000,-10.0,10.0);
  TH2D *pXZProd2 = new TH2D("pXZProd2","Lepton Beam Px Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH2D *pYZProd2 = new TH2D("pYZProd2","Lepton Beam Py Vs Vertex z",5000,-500.,500.,10000,-10.0,10.0);
  TH1D *pZ1 = new TH1D("pZ1","Hadron Beam Pz",10000,0.0,280.0);
  TH1D *pZ2 = new TH1D("pZ2","Lepton Beam Pz",10000,-20.,0.);

  TH1D *atan2PxPz1Hist = new TH1D("atan2PxPz1","",2500,0.0,0.05);
  //TH1D *atan2PxPz1Hist = new TH1D("atan2PxPz1","",500,0.24,0.26);
  TH1D *atan2PyPz1Hist = new TH1D("atan2PyPz1","",2500,-0.01,0.01);
  TH1D *atan2PyPtot1Hist = new TH1D("atan2PyPtot1","",500,-0.001,0.001);

  TH1D *vtxX = new TH1D("vtxX","Vertex x;[mm]",5000,-5.0,5.0);
  TH1D *vtxY = new TH1D("vtxY","Vertex y;[mm]",5000,-5.0,5.0);
  TH1D *vtxZ = new TH1D("vtxZ","Vertex z;[mm]",5000,-500.0,500.0);
  TH1D *vtxT = new TH1D("vtxT","Time;[mm]",5000,-500.0,500.0);
  TH2D *vtxYvsX = new TH2D("vtxYvsX","Vertex Y vs X;X [mm];Y [mm]",5000,-5.0,5.0,5000,-5.0,5.0);
  TH2D *vtxXvsT = new TH2D("vtxXvsT","Vertex X vs T;T [mm];X [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
  TH2D *vtxXvsZ = new TH2D("vtxXvsZ","Vertex X vs Z;Z [mm];X [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
  TH2D *vtxYvsZ = new TH2D("vtxYvsZ","Vertex Y vs Z;Z [mm];Y [mm]",5000,-500.0,500.0,5000,-5.0,5.0);
  TH2D *vtxTvsZ = new TH2D("vtxTvsZ","Interaction Time Vs Z-vertex;Z [mm];T [mm]",5000,-500.0,500.0,5000,-500.0,500.0);
  TH2D *vtxXvsTZSum = new TH2D("vtxXvsTZSum","Vertex X vs T+Z;T+Z [mm];X [mm]",5000,-500.,500.,5000,-5.,5.);
  TH2D *vtxXvsTZDiff = new TH2D("vtxXvsTZDiff","Vertex X vs T-Z;T-Z [mm];X [mm]",5000,-500.,500.,5000,-5.,5.);
  
  TH2D *lepVsHadPartZ = new TH2D("lepVsHadPartZ","Intrabunch Z Positions of Colliding Leptons Vs Hadrons",5000,-500.0,500.0,5000,-500.0,500.0);

  // Particle Quantities
  TH1D *partPtHist = new TH1D("partPt","Final State Particle Pt",500,0.,50.);
  TH1D *partEtaHist = new TH1D("partEta","Final State Particle Eta",400,-10.,10.);
  TH1D *partPhiHist = new TH1D("partPhi","Final State Particle Phi",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPtVsEtaHist = new TH2D("partPtVsEta","Final State Particle Pt Vs Eta",400,-10.,10.,500,0.,50.);
  TH2D *partPhiVsEtaHist = new TH2D("partPhiVsEta","Final State Particle Phi Vs Eta",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partStatusVsEtaHist = new TH2D("partStatusVsEta","Particle Status Code Vs Eta",400,-10.,10.,200,0.,200.);

  TH2D *partPtVsEtaBRHist = new TH2D("partPtVsEtaBR","Final State Particle Pt Vs Eta: Beam Remnant",400,-10.,10.,500,0.,50.);
  TH2D *partPhiVsEtaBRHist = new TH2D("partPhiVsEtaBR","Final State Particle Phi Vs Eta: Beam Remnant",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPtVsEtaHadHist = new TH2D("partPtVsEtaHad","Final State Particle Pt Vs Eta: Hadronization",400,-10.,10.,500,0.,50.);
  TH2D *partPhiVsEtaHadHist = new TH2D("partPhiVsEtaHad","Final State Particle Phi Vs Eta: Hadronization",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPtVsEtaDecayHist = new TH2D("partPtVsEtaDecay","Final State Particle Pt Vs Eta: Decay",400,-10.,10.,500,0.,50.);
  TH2D *partPhiVsEtaDecayHist = new TH2D("partPhiVsEtaDecay","Final State Particle Phi Vs Eta: Decay",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  // Hi Pt
  TH1D *partPtHiHist = new TH1D("partPtHi","Final State Particle Pt (Pt > 1 GeV)",500,0.,50.);
  TH1D *partEtaHiHist = new TH1D("partEtaHi","Final State Particle Eta (Pt > 1 GeV)",400,-10.,10.);
  TH1D *partPhiHiHist = new TH1D("partPhiHi","Final State Particle Phi (Pt > 1 GeV)",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsEtaHiHist = new TH2D("partPhiVsEtaHi","Final State Particle Phi Vs Eta (Pt > 1 GeV)",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partStatusVsEtaHiHist = new TH2D("partStatusVsEtaHi","Particle Status Code Vs Eta (Pt > 1 GeV)",400,-10.,10.,200,0.,200.);

  // 2->1 Outgoing Quark ("Jet") Quantities
  TH2D *pJetEtaVsP = new TH2D("pJetEtaVsP","Parton Eta Vs Momentum",200,0.,200.,100,-5.,5.);

  TH2D *pJetQ2VsP = new TH2D("pJetQ2VsP","Q2 Vs Parton Momentum",200,0.,200.,100,0.,5.);

  TH1D *pJetPtHist = new TH1D("pJetPt","Parton Pt",500,0.,50.);
  TH1D *pJetEtaHist = new TH1D("pJetEta","Parton Eta",100,-5.,5.);
  TH1D *pJetPhiHist = new TH1D("pJetPhi","Parton Phi",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *pJetPtVsEtaHist = new TH2D("pJetPtVsEta","Parton Pt Vs Eta",400,-10.,10.,500,0.,50.);
  TH2D *pJetPhiVsEtaHist = new TH2D("pJetPhiVsEta","Parton Phi Vs Eta",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());

  // "Jet" - Electron Quantities
  TH1D *pJetElecPhiHist = new TH1D("pJetElecPhi","Hard Parton  - Electron Phi - Pi",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *pJetKtElecPhiHist = new TH1D("pJetKtElecPhi","Parton - Electron Phi - Pi",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *pJetKtElecPhiAbsHist = new TH1D("pJetKtElecPhiAbs","|Parton - Electron Phi - Pi|",100,-1.0*TMath::Pi(),TMath::Pi());

  // True Jet Quantities
  TH1D *jetPtHist = new TH1D("jetPt","Jet Pt",500,0.,50.);
  TH1D *jetEtaHist = new TH1D("jetEta","Jet Eta",100,-5.,5.);
  TH1D *jetPhiHist = new TH1D("jetPhi","Jet Phi",100,-1.0*TMath::Pi(),TMath::Pi());

TH2D *jetPtVsPtNoCutHist = new TH2D("jetPtVsPtNoCut","Jet Pt Vs Parton Pt",500,0.,50.,500,0.,50.);

  TH2D *jetPtVsEtaHist = new TH2D("jetPtVsEta","Jet Pt Vs Eta",100,-5.,5.,500,0.,50.);
  TH2D *jetPhiVsEtaHist = new TH2D("jetPhiVsEta","Jet Phi Vs Eta",100,-5.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *jetPtVsPtHist = new TH2D("jetPtVsPt","Jet Pt Vs Parton Pt",500,0.,50.,500,0.,50.);
  TH2D *jetEtaVsEtaHist = new TH2D("jetEtaVsEta","Jet Eta Vs Parton Eta",100,-5.,5.,100,-5.,5.);
  TH2D *jetPhiVsPhiHist = new TH2D("jetPhiVsPhi","Jet Phi Vs Parton Phi",100,-1.0*TMath::Pi(),TMath::Pi(),100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *jetPhiVsPhiCutHist = new TH2D("jetPhiVsPhiCut","Jet Phi Vs Parton Phi (Jet eta < 2.5)",100,-1.0*TMath::Pi(),TMath::Pi(),100,-1.0*TMath::Pi(),TMath::Pi());

  /*
  TH2D *jetEventQ2VsXHist = new TH2D("jetEventQ2VsX","",12,-6.,0.,10,0.,5.);
  TH1D *jetEventXHist = new TH1D("jetEventXHist","",30,-6.,0.);

  TH2D *jetPhiVsEtaQ2xHist[35];
  for(int i=0; i<35; i++)
    {
      jetPhiVsEtaQ2xHist[i] = new TH2D(Form("jetPhiVsEtaQ2x_%d",i),"",100,-5.,5.,100,-1.0*TMath::Pi(),TMath::Pi());
    }

  TH2D *jetQ2VsXEtaHist[16];
  TH2D *jetPhiVsPhiEtaHist[16];
  for(int i=0; i<16; i++)
    {
      jetQ2VsXEtaHist[i] = new TH2D(Form("jetQ2VsXEta_%d",i),"",24,-6.,0.,20,0.,5.);
      jetPhiVsPhiEtaHist[i] = new TH2D(Form("jetPhiVsPhiEta_%d",i),"",100,-1.0*TMath::Pi(),TMath::Pi(),100,-1.0*TMath::Pi(),TMath::Pi());
    }
  */


  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC3::Pythia8ToHepMC3 topHepMC;

  // Specify file where HepMC events will be stored.
  //HepMC3::WriterAscii ascii_io(hepmcOut); // Write in HepMC3 Format
  HepMC3::WriterAsciiHepMC2 ascii_io(hepmcOut); // Write in HepMC2 Format for Delphes


  // Set Up Pythia Event
  Pythia8::Pythia p8;
  Pythia8::Event &event = p8.event;

  // A class to generate beam parameters according to own parametrization.
  BeamShapePtr myBeamShape = make_shared<eicBeamShape>(config,crab,hadE,lepE,xing);

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

      // Create FastJet Particle Containers
      vector<PseudoJet> particlesNoCut;
      //vector<PseudoJet> particlesCut;

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

      atan2PxPz1Hist->Fill(TMath::ATan2(p8.process[1].px(),p8.process[1].pz()));
      atan2PyPz1Hist->Fill(TMath::ATan2(p8.process[1].py(),p8.process[1].pz()));
      atan2PyPtot1Hist->Fill(TMath::ATan2(p8.process[1].py(),p8.process[1].pAbs()));
      
      vtxX->Fill(p8.process[0].xProd());
      vtxY->Fill(p8.process[0].yProd());
      vtxZ->Fill(p8.process[0].zProd());
      vtxT->Fill(p8.process[0].tProd());
      vtxYvsX->Fill(p8.process[0].xProd(),p8.process[0].yProd());
      vtxXvsT->Fill(p8.process[0].tProd(),p8.process[0].xProd());
      vtxXvsZ->Fill(p8.process[0].zProd(),p8.process[0].xProd());
      vtxYvsZ->Fill(p8.process[0].zProd(),p8.process[0].yProd());
      vtxTvsZ->Fill(p8.process[0].zProd(),p8.process[0].tProd());

      vtxXvsTZSum->Fill(p8.process[0].tProd()+p8.process[0].zProd(),p8.process[0].xProd());
      vtxXvsTZDiff->Fill(p8.process[0].tProd()-p8.process[0].zProd(),p8.process[0].xProd());

      double hadZ = p8.process[0].zProd() - TMath::Cos(0.0125)*p8.process[0].tProd();
      double lepZ = p8.process[0].zProd() + TMath::Cos(0.0125)*p8.process[0].tProd();
      lepVsHadPartZ->Fill(hadZ,lepZ);


      // Four-momenta of proton, electron, virtual photon/Z^0/W^+-.
      Pythia8::Vec4 pProton = event[1].p();
      Pythia8::Vec4 peIn    = event[4].p();
      //Pythia8::Vec4 peInAlt = event[2].p();
      Pythia8::Vec4 peOut   = event[6].p();
      Pythia8::Vec4 pPhoton = peIn - peOut;
      //Pythia8::Vec4 pPhotonAlt = peInAlt - peOut;

      // Q2, W2, Bjorken x, y.
      double Q2    = - pPhoton.m2Calc();
      //double Q2Alt = -pPhotonAlt.m2Calc();
      double W2    = (pProton + pPhoton).m2Calc();
      double x     = Q2 / (2. * pProton * pPhoton);
      //double xAlt  = Q2Alt / (2. * pProton * pPhotonAlt);
      double y     = (pProton * pPhoton) / (pProton * peIn);
      //double yAlt  = (pProton * pPhotonAlt) / (pProton * peInAlt);
      double jetPt = event[5].pT();
      double jetEta = event[5].eta();
      double jetPhi = event[5].phi();
      double jetP = event[5].pAbs();
      double elecPhi = event[6].phi();
      double jetKtPt = event[8].pT();
      double jetKtEta = event[8].eta();
      double jetKtPhi = event[8].phi();

      //if(Q2 < 10.0) cout << "Low Q2: " << Q2 << endl;

      //cout << setprecision(10) << eCMnom << " " << eCMnow << " " << Q2 << " " << x << " " << y << " " <<  p8.info.x1pdf() << " " << p8.info.x2pdf() << " " << p8.info.Q2Fac() << " " << p8.info.Q2Ren() << " " << p8.info.pT2Hat() << endl;
      //cout << setprecision(10) << Q2Alt << " " << xAlt << " " << yAlt << endl;
      //cout << setprecision(10) << eCMnom*eCMnom*x*y << " " << eCMnow*eCMnow*x*y << endl;
      //cout << endl;

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

	  double px = p8.event[i].px();
	  double py = p8.event[i].py();
	  double pz = p8.event[i].pz();
	  double E = p8.event[i].e();

	  //cout << i << " " << p8.event[i].id() << " " << partPt << " " << partEta << " " << partPhi << endl;

	  if(partFin) //partEta>-10.0 && partEta<10.0
	    {
	      partPtHist->Fill(partPt);
	      partPhiHist->Fill(partPhi);
	      partEtaHist->Fill(partEta);

	      partPtVsEtaHist->Fill(partEta,partPt);
	      partPhiVsEtaHist->Fill(partEta,partPhi);

	      partStatusVsEtaHist->Fill(partEta,p8.event[i].status());

	      if(p8.event[i].status() == 63)
		{
		  partPtVsEtaBRHist->Fill(partEta,partPt);
		  partPhiVsEtaBRHist->Fill(partEta,partPhi);
		}
	      if(p8.event[i].status() == 83 || p8.event[i].status() == 84)
		{
		  partPtVsEtaHadHist->Fill(partEta,partPt);
		  partPhiVsEtaHadHist->Fill(partEta,partPhi);
		}
	      if(p8.event[i].status() == 91)
		{
		  partPtVsEtaDecayHist->Fill(partEta,partPt);
		  partPhiVsEtaDecayHist->Fill(partEta,partPhi);
		}

	      if(partPt > 1.0)
		{
		  partPtHiHist->Fill(partPt);
		  partPhiHiHist->Fill(partPhi);
		  partEtaHiHist->Fill(partEta);
		  
		  partPhiVsEtaHiHist->Fill(partEta,partPhi);

		  partStatusVsEtaHiHist->Fill(partEta,p8.event[i].status());
		}
	    }

	  // Populate FastJet
	  if(partFin && partEta>-4.0 && partEta<4.0 && y<0.95 && y>0.01 && i > 7)
	    {
	      // No pT Cuts Particles for FastJet
	      fastjet::PseudoJet pNo(px,py,pz,E);
	      pNo.set_user_index(i);
	      particlesNoCut.push_back(pNo);
	    }
	}
      //cout << endl;


      if(jetP > 0.0)
	{
	  pJetEtaVsP->Fill(jetP,jetEta);

	  pJetQ2VsP->Fill(jetP,std::log10(Q2));

	  pJetPtHist->Fill(jetPt);
	  pJetEtaHist->Fill(jetEta);
	  pJetPhiHist->Fill(jetPhi);

	  pJetPtVsEtaHist->Fill(jetKtEta,jetKtPt);
	  pJetPhiVsEtaHist->Fill(jetKtEta,jetKtPhi);

	  pJetElecPhiHist->Fill(jetPhi - elecPhi - TMath::Pi());
	  pJetKtElecPhiHist->Fill(jetKtPhi - elecPhi - TMath::Pi());
	  pJetKtElecPhiAbsHist->Fill(TMath::Abs(jetKtPhi - elecPhi - TMath::Pi()));

	  //jetEventQ2VsXHist->Fill(std::log10(x),std::log10(Q2));
	  //jetEventXHist->Fill(std::log10(x));
	  //cout << "x = " << std::log10(x) << " Q2 = " << std::log10(Q2) << " Bin = " << jetEventQ2VsXHist->FindBin(std::log10(x),std::log10(Q2)) << endl;
	  //cout << "x = " << std::log10(x) << " Bin = " << jetEventXHist->FindBin(std::log10(x)) << endl;
	}

      // Cluster and Analyze Jets
      // Loop Over Jet Radii
      double R[2] = {1.0,0.4};
      for(int rad=0; rad<1; rad++)
	{
	  // Define Jet
	  JetDefinition jet_def_akt(antikt_algorithm,R[rad]);

	  // Cluster in Lab Frame
	  ClusterSequence cs_akt_lab_no(particlesNoCut, jet_def_akt);

	  // Min Jet Pt
	  double ptmin = 1.0;

	  // List of Jets
	  vector<PseudoJet> jets_akt_lab_no = sorted_by_pt(cs_akt_lab_no.inclusive_jets(ptmin));

	  // Loop Over Jets and Analyze
	  for(unsigned int jn=0; jn<jets_akt_lab_no.size(); jn++)
	    {
	      double fJetPt = jets_akt_lab_no[jn].pt();
	      double fJetEta = jets_akt_lab_no[jn].eta();
	      double fJetPhi = TVector2::Phi_mpi_pi(jets_akt_lab_no[jn].phi());

	      jetPtHist->Fill(fJetPt);
	      jetEtaHist->Fill(fJetEta);
	      jetPhiHist->Fill(fJetPhi);

	      jetPtVsPtNoCutHist->Fill(jetKtPt,fJetPt);

	      if(fJetPt > 5.0)
		{
		  jetPtVsEtaHist->Fill(fJetEta,fJetPt);
		  jetPhiVsEtaHist->Fill(fJetEta,fJetPhi);

		  jetPtVsPtHist->Fill(jetKtPt,fJetPt);
		  jetEtaVsEtaHist->Fill(jetKtEta,fJetEta);
		  jetPhiVsPhiHist->Fill(jetKtPhi,fJetPhi);

		  if(fJetEta < 2.5) jetPhiVsPhiCutHist->Fill(jetKtPhi,fJetPhi);

		  /*
		  int xQ2Index = xQ2Bin(std::log10(x),std::log10(Q2));

		  if(xQ2Index > -1) 
		    jetPhiVsEtaQ2xHist[xQ2Index]->Fill(fJetEta,fJetPhi);

		  int etaIndex = etaBin(fJetEta);
		  if(etaIndex > -1)
		    {
		      jetQ2VsXEtaHist[etaIndex]->Fill(std::log10(x),std::log10(Q2));
		      jetPhiVsPhiEtaHist[etaIndex]->Fill(jetKtPhi,fJetPhi);
		    }
		  */
		}

	      //if(jets_akt_lab_no[jn].pt() < 5.0) continue; // Relax Pt cut to see sub-threshold behavior
	      //if(jets_akt_lab_no[jn].eta() > 3.5-R[rad]) continue;
	      
	      //vector<PseudoJet> constituents = jets_akt_lab_no[jn].constituents();
	    }
	}

      // Construct new empty HepMC event and fill it.
      HepMC3::GenEvent hepmcevt;
      topHepMC.fill_next_event( p8, &hepmcevt );
      
      // Write the HepMC event to file.
      ascii_io.write_event(hepmcevt);
    }

  // List Statistics  
  p8.stat();
  
  // Write and Close Root File
  ofile->Write();
  ofile->Close();
  // Done.
  return 0;

					  
}
