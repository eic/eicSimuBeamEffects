// Create pythia eP collisions and save to HepMC
// based on example/main36.cc

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
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
#include "TMath.h"
#include "TVector2.h"

using namespace Pythia8;
using namespace fastjet;


int main(int argc, char* argv[])
{
  //const char* steerFile = 

  if(argc != 3)
    {
      cerr << "Wrong number of arguments" << endl;
      cerr << "program.exe steer out.hist.root" << endl;
      exit(EXIT_FAILURE);
    }

  const char* rootOut = argv[2];

  cout << "Steering File = " << argv[1] << endl;
  cout << "Root Output = " << rootOut << endl;

  // Open Root File
  TFile *ofile = TFile::Open(rootOut,"recreate");

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
  TH1D *partEtaHist = new TH1D("partEta","",400,-10.,10.);
  TH1D *partPhiHist = new TH1D("partPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsEtaHist = new TH2D("partPhiVsEta","",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());
  TH2D *partEtaVsQ2Hist = new TH2D("partEtaVsQ2","",50,0.,5.,400,-10.,10.);
  TH2D *partPhiVsQ2Hist = new TH2D("partPhiVsQ2","",50,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partStatusVsEtaHist = new TH2D("partStatusVsEta","",400,-10.,10.,200,0.,200.);

  // Hi Pt
  TH1D *partPtHiHist = new TH1D("partPtHi","",500,0.,50.);
  TH1D *partEtaHiHist = new TH1D("partEtaHi","",400,-10.,10.);
  TH1D *partPhiHiHist = new TH1D("partPhiHi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partPhiVsEtaHiHist = new TH2D("partPhiVsEtaHi","",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());
  TH2D *partEtaVsQ2HiHist = new TH2D("partEtaVsQ2Hi","",50,0.,5.,400,-10.,10.);
  TH2D *partPhiVsQ2HiHist = new TH2D("partPhiVsQ2Hi","",50,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *partStatusVsEtaHiHist = new TH2D("partStatusVsEtaHi","",400,-10.,10.,200,0.,200.);

  // "Jet" Quantities
  TH2D *pJetEtaVsP = new TH2D("pJetEtaVsP","",200,0.,200.,100,-5.,5.);

  TH2D *pJetQ2VsP = new TH2D("pJetQ2VsP","",200,0.,200.,100,0.,5.);

  TH1D *pJetPtHist = new TH1D("pJetPt","",500,0.,50.);
  TH1D *pJetEtaHist = new TH1D("pJetEta","",100,-5.,5.);
  TH1D *pJetPhiHist = new TH1D("pJetPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *pJetPhiVsEtaHist = new TH2D("pJetPhiVsEta","",400,-10.,10.,100,-1.0*TMath::Pi(),TMath::Pi());
  TH2D *pJetEtaVsQ2Hist = new TH2D("pJetEtaVsQ2","",50,0.,5.,400,-10.,10.);
  TH2D *pJetPhiVsQ2Hist = new TH2D("pJetPhiVsQ2","",50,0.,5.,100,-1.0*TMath::Pi(),TMath::Pi());

  // "Jet" - Electron Quantities
  TH1D *pJetElecPhiHist = new TH1D("pJetElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *pJetKtElecPhiHist = new TH1D("pJetKtElecPhi","",100,-1.0*TMath::Pi(),TMath::Pi());
  TH1D *pJetKtElecPhiAbsHist = new TH1D("pJetKtElecPhiAbs","",100,-1.0*TMath::Pi(),TMath::Pi());

  // True Jet Quantities
  TH1D *jetPtHist = new TH1D("jetPt","",500,0.,50.);
  TH1D *jetEtaHist = new TH1D("jetEta","",100,-5.,5.);
  TH1D *jetPhiHist = new TH1D("jetPhi","",100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *jetEtaVsEtaHist = new TH2D("jetEtaVsEta","",100,-5.,5.,100,-5.,5.);
  TH2D *jetPhiVsPhiHist = new TH2D("jetPhiVsPhi","",100,-1.0*TMath::Pi(),TMath::Pi(),100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *jetPhiVsPhiCutHist = new TH2D("jetPhiVsPhiCut","",100,-1.0*TMath::Pi(),TMath::Pi(),100,-1.0*TMath::Pi(),TMath::Pi());

  TH2D *jetPhiVsEtaHist = new TH2D("jetPhiVsEta","",100,-5.,5.,100,-1.0*TMath::Pi(),TMath::Pi());


  // Set Up Pythia Event
  Pythia8::Pythia p8;
  Pythia8::Event &event = p8.event;

  // A class to generate beam parameters according to own parametrization.
  BeamShapePtr myBeamShape = make_shared<eicBeamShape>(275.0,18.0,0.025,0.0030,-0.0015);

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

      //if(Q2 < 10.0) cout << "Low Q2: " << Q2 << endl;

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

	  if(partFin && partEta>-10.0 && partEta<10.0 && y<0.95 && y>0.01 && i > 7)
	    {
	      partPtHist->Fill(partPt);
	      partPhiHist->Fill(partPhi);
	      partEtaHist->Fill(partEta);

	      partPhiVsEtaHist->Fill(partEta,partPhi);
	      partEtaVsQ2Hist->Fill(std::log10(Q2),partEta);
	      partPhiVsQ2Hist->Fill(std::log10(Q2),partPhi);

	      partStatusVsEtaHist->Fill(partEta,p8.event[i].status());

	      if(partPt > 1.0)
		{
		  partPtHiHist->Fill(partPt);
		  partPhiHiHist->Fill(partPhi);
		  partEtaHiHist->Fill(partEta);
		  
		  partPhiVsEtaHiHist->Fill(partEta,partPhi);
		  partEtaVsQ2HiHist->Fill(std::log10(Q2),partEta);
		  partPhiVsQ2HiHist->Fill(std::log10(Q2),partPhi);

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

      if(y < 0.95 && y > 0.01 && jetP > 0.0)
	{
	  pJetEtaVsP->Fill(jetP,jetEta);

	  pJetQ2VsP->Fill(jetP,std::log10(Q2));

	  pJetPtHist->Fill(jetPt);
	  pJetEtaHist->Fill(jetEta);
	  pJetPhiHist->Fill(jetPhi);

	  pJetPhiVsEtaHist->Fill(jetEta,jetPhi);
	  pJetEtaVsQ2Hist->Fill(std::log10(Q2),jetEta);
	  pJetPhiVsQ2Hist->Fill(std::log10(Q2),jetPhi);

	  pJetElecPhiHist->Fill(jetPhi - elecPhi - TMath::Pi());
	  pJetKtElecPhiHist->Fill(jetKtPhi - elecPhi - TMath::Pi());
	  pJetKtElecPhiAbsHist->Fill(TMath::Abs(jetKtPhi - elecPhi - TMath::Pi()));
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
	      double fJetEta = jets_akt_lab_no[jn].eta();
	      double fJetPhi = TVector2::Phi_mpi_pi(jets_akt_lab_no[jn].phi());

	      jetPtHist->Fill(jets_akt_lab_no[jn].pt());
	      jetEtaHist->Fill(fJetEta);
	      jetPhiHist->Fill(fJetPhi);

	      if(jets_akt_lab_no[jn].pt() > 5.0)
		{
		  jetEtaVsEtaHist->Fill(jetEta,fJetEta);
		  jetPhiVsPhiHist->Fill(jetPhi,fJetPhi);

		  jetPhiVsEtaHist->Fill(fJetEta,fJetPhi);

		  if(fJetEta < 3.5) jetPhiVsPhiCutHist->Fill(jetPhi,fJetPhi);
		}

	      //if(jets_akt_lab_no[jn].pt() < 5.0) continue; // Relax Pt cut to see sub-threshold behavior
	      //if(jets_akt_lab_no[jn].eta() > 3.5-R[rad]) continue;
	      
	      //vector<PseudoJet> constituents = jets_akt_lab_no[jn].constituents();
	    }
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
