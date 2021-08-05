// Create pythia eP collisions and save to HepMC
// based on example/main36.cc

#include "Pythia8/Pythia.h"
#include "Pythia8/BeamShape.h"
#include "Pythia8Plugins/HepMC3.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "eicBeamShape.h"


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


int main(int argc, char* argv[])
{
  if(argc != 8)
    {
      cerr << "Wrong number of arguments" << endl;
      cerr << "program.exe steer configuration hadronE leptonE xangle out.hist.root out.hepmc" << endl;
      exit(EXIT_FAILURE);
    }

  const int config = atoi(argv[2]);
  const double hadE = atof(argv[3]);
  const double lepE = atof(argv[4]);
  const double xing = atof(argv[5]);
  const char* rootOut = argv[6];
  const char* hepmcOut = argv[7];

  cout << "Steering File = " << argv[1] << endl;
  if(config == 1) cout << "Configuration = High Divergence" << endl;
  if(config == 2) cout << "Configuration = High Acceptance" << endl; 
  if(config == 3) cout << "Ion Running" << endl;
  cout << "Hadron Energy = " << hadE << endl;
  cout << "Lepton Energy = " << lepE << endl;
  cout << "Beam Crossing Angle (Radians) = " << xing << endl;
  cout << "Root Output = " << rootOut << endl;
  cout << "HepMC Output = " << hepmcOut << endl;

  // Open Root File
  TFile *ofile = TFile::Open(rootOut,"recreate");

  // Histos
  TH1D *q2Hist = new TH1D("q2Hist","",100,0.,5.);
  TH1D *xHist = new TH1D("xHist","",100,-4.,1.);
  TH1D *yHist = new TH1D("yHist","",1000,-5.,0.);
  TH2D *phaseSpaceHist = new TH2D("phaseSpaceHist","",100,-4.,1.,100,0.,5.);


  // Interface for conversion from Pythia8::Event to HepMC event.
  HepMC3::Pythia8ToHepMC3 topHepMC;

  // Specify file where HepMC events will be stored.
  //HepMC3::WriterAscii ascii_io(hepmcOut); // Write in HepMC3 Format
  HepMC3::WriterAsciiHepMC2 ascii_io(hepmcOut); // Write in HepMC2 Format for Delphes

  // Set Up Pythia Event
  Pythia8::Pythia p8;
  Pythia8::Event &event = p8.event;

  // A class to generate beam parameters according to own parametrization.
  BeamShapePtr myBeamShape = make_shared<eicBeamShape>(config,hadE,lepE,xing);

  // Hand pointer to Pythia.
  // If you comment this out you get internal Gaussian-style implementation.
  p8.setBeamShapePtr(myBeamShape);

  // Read in Steering File & Make Other Settings
  p8.readFile(argv[1]);
  p8.readString("Main:timesAllowErrors = 10000"); // allow more errors, eP is brittle

  // Initialize Pythia
  p8.init();

  // Run
  int nevents = p8.mode("Main:numberOfEvents");
  if(nevents > 10000)
    {
      cout << "Limit nEvents to < 10000 due to size of HepMC file" << endl;
      exit(EXIT_FAILURE);
    }
  for(int ev=0; ev<nevents; ev++)
    {
      if(!p8.next()) continue;
      // p8.event.list();

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

      // Event Level
      q2Hist->Fill(std::log10(Q2));
      xHist->Fill(std::log10(x));
      yHist->Fill(std::log10(y));
      phaseSpaceHist->Fill(std::log10(x),std::log10(Q2));

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
