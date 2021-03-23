TFile *fa=0;
TFile *fb=0;
TFile *fc=0;
TFile *fd=0;
TFile *fe=0;
TFile *ff=0;

TCanvas *c[30];

//TString outFileName;

int plotBeamShape()
{

  //gStyle->SetPalette(1,0);
  //gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  //gStyle->SetOptFit(1111);

  TString rootHistFnameA;
  TString rootHistFnameB;
  TString rootHistFnameC;
  TString rootHistFnameD;
  TString rootHistFnameE;
  TString rootHistFnameF;

  rootHistFnameA="rootFiles/testDefault.hist.root";
  rootHistFnameB="rootFiles/testXing.hist.root";
  rootHistFnameC="rootFiles/testCrab.hist.root";
  rootHistFnameD="rootFiles/testDiv.hist.root";
  rootHistFnameE="rootFiles/testCrabDiv.hist.root";
  rootHistFnameF="rootFiles/testAll.hist.root";
  
  fa=new TFile(rootHistFnameA);
  assert(fa->IsOpen());

  fb=new TFile(rootHistFnameB);
  assert(fb->IsOpen());

  fc=new TFile(rootHistFnameC);
  assert(fc->IsOpen());

  fd=new TFile(rootHistFnameD);
  assert(fd->IsOpen());

  fe=new TFile(rootHistFnameE);
  assert(fe->IsOpen());

  ff=new TFile(rootHistFnameF);
  assert(ff->IsOpen());

  //outFileName = outputPostscriptFile;

  //c=new TCanvas("c1","Canvas",800,600);

  //printf("To put all plots into pdf, type plots2pdf()\n");
  //printf("To view individual page, type eachPage(#,0), where #=1-5\n");

  return 0;
}


void plotShapes()
{
  c[0]=new TCanvas("c0","CM Energy Difference",800,600);
  c[1]=new TCanvas("c1","Beam Z Momentum",1200,600);
  c[2]=new TCanvas("c2","Z-Vertex",800,600);
  c[3]=new TCanvas("c3","Hadron Beam Py Vs Px",800,600);
  c[4]=new TCanvas("c4","Lepton Beam Py Vs Px",800,600);
  c[5]=new TCanvas("c5","Hadron Beam XYZ",800,600);
  c[6]=new TCanvas("c6","Lepton Beam XYZ",800,600);
  c[7]=new TCanvas("c7","Particle Kinematics",800,600);

  c[0]->Clear();
  c[0]->Divide(1,1);

  c[1]->Clear();
  c[1]->Divide(2,1);

  c[2]->Clear();
  c[2]->Divide(1,1);

  c[3]->Clear();
  c[3]->Divide(2,2);

  c[4]->Clear();
  c[4]->Divide(2,2);

  c[5]->Clear();
  c[5]->Divide(2,2);

  c[6]->Clear();
  c[6]->Divide(2,2);

  c[7]->Clear();
  c[7]->Divide(2,2);

  TH1D *hVtx=(TH1D *)fb->Get("vtxZ");

  TH1D *hECM[6];
  TH1D *hPZHad[6];
  TH1D *hPZLep[6];
  TH2D *hPXYHad[6];
  TH2D *hPXYLep[6];
  TH2D *hPXYZHad[2][6];
  TH2D *hPXYZLep[2][6];
  TH1D *hPartKin[4][6];

  // Center of Mass
  hECM[0]=(TH1D *)fa->Get("eCM"); // Default
  hECM[1]=(TH1D *)fb->Get("eCM"); // Xing
  hECM[2]=(TH1D *)fc->Get("eCM"); // Crab
  hECM[3]=(TH1D *)fd->Get("eCM"); // Divergence
  hECM[4]=(TH1D *)fe->Get("eCM"); // Crab + Divergence
  hECM[5]=(TH1D *)ff->Get("eCM"); // All

  // Pz
  hPZHad[0]=(TH1D *)fa->Get("pZ1"); // Default
  hPZHad[1]=(TH1D *)fb->Get("pZ1"); // Xing
  hPZHad[2]=(TH1D *)fc->Get("pZ1"); // Crab
  hPZHad[3]=(TH1D *)fd->Get("pZ1"); // Divergence
  hPZHad[4]=(TH1D *)fe->Get("pZ1"); // Crab + Divergence
  hPZHad[5]=(TH1D *)ff->Get("pZ1"); // All

  hPZLep[0]=(TH1D *)fa->Get("pZ2"); // Default
  hPZLep[1]=(TH1D *)fb->Get("pZ2"); // Xing
  hPZLep[2]=(TH1D *)fc->Get("pZ2"); // Crab
  hPZLep[3]=(TH1D *)fd->Get("pZ2"); // Divergence
  hPZLep[4]=(TH1D *)fe->Get("pZ2"); // Crab + Divergence
  hPZLep[5]=(TH1D *)ff->Get("pZ2"); // All

  // Py Vs Px
  hPXYHad[0]=(TH2D *)fa->Get("pXY1"); // Default
  hPXYHad[1]=(TH2D *)fb->Get("pXY1"); // Xing
  hPXYHad[2]=(TH2D *)fc->Get("pXY1"); // Crab
  hPXYHad[3]=(TH2D *)fd->Get("pXY1"); // Divergence
  hPXYHad[4]=(TH2D *)fe->Get("pXY1"); // Crab + Divergence
  hPXYHad[5]=(TH2D *)ff->Get("pXY1"); // All

  hPXYLep[0]=(TH2D *)fa->Get("pXY2"); // Default
  hPXYLep[1]=(TH2D *)fb->Get("pXY2"); // Xing
  hPXYLep[2]=(TH2D *)fc->Get("pXY2"); // Crab
  hPXYLep[3]=(TH2D *)fd->Get("pXY2"); // Divergence
  hPXYLep[4]=(TH2D *)fe->Get("pXY2"); // Crab + Divergence
  hPXYLep[5]=(TH2D *)ff->Get("pXY2"); // All

  // Px Vs Z and Py Vs Z
  hPXYZHad[0][0]=(TH2D *)fa->Get("pXZProd1"); // Default
  hPXYZHad[0][1]=(TH2D *)fb->Get("pXZProd1"); // Xing
  hPXYZHad[0][2]=(TH2D *)fc->Get("pXZProd1"); // Crab
  hPXYZHad[0][3]=(TH2D *)fd->Get("pXZProd1"); // Divergence
  hPXYZHad[0][4]=(TH2D *)fe->Get("pXZProd1"); // Crab + Divergence
  hPXYZHad[0][5]=(TH2D *)ff->Get("pXZProd1"); // All

  hPXYZHad[1][0]=(TH2D *)fa->Get("pYZProd1"); // Default
  hPXYZHad[1][1]=(TH2D *)fb->Get("pYZProd1"); // Xing
  hPXYZHad[1][2]=(TH2D *)fc->Get("pYZProd1"); // Crab
  hPXYZHad[1][3]=(TH2D *)fd->Get("pYZProd1"); // Divergence
  hPXYZHad[1][4]=(TH2D *)fe->Get("pYZProd1"); // Crab + Divergence
  hPXYZHad[1][5]=(TH2D *)ff->Get("pYZProd1"); // All

  hPXYZLep[0][0]=(TH2D *)fa->Get("pXZProd2"); // Default
  hPXYZLep[0][1]=(TH2D *)fb->Get("pXZProd2"); // Xing
  hPXYZLep[0][2]=(TH2D *)fc->Get("pXZProd2"); // Crab
  hPXYZLep[0][3]=(TH2D *)fd->Get("pXZProd2"); // Divergence
  hPXYZLep[0][4]=(TH2D *)fe->Get("pXZProd2"); // Crab + Divergence
  hPXYZLep[0][5]=(TH2D *)ff->Get("pXZProd2"); // All

  hPXYZLep[1][0]=(TH2D *)fa->Get("pYZProd2"); // Default
  hPXYZLep[1][1]=(TH2D *)fb->Get("pYZProd2"); // Xing
  hPXYZLep[1][2]=(TH2D *)fc->Get("pYZProd2"); // Crab
  hPXYZLep[1][3]=(TH2D *)fd->Get("pYZProd2"); // Divergence
  hPXYZLep[1][4]=(TH2D *)fe->Get("pYZProd2"); // Crab + Divergence
  hPXYZLep[1][5]=(TH2D *)ff->Get("pYZProd2"); // All

  // Particle Kinematics
  hPartKin[0][0]=(TH1D *)fa->Get("partEta");
  hPartKin[0][1]=(TH1D *)fb->Get("partEta");
  hPartKin[0][2]=(TH1D *)fc->Get("partEta");
  hPartKin[0][3]=(TH1D *)fd->Get("partEta");
  hPartKin[0][4]=(TH1D *)fe->Get("partEta");
  hPartKin[0][5]=(TH1D *)ff->Get("partEta");

  hPartKin[1][0]=(TH1D *)fa->Get("partEtaHi");
  hPartKin[1][1]=(TH1D *)fb->Get("partEtaHi");
  hPartKin[1][2]=(TH1D *)fc->Get("partEtaHi");
  hPartKin[1][3]=(TH1D *)fd->Get("partEtaHi");
  hPartKin[1][4]=(TH1D *)fe->Get("partEtaHi");
  hPartKin[1][5]=(TH1D *)ff->Get("partEtaHi");

  hPartKin[2][0]=(TH1D *)fa->Get("partPhi");
  hPartKin[2][1]=(TH1D *)fb->Get("partPhi");
  hPartKin[2][2]=(TH1D *)fc->Get("partPhi");
  hPartKin[2][3]=(TH1D *)fd->Get("partPhi");
  hPartKin[2][4]=(TH1D *)fe->Get("partPhi");
  hPartKin[2][5]=(TH1D *)ff->Get("partPhi");

  hPartKin[3][0]=(TH1D *)fa->Get("partPhiHi");
  hPartKin[3][1]=(TH1D *)fb->Get("partPhiHi");
  hPartKin[3][2]=(TH1D *)fc->Get("partPhiHi");
  hPartKin[3][3]=(TH1D *)fd->Get("partPhiHi");
  hPartKin[3][4]=(TH1D *)fe->Get("partPhiHi");
  hPartKin[3][5]=(TH1D *)ff->Get("partPhiHi");

  for(int i=0; i<6; i++)
    {
      hECM[i]->Rebin(10);
      hPZHad[i]->Rebin(10);
      hPZLep[i]->Rebin(10);
    }

  c[0]->cd(1);
  hECM[1]->Draw("HIST");
  hECM[4]->SetLineColor(kGreen+2);
  hECM[4]->Draw("HISTSAME");
  hECM[5]->SetLineColor(kRed);
  hECM[5]->Draw("HISTSAME");
  hECM[0]->SetLineColor(kBlack);
  hECM[0]->Draw("HISTSAME");
  hECM[1]->SetTitle("Nominal - Modified CM Energy;E_{Nom} - E_{Mod}");

  TLegend *leg0;
  leg0 = new TLegend(0.60,0.75,0.90,0.50);
  leg0->SetFillColor(0);
  leg0->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg0->AddEntry(hECM[0],"Default","l");
  leg0->AddEntry(hECM[1],"Beam Xing Only","l");
  leg0->AddEntry(hECM[4],"Crab + Divergence Only","l");
  leg0->AddEntry(hECM[5],"All Effects","l");
  leg0->Draw();

  c[1]->cd(1);
  hPZHad[1]->Draw("HIST");
  hPZHad[4]->SetLineColor(kGreen+2);
  hPZHad[4]->Draw("HISTSAME");
  hPZHad[5]->SetLineColor(kRed);
  hPZHad[5]->Draw("HISTSAME");
  hPZHad[0]->SetLineColor(kBlack);
  hPZHad[0]->Draw("HISTSAME");
  hPZHad[1]->SetTitle("Hadron Beam Z-Momentum;p_{z} Hadron");
  hPZHad[1]->GetXaxis()->SetRangeUser(274.6,275.4);

  TLegend *leg1;
  leg1 = new TLegend(0.60,0.75,0.90,0.50);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg1->AddEntry(hPZHad[0],"Default","l");
  leg1->AddEntry(hPZHad[1],"Beam Xing Only","l");
  leg1->AddEntry(hPZHad[4],"Crab + Divergence Only","l");
  leg1->AddEntry(hPZHad[5],"All Effects","l");
  leg1->Draw();

  c[1]->cd(2);
  hPZLep[1]->Draw("HIST");
  hPZLep[4]->SetLineColor(kGreen+2);
  hPZLep[4]->Draw("HISTSAME");
  hPZLep[5]->SetLineColor(kRed);
  hPZLep[5]->Draw("HISTSAME");
  hPZLep[0]->SetLineColor(kBlack);
  hPZLep[0]->Draw("HISTSAME");
  hPZLep[1]->SetTitle("Lepton Beam Z-Momentum;p_{z} Lepton");
  hPZLep[1]->GetXaxis()->SetRangeUser(-18.02,-17.98);

  c[2]->cd(1);
  hVtx->Draw("HIST");
  hVtx->SetTitle("Smeared Vertex Z Position;Z-Vertex");

  c[3]->cd(1);
  hPXYHad[1]->Draw("COLZ");
  hPXYHad[1]->SetTitle("Py Vs Px: Xing Angle Only;Px;Py");
  hPXYHad[1]->GetXaxis()->SetRangeUser(6.5,7.5);
  hPXYHad[1]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[3]->cd(2);
  hPXYHad[2]->Draw("COLZ");
  hPXYHad[2]->SetTitle("Py Vs Px: Crabbing Only;Px;Py");
  hPXYHad[2]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYHad[2]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[3]->cd(3);
  hPXYHad[3]->Draw("COLZ");
  hPXYHad[3]->SetTitle("Py Vs Px: Divergence Only;Px;Py");
  hPXYHad[3]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYHad[3]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[3]->cd(4);
  hPXYHad[5]->Draw("COLZ");
  hPXYHad[5]->SetTitle("Py Vs Px: All Effects;Px;Py");
  hPXYHad[5]->GetXaxis()->SetRangeUser(6.5,7.5);
  hPXYHad[5]->GetYaxis()->SetRangeUser(-0.5,0.5);

  c[4]->cd(1);
  hPXYLep[1]->Draw("COLZ");
  hPXYLep[1]->SetTitle("Py Vs Px: Xing Angle Only;Px;Py");
  hPXYLep[1]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYLep[1]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[4]->cd(2);
  hPXYLep[2]->Draw("COLZ");
  hPXYLep[2]->SetTitle("Py Vs Px: Crabbing Only;Px;Py");
  hPXYLep[2]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYLep[2]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[4]->cd(3);
  hPXYLep[3]->Draw("COLZ");
  hPXYLep[3]->SetTitle("Py Vs Px: Divergence Only;Px;Py");
  hPXYLep[3]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYLep[3]->GetYaxis()->SetRangeUser(-0.5,0.5);
  c[4]->cd(4);
  hPXYLep[5]->Draw("COLZ");
  hPXYLep[5]->SetTitle("Py Vs Px: All Effects;Px;Py");
  hPXYLep[5]->GetXaxis()->SetRangeUser(-0.5,0.5);
  hPXYLep[5]->GetYaxis()->SetRangeUser(-0.5,0.5);

  // Hadron Px Py Vs Z-Vertex
  c[5]->cd(1);
  hPXYZHad[0][2]->Draw("COLZ");
  hPXYZHad[0][2]->SetTitle("Hadron Px Vs Z-Vertex: Crabbing Only;Z-Vertex;Px");
  hPXYZHad[0][2]->GetYaxis()->SetRangeUser(-1.,1.);
  c[5]->cd(2);
  hPXYZHad[0][5]->Draw("COLZ");
  hPXYZHad[0][5]->SetTitle("Hadron Px Vs Z-Vertex: All Effects;Z-Vertex;Px");
  hPXYZHad[0][5]->GetYaxis()->SetRangeUser(6.,8.);
  c[5]->cd(3);
  hPXYZHad[1][2]->Draw("COLZ");
  hPXYZHad[1][2]->SetTitle("Hadron Py Vs Z-Vertex: Crabbing Only;Z-Vertex;Py");
  hPXYZHad[1][2]->GetYaxis()->SetRangeUser(-1.,1.);
  c[5]->cd(4);
  hPXYZHad[1][5]->Draw("COLZ");
  hPXYZHad[1][5]->SetTitle("Hadron Py Vs Z-Vertex: All Effects;Z-Vertex;Py");
  hPXYZHad[1][5]->GetYaxis()->SetRangeUser(-1.,1.);

  // Lepton Px Py Vs Z-Vertex
  c[6]->cd(1);
  hPXYZLep[0][2]->Draw("COLZ");
  hPXYZLep[0][2]->SetTitle("Lepton Px Vs Z-Vertex: Crabbing Only;Z-Vertex;Px");
  hPXYZLep[0][2]->GetYaxis()->SetRangeUser(-1.,1.);
  c[6]->cd(2);
  hPXYZLep[0][5]->Draw("COLZ");
  hPXYZLep[0][5]->SetTitle("Lepton Px Vs Z-Vertex: All Effects;Z-Vertex;Px");
  hPXYZLep[0][5]->GetYaxis()->SetRangeUser(-1.,1.);
  c[6]->cd(3);
  hPXYZLep[1][2]->Draw("COLZ");
  hPXYZLep[1][2]->SetTitle("Lepton Py Vs Z-Vertex: Crabbing Only;Z-Vertex;Py");
  hPXYZLep[1][2]->GetYaxis()->SetRangeUser(-1.,1.);
  c[6]->cd(4);
  hPXYZLep[1][5]->Draw("COLZ");
  hPXYZLep[1][5]->SetTitle("Lepton Py Vs Z-Vertex: All Effects;Z-Vertex;Py");
  hPXYZLep[1][5]->GetYaxis()->SetRangeUser(-1.,1.);

  // Particle Kinematics
  c[7]->cd(1);
  hPartKin[0][0]->Draw("HIST");
  hPartKin[0][0]->SetLineColor(kBlack);
  hPartKin[0][4]->Draw("HISTSAME");
  hPartKin[0][4]->SetLineColor(kGreen+2);
  hPartKin[0][5]->Draw("HISTSAME");
  hPartKin[0][5]->SetLineColor(kRed);
  hPartKin[0][0]->SetTitle("Particle Eta (pT > 0.0);Eta");
  hPartKin[0][0]->GetYaxis()->SetRangeUser(0,300000.0);
  c[7]->cd(2);
  hPartKin[1][0]->Draw("HIST");
  hPartKin[1][0]->SetLineColor(kBlack);
  hPartKin[1][4]->Draw("HISTSAME");
  hPartKin[1][4]->SetLineColor(kGreen+2);
  hPartKin[1][5]->Draw("HISTSAME");
  hPartKin[1][5]->SetLineColor(kRed);
  hPartKin[1][0]->SetTitle("Particle Eta (pT > 1.0);Eta");
  hPartKin[1][0]->GetYaxis()->SetRangeUser(0,40000.0);
  c[7]->cd(3);
  hPartKin[2][0]->Draw("HIST");
  hPartKin[2][0]->SetLineColor(kBlack);
  hPartKin[2][4]->Draw("HISTSAME");
  hPartKin[2][4]->SetLineColor(kGreen+2);
  hPartKin[2][5]->Draw("HISTSAME");
  hPartKin[2][5]->SetLineColor(kRed);
  hPartKin[2][0]->SetTitle("Particle Phi (pT > 0.0);Phi");
  hPartKin[2][0]->GetYaxis()->SetRangeUser(0,140000.0);

  TLegend *leg7;
  leg7 = new TLegend(0.30,0.55,0.70,0.15);
  leg7->SetFillColor(0);
  leg7->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg7->AddEntry(hPartKin[2][0],"Default","l");
  leg7->AddEntry(hPartKin[2][4],"Crab + Divergence Only","l");
  leg7->AddEntry(hPartKin[2][5],"All Effects","l");
  leg7->Draw();

  c[7]->cd(4);
  hPartKin[3][0]->Draw("HIST");
  hPartKin[3][0]->SetLineColor(kBlack);
  hPartKin[3][4]->Draw("HISTSAME");
  hPartKin[3][4]->SetLineColor(kGreen+2);
  hPartKin[3][5]->Draw("HISTSAME");
  hPartKin[3][5]->SetLineColor(kRed);
  hPartKin[3][0]->SetTitle("Particle Phi (pT > 1.0);Phi");
  hPartKin[3][0]->GetYaxis()->SetRangeUser(0,20000.0);
}