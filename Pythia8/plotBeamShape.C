TFile *fa=0;
TFile *fb=0;
TFile *fc=0;
TFile *fd=0;
TFile *fe=0;
TFile *ff=0;
TFile *wrk=0;

TCanvas *c[30];

//TString outFileName;

int plotBeamShape(TString path="rootFiles/testDefault_new.hist.root")
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
  TString rootHistFnameWrk;

  rootHistFnameA="rootFiles/testDefault_new.hist.root";
  rootHistFnameB="rootFiles/testXing.hist.root";
  rootHistFnameC="rootFiles/testCrab.hist.root";
  rootHistFnameD="rootFiles/testDiv.hist.root";
  rootHistFnameE="rootFiles/testCrabDiv.hist.root";
  rootHistFnameF="rootFiles/testAll_new.hist.root";
  rootHistFnameWrk=path;
  
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

  wrk=new TFile(rootHistFnameWrk);
  assert(wrk->IsOpen());

  //outFileName = outputPostscriptFile;

  //c=new TCanvas("c1","Canvas",800,600);

  //printf("To put all plots into pdf, type plots2pdf()\n");
  //printf("To view individual page, type eachPage(#,0), where #=1-5\n");

  return 0;
}


void plotSingleConfigBeamMom()
{
  c[0]=new TCanvas("c0","CM Energy and Beam Momentum Diff",800,600);
  c[1]=new TCanvas("c1","Hadron Beam Energy Deviations",800,600);
  c[2]=new TCanvas("c2","Lepton Beam Energy Deviations",800,600);

  c[0]->Clear();
  c[0]->Divide(2,2);

  c[1]->Clear();
  c[1]->Divide(2,2);

  c[2]->Clear();
  c[2]->Divide(2,2);

  TH1D *hA=(TH1D *)wrk->Get("eCM");
  TH1D *hB=(TH1D *)wrk->Get("pZ1");
  TH1D *hC=(TH1D *)wrk->Get("pZ2");

  TH2D *hD=(TH2D *)wrk->Get("pXY1");
  TH2D *hE=(TH2D *)wrk->Get("pXZProd1");
  TH2D *hF=(TH2D *)wrk->Get("pYZProd1");

  TH2D *hG=(TH2D *)wrk->Get("pXY2");
  TH2D *hH=(TH2D *)wrk->Get("pXZProd2");
  TH2D *hI=(TH2D *)wrk->Get("pYZProd2");

  c[0]->cd(1);
  hA->Rebin(5);
  hA->Draw("HIST");
  c[0]->cd(2);
  hB->Draw("HIST");
  hB->SetXTitle("Hadron Pz [GeV]");
  c[0]->cd(3);
  hC->Draw("HIST");
  hC->SetXTitle("Lepton Pz [GeV]");

  c[1]->cd(1);
  hD->Draw("COLZ");
  hD->GetXaxis()->SetRangeUser(6,8);
  hD->GetYaxis()->SetRangeUser(-1,1);
  hD->SetTitle("Hadron Beam Py vs Px;Px [GeV];Py [GeV]");
  gPad->SetLogz();
  c[1]->cd(2);
  hE->Draw("COLZ");
  hE->GetXaxis()->SetRangeUser(-300,300);
  hE->GetYaxis()->SetRangeUser(6,8);
  hE->SetTitle("Hadron Beam Px vs Z Vertex;Z-Vert [mm];Px [GeV]");
  gPad->SetLogz();
  c[1]->cd(3);
  hF->Draw("COLZ");
  hF->GetXaxis()->SetRangeUser(-300,300);
  hF->GetYaxis()->SetRangeUser(-1,1);
  hF->SetTitle("Hadron Beam Py vs Z Vertex;Z-Vert [mm];Py [GeV]");
  gPad->SetLogz();

  c[2]->cd(1);
  hG->Draw("COLZ");
  hG->GetXaxis()->SetRangeUser(-0.1,0.1);
  hG->GetYaxis()->SetRangeUser(-0.1,0.1);
  hG->SetTitle("Lepton Beam Py vs Px;Px [GeV];Py [GeV]");
  gPad->SetLogz();
  c[2]->cd(2);
  hH->Draw("COLZ");
  hH->GetXaxis()->SetRangeUser(-300,300);
  hH->GetYaxis()->SetRangeUser(-0.1,0.1);
  hH->SetTitle("Lepton Beam Px vs Z Vertex;Z-Vert [mm];Px [GeV]");
  gPad->SetLogz();
  c[2]->cd(3);
  hI->Draw("COLZ");
  hI->GetXaxis()->SetRangeUser(-300,300);
  hI->GetYaxis()->SetRangeUser(-0.1,0.1);
  hI->SetTitle("Lepton Beam Py vs Z Vertex;Z-Vert [mm];Py [GeV]");
  gPad->SetLogz();
}


void plotSingleConfigVertex()
{
  c[0]=new TCanvas("c0","Vertex Distributions",800,600);
  c[1]=new TCanvas("c1","Vertex Correlations",800,600);
  c[2]=new TCanvas("c2","Beam Sizes and X-TZ Correlations",800,600);
  c[3]=new TCanvas("c3","Beam Angles",800,600);

  c[0]->Clear();
  c[0]->Divide(2,2);

  c[1]->Clear();
  c[1]->Divide(2,2);

  c[2]->Clear();
  c[2]->Divide(2,2);

  c[3]->Clear();
  c[3]->Divide(2,2);

  TH1D *hA=(TH1D *)wrk->Get("vtxX");
  TH1D *hB=(TH1D *)wrk->Get("vtxY");
  TH1D *hC=(TH1D *)wrk->Get("vtxZ");
  TH1D *hD=(TH1D *)wrk->Get("vtxT");

  TH2D *hE=(TH2D *)wrk->Get("vtxYvsX");
  TH2D *hF=(TH2D *)wrk->Get("vtxXvsZ");
  TH2D *hG=(TH2D *)wrk->Get("vtxXvsT");
  TH2D *hH=(TH2D *)wrk->Get("vtxTvsZ");

  TH2D *hI=(TH2D *)wrk->Get("lepVsHadPartZ");
  TH2D *hJ=(TH2D *)wrk->Get("vtxXvsTZDiff");
  TH2D *hK=(TH2D *)wrk->Get("vtxXvsTZSum");

  TH1D *hL=(TH1D *)wrk->Get("atan2PxPz1");
  TH1D *hM=(TH1D *)wrk->Get("atan2PyPz1");
  TH1D *hN=(TH1D *)wrk->Get("atan2PxPz2");
  TH1D *hO=(TH1D *)wrk->Get("atan2PyPz2");

  c[0]->cd(1);
  hA->Draw("HIST");
  hA->GetXaxis()->SetRangeUser(-1,1);
  c[0]->cd(2);
  hB->Draw("HIST");
  hB->GetXaxis()->SetRangeUser(-0.1,0.1);
  c[0]->cd(3);
  hC->Draw("HIST");
  hC->GetXaxis()->SetRangeUser(-300,300);
  c[0]->cd(4);
  hD->Draw("HIST");
  hD->GetXaxis()->SetRangeUser(-300,300);

  c[1]->cd(1);
  hE->Draw("COLZ");
  hE->GetXaxis()->SetRangeUser(-1,1);
  hE->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();
  c[1]->cd(2);
  hF->Draw("COLZ");
  hF->GetXaxis()->SetRangeUser(-300,300);
  hF->GetYaxis()->SetRangeUser(-1,1);
  gPad->SetLogz();
  c[1]->cd(3);
  hG->Draw("COLZ");
  hG->GetXaxis()->SetRangeUser(-300,300);
  hG->GetYaxis()->SetRangeUser(-1,1);
  gPad->SetLogz();
  c[1]->cd(4);
  hH->Draw("COLZ");
  hH->GetXaxis()->SetRangeUser(-300,300);
  hH->GetYaxis()->SetRangeUser(-300,300);
  gPad->SetLogz();

  c[2]->cd(1);
  hI->Draw("COLZ");
  hI->SetTitle("Lepton Vs Hadron Bunch Z Length;Hadron Z [mm];Lepton Z [mm]");
  gPad->SetLogz();
  c[2]->cd(2);
  hJ->Draw("COLZ");
  hJ->GetYaxis()->SetRangeUser(-1,1);
  gPad->SetLogz();
  c[2]->cd(3);
  hK->Draw("COLZ");
  hK->GetXaxis()->SetRangeUser(-100,100);
  hK->GetYaxis()->SetRangeUser(-1,1);
  gPad->SetLogz();

  c[3]->cd(1);
  hL->Draw("HIST");
  hL->SetTitle("Hadron Beam Angular Deviation: X;Rad");
  c[3]->cd(2);
  hM->Draw("HIST");
  hM->SetTitle("Hadron Beam Angular Deviation: Y;Rad");
  hM->GetXaxis()->SetRangeUser(-0.002,0.002);
  c[3]->cd(3);
  hN->Draw("HIST");
  hN->SetTitle("Lepton Beam Angular Deviation: X;Rad");
  hN->GetXaxis()->SetRangeUser(-0.002,0.002);
  c[3]->cd(4);
  hO->Draw("HIST");
  hO->SetTitle("Lepton Beam Angular Deviation: Y;Rad");
  hO->GetXaxis()->SetRangeUser(-0.002,0.002);
}


void plotSingleConfigKinematics()
{
  c[0]=new TCanvas("c0","Final State Particle Kinematics",800,600);

  c[0]->Clear();
  c[0]->Divide(2,2);

  TH2D *hA=(TH2D *)wrk->Get("partPtVsEta");
  TH2D *hB=(TH2D *)wrk->Get("partPhiVsEta");
  TH2D *hC=(TH2D *)wrk->Get("partPhiVsEtaHi");

  c[0]->cd(1);
  hA->Draw("COLZ");
  hA->GetYaxis()->SetRangeUser(0,20);
  hA->SetXTitle("Eta"); hA->SetYTitle("p_{T}");
  gPad->SetLogz();
  c[0]->cd(2);
  hB->Draw("COLZ");
  hB->SetXTitle("Eta"); hB->SetYTitle("Phi");
  gPad->SetLogz();
  c[0]->cd(3);
  hC->Draw("COLZ");
  hC->SetXTitle("Eta"); hC->SetYTitle("Phi");
  gPad->SetLogz();
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


void plotPartKinematics()
{
  c[0]=new TCanvas("c0","Particle Phi Vs Eta",800,600);
  c[1]=new TCanvas("c1","Particle Phi and Eta Projections",800,600);
  c[2]=new TCanvas("c2","Particle Phi in Eta Slices",800,600);
  c[3]=new TCanvas("c3","Particle Status Vs Eta",800,600);
  c[4]=new TCanvas("c4","Particle Eta in Status Slices",800,600);
  c[5]=new TCanvas("c5","Particle Pt Vs Eta",800,600);
  //c[6]=new TCanvas("c6","Lepton Beam XYZ",800,600);
  //c[7]=new TCanvas("c7","Particle Kinematics",800,600);

  c[0]->Clear();
  c[0]->Divide(2,2);

  c[1]->Clear();
  c[1]->Divide(2,2);

  c[2]->Clear();
  c[2]->Divide(2,2);

  c[3]->Clear();
  c[3]->Divide(2,2);

  c[4]->Clear();
  c[4]->Divide(2,2);

  c[5]->Clear();
  c[5]->Divide(2,2);

  //c[6]->Clear();
  //c[6]->Divide(2,2);

  //c[7]->Clear();
  //c[7]->Divide(2,2);

  TH2D *hPhiVsEta[4];
  TH2D *hStatusVsEta[4];
  TH2D *hPtVsEta[2];

  // Phi Vs Eta
  hPhiVsEta[0]=(TH2D *)fa->Get("partPhiVsEta");
  hPhiVsEta[1]=(TH2D *)fa->Get("partPhiVsEtaHi");
  hPhiVsEta[2]=(TH2D *)ff->Get("partPhiVsEta");
  hPhiVsEta[3]=(TH2D *)ff->Get("partPhiVsEtaHi");

  // Status Vs Eta
  hStatusVsEta[0]=(TH2D *)fa->Get("partStatusVsEta");
  hStatusVsEta[1]=(TH2D *)fa->Get("partStatusVsEtaHi");
  hStatusVsEta[2]=(TH2D *)ff->Get("partStatusVsEta");
  hStatusVsEta[3]=(TH2D *)ff->Get("partStatusVsEtaHi");

  // Pt Vs Eta
  hPtVsEta[0]=(TH2D *)fa->Get("partPtVsEta");
  hPtVsEta[1]=(TH2D *)ff->Get("partPtVsEta");

  // Projections
  TH1D *pPhiFull[4];
  TH1D *pEtaFull[4];
  TH1D *pPhiEtaSlice[4][4];
  TH1D *pEtaStatusSlice[4][3];
  TH1D *pPtEtaSlice[2][2];
  for(int i=0; i<4; i++)
    {
      pPhiFull[i] = hPhiVsEta[i]->ProjectionY(Form("p0_%d",i),1,400);
      pEtaFull[i] = hPhiVsEta[i]->ProjectionX(Form("p1_%d",i),1,100);

      pPhiEtaSlice[i][0] = hPhiVsEta[i]->ProjectionY(Form("p2_%d",i),160,180);
      pPhiEtaSlice[i][1] = hPhiVsEta[i]->ProjectionY(Form("p3_%d",i),200,220);
      pPhiEtaSlice[i][2] = hPhiVsEta[i]->ProjectionY(Form("p4_%d",i),240,260);
      pPhiEtaSlice[i][3] = hPhiVsEta[i]->ProjectionY(Form("p5_%d",i),280,300);

      pEtaStatusSlice[i][0] = hStatusVsEta[i]->ProjectionX(Form("p6_%d",i),61,66);
      pEtaStatusSlice[i][1] = hStatusVsEta[i]->ProjectionX(Form("p7_%d",i),81,86);
      pEtaStatusSlice[i][2] = hStatusVsEta[i]->ProjectionX(Form("p8_%d",i),91,96);

      if(i<2)
	{
	  pPtEtaSlice[i][0] = hPtVsEta[i]->ProjectionY(Form("p9_%d",i),1,500);
	  pPtEtaSlice[i][1] = hPtVsEta[i]->ProjectionY(Form("p10_%d",i),1,270);
	}
    }


  // Titles
  char titles0[4][100];
  sprintf(titles0[0],"Default Particle Phi Vs Eta;Eta;Phi");
  sprintf(titles0[1],"Default Particle Phi Vs Eta (pT > 1 GeV);Eta;Phi");
  sprintf(titles0[2],"All Effects Particle Phi Vs Eta;Eta;Phi");
  sprintf(titles0[3],"All Effects Particle Phi Vs Eta (pT > 1 GeV);Eta;Phi");

  char titles1[4][100];
  sprintf(titles1[0],"Default Particle Status Vs Eta;Eta;Status");
  sprintf(titles1[1],"Default Particle Status Vs Eta (pT > 1 GeV);Eta;Status");
  sprintf(titles1[2],"All Effects Particle Status Vs Eta;Eta;Status");
  sprintf(titles1[3],"All Effects Particle Status Vs Eta (pT > 1 GeV);Eta;Status");

  // Lines
  TLine *l00a = new TLine(-2.0,-1.0*TMath::Pi(),-2.0,TMath::Pi());
  TLine *l00b = new TLine(-1.0,-1.0*TMath::Pi(),-1.0,TMath::Pi());
  TLine *l01a = new TLine(0.0,-1.0*TMath::Pi(),0.0,TMath::Pi());
  TLine *l01b = new TLine(1.0,-1.0*TMath::Pi(),1.0,TMath::Pi());
  TLine *l02a = new TLine(2.0,-1.0*TMath::Pi(),2.0,TMath::Pi());
  TLine *l02b = new TLine(3.0,-1.0*TMath::Pi(),3.0,TMath::Pi());
  TLine *l03a = new TLine(4.0,-1.0*TMath::Pi(),4.0,TMath::Pi());
  TLine *l03b = new TLine(5.0,-1.0*TMath::Pi(),5.0,TMath::Pi());

  for(int i=0; i<4; i++)
    {
      // Part Phi Vs Eta
      c[0]->cd(i+1);
      hPhiVsEta[i]->Draw("COLZ");
      hPhiVsEta[i]->SetTitle(Form("%s",titles0[i]));
      gPad->SetLogz();

      if(i<2)
	{
	  l00a->Draw("SAME");
	  l00b->Draw("SAME");
	  l00a->SetLineColor(kBlue);
	  l00b->SetLineColor(kBlue);
	  l01a->Draw("SAME");
	  l01b->Draw("SAME");
	  l01a->SetLineColor(kRed);
	  l01b->SetLineColor(kRed);
	  l02a->Draw("SAME");
	  l02b->Draw("SAME");
	  l02a->SetLineColor(kGreen+2);
	  l02b->SetLineColor(kGreen+2);
	  l03a->Draw("SAME");
	  l03b->Draw("SAME");
	  l03a->SetLineColor(kBlack);
	  l03b->SetLineColor(kBlack);
	}

      // Particle Status Vs Eta
      c[3]->cd(i+1);
      hStatusVsEta[i]->Draw("COLZ");
      hStatusVsEta[i]->SetTitle(Form("%s",titles1[i]));
      gPad->SetLogz();
    }

  // Eta Phi Projections
  c[1]->cd(1);
  pPhiFull[0]->Draw("HIST");
  pPhiFull[0]->SetLineColor(kBlue);
  pPhiFull[2]->Draw("HISTSAME");
  pPhiFull[2]->SetLineColor(kRed);
  pPhiFull[0]->SetTitle("Particle Phi: All Eta;Phi");
  pPhiFull[0]->GetYaxis()->SetRangeUser(0,600000);
  c[1]->cd(2);
  pEtaFull[0]->Draw("HIST");
  pEtaFull[0]->SetLineColor(kBlue);
  pEtaFull[2]->Draw("HISTSAME");
  pEtaFull[2]->SetLineColor(kRed);
  pEtaFull[0]->SetTitle("Particle Eta: All Phi;Eta");
  pEtaFull[0]->GetYaxis()->SetRangeUser(0,400000);
  c[1]->cd(3);
  pPhiFull[1]->Draw("HIST");
  pPhiFull[1]->SetLineColor(kBlue);
  pPhiFull[3]->Draw("HISTSAME");
  pPhiFull[3]->SetLineColor(kRed);
  pPhiFull[1]->SetTitle("Particle Phi (pT > 1 GeV): All Eta;Phi");
  pPhiFull[1]->GetYaxis()->SetRangeUser(0,250000);
  c[1]->cd(4);
  pEtaFull[1]->Draw("HIST");
  pEtaFull[1]->SetLineColor(kBlue);
  pEtaFull[3]->Draw("HISTSAME");
  pEtaFull[3]->SetLineColor(kRed);
  pEtaFull[1]->SetTitle("Particle Eta (pT > 1 GeV): All Phi;Eta");
  pEtaFull[1]->GetYaxis()->SetRangeUser(1,500000);
  gPad->SetLogy();
  
  TLegend *leg1;
  leg1 = new TLegend(0.15,0.85,0.5,0.70);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg1->AddEntry(pEtaFull[1],"Default","l");
  leg1->AddEntry(pEtaFull[3],"All Effects","l");
  leg1->Draw();
  
  // Particle Phi Eta Slices
  c[2]->cd(1);
  pPhiEtaSlice[2][0]->Draw("HIST");
  pPhiEtaSlice[2][0]->SetLineColor(kBlue);
  pPhiEtaSlice[2][1]->Draw("HISTSAME");
  pPhiEtaSlice[2][1]->SetLineColor(kRed);
  pPhiEtaSlice[2][2]->Draw("HISTSAME");
  pPhiEtaSlice[2][2]->SetLineColor(kGreen+2);
  pPhiEtaSlice[2][3]->Draw("HISTSAME");
  pPhiEtaSlice[2][3]->SetLineColor(kBlack);
  pPhiEtaSlice[2][0]->SetTitle("All Effects Particle Phi in Eta Slices;Phi");
  pPhiEtaSlice[2][0]->GetYaxis()->SetRangeUser(500,500000);
  gPad->SetLogy();

  TLegend *leg2a;
  leg2a = new TLegend(0.15,0.35,0.45,0.20);
  leg2a->SetFillColor(0);
  leg2a->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg2a->AddEntry(pPhiEtaSlice[2][0],"-2 < #eta < -1","l");
  leg2a->AddEntry(pPhiEtaSlice[2][1],"0 < #eta < 1","l");
  leg2a->Draw();

  TLegend *leg2b;
  leg2b = new TLegend(0.55,0.35,0.85,0.20);
  leg2b->SetFillColor(0);
  leg2b->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg2b->AddEntry(pPhiEtaSlice[2][2],"2 < #eta < 3","l");
  leg2b->AddEntry(pPhiEtaSlice[2][3],"4 < #eta < 5","l");
  leg2b->Draw();

  c[2]->cd(2);
  pPhiEtaSlice[3][0]->Draw("HIST");
  pPhiEtaSlice[3][0]->SetLineColor(kBlue);
  pPhiEtaSlice[3][1]->Draw("HISTSAME");
  pPhiEtaSlice[3][1]->SetLineColor(kRed);
  pPhiEtaSlice[3][2]->Draw("HISTSAME");
  pPhiEtaSlice[3][2]->SetLineColor(kGreen+2);
  pPhiEtaSlice[3][3]->Draw("HISTSAME");
  pPhiEtaSlice[3][3]->SetLineColor(kBlack);
  pPhiEtaSlice[3][0]->SetTitle("All Effects Particle Phi in Eta Slices (pT > 1 GeV);Phi");
  pPhiEtaSlice[3][0]->GetYaxis()->SetRangeUser(1,500000);
  gPad->SetLogy();

  // Particle Eta in Status Slices
  c[4]->cd(1);
  pEtaStatusSlice[0][0]->Draw("HIST");
  pEtaStatusSlice[0][0]->SetLineColor(kBlue);
  pEtaStatusSlice[0][1]->Draw("HISTSAME");
  pEtaStatusSlice[0][1]->SetLineColor(kRed);
  pEtaStatusSlice[0][2]->Draw("HISTSAME");
  pEtaStatusSlice[0][2]->SetLineColor(kGreen+2);
  pEtaStatusSlice[0][0]->SetTitle("Default Particle Eta: Status Breakdown;Eta");
  pEtaStatusSlice[0][0]->GetYaxis()->SetRangeUser(1,300000);
  gPad->SetLogy();
  c[4]->cd(2);
  pEtaStatusSlice[1][0]->Draw("HIST");
  pEtaStatusSlice[1][0]->SetLineColor(kBlue);
  pEtaStatusSlice[1][1]->Draw("HISTSAME");
  pEtaStatusSlice[1][1]->SetLineColor(kRed);
  pEtaStatusSlice[1][2]->Draw("HISTSAME");
  pEtaStatusSlice[1][2]->SetLineColor(kGreen+2);
  pEtaStatusSlice[1][0]->SetTitle("Default Particle Eta: Status Breakdown (pT > 1 GeV);Eta");
  pEtaStatusSlice[1][0]->GetYaxis()->SetRangeUser(1,30000);
  gPad->SetLogy();
  c[4]->cd(3);
  pEtaStatusSlice[2][0]->Draw("HIST");
  pEtaStatusSlice[2][0]->SetLineColor(kBlue);
  pEtaStatusSlice[2][1]->Draw("HISTSAME");
  pEtaStatusSlice[2][1]->SetLineColor(kRed);
  pEtaStatusSlice[2][2]->Draw("HISTSAME");
  pEtaStatusSlice[2][2]->SetLineColor(kGreen+2);
  pEtaStatusSlice[2][0]->SetTitle("All Effects Particle Eta: Status Breakdown;Eta");
  pEtaStatusSlice[2][0]->GetYaxis()->SetRangeUser(1,500000);
  gPad->SetLogy();
  c[4]->cd(4);
  pEtaStatusSlice[3][0]->Draw("HIST");
  pEtaStatusSlice[3][0]->SetLineColor(kBlue);
  pEtaStatusSlice[3][1]->Draw("HISTSAME");
  pEtaStatusSlice[3][1]->SetLineColor(kRed);
  pEtaStatusSlice[3][2]->Draw("HISTSAME");
  pEtaStatusSlice[3][2]->SetLineColor(kGreen+2);
  pEtaStatusSlice[3][0]->SetTitle("All Effects Particle Eta: Status Breakdown (pT > 1 GeV);Eta");
  pEtaStatusSlice[3][0]->GetYaxis()->SetRangeUser(1,300000);
  gPad->SetLogy();

  TLegend *leg4;
  leg4 = new TLegend(0.15,0.85,0.40,0.60);
  leg4->SetFillColor(0);
  leg4->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg4->AddEntry(pEtaStatusSlice[3][0],"Beam Remnant","l");
  leg4->AddEntry(pEtaStatusSlice[3][1],"Hadronization","l");
  leg4->AddEntry(pEtaStatusSlice[3][2],"Decay","l");
  leg4->Draw();

  // Particle Pt Vs Eta
  c[5]->cd(1);
  hPtVsEta[0]->Draw("COLZ");
  hPtVsEta[0]->SetTitle("Default Particle Pt Vs Eta;Eta;Pt");
  hPtVsEta[0]->GetYaxis()->SetRangeUser(0,20);
  gPad->SetLogz();
  c[5]->cd(2);
  hPtVsEta[1]->Draw("COLZ");
  hPtVsEta[1]->SetTitle("All Effects Particle Pt Vs Eta;Eta;Pt");
  hPtVsEta[1]->GetYaxis()->SetRangeUser(0,20);
  gPad->SetLogz();
  c[5]->cd(3);
  pPtEtaSlice[0][0]->SetLineColor(kBlue);
  pPtEtaSlice[0][1]->SetLineColor(kRed);
  pPtEtaSlice[1][0]->SetLineColor(kGreen+2);
  pPtEtaSlice[1][1]->SetLineColor(kBlack);
  pPtEtaSlice[0][0]->SetTitle("Particle Pt Projections;Pt");
  pPtEtaSlice[0][0]->GetXaxis()->SetRangeUser(0,30);
  pPtEtaSlice[0][0]->DrawCopy("HIST");
  pPtEtaSlice[0][1]->DrawCopy("HISTSAME");
  pPtEtaSlice[1][0]->DrawCopy("HISTSAME");
  pPtEtaSlice[1][1]->DrawCopy("HISTSAME");
  gPad->SetLogy();

  TLegend *leg5;
  leg5 = new TLegend(0.25,0.85,0.70,0.55);
  leg5->SetFillColor(0);
  leg5->SetBorderSize(0);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg5->AddEntry(pPtEtaSlice[0][0],"All #eta: Default","l");
  leg5->AddEntry(pPtEtaSlice[0][1],"#eta < 3.5: Default","l");
  leg5->AddEntry(pPtEtaSlice[1][0],"All #eta: All Effects","l");
  leg5->AddEntry(pPtEtaSlice[1][1],"#eta < 3.5: All Effects","l");
  leg5->Draw();

  c[5]->cd(4);
  pPtEtaSlice[0][0]->Draw("HIST");
  pPtEtaSlice[0][1]->Draw("HISTSAME");
  pPtEtaSlice[1][0]->Draw("HISTSAME");
  pPtEtaSlice[1][1]->Draw("HISTSAME");
  pPtEtaSlice[0][0]->SetTitle("Particle Pt Projections;Pt");
  pPtEtaSlice[0][0]->GetXaxis()->SetRangeUser(0,5);
  gPad->SetLogy();
}
