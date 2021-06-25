TFile *fa=0;
TFile *fb=0;
TFile *fc=0;
TFile *fd=0;

TFile *fw=0;
TFile *fx=0;
TFile *fy=0;
TFile *fz=0;

TCanvas *c[30];

//TString outFileName;

int plotConfigComp()
{

  //gStyle->SetPalette(1,0);
  gStyle->SetOptStat(0);
  //gStyle->SetOptStat(1111);
  //gStyle->SetOptFit(1111);

  TString rootHistFnameA;
  TString rootHistFnameB;
  TString rootHistFnameC;
  TString rootHistFnameD;

  TString rootHistFnameW;
  TString rootHistFnameX;
  TString rootHistFnameY;
  TString rootHistFnameZ;

  rootHistFnameA="rootFiles/headonTestJin/test_crossDivNrgCrab_25mRad_18x275_v2.hist.root";
  rootHistFnameB="rootFiles/headonTestJin/test_crossDivNrgCrab_35mRad_18x275_v2.hist.root";
  rootHistFnameC="rootFiles/headonTestJin/test_crossDivNrgCrab_25mRad_5x41_v2.hist.root";
  rootHistFnameD="rootFiles/headonTestJin/test_crossDivNrgCrab_35mRad_5x41_v2.hist.root";

  rootHistFnameW="rootFiles/headonTestJin/test_headon_25mRad_18x275_v1.hist.root";
  rootHistFnameX="rootFiles/headonTestJin/test_headon_35mRad_18x275_v1.hist.root";
  rootHistFnameY="rootFiles/headonTestJin/test_headon_25mRad_5x41_v1.hist.root";
  rootHistFnameZ="rootFiles/headonTestJin/test_headon_35mRad_5x41_v1.hist.root";
  
  fa=new TFile(rootHistFnameA);
  assert(fa->IsOpen());

  fb=new TFile(rootHistFnameB);
  assert(fb->IsOpen());

  fc=new TFile(rootHistFnameC);
  assert(fc->IsOpen());

  fd=new TFile(rootHistFnameD);
  assert(fd->IsOpen());

  fw=new TFile(rootHistFnameW);
  assert(fw->IsOpen());

  fx=new TFile(rootHistFnameX);
  assert(fx->IsOpen());

  fy=new TFile(rootHistFnameY);
  assert(fy->IsOpen());

  fz=new TFile(rootHistFnameZ);
  assert(fz->IsOpen());

  return 0;
}


void plotBeamMom()
{
  c[0]=new TCanvas("c0","Energy Shifts",800,600);
  c[1]=new TCanvas("c1","Hadron Beam Py Vs Px",800,600);
  c[2]=new TCanvas("c2","Lepton Beam Py Vs Px",800,600);
  c[3]=new TCanvas("c3","Hadron Beam Px Vs Z Position",800,600);
  c[4]=new TCanvas("c4","Lepton Beam Px Vs Z Position",800,600);

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

  TH1D *hCM[4];
  TH2D *hXY1[4];
  TH2D *hXY2[4];
  TH2D *hXZProd1[4];
  TH2D *hXZProd2[4];

  hCM[0]=(TH1D *)fa->Get("eCM"); // 18x275 25
  hCM[1]=(TH1D *)fb->Get("eCM"); // 18x275 35
  hCM[2]=(TH1D *)fc->Get("eCM"); // 5x41 25
  hCM[3]=(TH1D *)fd->Get("eCM"); // 5x41 35

  hXY1[0]=(TH2D *)fa->Get("pXY1"); // 18x275 25
  hXY1[1]=(TH2D *)fb->Get("pXY1"); // 18x275 35
  hXY1[2]=(TH2D *)fc->Get("pXY1"); // 5x41 25
  hXY1[3]=(TH2D *)fd->Get("pXY1"); // 5x41 35

  hXY2[0]=(TH2D *)fa->Get("pXY2"); // 18x275 25
  hXY2[1]=(TH2D *)fb->Get("pXY2"); // 18x275 35
  hXY2[2]=(TH2D *)fc->Get("pXY2"); // 5x41 25
  hXY2[3]=(TH2D *)fd->Get("pXY2"); // 5x41 35

  hXZProd1[0]=(TH2D *)fa->Get("pXZProd1"); // 18x275 25
  hXZProd1[1]=(TH2D *)fb->Get("pXZProd1"); // 18x275 35
  hXZProd1[2]=(TH2D *)fc->Get("pXZProd1"); // 5x41 25
  hXZProd1[3]=(TH2D *)fd->Get("pXZProd1"); // 5x41 35

  hXZProd2[0]=(TH2D *)fa->Get("pXZProd2"); // 18x275 25
  hXZProd2[1]=(TH2D *)fb->Get("pXZProd2"); // 18x275 35
  hXZProd2[2]=(TH2D *)fc->Get("pXZProd2"); // 5x41 25
  hXZProd2[3]=(TH2D *)fd->Get("pXZProd2"); // 5x41 35

  c[0]->cd(1);
  hCM[0]->Draw("HIST");
  hCM[0]->Rebin(4);
  hCM[0]->SetTitle("Modified-Nominal CM Energy: 18x275 25mRad;Energy Diff [GeV]");
  hCM[0]->GetXaxis()->SetRangeUser(-0.5,0.5);
  c[0]->cd(2);
  hCM[1]->Draw("HIST");
  hCM[1]->Rebin(4);
  hCM[1]->SetTitle("Modified-Nominal CM Energy: 18x275 35mRad;Energy Diff [GeV]");
  hCM[1]->GetXaxis()->SetRangeUser(-0.5,0.5);
  c[0]->cd(3);
  hCM[2]->Draw("HIST");
  hCM[2]->SetTitle("Modified-Nominal CM Energy: 5x41 25mRad;Energy Diff [GeV]");
  hCM[2]->GetXaxis()->SetRangeUser(-0.5,0.5);
  c[0]->cd(4);
  hCM[3]->Draw("HIST");
  hCM[3]->SetTitle("Modified-Nominal CM Energy: 5x41 35mRad;Energy Diff [GeV]");
  hCM[3]->GetXaxis()->SetRangeUser(-0.5,0.5);

  c[1]->cd(1);
  hXY1[0]->Draw("COLZ");
  hXY1[0]->SetTitle("Hadron Beam Y vs X Momentum: 18x275 25mRad;Px;Py");
  hXY1[0]->GetXaxis()->SetRangeUser(6,8);
  hXY1[0]->GetYaxis()->SetRangeUser(-0.5,0.5);
  gPad->SetLogz();
  c[1]->cd(2);
  hXY1[1]->Draw("COLZ");
  hXY1[1]->SetTitle("Hadron Beam Y vs X Momentum: 18x275 35mRad;Px;Py");
  hXY1[1]->GetXaxis()->SetRangeUser(8,10);
  hXY1[1]->GetYaxis()->SetRangeUser(-0.5,0.5);
  gPad->SetLogz();
  c[1]->cd(3);
  hXY1[2]->Draw("COLZ");
  hXY1[2]->SetTitle("Hadron Beam Y vs X Momentum: 5x41 25mRad;Px;Py");
  hXY1[2]->GetXaxis()->SetRangeUser(0,2);
  hXY1[2]->GetYaxis()->SetRangeUser(-0.5,0.5);
  gPad->SetLogz();
  c[1]->cd(4);
  hXY1[3]->Draw("COLZ");
  hXY1[3]->SetTitle("Hadron Beam Y vs X Momentum: 5x41 35mRad;Px;Py");
  hXY1[3]->GetXaxis()->SetRangeUser(0,2);
  hXY1[3]->GetYaxis()->SetRangeUser(-0.5,0.5);
  gPad->SetLogz();

  c[2]->cd(1);
  hXY2[0]->Draw("COLZ");
  hXY2[0]->SetTitle("Lepton Beam Y vs X Momentum: 18x275 25mRad;Px;Py");
  hXY2[0]->GetXaxis()->SetRangeUser(-0.1,0.1);
  hXY2[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();
  c[2]->cd(2);
  hXY2[1]->Draw("COLZ");
  hXY2[1]->SetTitle("Lepton Beam Y vs X Momentum: 18x275 35mRad;Px;Py");
  hXY2[1]->GetXaxis()->SetRangeUser(-0.1,0.1);
  hXY2[1]->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();
  c[2]->cd(3);
  hXY2[2]->Draw("COLZ");
  hXY2[2]->SetTitle("Lepton Beam Y vs X Momentum: 5x41 25mRad;Px;Py");
  hXY2[2]->GetXaxis()->SetRangeUser(-0.1,0.1);
  hXY2[2]->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();
  c[2]->cd(4);
  hXY2[3]->Draw("COLZ");
  hXY2[3]->SetTitle("Lepton Beam Y vs X Momentum: 5x41 35mRad;Px;Py");
  hXY2[3]->GetXaxis()->SetRangeUser(-0.1,0.1);
  hXY2[3]->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();

  c[3]->cd(1);
  hXZProd1[0]->Draw("COLZ");
  hXZProd1[0]->SetTitle("Hadron Beam X Momentum vs Z Vertex Position: 18x275 25mRad;Z Vert [mm];Px");
  hXZProd1[0]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZProd1[0]->GetYaxis()->SetRangeUser(6.,8.);
  gPad->SetLogz();
  c[3]->cd(2);
  hXZProd1[1]->Draw("COLZ");
  hXZProd1[1]->SetTitle("Hadron Beam X Momentum vs Z Vertex Position: 18x275 35mRad;Z Vert [mm];Px");
  hXZProd1[1]->GetYaxis()->SetRangeUser(8.,10.);
  hXZProd1[1]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();
  c[3]->cd(3);
  hXZProd1[2]->Draw("COLZ");
  hXZProd1[2]->SetTitle("Hadron Beam X Momentum vs Z Vertex Position: 5x41 25mRad;Z Vert [mm];Px");
  hXZProd1[2]->GetYaxis()->SetRangeUser(0,2);
  hXZProd1[2]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();
  c[3]->cd(4);
  hXZProd1[3]->Draw("COLZ");
  hXZProd1[3]->SetTitle("Hadron Beam X Momentum vs Z Vertex Position: 5x41 35mRad;Z Vert [mm];Px");
  hXZProd1[3]->GetYaxis()->SetRangeUser(0,2);
  hXZProd1[3]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();

  c[4]->cd(1);
  hXZProd2[0]->Draw("COLZ");
  hXZProd2[0]->SetTitle("Lepton Beam X Momentum vs Z Vertex Position: 18x275 25mRad;Z Vert [mm];Px");
  hXZProd2[0]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZProd2[0]->GetYaxis()->SetRangeUser(-0.1,0.1);
  gPad->SetLogz();
  c[4]->cd(2);
  hXZProd2[1]->Draw("COLZ");
  hXZProd2[1]->SetTitle("Lepton Beam X Momentum vs Z Vertex Position: 18x275 35mRad;Z Vert [mm];Px");
  hXZProd2[1]->GetYaxis()->SetRangeUser(-0.1,0.1);
  hXZProd2[1]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();
  c[4]->cd(3);
  hXZProd2[2]->Draw("COLZ");
  hXZProd2[2]->SetTitle("Lepton Beam X Momentum vs Z Vertex Position: 5x41 25mRad;Z Vert [mm];Px");
  hXZProd2[2]->GetYaxis()->SetRangeUser(-0.1,0.1);
  hXZProd2[2]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();
  c[4]->cd(4);
  hXZProd2[3]->Draw("COLZ");
  hXZProd2[3]->SetTitle("Lepton Beam X Momentum vs Z Vertex Position: 5x41 35mRad;Z Vert [mm];Px");
  hXZProd2[3]->GetYaxis()->SetRangeUser(-0.1,0.1);
  hXZProd2[3]->GetXaxis()->SetRangeUser(-200.,200.);
  gPad->SetLogz();
}


void plotVertex()
{
  c[0]=new TCanvas("c0","Bunch Z Size",800,600);
  c[1]=new TCanvas("c1","Vertex Distributions",800,600);
  c[2]=new TCanvas("c2","Vertex Y Vs X",800,600);
  c[3]=new TCanvas("c3","Vertex X Vs Z",800,600);
  c[4]=new TCanvas("c4","Vertex X Vs T",800,600);
  c[5]=new TCanvas("c5","Vertex X Vs T+Z",800,600);

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

  TH2D *hBZ[4];
  TH1D *hVert[4][4];
  TH2D *hXY[4];
  TH2D *hXZ[4];
  TH2D *hXT[4];
  TH2D *hXTZ[4];

  hBZ[0]=(TH2D *)fa->Get("lepVsHadPartZ"); // 18x275 25
  hBZ[1]=(TH2D *)fb->Get("lepVsHadPartZ"); // 18x275 35
  hBZ[2]=(TH2D *)fc->Get("lepVsHadPartZ"); // 5x41 25
  hBZ[3]=(TH2D *)fd->Get("lepVsHadPartZ"); // 5x41 35

  hVert[0][0]=(TH1D *)fa->Get("vtxX");
  hVert[0][1]=(TH1D *)fb->Get("vtxX");
  hVert[0][2]=(TH1D *)fc->Get("vtxX");
  hVert[0][3]=(TH1D *)fd->Get("vtxX");
  hVert[1][0]=(TH1D *)fa->Get("vtxY");
  hVert[1][1]=(TH1D *)fb->Get("vtxY");
  hVert[1][2]=(TH1D *)fc->Get("vtxY");
  hVert[1][3]=(TH1D *)fd->Get("vtxY");
  hVert[2][0]=(TH1D *)fa->Get("vtxZ");
  hVert[2][1]=(TH1D *)fb->Get("vtxZ");
  hVert[2][2]=(TH1D *)fc->Get("vtxZ");
  hVert[2][3]=(TH1D *)fd->Get("vtxZ");
  hVert[3][0]=(TH1D *)fa->Get("vtxT");
  hVert[3][1]=(TH1D *)fb->Get("vtxT");
  hVert[3][2]=(TH1D *)fc->Get("vtxT");
  hVert[3][3]=(TH1D *)fd->Get("vtxT");

  hXY[0]=(TH2D *)fa->Get("vtxYvsX"); // 18x275 25
  hXY[1]=(TH2D *)fb->Get("vtxYvsX"); // 18x275 35
  hXY[2]=(TH2D *)fc->Get("vtxYvsX"); // 5x41 25
  hXY[3]=(TH2D *)fd->Get("vtxYvsX"); // 5x41 35

  hXZ[0]=(TH2D *)fa->Get("vtxXvsZ"); // 18x275 25
  hXZ[1]=(TH2D *)fb->Get("vtxXvsZ"); // 18x275 35
  hXZ[2]=(TH2D *)fc->Get("vtxXvsZ"); // 5x41 25
  hXZ[3]=(TH2D *)fd->Get("vtxXvsZ"); // 5x41 35

  hXT[0]=(TH2D *)fa->Get("vtxXvsT"); // 18x275 25
  hXT[1]=(TH2D *)fb->Get("vtxXvsT"); // 18x275 35
  hXT[2]=(TH2D *)fc->Get("vtxXvsT"); // 5x41 25
  hXT[3]=(TH2D *)fd->Get("vtxXvsT"); // 5x41 35

  hXTZ[0]=(TH2D *)fa->Get("vtxXvsTZSum"); // 18x275 25
  hXTZ[1]=(TH2D *)fb->Get("vtxXvsTZSum"); // 18x275 35
  hXTZ[2]=(TH2D *)fc->Get("vtxXvsTZSum"); // 5x41 25
  hXTZ[3]=(TH2D *)fd->Get("vtxXvsTZSum"); // 5x41 35

  c[0]->cd(1);
  hBZ[0]->Draw("COLZ");
  hBZ[0]->SetTitle("Lepton Vs Hadron Bunch Z Extent: 18x275 25mRad;Hadron Z Extent;Lepton Z Extent");
  gPad->SetLogz();
  c[0]->cd(2);
  hBZ[1]->Draw("COLZ");
  hBZ[1]->SetTitle("Lepton Vs Hadron Bunch Z Extent: 18x275 35mRad;Hadron Z Extent;Lepton Z Extent");
  gPad->SetLogz();
  c[0]->cd(3);
  hBZ[2]->Draw("COLZ");
  hBZ[2]->SetTitle("Lepton Vs Hadron Bunch Z Extent: 5x41 25mRad;Hadron Z Extent;Lepton Z Extent");
  gPad->SetLogz();
  c[0]->cd(4);
  hBZ[3]->Draw("COLZ");
  hBZ[3]->SetTitle("Lepton Vs Hadron Bunch Z Extent: 5x41 35mRad;Hadron Z Extent;Lepton Z Extent");
  gPad->SetLogz();

  int colors[4]={600,632,418,800};
  for(int i=0; i<4; i++)
    {
      c[1]->cd(i+1);
      hVert[i][0]->Draw("HIST");
      hVert[i][0]->SetLineColor(colors[0]);
      if(i==0) hVert[i][0]->GetXaxis()->SetRangeUser(-1,1);
      if(i==1) hVert[i][0]->GetXaxis()->SetRangeUser(-0.2,0.2);
      if(i==2) hVert[i][0]->GetXaxis()->SetRangeUser(-100,100);
      if(i==3) hVert[i][0]->GetXaxis()->SetRangeUser(-100,100);

      if(i==1)
	{
	  TLegend *leg1;
	  leg1 = new TLegend(0.15,0.85,0.47,0.60);
	  leg1->SetFillColor(0);
	  leg1->SetBorderSize(1);
	  //leg22->SetTextFont(63);
	  //leg22->SetTextSize(30);
	  leg1->AddEntry(hVert[0][0],"18x275; 25 mRad","l");
	  leg1->AddEntry(hVert[0][1],"18x275; 35 mRad","l");
	  leg1->AddEntry(hVert[0][2],"5x41; 25 mRad","l");
	  leg1->AddEntry(hVert[0][3],"5x41; 35 mRad","l");
	  leg1->Draw();
	}
      
      for(int j=1; j<4; j++)
	{
	  hVert[i][j]->Draw("HISTSAME");
	  hVert[i][j]->SetLineColor(colors[j]);
	}
    }

  c[2]->cd(1);
  hXY[0]->Draw("COLZ");
  hXY[0]->SetTitle("Vertex Y Vs X Position: 18x275 25mRad;X [mm];Y [mm]");
  hXY[0]->GetXaxis()->SetRangeUser(-1.,1.);
  hXY[0]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[2]->cd(2);
  hXY[1]->Draw("COLZ");
  hXY[1]->SetTitle("Vertex Y Vs X Position: 18x275 35mRad;X [mm];Y [mm]");
  hXY[1]->GetXaxis()->SetRangeUser(-1.,1.);
  hXY[1]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[2]->cd(3);
  hXY[2]->Draw("COLZ");
  hXY[2]->SetTitle("Vertex Y Vs X Position: 5x41 25mRad;X [mm];Y [mm]");
  hXY[2]->GetXaxis()->SetRangeUser(-1.,1.);
  hXY[2]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[2]->cd(4);
  hXY[3]->Draw("COLZ");
  hXY[3]->SetTitle("Vertex Y Vs X Position: 5x41 35mRad;X [mm];Y [mm]");
  hXY[3]->GetXaxis()->SetRangeUser(-1.,1.);
  hXY[3]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();

  c[3]->cd(1);
  hXZ[0]->Draw("COLZ");
  hXZ[0]->SetTitle("Vertex X Vs Z Position: 18x275 25mRad;Z [mm];X [mm]");
  hXZ[0]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZ[0]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[3]->cd(2);
  hXZ[1]->Draw("COLZ");
  hXZ[1]->SetTitle("Vertex X Vs Z Position: 18x275 35mRad;Z [mm];X [mm]");
  hXZ[1]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZ[1]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[3]->cd(3);
  hXZ[2]->Draw("COLZ");
  hXZ[2]->SetTitle("Vertex X Vs Z Position: 5x41 25mRad;Z [mm];X [mm]");
  hXZ[2]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZ[2]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[3]->cd(4);
  hXZ[3]->Draw("COLZ");
  hXZ[3]->SetTitle("Vertex X Vs Z Position: 5x41 35mRad;Z [mm];X [mm]");
  hXZ[3]->GetXaxis()->SetRangeUser(-200.,200.);
  hXZ[3]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();

  c[4]->cd(1);
  hXT[0]->Draw("COLZ");
  hXT[0]->SetTitle("Vertex X Vs T Position: 18x275 25mRad;T [mm];X [mm]");
  hXT[0]->GetXaxis()->SetRangeUser(-200.,200.);
  hXT[0]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[4]->cd(2);
  hXT[1]->Draw("COLZ");
  hXT[1]->SetTitle("Vertex X Vs T Position: 18x275 35mRad;T [mm];X [mm]");
  hXT[1]->GetXaxis()->SetRangeUser(-200.,200.);
  hXT[1]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[4]->cd(3);
  hXT[2]->Draw("COLZ");
  hXT[2]->SetTitle("Vertex X Vs T Position: 5x41 25mRad;T [mm];X [mm]");
  hXT[2]->GetXaxis()->SetRangeUser(-200.,200.);
  hXT[2]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[4]->cd(4);
  hXT[3]->Draw("COLZ");
  hXT[3]->SetTitle("Vertex X Vs T Position: 5x41 35mRad;T [mm];X [mm]");
  hXT[3]->GetXaxis()->SetRangeUser(-200.,200.);
  hXT[3]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();

  c[5]->cd(1);
  hXTZ[0]->Draw("COLZ");
  hXTZ[0]->SetTitle("Vertex X Vs T+Z: 18x275 25mRad;T+Z [mm];X [mm]");
  hXTZ[0]->GetXaxis()->SetRangeUser(-200.,200.);
  hXTZ[0]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[5]->cd(2);
  hXTZ[1]->Draw("COLZ");
  hXTZ[1]->SetTitle("Vertex X Vs T+Z: 18x275 35mRad;T+Z [mm];X [mm]");
  hXTZ[1]->GetXaxis()->SetRangeUser(-200.,200.);
  hXTZ[1]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[5]->cd(3);
  hXTZ[2]->Draw("COLZ");
  hXTZ[2]->SetTitle("Vertex X Vs T+Z: 5x41 25mRad;T+Z [mm];X [mm]");
  hXTZ[2]->GetXaxis()->SetRangeUser(-200.,200.);
  hXTZ[2]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
  c[5]->cd(4);
  hXTZ[3]->Draw("COLZ");
  hXTZ[3]->SetTitle("Vertex X Vs T+Z: 5x41 35mRad;T+Z [mm];X [mm]");
  hXTZ[3]->GetXaxis()->SetRangeUser(-200.,200.);
  hXTZ[3]->GetYaxis()->SetRangeUser(-1.,1.);
  gPad->SetLogz();
}


void plotFinalState()
{
  c[0]=new TCanvas("c0","Final State Phi Vs Eta",800,600);
  c[1]=new TCanvas("c1","Final State Eta",800,600);
  c[2]=new TCanvas("c2","Final State Phi",800,600);
  c[3]=new TCanvas("c3","Final State Pt Vs Eta",800,600);
  c[4]=new TCanvas("c4","Final State Phi Vs Eta: Default",800,600);
  c[5]=new TCanvas("c5","Final State Pt Vs Eta: Default",800,600);

  c[0]->Clear();
  c[0]->Divide(2,2);

  c[1]->Clear();
  c[1]->Divide(1,1);

  c[2]->Clear();
  c[2]->Divide(1,1);

  c[3]->Clear();
  c[3]->Divide(2,2);

  c[4]->Clear();
  c[4]->Divide(2,2);

  c[5]->Clear();
  c[5]->Divide(2,2);

  TH2D *hPE[6];
  TH1D *hEta[6];
  TH1D *hPhi[6];
  TH2D *hPtE[6];

  hPE[0]=(TH2D *)fa->Get("partPhiVsEta"); // 18x275 25
  hPE[1]=(TH2D *)fb->Get("partPhiVsEta"); // 18x275 35
  hPE[2]=(TH2D *)fc->Get("partPhiVsEta"); // 5x41 25
  hPE[3]=(TH2D *)fd->Get("partPhiVsEta"); // 5x41 35
  hPE[4]=(TH2D *)fw->Get("partPhiVsEta");
  hPE[5]=(TH2D *)fy->Get("partPhiVsEta");

  hEta[0]=(TH1D *)fa->Get("partEta");
  hEta[1]=(TH1D *)fb->Get("partEta");
  hEta[2]=(TH1D *)fc->Get("partEta");
  hEta[3]=(TH1D *)fd->Get("partEta");
  hEta[4]=(TH1D *)fw->Get("partEta");
  hEta[5]=(TH1D *)fy->Get("partEta");

  hPhi[0]=(TH1D *)fa->Get("partPhi");
  hPhi[1]=(TH1D *)fb->Get("partPhi");
  hPhi[2]=(TH1D *)fc->Get("partPhi");
  hPhi[3]=(TH1D *)fd->Get("partPhi");
  hPhi[4]=(TH1D *)fw->Get("partPhi");
  hPhi[5]=(TH1D *)fy->Get("partPhi");

  hPtE[0]=(TH2D *)fa->Get("partPtVsEta"); // 18x275 25
  hPtE[1]=(TH2D *)fb->Get("partPtVsEta"); // 18x275 35
  hPtE[2]=(TH2D *)fc->Get("partPtVsEta"); // 5x41 25
  hPtE[3]=(TH2D *)fd->Get("partPtVsEta"); // 5x41 35
  hPtE[4]=(TH2D *)fw->Get("partPtVsEta");
  hPtE[5]=(TH2D *)fy->Get("partPtVsEta");

  c[0]->cd(1);
  hPE[0]->Draw("COLZ");
  hPE[0]->SetTitle("Final State Particle Phi Vs Eta: 18x275 25mRad;Eta;Phi");
  gPad->SetLogz();
  c[0]->cd(2);
  hPE[1]->Draw("COLZ");
  hPE[1]->SetTitle("Final State Particle Phi Vs Eta: 18x275 35mRad;Eta;Phi");
  gPad->SetLogz();
  c[0]->cd(3);
  hPE[2]->Draw("COLZ");
  hPE[2]->SetTitle("Final State Particle Phi Vs Eta: 5x41 25mRad;Eta;Phi");
  gPad->SetLogz();
  c[0]->cd(4);
  hPE[3]->Draw("COLZ");
  hPE[3]->SetTitle("Final State Particle Phi Vs Eta: 5x41 35mRad;Eta;Phi");
  gPad->SetLogz();

  c[1]->cd(1);
  hEta[0]->Draw("HIST");
  hEta[0]->SetLineColor(kBlue);
  hEta[1]->Draw("HISTSAME");
  hEta[1]->SetLineColor(kRed);
  hEta[2]->Draw("HISTSAME");
  hEta[2]->SetLineColor(kGreen+2);
  hEta[3]->Draw("HISTSAME");
  hEta[3]->SetLineColor(kOrange);
  hEta[4]->Draw("HISTSAME");
  hEta[4]->SetLineColor(kGray);
  hEta[5]->Draw("HISTSAME");
  hEta[5]->SetLineColor(kGray+2);
  hEta[0]->SetTitle("Final State Particle Eta;Eta");
  hEta[0]->GetXaxis()->SetRangeUser(-4.,8.);
  hEta[0]->GetYaxis()->SetRangeUser(0.,700000.);

  TLegend *leg1;
  leg1 = new TLegend(0.15,0.85,0.50,0.60);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg1->AddEntry(hEta[0],"18x275; 25 mRad","l");
  leg1->AddEntry(hEta[1],"18x275; 35 mRad","l");
  leg1->AddEntry(hEta[2],"5x41; 25 mRad","l");
  leg1->AddEntry(hEta[3],"5x41; 35 mRad","l");
  leg1->AddEntry(hEta[4],"18x275; Nominal","l");
  leg1->AddEntry(hEta[5],"5x41; Nominal","l");
  leg1->Draw();

  c[2]->cd(1);
  hPhi[0]->Draw("HIST");
  hPhi[0]->SetLineColor(kBlue);
  hPhi[1]->Draw("HISTSAME");
  hPhi[1]->SetLineColor(kRed);
  hPhi[2]->Draw("HISTSAME");
  hPhi[2]->SetLineColor(kGreen+2);
  hPhi[3]->Draw("HISTSAME");
  hPhi[3]->SetLineColor(kOrange);
  hPhi[4]->Draw("HISTSAME");
  hPhi[4]->SetLineColor(kGray);
  hPhi[5]->Draw("HISTSAME");
  hPhi[5]->SetLineColor(kGray+2);
  hPhi[0]->SetTitle("Final State Particle Phi;Phi");
  //hPhi[0]->GetXaxis()->SetRangeUser(-4.,8.);
  hPhi[0]->GetYaxis()->SetRangeUser(0.,1000000.);

  TLegend *leg2;
  leg2 = new TLegend(0.15,0.85,0.45,0.60);
  leg2->SetFillColor(0);
  leg2->SetBorderSize(1);
  //leg22->SetTextFont(63);
  //leg22->SetTextSize(30);
  leg2->AddEntry(hPhi[0],"18x275; 25 mRad","l");
  leg2->AddEntry(hPhi[1],"18x275; 35 mRad","l");
  leg2->AddEntry(hPhi[2],"5x41; 25 mRad","l");
  leg2->AddEntry(hPhi[3],"5x41; 35 mRad","l");
  leg2->AddEntry(hPhi[4],"18x275; Nominal","l");
  leg2->AddEntry(hPhi[5],"5x41; Nominal","l");
  leg2->Draw();

  c[3]->cd(1);
  hPtE[0]->Draw("COLZ");
  hPtE[0]->SetTitle("Final State Particle Pt Vs Eta: 18x275 25mRad;Eta;Pt");
  gPad->SetLogz();
  c[3]->cd(2);
  hPtE[1]->Draw("COLZ");
  hPtE[1]->SetTitle("Final State Particle Pt Vs Eta: 18x275 35mRad;Eta;Pt");
  gPad->SetLogz();
  c[3]->cd(3);
  hPtE[2]->Draw("COLZ");
  hPtE[2]->SetTitle("Final State Particle Pt Vs Eta: 5x41 25mRad;Eta;Pt");
  hPtE[2]->GetYaxis()->SetRangeUser(0,20);
  gPad->SetLogz();
  c[3]->cd(4);
  hPtE[3]->Draw("COLZ");
  hPtE[3]->SetTitle("Final State Particle Pt Vs Eta: 5x41 35mRad;Eta;Pt");
  hPtE[3]->GetYaxis()->SetRangeUser(0,20);
  gPad->SetLogz();

  c[4]->cd(1);
  hPE[4]->Draw("COLZ");
  hPE[4]->SetTitle("Final State Particle Phi Vs Eta: 18x275 Nominal;Eta;Phi");
  gPad->SetLogz();
  c[4]->cd(3);
  hPE[5]->Draw("COLZ");
  hPE[5]->SetTitle("Final State Particle Phi Vs Eta: 5x41 Nominal;Eta;Phi");
  gPad->SetLogz();

  c[5]->cd(1);
  hPtE[4]->Draw("COLZ");
  hPtE[4]->SetTitle("Final State Particle Pt Vs Eta: 18x275 Nominal;Eta;Pt");
  gPad->SetLogz();
  c[5]->cd(3);
  hPtE[5]->Draw("COLZ");
  hPtE[5]->SetTitle("Final State Particle Pt Vs Eta: 18x275 Nominal;Eta;Pt");
  hPtE[5]->GetYaxis()->SetRangeUser(0,20);
  gPad->SetLogz();

}
