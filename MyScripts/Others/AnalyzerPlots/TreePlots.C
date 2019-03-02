void Save(TString path, TCanvas *c, TH1D* h, TString x, Double_t labelSize, Double_t titleSize) {

  h->Draw();
  h->GetXaxis()->SetTitle(x);
  h->SetFillColor(38);
  h->SetTitleSize(titleSize, "t");
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");
}

void Save(TString path, TCanvas *c, TH2D* h, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  h->Draw("colz");
  h->GetXaxis()->SetTitle(x);
  h->GetYaxis()->SetTitle(y);
  h->SetTitleSize(titleSize, "t");
  h->GetXaxis()->SetTitleSize(labelSize);
  h->GetYaxis()->SetTitleSize(labelSize);
  h->GetXaxis()->SetLabelSize(labelSize);
  h->GetYaxis()->SetLabelSize(labelSize);
  h->GetXaxis()->SetTitleOffset(1.4);
  h->GetYaxis()->SetTitleOffset(1.4);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");
}

void TreePlots(TString dir, TString histo1) {

  // dir = output dir, histo1 = histo to do cosmetics on

  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  TString path = "";

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("2016"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/2016/";
    else if (histo1.Contains("2017"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/2017/";
    else if (histo1.Contains("Capped"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/K3piCapped/";
    else if (!histo1.Contains("Capped") && histo1.Contains("K3pi"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/K3pi/";

  //path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/2/HeavyNeutrino/";
  }
  
  TFile *f = TFile::Open(histo1);

  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    return;
  }

  TTree* tree = (TTree*)f->Get("Passed");
  Double_t Weight;
  Double_t CHODTime1;
  Double_t CHODTime2;
  Double_t CDA;
  Double_t Zvertex;
  Double_t CDALine;
  Double_t ZCDALine;
  Double_t BeamlineDist;
  Double_t xSR;
  Double_t ySR;
  Double_t MuEoP;
  Double_t PiEoP;
  Double_t R;
  Double_t energyPi;
  Double_t energyMu;
  Double_t invMass;
  Double_t L0TPTime;
  TVector3 *Mom1 = new TVector3();
  TVector3 *Mom2 = new TVector3();
  TVector3 *TotMom = new TVector3();
  TVector3 *Vertex = new TVector3();
  TVector3 *threeMomPi = new TVector3();
  TVector3 *threeMomMu = new TVector3();
  Bool_t Target;
  Bool_t K3pi;
  Bool_t autoPass;
  Int_t Assoc;
  //TRecoCedarCandidate *KTAGcand;
  
  tree->SetBranchAddress("Weight", &Weight);
  tree->SetBranchAddress("CHODTime1", &CHODTime1);
  tree->SetBranchAddress("CHODTime2", &CHODTime2);
  tree->SetBranchAddress("CDA", &CDA);
  tree->SetBranchAddress("Zvertex", &Zvertex);
  tree->SetBranchAddress("CDALine", &CDALine);
  tree->SetBranchAddress("ZCDALine", &ZCDALine);
  tree->SetBranchAddress("BeamlineDist", &BeamlineDist);
  tree->SetBranchAddress("xSR", &xSR);
  tree->SetBranchAddress("ySR", &ySR);
  tree->SetBranchAddress("MuEoP", &MuEoP);
  tree->SetBranchAddress("PiEoP", &PiEoP);
  tree->SetBranchAddress("R", &R);
  tree->SetBranchAddress("energyPi", &energyPi);
  tree->SetBranchAddress("energyMu", &energyMu);
  tree->SetBranchAddress("invMass", &invMass);
  tree->SetBranchAddress("L0TPTime", &L0TPTime);
  tree->SetBranchAddress("Mom1", &Mom1);
  tree->SetBranchAddress("Mom2", &Mom2);
  tree->SetBranchAddress("TotMom", &TotMom);
  tree->SetBranchAddress("Vertex", &Vertex);
  tree->SetBranchAddress("threeMomPi", &threeMomPi);
  tree->SetBranchAddress("threeMomMu", &threeMomMu);
  tree->SetBranchAddress("Target", &Target);
  tree->SetBranchAddress("K3pi", &K3pi);
  tree->SetBranchAddress("autoPass", &autoPass);
  tree->SetBranchAddress("Assoc", &Assoc);
  //tree->SetBranchAddress("KTAGcand", KTAGcand);
    
  TH1D *hDistTime = new TH1D("hDistTime", "Parasitic background studies", 100, 0., 1000.);
  TH1D *hTimeTime = new TH1D("hTimeTime", "Combinatorial background studies", 100, -15., 15.);
  TH1D *hZTime = new TH1D("hZTime", "Prompt background studies", 500, 100., 190.);
  TH2D *hDistvsMassTime = new TH2D("hDistvsMassTime", "Parasitic background studies", 200, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRTime = new TH2D("hSRTime", "Signal region", 500, -50., 50., 50, 0., 0.1);

  TH1D *hDistSR = new TH1D("hDistSR", "Parasitic background studies", 100, 0., 1000.);
  TH1D *hTimeSR = new TH1D("hTimeSR", "Combinatorial background studies", 100, -15., 15.);
  TH1D *hZSR = new TH1D("hZSR", "Prompt background studies", 500, 100., 190.);
  TH2D *hDistvsMassSR = new TH2D("hDistvsMassSR", "Parasitic background studies", 200, 0.2, 2., 50, 0., 1000.);
  TH2D *hSRSR = new TH2D("hSRSR", "Signal region", 500, -50., 50., 50, 0., 0.1);
  
  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    if ((ZCDALine < -10000.) || (ZCDALine > 35000.) || ((ZCDALine >= -10000. && ZCDALine <= 35000.) && CDALine > 40.)) {
      hDistSR->Fill(BeamlineDist, Weight);
      hTimeSR->Fill(CHODTime1-CHODTime2, Weight);
      hZSR->Fill(Zvertex/1000., Weight);
      hDistvsMassSR->Fill(invMass/1000., BeamlineDist, Weight);
      hSRSR->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if ((CHODTime1-CHODTime2 < -3. && CHODTime1-CHODTime2 > -5.) || (CHODTime1-CHODTime2 > 3. && CHODTime1-CHODTime2 < 5.)) {
      hDistTime->Fill(BeamlineDist, Weight);
      hTimeTime->Fill(CHODTime1-CHODTime2, Weight);
      hZTime->Fill(Zvertex/1000., Weight);
      hDistvsMassTime->Fill(invMass/1000., BeamlineDist, Weight);
      hSRTime->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
  }

  Save(path, c, hDistSR, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimeSR, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZSR, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hDistvsMassSR, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRSR, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path, c, hDistTime, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hTimeTime, "Track time difference [ns]", labelSize, titleSize);
  Save(path, c, hZTime, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path, c, hDistvsMassTime, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path, c, hSRTime, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  tree->ResetBranchAddresses();
}
