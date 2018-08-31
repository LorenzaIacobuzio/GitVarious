Double_t labelSize = 0.05;
Double_t titleSize = 0.07;

void TH1Cosmetics(TH1* h1, Double_t labelSize, Double_t titleSize) {

  h1->SetFillColor(38);
  h1->SetTitleSize(titleSize, "t");
  h1->GetXaxis()->SetTitleSize(labelSize);
  h1->GetXaxis()->SetLabelSize(labelSize);
  h1->GetYaxis()->SetTitleSize(labelSize);
  h1->GetYaxis()->SetLabelSize(labelSize);
  h1->GetXaxis()->SetTitleOffset(1.4);
  h1->GetYaxis()->SetTitleOffset(1.4);
  gPad->Update();
  gStyle->SetStatW(0.2);
  gStyle->SetStatH(0.2);
  gStyle->SetOptStat(0);
  gPad->SetLogy(0);
  gPad->Update();
}

void Superimpose(TTree *t, TFile *f, TCanvas *c, TString path, TString variable, TString conditions, Int_t n, TString histoTitle, TString axisTitle, TString imageTitle, TString gaiaTitle) {

  t->Draw(variable, conditions, "", n, 0);
  TH1D *hMio = (TH1D*)gPad->GetPrimitive("h");
  TH1D* hGaia = (TH1D*)f->Get(gaiaTitle);
  hGaia->Scale(hMio->Integral()/hGaia->Integral());
  TH1Cosmetics(hMio, labelSize, titleSize);
  hMio->GetXaxis()->SetTitle(axisTitle);
  hMio->SetTitle(histoTitle);
  hMio->Draw();
  hGaia->Draw("sames");  
  gStyle->SetOptStat(0);
  gPad->Update();
  auto legend = new TLegend(0.71, 0.72, 0.98, 0.93);                                           
  legend->AddEntry(hGaia, "Toy MC (G. Lanfranchi)");                                           
  legend->AddEntry(hMio, "Full MC (L. Iacobuzio)");                                            
  legend->Draw();                                                                              
  c->Update();             
  c->SaveAs(path + imageTitle + "_comp.pdf");                                              
  c->SaveAs(path + imageTitle + "_comp.png");
}

void ToyComp(TString dir, TString histo1, Int_t n, TString prodCond) {

  TCanvas *c = new TCanvas();

  c->SetBottomMargin(0.2);
  
  // FIRST STEP HISTOS
    
  // One value plots for selection analyzer

  TString path = "";
  
  if (dir != "")
    path = dir;
  else
    path = "/home/li/Dropbox/PhD/Talks and papers/Notes/MCnote/images/Plots/2/HeavyNeutrinoScan/ToyMC/";

  TFile *f = TFile::Open(histo1);
  TTree* t = (TTree*)f->Get("MC");
  TFile *f1 = TFile::Open("/home/li/Desktop/NewHistos/Lorenza_D2Nmu_Npimu_comparison.root");
  TFile *f2 = TFile::Open("/home/li/Desktop/NewHistos/Lorenza_D2NKmu_Npimu_toys.root");
  TString prodCondition = "";

  if (prodCond.Contains("target"))
    prodCondition = " && fKineParts.fProdPos.Z()/1000. >= -0.4 && fKineParts.fProdPos.Z()/1000. <= 0.4";
  else if (prodCond.Contains("tax"))
    prodCondition = " && fKineParts.fProdPos.Z()/1000. >= 20. && fKineParts.fProdPos.Z()/1000. <= 27.";
  else
    prodCondition = "";
  
  // DS pN
  
  Superimpose(t, f1, c, path, "Generated.fKineParts.fMomGigaTrackerExit.Z()/1000. >> h(50, 0., 300.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"DS\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "N momentum from D_{S} #rightarrow N#mu, good only", "P [GeV/c]", "DS/pN", "hN_p_prodinacc;1");

  // DS ptN
    
  Superimpose(t, f1, c, path, "Generated.fKineParts.fMomGigaTrackerExit.T()/1000. >> h(20, 0., 2.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"DS\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "N transverse momentum from D_{S} #rightarrow N#mu, good only", "P_{t} [GeV/c]", "DS/ptN", "hN_p_prodinacc;2");

  // DS pMu
  
  Superimpose(t, f1, c, path, "Generated.fKineParts.fMomSpectrometerEntry.X()/1000. >> h(50, 0., 300.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"DS\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "Muon daughter momentum from D_{S} #rightarrow N#mu, good only", "P [GeV/c]", "DS/pMu", "hProd_p_prodinacc;1");
  
  // DS ptMu
  
  Superimpose(t, f1, c, path, "Generated.fKineParts.fMomSpectrometerEntry.Y()/1000. >> h(20, 0., 2.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"DS\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "Muon daughter transverse momentum from D_{S} #rightarrow N#mu, good only", "P_{t} [GeV/c]", "DS/ptMu", "hProd_pt_prodinacc;1");
  
  // D0 pN
  
  Superimpose(t, f2, c, path, "Generated.fKineParts.fMomMNP33Entry.Z()/1000. >> h(50, 0., 300.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"D0\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "N momentum from D^{ 0} #rightarrow K#muN, good only", "P [GeV/c]", "D0/pN", "hN_p_prodinacc;1");

  // D0 ptN

  Superimpose(t, f2, c, path, "Generated.fKineParts.fMomMNP33Entry.T()/1000. >> h(20, 0., 2.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"D0\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "N transverse momentum from D^{ 0} #rightarrow K#muN, good only", "P_{t} [GeV/c]", "D0/ptN", "hN_p_prodinacc;2");

  // D0 pMu
  
  Superimpose(t, f2, c, path, "Generated.fKineParts.fMomMNP33Exit.X()/1000. >> h(50, 0., 300.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"D0\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "Muon daughter momentum from D^{ 0} #rightarrow K#muN, good only", "P [GeV/c]", "D0/pMu", "hProd_p_prodinacc;1");

  // D0 ptMu
  
  Superimpose(t, f2, c, path, "Generated.fKineParts.fMomMNP33Exit.Y()/1000. >> h(20, 0., 2.)", "Generated.fKineParts.fPDGcode == 999 && Generated.fKineParts.fParticleName.Contains(\"D0\") && Generated.fKineParts.fParticleName.Contains(\"mu\") && Generated.fKineParts.fEndProcessName.Contains(\"good\")" + prodCondition, n, "Muon daughter transverse momentum from D^{ 0} #rightarrow K#muN, good only", "P_{t} [GeV/c]", "D0/ptMu", "hProd_pt_prodinacc;1");
}
