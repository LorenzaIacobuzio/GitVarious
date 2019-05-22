Double_t Par[6] = {-0.0600938, 0.909681, -2.16879, 2.03637, -0.806702, 0.112638};
Double_t minSigma = 0.00206844/4.;

void Save(TString path, TCanvas *c, TH1D* h, TString x, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();
  
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

  if (name.Contains("TimeSR") || name.Contains("ZSR"))
    gPad->SetLogy();
  else
    gPad->SetLogy(0);
  
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");

  return;
}

void Save(TString path, TCanvas *c, TH2D* h, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();
  
  if (h->GetEntries() < 100)
    h->Draw("text");
  else
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

  return;
}

void Save(TString path, TCanvas *c, TGraph* g, TString name, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  if (name.Contains("gMassMC"))
    g->Draw("AP*");
  else
    g->Draw("APL");
  
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.4);
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + name + ".pdf");
  c->SaveAs(path + name + ".png");

  return;
}
	  
TGraph* WindowScanFixedSigma(TH1D* h, Double_t massSigma, Double_t sigmaStep) {

  Int_t firstBin = h->FindFirstBinAbove(0.,1);
  Double_t firstBinValue = h->GetBinCenter(firstBin)-h->GetBinWidth(firstBin)/2.;
  Int_t lastBin = h->FindLastBinAbove(0.,1);
  Double_t lastBinValue = h->GetBinCenter(lastBin)-h->GetBinWidth(lastBin)/2.;
  Int_t counter = 0;
  TGraph *g = new TGraph();
  
  for (Double_t massHyp = firstBinValue-2.*sigmaStep*massSigma; massHyp <= lastBinValue+2.*sigmaStep*massSigma; massHyp += sigmaStep*massSigma) {
    Double_t nBkg = 0.;
    for (Double_t i = massHyp-massSigma; i <= massHyp+massSigma; i += massSigma) {
      nBkg += h->GetBinContent(h->FindBin(i+massSigma*0.000001));
    }
    g->SetPoint(counter, massHyp, nBkg);
    counter++;
  }
  
  return g;	 
}

TGraph* WindowScanNoFixedSigma(TH1D* h, Int_t howManySigmas) {

  Int_t firstBin = h->FindFirstBinAbove(0.,1);
  Double_t firstBinValue = h->GetBinCenter(firstBin)-h->GetBinWidth(firstBin)/2.;
  Int_t lastBin = h->FindLastBinAbove(0.,1);
  Double_t lastBinValue = h->GetBinCenter(lastBin)-h->GetBinWidth(lastBin)/2.;
  Int_t counter = 0;
  TGraph *res = new TGraph();

  for (Double_t massHyp = firstBinValue-minSigma; massHyp <= lastBinValue+minSigma; massHyp += minSigma) {
    Double_t sigma = Par[0]*TMath::Power(massHyp, 5.) + Par[1]*TMath::Power(massHyp, 4.) + Par[2]*TMath::Power(massHyp, 3.) + Par[3]*TMath::Power(massHyp, 2.) + Par[4]*TMath::Power(massHyp, 1.) + Par[5];
    TAxis *axis = h->GetXaxis();
    Int_t bmin = axis->FindBin(massHyp-howManySigmas*sigma);
    Int_t bmax = axis->FindBin(massHyp+howManySigmas*sigma);
    Double_t nBkg = h->Integral(bmin, bmax);
    nBkg -= h->GetBinContent(bmin)*(massHyp-howManySigmas*sigma-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    nBkg -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-(massHyp+howManySigmas*sigma))/axis->GetBinWidth(bmax);
    res->SetPoint(counter, massHyp, nBkg);
    counter++;
  }
  
  return res;
}
  
void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, Double_t &counterComb, Double_t &counterPar, Double_t &counterPrompt, Double_t &counterPromptPosNeg, Double_t &counterPromptPosNegFV) {
  
  TString path = "";
  TString analyzer = "";
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t ZCDALineMax = 35000.;
  Double_t ZCDALineMin = -10000.;
  Double_t CDALineMax = 40.;
  Double_t Time1Max = -2.;
  Double_t Time1Min = -4.;
  Double_t Time2Max = 4.;
  Double_t Time2Min = 2.;
  Double_t CDAMax = 10.;
  Double_t ZVertexMax = 120000.;
  Double_t ZVertexMin = 100000.;
  Double_t BeamdistMax = 150.;
  Double_t BeamdistMin = 100.;
  Double_t ZEndFV = 180000.;
  Double_t xGTKMin = -50.;
  Double_t xGTKMax = 50.;
  Double_t yGTKMin = -40;
  Double_t yGTKMax = 40;
  Double_t massSigma = 0.0045; //GeV
  Double_t sigmaStep = 0.25; // step for mass window center
  Double_t massMin = 0.2;
  Double_t massMax = 2.;
  Double_t massStep = 0.01;
  Int_t nMassBins = (massMax - massMin)/massSigma;
  Double_t binWidth = (massMax - massMin)/nMassBins;
  
  if (an.Contains("Pos"))
    analyzer = "Pos";
  else if (an.Contains("Neg"))
    analyzer = "Neg";
  else if (!an.Contains("Pos") && !an.Contains("Neg"))
    analyzer = "Zero";
  
  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("2016"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2016/" + analyzer + "/";
    else if (histo1.Contains("2017"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2017/" + analyzer + "/";
    else if (histo1.Contains("2018"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/2018/" + analyzer + "/";    
    else if (histo1.Contains("Data"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/" + analyzer + "/";
    else
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/MC/" + analyzer + "/";
  }
  
  TFile *f = TFile::Open(histo1);
  
  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    _exit(1);
  }

  TString treeName = an + "Passed";
  TTree* tree = (TTree*)f->Get(treeName);
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
  Double_t xGTK31;
  Double_t yGTK31;
  Double_t xGTK32;
  Double_t yGTK32;
  Double_t BeamCDA1;
  Double_t BeamCDA2;
  Double_t EtotLKr;
  TVector3 *Mom1 = new TVector3();
  TVector3 *Mom2 = new TVector3();
  TVector3 *TotMom = new TVector3();
  TVector3 *Vertex = new TVector3();
  TVector3 *Pos1 = new TVector3();
  TVector3 *Pos2 = new TVector3();
  TVector3 *BeamVtx1 = new TVector3();
  TVector3 *BeamVtx2 = new TVector3();
  Bool_t Target;
  Bool_t K3pi;
  Bool_t autoPass;
  Int_t Assoc;
  Bool_t CHANTIAssoc1;
  Bool_t CHANTIAssoc2;
  Int_t Charge1;
  Int_t Charge2;
  Int_t nSec;
  Int_t nCHOD;
  Int_t nNewCHOD;
  Int_t nLKr;
  std::vector<Double_t> *RICHInfo1 = new std::vector<Double_t>;
  std::vector<Double_t> *RICHInfo2 = new std::vector<Double_t>;
  
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
  tree->SetBranchAddress("xGTK31", &xGTK31);
  tree->SetBranchAddress("yGTK31", &yGTK31);
  tree->SetBranchAddress("xGTK32", &xGTK32);
  tree->SetBranchAddress("yGTK32", &yGTK32);
  tree->SetBranchAddress("BeamCDA1", &BeamCDA1);
  tree->SetBranchAddress("BeamCDA2", &BeamCDA2);
  tree->SetBranchAddress("EtotLKr", &EtotLKr);
  tree->SetBranchAddress("Mom1", &Mom1);
  tree->SetBranchAddress("Mom2", &Mom2);
  tree->SetBranchAddress("TotMom", &TotMom);
  tree->SetBranchAddress("Vertex", &Vertex);
  tree->SetBranchAddress("Pos1", &Pos1);
  tree->SetBranchAddress("Pos2", &Pos2);
  tree->SetBranchAddress("BeamVtx1", &BeamVtx1);
  tree->SetBranchAddress("BeamVtx2", &BeamVtx2);
  tree->SetBranchAddress("Target", &Target);
  tree->SetBranchAddress("K3pi", &K3pi);
  tree->SetBranchAddress("autoPass", &autoPass);
  tree->SetBranchAddress("Assoc", &Assoc);
  tree->SetBranchAddress("CHANTIAssoc1", &CHANTIAssoc1);
  tree->SetBranchAddress("CHANTIAssoc2", &CHANTIAssoc2);
  tree->SetBranchAddress("Charge1", &Charge1);
  tree->SetBranchAddress("Charge2", &Charge2);
  tree->SetBranchAddress("nSec", &nSec);
  tree->SetBranchAddress("nCHOD", &nCHOD);
  tree->SetBranchAddress("nNewCHOD", &nNewCHOD);
  tree->SetBranchAddress("nLKr", &nLKr);
  tree->SetBranchAddress("RICHInfo1", &RICHInfo1);
  tree->SetBranchAddress("RICHInfo2", &RICHInfo2);

  // Combinatorial
    
  TH1D *hDistComb = new TH1D("hDistComb", "Combinatorial background studies", 100, 0., 1000.);
  TH1D *hTimeComb = new TH1D("hTimeComb", "Combinatorial background studies", 100, -15., 15.);
  TH1D *hZComb = new TH1D("hZComb", "Combinatorial background studies", 500, 100., 190.);
  TH1D *hCDAComb = new TH1D("hCDAComb", "Combinatorial background studies", 100, 0., 0.5);
  TH1D *hInvMassComb = new TH1D("hInvMassComb", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombKolm = new TH1D("hInvMassCombKolm", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSR = new TH1D("hInvMassCombSR", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSREnriched = new TH1D("hInvMassCombSREnriched", "Enriched combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hMomPiPosComb = new TH1D("hMomPiPosComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuPosComb = new TH1D("hMomMuPosComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomPiNegComb = new TH1D("hMomPiNegComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuNegComb = new TH1D("hMomMuNegComb", "Combinatorial background studies", 100, 0., 200.);
  TH2D *hGTK3IDPiPosComb = new TH2D("hGTK3IDPiPosComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDMuPosComb = new TH2D("hGTK3IDMuPosComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDPiNegComb = new TH2D("hGTK3IDPiNegComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDMuNegComb = new TH2D("hGTK3IDMuNegComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hDistvsMassComb = new TH2D("hDistvsMassComb", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  TH2D *hGTK3XYPiPosComb = new TH2D("hGTK3XYPiPosComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYMuPosComb = new TH2D("hGTK3XYMuPosComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYPiNegComb = new TH2D("hGTK3XYPiNegComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYMuNegComb = new TH2D("hGTK3XYMuNegComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hSRComb = new TH2D("hSRComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalComb = new TH2D("hSRFinalComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalCombEnriched = new TH2D("hSRFinalCombEnriched", "Enriched signal region", 500, -50., 50., 50, 0., 0.1);
  TGraph *gBkg1SigmaComb = new TGraph();
  gBkg1SigmaComb->SetNameTitle("gBkg1SigmaComb", "Combinatorial background studies");
  TGraph *gBkg2SigmaComb = new TGraph();
  gBkg2SigmaComb->SetNameTitle("gBkg2SigmaComb", "Combinatorial background studies");
  
  // Prompt

  TH1D *hDistPrompt = new TH1D("hDistPrompt", "Prompt background studies", 100, 0., 1000.);
  TH1D *hTimePrompt = new TH1D("hTimePrompt", "Prompt background studies", 100, -15., 15.);
  TH1D *hZPrompt = new TH1D("hZPrompt", "Prompt background studies", 500, 100., 190.);
  TH1D *hCDAPrompt = new TH1D("hCDAPrompt", "Prompt background studies", 100, 0., 0.5);
  TH1D *hInvMassPrompt = new TH1D("hInvMassPrompt", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptSR = new TH1D("hInvMassPromptSR", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hMomPiPrompt = new TH1D("hMomPiPrompt", "Prompt background studies", 100, 0., 200.);
  TH1D *hMomMuPrompt = new TH1D("hMomMuPrompt", "Prompt background studies", 100, 0., 200.);
  TH2D *hDistvsMassPrompt = new TH2D("hDistvsMassPrompt", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  TH2D *hSRPrompt = new TH2D("hSRPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPrompt = new TH2D("hSRFinalPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TGraph *gBkg1SigmaPrompt = new TGraph();
  gBkg1SigmaComb->SetNameTitle("gBkg1SigmaPrompt", "Prompt background studies");
  TGraph *gBkg2SigmaPrompt = new TGraph();
  gBkg2SigmaComb->SetNameTitle("gBkg2SigmaPrompt", "Prompt background studies");

  // Parasitic

  TH1D *hDistPar = new TH1D("hDistPar", "Parasitic background studies", 100, 0., 1000.);
  TH1D *hTimePar = new TH1D("hTimePar", "Parasitic background studies", 100, -15., 15.);
  TH1D *hZPar = new TH1D("hZPar", "Parasitic background studies", 500, 100., 190.);
  TH1D *hCDAPar = new TH1D("hCDAPar", "Parasitic background studies", 100, 0., 0.5);
  TH1D *hInvMassPar = new TH1D("hInvMassPar", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParSR = new TH1D("hInvMassParSR", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hMomPiPar = new TH1D("hMomPiPar", "Parasitic background studies", 100, 0., 200.);
  TH1D *hMomMuPar = new TH1D("hMomMuPar", "Parasitic background studies", 100, 0., 200.);
  TH2D *hDistvsMassPar = new TH2D("hDistvsMassPar", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  TH2D *hDistvsMassParSR = new TH2D("hDistvsMassParSR", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  TH2D *hSRPar = new TH2D("hSRPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPar = new TH2D("hSRFinalPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TGraph *gBkg1SigmaPar = new TGraph();
  gBkg1SigmaPar->SetNameTitle("gBkg1SigmaPar", "Parasitic background studies");
  TGraph *gBkg2SigmaPar = new TGraph();
  gBkg2SigmaPar->SetNameTitle("gBkg2SigmaPar", "Parasitic background studies");

  // SR

  TH1D *hDistSR = new TH1D("hDistSR", "Signal region background studies", 100, 0., 1000.);
  TH1D *hTimeSR = new TH1D("hTimeSR", "Signal region background studies", 100, -15., 15.);
  TH1D *hZSR = new TH1D("hZSR", "Signal region background studies", 500, 100., 190.);
  TH1D *hCDASR = new TH1D("hCDASR", "Signal region background studies", 100, 0., 0.5);
  TH1D *hInvMassSR = new TH1D("hInvMassSR", "Signal region background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hMomPiSR = new TH1D("hMomPiSR", "Signal region background studies", 100, 0., 200.);
  TH1D *hMomMuSR = new TH1D("hMomMuSR", "Signal region background studies", 100, 0., 200.);
  TH2D *hDistvsMassSR = new TH2D("hDistvsMassSR", "Signal region background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  TH2D *hSRSR = new TH2D("hSRSR", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // MC only

  TH1D *hInvMassMC = new TH1D("hInvMassMC", "MC studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TGraph *gMassMC = new TGraph();
  gMassMC->SetNameTitle("gMassMC", "MC studies");
  
  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    
    if (i%1000 == 0)
      cout<<"Processing event n."<<i<<endl;

    // MC studies

    if (!histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data") && !an.Contains("Pos") && !an.Contains("Neg")) {
      hInvMassMC->Fill(invMass/1000.);
    }
    
    // 0-charge
    
    if ((ZCDALine < ZCDALineMin || ZCDALine > ZCDALineMax || (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine > CDALineMax)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events outside blinded region
      hDistSR->Fill(BeamlineDist, Weight);
      hTimeSR->Fill(CHODTime1-CHODTime2, Weight);
      hZSR->Fill(Zvertex/1000., Weight);
      hCDASR->Fill(CDALine/1000., Weight);
      hInvMassSR->Fill(invMass/1000., Weight);
      hMomPiSR->Fill(Mom1->Mag()/1000., Weight);
      hMomMuSR->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassSR->Fill(invMass/1000., BeamlineDist, Weight);
      hSRSR->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside time sidebands
      hDistComb->Fill(BeamlineDist, Weight);
      hTimeComb->Fill(CHODTime1-CHODTime2, Weight);
      hZComb->Fill(Zvertex/1000., Weight);
      hCDAComb->Fill(CDALine/1000., Weight);
      hInvMassComb->Fill(invMass/1000., Weight);            
      hDistvsMassComb->Fill(invMass/1000., BeamlineDist, Weight);
      hSRComb->Fill(ZCDALine/1000., CDALine/1000., Weight);

      // Study spike at 75 GeV

      if (Assoc == 1 && Charge1 == 1) { // mu+pi-
	hMomMuPosComb->Fill(Mom1->Mag()/1000., Weight);
	hMomPiNegComb->Fill(Mom2->Mag()/1000., Weight);
      }
      else if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	hMomMuNegComb->Fill(Mom1->Mag()/1000., Weight);
	hMomPiPosComb->Fill(Mom2->Mag()/1000., Weight);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	hMomPiPosComb->Fill(Mom1->Mag()/1000., Weight);
	hMomMuNegComb->Fill(Mom2->Mag()/1000., Weight);
      }
      else if (Assoc ==	2 && Charge1 == -1) { // pi-mu+
	hMomPiNegComb->Fill(Mom1->Mag()/1000., Weight);
	hMomMuPosComb->Fill(Mom2->Mag()/1000., Weight);
      }

      // Select momentum in [70-80] GeV to study spike

      if (Assoc == 1 && Charge1 == 1) { // mu+pi-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYMuPosComb->Fill(xGTK31, yGTK31);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYPiNegComb->Fill(xGTK32, yGTK32);
      }
      else if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYMuNegComb->Fill(xGTK31, yGTK31);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYPiPosComb->Fill(xGTK32, yGTK32);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYPiPosComb->Fill(xGTK31, yGTK31);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYMuNegComb->Fill(xGTK32, yGTK32);
      }
      else if (Assoc == 2 && Charge1 == -1) { // pi-mu+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYPiNegComb->Fill(xGTK31, yGTK32);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYMuPosComb->Fill(xGTK32, yGTK32);
      }

      // Select momentum in [70-80] GeV and xy @ GTK3 to study spike
      
      if (RICHInfo1->at(0) == 99)
	RICHInfo1->at(0) = 5;

      if (RICHInfo2->at(0) == 99)
	RICHInfo2->at(0) = 5;

      if (Assoc == 1 && Charge1 == 1) { // mu+pi-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDMuPosComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDPiNegComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDMuNegComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDPiPosComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
      if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDPiPosComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDMuNegComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
      if (Assoc == 2 && Charge1 == -1) { // pi-mu+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDPiNegComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDMuPosComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
    }
    
    if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside vertex sidebands
      hDistPrompt->Fill(BeamlineDist, Weight);
      hTimePrompt->Fill(CHODTime1-CHODTime2, Weight);
      hZPrompt->Fill(Zvertex/1000., Weight);
      hCDAPrompt->Fill(CDALine/1000., Weight);
      hInvMassPrompt->Fill(invMass/1000., Weight);
      hMomPiPrompt->Fill(Mom1->Mag()/1000., Weight);
      hMomMuPrompt->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassPrompt->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // all events inside beam distance sidebands
      hDistPar->Fill(BeamlineDist, Weight);
      hTimePar->Fill(CHODTime1-CHODTime2, Weight);
      hZPar->Fill(Zvertex/1000., Weight);
      hCDAPar->Fill(CDALine/1000., Weight);
      hInvMassPar->Fill(invMass/1000., Weight);
      hMomPiPar->Fill(Mom1->Mag()/1000., Weight);
      hMomMuPar->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassPar->Fill(invMass/1000., BeamlineDist, Weight);
      hSRPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    if (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine <= CDALineMax) { // all events inside blinded region...
      if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < 100. && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and time sidebands (enriched sample)
	hSRFinalCombEnriched->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSREnriched->Fill(invMass/1000., Weight);
      }
      if (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)) && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and time sidebands (not enriched sample to study spike at 75 GeV)
	counterComb++;
	hSRFinalComb->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSR->Fill(invMass/1000., Weight);
      }
      if (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and beam distance sidebands
	hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPar++;
	hInvMassParSR->Fill(invMass/1000., Weight);
	hDistvsMassParSR->Fill(invMass/1000., BeamlineDist, Weight);
      }
      if (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax && !an.Contains("Pos") && !an.Contains("Neg")) { // ...and vertex sidebands (0-charge)
	hSRFinalPrompt->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterPrompt++;
	hInvMassPromptSR->Fill(invMass/1000., Weight);	
      }
      if ((an.Contains("Pos") || an.Contains("Neg")) && Zvertex >= ZVertexMin && Zvertex <= ZVertexMax && CDA < CDAMax) { // ...and vertex sidebands (2-charge)
	counterPromptPosNeg++;
      }
      if ((an.Contains("Pos") || an.Contains("Neg")) && Zvertex >= ZVertexMax && Zvertex <= ZEndFV && CDA < CDAMax) { // ...and outside vertex sidebands (2-charge)
	counterPromptPosNegFV++;
      }
    }

    if (i>10000)
      break;

  }

  if (!an.Contains("Pos") && !an.Contains("Neg")) {
    
    // Sigma vs mass (MC studies)
    
    if (!histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data")) { 
      Int_t counter = 0;
      Double_t sigma = 0.;
      Double_t sigmaMin = 999.;
      hInvMassMC->Draw();
      for (Double_t mass = massMin+massStep; mass <= massMax-massStep; mass += massStep) {
	TF1 *f1 = new TF1("f1", "gaus", mass-massStep, mass+massStep);
	hInvMassMC->Fit("f1", "R");
	sigma = f1->GetParameter(2);
	gMassMC->SetPoint(counter, mass, sigma);
	if (sigma < sigmaMin && sigma != 0.)
	  sigmaMin = sigma;
	counter++;
      }

      gMassMC->Fit("pol5");
      TF1 *f = gMassMC->GetFunction("pol5");
      cout<<"MC signal mass resolution, fit parameters: "<<f->GetParameter(0)<<" "<<f->GetParameter(1)<<" "<<f->GetParameter(2)<<" "<<f->GetParameter(3)<<" "<<f->GetParameter(4)<<" "<<f->GetParameter(5)<<" and minimum sigma over whole range: "<<sigmaMin<<" "<<endl;
    }
    
    cout<<"Kolmogorov test for SR sample: "<<hInvMassCombSR->KolmogorovTest(hInvMassComb)<<" and enriched SR sample: "<<hInvMassCombSREnriched->KolmogorovTest(hInvMassComb)<<endl;
    hInvMassCombKolm = (TH1D*)hInvMassComb->Clone();
    hInvMassCombKolm->SetName("hInvMassCombKolm");
    hInvMassCombKolm->Scale(hInvMassCombSR->Integral()/hInvMassCombKolm->Integral());
    /*    
    TH1D *h = new TH1D("h","h",3,0.1,0.7);
    h->Fill(0.4);
    h->Fill(0.4);
    h->Fill(0.6);
    TGraph *mass = new TGraph(5);
    mass->SetPoint(0,0.2,0.2);
    mass->SetPoint(1,0.3,0.2);
    mass->SetPoint(2,0.4,0.2);
    mass->SetPoint(3,0.5,0.1);
    mass->SetPoint(4,0.6,0.1);
    TGraph *g = new TGraph();
    g = WindowScanNoFixedSigma(h,mass,1);
    */
    gBkg1SigmaComb = WindowScanNoFixedSigma(hInvMassCombKolm, 2.);
    gBkg2SigmaComb = WindowScanNoFixedSigma(hInvMassCombKolm, 2.);
  }
      
  Save(path + "SR/", c, hDistSR, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "SR/", c, hTimeSR, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "SR/", c, hZSR, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "SR/", c, hCDASR, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "SR/", c, hInvMassSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "SR/", c, hMomPiSR, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "SR/", c, hMomMuSR, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "SR/", c, hDistvsMassSR, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "SR/", c, hSRSR, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path + "Comb/", c, hDistComb, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hTimeComb, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "Comb/", c, hZComb, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "Comb/", c, hCDAComb, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Comb/", c, hInvMassComb, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Comb/", c, hInvMassCombKolm, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Comb/", c, hInvMassCombSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Comb/", c, hInvMassCombSREnriched, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Comb/", c, hMomPiPosComb, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Comb/", c, hMomMuPosComb, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Comb/", c, hMomPiNegComb, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Comb/", c, hMomMuNegComb, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3IDPiPosComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3IDMuPosComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3IDPiNegComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3IDMuNegComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hDistvsMassComb, "Reconstructed HNL mass [GeV/c^{2}]", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3XYPiPosComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3XYMuPosComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3XYPiNegComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hGTK3XYMuNegComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "Comb/", c, hSRComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Comb/", c, hSRFinalComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Comb/", c, hSRFinalCombEnriched, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  //Save(path + "Comb/", c, gBkg1SigmaComb, "gBkg1SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  //Save(path + "Comb/", c, gBkg2SigmaComb, "gBkg2SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  
  Save(path + "Prompt/", c, hDistPrompt, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Prompt/", c, hTimePrompt, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "Prompt/", c, hZPrompt, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "Prompt/", c, hCDAPrompt, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Prompt/", c, hInvMassPrompt, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Prompt/", c, hInvMassPromptSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Prompt/", c, hMomPiPrompt, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Prompt/", c, hMomMuPrompt, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Prompt/", c, hDistvsMassPrompt, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Prompt/", c, hSRPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Prompt/", c, hSRFinalPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  //Save(path + "Prompt/", c, gBkg1SigmaPrompt, "gBkg1SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  //Save(path + "Prompt/", c, gBkg2SigmaPrompt, "gBkg2SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  Save(path + "Par/", c, hDistPar, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Par/", c, hTimePar, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "Par/", c, hZPar, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "Par/", c, hCDAPar, "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Par/", c, hInvMassPar, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Par/", c, hInvMassParSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "Par/", c, hMomPiPar, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Par/", c, hMomMuPar, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "Par/", c, hDistvsMassPar, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Par/", c, hDistvsMassParSR, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "Par/", c, hSRPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "Par/", c, hSRFinalPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  //Save(path + "Par/", c, gBkg1SigmaPar, "gBkg1SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  //Save(path + "Par/", c, gBkg2SigmaPar, "gBkg2SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  if (!an.Contains("Pos") && !an.Contains("Neg") && !histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data")) {
    Save(path + "MC/", c, hInvMassMC, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
    Save(path + "MC/", c, gMassMC, "gMassMC", "Reconstructed invariant mass [GeV/c^{2}]", "Sigma [GeV/c^{2}]", labelSize, titleSize);
  }
  
  tree->ResetBranchAddresses();
  
  return;
}

void TreePlots(TString dir, TString histo1) {

  // dir = output dir, histo1 = histo to do cosmetics on

  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t counterComb = 0;
  Double_t counterPar = 0;
  Double_t counterPrompt = 0;
  Double_t counterPromptFV = 0;
  Double_t counterPromptPosNeg = 0;
  Double_t counterPromptPosNegFV = 0;
  
  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  if (histo1.Contains("2016") || histo1.Contains("2017") || histo1.Contains("2018") || histo1.Contains("Data")) {
    Analyzer(dir, histo1, "HeavyNeutrino", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
    //Analyzer(dir, histo1, "HeavyNeutrinoPos", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
    //Analyzer(dir, histo1, "HeavyNeutrinoNeg", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
    //counterPromptFV = counterPrompt*counterPromptPosNegFV/counterPromptPosNeg;
  }
  else
    Analyzer(dir, histo1, "HeavyNeutrino", c, counterComb, counterPar, counterPrompt, counterPromptPosNeg, counterPromptPosNegFV);
  
  cout<<"Number of events in SR and time sidebands: "<<counterComb<<", and beamdist sidebands: "<<counterPar<<", and Z sidebands (0-charge): "<<counterPrompt<<", and Z sidebands (2-charge): "<<counterPromptPosNeg<<", and FV (2-charge): "<<counterPromptPosNegFV<<", and FV (0-charge): "<<counterPromptFV<<endl;

  _exit(0);
}
