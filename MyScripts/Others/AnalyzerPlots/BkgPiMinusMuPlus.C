// Variables for sigma vs mass fit

Double_t Par[6] = {0., 0., 0., 0., 0., 0.};
Double_t minSigma = 0.;
Double_t massMin = 0.25;
Double_t massMax = 1.96;
Double_t massStep = 0.01;
Int_t nMassBins = (massMax - massMin)/massStep;
Double_t binWidth = (massMax - massMin)/(nMassBins);

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
    g->Draw("AL");
  
  g->GetXaxis()->SetTitleOffset(1.4);
  g->GetYaxis()->SetTitleOffset(1.5);
  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
  g->GetXaxis()->SetTitleSize(labelSize);
  g->GetYaxis()->SetTitleSize(labelSize);
  g->GetXaxis()->SetLabelSize(labelSize);
  g->GetYaxis()->SetLabelSize(labelSize);
  g->SetLineWidth(2);
  g->SetLineColor(kRed);
  gStyle->SetOptStat(0);
  gPad->Update();
  c->SaveAs(path + name + ".pdf");
  c->SaveAs(path + name + ".png");

  return;
}

// Scan on HNL mass for expected number of background events, with window around mass hypothesis at 1 or 2 sigma (fixed sigma)

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

// Scan on HNL mass for expected number of background events, with window around mass hypothesis at 1 or 2 sigma (non-fixed sigma: each MC mass hypothesis is fitted, sigma is plotted vs mass and fitted. Fit parameters are used to compute sigma for each mass hypothesis in this scan)

TGraph* WindowScanNoFixedSigma(TH1D* h, Int_t howManySigmas) {

  Int_t firstBin = h->FindFirstBinAbove(0.,1);
  Double_t firstBinValue = h->GetBinCenter(firstBin)-h->GetBinWidth(firstBin)/2.;
  Int_t lastBin = h->FindLastBinAbove(0.,1);
  Double_t lastBinValue = h->GetBinCenter(lastBin)-h->GetBinWidth(lastBin)/2.;
  Int_t counter = 0;
  TGraph *res = new TGraph();

  for (Double_t massHyp = massMin-minSigma; massHyp <= massMax+minSigma; massHyp += minSigma) {
    Double_t sigma = Par[5]*TMath::Power(massHyp, 5.) + Par[4]*TMath::Power(massHyp, 4.) + Par[3]*TMath::Power(massHyp, 3.) + Par[2]*TMath::Power(massHyp, 2.) + Par[1]*TMath::Power(massHyp, 1.) + Par[0];
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

// Sum graphs for total expected bkg

void SumGraphs(TGraph* g1, TGraph* g2, TGraph &g) {

  if (g1->GetN() != g2->GetN()) {
    cout<<g1->GetName()<<" and "<<g2->GetName()<<" do not have the same number of points!"<<endl;
    _exit(1);
  }
  
  for (Int_t i = 0; i < g1->GetN(); i++) {
    Double_t X = 0.;
    Double_t Y = 0.;
    Double_t Ytot = 0.;
    
    g1->GetPoint(i, X, Y);
    Ytot += Y;
    g2->GetPoint(i, X, Y);
    Ytot += Y;
    g.SetPoint(i, X, Ytot);
    X = 0.;
    Y = 0.;
    Ytot = 0.;
  }

  return;
}

void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, Double_t &counterComb, Double_t &counterPar, Double_t &counterPromptSBZ, Double_t &counterPromptFVZ, Double_t &counterPromptSBN, Double_t &counterPromptFVN, TGraph* gBkg1SigmaBuffer, TGraph* gBkg2SigmaBuffer, TGraph* gBkg1SigmaTot, TGraph* gBkg2SigmaTot) {

  // Generic variables
  
  TString path = "";
  TString analyzer = "";
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t ZCDALineMax = 40000.;//35000.;
  Double_t ZCDALineMin = -15000.;//-10000.;
  Double_t CDALineMax = 60.;//40.;
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
  Double_t yGTKMin = -30.;
  Double_t yGTKMax = 30.;
  Double_t massSigma = 0.0045; //GeV (fixed from 1 GeV mass);
  Double_t sigmaStep = 0.25; // step for mass window center

  // Finding paths
  
  if (an.Contains("Neg"))
    analyzer = "Neg";
  else if (!an.Contains("Pos") && !an.Contains("Neg"))
    analyzer = "Zero";
  
  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("Data"))
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/" + analyzer + "/";
    else
      path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/MC/" + analyzer + "/";
  }
  
  TFile *f = TFile::Open(histo1);
  
  if (f == 0) {
    cout << "Error: cannot open " << histo1 << endl;
    _exit(1);
  }

  // Variables for tree
  
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
  Double_t trueMass;
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
  Bool_t isControlTrigger;
  Int_t Charge1;
  Int_t Charge2;
  Int_t nSec;
  Int_t nCHOD;
  Int_t nNewCHOD;
  Int_t nLKr;
  std::vector<Double_t> *RICHInfo1 = new std::vector<Double_t>;
  std::vector<Double_t> *RICHInfo2 = new std::vector<Double_t>;
  std::string primitive;
  
  // Setting branches
  
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
  tree->SetBranchAddress("trueMass", &trueMass);
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
  tree->SetBranchAddress("isControlTrigger", &isControlTrigger);
  tree->SetBranchAddress("Charge1", &Charge1);
  tree->SetBranchAddress("Charge2", &Charge2);
  tree->SetBranchAddress("nSec", &nSec);
  tree->SetBranchAddress("nCHOD", &nCHOD);
  tree->SetBranchAddress("nNewCHOD", &nNewCHOD);
  tree->SetBranchAddress("nLKr", &nLKr);
  tree->SetBranchAddress("RICHInfo1", &RICHInfo1);
  tree->SetBranchAddress("RICHInfo2", &RICHInfo2);
  tree->SetBranchAddress("primitive", &primitive);

  // Setting histograms - splitting analysis for pi+mu- and pi-mu+ cases
  
  // Combinatorial - pi+mu-
    
  TH1D *hInvMassComb = new TH1D("hInvMassComb", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombKolm = new TH1D("hInvMassCombKolm", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSR = new TH1D("hInvMassCombSR", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  
  TH1D *hMomPiNegComb = new TH1D("hMomPiPosComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuPosComb = new TH1D("hMomMuPosComb", "Combinatorial background studies", 100, 0., 200.);
  TH2D *hGTK3XYPiNegComb = new TH2D("hGTK3XYPiNegComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYMuPosComb = new TH2D("hGTK3XYMuPosComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3IDPiNegComb = new TH2D("hGTK3IDPiNegComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDMuPosComb = new TH2D("hGTK3IDMuPosComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);

  TH2D *hSRComb = new TH2D("hSRComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalComb = new TH2D("hSRFinalComb", "Signal region", 500, -50., 50., 50, 0., 0.1);

  TGraph *gBkg1SigmaComb = new TGraph();
  gBkg1SigmaComb->SetNameTitle("PiMinusMuPlus/gBkg1SigmaComb", "Combinatorial background studies");
  TGraph *gBkg2SigmaComb = new TGraph();
  gBkg2SigmaComb->SetNameTitle("PiMinusMuPlus/gBkg2SigmaComb", "Combinatorial background studies");

  // Prompt - pi+mu-

  TH2D *hZTimePrompt = new TH2D("hZTimePrompt", "Muon-induced background studies", 500, 100., 190., 100, -15., 15.);

  TH1D *hInvMassPrompt = new TH1D("hInvMassPrompt", "Muon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptSR = new TH1D("hInvMassPromptSR", "Muon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptKolm = new TH1D("hInvMassPromptKolm", "Muon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  
  TH2D *hMom1VsMom2PeakPrompt = new TH2D("hMom1VsMom2PeakPrompt", "Muon-induced background studies", 100, 0., 200., 100, 0., 200.);
  TH2D *hMom1VsMom2NoPeakPrompt = new TH2D("hMom1VsMom2NoPeakPrompt", "Muon-induced background studies", 100, 0., 200., 100, 0., 200.);

  TH2D *hInvMassVsZPrompt = new TH2D("hInvMassVsZPrompt", "Muon-induced background studies", 300, 90., 120., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassVsZq2Prompt = new TH2D("hInvMassVsZq2Prompt", "Muon-induced background studies", 100, 90., 190., 50, 0.25, 1.);
  TH2D *hInvMassVsPq2Prompt = new TH2D("hInvMassVsPq2Prompt", "Muon-induced background studies", 100, 0., 300., 50, 0.25, 1.);
  TH2D *hInvMassVsDistq2Prompt = new TH2D("hInvMassVsDistq2Prompt", "Muon-induced background studies", 100, 0., 300., 50, 0.25, 1.);
  TH2D *hInvMassVsRICHq2Prompt = new TH2D("hInvMassVsRICHq2Prompt", "Muon-induced background studies", 6, -0.5, 5.5, 50, 0.25, 1.);
  TH2D *hSRPrompt = new TH2D("hSRPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPrompt = new TH2D("hSRFinalPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  
  TGraph *gBkg1SigmaPrompt = new TGraph();
  gBkg1SigmaPrompt->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPrompt", "Muon-induced background studies");
  TGraph *gBkg2SigmaPrompt = new TGraph();
  gBkg2SigmaPrompt->SetNameTitle("PiMinusMuPlus/gBkg2SigmaPrompt", "Muon-induced background studies");

  // Parasitic - pi+mu-

  TH1D *hRPar = new TH1D("hRPar", "Kaon-induced background studies", 1000, 0., 2000.);
  TH1D *hInvMassPar = new TH1D("hInvMassPar", "Kaon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParSR = new TH1D("hInvMassParSR", "Kaon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParKolm = new TH1D("hInvMassParKolm", "Kaon-induced background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassVsCDAFinalPar = new TH2D("hInvMassVsCDAPar", "Kaon-induced background studies", 100, 0., 1000., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassVsRFinalPar = new TH2D("hInvMassVsRFinalPar", "Kaon-induced background studies", 100, 0., 1000., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);

  TH2D *hSRPar = new TH2D("hSRPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPar = new TH2D("hSRFinalPar", "Signal region", 500, -50., 50., 50, 0., 0.1);

  TGraph *gBkg1SigmaPar = new TGraph();
  gBkg1SigmaPar->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPar", "Kaon-induced background studies");
  TGraph *gBkg2SigmaPar = new TGraph();
  gBkg2SigmaPar->SetNameTitle("PiMinusMuPlus/gBkg2SigmaPar", "Kaon-induced background studies");

  // SR - pi+mu-

  TH1D *hCDASR = new TH1D("hCDASR", "Signal region background studies", 100, 0., 1000.);
  TH1D *hTimeSR = new TH1D("hTimeSR", "Signal region background studies", 100, -15., 15.);
  TH1D *hZSR = new TH1D("hZSR", "Signal region background studies", 500, 100., 190.);
  TH1D *hZSRq2 = new TH1D("hZSRq2", "Signal region background studies", 500, 100., 190.);

  TH1D *hInvMassSR = new TH1D("hInvMassSR", "Signal region background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassVsCDASR = new TH2D("hInvMassVsCDASR", "Signal region background studies", 100, 0., 1000., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);

  TH2D *hSRSR = new TH2D("hSRSR", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // MC only - pi+mu-

  TH1D *hMass[172];
  for (Int_t i = 0; i < 172; i++) {
    hMass[i] = new TH1D(Form("hMass%d", i), Form("hMass%d", i), 100, (i*10.+245.)/1000., (i*10.+255.)/1000.);
  }
  TH2D *hInvMassRecoVsMC = new TH2D("hInvMassRecoVsMC", "MC studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TGraph *gMassMC = new TGraph();
  gMassMC->SetNameTitle("PiMinusMuPlus/gMassMC", "MC studies");

  // Scan tree
  
  Bool_t Zero = (!an.Contains("Pos") && !an.Contains("Neg"));
  Bool_t MC = (!histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data"));

  if (Zero && !MC)
    cout<<"**********Processing 0-charge events**********"<<endl;
  else if (an.Contains("Neg"))
    cout<<"**********Processing neg-charge events**********"<<endl;
  else if (Zero && MC)
    cout<<"**********Processing MC events**********"<<endl;

  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    // Booleans for bkg studies
    
    Bool_t outsideSR = ((ZCDALine < ZCDALineMin || ZCDALine > ZCDALineMax || (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine > CDALineMax)));
    Bool_t CDAIn = (CDA < CDAMax);
    Bool_t Comb = (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)));
    Bool_t Prompt = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t SB = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t FV = (Zvertex >= ZVertexMax && Zvertex <= ZEndFV);
    Bool_t Par = (BeamCDA1 <= 50. && BeamCDA2 <= 50.);
    Bool_t PiMinusMuPlus = ((Assoc == 2 && Charge1 == -1) || (Assoc == 1 && Charge1 == 1));
    Bool_t Spike1 = ((Assoc == 1 && Charge1 == -1) && (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) && (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax));
    Bool_t Spike2 = ((Assoc == 2 && Charge1 == 1) && (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) && (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax));
    Bool_t noSpike = (!Spike1 && !Spike2);

    // MC studies on reco mass vs true mass
    
    if (MC && Zero && noSpike && PiMinusMuPlus) {
      hInvMassRecoVsMC->Fill(trueMass, invMass/1000.);
      Int_t massIndex = std::round((trueMass*1000-250)/10);
      hMass[massIndex]->Fill(invMass/1000.);
    }
    
    // 1 - Events outside blinded region, to show SB for all bkgs
    
    if (outsideSR && Zero && noSpike && CDAIn && PiMinusMuPlus) { // all events outside blinded region
      hCDASR->Fill(BeamCDA1);
      hCDASR->Fill(BeamCDA2);
      hZSR->Fill(Zvertex/1000.);
      hTimeSR->Fill(CHODTime1-CHODTime2);
      hInvMassSR->Fill(invMass/1000.);
      hSRSR->Fill(ZCDALine/1000., CDALine/1000.);
      hInvMassVsCDASR->Fill(BeamCDA1, invMass/1000.);
      hInvMassVsCDASR->Fill(BeamCDA2, invMass/1000.);
    }
    
    // 2 - Combinatorial SB
    
    if (Comb && Zero && CDAIn && PiMinusMuPlus) { // all events inside time sidebands
      
     // Study spike at 75 GeV                                                                                       

      if (Assoc == 1 && Charge1 == 1) { // mu+pi-                                                                    
	hMomMuPosComb->Fill(Mom1->Mag()/1000.);
        hMomPiNegComb->Fill(Mom2->Mag()/1000.);
      }
      else if (Assoc == 2 && Charge1 == -1) { // pi-mu+                                                              
        hMomPiNegComb->Fill(Mom1->Mag()/1000.);
	hMomMuPosComb->Fill(Mom2->Mag()/1000.);
      }

      // Select momentum in [70-80] GeV to study spike                                                               

      if (Assoc == 1 && Charge1 == 1) { // mu+pi-                                                                    
        if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
          hGTK3XYMuPosComb->Fill(xGTK31, yGTK31);
        if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
          hGTK3XYPiNegComb->Fill(xGTK32, yGTK32);
      }
      else if (Assoc == 2 && Charge1 == -1) { // pi-mu+                                                              
        if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
          hGTK3XYPiNegComb->Fill(xGTK31, yGTK31);
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
      else if (Assoc == 2 && Charge1 == -1) { // pi-mu+                                                              
        if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
          if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
            hGTK3IDPiNegComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
        }
        if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
          if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)
            hGTK3IDMuPosComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
        }
      }
    
      if (noSpike) { // all events inside time sidebands and no spike
	hInvMassComb->Fill(invMass/1000.);            
	hSRComb->Fill(ZCDALine/1000., CDALine/1000.);
      }
    }
    
    // 3 - Prompt SB
    
    if (Prompt && Zero && noSpike && CDAIn && PiMinusMuPlus) { // all events inside vertex sidebands
      hZTimePrompt->Fill(Zvertex/1000., CHODTime1-CHODTime2);
      hSRPrompt->Fill(ZCDALine/1000., CDALine/1000.);
    }

    if (FV && an.Contains("Neg") && noSpike && CDAIn) { // prompt with q = 2 in FV
      hInvMassPrompt->Fill(invMass/1000.);
      hInvMassVsZq2Prompt->Fill(Zvertex/1000., invMass/1000.);
      hInvMassVsPq2Prompt->Fill(TotMom->Mag()/1000., invMass/1000.);
      hInvMassVsDistq2Prompt->Fill(BeamlineDist, invMass/1000.);
      
      if (RICHInfo1->at(0) == 99)
        RICHInfo1->at(0) = 5;
      if (RICHInfo2->at(0) == 99)
        RICHInfo2->at(0) = 5;
      
      hInvMassVsRICHq2Prompt->Fill(RICHInfo1->at(0), invMass/1000.);
      hInvMassVsRICHq2Prompt->Fill(RICHInfo2->at(0), invMass/1000.);

      if (invMass > 450. && invMass < 550.)
	hMom1VsMom2PeakPrompt->Fill(Mom2->Mag()/1000., Mom1->Mag()/1000.);
      else
	hMom1VsMom2NoPeakPrompt->Fill(Mom2->Mag()/1000., Mom1->Mag()/1000.);
    }
    
    // 4 - Parasitic SB
    
    if (Par && Zero && noSpike && CDAIn && PiMinusMuPlus) { // all events inside beam CDA sidebands
      hInvMassPar->Fill(invMass/1000.);
      hSRPar->Fill(ZCDALine/1000., CDALine/1000.);
      hRPar->Fill(R);
    }
    
    // 5 - Events inside blinded region
    
    if (!outsideSR && noSpike) { // all events inside blinded region...
      if (Comb && Zero && !MC && CDAIn && PiMinusMuPlus) { // ...and time sidebands (not enriched sample to study spike at 75 GeV)
	counterComb++;
	hSRFinalComb->Fill(ZCDALine/1000., CDALine/1000.);
	hInvMassCombSR->Fill(invMass/1000.);
      }
      if (Par && Zero && !MC && CDAIn && PiMinusMuPlus) { // ...and beam distance sidebands
	counterPar++;
	hSRFinalPar->Fill(ZCDALine/1000., CDALine/1000.);
	hInvMassParSR->Fill(invMass/1000.);
	hInvMassVsCDAFinalPar->Fill(BeamCDA1, invMass/1000.);
	hInvMassVsCDAFinalPar->Fill(BeamCDA2, invMass/1000.);
	hInvMassVsRFinalPar->Fill(R, invMass/1000.);
      }
      if (Prompt && Zero && !MC && CDAIn && PiMinusMuPlus) { // ...and vertex sidebands (0-charge)
	counterPromptSBZ++;
      }
      if (SB && an.Contains("Neg") && CDAIn) { // ...and vertex sidebands (pos-charge)
	counterPromptSBN++;
      }
      if (FV && an.Contains("Neg") && CDAIn) { // ...and in FV (pos-charge)
	counterPromptFVN++;
	hSRFinalPrompt->Fill(ZCDALine/1000., CDALine/1000.);
	hInvMassPromptSR->Fill(invMass/1000.);		
      }
    }
  }

  // Mass hypothesis scan

  // Sigma vs mass (MC studies)
    
  if (MC && Zero) { 
    Double_t sigma = 0.;
    Double_t sigmaMin = 999.;

    for (Int_t i = 0; i < 172; i++) {
      hMass[i]->Draw();
    }
    for (Int_t i = 0; i < 172; i++) {
      TF1 *f1 = new TF1("f1", "gaus", (i*10.+245.)/1000., (i*10.+255.)/1000.);
      hMass[i]->Fit("f1", "Rq");
      sigma = f1->GetParameter(2);
      gMassMC->SetPoint(i, (i*10.+250.)/1000., sigma);
      if (sigma < sigmaMin && sigma != 0.)
        sigmaMin = sigma;
    }

    minSigma = sigmaMin;
    gMassMC->Fit("pol5");
    TF1 *f = gMassMC->GetFunction("pol5");

    for (Int_t j = 0; j < 6; j++)
      Par[j] = f->GetParameter(j);
  }

  // Kolmogorov test and mass scan for data

  if (!MC) {   
    if (Zero) {
	
      // A - Combinatorial

      cout<<"Combinatorial bkg - Kolmogorov test for SR sample: "<<hInvMassCombSR->KolmogorovTest(hInvMassComb)<<endl;
      hInvMassCombKolm = (TH1D*)hInvMassComb->Clone();
      hInvMassCombKolm->SetName("hInvMassCombKolm");
      
      hInvMassCombKolm->Scale(hInvMassCombSR->Integral()/hInvMassCombKolm->Integral());
      gBkg1SigmaComb = WindowScanNoFixedSigma(hInvMassCombKolm, 1.);
      gBkg2SigmaComb = WindowScanNoFixedSigma(hInvMassCombKolm, 2.);
      gBkg1SigmaComb->SetNameTitle("PiMinusMuPlus/gBkg1SigmaComb", "Expected combinatorial background events");
      gBkg2SigmaComb->SetNameTitle("PiMinusMuPlus/gBkg1SigmaComb", "Expected combinatorial background events");

      // B - Parasitic

      cout<<"Parasitic bkg - Kolmogorov test for SR sample: "<<hInvMassParSR->KolmogorovTest(hInvMassPar)<<endl;
      hInvMassParKolm = (TH1D*)hInvMassPar->Clone();
      hInvMassParKolm->SetName("hInvMassParKolm");
      hInvMassParKolm->Scale(hInvMassParSR->Integral()/hInvMassParKolm->Integral());
      TF1 *f = new TF1("f", "gaus", 0.25, 0.3);
      hInvMassParKolm->Fit("f", "Rq");
      cout<<"Fit double peak in parasitic: "<<f->GetParameter(1)<<" +- "<<f->GetParameter(2)<<endl;

      gBkg1SigmaPar = WindowScanNoFixedSigma(hInvMassParKolm, 1.);
      gBkg2SigmaPar = WindowScanNoFixedSigma(hInvMassParKolm, 2.);

      gBkg1SigmaPar->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPar", "Expected kaon-induced background events");
      gBkg2SigmaPar->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPar", "Expected kaon-induced background events");

      SumGraphs(gBkg1SigmaComb, gBkg1SigmaPar, *gBkg1SigmaBuffer);
      SumGraphs(gBkg2SigmaComb, gBkg2SigmaPar, *gBkg2SigmaBuffer);

      *gBkg1SigmaTot = *gBkg1SigmaBuffer;
      *gBkg2SigmaTot = *gBkg2SigmaBuffer;
    }
    /*
    if (an.Contains("Neg")) {

      // C - Prompt

      hInvMassPromptSR->Scale(counterPromptSBZ/counterPromptSBN);
      hInvMassPrompt->Scale(counterPromptSBZ/counterPromptSBN);
      cout<<"Prompt bkg - Kolmogorov test for SR sample: "<<hInvMassPromptSR->KolmogorovTest(hInvMassPrompt)<<endl;
      hInvMassPromptKolm = (TH1D*)hInvMassPrompt->Clone();
      hInvMassPromptKolm->SetName("hInvMassPromptKolm");
      hInvMassPromptKolm->Scale(hInvMassPromptSR->Integral()/hInvMassPromptKolm->Integral());
      gBkg1SigmaPrompt = WindowScanNoFixedSigma(hInvMassPromptKolm, 1.);
      gBkg2SigmaPrompt = WindowScanNoFixedSigma(hInvMassPromptKolm, 2.);
      gBkg1SigmaPrompt->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPrompt", "Expected muon-induced background events");
      gBkg2SigmaPrompt->SetNameTitle("PiMinusMuPlus/gBkg1SigmaPrompt", "Expected muon-induced background events");

      SumGraphs(gBkg1SigmaPrompt, gBkg1SigmaBuffer, *gBkg1SigmaTot);
      SumGraphs(gBkg2SigmaPrompt, gBkg2SigmaBuffer, *gBkg2SigmaTot);
      
      for (Int_t j = 0; j < 6; j++) 
	Par[j] = 0.;
      
      minSigma = 0.;
    }
    */
  }

  // Saving histograms

  Save(path + "PiMinusMuPlus/SR/", c, hCDASR, "Track-beamline CDA [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hTimeSR, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hZSR, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hZSRq2, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hInvMassSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hSRSR, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/SR/", c, hInvMassVsCDASR, "Track-beamline CDA [mm]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);

  Save(path + "PiMinusMuPlus/Comb/", c, hInvMassComb, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hInvMassCombKolm, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hInvMassCombSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hMomPiNegComb, "Track momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hMomMuPosComb, "Track momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hGTK3IDPiNegComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hGTK3IDMuPosComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hGTK3XYPiNegComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hGTK3XYMuPosComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hSRComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, hSRFinalComb, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, gBkg1SigmaComb, "gBkg1SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Comb/", c, gBkg2SigmaComb, "gBkg2SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassPrompt, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassPromptSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassPromptKolm, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hMom1VsMom2PeakPrompt, "Momentum of track1 [GeV/c]", "Momentum of track2 [GeV/c]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hMom1VsMom2NoPeakPrompt, "Momentum of track1 [GeV/c]", "Momentum of track2 [GeV/c]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hZTimePrompt, "Z coordinate of vertex [m]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassVsZPrompt, "Z coordinate of vertex [m]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassVsZq2Prompt, "Z coordinate of vertex [m]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassVsPq2Prompt, "Modulus of mother momentum [GeV/c]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassVsDistq2Prompt, "Vertex-beam distance [mm]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hInvMassVsRICHq2Prompt, "RICH hypothesis", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hSRPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, hSRFinalPrompt, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, gBkg1SigmaPrompt, "gBkg1SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Prompt/", c, gBkg2SigmaPrompt, "gBkg2SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  Save(path + "PiMinusMuPlus/Par/", c, hRPar, "Two-track distance at CH1 [mm]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hInvMassPar, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hInvMassParSR, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hInvMassParKolm, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hSRPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hSRFinalPar, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hInvMassVsCDAFinalPar, "Track-beamline CDA [mm]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, hInvMassVsRFinalPar, "Two-track distance at CH1 [mm]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, gBkg1SigmaPar, "gBkg1SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiMinusMuPlus/Par/", c, gBkg2SigmaPar, "gBkg2SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  if (!an.Contains("Pos")  && !histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data")) {
    Save(path + "PiMinusMuPlus/MC/", c, hInvMassRecoVsMC, "MC true invariant mass [GeV/c^{2}]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
    Save(path + "PiMinusMuPlus/MC/", c, gMassMC, "gMassMC", "Reconstructed invariant mass [GeV/c^{2}]", "Mass resolution [GeV/c^{2}]", labelSize, titleSize);
  }

  tree->ResetBranchAddresses();
    
  return;
}

void BkgPiMinusMuPlus(TString dir, TString histo1, TString histo2) {

  // dir = output dir, histo1 = Data, histo2 = MC

  gErrorIgnoreLevel = kFatal;

  // Counters for pos/zero charge studies
  
  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t counterComb = 0;
  Double_t counterPar = 0;
  Double_t counterPromptSBZ = 0;
  Double_t counterPromptFVZ = 0;
  Double_t counterPromptSBN = 0;
  Double_t counterPromptFVN = 0;
  Double_t N = 0.;
  TString path = "/Users/lorenza/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/Zero/PiMinusMuPlus/";

  // Total expected bkg
  
  TGraph *gBkg1SigmaBuffer = new TGraph();
  gBkg1SigmaBuffer->SetNameTitle("gBkg1SigmaBuffer", "Total expected background");
  TGraph *gBkg2SigmaBuffer = new TGraph();
  gBkg2SigmaBuffer->SetNameTitle("gBkg2SigmaBuffer", "Total expected background");
  TGraph *gBkg1SigmaTot = new TGraph();
  gBkg1SigmaTot->SetNameTitle("gBkg1SigmaTot", "Total expected background");
  TGraph *gBkg2SigmaTot = new TGraph();
  gBkg2SigmaTot->SetNameTitle("gBkg2SigmaTot", "Total expected background");

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  Analyzer(dir, histo2, "HeavyNeutrino", c, counterComb, counterPar, counterPromptSBZ, counterPromptFVZ, counterPromptSBN, counterPromptFVN, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot);
  Analyzer(dir, histo1, "HeavyNeutrino", c, counterComb, counterPar, counterPromptSBZ, counterPromptFVZ, counterPromptSBN, counterPromptFVN, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot);
  Analyzer(dir, histo1, "HeavyNeutrinoNeg", c, counterComb, counterPar, counterPromptSBZ, counterPromptFVZ, counterPromptSBN, counterPromptFVN, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot);
  
  Save(path + "Total/", c, gBkg1SigmaTot, "gBkg1SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "Total/", c, gBkg2SigmaTot, "gBkg2SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  if (counterPromptSBN != 0.)
    N = counterPromptFVN*counterPromptSBZ/counterPromptSBN;
  
  counterPromptFVZ = N;
  
  cout<<"PI+MU-: Number of events in SR and time sidebands: "<<counterComb<<", and beamdist sidebands: "<<counterPar<<", and Z sidebands (0-charge): "<<counterPromptSBZ<<", and Z sidebands (neg-charge): "<<counterPromptSBN<<", and FV (neg-charge): "<<counterPromptFVN<<", and FV (0-charge): "<<counterPromptFVZ<<endl;

  _exit(0);
}
