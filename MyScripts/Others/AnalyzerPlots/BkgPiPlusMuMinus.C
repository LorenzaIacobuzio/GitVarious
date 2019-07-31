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
  /*
  if (name.Contains("ZPrompt"))
    gPad->SetLogy();
  else
    gPad->SetLogy(0);
  */
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
  g->GetYaxis()->SetTitleOffset(1.5);
  g->GetXaxis()->SetTitle(x);
  g->GetYaxis()->SetTitle(y);
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
    //cout<<massHyp<<" "<<massMin-minSigma<<" "<<massMax+minSigma<<" "<<minSigma<<" "<<sigma<<endl;
  }

  return res;
}
/*

// Compute CDA between two lines

void ComputeCDA(TVector3 p1, TVector3 p2, TVector3 dir1, TVector3 dir2, TVector3 &vertex, Double_t &CDA) {

Double_t A = dir1*dir1;
Double_t B = dir2*dir2;
Double_t C = dir1*dir2;
Double_t det = C*C - A*B;
if (det == 0.) { // the two lines are parallel
vertex.SetXYZ(0., 0., 0.);
CDA = -999;
}
else {
TVector3 R12 = p1 - p2;
double D = R12*dir1;
double E = R12*dir2;
double T1 = (B*D - C*E)/det;
double T2 = (C*D - A*E)/det;
TVector3 Q1 = p1 + T1*dir1;
TVector3 Q2 = p2 + T2*dir2;
vertex = 0.5*(Q1+Q2);
CDA = (Q1-Q2).Mag();
}
}
*/

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

void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, Double_t &counterCombPiPlusMuMinus, Double_t &counterParPiPlusMuMinus, Double_t &counterPromptSBZPiPlusMuMinus, Double_t &counterPromptFVZPiPlusMuMinus, Double_t &counterPromptSBPPiPlusMuMinus, Double_t &counterPromptFVPPiPlusMuMinus, TGraph* gBkg1SigmaBufferPiPlusMuMinus, TGraph* gBkg2SigmaBufferPiPlusMuMinus, TGraph* gBkg1SigmaTotPiPlusMuMinus, TGraph* gBkg2SigmaTotPiPlusMuMinus) {

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
  
  if (an.Contains("Pos"))
    analyzer = "Pos";
  else if (!an.Contains("Pos") && !an.Contains("Neg"))
    analyzer = "Zero";
  
  if (dir != "")
    path = dir;
  else {
    if (histo1.Contains("Data"))
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/" + analyzer + "/";
    else
      path = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/MC/" + analyzer + "/";
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
  Int_t Charge1;
  Int_t Charge2;
  Int_t nSec;
  Int_t nCHOD;
  Int_t nNewCHOD;
  Int_t nLKr;
  std::vector<Double_t> *RICHInfo1 = new std::vector<Double_t>;
  std::vector<Double_t> *RICHInfo2 = new std::vector<Double_t>;

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
  tree->SetBranchAddress("Charge1", &Charge1);
  tree->SetBranchAddress("Charge2", &Charge2);
  tree->SetBranchAddress("nSec", &nSec);
  tree->SetBranchAddress("nCHOD", &nCHOD);
  tree->SetBranchAddress("nNewCHOD", &nNewCHOD);
  tree->SetBranchAddress("nLKr", &nLKr);
  tree->SetBranchAddress("RICHInfo1", &RICHInfo1);
  tree->SetBranchAddress("RICHInfo2", &RICHInfo2);

  // Setting histograms - splitting analysis for pi+mu- and pi-mu+ cases
  
  // Combinatorial - pi+mu-
    
  TH1D *hTimeCombPiPlusMuMinus = new TH1D("hTimeComb", "Combinatorial background studies", 100, -15., 15.);

  TH1D *hInvMassCombPiPlusMuMinus = new TH1D("hInvMassComb", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombKolmPiPlusMuMinus = new TH1D("hInvMassCombKolm", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSRPiPlusMuMinus = new TH1D("hInvMassCombSR", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSREnrichedPiPlusMuMinus = new TH1D("hInvMassCombSREnriched", "Enriched combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  
  TH1D *hMomPiPosCombPiPlusMuMinus = new TH1D("hMomPiPosComb", "Combinatorial background studies", 100, 0., 200.);
  TH1D *hMomMuNegCombPiPlusMuMinus = new TH1D("hMomMuNegComb", "Combinatorial background studies", 100, 0., 200.);

  TH1D *hInvMassCombExtra1PiPlusMuMinus = new TH1D("hInvMassCombExtra1", "Combinatorial background studies - K^{+}#rightarrow #pi^{+} + #mu^{-} halo", 1000, -1., 2.); // K->pi + mu halo (peak at m_pi0^2)
  TH1D *hInvMassCombExtra2PiPlusMuMinus = new TH1D("hInvMassCombExtra2", "Combinatorial background studies - K^{+}#rightarrow #mu^{+} + beam #pi^{-}", 1000, -1., 2.); // K->mu + beam pi (peak at 0)
  TH1D *hInvMassCombExtra3PiPlusMuMinus = new TH1D("hInvMassCombExtra3", "Combinatorial background studies - K^{+}#rightarrow #pi^{+}#pi^{+}#pi^{-} (lost+decayed)", 1000, -1., 4.); // K->3pi (pi lost +  pi->mu) (peak at m_pi^2)
  TH1D *hInvMassCombExtra4PiPlusMuMinus = new TH1D("hInvMassCombExtra4", "Combinatorial background studies - #pi^{+}#rightarrow #mu^{+} + beam #pi^{-}", 1000, -1., 2.); // pi->mu + beam pi (peak at 0)

  TH2D *hGTK3XYPiPosCombPiPlusMuMinus = new TH2D("hGTK3XYPiPosComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYMuNegCombPiPlusMuMinus = new TH2D("hGTK3XYMuNegComb", "Combinatorial background studies", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3IDPiPosCombPiPlusMuMinus = new TH2D("hGTK3IDPiPosComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDMuNegCombPiPlusMuMinus = new TH2D("hGTK3IDMuNegComb", "Combinatorial background studies", 6, -0.5, 5.5, 100, 0., 250.);

  TH2D *hSRCombPiPlusMuMinus = new TH2D("hSRComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalCombPiPlusMuMinus = new TH2D("hSRFinalComb", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalCombEnrichedPiPlusMuMinus = new TH2D("hSRFinalCombEnriched", "Enriched signal region", 500, -50., 50., 50, 0., 0.1);

  TGraph *gBkg1SigmaCombPiPlusMuMinus = new TGraph();
  gBkg1SigmaCombPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg1SigmaComb", "Combinatorial background studies");
  TGraph *gBkg2SigmaCombPiPlusMuMinus = new TGraph();
  gBkg2SigmaCombPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg2SigmaComb", "Combinatorial background studies");

  // Prompt - pi+mu-

  TH1D *hZPromptPiPlusMuMinus = new TH1D("hZPrompt", "Prompt background studies", 500, 100., 190.);
  TH1D *hZPromptPPiPlusMuMinus = new TH1D("hZPromptP", "Prompt background studies", 500, 100., 190.);
  TH2D *hZTimePromptPiPlusMuMinus = new TH2D("hZTimePrompt", "Prompt background studies", 500, 100., 190., 100, -15., 15.);
  TH2D *hZTimePromptPPiPlusMuMinus = new TH2D("hZTimePromptP", "Prompt background studies", 500, 100., 190., 100, -15., 15.);

  TH1D *hInvMassPromptPiPlusMuMinus = new TH1D("hInvMassPrompt", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptSRPiPlusMuMinus = new TH1D("hInvMassPromptSR", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptKolmPiPlusMuMinus = new TH1D("hInvMassPromptKolm", "Prompt background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  
  TH1D *hMomPiPromptPiPlusMuMinus = new TH1D("hMomPiPrompt", "Prompt background studies", 100, 0., 200.);
  TH1D *hMomMuPromptPiPlusMuMinus = new TH1D("hMomMuPrompt", "Prompt background studies", 100, 0., 200.);

  TH2D *hInvMassVsZPromptPiPlusMuMinus = new TH2D("hInvMassVsZPrompt", "Prompt background studies", 300, 90., 120., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassVsZq2PromptPiPlusMuMinus = new TH2D("hInvMassVsZq2Prompt", "Prompt background studies", 100, 90., 190., 50, 0.25, 1.);
  TH2D *hInvMassVsPq2PromptPiPlusMuMinus = new TH2D("hInvMassVsPq2Prompt", "Prompt background studies", 100, 0., 300., 50, 0.25, 1.);
  TH2D *hInvMassVsDistq2PromptPiPlusMuMinus = new TH2D("hInvMassVsDistq2Prompt", "Prompt background studies", 100, 0., 300., 50, 0.25, 1.);
  TH2D *hSRPromptPiPlusMuMinus = new TH2D("hSRPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalPromptPiPlusMuMinus = new TH2D("hSRFinalPrompt", "Signal region", 500, -50., 50., 50, 0., 0.1);
  
  TGraph *gBkg1SigmaPromptPiPlusMuMinus = new TGraph();
  gBkg1SigmaPromptPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPrompt", "Prompt background studies");
  TGraph *gBkg2SigmaPromptPiPlusMuMinus = new TGraph();
  gBkg2SigmaPromptPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg2SigmaPrompt", "Prompt background studies");

  // Parasitic - pi+mu-

  TH1D *hDistParPiPlusMuMinus = new TH1D("hDistPar", "Parasitic background studies", 100, 0., 1000.);
  TH2D *hDistvsMassParPiPlusMuMinus = new TH2D("hDistvsMassPar", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);
  
  TH1D *hInvMassParPiPlusMuMinus = new TH1D("hInvMassPar", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParSRPiPlusMuMinus = new TH1D("hInvMassParSR", "Combinatorial background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParKolmPiPlusMuMinus = new TH1D("hInvMassParKolm", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);

  TH1D *hMomPiParPiPlusMuMinus = new TH1D("hMomPiPar", "Parasitic background studies", 100, 0., 200.);
  TH1D *hMomMuParPiPlusMuMinus = new TH1D("hMomMuPar", "Parasitic background studies", 100, 0., 200.);

  TH2D *hSRParPiPlusMuMinus = new TH2D("hSRPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hSRFinalParPiPlusMuMinus = new TH2D("hSRFinalPar", "Signal region", 500, -50., 50., 50, 0., 0.1);
  TH2D *hDistvsMassParSRPiPlusMuMinus = new TH2D("hDistvsMassParSR", "Parasitic background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);

  TGraph *gBkg1SigmaParPiPlusMuMinus = new TGraph();
  gBkg1SigmaParPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPar", "Parasitic background studies");
  TGraph *gBkg2SigmaParPiPlusMuMinus = new TGraph();
  gBkg2SigmaParPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gBkg2SigmaPar", "Parasitic background studies");

  // SR - pi+mu-

  TH1D *hDistSRPiPlusMuMinus = new TH1D("hDistSR", "Signal region background studies", 100, 0., 1000.);
  TH1D *hCDASRPiPlusMuMinus = new TH1D("hCDASR", "Signal region background studies", 100, 0., 1000.);
  TH1D *hTimeSRPiPlusMuMinus = new TH1D("hTimeSR", "Signal region background studies", 100, -15., 15.);
  TH1D *hZSRPiPlusMuMinus = new TH1D("hZSR", "Signal region background studies", 500, 100., 190.);
  TH2D *hDistvsMassSRPiPlusMuMinus = new TH2D("hDistvsMassSR", "Signal region background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., 50, 0., 1000.);

  TH1D *hInvMassSRPiPlusMuMinus = new TH1D("hInvMassSR", "Signal region background studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);

  TH1D *hMomPiSRPiPlusMuMinus = new TH1D("hMomPiSR", "Signal region background studies", 100, 0., 200.);
  TH1D *hMomMuSRPiPlusMuMinus = new TH1D("hMomMuSR", "Signal region background studies", 100, 0., 200.);

  TH2D *hSRSRPiPlusMuMinus = new TH2D("hSRSR", "Signal region", 500, -50., 50., 50, 0., 0.1);

  // MC only - pi+mu-

  TH1D *hInvMassMCPiPlusMuMinus = new TH1D("hInvMassMC", "MC studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH2D *hInvMassRecoVsMCPiPlusMuMinus = new TH2D("hInvMassRecoVsMCPiPlusMuMinus", "MC studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TGraph *gMassMCPiPlusMuMinus = new TGraph();
  gMassMCPiPlusMuMinus->SetNameTitle("PiPlusMuMinus/gMassMC", "MC studies");

  // Scan tree
  
  Bool_t Zero = (!an.Contains("Pos") && !an.Contains("Neg"));
  Bool_t MC = (!histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data"));

  if (Zero && !MC)
    cout<<"**********Processing 0-charge events**********"<<endl;
  else if (an.Contains("Pos"))
    cout<<"**********Processing pos-charge events**********"<<endl;
  else if (Zero && MC)
    cout<<"**********Processing MC events**********"<<endl;

  for(Int_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);

    // Booleans for bkg studies
    
    Bool_t outsideSR = ((ZCDALine < ZCDALineMin || ZCDALine > ZCDALineMax || (ZCDALine >= ZCDALineMin && ZCDALine <= ZCDALineMax && CDALine > CDALineMax)));
    Bool_t CDAIn = (CDA < CDAMax);
    Bool_t CDAInEnriched = (CDA < 100.);
    Bool_t Comb = (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)));
    Bool_t Prompt = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t SB = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t FV = (Zvertex >= ZVertexMax && Zvertex <= ZEndFV);
    //Bool_t Par = (BeamlineDist >= BeamdistMin && BeamlineDist <= BeamdistMax);
    Bool_t Par = (BeamCDA1 <= 50. && BeamCDA2 <= 50.);
    Bool_t PiPlusMuMinus = ((Assoc == 2 && Charge1 == 1) || (Assoc == 1 && Charge1 == -1));
    Bool_t Spike1 = ((Assoc == 1 && Charge1 == -1) && (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) && (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax));
    Bool_t Spike2 = ((Assoc == 2 && Charge1 == 1) && (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) && (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax));
    Bool_t noSpike = (!Spike1 && !Spike2);

    // MC studies on reco mass vs true mass
    
    if (MC && Zero && noSpike && PiPlusMuMinus) {
      hInvMassMCPiPlusMuMinus->Fill(invMass/1000.);
      hInvMassRecoVsMCPiPlusMuMinus->Fill(trueMass, invMass/1000.);
    }
    
    // 1 - Events outside blinded region, to show SB for all bkgs
    
    if (outsideSR && Zero && noSpike && CDAIn && PiPlusMuMinus) { // all events outside blinded region
      hDistSRPiPlusMuMinus->Fill(BeamlineDist, Weight);
      hCDASRPiPlusMuMinus->Fill(BeamCDA1);//, Weight);
      hCDASRPiPlusMuMinus->Fill(BeamCDA2);//, Weight);
      hTimeSRPiPlusMuMinus->Fill(CHODTime1-CHODTime2, Weight);
      hZSRPiPlusMuMinus->Fill(Zvertex/1000., Weight);
      hInvMassSRPiPlusMuMinus->Fill(invMass/1000., Weight);
      hMomPiSRPiPlusMuMinus->Fill(Mom1->Mag()/1000., Weight);
      hMomMuSRPiPlusMuMinus->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassSRPiPlusMuMinus->Fill(invMass/1000., BeamlineDist, Weight);
      hSRSRPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    // 2 - Combinatorial SB
    
    if (Comb && Zero && CDAIn && PiPlusMuMinus) { // all events inside time sidebands
      
      // Study spike at 75 GeV
      
      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	hMomMuNegCombPiPlusMuMinus->Fill(Mom1->Mag()/1000., Weight);
	hMomPiPosCombPiPlusMuMinus->Fill(Mom2->Mag()/1000., Weight);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	hMomPiPosCombPiPlusMuMinus->Fill(Mom1->Mag()/1000., Weight);
	hMomMuNegCombPiPlusMuMinus->Fill(Mom2->Mag()/1000., Weight);
      }
      
      // Select momentum in [70-80] GeV to study spike
      
      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYMuNegCombPiPlusMuMinus->Fill(xGTK31, yGTK31, Weight);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYPiPosCombPiPlusMuMinus->Fill(xGTK32, yGTK32, Weight);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYPiPosCombPiPlusMuMinus->Fill(xGTK31, yGTK31, Weight);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYMuNegCombPiPlusMuMinus->Fill(xGTK32, yGTK32, Weight);
      }

      // Select momentum in [70-80] GeV and xy @ GTK3 to study spike
      
      if (RICHInfo1->at(0) == 99)
	RICHInfo1->at(0) = 5;

      if (RICHInfo2->at(0) == 99)
	RICHInfo2->at(0) = 5;

      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDMuNegCombPiPlusMuMinus->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDPiPosCombPiPlusMuMinus->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDPiPosCombPiPlusMuMinus->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDMuNegCombPiPlusMuMinus->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
    
      if (noSpike) { // all events inside time sidebands and no spike
	hTimeCombPiPlusMuMinus->Fill(CHODTime1-CHODTime2, Weight);
	hInvMassCombPiPlusMuMinus->Fill(invMass/1000., Weight);            
	hSRCombPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
      }
    
      // Extra combinatorial bkg
      
      if (outsideSR && noSpike && PiPlusMuMinus) { // all events outside blinded region
	/*
	  TVector3 vertex;
	  Double_t CDAcomp;
	  TVector3 momK(0., 0., 1.2E-3);
	  TVector3 posK(0., 0., 102000.);
	  //TVector3 mom1(Mom1.Px(), Mom1.Py(), Mom1.Pz());
	  //TVector3 mom2(Mom2.Px(), Mom2.Py(), Mom2.Pz());
	  //TVector3 pos1(Pos1.x(), Pos1.y(), Pos1.z());
	  //TVector3 pos2(Pos2.x(), Pos2.y(), Pos2.z());
	  ComputeCDA(*Pos1, posK, *Mom1, momK, vertex, CDAcomp);
	*/
	Double_t muMass = 105.66;
	Double_t piMass = 139.57;
	Double_t kMass = 493.68;
	Double_t kMom = 75000.;
	TVector3 *threeMomK = new TVector3(0., 0., 75000.);
	Double_t kE = TMath::Sqrt(kMom*kMom + kMass*kMass);
	Double_t piE = TMath::Sqrt(kMom*kMom + piMass*piMass);
	Double_t hypMuE1 = TMath::Sqrt(Mom1->Px()*Mom1->Px() + Mom1->Py()*Mom1->Py() + Mom1->Pz()*Mom1->Pz() + muMass*muMass);
	Double_t hypMuE2 = TMath::Sqrt(Mom2->Px()*Mom2->Px() + Mom2->Py()*Mom2->Py() + Mom2->Pz()*Mom2->Pz() + muMass*muMass);
	Double_t hypPiE1 = TMath::Sqrt(Mom1->Px()*Mom1->Px() + Mom1->Py()*Mom1->Py() + Mom1->Pz()*Mom1->Pz() + piMass*piMass);
	Double_t hypPiE2 = TMath::Sqrt(Mom2->Px()*Mom2->Px() + Mom2->Py()*Mom2->Py() + Mom2->Pz()*Mom2->Pz() + piMass*piMass);

	if (Charge1 == 1.)
	  hInvMassCombExtra1PiPlusMuMinus->Fill(((kE + hypMuE1)*(kE + hypMuE1) - (*threeMomK + *Mom1).Mag2())/1.E6);
	else if (Charge2 == 1.)
	  hInvMassCombExtra1PiPlusMuMinus->Fill(((kE + hypMuE2)*(kE + hypMuE2) - (*threeMomK + *Mom2).Mag2())/1.E6);

	if (Charge1 == 1.)
	  hInvMassCombExtra2PiPlusMuMinus->Fill(((kE + hypPiE1)*(kE + hypPiE1) - (*threeMomK + *Mom1).Mag2())/1.E6);
	else if	(Charge2 == 1.)
	  hInvMassCombExtra2PiPlusMuMinus->Fill(((kE + hypPiE2)*(kE + hypPiE2) - (*threeMomK + *Mom2).Mag2())/1.E6);

	if (Charge1 == 1.)
	  hInvMassCombExtra4PiPlusMuMinus->Fill(((piE + hypMuE1)*(piE + hypMuE1) - (*threeMomK + *Mom1).Mag2())/1.E6);
	else if (Charge2 == 1.)
	  hInvMassCombExtra4PiPlusMuMinus->Fill(((piE + hypMuE2)*(piE + hypMuE2) - (*threeMomK + *Mom2).Mag2())/1.E6);
	
	hInvMassCombExtra3PiPlusMuMinus->Fill(((kE + hypPiE1 + hypPiE2)*(kE + hypPiE1 + hypPiE2) - (*threeMomK + *Mom1 + *Mom2).Mag2())/1.E6);
      }
    }
    
    // 3 - Prompt SB
    
    if (Prompt && Zero && noSpike && CDAIn && PiPlusMuMinus) { // all events inside vertex sidebands
      hZPromptPiPlusMuMinus->Fill(Zvertex/1000., Weight);
      hZTimePromptPiPlusMuMinus->Fill(Zvertex/1000., CHODTime1-CHODTime2, Weight);
      hMomPiPromptPiPlusMuMinus->Fill(Mom1->Mag()/1000., Weight);
      hMomMuPromptPiPlusMuMinus->Fill(Mom2->Mag()/1000., Weight);
      hSRPromptPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }

    if (FV && an.Contains("Pos") && noSpike && CDAIn && PiPlusMuMinus) { // prompt with q = 2 in FV
      hInvMassPromptPiPlusMuMinus->Fill(invMass/1000., Weight);
      hInvMassVsZq2PromptPiPlusMuMinus->Fill(Zvertex/1000., invMass/1000.);
      hInvMassVsPq2PromptPiPlusMuMinus->Fill(TotMom->Mag()/1000., invMass/1000.);
      hInvMassVsDistq2PromptPiPlusMuMinus->Fill(BeamlineDist, invMass/1000.);
    }

    if (Prompt && an.Contains("Pos") && noSpike && CDAIn && PiPlusMuMinus) { // prompt with q = 2 in SB
      hZPromptPPiPlusMuMinus->Fill(Zvertex/1000., Weight);
      hZTimePromptPPiPlusMuMinus->Fill(Zvertex/1000., CHODTime1-CHODTime2, Weight);
    }
    
    // 4 - Parasitic SB
    
    if (Par && Zero && noSpike && CDAIn && PiPlusMuMinus) { // all events inside beam distance sidebands
      hDistParPiPlusMuMinus->Fill(BeamlineDist, Weight);
      hInvMassParPiPlusMuMinus->Fill(invMass/1000., Weight);
      hMomPiParPiPlusMuMinus->Fill(Mom1->Mag()/1000., Weight);
      hMomMuParPiPlusMuMinus->Fill(Mom2->Mag()/1000., Weight);
      hDistvsMassParPiPlusMuMinus->Fill(invMass/1000., BeamlineDist, Weight);
      hSRParPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
    }
    
    // 5 - Events inside blinded region
    
    if (!outsideSR && noSpike && PiPlusMuMinus) { // all events inside blinded region...
      
      if (Comb && CDAInEnriched && Zero && !MC) { // ...and time sidebands (enriched sample)
	hSRFinalCombEnrichedPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSREnrichedPiPlusMuMinus->Fill(invMass/1000., Weight);
      }
      if (Comb && Zero && !MC && CDAIn) { // ...and time sidebands (not enriched sample to study spike at 75 GeV)
	counterCombPiPlusMuMinus++;
	hSRFinalCombPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassCombSRPiPlusMuMinus->Fill(invMass/1000., Weight);
      }
      if (Par && Zero && !MC && CDAIn) { // ...and beam distance sidebands
	hSRFinalParPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
	counterParPiPlusMuMinus++;
	hInvMassParSRPiPlusMuMinus->Fill(invMass/1000., Weight);
	hDistvsMassParSRPiPlusMuMinus->Fill(invMass/1000., BeamlineDist, Weight);
      }
      if (Prompt && Zero && !MC && CDAIn) { // ...and vertex sidebands (0-charge)
	counterPromptSBZPiPlusMuMinus++;
      }
      if (SB && an.Contains("Pos") && CDAIn) { // ...and vertex sidebands (pos-charge)
	counterPromptSBPPiPlusMuMinus++;
      }
      if (FV && an.Contains("Pos") && CDAIn) { // ...and in FV (pos-charge)
	hSRFinalPromptPiPlusMuMinus->Fill(ZCDALine/1000., CDALine/1000., Weight);
	hInvMassPromptSRPiPlusMuMinus->Fill(invMass/1000., Weight);		
	counterPromptFVPPiPlusMuMinus++;
      }
    }
  }

  // Mass hypothesis scan

  // Sigma vs mass (MC studies)
    
  if (MC && Zero) { 
    Int_t counter = 0;
    Double_t sigma = 0.;
    Double_t sigmaMin = 999.;

    hInvMassMCPiPlusMuMinus->Draw();

    for (Double_t mass = massMin; mass < massMax; mass += massStep) {
      TF1 *f1 = new TF1("f1", "gaus", mass-massStep, mass+massStep);
      hInvMassMCPiPlusMuMinus->Fit("f1", "Rq");
      sigma = f1->GetParameter(2);
      if (sigma < 0.001 || sigma > 0.01)
	gMassMCPiPlusMuMinus->SetPoint(counter, mass, 0.0055);
      else
	gMassMCPiPlusMuMinus->SetPoint(counter, mass, sigma);
      if (sigma < sigmaMin && sigma != 0.)
	sigmaMin = sigma;
      counter++;
    }

    minSigma = sigmaMin;
    gMassMCPiPlusMuMinus->Fit("pol5");
    TF1 *f = gMassMCPiPlusMuMinus->GetFunction("pol5");

    for (Int_t j = 0; j < 6; j++)
      Par[j] = f->GetParameter(j);
  }

  // Kolmogorov test and mass scan for data

  if (!MC) {   
    if (Zero) {
	
      // A - Combinatorial

      cout<<"Combinatorial bkg - Kolmogorov test for SR sample: "<<hInvMassCombSRPiPlusMuMinus->KolmogorovTest(hInvMassCombPiPlusMuMinus)<<" and enriched SR sample: "<<hInvMassCombSREnrichedPiPlusMuMinus->KolmogorovTest(hInvMassCombPiPlusMuMinus)<<endl;
      hInvMassCombKolmPiPlusMuMinus = (TH1D*)hInvMassCombPiPlusMuMinus->Clone();
      hInvMassCombKolmPiPlusMuMinus->SetName("hInvMassCombKolm");
      
      hInvMassCombKolmPiPlusMuMinus->Scale(hInvMassCombSRPiPlusMuMinus->Integral()/hInvMassCombKolmPiPlusMuMinus->Integral());
      gBkg1SigmaCombPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassCombKolmPiPlusMuMinus, 1.);
      gBkg2SigmaCombPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassCombKolmPiPlusMuMinus, 2.);
      
      // B - Parasitic

      cout<<"Parasitic bkg - Kolmogorov test for SR sample: "<<hInvMassParSRPiPlusMuMinus->KolmogorovTest(hInvMassParPiPlusMuMinus)<<endl;
      hInvMassParKolmPiPlusMuMinus = (TH1D*)hInvMassParPiPlusMuMinus->Clone();
      hInvMassParKolmPiPlusMuMinus->SetName("hInvMassParKolm");
      hInvMassParKolmPiPlusMuMinus->Scale(hInvMassParSRPiPlusMuMinus->Integral()/hInvMassParKolmPiPlusMuMinus->Integral());
      gBkg1SigmaParPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassParKolmPiPlusMuMinus, 1.);
      gBkg2SigmaParPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassParKolmPiPlusMuMinus, 2.);

      SumGraphs(gBkg1SigmaCombPiPlusMuMinus, gBkg1SigmaParPiPlusMuMinus, *gBkg1SigmaBufferPiPlusMuMinus);
      SumGraphs(gBkg2SigmaCombPiPlusMuMinus, gBkg2SigmaParPiPlusMuMinus, *gBkg2SigmaBufferPiPlusMuMinus);

      //REMOVE
      *gBkg1SigmaTotPiPlusMuMinus = *gBkg1SigmaBufferPiPlusMuMinus;
      *gBkg2SigmaTotPiPlusMuMinus = *gBkg2SigmaBufferPiPlusMuMinus;
    }

    if (an.Contains("Pos")) {

      // C - Prompt

      hInvMassPromptSRPiPlusMuMinus->Scale(counterPromptSBZPiPlusMuMinus/counterPromptSBPPiPlusMuMinus);
      cout<<"Prompt bkg - Kolmogorov test for SR sample: "<<hInvMassPromptSRPiPlusMuMinus->KolmogorovTest(hInvMassPromptPiPlusMuMinus)<<endl;
      hInvMassPromptKolmPiPlusMuMinus = (TH1D*)hInvMassPromptPiPlusMuMinus->Clone();
      hInvMassPromptKolmPiPlusMuMinus->SetName("hInvMassPromptKolm");
      hInvMassPromptKolmPiPlusMuMinus->Scale(hInvMassPromptSRPiPlusMuMinus->Integral()/hInvMassPromptKolmPiPlusMuMinus->Integral());
      gBkg1SigmaPromptPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassPromptKolmPiPlusMuMinus, 1.);
      gBkg2SigmaPromptPiPlusMuMinus = WindowScanNoFixedSigma(hInvMassPromptKolmPiPlusMuMinus, 2.);
	
      SumGraphs(gBkg1SigmaPromptPiPlusMuMinus, gBkg1SigmaBufferPiPlusMuMinus, *gBkg1SigmaTotPiPlusMuMinus);
      SumGraphs(gBkg2SigmaPromptPiPlusMuMinus, gBkg2SigmaBufferPiPlusMuMinus, *gBkg2SigmaTotPiPlusMuMinus);
      
      for (Int_t j = 0; j < 6; j++) 
	Par[j] = 0.;
      
      minSigma = 0.;
    }
  }

  // Saving histograms

  Save(path + "PiPlusMuMinus/SR/", c, hDistSRPiPlusMuMinus, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hCDASRPiPlusMuMinus, "Track-beamline CDA [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hTimeSRPiPlusMuMinus, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZSRPiPlusMuMinus, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hInvMassSRPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hMomPiSRPiPlusMuMinus, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hMomMuSRPiPlusMuMinus, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hDistvsMassSRPiPlusMuMinus, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hSRSRPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);

  Save(path + "PiPlusMuMinus/Comb/", c, hTimeCombPiPlusMuMinus, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombExtra1PiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombExtra2PiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombExtra3PiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombExtra4PiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombKolmPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombSRPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombSREnrichedPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hMomPiPosCombPiPlusMuMinus, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hMomMuNegCombPiPlusMuMinus, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3IDPiPosCombPiPlusMuMinus, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3IDMuNegCombPiPlusMuMinus, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3XYPiPosCombPiPlusMuMinus, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3XYMuNegCombPiPlusMuMinus, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hSRCombPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hSRFinalCombPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hSRFinalCombEnrichedPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, gBkg1SigmaCombPiPlusMuMinus, "gBkg1SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, gBkg2SigmaCombPiPlusMuMinus, "gBkg2SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  
  Save(path + "PiPlusMuMinus/Prompt/", c, hZPromptPiPlusMuMinus, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZPromptPPiPlusMuMinus, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPromptPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPromptSRPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPromptKolmPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hMomPiPromptPiPlusMuMinus, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hMomMuPromptPiPlusMuMinus, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZTimePromptPiPlusMuMinus, "Z coordinate of vertex [m]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZTimePromptPPiPlusMuMinus, "Z coordinate of vertex [m]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassVsZPromptPiPlusMuMinus, "Z coordinate of vertex [m]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassVsZq2PromptPiPlusMuMinus, "Z coordinate of vertex [m]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassVsPq2PromptPiPlusMuMinus, "Modulus of mother momentum [GeV/c]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassVsDistq2PromptPiPlusMuMinus, "Vertex-beam distance [mm]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hSRPromptPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hSRFinalPromptPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, gBkg1SigmaPromptPiPlusMuMinus, "gBkg1SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, gBkg2SigmaPromptPiPlusMuMinus, "gBkg2SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  Save(path + "PiPlusMuMinus/Par/", c, hDistParPiPlusMuMinus, "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassParPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassParSRPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassParKolmPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hMomPiParPiPlusMuMinus, "Pion momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hMomMuParPiPlusMuMinus, "Muon momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hDistvsMassParPiPlusMuMinus, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hDistvsMassParSRPiPlusMuMinus, "Reconstructed HNL mass", "Vertex-beamline distance [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hSRParPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hSRFinalParPiPlusMuMinus, "Z of CDA of mother wrt target-TAX line [m]", "CDA of mother wrt target-TAX line [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, gBkg1SigmaParPiPlusMuMinus, "gBkg1SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, gBkg2SigmaParPiPlusMuMinus, "gBkg2SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  if (!an.Contains("Pos")  && !histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data")) {
    Save(path + "PiPlusMuMinus/MC/", c, hInvMassMCPiPlusMuMinus, "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
    Save(path + "PiPlusMuMinus/MC/", c, hInvMassRecoVsMCPiPlusMuMinus, "MC true invariant mass [GeV/c^{2}]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
    Save(path + "PiPlusMuMinus/MC/", c, gMassMCPiPlusMuMinus, "gMassMC", "Reconstructed invariant mass [GeV/c^{2}]", "Sigma [GeV/c^{2}]", labelSize, titleSize);
  }

  tree->ResetBranchAddresses();
    
  return;
}

void BkgPiPlusMuMinus(TString dir, TString histo1, TString histo2) {

  // dir = output dir, histo1 = Data, histo2 = MC

  gErrorIgnoreLevel = kFatal;

  // Counters for pos/zero charge studies
  
  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  Double_t counterCombPiPlusMuMinus = 0;
  Double_t counterParPiPlusMuMinus = 0;
  Double_t counterPromptSBZPiPlusMuMinus = 0;
  Double_t counterPromptFVZPiPlusMuMinus = 0;
  Double_t counterPromptSBPPiPlusMuMinus = 0;
  Double_t counterPromptFVPPiPlusMuMinus = 0;
  Double_t PPiPlusMuMinus = 0.;
  TString pathPiPlusMuMinus = "/home/li/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/Zero/PiPlusMuMinus/";

  // Total expected bkg
  
  TGraph *gBkg1SigmaBufferPiPlusMuMinus = new TGraph();
  gBkg1SigmaBufferPiPlusMuMinus->SetNameTitle("gBkg1SigmaBuffer", "Total expected background");
  TGraph *gBkg2SigmaBufferPiPlusMuMinus = new TGraph();
  gBkg2SigmaBufferPiPlusMuMinus->SetNameTitle("gBkg2SigmaBuffer", "Total expected background");
  TGraph *gBkg1SigmaTotPiPlusMuMinus = new TGraph();
  gBkg1SigmaTotPiPlusMuMinus->SetNameTitle("gBkg1SigmaTot", "Total expected background");
  TGraph *gBkg2SigmaTotPiPlusMuMinus = new TGraph();
  gBkg2SigmaTotPiPlusMuMinus->SetNameTitle("gBkg2SigmaTot", "Total expected background");

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  Analyzer(dir, histo2, "HeavyNeutrino", c, counterCombPiPlusMuMinus, counterParPiPlusMuMinus, counterPromptSBZPiPlusMuMinus, counterPromptFVZPiPlusMuMinus, counterPromptSBPPiPlusMuMinus, counterPromptFVPPiPlusMuMinus, gBkg1SigmaBufferPiPlusMuMinus, gBkg2SigmaBufferPiPlusMuMinus, gBkg1SigmaTotPiPlusMuMinus, gBkg2SigmaTotPiPlusMuMinus);
  Analyzer(dir, histo1, "HeavyNeutrino", c, counterCombPiPlusMuMinus, counterParPiPlusMuMinus, counterPromptSBZPiPlusMuMinus, counterPromptFVZPiPlusMuMinus, counterPromptSBPPiPlusMuMinus, counterPromptFVPPiPlusMuMinus, gBkg1SigmaBufferPiPlusMuMinus, gBkg2SigmaBufferPiPlusMuMinus, gBkg1SigmaTotPiPlusMuMinus, gBkg2SigmaTotPiPlusMuMinus);
  Analyzer(dir, histo1, "HeavyNeutrinoPos", c, counterCombPiPlusMuMinus, counterParPiPlusMuMinus, counterPromptSBZPiPlusMuMinus, counterPromptFVZPiPlusMuMinus, counterPromptSBPPiPlusMuMinus, counterPromptFVPPiPlusMuMinus, gBkg1SigmaBufferPiPlusMuMinus, gBkg2SigmaBufferPiPlusMuMinus, gBkg1SigmaTotPiPlusMuMinus, gBkg2SigmaTotPiPlusMuMinus);
  
  Save(pathPiPlusMuMinus + "Total/", c, gBkg1SigmaTotPiPlusMuMinus, "gBkg1SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(pathPiPlusMuMinus + "Total/", c, gBkg2SigmaTotPiPlusMuMinus, "gBkg2SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  if (counterPromptSBPPiPlusMuMinus != 0.)
    PPiPlusMuMinus = counterPromptFVPPiPlusMuMinus*counterPromptSBZPiPlusMuMinus/counterPromptSBPPiPlusMuMinus;
  
  counterPromptFVZPiPlusMuMinus = PPiPlusMuMinus;
  
  cout<<"PI+MU-: Number of events in SR and time sidebands: "<<counterCombPiPlusMuMinus<<", and beamdist sidebands: "<<counterParPiPlusMuMinus<<", and Z sidebands (0-charge): "<<counterPromptSBZPiPlusMuMinus<<", and Z sidebands (pos-charge): "<<counterPromptSBPPiPlusMuMinus<<", and FV (pos-charge): "<<counterPromptFVPPiPlusMuMinus<<", and FV (0-charge): "<<counterPromptFVZPiPlusMuMinus<<endl;

  _exit(0);
}
