// Variables for sigma vs mass fit

Double_t Par[6] = {0., 0., 0., 0., 0., 0.};
Double_t minSigma = 0.;
Double_t massMin = 0.25;
Double_t massMinForBkg = 0.5;
Double_t massMax = 1.96;
Double_t massStep = 0.01;
Int_t nMassBins = (massMax - massMin)/massStep;
Int_t nMassBinsForBkg = (massMax - massMinForBkg)/massStep;
Double_t binWidth = (massMax - massMin)/(nMassBins);

void Save(TString path, TCanvas *c, TH1D* h, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();
  
  h->Draw("hist");
  h->GetXaxis()->SetTitle(x);
  h->GetYaxis()->SetTitle(y);
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
  h->Write();
  c->SaveAs(path + h->GetName() + ".pdf");
  c->SaveAs(path + h->GetName() + ".png");

  return;
}

void Save(TString path, TCanvas *c, TH1D* h, TString x, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();
  
  h->Draw("hist");
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

  return;
}

void Save(TString path, TCanvas *c, TH2D* h, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  TString name = h->GetName();

  h->Draw("colz");

  if (h->GetEntries() < 100)
    gStyle->SetPalette(1);
  else
    gStyle->SetPalette(kBird);

  if(name.Contains("InvMassReco")) {
    h->GetXaxis()->SetRangeUser(0.5, 2.2);
    h->GetYaxis()->SetRangeUser(0.5, 2.2);
  }
  
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

  if (name.Contains("gMassMC")) {
    g->Draw("AP*");
    g->GetXaxis()->SetRangeUser(0.5, 2.2);
    g->GetYaxis()->SetRangeUser(0.001, 0.0075);
  }
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

void Save(TString path, TCanvas *c, TGraphAsymmErrors* g, TString name, TString x, TString y, Double_t labelSize, Double_t titleSize) {

  if (name.Contains("gMassMC"))
    g->Draw("AP*");
  else {
    g->SetFillStyle(3001);
    g->SetFillColor(kBlue);
    g->Draw("AL3");
  }
  
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

// Scan on HNL mass for expected number of background events, with window around mass hypothesis at 1 or 2 sigma (non-fixed sigma: each MC mass hypothesis is fitted, sigma is plotted vs mass and fitted. Fit parameters are used to compute sigma for each mass hypothesis in this scan)

TGraphAsymmErrors* WindowScanNoFixedSigma(TH1D* h, Int_t howManySigmas) {

  Int_t firstBin = h->FindFirstBinAbove(0.,1);
  Double_t firstBinValue = h->GetBinCenter(firstBin)-h->GetBinWidth(firstBin)/2.;
  Int_t lastBin = h->FindLastBinAbove(0.,1);
  Double_t lastBinValue = h->GetBinCenter(lastBin)-h->GetBinWidth(lastBin)/2.;
  Int_t counter = 0;
  TGraphAsymmErrors *res = new TGraphAsymmErrors();
  TString name = h->GetName();
  
  for (Double_t massHyp = massMinForBkg-minSigma; massHyp <= massMax+minSigma; massHyp += minSigma) {
    Double_t sigma = Par[5]*TMath::Power(massHyp, 5.) + Par[4]*TMath::Power(massHyp, 4.) + Par[3]*TMath::Power(massHyp, 3.) + Par[2]*TMath::Power(massHyp, 2.) + Par[1]*TMath::Power(massHyp, 1.) + Par[0];
    TAxis *axis = h->GetXaxis();
    Int_t bmin = axis->FindBin(massHyp-howManySigmas*sigma);
    Int_t bmax = axis->FindBin(massHyp+howManySigmas*sigma);
    Double_t nBkg = h->Integral(bmin, bmax);
    nBkg -= h->GetBinContent(bmin)*(massHyp-howManySigmas*sigma-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    nBkg -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-(massHyp+howManySigmas*sigma))/axis->GetBinWidth(bmax);
    res->SetPoint(counter, massHyp, nBkg);

    if (name.Contains("omb") || name.Contains("otal"))
      res->SetPointError(counter, 0., 0., nBkg*0.5, nBkg*0.5);
    else
      res->SetPointError(counter, 0., 0., 0., 0.);

    counter++;
  }

  return res;
}

// Sum graphs for total expected bkg

void SumGraphs(TGraphAsymmErrors* g1, TGraphAsymmErrors* g2, TGraphAsymmErrors &g) {

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
    g.SetPointError(i, 0., 0., g1->GetErrorYlow(i)+g2->GetErrorYlow(i), g1->GetErrorYhigh(i)+g2->GetErrorYhigh(i));
    X = 0.;
    Y = 0.;
    Ytot = 0.;
  }

  return;
}

void Analyzer(TString dir, TString histo1, TString an, TCanvas* c, TGraphAsymmErrors* gBkg1SigmaBuffer, TGraphAsymmErrors* gBkg2SigmaBuffer, TGraphAsymmErrors* gBkg1SigmaTot, TGraphAsymmErrors* gBkg2SigmaTot, TFile *of) {

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
      path = "/mnt/c/Users/loren/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/" + analyzer + "/";
    else
      path = "/mnt/c/Users/loren/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/MC/" + analyzer + "/";
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
    
  TH1D *hInvMassComb = new TH1D("hInvMassComb", "Events in sidebands", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombKolm = new TH1D("hInvMassCombKolm", "Events in sidebands", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassCombSR = new TH1D("hInvMassCombSR", "Events in sidebands and blinded region", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);
  
  TH1D *hMomPiPosComb = new TH1D("hMomPiPosComb", "Events in sidebands", 100, 0., 200.);
  TH1D *hMomMuNegComb = new TH1D("hMomMuNegComb", "Events in sidebands", 100, 0., 200.);
  TH2D *hGTK3XYPiPosComb = new TH2D("hGTK3XYPiPosComb", "Events in sidebands", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3XYMuNegComb = new TH2D("hGTK3XYMuNegComb", "Events in sidebands", 100, -50., 50., 100, -50., 50.);
  TH2D *hGTK3IDPiPosComb = new TH2D("hGTK3IDPiPosComb", "Events in sidebands", 6, -0.5, 5.5, 100, 0., 250.);
  TH2D *hGTK3IDMuNegComb = new TH2D("hGTK3IDMuNegComb", "Events in sidebands", 6, -0.5, 5.5, 100, 0., 250.);

  TGraphAsymmErrors *gBkg1SigmaComb = new TGraphAsymmErrors();
  gBkg1SigmaComb->SetNameTitle("PiPlusMuMinus/gBkg1SigmaComb", "Expected combinatorial background events");
  TGraphAsymmErrors *gBkg2SigmaComb = new TGraphAsymmErrors();			
  gBkg2SigmaComb->SetNameTitle("PiPlusMuMinus/gBkg2SigmaComb", "Expected combinatorial background events");

  // Prompt - pi+mu-

  TH1D *hInvMassPrompt = new TH1D("hInvMassPrompt", "Events in sidebands", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptSR = new TH1D("hInvMassPromptSR", "Events in sidebands and blinded region", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassPromptKolm = new TH1D("hInvMassPromptKolm", "Events in sidebands", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2.);

  TH2D *hZdT = new TH2D("hZdT", "Events in sideband and BR", 100, -20., 20., 200, 90., 190.);
  TH2D *hZR = new TH2D("hZR", "Events in sideband and BR", 200, 100., 800., 200, 90., 190.);
  TH2D *hZM = new TH2D("hZM", "Events in sideband and BR", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2., 200, 90., 190.);
  TH2D *hdTM = new TH2D("hdTM", "Events in sideband and BR", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2., 100, -20., 20.);
  TH2D *hRM = new TH2D("hRM", "Events in sideband and BR", nMassBinsForBkg+1, massMinForBkg-binWidth/2., massMax+binWidth/2., 200, 100., 800.);
  TH2D *hdTR = new TH2D("hdTR", "Events in sideband and BR", 200, 100., 800., 100, -20., 20.);
  TGraphAsymmErrors *gBkg1SigmaPrompt = new TGraphAsymmErrors();
  gBkg1SigmaPrompt->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPrompt", "Expected muon-induced background events");
  TGraphAsymmErrors *gBkg2SigmaPrompt = new TGraphAsymmErrors();			   
  gBkg2SigmaPrompt->SetNameTitle("PiPlusMuMinus/gBkg2SigmaPrompt", "Expected muon-induced background events");

  // Parasitic - pi+mu-
  /*
  TH1D *hInvMassPar = new TH1D("hInvMassPar", "Events in sidebands", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParSR = new TH1D("hInvMassParSR", "Events in sidebands and blinded region", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TH1D *hInvMassParKolm = new TH1D("hInvMassParKolm", "Events in sidebands", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);

  TGraphAsymmErrors *gBkg1SigmaPar = new TGraphAsymmErrors();
  gBkg1SigmaPar->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPar", "Expected kaon-induced background events");
  TGraphAsymmErrors *gBkg2SigmaPar = new TGraphAsymmErrors();			      
  gBkg2SigmaPar->SetNameTitle("PiPlusMuMinus/gBkg2SigmaPar", "Expected kaon-induced background events");
  */
  // SR - pi+mu-

  TH1D *hCDASR = new TH1D("hCDASR", "Events outside blinded region", 100, 0., 1000.);
  TH1D *hTimeSR = new TH1D("hTimeSR", "Events outside blinded region", 100, -15., 15.);
  TH1D *hZSR = new TH1D("hZSR", "Events outside blinded region", 500, 100., 190.);
  TH2D *hZTimeSR = new TH2D("hZTimeSR", "Events in sidebands", 500, 100., 190., 100, -15., 15.);
  TH1D *hZTimeSR1 = new TH1D("hZTimeSR1", "Events in sidebands", 500, 100., 190.);
  TH1D *hZTimeSR2 = new TH1D("hZTimeSR2", "Events in sidebands", 100, -15., 15.);
  TH1D *hZTimeSR4 = new TH1D("hZTimeSR4", "Events in sidebands", 500, 100., 190.);
  TH2D *hZTimeSRBR = new TH2D("hZTimeSRBR", "Events in sidebands", 500, 100., 190., 100, -15., 15.);
  TH1D *hZTimeSRBR1 = new TH1D("hZTimeSRBR1", "Events in sidebands", 500, 100., 190.);
  TH1D *hZTimeSRBR4 = new TH1D("hZTimeSRBR4", "Events in sidebands", 500, 100., 190.);
  TH2D *hZRadSR = new TH2D("hZRadSR", "Events in sidebands", 150, 100., 190., 150, 90., 1000.);
  TH2D *hZRadSRBR = new TH2D("hZRadSRBR", "Events in sidebands", 150, 100., 190., 150, 90., 1000.);

  // MC only - pi+mu-

  TH1D *hMass[172];
  for (Int_t i = 0; i < 172; i++) {
    hMass[i] = new TH1D(Form("hMass%d", i), Form("hMass%d", i), 100, (i*10.+245.)/1000., (i*10.+255.)/1000.);
  }
  TH2D *hInvMassRecoVsMC = new TH2D("hInvMassRecoVsMC", "Mass resolution studies", nMassBins+1, massMin-binWidth/2., massMax+binWidth/2., nMassBins+1, massMin-binWidth/2., massMax+binWidth/2.);
  TGraph *gMassMC = new TGraph();
  gMassMC->SetNameTitle("PiPlusMuMinus/gMassMC", "Mass resolution studies");

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
    Bool_t Comb = (((CHODTime1-CHODTime2 < Time1Max && CHODTime1-CHODTime2 > Time1Min) || (CHODTime1-CHODTime2 > Time2Min && CHODTime1-CHODTime2 < Time2Max)));
    Bool_t Prompt = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t SB = (Zvertex >= ZVertexMin && Zvertex <= ZVertexMax);
    Bool_t FV = (Zvertex >= ZVertexMax && Zvertex <= ZEndFV);
    Bool_t Par = (BeamCDA1 <= 50. && BeamCDA2 <= 50.);
    Bool_t PiPlusMuMinus = ((Assoc == 2 && Charge1 == 1) || (Assoc == 1 && Charge1 == -1));
    Bool_t Spike1 = ((Assoc == 1 && Charge1 == -1) && (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) && (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax));
    Bool_t Spike2 = ((Assoc == 2 && Charge1 == 1) && (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) && (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax));
    Bool_t noSpike = (!Spike1 && !Spike2);
    Bool_t LAVcomp = (BeamlineDist > 500.);
    //Bool_t Kcomp = (Zvertex > 107000. && Zvertex < 130000. && BeamlineDist > 200. && BeamlineDist < 300.);
    Bool_t noExtra = (!LAVcomp);// && !Kcomp);
    
    // MC studies on reco mass vs true mass
    
    if (MC && Zero && noSpike && noExtra && PiPlusMuMinus) {
      hInvMassRecoVsMC->Fill(trueMass, invMass/1000.);
      Int_t massIndex = std::round((trueMass*1000-250)/10);
      hMass[massIndex]->Fill(invMass/1000.);
    }
    
    // 1 - Events outside blinded region, to show SB for all bkgs
    
    if (outsideSR && Zero && noSpike && noExtra && CDAIn && PiPlusMuMinus) { // all events outside blinded region
      hCDASR->Fill(BeamCDA1, Weight);
      hCDASR->Fill(BeamCDA2, Weight);
      hTimeSR->Fill(CHODTime1-CHODTime2, Weight);
    }
    if (outsideSR && Zero && noSpike && noExtra && CDAIn && PiPlusMuMinus) {
      hZSR->Fill(Zvertex/1000., Weight);
      hZTimeSR->Fill(Zvertex/1000., CHODTime1-CHODTime2, Weight);
      if (Zvertex > 107000. && Zvertex < 130000. && BeamlineDist > 200. && BeamlineDist < 300.)
	hZTimeSR2->Fill(CHODTime1-CHODTime2, Weight);
      if ((CHODTime1-CHODTime2) > -2 && (CHODTime1-CHODTime2) < 2) {
	hZTimeSR1->Fill(Zvertex/1000., Weight);
	hZRadSR->Fill(Zvertex/1000., BeamlineDist, Weight);
      }
      if (((CHODTime1-CHODTime2) > 2.5 && (CHODTime1-CHODTime2) < 5.) || ((CHODTime1-CHODTime2) > -5. && (CHODTime1-CHODTime2) < -2.5))
	hZTimeSR4->Fill(Zvertex/1000., Weight);
    }

    // 2 - Combinatorial SB
    
    if (Comb && Zero && CDAIn && noExtra && PiPlusMuMinus) { // all events inside time sidebands
      
      // Study spike at 75 GeV
      
      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	hMomMuNegComb->Fill(Mom1->Mag()/1000., Weight);
	hMomPiPosComb->Fill(Mom2->Mag()/1000., Weight);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	hMomPiPosComb->Fill(Mom1->Mag()/1000., Weight);
	hMomMuNegComb->Fill(Mom2->Mag()/1000., Weight);
      }
      
      // Select momentum in [70-80] GeV to study spike
      
      if (Assoc == 1 && Charge1 == -1) { // mu-pi+
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYMuNegComb->Fill(xGTK31, yGTK31, Weight);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYPiPosComb->Fill(xGTK32, yGTK32, Weight);
      }
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.)
	  hGTK3XYPiPosComb->Fill(xGTK31, yGTK31, Weight);
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.)
	  hGTK3XYMuNegComb->Fill(xGTK32, yGTK32, Weight);
      }

      // Select momentum in [70-80] GeV and xy @ GTK3 to study spike
      
      if (RICHInfo1->at(0) == 99)
	RICHInfo1->at(0) = 5;

      if (RICHInfo2->at(0) == 99)
	RICHInfo2->at(0) = 5;

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
      else if (Assoc == 2 && Charge1 == 1) { // pi+mu-
	if (Mom1->Mag()/1000. >= 70. && Mom1->Mag()/1000. <= 80.) {
	  if (xGTK31 >= xGTKMin && xGTK31 <= xGTKMax && yGTK31 >= yGTKMin && yGTK31 <= yGTKMax)
	    hGTK3IDPiPosComb->Fill(RICHInfo1->at(0), RICHInfo1->at(5));
	}
	if (Mom2->Mag()/1000. >= 70. && Mom2->Mag()/1000. <= 80.) {
	  if (xGTK32 >= xGTKMin && xGTK32 <= xGTKMax && yGTK32 >= yGTKMin && yGTK32 <= yGTKMax)	  
	    hGTK3IDMuNegComb->Fill(RICHInfo2->at(0), RICHInfo2->at(5));
	}
      }
    
      if (noSpike) { // all events inside time sidebands and no spike
	hInvMassComb->Fill(invMass/1000., Weight);            
      }
    }
    
    // 3 - Prompt SB
    
    if (FV && Zero && noSpike && noExtra && CDAIn && PiPlusMuMinus) {
      hInvMassPrompt->Fill(invMass/1000., Weight);
    }
    
    // 4 - Parasitic SB
    /*    
    if (Par && Zero && noSpike && noExtra && CDAIn && PiPlusMuMinus) { // all events inside beam distance sidebands
      hInvMassPar->Fill(invMass/1000., Weight);
    }
    */
    // 5 - Events inside blinded region
    
    if (!outsideSR && noSpike && noExtra) { // all events inside blinded region...
      if (Comb && Zero && !MC && CDAIn && PiPlusMuMinus) { // ...and time sidebands (not enriched sample to study spike at 75 GeV)
	hInvMassCombSR->Fill(invMass/1000., Weight);
      }
      /*
      if (Par && Zero && !MC && CDAIn && PiPlusMuMinus) { // ...and beam distance sidebands
	hInvMassParSR->Fill(invMass/1000., Weight);
      }
      */
      if (SB && Zero && CDAIn && PiPlusMuMinus) {
	hInvMassPromptSR->Fill(invMass/1000., Weight);		
	hZTimeSRBR->Fill(Zvertex/1000., CHODTime1-CHODTime2, Weight);
	hZdT->Fill(CHODTime1-CHODTime2, Zvertex/1000., Weight);
	hZR->Fill(BeamlineDist, Zvertex/1000., Weight);
	hZM->Fill(invMass/1000., Zvertex/1000., Weight);
	hdTM->Fill(invMass/1000., CHODTime1-CHODTime2, Weight);
	hdTR->Fill(BeamlineDist, CHODTime1-CHODTime2, Weight);
	hRM->Fill(invMass/1000., BeamlineDist, Weight);
	if ((CHODTime1-CHODTime2) > -2 && (CHODTime1-CHODTime2) < 2) {
	  hZTimeSRBR1->Fill(Zvertex/1000., Weight);
	  hZRadSRBR->Fill(Zvertex/1000., BeamlineDist, Weight);
	}
	if (((CHODTime1-CHODTime2) > 2.5 && (CHODTime1-CHODTime2) < 5.) || ((CHODTime1-CHODTime2) > -5. && (CHODTime1-CHODTime2) < -2.5))
	  hZTimeSRBR4->Fill(Zvertex/1000., Weight);
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
      TF1 *f2 = new TF1("f2", "gaus", (i*10.+249.)/1000., (i*10.+251.)/1000.);
      if ((i*10.+245.)/1000. < 1.5) {
	hMass[i]->Fit("f1", "Rq");
	sigma = f1->GetParameter(2);
      }
      else {
	hMass[i]->Fit("f2", "Rq");
	sigma = f2->GetParameter(2);
      }
      gMassMC->SetPoint(i, (i*10.+250.)/1000., sigma);  
      if (sigma < sigmaMin && sigma != 0.)                 
        sigmaMin = sigma;                                  
    }
    minSigma = sigmaMin;
    gMassMC->Fit("pol1");
    TF1 *f = gMassMC->GetFunction("pol1");

    for (Int_t j = 0; j < 2; j++) {
      Par[j] = f->GetParameter(j);
    }
  }

  // Extrapolation of n bkg evt for prompt bkg

  Double_t N = 0.;
  
  if (Zero && !MC) {
    TH1D *h1 = new TH1D("h1", "Events per bin (size=0.18m)", 5, -0.5, 4.5);
    for (Int_t i = 1; i < hZTimeSR1->GetNbinsX(); i++){
      if (hZTimeSR1->GetBinCenter(i) > 110. && hZTimeSR1->GetBinCenter(i) < 120.) {
	h1->Fill(hZTimeSR1->GetBinContent(i));
      }
    }
    
    TF1 *p = new TF1("p", "[0]*TMath::Poisson(x, [1])");
    p->SetParameter(0, 1.);
    p->SetParameter(1, 1.);
    h1->Fit("p", "q");

    Double_t mean = p->GetParameter(1);
    Double_t conv = 1./hZTimeSR1->GetBinWidth(1);
    Double_t nSB = mean*conv*20.;
    Double_t nFV = mean*conv*60.;
    Double_t ratio = nSB/nFV;
    Double_t n = 2.3;

    N = n*nFV/nSB;
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
      gBkg1SigmaComb->SetNameTitle("PiPlusMuMinus/gBkg1SigmaComb", "Expected combinatorial background events");
      gBkg2SigmaComb->SetNameTitle("PiPlusMuMinus/gBkg1SigmaComb", "Expected combinatorial background events");
      
      // B - Parasitic
      /*
      cout<<"Parasitic bkg - Kolmogorov test for SR sample: "<<hInvMassParSR->KolmogorovTest(hInvMassPar)<<endl;
      hInvMassParKolm = (TH1D*)hInvMassPar->Clone();
      hInvMassParKolm->SetName("hInvMassParKolm");
      hInvMassParKolm->Scale(hInvMassParSR->Integral()/hInvMassParKolm->Integral());
      gBkg1SigmaPar = WindowScanNoFixedSigma(hInvMassParKolm, 1.);
      gBkg2SigmaPar = WindowScanNoFixedSigma(hInvMassParKolm, 2.);
      gBkg1SigmaPar->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPar", "Expected kaon-induced background events");
      gBkg2SigmaPar->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPar", "Expected kaon-induced background events");

      SumGraphs(gBkg1SigmaComb, gBkg1SigmaPar, *gBkg1SigmaBuffer);
      SumGraphs(gBkg2SigmaComb, gBkg2SigmaPar, *gBkg2SigmaBuffer);
      */
    }

    if (Zero) {

      // C - Prompt

      cout<<"Prompt bkg - Kolmogorov test for SR sample: "<<hInvMassPromptSR->KolmogorovTest(hInvMassPrompt)<<endl;
      hInvMassPromptKolm = (TH1D*)hInvMassPrompt->Clone();
      hInvMassPromptKolm->SetName("hInvMassPromptKolm");
      hInvMassPromptKolm->Scale(N/hInvMassPromptKolm->Integral());
      gBkg1SigmaPrompt = WindowScanNoFixedSigma(hInvMassPromptKolm, 1.);
      gBkg2SigmaPrompt = WindowScanNoFixedSigma(hInvMassPromptKolm, 2.);
      gBkg1SigmaPrompt->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPrompt", "Expected muon-induced background events");
      gBkg2SigmaPrompt->SetNameTitle("PiPlusMuMinus/gBkg1SigmaPrompt", "Expected muon-induced background events");
	
      SumGraphs(gBkg1SigmaComb, gBkg1SigmaPrompt, *gBkg1SigmaTot);
      SumGraphs(gBkg2SigmaComb, gBkg2SigmaPrompt, *gBkg2SigmaTot);
      
      for (Int_t j = 0; j < 2; j++) 
	Par[j] = 0.;
      
      minSigma = 0.;
    }
  }

  // Saving histograms
  
  Save(path + "PiPlusMuMinus/SR/", c, hCDASR, "Track-beamline CDA [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hTimeSR, "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZSR, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSR, "Z coordinate of vertex [m]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSR1, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSRBR1, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSR2, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSR4, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZTimeSRBR4, "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZRadSR, "Z coordinate of vertex [m]", "Vertex radial position [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/SR/", c, hZRadSRBR, "Z coordinate of vertex [m]", "Vertex radial position [mm]", labelSize, titleSize);

  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassComb, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombKolm, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hInvMassCombSR, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hMomPiPosComb, "Track momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hMomMuNegComb, "Track momentum [GeV/c]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3IDPiPosComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3IDMuNegComb, "RICH hypothesis", "RICH radius [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3XYPiPosComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, hGTK3XYMuNegComb, "X at GTK3 [mm]", "Y at GTK3 [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, gBkg1SigmaComb, "gBkg1SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Comb/", c, gBkg2SigmaComb, "gBkg2SigmaComb", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPrompt, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPromptSR, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hInvMassPromptKolm, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZdT, "Track time difference [ns]", "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZR, "Vertex radial position [mm]", "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hZM, "Reconstructed invariant mass [GeV/c^{2}]", "Z coordinate of vertex [m]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hdTM, "Reconstructed invariant mass [GeV/c^{2}]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hdTR, "Vertex radial position [mm]", "Track time difference [ns]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, hRM, "Reconstructed invariant mass [GeV/c^{2}]", "Vertex radial position [mm]", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, gBkg1SigmaPrompt, "gBkg1SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Prompt/", c, gBkg2SigmaPrompt, "gBkg2SigmaPrompt", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  /*
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassPar, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassParSR, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, hInvMassParKolm, "Reconstructed invariant mass [GeV/c^{2}]", "Events/10 MeV", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, gBkg1SigmaPar, "gBkg1SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "PiPlusMuMinus/Par/", c, gBkg2SigmaPar, "gBkg2SigmaPar", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  */
  if (!an.Contains("Pos")  && !histo1.Contains("2016") && !histo1.Contains("2017") && !histo1.Contains("2018") && !histo1.Contains("Data")) {
    Save(path + "PiPlusMuMinus/MC/", c, hInvMassRecoVsMC, "N true invariant mass [GeV/c^{2}]", "Reconstructed invariant mass [GeV/c^{2}]", labelSize, titleSize);
    Save(path + "PiPlusMuMinus/MC/", c, gMassMC, "gMassMC", "Reconstructed invariant mass [GeV/c^{2}]", "Mass resolution [GeV/c^{2}]", labelSize, titleSize);
  }
  /*  
  if (analyzer.Contains("Zero")) {
    if (MC)
      of->cd("Zero/MC/");
    else
      of->cd("Zero/Data/");
    
    hZTimeSR->Write();
    hZTimeSR1->Write();
    hZTimeSR4->Write();
    hZRadSR->Write();
    hZTimeSRBR->Write();
    hZTimeSRBR1->Write();
    hZTimeSRBR4->Write();
    hZRadSRBR->Write();
  }
  */
  tree->ResetBranchAddresses();
    
  return;
}

void BkgPiPlusMuMinus(TString dir, TString histo1, TString histo2) {

  // dir = output dir, histo1 = Data, histo2 = MC, histo3 = control

  gErrorIgnoreLevel = kFatal;

  // Counters for pos/zero charge studies
  
  TCanvas *c = new TCanvas();  
  Double_t labelSize = 0.05;
  Double_t titleSize = 0.07;  
  TString path = "/mnt/c/Users/loren/cernbox/PhD/TalksAndPapers/Notes/MCnote/images/Plots/Data/All/Zero/PiPlusMuMinus/";

  TFile *of = new TFile("/home/lorenza/GitVarious/MyScripts/Others/AnalyzerPlots/bkg.root", "RECREATE");
  /*
  of->cd();
  gDirectory->mkdir("Zero");
  of->cd("Zero");
  gDirectory->mkdir("MC");
  gDirectory->mkdir("Data");
  of->cd();
  */
  // Total expected bkg
  
  TGraphAsymmErrors *gBkg1SigmaBuffer = new TGraphAsymmErrors();
  gBkg1SigmaBuffer->SetNameTitle("gBkg1SigmaBuffer", "Total expected background");
  TGraphAsymmErrors *gBkg2SigmaBuffer = new TGraphAsymmErrors();
  gBkg2SigmaBuffer->SetNameTitle("gBkg2SigmaBuffer", "Total expected background");
  TGraphAsymmErrors *gBkg1SigmaTot = new TGraphAsymmErrors();
  gBkg1SigmaTot->SetNameTitle("gBkg1SigmaTot", "Total expected background");
  TGraphAsymmErrors *gBkg2SigmaTot = new TGraphAsymmErrors();
  gBkg2SigmaTot->SetNameTitle("gBkg2SigmaTot", "Total expected background");

  c->SetRightMargin(0.2);
  c->SetLeftMargin(0.2);
  c->SetBottomMargin(0.25);
  c->SetTopMargin(0.15);
  c->SetGrid();
  c->RedrawAxis();

  Analyzer(dir, histo2, "HeavyNeutrino", c, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot, of);
  Analyzer(dir, histo1, "HeavyNeutrino", c, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot, of);
  //Analyzer(dir, histo1, "HeavyNeutrinoPos", c, gBkg1SigmaBuffer, gBkg2SigmaBuffer, gBkg1SigmaTot, gBkg2SigmaTot, of);
  
  Save(path + "Total/", c, gBkg1SigmaTot, "gBkg1SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);
  Save(path + "Total/", c, gBkg2SigmaTot, "gBkg2SigmaTot", "N mass [GeV/c^{2}]", "N_{exp}", labelSize, titleSize);

  of->Close();
  _exit(0);
}
